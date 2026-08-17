package math.linalg;

import java.util.Random;
import java.util.stream.IntStream;

/**
 * Principal Component Analysis restricted to the leading components, for the case where
 * {@code no_dims} is much smaller than the number of features.
 *
 * <p>{@link JacobiPCA} decomposes the whole matrix and therefore costs {@code O(m n^2)} for the QR
 * plus {@code O(sweeps n^3)} for the Jacobi sweeps, computing all {@code n} components even when only
 * two or three are wanted. This class instead runs a randomized subspace iteration: it only ever multiplies
 * the data by a thin {@code n x l} block, with {@code l = no_dims + oversampling}, so the cost is
 * {@code O(iterations * m * n * l)} and no {@code n x n} matrix is ever formed. The centered matrix
 * is not materialized either.</p>
 *
 * <p>Input convention, output convention and argument validation are identical to
 * {@link JacobiPCA}: {@code matrix} is row-major with one sample per row, the returned array holds
 * the centered samples projected onto the leading {@code no_dims} principal directions, and each
 * component's sign is canonicalized so that its entry of largest magnitude is positive. The two
 * classes are therefore interchangeable at a call site.</p>
 *
 * <p>The result is an approximation. The iteration stops once the leading {@code no_dims} Ritz
 * values are stable to {@code tolerance} relative to the largest one; if that does not happen within
 * {@code maxIterations}, {@link #converged()} returns {@code false} and the caller should fall back
 * to the exact {@link JacobiPCA}. Where the spectrum is degenerate the Ritz values still converge
 * while the directions inside the tied subspace remain arbitrary - that is a property of the data,
 * not of this method, and an exact solver picks an equally arbitrary direction there.</p>
 *
 * <p>Output is deterministic: the starting block is drawn from a fixed seed and every accumulation
 * runs in a fixed order within a single thread, so parallel execution does not change the result.</p>
 */
public final class TruncatedPCA {

    private static final int    DEFAULT_OVERSAMPLING   = 6;
    private static final double DEFAULT_TOLERANCE      = 1e-9;
    private static final int    DEFAULT_MAX_ITERATIONS = 30;
    private static final long   DEFAULT_SEED           = 20260815L;

    /** Go parallel from roughly 65k elements, matching {@link FlatParallelJacobiSVD}. */
    private static final long PARALLEL_THRESHOLD = 1L << 16;

    private final int    oversampling;
    private final double tolerance;
    private final int    maxIterations;
    private final long   seed;

    // fitted state
    private double[]   mean;            // per-feature mean, length n
    private double[][] components;      // [no_dims][n], row k = k-th principal direction
    private double[]   singularValues;  // [no_dims], descending
    private boolean    converged;
    private int        iterations;
    private int m, n;                   // samples, features of the last fit

    public TruncatedPCA() {
        this(DEFAULT_OVERSAMPLING, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATIONS, DEFAULT_SEED);
    }

    /**
     * @param oversampling  extra dimensions carried along beyond {@code no_dims}; more of them make
     *                      the leading components converge in fewer iterations
     * @param tolerance     relative tolerance for the stability test of the leading Ritz values
     * @param maxIterations hard cap; reaching it leaves {@link #converged()} at {@code false}
     * @param seed          seed of the starting block, fixed so that results are reproducible
     */
    public TruncatedPCA(int oversampling, double tolerance, int maxIterations, long seed) {
        if (oversampling < 0) throw new IllegalArgumentException("oversampling must be >= 0");
        if (!(tolerance >= 0.0)) throw new IllegalArgumentException("tolerance must be >= 0");
        if (maxIterations < 1) throw new IllegalArgumentException("maxIterations must be >= 1");
        this.oversampling = oversampling;
        this.tolerance = tolerance;
        this.maxIterations = maxIterations;
        this.seed = seed;
    }

    /**
     * An instance that runs a fixed number of subspace iterations instead of testing the Ritz values
     * for stability.
     *
     * <p>Use this when many components are wanted. The stability test asks every kept direction to
     * settle, and the tail of a spectrum with dozens of components is nearly tied, so the test is
     * never satisfied and only the iteration cap ever ends the loop - the test then costs work
     * without deciding anything. A fixed count says the same thing directly. It is also the right
     * answer when the caller needs the retained <em>subspace</em> rather than the individual
     * directions inside it, because a nearly tied tail leaves those directions undetermined anyway
     * while the subspace they span is already accurate.</p>
     *
     * <p>{@link #converged()} then reports whether the Ritz values stopped moving altogether before
     * the count was exhausted, which is a stronger statement than the tolerance test and is not
     * normally reached; a caller of this factory has already decided that {@code iterations} is
     * enough and does not need to consult it.</p>
     *
     * @param oversampling extra dimensions carried along beyond {@code no_dims}
     * @param iterations   number of subspace iterations to run
     * @return a {@code TruncatedPCA} that always runs exactly {@code iterations} iterations
     */
    public static TruncatedPCA fixedIterations(int oversampling, int iterations) {
        return new TruncatedPCA(oversampling, 0.0, iterations, DEFAULT_SEED);
    }

    /**
     * Fits on {@code matrix} (row-major, row = sample) and returns the samples projected onto the
     * first {@code no_dims} principal components. Check {@link #converged()} afterwards.
     *
     * @param matrix  row-major data, m samples x n features (rectangular, non-empty)
     * @param no_dims number of components to keep, 1 &lt;= no_dims &lt;= min(m, n)
     * @return {@code double[m][no_dims]} of projected coordinates
     */
    public double[][] pca(double[][] matrix, int no_dims) {
        validate(matrix, no_dims);

        m = matrix.length;
        n = matrix[0].length;
        converged = false;
        iterations = 0;

        // 1) per-feature mean; the centered matrix stays implicit
        mean = new double[n];
        for (double[] row : matrix)
            for (int j = 0; j < n; j++) mean[j] += row[j];
        for (int j = 0; j < n; j++) mean[j] /= m;

        final int l = Math.min(Math.min(n, m), no_dims + oversampling);
        final boolean parallel = (long) m * n >= PARALLEL_THRESHOLD;

        // 2) starting block, l x n, row k = one basis vector of the search subspace
        double[][] v = startingBlock(l, n);
        orthonormalizeRows(v);

        // 3) subspace iteration: v <- Ac^T (Ac v^T), with the Ritz values of Ac v^T as the
        //    convergence indicator. y and v are always consistent when the loop is left.
        double[][] y = null;
        double[] previous = null;
        FlatParallelJacobiSVD.Result ritz = null;

        for (int it = 1; it <= maxIterations; it++) {
            y = project(matrix, v, parallel);                  // m x l
            ritz = decomposeGram(y, l);                        // eigen decomposition of y^T y
            iterations = it;

            double[] current = ritz.sigma;
            if (previous != null && resolved(previous, current, no_dims)) {
                converged = true;
                break;
            }
            previous = current.clone();
            if (it == maxIterations) break;

            v = backProject(matrix, y, parallel);              // l x n
            orthonormalizeRows(v);
        }

        if (!ritz.converged) {
            // the small l x l decomposition itself failed; the result cannot be trusted
            converged = false;
        }

        // 4) Rayleigh-Ritz: rotate the subspace into the principal directions. y^T y is symmetric
        //    positive semi-definite, so its singular vectors are its eigenvectors and its singular
        //    values are the squared singular values of the data restricted to the subspace.
        components = new double[no_dims][n];
        singularValues = new double[no_dims];
        for (int k = 0; k < no_dims; k++) {
            singularValues[k] = Math.sqrt(Math.max(0.0, ritz.sigma[k]));
            double[] ck = components[k];
            for (int p = 0; p < l; p++) {
                double w = ritz.V[k * l + p];                  // entry p of eigenvector k
                if (w == 0.0) continue;
                double[] vp = v[p];
                for (int j = 0; j < n; j++) ck[j] += w * vp[j];
            }
        }

        double[][] projected = new double[m][no_dims];
        final double[][] fy = y;
        final FlatParallelJacobiSVD.Result fritz = ritz;
        final int fl = l, fno = no_dims;
        forEachRow(m, parallel, i -> {
            double[] yi = fy[i];
            double[] out = projected[i];
            for (int k = 0; k < fno; k++) {
                double s = 0.0;
                for (int p = 0; p < fl; p++) s += yi[p] * fritz.V[k * fl + p];
                out[k] = s;
            }
        });

        // 5) canonicalize signs exactly as JacobiPCA does, and carry the flip into the projection
        //    so that components and coordinates stay consistent
        for (int k = 0; k < no_dims; k++) {
            if (needsFlip(components[k])) {
                double[] ck = components[k];
                for (int j = 0; j < n; j++) ck[j] = -ck[j];
                for (int i = 0; i < m; i++) projected[i][k] = -projected[i][k];
            }
        }
        return projected;
    }

    /**
     * Whether the last {@link #pca} call reached its tolerance before {@code maxIterations}.
     *
     * @return {@code true} if the last fit converged
     */
    public boolean converged() { return converged; }

    /**
     * Iterations the last {@link #pca} call needed.
     *
     * @return the number of subspace iterations performed by the last fit
     */
    public int getIterations() { return iterations; }

    /**
     * Per-feature mean used for centering (length n).
     *
     * @return the per-feature mean of the last fit
     */
    public double[] getMean() { return mean; }

    /**
     * Kept principal directions, {@code [no_dims][n]}, row k = k-th component.
     *
     * @return the principal directions of the last fit
     */
    public double[][] getComponents() { return components; }

    /**
     * Singular values of the centered matrix for the kept components (descending).
     *
     * @return the singular values of the last fit
     */
    public double[] getSingularValues() { return singularValues; }

    // ------------------------------------------------------------------

    private double[][] startingBlock(int l, int n) {
        Random random = new Random(seed);
        double[][] v = new double[l][n];
        for (int k = 0; k < l; k++)
            for (int j = 0; j < n; j++) v[k][j] = random.nextGaussian();
        return v;
    }

    /** (A - 1 mean^T) v^T, m x l */
    private double[][] project(double[][] a, double[][] v, boolean parallel) {
        final int l = v.length;
        final double[][] out = new double[m][l];
        final double[] mu = mean;
        final int cols = n;
        forEachRow(m, parallel, i -> {
            double[] ai = a[i];
            double[] oi = out[i];
            for (int k = 0; k < l; k++) {
                double[] vk = v[k];
                double s = 0.0;
                for (int j = 0; j < cols; j++) s += (ai[j] - mu[j]) * vk[j];
                oi[k] = s;
            }
        });
        return out;
    }

    /** (A - 1 mean^T)^T y, l x n */
    private double[][] backProject(double[][] a, double[][] y, boolean parallel) {
        final int l = y[0].length;
        final double[][] out = new double[l][n];
        final double[] mu = mean;
        final int rows = m, cols = n;
        forEachRow(l, parallel, k -> {
            double[] ok = out[k];
            for (int i = 0; i < rows; i++) {
                double yik = y[i][k];
                if (yik == 0.0) continue;
                double[] ai = a[i];
                for (int j = 0; j < cols; j++) ok[j] += yik * (ai[j] - mu[j]);
            }
        });
        return out;
    }

    /** Eigen decomposition of the l x l Gram matrix y^T y, via the existing Jacobi SVD. */
    private static FlatParallelJacobiSVD.Result decomposeGram(double[][] y, int l) {
        final int rows = y.length;
        double[] gram = new double[l * l];
        for (int p = 0; p < l; p++) {
            for (int q = p; q < l; q++) {
                double s = 0.0;
                for (int i = 0; i < rows; i++) s += y[i][p] * y[i][q];
                gram[q * l + p] = s;
                gram[p * l + q] = s;
            }
        }
        return new FlatParallelJacobiSVD().decomposeInPlace(gram, l, l);
    }

    /**
     * Whether the leading {@code keep} directions are as resolved as they are going to get.
     * <p>
     * Measured against the largest Ritz value, so a direction that carries very little variance is
     * not asked for more relative accuracy than the leading one. Such a direction may well fail
     * this test for many iterations, because where the spectrum is nearly tied the data does not
     * determine the direction at all - deciding whether that matters is left to the caller, which
     * knows what the components are for.
     */
    private boolean resolved(double[] previous, double[] current, int keep) {
        double scale = current[0];
        if (!(scale > 0.0)) return true;               // all-zero data: nothing left to converge to
        for (int k = 0; k < keep; k++) {
            if (Math.abs(current[k] - previous[k]) > tolerance * scale) return false;
        }
        return true;
    }

    /** modified Gram-Schmidt over the rows, run twice for stability */
    private static void orthonormalizeRows(double[][] v) {
        final int cols = v[0].length;
        for (int pass = 0; pass < 2; pass++) {
            for (int k = 0; k < v.length; k++) {
                double[] vk = v[k];
                for (int p = 0; p < k; p++) {
                    double[] vp = v[p];
                    double dot = 0.0;
                    for (int j = 0; j < cols; j++) dot += vk[j] * vp[j];
                    for (int j = 0; j < cols; j++) vk[j] -= dot * vp[j];
                }
                double norm = 0.0;
                for (int j = 0; j < cols; j++) norm += vk[j] * vk[j];
                norm = Math.sqrt(norm);
                if (norm > 0.0 && norm * 0.0 == 0.0) {
                    for (int j = 0; j < cols; j++) vk[j] /= norm;
                } else {
                    // linearly dependent direction: the subspace has less rank than l. Zeroing the
                    // row makes it contribute a zero Ritz value, which sorts last and is therefore
                    // never among the kept components.
                    java.util.Arrays.fill(vk, 0.0);
                }
            }
        }
    }

    /** {@code true} if the entry of largest magnitude is negative (the svd_flip convention). */
    private static boolean needsFlip(double[] component) {
        int argmax = 0;
        double best = -1.0;
        for (int j = 0; j < component.length; j++) {
            double abs = Math.abs(component[j]);
            if (abs > best) { best = abs; argmax = j; }
        }
        return component[argmax] < 0.0;
    }

    private static void forEachRow(int count, boolean parallel, java.util.function.IntConsumer body) {
        IntStream range = IntStream.range(0, count);
        if (parallel) range = range.parallel();
        range.forEach(body);
    }

    private static void validate(double[][] matrix, int no_dims) {
        if (matrix == null || matrix.length == 0 || matrix[0].length == 0)
            throw new IllegalArgumentException("empty matrix");
        int n = matrix[0].length;
        for (double[] row : matrix)
            if (row.length != n)
                throw new IllegalArgumentException("matrix must be rectangular");
        int maxComp = Math.min(matrix.length, n);
        if (no_dims < 1 || no_dims > maxComp)
            throw new IllegalArgumentException(
                "no_dims must be in [1, min(m,n)] = [1, " + maxComp + "], got " + no_dims);
    }
}
