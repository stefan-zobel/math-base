package math.linalg;

/**
 * Principal Component Analysis on a <b>covariance matrix</b> (or a correlation,
 * kernel or graph Laplacian matrix), built on {@link SymmetricJacobiEigen}.
 * The principal directions are the eigenvectors of the input, the explained
 * variances are its eigenvalues, both in descending order.
 *
 * <p><b>Why there is a second PCA implementation.</b> {@link JacobiPCA} is the
 * one to use whenever the raw data matrix is at hand: it decomposes the
 * centered data directly and is the numerically better route (see below). What
 * it cannot do is start from a second-moment matrix, and that is a real gap:</p>
 * <ul>
 *   <li>the data was streamed and discarded, and only a running covariance
 *       survives (see {@code math.probe.CovarianceAccumulator});</li>
 *   <li>the data set does not fit in memory, but an {@code n x n} matrix does;</li>
 *   <li>the covariance is the model, not a summary of samples - factor models,
 *       risk matrices, a shrinkage estimator, a covariance published without the
 *       underlying observations;</li>
 *   <li>the matrix is not a covariance at all but is still to be analyzed
 *       spectrally: a correlation matrix, a kernel/Gram matrix (kernel PCA), a
 *       graph Laplacian (spectral clustering, Fiedler vector).</li>
 * </ul>
 *
 * <p><b>What it costs.</b> Forming {@code C = A_c^T A_c / (m-1)} squares the
 * condition number of the data. Small eigenvalues therefore drown at
 * {@code sqrt(eps) ~ 1.5e-8} relative to the largest one, where the SVD route of
 * {@code JacobiPCA} resolves them down to {@code eps ~ 2.2e-16}. That is eight
 * digits, and it is given away precisely by the property one-sided Jacobi was
 * chosen for, namely high relative accuracy on small singular values. For the
 * leading components of a spectrum with well separated variances the loss is
 * invisible; for a near degenerate spectrum, an ill conditioned data set, or
 * components far down the tail it is not. This class is therefore about
 * <em>reach</em>, not about speed - it is not measurably faster than
 * {@code JacobiPCA} on data that both can handle.</p>
 *
 * <p><b>What the eigensolver adds over the SVD here.</b> A covariance matrix is
 * positive semidefinite in exact arithmetic only. Estimators that shrink,
 * accumulate in a stream, or handle missing values pairwise routinely produce a
 * matrix with slightly negative eigenvalues. An SVD would report {@code |lambda|}
 * and hide the defect; {@link #getEigenvalues()} and {@link #getMinEigenvalue()}
 * show it. Whitening (see {@link #whiten(double[][])}) is exactly the operation
 * that would silently produce nonsense on such an input, and it is refused here
 * instead.</p>
 *
 * <p>Input convention: {@code cov} is flat and column-major
 * ({@code C[j*n + i]}), the layout used throughout {@code math.linalg}. For a
 * symmetric matrix column-major and row-major coincide, so the caller does not
 * have to care. Component signs follow the same convention as
 * {@link JacobiPCA}: each direction is flipped so that its entry of largest
 * magnitude is positive, which makes the output of the two implementations
 * directly comparable.</p>
 *
 * @since 1.5.1
 */
public final class CovariancePCA {

    /** Relative tolerance on the asymmetry of the input matrix. */
    private static final double SYMMETRY_TOL = 1.0e-10;

    private final SymmetricJacobiEigen eig;

    // fitted state
    private double[][] components;        // [noDims][n], row k = k-th principal direction
    private double[]   explainedVariance; // [noDims], the leading eigenvalues
    private double[]   eigenvalues;       // [n], all of them, descending
    private double[]   mean;              // [n] or null if none was supplied
    private double     totalVariance;     // trace of the input
    private int        n;

    /** Creates a PCA using a default-configured {@link SymmetricJacobiEigen}. */
    public CovariancePCA() {
        this(new SymmetricJacobiEigen());
    }

    /**
     * Creates a PCA using a specific eigensolver configuration.
     *
     * @param eig the symmetric eigensolver to decompose the input with
     */
    public CovariancePCA(SymmetricJacobiEigen eig) {
        this.eig = eig;
    }

    /**
     * Fits the PCA on a covariance matrix without a mean vector.
     * {@link #transform(double[][])} then assumes its input is already centered.
     *
     * @param cov     symmetric {@code n x n} matrix, flat column-major, length
     *                {@code n*n}; left unchanged
     * @param n       the dimension, i.e. the number of features
     * @param noDims  number of components to keep, {@code 1 <= noDims <= n}
     * @return this instance, fitted
     * @throws IllegalArgumentException if the arguments are inconsistent or the
     *                                  matrix is not symmetric
     * @throws IllegalStateException    if the eigensolver did not converge
     */
    public CovariancePCA fit(double[] cov, int n, int noDims) {
        return fit(cov, n, null, noDims);
    }

    /**
     * Fits the PCA on a covariance matrix together with the mean vector of the
     * data it came from, which {@link #transform(double[][])} needs in order to
     * center new samples.
     *
     * @param cov     symmetric {@code n x n} matrix, flat column-major, length
     *                {@code n*n}; left unchanged
     * @param n       the dimension, i.e. the number of features
     * @param mean    per-feature mean of length {@code n}, or {@code null} if
     *                samples handed to {@code transform} are already centered;
     *                copied defensively
     * @param noDims  number of components to keep, {@code 1 <= noDims <= n}
     * @return this instance, fitted
     * @throws IllegalArgumentException if the arguments are inconsistent or the
     *                                  matrix is not symmetric
     * @throws IllegalStateException    if the eigensolver did not converge
     */
    public CovariancePCA fit(double[] cov, int n, double[] mean, int noDims) {
        validate(cov, n, mean, noDims);

        SymmetricJacobiEigen.Result r = eig.decompose(cov, n);
        if (!r.converged) {
            throw new IllegalStateException(
                    "Jacobi eigensolver did not converge; PCA result would be unreliable");
        }

        this.n = n;
        this.mean = (mean == null) ? null : mean.clone();
        this.eigenvalues = r.lambda;                 // freshly allocated by the solver
        this.explainedVariance = new double[noDims];
        this.components = new double[noDims][n];
        for (int k = 0; k < noDims; k++) {
            explainedVariance[k] = r.lambda[k];
            System.arraycopy(r.V, k * n, components[k], 0, n);   // column k of V
        }
        canonicalizeSigns(components);

        double trace = 0.0;
        for (int j = 0; j < n; j++) {
            trace += cov[j * n + j];
        }
        this.totalVariance = trace;
        return this;
    }

    /**
     * Projects samples onto the kept components. Samples are centered with the
     * mean vector passed to {@code fit}; if none was passed, they are used as
     * they are.
     *
     * @param samples row-major samples, {@code n} features each
     * @return projected coordinates, {@code [samples.length][noDims]}
     * @throws IllegalStateException    if {@code fit} has not been called
     * @throws IllegalArgumentException if a sample has the wrong length
     */
    public double[][] transform(double[][] samples) {
        requireFitted();
        final int noDims = components.length;
        double[][] out = new double[samples.length][noDims];
        for (int i = 0; i < samples.length; i++) {
            double[] s = samples[i];
            if (s.length != n) {
                throw new IllegalArgumentException(
                        "expected " + n + " features, got " + s.length);
            }
            for (int k = 0; k < noDims; k++) {
                double[] ck = components[k];
                double dot = 0.0;
                for (int j = 0; j < n; j++) {
                    dot += (mean == null ? s[j] : s[j] - mean[j]) * ck[j];
                }
                out[i][k] = dot;
            }
        }
        return out;
    }

    /**
     * Projects samples like {@link #transform(double[][])} and additionally
     * divides each coordinate by {@code sqrt(lambda_k)}, so that the transformed
     * data has unit variance along every kept component (PCA whitening, the
     * building block of the Mahalanobis distance).
     *
     * @param samples row-major samples, {@code n} features each
     * @return whitened coordinates, {@code [samples.length][noDims]}
     * @throws IllegalStateException    if {@code fit} has not been called, or if
     *                                  a kept eigenvalue is not strictly
     *                                  positive - the input is then rank
     *                                  deficient or not positive definite and
     *                                  whitening is undefined
     * @throws IllegalArgumentException if a sample has the wrong length
     */
    public double[][] whiten(double[][] samples) {
        requireFitted();
        final int noDims = components.length;
        double[] inv = new double[noDims];
        for (int k = 0; k < noDims; k++) {
            double lambda = explainedVariance[k];
            if (!(lambda > 0.0)) {
                throw new IllegalStateException(
                        "cannot whiten: eigenvalue " + k + " is " + lambda
                                + " (input is rank deficient or not positive definite)");
            }
            inv[k] = 1.0 / Math.sqrt(lambda);
        }
        double[][] out = transform(samples);
        for (int i = 0; i < out.length; i++) {
            for (int k = 0; k < noDims; k++) {
                out[i][k] *= inv[k];
            }
        }
        return out;
    }

    /**
     * Kept principal directions, {@code [noDims][n]}, row {@code k} being the
     * {@code k}-th direction.
     *
     * @return the component matrix of the last fit
     */
    public double[][] getComponents() {
        return components;
    }

    /**
     * Variance explained by each kept component, i.e. the leading eigenvalues in
     * descending order.
     *
     * @return explained variance per kept component
     */
    public double[] getExplainedVariance() {
        return explainedVariance;
    }

    /**
     * Fraction of the total variance explained by each kept component,
     * {@code lambda_k / trace}. Zero throughout if the input has zero trace.
     *
     * @return explained variance ratio per kept component
     */
    public double[] getExplainedVarianceRatio() {
        requireFitted();
        double[] ratio = new double[explainedVariance.length];
        if (totalVariance != 0.0) {
            for (int k = 0; k < ratio.length; k++) {
                ratio[k] = explainedVariance[k] / totalVariance;
            }
        }
        return ratio;
    }

    /**
     * All {@code n} eigenvalues of the input in descending order, not just the
     * kept ones. Useful for choosing a cutoff and for inspecting the tail of the
     * spectrum.
     *
     * @return the full eigenvalue spectrum
     */
    public double[] getEigenvalues() {
        return eigenvalues;
    }

    /**
     * The smallest eigenvalue of the input. A value that is negative by more
     * than rounding means the matrix is not positive semidefinite and hence not
     * a valid covariance matrix.
     *
     * @return the smallest eigenvalue
     * @throws IllegalStateException if {@code fit} has not been called
     */
    public double getMinEigenvalue() {
        requireFitted();
        return eigenvalues[n - 1];
    }

    /**
     * Total variance of the input, i.e. its trace, which equals the sum of all
     * eigenvalues.
     *
     * @return the trace of the fitted matrix
     */
    public double getTotalVariance() {
        return totalVariance;
    }

    /**
     * The mean vector passed to {@code fit}, or {@code null} if none was given.
     *
     * @return a fresh copy of the mean vector, or {@code null}
     */
    public double[] getMean() {
        return (mean == null) ? null : mean.clone();
    }

    // ------------------------------------------------------------------

    private void requireFitted() {
        if (components == null) {
            throw new IllegalStateException("call fit(...) first");
        }
    }

    /**
     * Flips each direction so that its entry of largest magnitude is positive.
     * Same convention as {@code JacobiPCA}, so that the two implementations
     * produce comparable output for the same data.
     */
    private static void canonicalizeSigns(double[][] components) {
        for (int k = 0; k < components.length; k++) {
            double[] comp = components[k];
            int argmax = 0;
            double best = -1.0;
            for (int j = 0; j < comp.length; j++) {
                double abs = Math.abs(comp[j]);
                if (abs > best) { best = abs; argmax = j; }
            }
            if (comp[argmax] < 0.0) {
                for (int j = 0; j < comp.length; j++) comp[j] = -comp[j];
            }
        }
    }

    private static void validate(double[] cov, int n, double[] mean, int noDims) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1, got " + n);
        }
        if (cov == null || cov.length != n * n) {
            throw new IllegalArgumentException(
                    "cov must have length n*n = " + (n * n)
                            + ", got " + (cov == null ? "null" : Integer.toString(cov.length)));
        }
        if (noDims < 1 || noDims > n) {
            throw new IllegalArgumentException(
                    "noDims must be in [1, " + n + "], got " + noDims);
        }
        if (mean != null && mean.length != n) {
            throw new IllegalArgumentException(
                    "mean must have length " + n + ", got " + mean.length);
        }
        double maxAbs = 0.0;
        for (int idx = 0; idx < cov.length; idx++) {
            maxAbs = Math.max(maxAbs, Math.abs(cov[idx]));
        }
        double asym = SymmetricJacobiEigen.symmetryError(cov, n);
        if (asym > SYMMETRY_TOL * Math.max(1.0, maxAbs)) {
            throw new IllegalArgumentException(
                    "input is not symmetric: max|C - C^T| = " + asym);
        }
    }

    // ------------------------------------------------------------------
    // Demo / self-check
    // ------------------------------------------------------------------

    /**
     * Demo and self-check; prints residuals and an overall verdict.
     *
     * @param args ignored
     */
    public static void main(String[] args) {
        boolean ok = true;
        ok &= streamedVersusRawData();
        ok &= correlationInput();
        ok &= rankDeficientInput();
        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }

    /**
     * The covariance route must reproduce what {@code JacobiPCA} computes from
     * the raw data, even though the covariance was accumulated in a stream and
     * the data itself was never kept.
     */
    private static boolean streamedVersusRawData() {
        final int m = 500, n = 8, k = 3;
        double[][] data = makeData(m, n, 111L);

        // stream the rows through the accumulator, keeping nothing
        math.probe.CovarianceAccumulator acc = new math.probe.CovarianceAccumulator(n);
        acc.addAll(data);
        double[] cov = acc.covariance();

        CovariancePCA cpca = new CovariancePCA().fit(cov, n, acc.mean(), k);
        JacobiPCA rpca = new JacobiPCA();
        rpca.pca(data, k);

        double compErr = maxDiff(cpca.getComponents(), rpca.getComponents());
        double varErr = maxRelDiff(cpca.getExplainedVariance(), rpca.getExplainedVariance());
        double meanErr = maxRelDiff(acc.mean(), rpca.getMean());
        double orthoErr = orthonormalityError(cpca.getComponents());

        // whitened data must have unit variance along every kept component
        double[][] w = cpca.whiten(data);
        double whitenErr = 0.0;
        for (int c = 0; c < k; c++) {
            double s2 = 0.0;
            for (int i = 0; i < m; i++) s2 += w[i][c] * w[i][c];
            whitenErr = Math.max(whitenErr, Math.abs(s2 / (m - 1) - 1.0));
        }

        double ratioSum = 0.0;
        for (double r : cpca.getExplainedVarianceRatio()) ratioSum += r;

        System.out.println("=== Streamed covariance vs. raw data (500 x 8, 3 components) ===");
        System.out.printf("  max|component_cov - component_raw|  = %.3e%n", compErr);
        System.out.printf("  max rel. diff of explained variance = %.3e%n", varErr);
        System.out.printf("  max rel. diff of the mean vector    = %.3e%n", meanErr);
        System.out.printf("  max |C C^T - I|                     = %.3e%n", orthoErr);
        System.out.printf("  max |var(whitened) - 1|             = %.3e%n", whitenErr);
        System.out.printf("  explained variance ratio (sum)      = %.5f%n", ratioSum);
        System.out.println();

        return compErr < 1e-6 && varErr < 1e-9 && meanErr < 1e-12
                && orthoErr < 1e-12 && whitenErr < 1e-12
                && ratioSum > 0.0 && ratioSum <= 1.0 + 1e-12;
    }

    /** A correlation matrix has unit diagonal, so its eigenvalues must sum to n. */
    private static boolean correlationInput() {
        final int m = 300, n = 6, k = 2;
        double[][] data = makeData(m, n, 222L);
        math.probe.CovarianceAccumulator acc = new math.probe.CovarianceAccumulator(n);
        acc.addAll(data);

        CovariancePCA pca = new CovariancePCA().fit(acc.correlation(), n, k);
        double sum = 0.0;
        for (double lambda : pca.getEigenvalues()) sum += lambda;

        double orthoErr = orthonormalityError(pca.getComponents());
        System.out.println("=== Correlation matrix input (300 x 6, 2 components) ===");
        System.out.printf("  |sum(lambda) - n|                   = %.3e%n", Math.abs(sum - n));
        System.out.printf("  |trace - n|                         = %.3e%n",
                Math.abs(pca.getTotalVariance() - n));
        System.out.printf("  max |C C^T - I|                     = %.3e%n", orthoErr);
        System.out.println();

        return Math.abs(sum - n) < 1e-12 * n && orthoErr < 1e-12;
    }

    /**
     * A feature that is an exact linear combination of two others makes the
     * covariance singular. The eigenvector basis stays complete and orthonormal
     * - an SVD would hand out a zero column instead - and the defect shows up in
     * the smallest eigenvalue. Whether rounding puts that eigenvalue just above
     * or just below zero is a coin flip, which is why the whitening guard cannot
     * be the only line of defence: the amplification factor
     * {@code 1/sqrt(lambda_min)} is the number to look at.
     */
    private static boolean rankDeficientInput() {
        final int m = 200, n = 4;
        double[][] base = makeData(m, n - 1, 333L);
        double[][] data = new double[m][n];
        for (int i = 0; i < m; i++) {
            System.arraycopy(base[i], 0, data[i], 0, n - 1);
            data[i][n - 1] = base[i][0] + base[i][1];    // linearly dependent
        }
        math.probe.CovarianceAccumulator acc = new math.probe.CovarianceAccumulator(n);
        acc.addAll(data);

        CovariancePCA pca = new CovariancePCA().fit(acc.covariance(), n, n);
        double minLambda = pca.getMinEigenvalue();
        double maxLambda = pca.getEigenvalues()[0];
        double orthoErr = orthonormalityError(pca.getComponents());

        // an exactly singular input must be caught by the whitening guard
        boolean refused = false;
        try {
            new CovariancePCA().fit(new double[] { 4.0, 0.0, 0.0, 0.0 }, 2, 2)
                    .whiten(new double[][] { { 1.0, 1.0 } });
        } catch (IllegalStateException expected) {
            refused = true;
        }

        System.out.println("=== Rank deficient covariance (200 x 4, rank 3) ===");
        System.out.printf("  lambda_min / lambda_max             = %.3e%n", minLambda / maxLambda);
        System.out.printf("  1 / sqrt(|lambda_min|)              = %.3e%n",
                1.0 / Math.sqrt(Math.abs(minLambda)));
        System.out.printf("  max |C C^T - I|                     = %.3e%n", orthoErr);
        System.out.println("  whitening of an exact zero refused  = " + refused);
        System.out.println();

        return Math.abs(minLambda) < 1e-12 * maxLambda && orthoErr < 1e-12 && refused;
    }

    /** max |C C^T - I| for components C ([k][n]). */
    private static double orthonormalityError(double[][] c) {
        double max = 0.0;
        for (int p = 0; p < c.length; p++) {
            for (int q = 0; q < c.length; q++) {
                double dot = 0.0;
                for (int j = 0; j < c[p].length; j++) dot += c[p][j] * c[q][j];
                max = Math.max(max, Math.abs(dot - (p == q ? 1.0 : 0.0)));
            }
        }
        return max;
    }

    private static double maxDiff(double[][] a, double[][] b) {
        double max = 0.0;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                max = Math.max(max, Math.abs(a[i][j] - b[i][j]));
            }
        }
        return max;
    }

    private static double maxRelDiff(double[] a, double[] b) {
        double max = 0.0;
        for (int i = 0; i < a.length; i++) {
            double scale = Math.max(1e-300, Math.abs(a[i]));
            max = Math.max(max, Math.abs(a[i] - b[i]) / scale);
        }
        return max;
    }

    /** Deterministic synthetic data with three dominant, well separated directions. */
    private static double[][] makeData(int m, int n, long seed) {
        double[][] d = new double[m][n];
        long lcg = seed;
        for (int i = 0; i < m; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f1 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f2 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f3 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            for (int j = 0; j < n; j++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double noise = (((lcg >>> 11) / (double) (1L << 53)) - 0.5) * 0.02;
                d[i][j] = 6.0 * f1 * Math.cos(j) + 3.0 * f2 * Math.sin(j)
                        + 1.5 * f3 * Math.cos(2 * j) + 10.0 + noise;
            }
        }
        return d;
    }
}
