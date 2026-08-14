package math.linalg;

/**
 * Principal Component Analysis built on {@link FlatParallelJacobiSVD}, without any
 * external dependency.
 *
 * <p>Input convention: {@code double[][] matrix} is <b>row-major</b>, i.e. each
 * {@code matrix[i]} is one sample and each column is one feature. With m = number of
 * samples and n = number of features:</p>
 * <ul>
 *   <li>the data is mean-centered per feature (column),</li>
 *   <li>the top {@code no_dims} principal directions are the leading right singular
 *       vectors of the centered matrix,</li>
 *   <li>each centered sample is projected onto those directions.</li>
 * </ul>
 *
 * <p>Both shapes are handled:</p>
 * <ul>
 *   <li><b>m &gt;= n</b> (tall/typical): SVD of the centered matrix A directly; the
 *       principal directions are the columns of V.</li>
 *   <li><b>m &lt; n</b> (wide, more features than samples): the underlying Jacobi SVD
 *       needs rows &gt;= cols, so we decompose A^T instead. The principal directions of
 *       A are then the <i>left</i> singular vectors of A^T (the columns of U). There are
 *       at most min(m, n) of them.</li>
 * </ul>
 *
 * <p>Singular vectors are only defined up to sign. To make the output deterministic and
 * independent of the solver path, each component's sign is canonicalized so that its
 * entry of largest magnitude is positive (the "svd_flip" convention). This differs from
 * an arbitrary blanket sign flip but is irrelevant for downstream use such as t-SNE.</p>
 * 
 * @since 1.5.1
 */
public final class JacobiPCA {

    private final FlatParallelJacobiSVD svd;

    // fitted state
    private double[]   mean;            // per-feature mean, length n
    private double[][] components;      // [no_dims][n], row k = k-th principal direction
    private double[]   singularValues;  // [no_dims], descending
    private int m, n;                   // samples, features of the last fit

    public JacobiPCA() {
        this.svd = new FlatParallelJacobiSVD();
    }

    public JacobiPCA(FlatParallelJacobiSVD svd) {
        this.svd = svd;
    }

    /**
     * Fits the PCA on {@code matrix} (row-major, row = sample) and returns the samples
     * projected onto the first {@code no_dims} principal components.
     *
     * @param matrix  row-major data, m samples x n features (rectangular, non-empty)
     * @param no_dims number of components to keep, 1 &lt;= no_dims &lt;= min(m, n)
     * @return {@code double[m][no_dims]} of projected coordinates
     */
    public double[][] pca(double[][] matrix, int no_dims) {
        validate(matrix, no_dims);

        m = matrix.length;
        n = matrix[0].length;

        // 1) per-feature mean
        mean = new double[n];
        for (double[] row : matrix)
            for (int j = 0; j < n; j++) mean[j] += row[j];
        for (int j = 0; j < n; j++) mean[j] /= m;

        // 2) centered copy (reused for packing and projection)
        double[][] centered = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                centered[i][j] = matrix[i][j] - mean[j];

        // 3) SVD -> principal directions, picking the right path for the shape
        components     = new double[no_dims][n];
        singularValues = new double[no_dims];

        if (m > n) {
            // Tall case: QR-precondition. A = Q R with Q having orthonormal columns, so the
            // right singular vectors and singular values of A equal those of the small n x n R
            // (A = Q R = (Q U_R) Sigma V^T). We only need V and sigma for PCA -> no need for Q.
            // Cost drops from O(sweeps*m*n^2) to O(m*n^2) [QR] + O(sweeps*n^3) [Jacobi on R],
            // and the Jacobi part no longer touches an m x n array.
            double[] flat = new double[m * n];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    flat[j * m + i] = centered[i][j];   // column-major: flat[j*m + i]

            double[] R = qrToR(flat, m, n);             // consumes flat, returns n x n (col-major)
            FlatParallelJacobiSVD.Result r = svd.decomposeInPlace(R, n, n);
            requireConverged(r);
            for (int k = 0; k < no_dims; k++) {
                singularValues[k] = r.sigma[k];
                for (int j = 0; j < n; j++)
                    components[k][j] = r.V[k * n + j];   // column k of V (right singular vector)
            }
        } else if (m == n) {
            // Square case: QR would add O(n^3) without shrinking the problem -> SVD directly.
            double[] flat = new double[m * n];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    flat[j * m + i] = centered[i][j];

            FlatParallelJacobiSVD.Result r = svd.decomposeInPlace(flat, m, n);
            requireConverged(r);
            for (int k = 0; k < no_dims; k++) {
                singularValues[k] = r.sigma[k];
                for (int j = 0; j < n; j++)
                    components[k][j] = r.V[k * n + j];   // column k of V
            }
        } else {
            // Wide case (m < n): decompose B = A^T (n x m, rows >= cols since n > m).
            // A = B^T = V_B Sigma U_B^T  =>  right vectors of A = columns of U_B (length n).
            // QR would not help here: we need B's LEFT vectors, which need Q of B's QR.
            double[] bflat = new double[n * m];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    bflat[i * n + j] = centered[i][j];   // B[j][i]=A[i][j], col-major over n rows

            FlatParallelJacobiSVD.Result r = svd.decomposeInPlace(bflat, n, m);
            requireConverged(r);
            for (int k = 0; k < no_dims; k++) {
                singularValues[k] = r.sigma[k];
                for (int j = 0; j < n; j++)
                    components[k][j] = r.U[k * n + j];   // column k of U_B (length n)
            }
        }

        // 4) canonicalize signs (largest-magnitude entry of each component positive)
        canonicalizeSigns(components);

        // 5) project the centered samples onto the kept components
        double[][] projected = new double[m][no_dims];
        for (int i = 0; i < m; i++) {
            double[] ci = centered[i];
            for (int k = 0; k < no_dims; k++) {
                double[] ck = components[k];
                double dot = 0.0;
                for (int j = 0; j < n; j++) dot += ci[j] * ck[j];
                projected[i][k] = dot;
            }
        }
        return projected;
    }

    /**
     * Projects new, unseen samples (row-major, n features each) onto the fitted
     * components. {@link #pca} must have been called first.
     *
     * @param samples row-major samples to project, each row with {@code n} features
     * @return projected coordinates with shape {@code [samples.length][no_dims]}
     */
    public double[][] transform(double[][] samples) {
        if (components == null) throw new IllegalStateException("call pca(...) first");
        int rows = samples.length;
        int no_dims = components.length;
        double[][] out = new double[rows][no_dims];
        for (int i = 0; i < rows; i++) {
            if (samples[i].length != n)
                throw new IllegalArgumentException("unexpected sample length");
            for (int k = 0; k < no_dims; k++) {
                double[] ck = components[k];
                double dot = 0.0;
                for (int j = 0; j < n; j++) dot += (samples[i][j] - mean[j]) * ck[j];
                out[i][k] = dot;
            }
        }
        return out;
    }

    /**
     * Per-feature mean used for centering (length n).
     *
     * @return mean vector of the last fitted dataset
     */
    public double[] getMean() { return mean; }

    /**
     * Kept principal directions, {@code [no_dims][n]}, row k = k-th component.
     *
     * @return principal component matrix of the last fit
     */
    public double[][] getComponents() { return components; }

    /**
     * Singular values of the centered matrix for the kept components (descending).
     *
     * @return singular values for the kept components
     */
    public double[] getSingularValues() { return singularValues; }

    /**
     * Variance explained by each kept component: sigma_k^2 / (m - 1).
     * (Eigenvalues of the sample covariance matrix.)
     *
     * @return explained variance per kept component
     */
    public double[] getExplainedVariance() {
        double denom = Math.max(1, m - 1);
        double[] ev = new double[singularValues.length];
        for (int k = 0; k < ev.length; k++)
            ev[k] = singularValues[k] * singularValues[k] / denom;
        return ev;
    }

    // ------------------------------------------------------------------

    /**
     * Householder QR of a column-major m x n matrix (m >= n), in place; returns the
     * n x n upper-triangular factor R (column-major). Q is not formed. The input array
     * {@code a} is destroyed. For PCA only R is needed, since A = Q R shares A's right
     * singular vectors and singular values.
     */
    private static double[] qrToR(double[] a, int m, int n) {
        for (int k = 0; k < n; k++) {
            int kk = k * m + k;
            double x0 = a[kk];

            // overflow-safe norm of the sub-column rows k..m-1 (BLAS dnrm2 style)
            double norm = safeNorm(a, k * m + k, m - k);
            if (norm == 0.0) continue;                 // zero column -> R[k][k] = 0

            double alpha = (x0 >= 0.0) ? -norm : norm; // reflect away from x0 for stability
            double v0    = x0 - alpha;                 // Householder vector: v = (v0, x_{k+1..})
            // vnorm2 = v0^2 + sum_{i>k} x_i^2 = 2*norm*(norm + |x0|), avoiding norm^2 overflow
            double vnorm2 = 2.0 * norm * (norm + Math.abs(x0));
            if (vnorm2 == 0.0) { a[kk] = alpha; continue; }
            double beta = 2.0 / vnorm2;

            // apply H = I - beta v v^T to the trailing columns j = k+1..n-1
            for (int j = k + 1; j < n; j++) {
                int jk = j * m + k;
                double s = v0 * a[jk];
                for (int i = k + 1; i < m; i++) s += a[k * m + i] * a[j * m + i];
                s *= beta;
                a[jk] -= s * v0;
                for (int i = k + 1; i < m; i++) a[j * m + i] -= s * a[k * m + i];
            }
            a[kk] = alpha;                             // R diagonal
        }

        // extract upper-triangular R (n x n), column-major
        double[] R = new double[n * n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i <= j; i++)
                R[j * n + i] = a[j * m + i];
        return R;
    }

    /** ||a[off..off+len)||_2 without over-/underflow (BLAS dnrm2 style). */
    private static double safeNorm(double[] a, int off, int len) {
        double scale = 0.0, ssq = 1.0;
        for (int i = 0; i < len; i++) {
            double x = a[off + i];
            if (x != 0.0) {
                double ax = Math.abs(x);
                if (scale < ax) { double r = scale / ax; ssq = 1.0 + ssq * r * r; scale = ax; }
                else            { double r = ax / scale; ssq += r * r; }
            }
        }
        return scale * Math.sqrt(ssq);
    }

    private static void requireConverged(FlatParallelJacobiSVD.Result r) {
        if (!r.converged)
            throw new IllegalStateException(
                "Jacobi SVD did not converge; PCA result would be unreliable");
    }

    private static void canonicalizeSigns(double[][] components) {
        for (double[] comp : components) {
            int argmax = 0;
            double best = -1.0;
            for (int j = 0; j < comp.length; j++) {
                double abs = Math.abs(comp[j]);
                if (abs > best) { best = abs; argmax = j; }
            }
            if (comp[argmax] < 0.0)
                for (int j = 0; j < comp.length; j++) comp[j] = -comp[j];
        }
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

    // ------------------------------------------------------------------
    // Demo / self-check
    // ------------------------------------------------------------------

    public static void main(String[] args) {
        // Tall case: 500 samples x 8 features, correlated structure.
        runDemo("Tall (m >= n): 500 x 8", makeData(500, 8, 111L), 3);
        // Wide case: 6 samples x 40 features -> exercises the A^T path.
        runDemo("Wide (m < n):  6 x 40",   makeData(6, 40, 222L),   3);
    }

    private static void runDemo(String label, double[][] data, int noDims) {
        JacobiPCA pca = new JacobiPCA();
        double[][] y = pca.pca(data, noDims);

        System.out.println("=== " + label + " ===");
        System.out.printf("  projected shape          = %d x %d%n", y.length, y[0].length);
        double[] sv = pca.getSingularValues();
        System.out.print  ("  singular values          = [");
        for (int k = 0; k < sv.length; k++)
            System.out.printf("%s%.5g", k == 0 ? "" : ", ", sv[k]);
        System.out.println("]");
        System.out.printf("  components orthonormal    = %.3e (max |C C^T - I|)%n",
                componentOrthoError(pca.getComponents()));
        // components must be sorted by singular value (variance) descending
        boolean desc = true;
        for (int k = 0; k + 1 < sv.length; k++) if (sv[k] < sv[k + 1]) desc = false;
        System.out.println("  singular values sorted    = " + desc);
        System.out.println();
    }

    /** max |C C^T - I| for components C ([k][n]); each row should be unit and mutually orthogonal. */
    private static double componentOrthoError(double[][] c) {
        double max = 0;
        for (int p = 0; p < c.length; p++)
            for (int q = 0; q < c.length; q++) {
                double dot = 0;
                for (int j = 0; j < c[p].length; j++) dot += c[p][j] * c[q][j];
                max = Math.max(max, Math.abs(dot - (p == q ? 1.0 : 0.0)));
            }
        return max;
    }

    /** Deterministic synthetic data with a few dominant directions. */
    private static double[][] makeData(int m, int n, long seed) {
        double[][] d = new double[m][n];
        long lcg = seed;
        for (int i = 0; i < m; i++) {
            // two latent factors drive most of the variance
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f1 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f2 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            for (int j = 0; j < n; j++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double noise = (((lcg >>> 11) / (double) (1L << 53)) - 0.5) * 0.05;
                d[i][j] = 3.0 * f1 * Math.cos(j) + 2.0 * f2 * Math.sin(j) + noise;
            }
        }
        return d;
    }
}
