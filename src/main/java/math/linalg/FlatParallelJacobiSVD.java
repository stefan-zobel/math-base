package math.linalg;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Stream;

/**
 * One-sided Jacobi SVD (Hestenes method) with
 *  - flat, column-major storage (BLAS/LAPACK layout: A[j*m + i]),
 *  - parallel round-robin ordering (Brent-Luk): per round n/2 column-disjoint
 *    rotations that can run concurrently without conflict.
 *
 * Decomposes A ({@code m x n}, {@code m >= n}) into {@code A = U * diag(sigma) * V^T}.
 * 
 * @since 1.5.1
 */
public final class FlatParallelJacobiSVD {

    /** Result; U, V stored column-major. */
    public static final class Result {
        public final double[] U;      // m x n, U[j*m + i]
        public final double[] V;      // n x n, V[j*n + i]
        public final double[] sigma;  // n, descending
        public final int m, n;
        public final boolean converged; // false if maxSweeps was exhausted (result may be inaccurate)
        Result(double[] U, double[] V, double[] sigma, int m, int n, boolean converged) {
            this.U = U; this.V = V; this.sigma = sigma; this.m = m; this.n = n;
            this.converged = converged;
        }
    }

    private final double eps;               // relative orthogonality threshold
    private final int    maxSweeps;         // safety upper bound
    private final long   parallelThreshold; // go parallel once m*n >= this size

    public FlatParallelJacobiSVD() {
        this(1e-15, 60, 1L << 16);          // default: parallel from ~65k elements
    }

    public FlatParallelJacobiSVD(double eps, int maxSweeps, long parallelThreshold) {
        this.eps = eps;
        this.maxSweeps = maxSweeps;
        this.parallelThreshold = parallelThreshold;
    }

    /**
     * @param A column-major packed, length m*n, element (i,j) = A[j*m + i]. Left unchanged.
     * @param m number of rows
     * @param n number of columns
     * @return decomposition result containing {@code U}, {@code V}, singular values, and convergence flag
     */
    public Result decompose(double[] A, int m, int n) {
        if (m < n) throw new IllegalArgumentException("expected m >= n");
        if (A.length != (long) m * n) throw new IllegalArgumentException("length != m*n");
        return decomposeInPlace(A.clone(), m, n);
    }

    /**
     * In-place variant: {@code a} is used as the working array and is overwritten (it
     * becomes U*Sigma / then U). {@code Result.U} aliases {@code a}, so do not reuse the
     * passed array afterwards. Saves one m*n copy versus {@link #decompose}.
     *
     * @param a column-major packed, length m*n, element (i,j) = a[j*m + i]. Destroyed.
     * @param m number of rows
     * @param n number of columns
     * @return decomposition result containing {@code U}, {@code V}, singular values, and convergence flag
     */
    public Result decomposeInPlace(double[] a, int m, int n) {
        if (m < n) throw new IllegalArgumentException("expected m >= n");
        if (a.length != (long) m * n) throw new IllegalArgumentException("length != m*n");

        final double[] v = new double[n * n];
        for (int j = 0; j < n; j++) v[j * n + j] = 1.0;   // V = I

        // --- global scaling into a safe range (largest column norm ~ 1) ---
        // Prevents over-/underflow when squaring in rotate(). Scalar scaling is
        // exactly compatible with the SVD: SVD(c*A) = U * (c*Sigma) * V^T.
        double maxNorm = 0.0;
        for (int j = 0; j < n; j++)
            maxNorm = Math.max(maxNorm, safeColumnNorm(a, j * m, m));

        double scale = 1.0;
        if (maxNorm > 0.0) {
            // exact: only shift the exponent (power of two -> no rounding error)
            scale = Math.scalb(1.0, -Math.getExponent(maxNorm));
            if (scale != 1.0)
                for (int idx = 0; idx < a.length; idx++) a[idx] *= scale;
        }

        // Round-robin schedule (padded up to the next even n).
        final int nEven = (n % 2 == 0) ? n : n + 1;
        final int[][][] schedule = roundRobin(nEven);
        final boolean parallel = (long) m * n >= parallelThreshold;

        boolean didConverge = false;
        for (int sweep = 0; sweep < maxSweeps; sweep++) {
            final AtomicBoolean converged = new AtomicBoolean(true);

            for (int[][] round : schedule) {
                // All pairs of ONE round are column-disjoint -> safe in parallel.
                Stream<int[]> s = Arrays.stream(round);
                if (parallel) s = s.parallel();
                s.forEach(pair -> {
                    int p = pair[0], q = pair[1];
                    if (p >= n || q >= n) return;          // dummy column -> skip
                    if (rotate(a, m, v, n, p, q)) converged.set(false);
                });
                // End of the stream acts as a barrier between rounds.
            }
            if (converged.get()) { didConverge = true; break; }
        }
        // Convergence guard: never return a silently non-converged result unnoticed.
        if (!didConverge)
            System.err.println("[FlatParallelJacobiSVD] WARNING: not converged after "
                    + maxSweeps + " sweeps (m=" + m + ", n=" + n
                    + "); result may be inaccurate. Consider raising maxSweeps.");

        // Singular values = column norms; U = normalized columns (in-place in a).
        final double[] sigma = new double[n];
        for (int j = 0; j < n; j++) {
            int ja = j * m;
            double norm = safeColumnNorm(a, ja, m);   // overflow-safe instead of norm += x*x
            sigma[j] = norm;
            if (norm > 0) {
                double inv = 1.0 / norm;
                for (int i = 0; i < m; i++) a[ja + i] *= inv;
            }
        }

        // undo the global scaling: we decomposed scale*A -> sigma_A = sigma/scale
        if (scale != 1.0) {
            double inv = 1.0 / scale;
            for (int j = 0; j < n; j++) sigma[j] *= inv;
        }

        sortDescending(sigma, a, m, v, n);
        return new Result(a, v, sigma, m, n, didConverge);
    }

    /**
     * A single Jacobi rotation on column pair (p, q). Returns true if a rotation was
     * actually applied (the pair was not yet orthogonal enough).
     */
    private boolean rotate(double[] a, int m, double[] v, int n, int p, int q) {
        final int pa = p * m, qa = q * m;

        // --- dot products with 4 accumulators each (ILP instead of one latency chain) ---
        double a0=0,a1=0,a2=0,a3=0;   // app parts
        double b0=0,b1=0,b2=0,b3=0;   // aqq parts
        double c0=0,c1=0,c2=0,c3=0;   // apq parts
        int i = 0;
        int limit = m - (m & 3);      // largest multiple of 4
        for (; i < limit; i += 4) {
            double x0=a[pa+i],   y0=a[qa+i];
            double x1=a[pa+i+1], y1=a[qa+i+1];
            double x2=a[pa+i+2], y2=a[qa+i+2];
            double x3=a[pa+i+3], y3=a[qa+i+3];
            a0+=x0*x0; a1+=x1*x1; a2+=x2*x2; a3+=x3*x3;
            b0+=y0*y0; b1+=y1*y1; b2+=y2*y2; b3+=y3*y3;
            c0+=x0*y0; c1+=x1*y1; c2+=x2*y2; c3+=x3*y3;
        }
        double app=(a0+a1)+(a2+a3);
        double aqq=(b0+b1)+(b2+b3);
        double apq=(c0+c1)+(c2+c3);
        for (; i < m; i++) {          // remainder (< 4 elements)
            double x=a[pa+i], y=a[qa+i];
            app+=x*x; aqq+=y*y; apq+=x*y;
        }

        if (Math.abs(apq) <= eps * Math.sqrt(app * aqq)) return false;

        double tau = (aqq - app) / (2.0 * apq);
        // sign(tau) has to be +1 at tau == 0, which is what Math.signum does NOT
        // give: it returns 0.0 there, so t, s and the whole rotation collapse to
        // the identity. Two columns of exactly equal norm - duplicated columns,
        // or any symmetric pair - then never get orthogonalized and the sweeps
        // run to the limit without converging. Same form as
        // SymmetricJacobiEigen, including its overflow guard.
        double at  = Math.abs(tau);
        double t   = (at > 1e150) ? 0.5 / tau
                : ((tau >= 0.0) ? 1.0 : -1.0) / (at + Math.sqrt(1.0 + tau * tau));
        double c   = 1.0 / Math.sqrt(1.0 + t * t);
        double s   = t * c;

        // --- update columns p,q of A (unroll x4, mainly to cut loop overhead) ---
        i = 0;
        for (; i < limit; i += 4) {
            double x0=a[pa+i],   y0=a[qa+i];
            double x1=a[pa+i+1], y1=a[qa+i+1];
            double x2=a[pa+i+2], y2=a[qa+i+2];
            double x3=a[pa+i+3], y3=a[qa+i+3];
            a[pa+i]  = c*x0 - s*y0;  a[qa+i]  = s*x0 + c*y0;
            a[pa+i+1]= c*x1 - s*y1;  a[qa+i+1]= s*x1 + c*y1;
            a[pa+i+2]= c*x2 - s*y2;  a[qa+i+2]= s*x2 + c*y2;
            a[pa+i+3]= c*x3 - s*y3;  a[qa+i+3]= s*x3 + c*y3;
        }
        for (; i < m; i++) {
            double x=a[pa+i], y=a[qa+i];
            a[pa+i]=c*x - s*y; a[qa+i]=s*x + c*y;
        }

        // accumulate V -- n is small, unrolling rarely pays off here.
        final int pv = p * n, qv = q * n;
        for (int k = 0; k < n; k++) {
            double x = v[pv+k], y = v[qv+k];
            v[pv+k] = c*x - s*y;
            v[qv+k] = s*x + c*y;
        }
        return true;
   }

    /**
     * Circle method: for even nEven -> (nEven-1) rounds with nEven/2 column-disjoint
     * pairs each; every pair from 0..nEven-1 appears exactly once.
     */
    private static int[][][] roundRobin(int nEven) {
        int[] arr = new int[nEven];
        for (int i = 0; i < nEven; i++) arr[i] = i;

        int rounds = nEven - 1, perRound = nEven / 2;
        int[][][] schedule = new int[rounds][perRound][2];

        for (int r = 0; r < rounds; r++) {
            for (int i = 0; i < perRound; i++) {
                int x = arr[i], y = arr[nEven - 1 - i];
                schedule[r][i][0] = Math.min(x, y);
                schedule[r][i][1] = Math.max(x, y);
            }
            // rotate arr[1..] by one, keep arr[0] fixed.
            int last = arr[nEven - 1];
            for (int i = nEven - 1; i > 1; i--) arr[i] = arr[i - 1];
            arr[1] = last;
        }
        return schedule;
    }

    /** Sort descending; carry the U and V column blocks along physically. */
    private static void sortDescending(double[] sigma, double[] u, int m, double[] v, int n) {
        for (int i = 0; i < n - 1; i++) {
            int max = i;
            for (int j = i + 1; j < n; j++) if (sigma[j] > sigma[max]) max = j;
            if (max != i) {
                double t = sigma[i]; sigma[i] = sigma[max]; sigma[max] = t;
                swapBlock(u, i * m, max * m, m);
                swapBlock(v, i * n, max * n, n);
            }
        }
    }

    private static void swapBlock(double[] arr, int offA, int offB, int len) {
        for (int k = 0; k < len; k++) {
            double t = arr[offA + k]; arr[offA + k] = arr[offB + k]; arr[offB + k] = t;
        }
    }

    /**
     * ||column||_2 without over-/underflow (BLAS dnrm2 style with a running scale).
     * @param off start offset of the column in the flat array, m elements from there.
     */
    private static double safeColumnNorm(double[] a, int off, int m) {
        double scale = 0.0, ssq = 1.0;
        for (int i = 0; i < m; i++) {
            double x = a[off + i];
            if (x != 0.0) {
                double ax = Math.abs(x);
                if (scale < ax) {
                    double r = scale / ax;
                    ssq = 1.0 + ssq * r * r;
                    scale = ax;
                } else {
                    double r = ax / scale;
                    ssq += r * r;
                }
            }
        }
        return scale * Math.sqrt(ssq);
    }

    // ---------------------------------------------------------------------
    // Reconstruction & verification
    // ---------------------------------------------------------------------

    /**
     * Reconstructs {@code A' = U * diag(sigma) * V^T}, column-major.
     *
     * @param r decomposition result
     * @return reconstructed matrix in column-major layout
     */
    public static double[] reconstruct(Result r) {
        int m = r.m, n = r.n;
        double[] A = new double[m * n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++) {
                double sum = 0;
                for (int k = 0; k < n; k++)
                    sum += r.U[k * m + i] * r.sigma[k] * r.V[k * n + j];
                //          U_{ik}            sigma_k     V_{jk} = V[k*n + j]
                A[j * m + i] = sum;
            }
        return A;
    }

    /**
     * Computes {@code max |A - U S V^T|}.
     *
     * @param A original matrix in column-major layout
     * @param r decomposition result
     * @return maximum absolute reconstruction error
     */
    public static double reconstructionError(double[] A, Result r) {
        double[] rec = reconstruct(r);
        double max = 0;
        for (int idx = 0; idx < A.length; idx++)
            max = Math.max(max, Math.abs(A[idx] - rec[idx]));
        return max;
    }

    /**
     * Computes {@code max |(Q^T Q) - I|} for {@code Q} (rows x cols), column-major {@code Q[j*rows + i]}.
     *
     * @param Q matrix in column-major layout
     * @param rows number of rows of {@code Q}
     * @param cols number of columns of {@code Q}
     * @return maximum absolute orthonormality error
     */
    public static double orthonormalityError(double[] Q, int rows, int cols) {
        double max = 0;
        for (int p = 0; p < cols; p++)
            for (int q = 0; q < cols; q++) {
                double dot = 0;
                for (int i = 0; i < rows; i++)
                    dot += Q[p * rows + i] * Q[q * rows + i];
                max = Math.max(max, Math.abs(dot - (p == q ? 1.0 : 0.0)));
            }
        return max;
    }

    // ---------------------------------------------------------------------
    // Demo
    // ---------------------------------------------------------------------

    public static void main(String[] args) {
        // 1) Small fixed example (m=4, n=3), written row-wise...
        double[][] rows = {
            { 1,  2,  3},
            { 4,  5,  6},
            { 7,  8, 10},
            { 2,  0,  1},
        };
        int m = rows.length, n = rows[0].length;
        double[] A = new double[m * n];                 // ...pack column-major
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A[j * m + i] = rows[i][j];

        report("Small 4x3 example", A, m, n);

        // 2) Larger, deterministically generated example -> runs in parallel.
        int M = 400, N = 250;
        double[] B = new double[M * N];
        long lcg = 123456789L;                          // simple LCG, no Math.random
        for (int idx = 0; idx < B.length; idx++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            B[idx] = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
        }
        report("Large " + M + "x" + N + " example", B, M, N);
    }

    private static void report(String label, double[] A, int m, int n) {
        FlatParallelJacobiSVD svd = new FlatParallelJacobiSVD();
        long t0 = System.nanoTime();
        Result r = svd.decompose(A, m, n);
        double ms = (System.nanoTime() - t0) / 1e6;

        double recErr = reconstructionError(A, r);
        double uErr   = orthonormalityError(r.U, m, n);
        double vErr   = orthonormalityError(r.V, n, n);

        System.out.println("=== " + label + " ===");
        System.out.printf("  time                     = %.2f ms%n", ms);
        System.out.printf("  converged                = %b%n", r.converged);
        System.out.printf("  max|A - U S V^T|         = %.3e%n", recErr);
        System.out.printf("  max|U^T U - I|           = %.3e%n", uErr);
        System.out.printf("  max|V^T V - I|           = %.3e%n", vErr);
        int show = Math.min(6, n);
        System.out.print  ("  sigma[0.." + (show - 1) + "]           = [");
        for (int j = 0; j < show; j++)
            System.out.printf("%s%.5g", j == 0 ? "" : ", ", r.sigma[j]);
        System.out.println(n > show ? ", ...]" : "]");
        boolean ok = recErr < 1e-9 && uErr < 1e-9 && vErr < 1e-9;
        System.out.println("  " + (ok ? ">>> OK" : ">>> FAILED") + "\n");
    }
}
