package math.linalg;

import java.util.stream.IntStream;

/**
 * Two-sided (classical) Jacobi eigenvalue decomposition for real symmetric
 * matrices with
 * <ul>
 * <li>flat, column-major storage (BLAS/LAPACK layout: {@code A[j*n + i]}) -
 * note that for a symmetric matrix row-major and column-major storage are
 * identical, so the caller does not have to care,</li>
 * <li>parallel round-robin ordering (Brent-Luk): per round {@code n/2}
 * index-disjoint rotations, so that every column belongs to exactly one of
 * them. Each rotation can therefore apply both sides of its two-sided update on
 * the two columns it owns, which lets a whole round run as a single parallel
 * dispatch without a barrier.</li>
 * </ul>
 *
 * Decomposes a symmetric {@code n x n} matrix {@code A} into
 * {@code A = V * diag(lambda) * V^T} with an orthogonal {@code V}. Column
 * {@code j} of {@code V} is the eigenvector belonging to {@code lambda[j]}, and
 * the eigenvalues are returned in <em>descending</em> order.
 * <p>
 * The method is backward stable and computes small eigenvalues with high
 * relative accuracy, but for large {@code n} it is slower than a tridiagonal
 * QL/QR approach. It is a good fit for small to medium sized matrices and for
 * cases where accuracy matters more than raw speed.
 *
 * @since 1.5.1
 */
public final class SymmetricJacobiEigen {

    /** Result; {@code V} stored column-major. */
    public static final class Result {
        /** {@code n x n} orthogonal eigenvector matrix, {@code V[j*n + i]}. */
        public final double[] V;
        /** {@code n} eigenvalues in descending order. */
        public final double[] lambda;
        /** the dimension of the decomposed matrix. */
        public final int n;
        /** false if maxSweeps was exhausted (result may be inaccurate). */
        public final boolean converged;

        Result(double[] V, double[] lambda, int n, boolean converged) {
            this.V = V; this.lambda = lambda; this.n = n;
            this.converged = converged;
        }

        /**
         * The eigenvector belonging to {@code lambda[j]} as a freshly allocated
         * array of length {@code n}.
         *
         * @param j index of the eigenvalue / eigenvector
         * @return a copy of column {@code j} of {@code V}
         */
        public double[] eigenvector(int j) {
            if (j < 0 || j >= n) throw new IndexOutOfBoundsException("j: " + j);
            double[] col = new double[n];
            System.arraycopy(V, j * n, col, 0, n);
            return col;
        }
    }

    private final double eps;               // relative off-diagonal threshold
    private final int    maxSweeps;         // safety upper bound
    private final long   parallelThreshold; // go parallel once n*n >= this size

    public SymmetricJacobiEigen() {
        this(1e-15, 60, 1L << 16);          // default: parallel from ~65k elements
    }

    public SymmetricJacobiEigen(double eps, int maxSweeps, long parallelThreshold) {
        this.eps = eps;
        this.maxSweeps = maxSweeps;
        this.parallelThreshold = parallelThreshold;
    }

    /**
     * @param A column-major packed, length n*n, element (i,j) = A[j*n + i]. Left unchanged.
     * @param n the dimension of the (symmetric) matrix
     * @return decomposition result containing {@code V}, the eigenvalues, and the convergence flag
     */
    public Result decompose(double[] A, int n) {
        checkArgs(A, n);
        return decomposeInPlace(A.clone(), n);
    }

    /**
     * In-place variant: {@code a} is used as the working array and is overwritten (it
     * becomes the diagonalized matrix). Saves one n*n copy versus {@link #decompose}.
     * <p>
     * A matrix that is only symmetric up to rounding is silently symmetrized; for a
     * matrix that is exactly symmetric this is a no-op. Use
     * {@link #symmetryError(double[], int)} beforehand if a badly non-symmetric input
     * should be treated as a caller error.
     *
     * @param a column-major packed, length n*n, element (i,j) = a[j*n + i]. Destroyed.
     * @param n the dimension of the (symmetric) matrix
     * @return decomposition result containing {@code V}, the eigenvalues, and the convergence flag
     */
    public Result decomposeInPlace(double[] a, int n) {
        checkArgs(a, n);

        final double[] v = new double[n * n];
        for (int j = 0; j < n; j++) v[j * n + j] = 1.0;   // V = I

        if (n == 1) {
            return new Result(v, new double[] { a[0] }, 1, true);
        }

        // Both passes below assume that both triangles agree exactly.
        symmetrize(a, n);

        // --- global scaling into a safe range (largest element ~ 1) ---
        // Guards the (aqq - app) difference in prepare() against overflow. Scalar
        // scaling is exactly compatible: eig(c*A) = c * eig(A), same eigenvectors.
        double maxAbs = 0.0;
        for (int idx = 0; idx < a.length; idx++) maxAbs = Math.max(maxAbs, Math.abs(a[idx]));

        double scale = 1.0;
        if (maxAbs > 0.0 && maxAbs < Double.POSITIVE_INFINITY) {
            // exact: only shift the exponent (power of two -> no rounding error)
            scale = Math.scalb(1.0, -Math.getExponent(maxAbs));
            if (scale != 1.0)
                for (int idx = 0; idx < a.length; idx++) a[idx] *= scale;
        }

        // Round-robin schedule (padded up to the next even n).
        final int nEven = (n % 2 == 0) ? n : n + 1;
        final int[][][] schedule = roundRobin(nEven);
        final boolean parallel = (long) n * n >= parallelThreshold;

        // Per-round scratch for the rotations of the current round.
        final int maxPairs = nEven / 2;
        final int[]    pp = new int[maxPairs];
        final int[]    qq = new int[maxPairs];
        final double[] cc = new double[maxPairs];
        final double[] ss = new double[maxPairs];
        final double[] dp = new double[maxPairs];   // exact new a(p,p)
        final double[] dq = new double[maxPairs];   // exact new a(q,q)
        // Per-slot view of the same round: the two columns owned by slot s and the
        // index of its rotation in the compact arrays above (-1 if it was skipped).
        final int[] colP = new int[maxPairs];
        final int[] colQ = new int[maxPairs];
        final int[] act  = new int[maxPairs];

        boolean didConverge = false;
        for (int sweep = 0; sweep < maxSweeps; sweep++) {
            boolean rotated = false;

            for (int[][] round : schedule) {
                // All pairs of ONE round are index-disjoint. A two-sided rotation
                // touches rows p,q AND columns p,q, so two disjoint pairs still
                // collide in the elements (p_i, p_j). The update is therefore split
                // into a right pass (columns) and a left pass (rows) - but both are
                // driven by the slot that owns the columns, so no barrier is needed
                // between them; see applyRound.
                final int cnt = prepare(a, n, round, pp, qq, cc, ss, dp, dq, colP, colQ, act);
                if (cnt == 0) continue;
                rotated = true;

                applyRound(a, v, n, cnt, pp, qq, cc, ss, colP, colQ, act, maxPairs, parallel);

                // Only rotation (p,q) touches the 2x2 block (p,q) - restore it from
                // the closed form so that the annihilated entries are exactly zero.
                for (int i = 0; i < cnt; i++) {
                    final int p = pp[i], q = qq[i];
                    a[p * n + p] = dp[i];
                    a[q * n + q] = dq[i];
                    a[q * n + p] = 0.0;
                    a[p * n + q] = 0.0;
                }
            }
            if (!rotated) { didConverge = true; break; }
        }
        // Convergence guard: never return a silently non-converged result unnoticed.
        if (!didConverge)
            System.err.println("[SymmetricJacobiEigen] WARNING: not converged after "
                    + maxSweeps + " sweeps (n=" + n
                    + "); result may be inaccurate. Consider raising maxSweeps.");

        final double[] lambda = new double[n];
        final double inv = (scale != 1.0) ? 1.0 / scale : 1.0;  // undo the global scaling
        for (int j = 0; j < n; j++) lambda[j] = a[j * n + j] * inv;

        sortDescending(lambda, v, n);
        return new Result(v, lambda, n, didConverge);
    }

    /**
     * Computes the rotations of one round from the current 2x2 diagonal blocks and
     * collects them in the scratch arrays. Pairs that are already orthogonal enough
     * are skipped. Returns the number of rotations to apply.
     * <p>
     * Fills two views of the same round: the compact arrays {@code pp/qq/cc/ss/dp/dq}
     * hold the {@code cnt} rotations that actually have to be applied, while
     * {@code colP/colQ/act} keep one entry per slot of the schedule. A slot whose
     * pair was skipped still owns its columns - they need the left rotations of the
     * other pairs - which is what {@code act[s] < 0} encodes.
     */
    private int prepare(double[] a, int n, int[][] round,
                        int[] pp, int[] qq, double[] cc, double[] ss, double[] dp, double[] dq,
                        int[] colP, int[] colQ, int[] act) {
        int cnt = 0;
        for (int slot = 0; slot < round.length; slot++) {
            final int p = round[slot][0], q = round[slot][1];
            colP[slot] = (p < n) ? p : -1;                  // dummy index -> no column
            colQ[slot] = (q < n) ? q : -1;
            act[slot] = -1;
            if (p >= n || q >= n) continue;                 // dummy index -> skip
            final double apq = a[q * n + p];
            if (apq == 0.0) continue;
            final double app = a[p * n + p];
            final double aqq = a[q * n + q];
            // relative criterion, factored to avoid overflow in |app| * |aqq|
            if (Math.abs(apq) <= eps * Math.sqrt(Math.abs(app)) * Math.sqrt(Math.abs(aqq))) continue;

            // t = tan(phi) as the root of t^2 + 2*theta*t - 1 = 0 with the smaller
            // magnitude (the rotation closer to the identity).
            final double theta = (aqq - app) / (2.0 * apq);
            final double at = Math.abs(theta);
            final double t = (at > 1e150) ? 0.5 / theta     // 1 + theta^2 would overflow
                    : ((theta >= 0.0) ? 1.0 : -1.0) / (at + Math.sqrt(1.0 + theta * theta));
            final double c = 1.0 / Math.sqrt(1.0 + t * t);
            final double s = t * c;
            if (s == 0.0) continue;                         // underflow: rotation is the identity

            act[slot] = cnt;
            pp[cnt] = p; qq[cnt] = q; cc[cnt] = c; ss[cnt] = s;
            dp[cnt] = app - t * apq;
            dq[cnt] = aqq + t * apq;
            cnt++;
        }
        return cnt;
    }

    /**
     * Applies one whole round, {@code A <- J^T*A*J} and {@code V <- V*J}, in a
     * single parallel pass over the slots of the schedule.
     * <p>
     * The rotations of a round are index-disjoint, so every column of the matrix
     * belongs to exactly one slot. Slot s therefore owns its two columns outright
     * and performs both halves of the two-sided update on them without ever
     * touching anyone else's memory:
     * <ol>
     * <li>right pass - mix the two owned columns of {@code A} (and of {@code V})
     *     against each other, which is {@code A <- A*J} restricted to those two
     *     columns;</li>
     * <li>left pass - apply <em>all</em> rotations of the round to the entries
     *     <em>inside</em> each owned column, which is {@code A <- J^T*A} restricted
     *     to the same two columns.</li>
     * </ol>
     * The left pass only ever mixes entries within one column, and the right pass
     * only needs the pre-round values of the two owned columns, so no barrier is
     * required between the two - which is what lets them be fused. Fusing them
     * halves the number of parallel dispatches per sweep and, more importantly,
     * keeps the two columns (2*n doubles) in cache between the passes instead of
     * streaming the whole matrix through memory a second time.
     * <p>
     * A slot whose own pair was skipped as already converged ({@code act[s] < 0})
     * still has to run the left pass on its columns.
     */
    private static void applyRound(double[] a, double[] v, int n, int cnt,
                                   int[] pp, int[] qq, double[] cc, double[] ss,
                                   int[] colP, int[] colQ, int[] act, int slots, boolean parallel) {
        IntStream stream = IntStream.range(0, slots);
        if (parallel) stream = stream.parallel();
        stream.forEach(slot -> {
            final int p = colP[slot], q = colQ[slot];
            final int r = act[slot];
            if (r >= 0) {
                final int pOff = p * n, qOff = q * n;
                final double c = cc[r], sn = ss[r];
                rotateColumns(a, pOff, qOff, n, c, sn);
                rotateColumns(v, pOff, qOff, n, c, sn);
            }
            if (p >= 0) rotateWithinColumn(a, p * n, cnt, pp, qq, cc, ss);
            if (q >= 0) rotateWithinColumn(a, q * n, cnt, pp, qq, cc, ss);
        });
    }

    /**
     * Left pass for a single column: applies every rotation of the round to the
     * entry pair (p,q) inside that column. Row p of a column-major matrix has
     * stride n, so doing this one whole column at a time - rather than one pair at
     * a time - keeps all accesses inside a single contiguous block.
     */
    private static void rotateWithinColumn(double[] a, int off, int cnt,
                                           int[] pp, int[] qq, double[] cc, double[] ss) {
        for (int i = 0; i < cnt; i++) {
            final int ip = off + pp[i], iq = off + qq[i];
            final double x = a[ip], y = a[iq];
            final double c = cc[i], s = ss[i];
            a[ip] = c * x - s * y;
            a[iq] = s * x + c * y;
        }
    }

    /** new col p = c*col p - s*col q, new col q = s*col p + c*col q. */
    private static void rotateColumns(double[] m, int pOff, int qOff, int n, double c, double s) {
        int k = 0;
        final int limit = n - (n & 3);      // largest multiple of 4
        for (; k < limit; k += 4) {
            double x0=m[pOff+k],   y0=m[qOff+k];
            double x1=m[pOff+k+1], y1=m[qOff+k+1];
            double x2=m[pOff+k+2], y2=m[qOff+k+2];
            double x3=m[pOff+k+3], y3=m[qOff+k+3];
            m[pOff+k]  = c*x0 - s*y0;  m[qOff+k]  = s*x0 + c*y0;
            m[pOff+k+1]= c*x1 - s*y1;  m[qOff+k+1]= s*x1 + c*y1;
            m[pOff+k+2]= c*x2 - s*y2;  m[qOff+k+2]= s*x2 + c*y2;
            m[pOff+k+3]= c*x3 - s*y3;  m[qOff+k+3]= s*x3 + c*y3;
        }
        for (; k < n; k++) {                // remainder (< 4 elements)
            double x = m[pOff+k], y = m[qOff+k];
            m[pOff+k] = c*x - s*y;
            m[qOff+k] = s*x + c*y;
        }
    }

    /**
     * Circle method: for even nEven -> (nEven-1) rounds with nEven/2 index-disjoint
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

    /** Averages both triangles; exact (and thus a no-op) for symmetric input. */
    private static void symmetrize(double[] a, int n) {
        for (int j = 0; j < n; j++)
            for (int i = j + 1; i < n; i++) {
                double mean = 0.5 * (a[j * n + i] + a[i * n + j]);
                a[j * n + i] = mean;
                a[i * n + j] = mean;
            }
    }

    /** Sort descending; carry the V column blocks along physically. */
    private static void sortDescending(double[] lambda, double[] v, int n) {
        for (int i = 0; i < n - 1; i++) {
            int max = i;
            for (int j = i + 1; j < n; j++) if (lambda[j] > lambda[max]) max = j;
            if (max != i) {
                double t = lambda[i]; lambda[i] = lambda[max]; lambda[max] = t;
                swapBlock(v, i * n, max * n, n);
            }
        }
    }

    private static void swapBlock(double[] arr, int offA, int offB, int len) {
        for (int k = 0; k < len; k++) {
            double t = arr[offA + k]; arr[offA + k] = arr[offB + k]; arr[offB + k] = t;
        }
    }

    private static void checkArgs(double[] a, int n) {
        if (n < 1) throw new IllegalArgumentException("expected n >= 1, got " + n);
        if (a.length != (long) n * n) throw new IllegalArgumentException("length != n*n");
    }

    // ---------------------------------------------------------------------
    // Reconstruction & verification
    // ---------------------------------------------------------------------

    /**
     * Reconstructs {@code A' = V * diag(lambda) * V^T}, column-major.
     *
     * @param r decomposition result
     * @return reconstructed matrix in column-major layout
     */
    public static double[] reconstruct(Result r) {
        final int n = r.n;
        double[] A = new double[n * n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int k = 0; k < n; k++)
                    sum += r.V[k * n + i] * r.lambda[k] * r.V[k * n + j];
                //          V_{ik}            lambda_k    V_{jk} = V[k*n + j]
                A[j * n + i] = sum;
            }
        return A;
    }

    /**
     * Computes {@code max |A - V L V^T|}.
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
     * Computes {@code max |A*v_j - lambda_j*v_j|} over all eigenpairs, i.e. how well
     * the eigenvalue equation itself is satisfied.
     *
     * @param A original matrix in column-major layout
     * @param r decomposition result
     * @return maximum absolute residual of the eigenvalue equation
     */
    public static double residualError(double[] A, Result r) {
        final int n = r.n;
        double max = 0;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int k = 0; k < n; k++)
                    sum += A[k * n + i] * r.V[j * n + k];
                max = Math.max(max, Math.abs(sum - r.lambda[j] * r.V[j * n + i]));
            }
        return max;
    }

    /**
     * Computes {@code max |A - A^T|}, i.e. how far the input deviates from being
     * symmetric.
     *
     * @param A matrix in column-major layout
     * @param n the dimension of the matrix
     * @return maximum absolute deviation from symmetry
     */
    public static double symmetryError(double[] A, int n) {
        double max = 0;
        for (int j = 0; j < n; j++)
            for (int i = j + 1; i < n; i++)
                max = Math.max(max, Math.abs(A[j * n + i] - A[i * n + j]));
        return max;
    }

    // ---------------------------------------------------------------------
    // Demo
    // ---------------------------------------------------------------------

    public static void main(String[] args) {
        // 1) Small fixed example (n=4), written row-wise (same thing here)...
        double[][] rows = {
            { 4,  1, -2,  2},
            { 1,  2,  0,  1},
            {-2,  0,  3, -2},
            { 2,  1, -2, -1},
        };
        int n = rows.length;
        double[] A = new double[n * n];                 // ...packed column-major
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[j * n + i] = rows[i][j];

        report("Small 4x4 example", A, n);

        // 2) Larger, deterministically generated example -> runs in parallel.
        int N = 300;
        double[] B = new double[N * N];
        long lcg = 123456789L;                          // simple LCG, no Math.random
        for (int j = 0; j < N; j++)
            for (int i = j; i < N; i++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double x = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
                B[j * N + i] = x;
                B[i * N + j] = x;
            }
        report("Large " + N + "x" + N + " example", B, N);
    }

    private static void report(String label, double[] A, int n) {
        SymmetricJacobiEigen eig = new SymmetricJacobiEigen();
        long t0 = System.nanoTime();
        Result r = eig.decompose(A, n);
        double ms = (System.nanoTime() - t0) / 1e6;

        double recErr = reconstructionError(A, r);
        double resErr = residualError(A, r);
        double vErr   = FlatParallelJacobiSVD.orthonormalityError(r.V, n, n);

        double trace = 0.0, sumLambda = 0.0;
        for (int i = 0; i < n; i++) { trace += A[i * n + i]; sumLambda += r.lambda[i]; }

        System.out.println("=== " + label + " ===");
        System.out.printf("  time                     = %.2f ms%n", ms);
        System.out.printf("  converged                = %b%n", r.converged);
        System.out.printf("  max|A - V L V^T|         = %.3e%n", recErr);
        System.out.printf("  max|A v - lambda v|      = %.3e%n", resErr);
        System.out.printf("  max|V^T V - I|           = %.3e%n", vErr);
        System.out.printf("  |trace - sum(lambda)|    = %.3e%n", Math.abs(trace - sumLambda));
        int show = Math.min(6, n);
        System.out.print  ("  lambda[0.." + (show - 1) + "]          = [");
        for (int j = 0; j < show; j++)
            System.out.printf("%s%.5g", j == 0 ? "" : ", ", r.lambda[j]);
        System.out.println(n > show ? ", ...]" : "]");
        boolean ok = recErr < 1e-9 && resErr < 1e-9 && vErr < 1e-9;
        System.out.println("  " + (ok ? ">>> OK" : ">>> FAILED") + "\n");
    }
}
