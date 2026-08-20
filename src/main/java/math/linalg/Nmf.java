package math.linalg;

import java.util.Arrays;

import math.gemm.Dgemm;
import math.gemm.Trans;
import math.rng.SplitMix64;

/**
 * Non-negative matrix factorization: approximates a non-negative {@code X}
 * ({@code m x n}) by a product of two non-negative factors {@code W} ({@code m x k})
 * and {@code H} ({@code k x n}) of a chosen rank {@code k}, minimizing
 *
 * <pre>
 * 0.5 ||X - W H||_F^2  +  penalty(W)  +  penalty(H)
 * </pre>
 *
 * <p>It is the additive counterpart to {@link JacobiPCA}: principal components are
 * orthogonal and carry signs, so a sample is a sum of positive and negative
 * contributions that cancel; the factors here are non-negative throughout, so a
 * sample is a sum of parts that only add. Spectral unmixing, topic extraction and
 * source separation want the second.</p>
 *
 * <p>The engine is hierarchical alternating least squares (HALS): one rank-one
 * factor at a time, each coordinate minimized in closed form, which is exact
 * because the subproblem for a single entry is scalar and non-negativity is
 * therefore a {@code max(0, .)}. The two {@code O(m n k)} products run through
 * {@link Dgemm} and inherit its parallelism; nothing else in a sweep exceeds
 * {@code O(k^2 (m + n))}.</p>
 *
 * <p>The problem is <em>not</em> convex, so what is returned is a stationary point,
 * not a global minimum, and it depends on the starting point. Two starts are
 * offered and both are reproducible: {@link Init#RANDOM} scales its entries to the
 * data so that {@code W H} begins at the right magnitude, and {@link Init#NNDSVD}
 * is derived from the leading singular triplets and does not depend on a seed at
 * all.</p>
 *
 * <p>The factorization is only determined up to a positive diagonal
 * {@code W H = (W D)(D^-1 H)} and up to a permutation of the components. Both are
 * pinned on the way out: the columns of {@code W} come back with unit Euclidean
 * norm and the components are ordered by descending
 * {@link Result#componentEnergy}, so two runs are comparable. A penalty needs the
 * same normalization <em>during</em> the run, since otherwise an L1 penalty on one
 * factor is simply evaded by inflating the other; the class does that on its own
 * whenever a penalty is active.</p>
 *
 * @see <a href="https://en.wikipedia.org/wiki/Non-negative_matrix_factorization">Non-negative
 *      matrix factorization</a>
 * @since 1.5.2
 */
public final class Nmf {

    /** Default relative tolerance on the projected gradient. */
    private static final double DEFAULT_TOLERANCE = 1.0e-5;

    /** Default cap on the number of sweeps. */
    private static final int DEFAULT_MAX_ITERATIONS = 500;

    /** Default seed of the random start, fixed so that results are reproducible. */
    private static final long DEFAULT_SEED = 20260821L;

    /**
     * How the factors are seeded before the first sweep.
     *
     * @since 1.5.2
     */
    public enum Init {

        /**
         * Uniform entries scaled by {@code sqrt(mean(X) / k)}, so that the initial
         * product already has the magnitude of the data. Reproducible through the
         * seed handed to the constructor.
         */
        RANDOM,

        /**
         * Non-negative double singular value decomposition (Boutsidis and
         * Gallopoulos): the leading {@code k} singular triplets, each split into its
         * positive and its negative part with the heavier one kept. Deterministic and
         * markedly better than a random start at short iteration budgets, at the cost
         * of one {@link FlatParallelJacobiSVD} of the whole matrix.
         */
        NNDSVD
    }

    /**
     * Why the run ended.
     *
     * @since 1.5.2
     */
    public enum Status {

        /** The projected gradient fell to the requested fraction of its initial value. */
        TOLERANCE_REACHED(true),

        /** The sweep cap was exhausted first; the factors are usable but not converged. */
        TOO_MANY_ITERATIONS(false);

        private final boolean success;

        Status(boolean success) {
            this.success = success;
        }

        /**
         * Whether this status denotes a converged run.
         *
         * @return {@code true} for {@link #TOLERANCE_REACHED}
         */
        public boolean isSuccess() {
            return success;
        }
    }

    /**
     * An elastic net penalty on one factor, {@code lambda * (alpha * sum|a| +
     * 0.5 * (1 - alpha) * sum a^2)}, in the vocabulary of {@link Lasso}. Since the
     * factors are non-negative the L1 part is simply their sum, and it drives
     * entries to exactly zero.
     *
     * @since 1.5.2
     */
    public static final class Penalty {

        /** Overall strength; {@code 0} switches the penalty off. */
        public final double lambda;

        /** Mixing: {@code 1} is pure L1, {@code 0} is pure L2. */
        public final double alpha;

        private Penalty(double lambda, double alpha) {
            this.lambda = lambda;
            this.alpha = alpha;
        }

        /**
         * No penalty at all.
         *
         * @return a penalty of strength zero
         */
        public static Penalty none() {
            return new Penalty(0.0, 1.0);
        }

        /**
         * Pure L1, which makes the factor sparse with exact zeros.
         *
         * @param lambda strength, {@code >= 0} and finite
         * @return the penalty
         */
        public static Penalty lasso(double lambda) {
            return elasticNet(lambda, 1.0);
        }

        /**
         * Pure L2, which shrinks the factor without creating zeros.
         *
         * @param lambda strength, {@code >= 0} and finite
         * @return the penalty
         */
        public static Penalty ridge(double lambda) {
            return elasticNet(lambda, 0.0);
        }

        /**
         * A mixture of the two.
         *
         * @param lambda strength, {@code >= 0} and finite
         * @param alpha  mixing in {@code [0, 1]}, {@code 1} being pure L1
         * @return the penalty
         */
        public static Penalty elasticNet(double lambda, double alpha) {
            if (!(lambda >= 0.0) || Double.isInfinite(lambda)) {
                throw new IllegalArgumentException("lambda must be finite and >= 0 : " + lambda);
            }
            if (!(alpha >= 0.0) || !(alpha <= 1.0)) {
                throw new IllegalArgumentException("alpha must lie in [0, 1] : " + alpha);
            }
            return new Penalty(lambda, alpha);
        }

        boolean isActive() {
            return lambda > 0.0;
        }

        double l1() {
            return lambda * alpha;
        }

        double l2() {
            return lambda * (1.0 - alpha);
        }

        @Override
        public String toString() {
            if (!isActive()) {
                return "Penalty[none]";
            }
            return "Penalty[lambda=" + lambda + ", alpha=" + alpha + "]";
        }
    }

    /**
     * The factorization and what is known about it. Both factors are stored flat and
     * column-major, in the layout the rest of {@code math.linalg} uses.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** {@code m x k}, entry {@code (i, r)} at {@code W[r * m + i]}, columns of unit norm. */
        public final double[] W;

        /** {@code k x n}, entry {@code (r, j)} at {@code H[j * k + r]}. */
        public final double[] H;

        /** Rows of {@code X}. */
        public final int m;

        /** Columns of {@code X}. */
        public final int n;

        /** Rank of the factorization. */
        public final int k;

        /** {@code ||X - W H||_F / ||X||_F}, recomputed from the factors, not accumulated. */
        public final double relativeError;

        /** {@code 0.5 ||X - W H||_F^2} plus the penalties, the quantity actually minimized. */
        public final double objective;

        /**
         * {@code ||H_r||_2} per component, descending. Since the columns of {@code W}
         * have unit norm this is the weight the component carries in the product, the
         * counterpart of a singular value.
         */
        public final double[] componentEnergy;

        /** Number of components that collapsed to zero and contribute nothing. */
        public final int deadComponents;

        /** Norm of the projected gradient at the returned factors. */
        public final double projectedGradient;

        /**
         * {@link #projectedGradient} divided by its value at the starting point, which
         * is the quantity compared against the tolerance. Invariant under
         * {@code X -> c X} as long as no penalty is active; a penalty carries units and
         * breaks that.
         */
        public final double projectedGradientRatio;

        /** Completed sweeps. */
        public final int iterations;

        /** Why the run ended. */
        public final Status status;

        /** Shorthand for {@code status.isSuccess()}. */
        public final boolean converged;

        Result(double[] W, double[] H, int m, int n, int k, double relativeError, double objective,
                double[] componentEnergy, int deadComponents, double projectedGradient,
                double projectedGradientRatio, int iterations, Status status) {
            this.W = W;
            this.H = H;
            this.m = m;
            this.n = n;
            this.k = k;
            this.relativeError = relativeError;
            this.objective = objective;
            this.componentEnergy = componentEnergy;
            this.deadComponents = deadComponents;
            this.projectedGradient = projectedGradient;
            this.projectedGradientRatio = projectedGradientRatio;
            this.iterations = iterations;
            this.status = status;
            this.converged = status.isSuccess();
        }

        @Override
        public String toString() {
            return "Nmf.Result[" + m + "x" + n + ", k=" + k + ", relativeError=" + relativeError
                    + ", iterations=" + iterations + ", status=" + status + ", dead=" + deadComponents
                    + ", pgRatio=" + projectedGradientRatio + "]";
        }
    }

    private final double tolerance;
    private final int maxIterations;
    private final Init init;
    private final long seed;
    private final Penalty penaltyW;
    private final Penalty penaltyH;

    /**
     * An instance with the calibrated defaults: relative tolerance {@code 1e-5} on
     * the projected gradient, at most {@code 500} sweeps, a scaled random start with
     * a fixed seed and no penalty.
     */
    public Nmf() {
        this(DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATIONS, Init.RANDOM, DEFAULT_SEED, Penalty.none(),
                Penalty.none());
    }

    /**
     * @param tolerance     stop once the projected gradient has fallen to this
     *                      fraction of its value at the starting point. The relative
     *                      reconstruction error then lands within roughly this
     *                      fraction of what an exhaustive run would reach.
     * @param maxIterations cap on the sweeps; reaching it leaves
     *                      {@link Result#converged} at {@code false}
     * @param init          how to seed the factors
     * @param seed          seed of {@link Init#RANDOM}, ignored by {@link Init#NNDSVD}
     * @param penaltyW      penalty on {@code W}, never {@code null}
     * @param penaltyH      penalty on {@code H}, never {@code null}
     */
    public Nmf(double tolerance, int maxIterations, Init init, long seed, Penalty penaltyW,
            Penalty penaltyH) {
        if (!(tolerance >= 0.0)) {
            throw new IllegalArgumentException("tolerance must be >= 0 : " + tolerance);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException("maxIterations must be >= 1 : " + maxIterations);
        }
        if (init == null) {
            throw new IllegalArgumentException("init must not be null");
        }
        if (penaltyW == null || penaltyH == null) {
            throw new IllegalArgumentException("penalties must not be null");
        }
        this.tolerance = tolerance;
        this.maxIterations = maxIterations;
        this.init = init;
        this.seed = seed;
        this.penaltyW = penaltyW;
        this.penaltyH = penaltyH;
    }

    /**
     * Factorizes {@code X}.
     *
     * @param X column-major packed, length {@code m * n}, entry {@code (i, j)} at
     *          {@code X[j * m + i]}; every entry finite and {@code >= 0}. Left
     *          unchanged.
     * @param m number of rows, {@code >= 1}
     * @param n number of columns, {@code >= 1}
     * @param k rank, {@code 1 <= k <= min(m, n)}
     * @return the factorization; check {@link Result#converged}
     */
    public Result factor(double[] X, int m, int n, int k) {
        check(X, m, n, k);
        return run(X, m, n, k);
    }

    /**
     * Factorizes {@code X}, taking the shape from the matrix.
     *
     * @param X the data, every entry finite and {@code >= 0}; left unchanged
     * @param k rank, {@code 1 <= k <= min(rows, columns)}
     * @return the factorization; check {@link Result#converged}
     */
    public Result factor(DMatrix X, int k) {
        if (X == null) {
            throw new IllegalArgumentException("X must not be null");
        }
        return factor(X.getArrayUnsafe(), X.numRows(), X.numColumns(), k);
    }

    // ------------------------------------------------------------------ loop --

    private Result run(double[] x, int m, int n, int k) {
        double[] w = new double[m * k];
        double[] h = new double[k * n];
        initialize(x, w, h, m, n, k);

        final double l1H = penaltyH.l1();
        final double l2H = penaltyH.l2();
        final double l1W = penaltyW.l1();
        final double l2W = penaltyW.l2();
        // An L1 penalty on one factor is meaningless while the other is free to
        // absorb any scale, so the normalization moves into the loop as soon as a
        // penalty is active. Without one it is a no-op on the trajectory and is
        // deferred to the end, which keeps the stopping measure scale invariant.
        final boolean normalizeEverySweep = penaltyW.isActive() || penaltyH.isActive();

        double[] wtw = new double[k * k];
        double[] wtx = new double[k * n];
        double[] hht = new double[k * k];
        double[] xht = new double[m * k];

        gramW(w, m, k, wtw);
        crossLeft(w, x, m, n, k, wtx);
        gramH(h, k, n, hht);
        crossRight(x, h, m, n, k, xht);
        final double pg0 = projectedGradient(w, h, wtw, wtx, hht, xht, m, n, k, l1H, l2H, l1W, l2W);

        int iterations = 0;
        double pg = pg0;
        Status status = Status.TOO_MANY_ITERATIONS;
        if (!(pg0 > 0.0)) {
            // already stationary, which is also the answer for an all zero X
            status = Status.TOLERANCE_REACHED;
        } else {
            for (int it = 1; it <= maxIterations; it++) {
                sweepH(h, wtw, wtx, n, k, l1H, l2H);
                if (normalizeEverySweep) {
                    normalize(w, h, m, n, k);
                }
                gramH(h, k, n, hht);
                crossRight(x, h, m, n, k, xht);
                sweepW(w, hht, xht, m, k, l1W, l2W);
                gramW(w, m, k, wtw);
                crossLeft(w, x, m, n, k, wtx);
                iterations = it;
                pg = projectedGradient(w, h, wtw, wtx, hht, xht, m, n, k, l1H, l2H, l1W, l2W);
                if (pg <= tolerance * pg0) {
                    status = Status.TOLERANCE_REACHED;
                    break;
                }
            }
        }

        normalize(w, h, m, n, k);
        double[] energy = orderByEnergy(w, h, m, n, k);
        int dead = 0;
        for (int r = 0; r < k; r++) {
            double columnNorm = 0.0;
            for (int i = 0; i < m; i++) {
                columnNorm += w[r * m + i] * w[r * m + i];
            }
            if (!(energy[r] > 0.0) || !(columnNorm > 0.0)) {
                dead++;
            }
        }
        double normSq = normSquared(x);
        double residual = residualSquared(x, w, h, m, n, k);
        double relative = (normSq > 0.0) ? Math.sqrt(residual / normSq) : 0.0;
        double objective = 0.5 * residual + penaltyValue(w, penaltyW) + penaltyValue(h, penaltyH);
        return new Result(w, h, m, n, k, relative, objective, energy, dead, pg, pg / pg0, iterations,
                status);
    }

    // ---------------------------------------------------------------- sweeps --

    /**
     * {@code H[r][j] <- max(0, (s + d h - l1) / (d + l2))} with
     * {@code s = (W'X - W'W H)[r][j]} and {@code d = (W'W)[r][r]}, which is the exact
     * minimizer of the objective in that single entry.
     */
    private static void sweepH(double[] h, double[] wtw, double[] wtx, int n, int k, double l1,
            double l2) {
        for (int r = 0; r < k; r++) {
            double d = wtw[r * k + r];
            if (!(d > 0.0)) {
                continue;
            }
            double den = d + l2;
            for (int j = 0; j < n; j++) {
                int base = j * k;
                double s = wtx[base + r];
                for (int c = 0; c < k; c++) {
                    s -= wtw[c * k + r] * h[base + c];
                }
                double v = (s + d * h[base + r] - l1) / den;
                h[base + r] = (v > 0.0) ? v : 0.0;
            }
        }
    }

    /** The same for {@code W}, with {@code H H'} and {@code X H'} in place of the two. */
    private static void sweepW(double[] w, double[] hht, double[] xht, int m, int k, double l1,
            double l2) {
        for (int r = 0; r < k; r++) {
            double d = hht[r * k + r];
            if (!(d > 0.0)) {
                continue;
            }
            double den = d + l2;
            int cb = r * m;
            for (int i = 0; i < m; i++) {
                double s = xht[cb + i];
                for (int c = 0; c < k; c++) {
                    s -= hht[c * k + r] * w[c * m + i];
                }
                double v = (s + d * w[cb + i] - l1) / den;
                w[cb + i] = (v > 0.0) ? v : 0.0;
            }
        }
    }

    // --------------------------------------------------------------- kernels --

    /** {@code W'W}, {@code k x k}, symmetric. */
    private static void gramW(double[] w, int m, int k, double[] out) {
        for (int a = 0; a < k; a++) {
            for (int b = 0; b <= a; b++) {
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    s += w[a * m + i] * w[b * m + i];
                }
                out[a * k + b] = s;
                out[b * k + a] = s;
            }
        }
    }

    /** {@code H H'}, {@code k x k}, symmetric. */
    private static void gramH(double[] h, int k, int n, double[] out) {
        Arrays.fill(out, 0.0);
        for (int j = 0; j < n; j++) {
            int base = j * k;
            for (int a = 0; a < k; a++) {
                double va = h[base + a];
                if (va == 0.0) {
                    continue;
                }
                for (int b = 0; b <= a; b++) {
                    out[a * k + b] += va * h[base + b];
                }
            }
        }
        for (int a = 0; a < k; a++) {
            for (int b = 0; b < a; b++) {
                out[b * k + a] = out[a * k + b];
            }
        }
    }

    /** {@code W'X}, {@code k x n}. */
    private static void crossLeft(double[] w, double[] x, int m, int n, int k, double[] out) {
        Dgemm.dgemm(Trans.TRANS, Trans.NO_TRANS, k, n, m, 1.0, w, 0, m, x, 0, m, 0.0, out, 0, k);
    }

    /** {@code X H'}, {@code m x k}. */
    private static void crossRight(double[] x, double[] h, int m, int n, int k, double[] out) {
        Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m, k, n, 1.0, x, 0, m, h, 0, k, 0.0, out, 0, m);
    }

    // -------------------------------------------------------------- gradient --

    /**
     * Norm of the projected gradient of both factors: the plain gradient where an
     * entry is positive, only its negative part where the entry sits at zero. It
     * vanishes exactly at a stationary point of the constrained problem.
     */
    private static double projectedGradient(double[] w, double[] h, double[] wtw, double[] wtx,
            double[] hht, double[] xht, int m, int n, int k, double l1H, double l2H, double l1W,
            double l2W) {
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            int base = j * k;
            for (int r = 0; r < k; r++) {
                double g = -wtx[base + r] + l1H + l2H * h[base + r];
                for (int c = 0; c < k; c++) {
                    g += wtw[c * k + r] * h[base + c];
                }
                double p = (h[base + r] > 0.0) ? g : Math.min(0.0, g);
                s += p * p;
            }
        }
        for (int r = 0; r < k; r++) {
            int cb = r * m;
            for (int i = 0; i < m; i++) {
                double g = -xht[cb + i] + l1W + l2W * w[cb + i];
                for (int c = 0; c < k; c++) {
                    g += hht[c * k + r] * w[c * m + i];
                }
                double p = (w[cb + i] > 0.0) ? g : Math.min(0.0, g);
                s += p * p;
            }
        }
        return Math.sqrt(s);
    }

    // -------------------------------------------------------- scale and order --

    /** Scales the columns of {@code W} to unit norm and moves the factor into {@code H}. */
    private static void normalize(double[] w, double[] h, int m, int n, int k) {
        for (int r = 0; r < k; r++) {
            int cb = r * m;
            double s = 0.0;
            for (int i = 0; i < m; i++) {
                s += w[cb + i] * w[cb + i];
            }
            double nrm = Math.sqrt(s);
            if (!(nrm > 0.0)) {
                continue;
            }
            double inv = 1.0 / nrm;
            for (int i = 0; i < m; i++) {
                w[cb + i] *= inv;
            }
            for (int j = 0; j < n; j++) {
                h[j * k + r] *= nrm;
            }
        }
    }

    /**
     * Sorts the components by descending {@code ||H_r||} and permutes both factors
     * accordingly.
     *
     * @return the energies, descending
     */
    private static double[] orderByEnergy(double[] w, double[] h, int m, int n, int k) {
        double[] energy = new double[k];
        for (int j = 0; j < n; j++) {
            int base = j * k;
            for (int r = 0; r < k; r++) {
                energy[r] += h[base + r] * h[base + r];
            }
        }
        for (int r = 0; r < k; r++) {
            energy[r] = Math.sqrt(energy[r]);
        }
        // insertion sort by descending energy; k is small and this keeps the
        // permutation on primitive arrays
        int[] order = new int[k];
        for (int r = 0; r < k; r++) {
            order[r] = r;
        }
        for (int r = 1; r < k; r++) {
            int pick = order[r];
            double value = energy[pick];
            int p = r - 1;
            while (p >= 0 && energy[order[p]] < value) {
                order[p + 1] = order[p];
                p--;
            }
            order[p + 1] = pick;
        }
        boolean identity = true;
        for (int r = 0; r < k; r++) {
            if (order[r] != r) {
                identity = false;
                break;
            }
        }
        if (identity) {
            return energy;
        }
        double[] w2 = new double[w.length];
        double[] h2 = new double[h.length];
        double[] e2 = new double[k];
        for (int r = 0; r < k; r++) {
            int src = order[r];
            System.arraycopy(w, src * m, w2, r * m, m);
            for (int j = 0; j < n; j++) {
                h2[j * k + r] = h[j * k + src];
            }
            e2[r] = energy[src];
        }
        System.arraycopy(w2, 0, w, 0, w.length);
        System.arraycopy(h2, 0, h, 0, h.length);
        return e2;
    }

    // ---------------------------------------------------------------- starts --

    private void initialize(double[] x, double[] w, double[] h, int m, int n, int k) {
        if (init == Init.NNDSVD) {
            nndsvd(x, w, h, m, n, k);
            // A zero column of W would freeze that component for good, so anything the
            // splitting left empty falls back to the random start for those entries.
            repairEmptyComponents(x, w, h, m, n, k);
        } else {
            randomStart(x, w, h, m, n, k, seed);
        }
    }

    /** Uniform entries scaled to {@code sqrt(mean(X) / k)}. */
    private static void randomStart(double[] x, double[] w, double[] h, int m, int n, int k,
            long seed) {
        double mean = 0.0;
        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }
        mean /= x.length;
        double scale = Math.sqrt(mean / k);
        if (!(scale > 0.0)) {
            scale = 1.0;
        }
        SplitMix64 rng = new SplitMix64(seed);
        for (int i = 0; i < w.length; i++) {
            w[i] = scale * (rng.nextDouble() + 0.1);
        }
        for (int i = 0; i < h.length; i++) {
            h[i] = scale * (rng.nextDouble() + 0.1);
        }
    }

    /**
     * Non-negative double SVD. The leading triplet of a non-negative matrix is
     * non-negative by Perron-Frobenius and is taken as it stands; every later one is
     * split into its positive and its negative part and the heavier of the two is
     * kept.
     */
    private static void nndsvd(double[] x, double[] w, double[] h, int m, int n, int k) {
        boolean wide = m < n;
        double[] a;
        int rows;
        int cols;
        if (wide) {
            a = new double[n * m];
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    a[i * n + j] = x[j * m + i];
                }
            }
            rows = n;
            cols = m;
        } else {
            a = x.clone();
            rows = m;
            cols = n;
        }
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decomposeInPlace(a, rows, cols);

        double[] u = new double[m * k];
        double[] v = new double[n * k];
        for (int r = 0; r < k; r++) {
            for (int i = 0; i < rows; i++) {
                double value = svd.U[r * rows + i];
                if (wide) {
                    v[r * n + i] = value;
                } else {
                    u[r * m + i] = value;
                }
            }
            for (int i = 0; i < cols; i++) {
                double value = svd.V[r * cols + i];
                if (wide) {
                    u[r * m + i] = value;
                } else {
                    v[r * n + i] = value;
                }
            }
        }

        double s0 = Math.sqrt(svd.sigma[0]);
        for (int i = 0; i < m; i++) {
            w[i] = s0 * Math.abs(u[i]);
        }
        for (int j = 0; j < n; j++) {
            h[j * k] = s0 * Math.abs(v[j]);
        }

        double[] up = new double[m];
        double[] un = new double[m];
        double[] vp = new double[n];
        double[] vn = new double[n];
        for (int r = 1; r < k; r++) {
            double nup = 0.0;
            double nun = 0.0;
            for (int i = 0; i < m; i++) {
                double value = u[r * m + i];
                up[i] = (value > 0.0) ? value : 0.0;
                un[i] = (value < 0.0) ? -value : 0.0;
                nup += up[i] * up[i];
                nun += un[i] * un[i];
            }
            double nvp = 0.0;
            double nvn = 0.0;
            for (int j = 0; j < n; j++) {
                double value = v[r * n + j];
                vp[j] = (value > 0.0) ? value : 0.0;
                vn[j] = (value < 0.0) ? -value : 0.0;
                nvp += vp[j] * vp[j];
                nvn += vn[j] * vn[j];
            }
            nup = Math.sqrt(nup);
            nun = Math.sqrt(nun);
            nvp = Math.sqrt(nvp);
            nvn = Math.sqrt(nvn);
            double termP = nup * nvp;
            double termN = nun * nvn;
            double[] uu;
            double[] vv;
            double nu;
            double nv;
            double term;
            if (termP >= termN) {
                uu = up;
                vv = vp;
                nu = nup;
                nv = nvp;
                term = termP;
            } else {
                uu = un;
                vv = vn;
                nu = nun;
                nv = nvn;
                term = termN;
            }
            double scale = Math.sqrt(svd.sigma[r] * term);
            if (nu > 0.0 && nv > 0.0 && scale > 0.0) {
                for (int i = 0; i < m; i++) {
                    w[r * m + i] = scale * uu[i] / nu;
                }
                for (int j = 0; j < n; j++) {
                    h[j * k + r] = scale * vv[j] / nv;
                }
            }
        }
    }

    /** Re-seeds any component that came out of the start with a zero factor. */
    private static void repairEmptyComponents(double[] x, double[] w, double[] h, int m, int n, int k) {
        double mean = 0.0;
        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }
        mean /= x.length;
        double scale = Math.sqrt(mean / k);
        if (!(scale > 0.0)) {
            return;
        }
        SplitMix64 rng = new SplitMix64(DEFAULT_SEED);
        for (int r = 0; r < k; r++) {
            double sw = 0.0;
            for (int i = 0; i < m; i++) {
                sw += w[r * m + i];
            }
            double sh = 0.0;
            for (int j = 0; j < n; j++) {
                sh += h[j * k + r];
            }
            if (sw > 0.0 && sh > 0.0) {
                continue;
            }
            for (int i = 0; i < m; i++) {
                w[r * m + i] = scale * (rng.nextDouble() + 0.1);
            }
            for (int j = 0; j < n; j++) {
                h[j * k + r] = scale * (rng.nextDouble() + 0.1);
            }
        }
    }

    // ----------------------------------------------------------------- misc --

    private static double normSquared(double[] a) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) {
            s += a[i] * a[i];
        }
        return s;
    }

    /** {@code ||X - W H||_F^2}, formed from the factors rather than accumulated. */
    private static double residualSquared(double[] x, double[] w, double[] h, int m, int n, int k) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int hb = j * k;
            int xb = j * m;
            for (int i = 0; i < m; i++) {
                double fit = 0.0;
                for (int r = 0; r < k; r++) {
                    fit += w[r * m + i] * h[hb + r];
                }
                double d = x[xb + i] - fit;
                sum += d * d;
            }
        }
        return sum;
    }

    private static double penaltyValue(double[] a, Penalty p) {
        if (!p.isActive()) {
            return 0.0;
        }
        double l1 = p.l1();
        double l2 = p.l2();
        double sum = 0.0;
        double squares = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i];
            squares += a[i] * a[i];
        }
        return l1 * sum + 0.5 * l2 * squares;
    }

    /**
     * Multiplies the factors back out.
     *
     * @param r a factorization
     * @return {@code W H}, {@code m x n} column-major
     */
    public static double[] reconstruct(Result r) {
        if (r == null) {
            throw new IllegalArgumentException("r must not be null");
        }
        double[] out = new double[r.m * r.n];
        Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, r.m, r.n, r.k, 1.0, r.W, 0, r.m, r.H, 0, r.k, 0.0,
                out, 0, r.m);
        return out;
    }

    /**
     * Relative reconstruction error of a factorization against data of the caller's
     * choosing, which need not be the data it was fitted to.
     *
     * @param X column-major packed, length {@code r.m * r.n}
     * @param r a factorization
     * @return {@code ||X - W H||_F / ||X||_F}, and {@code 0} for an all zero {@code X}
     */
    public static double reconstructionError(double[] X, Result r) {
        if (r == null) {
            throw new IllegalArgumentException("r must not be null");
        }
        if (X == null || X.length != (long) r.m * r.n) {
            throw new IllegalArgumentException("length != m*n");
        }
        double den = normSquared(X);
        if (!(den > 0.0)) {
            return 0.0;
        }
        return Math.sqrt(residualSquared(X, r.W, r.H, r.m, r.n, r.k) / den);
    }

    private static void check(double[] x, int m, int n, int k) {
        if (x == null) {
            throw new IllegalArgumentException("X must not be null");
        }
        if (m < 1 || n < 1) {
            throw new IllegalArgumentException("m and n must be >= 1 : " + m + "x" + n);
        }
        if (x.length != (long) m * n) {
            throw new IllegalArgumentException("length != m*n : " + x.length + " != " + m + "*" + n);
        }
        if (k < 1 || k > Math.min(m, n)) {
            throw new IllegalArgumentException("k must lie in [1, min(m, n)] : " + k);
        }
        for (int i = 0; i < x.length; i++) {
            if (!(x[i] >= 0.0) || Double.isInfinite(x[i])) {
                throw new IllegalArgumentException(
                        "X[" + (i % m) + ", " + (i / m) + "] is not finite and >= 0 : " + x[i]);
            }
        }
    }

    /**
     * Self-check on parts-based data of a known rank.
     *
     * @param args ignored
     */
    public static void main(String[] args) {
        int m = 400;
        int n = 200;
        int k = 12;
        long state = 20260821L;
        double[] w0 = new double[m * k];
        double[] h0 = new double[k * n];
        for (int i = 0; i < w0.length; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            double u = (state >>> 11) * 0x1.0p-53;
            w0[i] = (((i % m) % k) == (i / m)) ? 1.0 + u : ((u < 0.1) ? u * 0.2 : 0.0);
        }
        for (int i = 0; i < h0.length; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            double u = (state >>> 11) * 0x1.0p-53;
            h0[i] = (u < 0.4) ? u * 2.0 : 0.0;
        }
        double[] x = new double[m * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double v = 0.0;
                for (int r = 0; r < k; r++) {
                    v += w0[r * m + i] * h0[j * k + r];
                }
                x[j * m + i] = v;
            }
        }

        Result exact = new Nmf(1.0e-8, 2000, Init.NNDSVD, 0L, Penalty.none(), Penalty.none()).factor(x, m,
                n, k);
        Result plain = new Nmf().factor(x, m, n, k);

        double negative = 0.0;
        for (int i = 0; i < exact.W.length; i++) {
            negative = Math.min(negative, exact.W[i]);
        }
        for (int i = 0; i < exact.H.length; i++) {
            negative = Math.min(negative, exact.H[i]);
        }
        double worstNorm = 0.0;
        for (int r = 0; r < k; r++) {
            double s = 0.0;
            for (int i = 0; i < m; i++) {
                s += exact.W[r * m + i] * exact.W[r * m + i];
            }
            worstNorm = Math.max(worstNorm, Math.abs(Math.sqrt(s) - 1.0));
        }
        boolean descending = true;
        for (int r = 1; r < k; r++) {
            if (exact.componentEnergy[r] > exact.componentEnergy[r - 1]) {
                descending = false;
            }
        }
        double roundTrip = reconstructionError(x, exact);

        System.out.println("Nmf self-check, X = " + m + "x" + n + " of exact rank " + k);
        System.out.println("  NNDSVD, tol 1e-8 : " + exact);
        System.out.println("  defaults         : " + plain);
        System.out.printf("  most negative entry            : %.3e%n", negative);
        System.out.printf("  worst | ||W_r|| - 1 |          : %.3e%n", worstNorm);
        System.out.printf("  reconstruction error           : %.3e%n", roundTrip);
        System.out.println("  component energies descending  : " + descending);

        boolean ok = negative == 0.0 && worstNorm < 1.0e-12 && descending && roundTrip < 1.0e-8
                && exact.converged && exact.deadComponents == 0
                && Math.abs(roundTrip - exact.relativeError) <= 1.0e-15;
        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }
}
