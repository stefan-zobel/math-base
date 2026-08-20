package math.linalg;

import java.util.Arrays;

import math.MathConsts;

/**
 * Active set core of {@link BoundedLeastSquares}, on flat column-major arrays.
 * <p>
 * Solves {@code min ||Ax - b||_2} subject to {@code l <= x <= u} by the method
 * of Stark and Parker (1995), which is Lawson and Hanson's non-negative least
 * squares generalized to two-sided bounds. Every variable is either held at one
 * of its bounds or free; the free block is solved as an unconstrained least
 * squares problem, and a variable that leaves its box on the way is put back on
 * the bound it crossed and frozen there.
 */
final class Bvls {

    /** the variable sits exactly on its lower bound */
    static final byte AT_LOWER = 0;

    /** the variable is free and the gradient vanishes in its direction */
    static final byte FREE = 1;

    /** the variable sits exactly on its upper bound */
    static final byte AT_UPPER = 2;

    /**
     * Absolute floor of the optimality test, in units of
     * {@code u * ||A||_F * ||r||}. The gradient is {@code w = A'r}, so
     * {@code |w_j| <= ||a_j|| * ||r||} and this is the scale on which "zero"
     * has to be read; a fixed absolute tolerance would not survive a rescaling
     * of the problem.
     */
    private static final double TOLERANCE_FACTOR = 100.0;

    /**
     * How far above the free set's own gradient a violation has to stand
     * before the variable is freed. The free variables carry {@code w_j == 0}
     * by construction, so whatever is left there is the noise of the inner
     * solve on this particular problem -- a better floor than any constant,
     * because it grows with the conditioning instead of being guessed for it.
     */
    private static final double FREE_FLOOR_FACTOR = 10.0;

    /** Outcome of {@link Bvls#solve}. */
    static final class Fit {

        final double[] x;
        final byte[] state;
        final int iterations;
        final boolean converged;
        final boolean rankDeficient;

        Fit(double[] x, byte[] state, int iterations, boolean converged, boolean rankDeficient) {
            this.x = x;
            this.state = state;
            this.iterations = iterations;
            this.converged = converged;
            this.rankDeficient = rankDeficient;
        }
    }

    private final double[] a;
    private final int m;
    private final int n;
    private final double[] b;
    private final double[] lower;
    private final double[] upper;
    private final double[] x;
    private final byte[] state;
    private final double[] residual;
    private final double[] gradient;
    private final double normA;
    private boolean rankDeficient;

    /**
     * Iteration cap, a backstop rather than a working limit. The measured
     * worst case over conditioned families up to {@code kappa = 1e15}, and
     * over degenerate designs with repeated columns and ties, is
     * {@code 2.7 * min(m, n)} outer iterations, so this leaves a factor of
     * three and a half.
     *
     * @param m
     *            number of rows
     * @param n
     *            number of columns
     * @return the default cap on outer iterations
     */
    static int defaultMaxIterations(int m, int n) {
        return Math.max(100, 10 * Math.min(m, n));
    }

    /**
     * Solves the bounded problem.
     *
     * @param a
     *            the {@code m x n} design, column-major, left unchanged
     * @param m
     *            number of rows
     * @param n
     *            number of columns
     * @param b
     *            the right hand side, length {@code m}, left unchanged
     * @param lower
     *            lower bounds, length {@code n}, may contain
     *            {@code -Infinity}
     * @param upper
     *            upper bounds, length {@code n}, may contain
     *            {@code +Infinity}
     * @param maxIterations
     *            cap on outer iterations
     * @return the solution together with which bounds ended up binding
     */
    static Fit solve(double[] a, int m, int n, double[] b, double[] lower, double[] upper, int maxIterations) {
        return new Bvls(a, m, n, b, lower, upper).run(maxIterations);
    }

    private Bvls(double[] a, int m, int n, double[] b, double[] lower, double[] upper) {
        // A and b are scaled by one common power of two, which leaves the
        // solution and the bounds untouched -- ||cAx - cb|| = |c| ||Ax - b|| --
        // and is exact, since only the exponent moves. Without it the whole
        // optimality test underflows away: at entries of 1e-200 both the
        // gradient A'r and the tolerance eps*||A||_F*||r|| are products of two
        // such numbers, so both collapse to zero and every comparison between
        // them holds vacuously. Measured: the same problem scaled by 1e-200
        // came back with a different active set.
        double magnitude = Math.max(maxAbs(a), maxAbs(b));
        double scale = (magnitude > 0.0) ? Math.scalb(1.0, -Math.getExponent(magnitude)) : 1.0;
        if (scale != 1.0) {
            a = scaled(a, scale);
            b = scaled(b, scale);
        }
        this.a = a;
        this.m = m;
        this.n = n;
        this.b = b;
        this.lower = lower;
        this.upper = upper;
        this.x = new double[n];
        this.state = new byte[n];
        this.residual = new double[m];
        this.gradient = new double[n];
        this.normA = frobeniusNorm(a);
    }

    private Fit run(int maxIterations) {
        initialize();
        // A coefficient that is unbounded on both sides has no bound to be held
        // at and starts free, so unlike in the non-negative case the free set
        // can be non-empty before the first violation is looked for. Without
        // this the loop below would find nothing to release, report
        // convergence and hand back the starting point: with every bound
        // infinite it returned all zeros instead of the unconstrained fit.
        refineFreeSet();
        boolean[] ineligible = new boolean[n];
        computeResidual();
        double previousObjective = sumOfSquares(residual);
        int iterations = 0;
        boolean converged = false;

        while (iterations < maxIterations) {
            iterations++;
            computeResidual();
            computeGradient();
            int candidate = mostViolating(ineligible);
            if (candidate < 0) {
                converged = true;
                break;
            }
            state[candidate] = FREE;
            refineFreeSet();
            computeResidual();
            double objective = sumOfSquares(residual);
            // Termination rests on the objective rather than on the active
            // set. Lawson and Hanson's rule -- clear the flags whenever the
            // freed variable stays free -- is not enough in floating point: at
            // a condition number of 1e12 a pair of variables was measured
            // taking turns indefinitely, each of them staying free while the
            // residual no longer moved. Clearing only after a strict decrease
            // bounds the whole search, because a pass that achieves nothing
            // can free each variable at most once before it runs out of
            // candidates.
            if (objective < previousObjective) {
                previousObjective = objective;
                Arrays.fill(ineligible, false);
            } else {
                ineligible[candidate] = true;
            }
        }
        return new Fit(x, state, iterations, converged, rankDeficient);
    }

    /**
     * Every variable starts on the bound nearer zero, which for
     * {@code l = 0, u = +inf} is exactly Lawson and Hanson's non-negative
     * least squares start. A variable that is bounded on neither side has no
     * bound to sit on and starts free.
     */
    private void initialize() {
        for (int j = 0; j < n; j++) {
            boolean lowerInfinite = Double.isInfinite(lower[j]);
            boolean upperInfinite = Double.isInfinite(upper[j]);
            if (lowerInfinite && upperInfinite) {
                state[j] = FREE;
                x[j] = 0.0;
            } else if (upperInfinite || (!lowerInfinite && Math.abs(lower[j]) <= Math.abs(upper[j]))) {
                state[j] = AT_LOWER;
                x[j] = lower[j];
            } else {
                state[j] = AT_UPPER;
                x[j] = upper[j];
            }
        }
    }

    /**
     * The bound variable whose gradient argues hardest for letting it move
     * inwards, or {@code -1} when none of them does, which is the optimality
     * condition of the bounded problem.
     */
    private int mostViolating(boolean[] ineligible) {
        int free = 0;
        double freeFloor = 0.0;
        for (int j = 0; j < n; j++) {
            if (state[j] == FREE) {
                free++;
                freeFloor = Math.max(freeFloor, Math.abs(gradient[j]));
            }
        }
        // the free block is factorized as an m x k matrix, so it cannot hold
        // more columns than there are rows; at k == m the residual is already
        // zero and there is nothing left to gain anyway
        if (free >= m) {
            return -1;
        }
        double tolerance = Math.max(TOLERANCE_FACTOR * MathConsts.MACH_EPS_DBL * normA * frobeniusNorm(residual),
                FREE_FLOOR_FACTOR * freeFloor);
        int best = -1;
        double bestViolation = tolerance;
        for (int j = 0; j < n; j++) {
            if (state[j] == FREE || ineligible[j] || lower[j] == upper[j]) {
                continue;
            }
            double violation = (state[j] == AT_LOWER) ? gradient[j] : -gradient[j];
            if (violation > bestViolation) {
                bestViolation = violation;
                best = j;
            }
        }
        return best;
    }

    /**
     * Solves the unconstrained problem over the free columns and, as long as
     * that answer leaves the box, walks as far towards it as feasibility
     * allows and freezes whatever reached a bound.
     */
    private void refineFreeSet() {
        while (true) {
            int k = 0;
            for (int j = 0; j < n; j++) {
                if (state[j] == FREE) {
                    k++;
                }
            }
            if (k == 0) {
                return;
            }
            int[] index = new int[k];
            int p = 0;
            for (int j = 0; j < n; j++) {
                if (state[j] == FREE) {
                    index[p++] = j;
                }
            }
            // column-major pays off here: a free column is contiguous
            double[] free = new double[m * k];
            for (int c = 0; c < k; c++) {
                System.arraycopy(a, index[c] * m, free, c * m, m);
            }
            double[] rhs = new double[m];
            System.arraycopy(b, 0, rhs, 0, m);
            for (int j = 0; j < n; j++) {
                if (state[j] != FREE && x[j] != 0.0) {
                    int column = j * m;
                    for (int i = 0; i < m; i++) {
                        rhs[i] -= a[column + i] * x[j];
                    }
                }
            }
            // decomposeInPlace consumes free, which is a scratch copy anyway
            FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decomposeInPlace(free, m, k);
            if (!svd.converged) {
                throw new RuntimeException("the singular value decomposition of the free block did not converge");
            }
            double singularTolerance = svd.sigma[0] * Math.max(m, k) * MathConsts.MACH_EPS_DBL;
            for (int i = 0; i < k; i++) {
                if (svd.sigma[i] <= singularTolerance) {
                    rankDeficient = true;
                    break;
                }
            }
            double[] z = SvdLeastSquares.solveTruncated(svd, rhs, singularTolerance);

            boolean inside = true;
            for (int c = 0; c < k; c++) {
                if (z[c] < lower[index[c]] || z[c] > upper[index[c]]) {
                    inside = false;
                    break;
                }
            }
            if (inside) {
                for (int c = 0; c < k; c++) {
                    x[index[c]] = z[c];
                }
                return;
            }
            double alpha = longestFeasibleStep(index, z);
            for (int c = 0; c < k; c++) {
                int j = index[c];
                x[j] = x[j] + alpha * (z[c] - x[j]);
            }
            if (!freezeWhatReachedABound(index) && !freezeTheWorstOffender(index, z)) {
                return;
            }
        }
    }

    /** The largest fraction of the way to {@code z} that stays inside the box. */
    private double longestFeasibleStep(int[] index, double[] z) {
        double alpha = 1.0;
        for (int c = 0; c < index.length; c++) {
            int j = index[c];
            double bound;
            if (z[c] < lower[j]) {
                bound = lower[j];
            } else if (z[c] > upper[j]) {
                bound = upper[j];
            } else {
                continue;
            }
            double direction = z[c] - x[j];
            if (direction == 0.0) {
                // already on the bound it is being pushed through
                return 0.0;
            }
            double step = (bound - x[j]) / direction;
            if (step < alpha) {
                alpha = step;
            }
        }
        return (alpha > 0.0) ? alpha : 0.0;
    }

    /**
     * Freezes every free variable that arrived at a bound, assigning the bound
     * itself rather than the value the step happened to land on. That
     * assignment is what makes {@code activeSet} exact.
     */
    private boolean freezeWhatReachedABound(int[] index) {
        boolean any = false;
        for (int c = 0; c < index.length; c++) {
            int j = index[c];
            if (x[j] <= lower[j]) {
                x[j] = lower[j];
                state[j] = AT_LOWER;
                any = true;
            } else if (x[j] >= upper[j]) {
                x[j] = upper[j];
                state[j] = AT_UPPER;
                any = true;
            }
        }
        return any;
    }

    /**
     * Fallback for the case where rounding leaves the step short of every
     * bound although the target was outside: freezing the variable that
     * overshot furthest guarantees the free set shrinks, so the loop above
     * cannot spin.
     */
    private boolean freezeTheWorstOffender(int[] index, double[] z) {
        int worst = -1;
        int worstColumn = -1;
        double worstExcess = 0.0;
        for (int c = 0; c < index.length; c++) {
            int j = index[c];
            double excess = Math.max(lower[j] - z[c], z[c] - upper[j]);
            if (excess > worstExcess) {
                worstExcess = excess;
                worst = j;
                worstColumn = c;
            }
        }
        if (worst < 0) {
            return false;
        }
        if (z[worstColumn] < lower[worst]) {
            x[worst] = lower[worst];
            state[worst] = AT_LOWER;
        } else {
            x[worst] = upper[worst];
            state[worst] = AT_UPPER;
        }
        return true;
    }

    private void computeResidual() {
        System.arraycopy(b, 0, residual, 0, m);
        for (int j = 0; j < n; j++) {
            double xj = x[j];
            if (xj == 0.0) {
                continue;
            }
            int column = j * m;
            for (int i = 0; i < m; i++) {
                residual[i] -= a[column + i] * xj;
            }
        }
    }

    /** {@code w = A'r}, the negative gradient of {@code 0.5 ||Ax - b||^2}. */
    private void computeGradient() {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            int column = j * m;
            for (int i = 0; i < m; i++) {
                sum += a[column + i] * residual[i];
            }
            gradient[j] = sum;
        }
    }

    private static double sumOfSquares(double[] v) {
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            sum += v[i] * v[i];
        }
        return sum;
    }

    private static double maxAbs(double[] v) {
        double max = 0.0;
        for (int i = 0; i < v.length; i++) {
            double abs = Math.abs(v[i]);
            if (abs > max) {
                max = abs;
            }
        }
        return max;
    }

    private static double[] scaled(double[] v, double factor) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            out[i] = v[i] * factor;
        }
        return out;
    }

    /** Scaled so that a design of huge or tiny entries does not overflow here. */
    private static double frobeniusNorm(double[] v) {
        double max = 0.0;
        for (int i = 0; i < v.length; i++) {
            double abs = Math.abs(v[i]);
            if (abs > max) {
                max = abs;
            }
        }
        if (max == 0.0 || !(max < Double.POSITIVE_INFINITY)) {
            return max;
        }
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            double scaled = v[i] / max;
            sum += scaled * scaled;
        }
        return max * Math.sqrt(sum);
    }
}
