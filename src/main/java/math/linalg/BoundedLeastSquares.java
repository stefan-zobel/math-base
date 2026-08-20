package math.linalg;

import java.util.Arrays;

/**
 * Least squares under hard bounds: {@code min ||Ax - b||_2} subject to
 * {@code l <= x <= u}, with {@code l = 0} and {@code u = +inf} giving
 * non-negative least squares.
 * <p>
 * Unlike {@link Ridge} and {@link Lasso}, which push coefficients towards zero
 * with a penalty, this is a hard restriction: a coefficient either satisfies it
 * or the solution is rejected. The active set method behind it is finite -- it
 * arrives at the exact minimizer after a bounded number of changes to the set
 * of binding constraints, rather than converging to it -- and a coefficient
 * held at a bound carries that bound exactly, bit for bit.
 * <p>
 * See <a href="https://en.wikipedia.org/wiki/Non-negative_least_squares">
 * Non-negative least squares</a>; the two-sided form is due to Stark and Parker
 * (1995), the non-negative one to Lawson and Hanson (1974).
 *
 * @since 1.5.2
 */
public final class BoundedLeastSquares {

    /**
     * Where a coefficient ended up, which is what distinguishes a binding
     * constraint from a value that merely came out small.
     *
     * @since 1.5.2
     */
    public enum Bound {

        /** The coefficient is held at its lower bound. */
        AT_LOWER,

        /** The coefficient is interior and unconstrained at the solution. */
        FREE,

        /** The coefficient is held at its upper bound. */
        AT_UPPER
    }

    /**
     * The outcome of a bounded fit. The solution is always feasible, whatever
     * the value of {@link #converged}.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** The coefficients, length {@code A.numColumns()}. */
        public final double[] solution;

        /** {@code b - A x}, length {@code A.numRows()}. */
        public final double[] residuals;

        /** The Euclidean norm of {@link #residuals}. */
        public final double residualNorm;

        /**
         * Which constraint binds for each coefficient. An entry of
         * {@link Bound#AT_LOWER} means the coefficient equals its lower bound
         * exactly, not approximately.
         */
        public final Bound[] activeSet;

        /** How many outer iterations were performed. */
        public final int iterations;

        /**
         * Whether the optimality condition was met. {@code false} means the
         * iteration cap ran out, and nothing else.
         */
        public final boolean converged;

        /**
         * Whether the free block was rank deficient at some point, in which
         * case the minimizer is not unique and this is one of several. The
         * solution reported is the one of smallest norm over the free set.
         */
        public final boolean rankDeficient;

        /** How many coefficients came back at a bound rather than free. */
        public final int atBound;

        Result(double[] solution, double[] residuals, double residualNorm, Bound[] activeSet, int iterations,
                boolean converged, boolean rankDeficient) {
            this.solution = solution;
            this.residuals = residuals;
            this.residualNorm = residualNorm;
            this.activeSet = activeSet;
            this.iterations = iterations;
            this.converged = converged;
            this.rankDeficient = rankDeficient;
            int bound = 0;
            for (int j = 0; j < activeSet.length; j++) {
                if (activeSet[j] != Bound.FREE) {
                    bound++;
                }
            }
            this.atBound = bound;
        }

        @Override
        public String toString() {
            return "BoundedLeastSquares.Result[||r|| = " + residualNorm + ", " + atBound + " of " + activeSet.length
                    + " at a bound, " + iterations + " iterations" + (converged ? "" : ", NOT converged")
                    + (rankDeficient ? ", rank deficient" : "") + "]";
        }
    }

    /**
     * Non-negative least squares, {@code min ||Ax - b||} subject to
     * {@code x >= 0}.
     * <p>
     * The restriction is sparsifying on its own: the non-negative orthant has
     * corners, so coefficients land at exactly zero without any penalty and
     * without a tuning parameter to choose.
     *
     * @param a
     *            the {@code m x n} design matrix, left unchanged
     * @param b
     *            the right hand side, {@code m x 1}, left unchanged
     * @return the fit
     * @throws IllegalArgumentException
     *             if the shapes disagree or an entry is not finite
     */
    public static Result nonNegative(DMatrix a, DMatrix b) {
        int n = (a == null) ? 0 : a.numColumns();
        double[] lower = new double[n];
        double[] upper = new double[n];
        Arrays.fill(upper, Double.POSITIVE_INFINITY);
        return bounded(a, b, lower, upper);
    }

    /**
     * Bounded variable least squares, {@code min ||Ax - b||} subject to
     * {@code l <= x <= u}.
     * <p>
     * Infinite bounds are allowed on either side. {@code lower[j] ==
     * upper[j]} pins that coefficient to the common value.
     *
     * @param a
     *            the {@code m x n} design matrix, left unchanged
     * @param b
     *            the right hand side, {@code m x 1}, left unchanged
     * @param lower
     *            lower bounds, length {@code n}, left unchanged
     * @param upper
     *            upper bounds, length {@code n}, left unchanged
     * @return the fit
     * @throws IllegalArgumentException
     *             if the shapes disagree, an entry of {@code a} or {@code b} is
     *             not finite, a bound is {@code NaN}, some
     *             {@code lower[j] > upper[j]}, or more coefficients are
     *             unbounded on both sides than {@code a} has rows
     */
    public static Result bounded(DMatrix a, DMatrix b, double[] lower, double[] upper) {
        check(a, b, lower, upper);
        int m = a.numRows();
        int n = a.numColumns();
        Bvls.Fit fit = Bvls.solve(a.getArrayUnsafe(), m, n, b.getArrayUnsafe(), lower, upper,
                Bvls.defaultMaxIterations(m, n));
        double[] residuals = new double[m];
        System.arraycopy(b.getArrayUnsafe(), 0, residuals, 0, m);
        double[] design = a.getArrayUnsafe();
        for (int j = 0; j < n; j++) {
            double xj = fit.x[j];
            if (xj == 0.0) {
                continue;
            }
            int column = j * m;
            for (int i = 0; i < m; i++) {
                residuals[i] -= design[column + i] * xj;
            }
        }
        Bound[] activeSet = new Bound[n];
        for (int j = 0; j < n; j++) {
            if (fit.state[j] == Bvls.AT_LOWER) {
                activeSet[j] = Bound.AT_LOWER;
            } else if (fit.state[j] == Bvls.AT_UPPER) {
                activeSet[j] = Bound.AT_UPPER;
            } else {
                activeSet[j] = Bound.FREE;
            }
        }
        return new Result(fit.x, residuals, euclidean(residuals), activeSet, fit.iterations, fit.converged,
                fit.rankDeficient);
    }

    private static void check(DMatrix a, DMatrix b, double[] lower, double[] upper) {
        if (a == null || b == null || lower == null || upper == null) {
            throw new IllegalArgumentException("null argument");
        }
        if (a.numRows() != b.numRows()) {
            throw new IllegalArgumentException(
                    "A.numRows != b.numRows : " + a.numRows() + " != " + b.numRows());
        }
        if (b.numColumns() != 1) {
            throw new IllegalArgumentException("b must have exactly one column, has " + b.numColumns());
        }
        int m = a.numRows();
        int n = a.numColumns();
        if (m < 1 || n < 1) {
            throw new IllegalArgumentException("A must not be empty : " + m + " x " + n);
        }
        if (lower.length != n || upper.length != n) {
            throw new IllegalArgumentException("bounds must have length " + n + " : lower " + lower.length
                    + ", upper " + upper.length);
        }
        checkFinite(a.getArrayUnsafe(), m, "A");
        checkFinite(b.getArrayUnsafe(), m, "b");
        int unbounded = 0;
        for (int j = 0; j < n; j++) {
            if (Double.isNaN(lower[j]) || Double.isNaN(upper[j])) {
                throw new IllegalArgumentException("bound " + j + " is NaN");
            }
            if (lower[j] > upper[j]) {
                throw new IllegalArgumentException(
                        "lower > upper at " + j + " : " + lower[j] + " > " + upper[j]);
            }
            if (lower[j] == Double.POSITIVE_INFINITY || upper[j] == Double.NEGATIVE_INFINITY) {
                throw new IllegalArgumentException("the box is empty at " + j + " : [" + lower[j] + ", "
                        + upper[j] + "]");
            }
            if (Double.isInfinite(lower[j]) && Double.isInfinite(upper[j])) {
                unbounded++;
            }
        }
        // such a coefficient has no bound to be held at, so it is free from the
        // start, and the free block is factorized as an m x k matrix
        if (unbounded > m) {
            throw new IllegalArgumentException("more coefficients are unbounded on both sides than A has rows : "
                    + unbounded + " > " + m);
        }
    }

    private static void checkFinite(double[] values, int rows, String name) {
        for (int i = 0; i < values.length; i++) {
            if (!isFinite(values[i])) {
                throw new IllegalArgumentException(
                        name + "[" + (i % rows) + ", " + (i / rows) + "] is not finite : " + values[i]);
            }
        }
    }

    private static boolean isFinite(double d) {
        return !Double.isNaN(d) && !Double.isInfinite(d);
    }

    private static double euclidean(double[] v) {
        double max = 0.0;
        for (int i = 0; i < v.length; i++) {
            double abs = Math.abs(v[i]);
            if (abs > max) {
                max = abs;
            }
        }
        if (max == 0.0) {
            return 0.0;
        }
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            double scaled = v[i] / max;
            sum += scaled * scaled;
        }
        return max * Math.sqrt(sum);
    }

    private BoundedLeastSquares() {
        throw new AssertionError();
    }
}
