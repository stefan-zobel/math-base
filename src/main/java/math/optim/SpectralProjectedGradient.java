package math.optim;

import java.util.Arrays;

import math.MathConsts;
import math.fun.DiffDMultiFunction;

/**
 * The spectral projected gradient method: minimization of a differentiable
 * function over a closed convex set, given only a {@link Projection} onto that
 * set. It is the answer to keeping a variance positive, a weight in
 * {@code [0, 1]} or a set of mixture proportions on the simplex, none of which
 * the unconstrained optimizers of this package can express.
 * <p>
 * Each iteration moves along {@code P(x - lambda*g) - x}, the projected
 * gradient direction, and searches on the segment from {@code x} to its end.
 * That segment stays inside the set because the set is convex and both of its
 * endpoints are feasible, so no projection is needed inside the line search.
 * The search stops when the projected gradient {@code P(x - g) - x} is small,
 * which is the first-order optimality condition of the constrained problem.
 * <p>
 * The method is of first order, so convergence is linear rather than
 * superlinear. What makes it competitive in practice is the Barzilai-Borwein
 * choice of {@code lambda} and a line search that tolerates an occasional
 * increase of the function value.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Barzilai-Borwein_method">Wikipedia
 * Barzilai-Borwein method</a>.
 *
 * @since 1.5.2
 */
public final class SpectralProjectedGradient {

    /** Tolerance on the projected gradient if none is given. */
    private static final double DEFAULT_KKT_TOLERANCE = Math.sqrt(MathConsts.MACH_EPS_DBL);
    /** Tolerance on the relative step if none is given. */
    private static final double DEFAULT_STEP_TOLERANCE = 1.0e-12;
    /** Constant part of the default iteration budget. */
    private static final int ITERATIONS_BASE = 1000;
    /** Part of the default iteration budget that grows with the dimension. */
    private static final int ITERATIONS_PER_DIMENSION = 100;
    /** How many past values the line search may look back over. */
    private static final int DEFAULT_MEMORY = 10;
    /** Fraction of the predicted decrease a step has to deliver. */
    private static final double ARMIJO = 1.0e-4;
    /** Below this step fraction the line search gives up. */
    private static final double MIN_ALPHA = 1.0e-30;
    /** Smallest spectral step, as in Birgin, Martinez and Raydan. */
    private static final double LAMBDA_MIN = 1.0e-30;
    /** Largest spectral step, as in Birgin, Martinez and Raydan. */
    private static final double LAMBDA_MAX = 1.0e30;

    private final double kktTolerance;
    private final double stepTolerance;
    private final int maxIterations;
    private final int memory;

    /**
     * Why the search stopped. Both {@link #KKT_TOLERANCE_REACHED} and
     * {@link #STEP_TOLERANCE_REACHED} are successful outcomes; the other three
     * mean the point returned is wherever the search happened to be.
     *
     * @since 1.5.2
     */
    public enum Status {

        /**
         * The projected gradient fell below the tolerance, which is the
         * first-order optimality condition of the constrained problem.
         */
        KKT_TOLERANCE_REACHED(true),

        /**
         * The step became too small to change the point meaningfully. The
         * projected gradient is still above its tolerance, so this is
         * convergence of the iterates rather than of the optimality measure.
         */
        STEP_TOLERANCE_REACHED(true),

        /** The iteration budget ran out. */
        TOO_MANY_ITERATIONS(false),

        /**
         * No step along the projected gradient direction reduced the function
         * value, or the direction was not one of descent, which a projection
         * onto a set that is not convex can cause.
         */
        LINE_SEARCH_FAILED(false),

        /** The gradient at the current point is not finite. */
        NOT_FINITE(false);

        private final boolean success;

        Status(boolean success) {
            this.success = success;
        }

        /**
         * Whether this outcome is one that the caller may use.
         *
         * @return {@code true} for the two tolerance outcomes
         */
        public boolean isSuccess() {
            return success;
        }
    }

    /**
     * The outcome of a minimization. The point is always feasible, whatever
     * the status, because every iterate is produced by the projection.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** The point the search settled on, always inside the set. */
        public final double[] point;

        /** The function value at {@link #point}. */
        public final double value;

        /**
         * The infinity norm of the projected gradient {@code P(x - g) - x} at
         * {@link #point}, which is zero exactly at a constrained stationary
         * point and is the quantity {@code kktTolerance} is compared against.
         */
        public final double projectedGradientNorm;

        /** How many iterations were performed. */
        public final int iterations;

        /** How often the function was evaluated. */
        public final int functionEvaluations;

        /** How often the gradient was evaluated. */
        public final int gradientEvaluations;

        /** Why the search stopped. */
        public final Status status;

        /** Shorthand for {@code status.isSuccess()}. */
        public final boolean converged;

        Result(double[] point, double value, double projectedGradientNorm, int iterations, int functionEvaluations,
                int gradientEvaluations, Status status) {
            this.point = point;
            this.value = value;
            this.projectedGradientNorm = projectedGradientNorm;
            this.iterations = iterations;
            this.functionEvaluations = functionEvaluations;
            this.gradientEvaluations = gradientEvaluations;
            this.status = status;
            this.converged = status.isSuccess();
        }
    }

    /**
     * Creates a minimizer with default settings.
     */
    public SpectralProjectedGradient() {
        this(DEFAULT_KKT_TOLERANCE, DEFAULT_STEP_TOLERANCE, 0, DEFAULT_MEMORY);
    }

    /**
     * Creates a minimizer.
     *
     * @param kktTolerance
     *            the search stops once the infinity norm of the projected
     *            gradient falls to this value; not negative
     * @param stepTolerance
     *            the search stops once the step falls to this value relative
     *            to the size of the current point; not negative
     * @param maxIterations
     *            maximal number of iterations; {@code 0} or less selects
     *            {@code 1000 + 100 * n}
     * @param memory
     *            how many past function values the line search compares
     *            against, at least one; {@code 1} makes the search monotone
     */
    public SpectralProjectedGradient(double kktTolerance, double stepTolerance, int maxIterations, int memory) {
        if (!(kktTolerance >= 0.0)) {
            throw new IllegalArgumentException("kktTolerance must be >= 0.0 : " + kktTolerance);
        }
        if (!(stepTolerance >= 0.0)) {
            throw new IllegalArgumentException("stepTolerance must be >= 0.0 : " + stepTolerance);
        }
        if (memory < 1) {
            throw new IllegalArgumentException("memory must be >= 1 : " + memory);
        }
        this.kktTolerance = kktTolerance;
        this.stepTolerance = stepTolerance;
        this.maxIterations = maxIterations;
        this.memory = memory;
    }

    /**
     * Minimizes {@code f} over the box {@code lower <= x <= upper}, the most
     * common case. Equivalent to passing {@link Projection#box(double[],
     * double[])}.
     *
     * @param f
     *            the function to minimize
     * @param start
     *            the starting point, not modified, projected into the box if
     *            it lies outside
     * @param lower
     *            the lower bounds, possibly infinite
     * @param upper
     *            the upper bounds, possibly infinite
     * @return the point the search has settled on and how it got there
     */
    public Result minimize(DiffDMultiFunction f, double[] start, double[] lower, double[] upper) {
        return minimize(f, start, Projection.box(lower, upper));
    }

    /**
     * Minimizes {@code f} over the set that {@code projection} projects onto.
     *
     * @param f
     *            the function to minimize
     * @param start
     *            the starting point, not modified, projected into the set if
     *            it lies outside
     * @param projection
     *            the projection onto the feasible set, which has to be closed
     *            and convex
     * @return the point the search has settled on and how it got there
     */
    public Result minimize(DiffDMultiFunction f, double[] start, Projection projection) {
        if (f == null) {
            throw new IllegalArgumentException("f is null");
        }
        if (start == null) {
            throw new IllegalArgumentException("start is null");
        }
        if (projection == null) {
            throw new IllegalArgumentException("projection is null");
        }
        int n = start.length;
        if (n < 1) {
            throw new IllegalArgumentException("start must have at least one component");
        }
        for (int i = 0; i < n; i++) {
            if (!isFinite(start[i])) {
                throw new IllegalArgumentException("start[" + i + "] is not finite : " + start[i]);
            }
        }

        double[] x = start.clone();
        projection.projectInto(x);
        for (int i = 0; i < n; i++) {
            if (!isFinite(x[i])) {
                throw new IllegalArgumentException(
                        "the projected starting point is not finite at index " + i + " : " + x[i]);
            }
        }

        double[] gradient = new double[n];
        double[] previousGradient = new double[n];
        double[] previousStep = new double[n];
        double[] direction = new double[n];
        double[] candidate = new double[n];
        double[] scratch = new double[n];

        double value = f.apply(x);
        int functionEvaluations = 1;
        if (!isFinite(value)) {
            throw new IllegalArgumentException("f is not finite at the projected starting point : " + value);
        }
        int gradientEvaluations = 0;
        int iterations = 0;
        int budget = (maxIterations > 0) ? maxIterations : defaultBudget(n);
        double lambda = 1.0;
        double kkt = Double.NaN;
        boolean stepIsSmall = false;
        Status status;

        // the last few accepted values, so that the line search can measure
        // against the worst of them rather than against the latest. Unfilled
        // slots never win the maximum
        double[] history = new double[memory];
        Arrays.fill(history, Double.NEGATIVE_INFINITY);
        history[0] = value;
        int nextSlot = 1 % memory;

        for (;;) {
            f.derivativeAt(x, gradient);
            gradientEvaluations++;
            if (!allFinite(gradient)) {
                status = Status.NOT_FINITE;
                break;
            }
            kkt = projectedGradientNorm(projection, x, gradient, scratch);
            if (kkt <= kktTolerance) {
                status = Status.KKT_TOLERANCE_REACHED;
                break;
            }
            if (stepIsSmall) {
                status = Status.STEP_TOLERANCE_REACHED;
                break;
            }
            if (iterations >= budget) {
                status = Status.TOO_MANY_ITERATIONS;
                break;
            }

            // The Barzilai-Borwein step: the one scalar multiple of the
            // identity that best satisfies the secant equation over the step
            // just taken. It carries the curvature that a plain gradient
            // method throws away and costs two dot products. A curvature that
            // is not positive says nothing, so the longest allowed step is
            // taken instead
            if (iterations == 0) {
                lambda = clamp(1.0 / kkt, LAMBDA_MIN, LAMBDA_MAX);
            } else {
                double stepDotYield = 0.0;
                double stepDotStep = 0.0;
                for (int i = 0; i < n; i++) {
                    double yield = gradient[i] - previousGradient[i];
                    stepDotYield += previousStep[i] * yield;
                    stepDotStep += previousStep[i] * previousStep[i];
                }
                lambda = (stepDotYield > 0.0) ? clamp(stepDotStep / stepDotYield, LAMBDA_MIN, LAMBDA_MAX)
                        : LAMBDA_MAX;
            }

            // the projected gradient direction, and the slope along it. For a
            // projection onto a convex set the slope is at most
            // -||direction||^2 / lambda, so anything else means the set the
            // projection describes is not convex
            for (int i = 0; i < n; i++) {
                candidate[i] = x[i] - lambda * gradient[i];
            }
            projection.projectInto(candidate);
            double slope = 0.0;
            for (int i = 0; i < n; i++) {
                double d = candidate[i] - x[i];
                direction[i] = d;
                slope += gradient[i] * d;
            }
            if (!(slope < 0.0)) {
                status = Status.LINE_SEARCH_FAILED;
                break;
            }

            // the first trial is the projected point itself rather than
            // x + 1.0 * direction, which is the same point only up to
            // rounding. Taken literally, a full step then lands exactly on
            // the boundary of the set, which is what makes a bound come out
            // as the bound and not as one ulp beside it
            // the worst of the last few values, not the latest one. Letting
            // the value rise for a while is what keeps the spectral step
            // usable: it is a step length, not a descent direction, and
            // forcing a decrease at every single iteration throws most of it
            // away again
            double reference = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < history.length; i++) {
                if (history[i] > reference) {
                    reference = history[i];
                }
            }

            double alpha = 1.0;
            double trialValue = f.apply(candidate);
            functionEvaluations++;
            // negated, so that a trial value of NaN backtracks
            boolean accepted = trialValue <= reference + ARMIJO * slope;
            while (!accepted) {
                // plain halving. Birgin, Martinez and Raydan backtrack by the
                // minimizer of the quadratic through f(x), the slope and the
                // trial value, safeguarded into [0.1, 0.9] of the current
                // step. Measured here it does not pay: on Rosenbrock in a box
                // it loses four to five orders of magnitude in the projected
                // gradient at equal budget for every dimension from four
                // upwards, and on a regularized logistic loss and a two-term
                // exponential fit it is a tie. The safeguard band is not the
                // reason -- the interpolated step lands below half the current
                // one almost every time, so the upper end of the band never
                // binds
                alpha *= 0.5;
                if (alpha < MIN_ALPHA) {
                    break;
                }
                for (int i = 0; i < n; i++) {
                    candidate[i] = x[i] + alpha * direction[i];
                }
                trialValue = f.apply(candidate);
                functionEvaluations++;
                accepted = trialValue <= reference + ARMIJO * alpha * slope;
            }
            if (!accepted) {
                status = Status.LINE_SEARCH_FAILED;
                break;
            }
            history[nextSlot] = trialValue;
            nextSlot++;
            if (nextSlot == history.length) {
                nextSlot = 0;
            }

            double largestStep = 0.0;
            double size = 0.0;
            for (int i = 0; i < n; i++) {
                double d = candidate[i] - x[i];
                previousStep[i] = d;
                double magnitude = Math.abs(d);
                if (magnitude > largestStep) {
                    largestStep = magnitude;
                }
                double a = Math.abs(x[i]);
                if (a > size) {
                    size = a;
                }
            }
            stepIsSmall = largestStep <= stepTolerance * (1.0 + size);

            System.arraycopy(gradient, 0, previousGradient, 0, n);
            System.arraycopy(candidate, 0, x, 0, n);
            value = trialValue;
            iterations++;
        }

        return new Result(x, value, kkt, iterations, functionEvaluations, gradientEvaluations, status);
    }

    /**
     * The infinity norm of {@code P(x - g) - x}, the first-order optimality
     * measure of the constrained problem. It deliberately uses a unit step
     * rather than the current spectral one: with the current step the measure
     * would shrink whenever the step does, and the search would stop early.
     */
    private static double projectedGradientNorm(Projection projection, double[] x, double[] gradient,
            double[] scratch) {
        for (int i = 0; i < x.length; i++) {
            scratch[i] = x[i] - gradient[i];
        }
        projection.projectInto(scratch);
        double norm = 0.0;
        for (int i = 0; i < x.length; i++) {
            double d = Math.abs(scratch[i] - x[i]);
            if (d > norm) {
                norm = d;
            }
        }
        return norm;
    }

    /** Confines {@code v} to {@code [lo, hi]}, sending {@code NaN} to {@code lo}. */
    private static double clamp(double v, double lo, double hi) {
        if (!(v > lo)) {
            return lo;
        }
        return (v < hi) ? v : hi;
    }

    private static int defaultBudget(int n) {
        long budget = ITERATIONS_BASE + (long) ITERATIONS_PER_DIMENSION * n;
        return (budget > Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) budget;
    }

    private static boolean allFinite(double[] v) {
        for (int i = 0; i < v.length; i++) {
            if (!isFinite(v[i])) {
                return false;
            }
        }
        return true;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
