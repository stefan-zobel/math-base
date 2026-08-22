package math.solve;

import math.MathConsts;
import math.fun.DVectorFunction;
import math.fun.DiffDVectorFunction;
import math.fun.NumericalDiffDVectorFunction;
import math.lapack.Dgesv;

/**
 * Solves a system of {@code n} nonlinear equations in {@code n} unknowns,
 * {@code F(x) = 0}, by Newton's method or by Broyden's, each of them globalized
 * by a backtracking line search on the merit function
 * {@code f(x) = 0.5 * F(x)'F(x)}.
 * <p>
 * This is the nonlinear counterpart of {@link LinearEquationsSolver} and the
 * multidimensional counterpart of {@link RootFinder}. A pure Newton step
 * {@code J(x) s = -F(x)} converges quadratically once it is close enough and
 * diverges cheerfully when it is not; the line search is what makes it usable
 * from an arbitrary starting point, at the price that the iteration can settle
 * into a local minimum of {@code f} where {@code F} is not zero. That outcome
 * is a fact about the system rather than a failure of the search, and it is
 * reported as {@link Status#SPURIOUS_MINIMUM} rather than hidden behind a
 * convergence flag.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Newton%27s_method#Systems_of_equations">
 *      Newton's method for systems of equations</a>
 * @see <a href="https://en.wikipedia.org/wiki/Broyden%27s_method">Broyden's
 *      method</a>
 * @since 1.5.2
 */
public final class NonlinearEquationsSolver {

    /**
     * Why the iteration stopped.
     * <p>
     * Only {@link #SOLVED} means the equations are satisfied. The others say
     * what stood in the way, and the point that comes back with them is still
     * the best one found -- {@link Result#residualNorm} says how good that is.
     *
     * @since 1.5.2
     */
    public enum Status {

        /**
         * The residual fell to or below the requested tolerance. This is the
         * only outcome that answers the question that was asked.
         */
        SOLVED(true),
        /**
         * The point stopped moving -- every component changed by less than the
         * step tolerance, relative to its own magnitude -- while the residual
         * is still above the residual tolerance. Usually this means the
         * requested residual tolerance is below what the arithmetic can deliver
         * for this system, and the point is as good as it gets; compare
         * {@link Result#residualNorm} against the residual at the starting
         * point before discarding it.
         */
        STEP_TOLERANCE_REACHED(false),
        /**
         * The iteration reached a local minimum of {@code 0.5 * F'F} at which
         * {@code F} is not zero. The scaled gradient of the merit function is
         * negligible there, so no further step of any length can reduce it:
         * the search is not stuck, it has arrived somewhere that is not a root.
         * More iterations do not help; a different starting point might, and a
         * system with no root at all ends here every time.
         */
        SPURIOUS_MINIMUM(false),
        /**
         * The Jacobian could not be used to produce a step, and neither could
         * the gradient of the merit function. On the {@code newton} path a
         * singular Jacobian on its own does not end the run -- the gradient
         * takes over -- so this means the gradient vanished as well, or the
         * Jacobian was not finite. On the {@code broyden} path it covers the
         * first Jacobian, or a restart, being impossible to invert.
         */
        SINGULAR_JACOBIAN(false),
        /**
         * The line search could not find a shorter step that reduces the merit
         * function, and on the {@code broyden} path rebuilding the Jacobian did
         * not change that. The direction was a descent direction, so this is
         * arithmetic rather than geometry: the residuals are not computed
         * accurately enough to see the decrease that is there.
         */
        LINE_SEARCH_STALLED(false),
        /** The iteration budget ran out with the residual still too large. */
        ITERATION_LIMIT(false),
        /**
         * The residuals at the starting point are not all finite, so there was
         * nothing to start from. A non-finite value encountered later is not
         * fatal -- the line search simply backtracks out of it.
         */
        NOT_FINITE(false);

        private final boolean success;

        private Status(boolean success) {
            this.success = success;
        }

        /**
         * Whether the equations are satisfied to the requested tolerance.
         *
         * @return {@code true} for {@link #SOLVED} and nothing else
         */
        public boolean isSuccess() {
            return success;
        }
    }

    /**
     * The outcome of a solve.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** The point reached, of length {@code n}. */
        public final double[] solution;
        /** The residuals {@code F} at {@link #solution}, of length {@code n}. */
        public final double[] residuals;
        /**
         * The largest absolute residual, {@code ||F||} in the infinity norm.
         * The infinity norm rather than the Euclidean one so that a residual
         * far too large to square still produces a number that can be compared.
         */
        public final double residualNorm;
        /** How many steps were taken. */
        public final int iterations;
        /**
         * How often {@code F} was evaluated, including the evaluations spent
         * inside the line search and inside an approximated Jacobian.
         */
        public final int functionEvaluations;
        /**
         * How often the Jacobian was asked for. On the derivative-free paths
         * each of these costs a further {@code n + 1} entries in
         * {@link #functionEvaluations}.
         */
        public final int jacobianEvaluations;
        /**
         * How often Broyden's approximation had to be thrown away and rebuilt
         * from a fresh Jacobian. Always {@code 0} on the {@code newton} path,
         * which keeps no approximation to throw away.
         */
        public final int restarts;
        /** Why the iteration stopped. */
        public final Status status;
        /** Shorthand for {@code status.isSuccess()}. */
        public final boolean converged;

        private Result(double[] solution, double[] residuals, double residualNorm, int iterations,
                int functionEvaluations, int jacobianEvaluations, int restarts, Status status) {
            this.solution = solution;
            this.residuals = residuals;
            this.residualNorm = residualNorm;
            this.iterations = iterations;
            this.functionEvaluations = functionEvaluations;
            this.jacobianEvaluations = jacobianEvaluations;
            this.restarts = restarts;
            this.status = status;
            this.converged = status.isSuccess();
        }
    }

    /** Default largest absolute residual that counts as a root. */
    public static final double DEFAULT_RESIDUAL_TOLERANCE = 1.0e-12;

    /** Default relative change below which the point counts as settled. */
    public static final double DEFAULT_STEP_TOLERANCE = 4.0 * MathConsts.MACH_EPS_DBL;

    /** Default budget of steps. */
    public static final int DEFAULT_MAX_ITERATIONS = 200;

    /**
     * The scaled gradient of the merit function below which a settled point is
     * a minimum of {@code 0.5 * F'F} rather than a stuck iteration.
     */
    private static final double GRADIENT_TOLERANCE = 1.0e-12;

    /**
     * A settled point is only called a spurious minimum when its residual is
     * still this much larger than the arithmetic could explain. Without it the
     * test misfires at a genuine root: there the gradient is negligible too,
     * because {@code F} itself is.
     */
    private static final double SPURIOUS_RESIDUAL_FACTOR = Math.sqrt(MathConsts.MACH_EPS_DBL);

    /** Armijo's sufficient decrease constant, the value everyone uses. */
    private static final double ARMIJO = 1.0e-4;

    /** The first step is at most this many times the size of the point. */
    private static final double STEP_BOUND_FACTOR = 100.0;

    private final double residualTolerance;
    private final double stepTolerance;
    private final int maxIterations;

    /**
     * Creates a solver with the default tolerances and budget.
     */
    public NonlinearEquationsSolver() {
        this(DEFAULT_RESIDUAL_TOLERANCE, DEFAULT_STEP_TOLERANCE, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a solver with the given tolerances and budget.
     *
     * @param residualTolerance
     *            the largest absolute residual that counts as a root, finite
     *            and greater than {@code 0}. It is <em>absolute</em>, so it
     *            lives in the units of {@code F}: multiplying every equation by
     *            a million leaves the root exactly where it was but raises the
     *            smallest reachable residual by the same factor, and the
     *            default then cannot be met by a system that meets it
     *            comfortably unscaled. Asking for less than the arithmetic can
     *            deliver does not produce a wrong answer, it produces
     *            {@link Status#STEP_TOLERANCE_REACHED} or
     *            {@link Status#LINE_SEARCH_STALLED} at a point that is
     *            nonetheless correct
     * @param stepTolerance
     *            the relative change below which the point counts as settled,
     *            finite and greater than {@code 0}. It is compared against
     *            {@code |dx_i| / max(|x_i|, 1)} taken over all components, so
     *            it is relative for a large component and absolute for a small
     *            one
     * @param maxIterations
     *            budget of steps, {@code 1} or greater. Line search
     *            evaluations do not count against it
     * @throws IllegalArgumentException
     *             if a tolerance is not finite and positive or if
     *             {@code maxIterations} is smaller than {@code 1}
     */
    public NonlinearEquationsSolver(double residualTolerance, double stepTolerance, int maxIterations) {
        if (!(residualTolerance > 0.0) || Double.isInfinite(residualTolerance)) {
            throw new IllegalArgumentException(
                    "residualTolerance must be finite and positive : " + residualTolerance);
        }
        if (!(stepTolerance > 0.0) || Double.isInfinite(stepTolerance)) {
            throw new IllegalArgumentException("stepTolerance must be finite and positive : " + stepTolerance);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException("maxIterations must be >= 1 : " + maxIterations);
        }
        this.residualTolerance = residualTolerance;
        this.stepTolerance = stepTolerance;
        this.maxIterations = maxIterations;
    }

    /**
     * Solves {@code F(x) = 0} by Newton's method, starting from
     * {@code initial}.
     * <p>
     * Each step solves {@code J(x) s = -F(x)} and then shortens {@code s} until
     * the merit function decreases. Where the Jacobian is numerically singular,
     * or where round-off makes the computed direction point uphill, the step of
     * that iteration is the gradient of the merit function instead: a worse
     * direction, but one that always exists, so a Jacobian that loses rank on
     * the way does not end the run.
     *
     * @param function
     *            the equations and their Jacobian, mapping {@code n} arguments
     *            to {@code n} results
     * @param initial
     *            the starting point, of length {@code n}, all finite (not
     *            modified)
     * @return the point reached, the residuals there, and why the iteration
     *         stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, or not finite
     */
    public Result newton(DiffDVectorFunction function, double[] initial) {
        checkArguments(function, initial);
        Run run = new Run(counting(function), initial, residualTolerance, stepTolerance);
        return newton(run);
    }

    /**
     * Solves {@code F(x) = 0} by Newton's method without being given a
     * derivative, approximating the Jacobian by forward differences instead.
     * <p>
     * Each Jacobian costs {@code n} further evaluations of {@code F} and is
     * accurate to about the square root of the machine epsilon, which usually
     * costs a step or two of the quadratic convergence near the root but rarely
     * the root itself. Note that Java picks between this method and
     * {@link #newton(DiffDVectorFunction, double[])} by the <em>static</em>
     * type of the argument, so a {@link DiffDVectorFunction} held in a variable
     * declared {@link DVectorFunction} lands here and its derivative goes
     * unused.
     *
     * @param function
     *            the equations, mapping {@code n} arguments to {@code n}
     *            results
     * @param initial
     *            the starting point, of length {@code n}, all finite (not
     *            modified)
     * @return the point reached, the residuals there, and why the iteration
     *         stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, or not finite
     */
    public Result newton(DVectorFunction function, double[] initial) {
        checkArguments(function, initial);
        Run run = new Run(approximating(function, initial.length), initial, residualTolerance, stepTolerance);
        return newton(run);
    }

    /**
     * Solves {@code F(x) = 0} by Broyden's method, starting from
     * {@code initial}.
     * <p>
     * The Jacobian is evaluated and inverted once, at the starting point, and
     * from then on the inverse is corrected by a rank one update that costs
     * {@code O(n^2)} rather than the {@code O(n^3)} of a factorization. That is
     * the whole point of the method and it only pays when {@code F} or its
     * derivative is expensive; on a cheap system Newton reaches the root in
     * fewer steps and usually in less time.
     * <p>
     * When the update degrades far enough that the line search can no longer
     * make progress, the approximation is thrown away and rebuilt from a fresh
     * Jacobian at the current point. Two failures in a row end the run.
     *
     * @param function
     *            the equations and their Jacobian, mapping {@code n} arguments
     *            to {@code n} results
     * @param initial
     *            the starting point, of length {@code n}, all finite (not
     *            modified)
     * @return the point reached, the residuals there, and why the iteration
     *         stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, or not finite
     */
    public Result broyden(DiffDVectorFunction function, double[] initial) {
        checkArguments(function, initial);
        Run run = new Run(counting(function), initial, residualTolerance, stepTolerance);
        return broyden(run);
    }

    /**
     * Solves {@code F(x) = 0} by Broyden's method without being given a
     * derivative, approximating the one Jacobian it needs by forward
     * differences.
     * <p>
     * This is the cheapest of the four entry points in evaluations of
     * {@code F}: one approximated Jacobian at the start, one further evaluation
     * per accepted step, and whatever the line search rejects along the way.
     * The same static-type caveat as in
     * {@link #newton(DVectorFunction, double[])} applies.
     *
     * @param function
     *            the equations, mapping {@code n} arguments to {@code n}
     *            results
     * @param initial
     *            the starting point, of length {@code n}, all finite (not
     *            modified)
     * @return the point reached, the residuals there, and why the iteration
     *         stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, or not finite
     */
    public Result broyden(DVectorFunction function, double[] initial) {
        checkArguments(function, initial);
        Run run = new Run(approximating(function, initial.length), initial, residualTolerance, stepTolerance);
        return broyden(run);
    }

    private Result newton(Run run) {
        int n = run.n;
        Status initial = run.start();
        if (initial != null) {
            return run.finish(initial, 0);
        }
        double[] jacobian = new double[n * n];
        double[] factor = new double[n * n];
        double[] gradient = new double[n];
        double[] direction = new double[n];
        int[] pivots = new int[n];
        double[] xNext = new double[n];
        double[] fNext = new double[n];

        for (int iteration = 1; iteration <= maxIterations; ++iteration) {
            run.jacobianAt(run.x, jacobian);
            transposeTimes(jacobian, run.fx, gradient, n);

            double slope = newtonDirection(jacobian, factor, pivots, run, direction, gradient);
            if (Double.isNaN(slope)) {
                // the gradient of the merit function has nothing left to offer
                // either, which at a point where F is not zero is the
                // definition of a spurious minimum rather than a mishap
                return run.finish(settled(run, Status.SINGULAR_JACOBIAN), iteration - 1);
            }

            boolean stepped = lineSearch(run, direction, slope, xNext, fNext);
            Status stopped = run.advance(stepped, xNext, fNext, direction);
            if (stopped != null) {
                return run.finish(stopped == Status.SOLVED ? Status.SOLVED : settled(run, stopped), iteration);
            }
        }
        return run.finish(settled(run, Status.ITERATION_LIMIT), maxIterations);
    }

    private Result broyden(Run run) {
        int n = run.n;
        Status initial = run.start();
        if (initial != null) {
            return run.finish(initial, 0);
        }
        double[] inverse = new double[n * n];
        double[] jacobian = new double[n * n];
        double[] factor = new double[n * n];
        int[] pivots = new int[n];
        double[] direction = new double[n];
        double[] xNext = new double[n];
        double[] fNext = new double[n];
        double[] change = new double[n];
        double[] inverseTimesChange = new double[n];
        double[] leftProduct = new double[n];

        if (!invertJacobianAt(run, jacobian, factor, pivots, inverse)) {
            return run.finish(settled(run, Status.SINGULAR_JACOBIAN), 0);
        }
        boolean rebuilt = false;
        boolean retried = false;

        for (int iteration = 1; iteration <= maxIterations; ++iteration) {
            timesNegated(inverse, run.fx, direction, n);
            // with an exact inverse the identity grad(f)'s = -F'F holds without
            // computing the gradient at all, and after an update it is still
            // the sufficient decrease the step ought to produce
            double slope = -2.0 * run.merit;

            boolean stepped = lineSearch(run, direction, slope, xNext, fNext);
            if (!stepped) {
                if (rebuilt) {
                    return run.finish(settled(run, Status.LINE_SEARCH_STALLED), iteration - 1);
                }
                if (!invertJacobianAt(run, jacobian, factor, pivots, inverse)) {
                    return run.finish(settled(run, Status.SINGULAR_JACOBIAN), iteration - 1);
                }
                run.restarts++;
                rebuilt = true;
                continue;
            }
            rebuilt = false;

            for (int i = 0; i < n; ++i) {
                change[i] = fNext[i] - run.fx[i];
            }
            Status stopped = run.advance(true, xNext, fNext, direction);
            if (stopped == Status.SOLVED) {
                return run.finish(Status.SOLVED, iteration);
            }
            if (stopped != null) {
                // the point stopped moving while the equations are still
                // unsatisfied. On this path that is as often a worn out
                // approximation as it is a real obstacle, so it is worth one
                // fresh Jacobian before the run is given up on -- but only one,
                // or a genuinely settled point would restart its way through
                // the whole budget
                if (retried || !invertJacobianAt(run, jacobian, factor, pivots, inverse)) {
                    return run.finish(settled(run, stopped), iteration);
                }
                run.restarts++;
                retried = true;
                rebuilt = true;
                continue;
            }
            update(inverse, run.step, change, inverseTimesChange, leftProduct, n);
        }
        return run.finish(settled(run, Status.ITERATION_LIMIT), maxIterations);
    }

    /**
     * The Newton direction, or the gradient direction where the Jacobian cannot
     * supply one. Returns the slope of the merit function along it, or
     * {@code Double.NaN} when even the fallback is unavailable, which can only
     * happen at the first iteration.
     */
    private static double newtonDirection(double[] jacobian, double[] factor, int[] pivots, Run run,
            double[] direction, double[] gradient) {
        int n = run.n;
        boolean usable = isFinite(jacobian);
        if (usable) {
            System.arraycopy(jacobian, 0, factor, 0, n * n);
            for (int i = 0; i < n; ++i) {
                direction[i] = -run.fx[i];
            }
            usable = Dgesv.dgesv(n, 1, factor, 0, n, pivots, 0, direction, 0, n)
                    && LinearEquationsSolver.diagonalRatio(factor, n, n, n) > LinearEquationsSolver
                            .rankThreshold(n, n)
                    && isFinite(direction);
        }
        double slope = usable ? dot(gradient, direction, n) : 0.0;
        if (usable && slope < 0.0) {
            return slope;
        }
        // the Jacobian is singular, or round-off on an ill conditioned one
        // turned the Newton step into an ascent direction. The gradient of the
        // merit function always descends where it does not vanish, and scaling
        // it to the size of the point keeps the line search from starting with
        // a step whose length has nothing to do with the problem
        double norm = twoNorm(gradient, n);
        if (!(norm > 0.0) || !isFinite(gradient)) {
            return Double.NaN;
        }
        double size = Math.max(twoNorm(run.x, n), n) / norm;
        for (int i = 0; i < n; ++i) {
            direction[i] = -size * gradient[i];
        }
        return -size * norm * norm;
    }

    /**
     * Evaluates the Jacobian at the current point and inverts it into
     * {@code inverse}. Returns {@code false} if it is not usable.
     */
    private static boolean invertJacobianAt(Run run, double[] jacobian, double[] factor, int[] pivots,
            double[] inverse) {
        int n = run.n;
        run.jacobianAt(run.x, jacobian);
        if (!isFinite(jacobian)) {
            return false;
        }
        System.arraycopy(jacobian, 0, factor, 0, n * n);
        for (int i = 0; i < n * n; ++i) {
            inverse[i] = 0.0;
        }
        for (int i = 0; i < n; ++i) {
            inverse[i * n + i] = 1.0;
        }
        return Dgesv.dgesv(n, n, factor, 0, n, pivots, 0, inverse, 0, n)
                && LinearEquationsSolver.diagonalRatio(factor, n, n, n) > LinearEquationsSolver.rankThreshold(n, n)
                && isFinite(inverse);
    }

    /**
     * Broyden's rank one correction of the inverse, the Sherman-Morrison form
     * of the update that makes {@code B} reproduce the change in {@code F} over
     * the step just taken:
     * {@code H += ((dx - H dy) (dx' H)) / (dx' H dy)}.
     */
    private static void update(double[] inverse, double[] step, double[] change, double[] inverseTimesChange,
            double[] leftProduct, int n) {
        times(inverse, change, inverseTimesChange, n);
        transposeTimes(inverse, step, leftProduct, n);
        double denominator = dot(step, inverseTimesChange, n);
        // the update is a division by a number that is a difference, so it can
        // be negligible relative to what it divides. Skipping leaves the
        // previous inverse in place, which is the better of the two
        double scale = twoNorm(step, n) * twoNorm(inverseTimesChange, n);
        if (!(Math.abs(denominator) > MathConsts.MACH_EPS_DBL * scale) || !isFinite(leftProduct)) {
            return;
        }
        for (int j = 0; j < n; ++j) {
            double right = leftProduct[j] / denominator;
            if (right == 0.0) {
                continue;
            }
            int column = j * n;
            for (int i = 0; i < n; ++i) {
                inverse[column + i] += (step[i] - inverseTimesChange[i]) * right;
            }
        }
    }

    /**
     * Backtracking line search on {@code 0.5 * F'F}, quadratic on the first
     * backtrack and cubic afterwards. Returns {@code false} when it ran below
     * the shortest step that could still move the point, in which case
     * {@code xNext} and {@code fNext} are not meaningful.
     */
    private static boolean lineSearch(Run run, double[] direction, double slope, double[] xNext, double[] fNext) {
        int n = run.n;
        double[] x = run.x;
        double meritOld = run.merit;

        double length = twoNorm(direction, n);
        if (length > run.stepBound) {
            double shrink = run.stepBound / length;
            for (int i = 0; i < n; ++i) {
                direction[i] *= shrink;
            }
            slope *= shrink;
        }

        double relative = 0.0;
        for (int i = 0; i < n; ++i) {
            double component = Math.abs(direction[i]) / Math.max(Math.abs(x[i]), 1.0);
            if (component > relative) {
                relative = component;
            }
        }
        double smallest = run.stepTolerance / relative;

        double lambda = 1.0;
        double previousLambda = 0.0;
        double previousMerit = 0.0;

        for (;;) {
            for (int i = 0; i < n; ++i) {
                xNext[i] = x[i] + lambda * direction[i];
            }
            run.valueAt(xNext, fNext);
            double merit = merit(fNext, n);
            if (merit <= meritOld + ARMIJO * lambda * slope) {
                run.trialMerit = merit;
                for (int i = 0; i < n; ++i) {
                    run.step[i] = lambda * direction[i];
                }
                return true;
            }
            if (lambda <= smallest) {
                return false;
            }

            double next;
            if (!Double.isFinite(merit)) {
                // no model can be fitted through a value that is not a number,
                // and the point of backtracking out of one is to leave quickly
                next = 0.1 * lambda;
            } else if (lambda == 1.0) {
                next = -slope / (2.0 * (merit - meritOld - slope));
            } else {
                double first = merit - meritOld - lambda * slope;
                double second = previousMerit - meritOld - previousLambda * slope;
                double a = (first / (lambda * lambda) - second / (previousLambda * previousLambda))
                        / (lambda - previousLambda);
                double b = (-previousLambda * first / (lambda * lambda)
                        + lambda * second / (previousLambda * previousLambda)) / (lambda - previousLambda);
                if (a == 0.0) {
                    next = -slope / (2.0 * b);
                } else {
                    double discriminant = b * b - 3.0 * a * slope;
                    if (discriminant < 0.0) {
                        next = 0.5 * lambda;
                    } else if (b <= 0.0) {
                        next = (-b + Math.sqrt(discriminant)) / (3.0 * a);
                    } else {
                        next = -slope / (b + Math.sqrt(discriminant));
                    }
                }
                if (next > 0.5 * lambda) {
                    next = 0.5 * lambda;
                }
            }
            previousLambda = lambda;
            previousMerit = merit;
            // negated so that a model that produced a NaN lands on the floor
            // rather than carrying it into the next round
            lambda = (next >= 0.1 * lambda) ? next : 0.1 * lambda;
        }
    }

    /**
     * Which of the two settled outcomes to report. The Jacobian is evaluated
     * once more for it, which happens only on a path that is not returning a
     * root anyway.
     */
    private static Status settled(Run run, Status stopped) {
        int n = run.n;
        double[] jacobian = new double[n * n];
        double[] gradient = new double[n];
        run.jacobianAt(run.x, jacobian);
        if (!isFinite(jacobian)) {
            return stopped;
        }
        transposeTimes(jacobian, run.fx, gradient, n);
        double denominator = Math.max(run.merit, 0.5 * n);
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double component = Math.abs(gradient[i]) * Math.max(Math.abs(run.x[i]), 1.0) / denominator;
            if (component > worst) {
                worst = component;
            }
        }
        boolean stillLarge = infinityNorm(run.fx, n) > SPURIOUS_RESIDUAL_FACTOR * Math.max(1.0, run.initialResidual);
        return (worst <= GRADIENT_TOLERANCE && stillLarge) ? Status.SPURIOUS_MINIMUM : stopped;
    }

    private static void checkArguments(DVectorFunction function, double[] initial) {
        if (function == null) {
            throw new IllegalArgumentException("function is null");
        }
        if (initial == null) {
            throw new IllegalArgumentException("initial is null");
        }
        if (initial.length == 0) {
            throw new IllegalArgumentException("initial is empty");
        }
        for (int i = 0; i < initial.length; ++i) {
            if (!Double.isFinite(initial[i])) {
                throw new IllegalArgumentException("initial[" + i + "] is not finite : " + initial[i]);
            }
        }
    }

    // --- the arithmetic, written out rather than taken from math.linalg.
    // VectorOps exists twice, and its Java 25 twin accumulates in lane order
    // and multiplies with fma, so a solver built on it would take a different
    // path on the two releases. These loops take the same one. ---

    private static double dot(double[] a, double[] b, int n) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    private static double twoNorm(double[] a, int n) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += a[i] * a[i];
        }
        return Math.sqrt(sum);
    }

    private static double infinityNorm(double[] a, int n) {
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double component = Math.abs(a[i]);
            if (Double.isNaN(component)) {
                return Double.NaN;
            }
            if (component > worst) {
                worst = component;
            }
        }
        return worst;
    }

    private static double merit(double[] f, int n) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += f[i] * f[i];
        }
        return 0.5 * sum;
    }

    private static boolean isFinite(double[] a) {
        for (int i = 0; i < a.length; ++i) {
            if (!Double.isFinite(a[i])) {
                return false;
            }
        }
        return true;
    }

    /** {@code out = A' v} for a column-major {@code n x n} matrix. */
    private static void transposeTimes(double[] a, double[] v, double[] out, int n) {
        for (int j = 0; j < n; ++j) {
            out[j] = dotColumn(a, j * n, v, n);
        }
    }

    private static double dotColumn(double[] a, int offset, double[] v, int n) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += a[offset + i] * v[i];
        }
        return sum;
    }

    /** {@code out = A v} for a column-major {@code n x n} matrix. */
    private static void times(double[] a, double[] v, double[] out, int n) {
        for (int i = 0; i < n; ++i) {
            out[i] = 0.0;
        }
        for (int j = 0; j < n; ++j) {
            double factor = v[j];
            if (factor == 0.0) {
                continue;
            }
            int column = j * n;
            for (int i = 0; i < n; ++i) {
                out[i] += a[column + i] * factor;
            }
        }
    }

    /** {@code out = -A v}, the Broyden step. */
    private static void timesNegated(double[] a, double[] v, double[] out, int n) {
        times(a, v, out, n);
        for (int i = 0; i < n; ++i) {
            out[i] = -out[i];
        }
    }

    /** What the two wrappers below have in common. */
    private interface Counted extends DiffDVectorFunction {

        /** How often the function underneath was evaluated. */
        int evaluations();
    }

    private static Counted counting(DiffDVectorFunction function) {
        return new Counting(function);
    }

    private static Counted approximating(DVectorFunction function, int n) {
        return new Approximating(function, n);
    }

    /** Routes every evaluation past the counter without changing any of them. */
    private static final class Counting implements Counted {

        private final DiffDVectorFunction function;
        int evaluations;

        Counting(DiffDVectorFunction function) {
            this.function = function;
        }

        @Override
        public void valueAt(double[] x, double[] values) {
            ++evaluations;
            function.valueAt(x, values);
        }

        @Override
        public void jacobianAt(double[] x, double[] jacobian) {
            function.jacobianAt(x, jacobian);
        }

        @Override
        public int evaluations() {
            return evaluations;
        }
    }

    /** The same for a function that has to be differenced. */
    private static final class Approximating extends NumericalDiffDVectorFunction implements Counted {

        private final DVectorFunction function;
        int evaluations;

        Approximating(DVectorFunction function, int n) {
            super(n);
            this.function = function;
        }

        @Override
        public void valueAt(double[] x, double[] values) {
            ++evaluations;
            function.valueAt(x, values);
        }

        @Override
        public int evaluations() {
            return evaluations;
        }
    }

    /** Everything one call to {@code newton} or {@code broyden} carries. */
    private static final class Run {

        final Counted function;
        final int n;
        final double residualTolerance;
        final double stepTolerance;
        final double[] x;
        final double[] fx;
        final double[] step;
        final double stepBound;
        double merit;
        double trialMerit;
        double initialResidual;
        int jacobianEvaluations;
        int restarts;

        Run(Counted function, double[] initial, double residualTolerance, double stepTolerance) {
            this.function = function;
            this.n = initial.length;
            this.residualTolerance = residualTolerance;
            this.stepTolerance = stepTolerance;
            this.x = initial.clone();
            this.fx = new double[n];
            this.step = new double[n];
            this.stepBound = STEP_BOUND_FACTOR * Math.max(twoNorm(this.x, n), n);
        }

        void valueAt(double[] point, double[] values) {
            function.valueAt(point, values);
        }

        void jacobianAt(double[] point, double[] jacobian) {
            ++jacobianEvaluations;
            function.jacobianAt(point, jacobian);
        }

        /**
         * Evaluates the starting point. Returns the status to report at once,
         * or {@code null} when there is work to do.
         */
        Status start() {
            valueAt(x, fx);
            if (!isFinite(fx)) {
                merit = merit(fx, n);
                initialResidual = infinityNorm(fx, n);
                return Status.NOT_FINITE;
            }
            merit = merit(fx, n);
            initialResidual = infinityNorm(fx, n);
            return (initialResidual <= residualTolerance) ? Status.SOLVED : null;
        }

        /**
         * Takes the step the line search found. Returns the status to stop on,
         * or {@code null} to keep going.
         */
        Status advance(boolean stepped, double[] xNext, double[] fNext, double[] direction) {
            if (!stepped) {
                return Status.LINE_SEARCH_STALLED;
            }
            double settledBy = 0.0;
            for (int i = 0; i < n; ++i) {
                double component = Math.abs(xNext[i] - x[i]) / Math.max(Math.abs(xNext[i]), 1.0);
                if (component > settledBy) {
                    settledBy = component;
                }
            }
            double before = merit;
            System.arraycopy(xNext, 0, x, 0, n);
            System.arraycopy(fNext, 0, fx, 0, n);
            merit = trialMerit;
            if (infinityNorm(fx, n) <= residualTolerance) {
                return Status.SOLVED;
            }
            // the line search accepts nothing that does not satisfy Armijo, and
            // Armijo asks for a strict decrease -- so an accepted step that
            // leaves the merit exactly where it was is one whose required
            // decrease underflowed. Nothing shorter can do better and nothing
            // longer is acceptable, so the budget would only be spent proving
            // it: a system with no root sits here for a hundred thousand
            // iterations and moves its point by five parts in ten thousand
            if (!(merit < before)) {
                return Status.LINE_SEARCH_STALLED;
            }
            return (settledBy <= stepTolerance) ? Status.STEP_TOLERANCE_REACHED : null;
        }

        Result finish(Status status, int iterations) {
            return new Result(x.clone(), fx.clone(), infinityNorm(fx, n), iterations, function.evaluations(),
                    jacobianEvaluations, restarts, status);
        }
    }

    /**
     * Solves two small systems and checks the residuals it reports.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        // (1) Rosenbrock as a system, whose only root is (1, 1)
        DiffDVectorFunction rosenbrock = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = 10.0 * (x[1] - x[0] * x[0]);
                f[1] = 1.0 - x[0];
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                j[0] = -20.0 * x[0];
                j[1] = -1.0;
                j[2] = 10.0;
                j[3] = 0.0;
            }
        };
        // (2) two circles meeting at (1, 0) and at (0, 1)
        DVectorFunction circles = new DVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = x[0] * x[0] + x[1] * x[1] - 1.0;
                f[1] = x[0] + x[1] - 1.0;
            }
        };

        NonlinearEquationsSolver solver = new NonlinearEquationsSolver();
        Result newton = solver.newton(rosenbrock, new double[] { -1.2, 1.0 });
        Result broyden = solver.broyden(rosenbrock, new double[] { -1.2, 1.0 });
        Result difference = solver.newton(circles, new double[] { 2.0, 3.0 });

        System.out.println("Rosenbrock, newton  : x = (" + newton.solution[0] + ", " + newton.solution[1]
                + "), ||F|| = " + newton.residualNorm + ", " + newton.status + " after " + newton.iterations
                + " steps");
        System.out.println("Rosenbrock, broyden : x = (" + broyden.solution[0] + ", " + broyden.solution[1]
                + "), ||F|| = " + broyden.residualNorm + ", " + broyden.status + " after " + broyden.iterations
                + " steps");
        System.out.println("Circles, differenced: x = (" + difference.solution[0] + ", " + difference.solution[1]
                + "), ||F|| = " + difference.residualNorm + ", " + difference.status + " after "
                + difference.iterations + " steps");

        boolean ok = newton.converged && broyden.converged && difference.converged
                && Math.abs(newton.solution[0] - 1.0) < 1.0e-9 && Math.abs(newton.solution[1] - 1.0) < 1.0e-9
                && Math.abs(broyden.solution[0] - 1.0) < 1.0e-9 && Math.abs(broyden.solution[1] - 1.0) < 1.0e-9
                && Math.abs(difference.solution[0] * difference.solution[1]) < 1.0e-9;
        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }
}
