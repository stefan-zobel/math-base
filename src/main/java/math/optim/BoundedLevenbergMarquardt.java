package math.optim;

import math.MathConsts;
import math.fun.DVectorFunction;
import math.fun.DiffDVectorFunction;
import math.linalg.BoundedLeastSquares;
import math.linalg.DMatrix;
import math.linalg.FlatParallelJacobiSVD;

/**
 * Nonlinear least squares under hard bounds: given {@code m} residuals over
 * {@code n} parameters, minimize the sum of their squares subject to
 * {@code l <= x <= u}.
 * <p>
 * Where {@link LevenbergMarquardt} damps the Gauss-Newton step and solves
 * {@code min ||J s + r||^2 + lambda ||D s||^2} for it, this one solves the same
 * subproblem with the box carried along,
 * {@code l - x <= s <= u - x}, through
 * {@link BoundedLeastSquares}, which is exact and terminates finitely. The
 * restriction therefore sits where the step is formed rather than around the
 * search, and <strong>every point the model is evaluated at is inside the
 * box</strong> -- which matters because a bound usually exists precisely because
 * the model is undefined beyond it, a variance that must stay positive or a
 * logarithm of a positive argument.
 * <p>
 * The parameters are scaled internally from the columns of the Jacobian, so a
 * model whose parameters differ by orders of magnitude does not need to be
 * rescaled by hand.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm">
 *      Levenberg-Marquardt algorithm</a>
 * @since 1.5.2
 */
public final class BoundedLevenbergMarquardt {

    /**
     * Why the search stopped. The three tolerance outcomes are successes; a
     * point that is feasible comes back in every case.
     *
     * @since 1.5.2
     */
    public enum Status {

        /** The relative decrease of the sum of squares is at most {@code valueTolerance}. */
        VALUE_TOLERANCE_REACHED(true),

        /** The relative change of the parameters is at most {@code stepTolerance}. */
        STEP_TOLERANCE_REACHED(true),

        /**
         * The scaled optimality measure is at most {@code kktTolerance}, which
         * is the first-order condition of the constrained problem. Switched off
         * by default, as in {@link LevenbergMarquardt} and in MINPACK's own
         * drivers.
         */
        KKT_TOLERANCE_REACHED(true),

        /**
         * The evaluation budget ran out. This says nothing about the point that
         * comes back with it, which is the best feasible one found.
         */
        TOO_MANY_EVALUATIONS(false),

        /**
         * The residuals or the Jacobian were not finite at a point inside the
         * box, so the search cannot continue. A trial point where only the
         * residuals fail is not this: it is treated as a step that was too long
         * and the damping is raised instead.
         */
        NOT_FINITE(false),

        /**
         * No step short enough to reduce the sum of squares could be found
         * before the damping overflowed. The model is not locally quadratic
         * where the search ended, or the residuals are computed too
         * inaccurately for the progress being asked of them.
         */
        STEP_SOLVE_FAILED(false);

        private final boolean success;

        private Status(boolean success) {
            this.success = success;
        }

        /**
         * Whether the returned parameters are a minimum to the accuracy that
         * was asked for.
         *
         * @return {@code true} for the three tolerance outcomes
         */
        public boolean isSuccess() {
            return success;
        }
    }

    /**
     * The outcome of a bounded fit. {@link #parameters} is feasible whatever
     * the status, because every iterate is.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** The parameters found, of length {@code n}, inside the box. */
        public final double[] parameters;

        /** The residuals at {@link #parameters}, of length {@code m}. */
        public final double[] residuals;

        /** The sum of the squared {@link #residuals}, the quantity minimized. */
        public final double sumOfSquares;

        /**
         * Which constraint binds for each parameter. An entry of
         * {@link BoundedLeastSquares.Bound#AT_LOWER} means the parameter equals
         * its lower bound exactly, not approximately, so it can be tested with
         * {@code ==}.
         */
        public final BoundedLeastSquares.Bound[] activeSet;

        /** How many parameters came back at a bound rather than free. */
        public final int atBound;

        /**
         * The largest amount by which the gradient wants to move a parameter in
         * a direction the box allows: {@code |g_j|} where the parameter is
         * free, and the part of {@code g_j} pointing into the box where it is
         * held. Zero at an exact solution.
         * <p>
         * This is in the units of the gradient and therefore grows with the
         * size of the data; {@link #scaledKktViolation} is the version that
         * does not.
         */
        public final double kktViolation;

        /**
         * {@link #kktViolation} divided by {@code ||J_j|| ||r||}, which is
         * MINPACK's measure of the gradient and is dimensionless. This is the
         * quantity {@code kktTolerance} is compared against.
         * <p>
         * It degenerates at both ends, which is why no tolerance is applied to
         * it by default. Where the residual vanishes at the solution it is the
         * cosine of an angle to the zero vector and comes out large on a
         * perfect fit; where a Jacobian column vanishes -- a parameter the
         * residuals do not depend on at the point reached -- the denominator
         * does instead. Read {@link #kktViolation} in either case.
         */
        public final double scaledKktViolation;

        /** How many steps were accepted. */
        public final int iterations;

        /** How often the residuals were evaluated. */
        public final int functionEvaluations;

        /** How often the Jacobian was evaluated. */
        public final int jacobianEvaluations;

        /** Why the search stopped. */
        public final Status status;

        /** Shorthand for {@code status.isSuccess()}. */
        public final boolean converged;

        /**
         * {@code m} minus the number of free parameters. A parameter held at a
         * bound is not estimated from the data and consumes no degree of
         * freedom, which is where this differs from
         * {@link LevenbergMarquardt.Result#degreesOfFreedom}.
         */
        public final int degreesOfFreedom;

        /**
         * The estimated variance-covariance matrix of {@link #parameters},
         * {@code n} by {@code n}, flat and column-major, or {@code null} where
         * it cannot be formed.
         * <p>
         * It covers the <em>free</em> parameters only; the row and column of a
         * parameter held at a bound are zero, because a parameter fixed by a
         * restriction has no sampling variance under this linearization. The
         * same two caveats as on
         * {@link LevenbergMarquardt.Result#covariance} apply on top of that:
         * the residuals must be independent draws of equal variance, and for a
         * model nonlinear in its parameters this is the linearization at the
         * solution, exact only in the limit of many observations. A bootstrap
         * over the residuals answers the same question without either
         * assumption, and {@code math.probe.Bootstrap} is in this library.
         * <p>
         * {@code null} when no parameter is free, when
         * {@link #degreesOfFreedom} is not positive, and when the free block of
         * the Jacobian is numerically rank deficient. None of the three means
         * the fit failed.
         */
        public final double[] covariance;

        /**
         * The square roots of the diagonal of {@link #covariance}, of length
         * {@code n}, or {@code null} whenever that is. A parameter at a bound
         * carries {@code 0.0}.
         */
        public final double[] standardErrors;

        Result(double[] parameters, double[] residuals, double sumOfSquares,
                BoundedLeastSquares.Bound[] activeSet, double kktViolation, double scaledKktViolation,
                int iterations, int functionEvaluations, int jacobianEvaluations, Status status,
                int degreesOfFreedom, double[] covariance) {
            this.parameters = parameters;
            this.residuals = residuals;
            this.sumOfSquares = sumOfSquares;
            this.activeSet = activeSet;
            this.kktViolation = kktViolation;
            this.scaledKktViolation = scaledKktViolation;
            this.iterations = iterations;
            this.functionEvaluations = functionEvaluations;
            this.jacobianEvaluations = jacobianEvaluations;
            this.status = status;
            this.converged = status.isSuccess();
            this.degreesOfFreedom = degreesOfFreedom;
            this.covariance = covariance;
            int bound = 0;
            for (int j = 0; j < activeSet.length; j++) {
                if (activeSet[j] != BoundedLeastSquares.Bound.FREE) {
                    bound++;
                }
            }
            this.atBound = bound;
            if (covariance == null) {
                this.standardErrors = null;
            } else {
                int n = parameters.length;
                double[] errors = new double[n];
                for (int j = 0; j < n; j++) {
                    errors[j] = Math.sqrt(covariance[j * n + j]);
                }
                this.standardErrors = errors;
            }
        }

        @Override
        public String toString() {
            return "BoundedLevenbergMarquardt.Result[ssq = " + sumOfSquares + ", " + atBound + " of "
                    + activeSet.length + " at a bound, " + iterations + " steps, " + functionEvaluations
                    + " evaluations, " + status + "]";
        }
    }

    /**
     * Default for both the step and the value tolerance. MINPACK recommends
     * nothing below the square root of the machine epsilon, because the sum of
     * squares is flat at its minimum and half the digits are gone by the time
     * it is reached.
     */
    private static final double DEFAULT_TOLERANCE = Math.sqrt(MathConsts.MACH_EPS_DBL);

    /**
     * Default initial damping. Dimensionless, because the damping multiplies
     * the squared column norms of the Jacobian rather than the step itself;
     * {@code 1e-3} is the customary choice for a starting point believed to be
     * near the minimum, and larger values buy a shorter and safer first step.
     */
    private static final double DEFAULT_INITIAL_DAMPING = 1.0e-3;

    /** Stands for the budget of {@code 100 * (n + 1)} that MINPACK defaults to. */
    private static final int EVALUATIONS_FROM_PROBLEM_SIZE = 0;

    /** Residuals are assumed accurate to machine precision unless said otherwise. */
    private static final double DEFAULT_FUNCTION_ACCURACY = 0.0;

    /** Beyond this the damping has stopped meaning anything. */
    private static final double MAX_DAMPING = 1.0e250;

    private final double valueTolerance;
    private final double stepTolerance;
    private final double kktTolerance;
    private final int maxEvaluations;
    private final double initialDamping;
    private final double functionAccuracy;

    /**
     * Creates a solver with the default stopping rules: both the value and the
     * step tolerance at the square root of the machine epsilon, no optimality
     * tolerance, and a budget of {@code 100 * (n + 1)} evaluations.
     */
    public BoundedLevenbergMarquardt() {
        this(DEFAULT_TOLERANCE, DEFAULT_TOLERANCE, 0.0, EVALUATIONS_FROM_PROBLEM_SIZE, DEFAULT_INITIAL_DAMPING);
    }

    /**
     * Creates a solver with explicit stopping rules. The search stops as soon
     * as any one of the three tolerances is met, so a tolerance of {@code 0}
     * switches that test off rather than demanding the impossible.
     *
     * @param valueTolerance
     *            relative decrease of the sum of squares that is accepted,
     *            {@code 0} or greater
     * @param stepTolerance
     *            relative change of the parameters that is accepted, {@code 0}
     *            or greater
     * @param kktTolerance
     *            scaled optimality measure that is accepted, {@code 0} or
     *            greater. {@code 0} switches the test off, which is the
     *            default and what MINPACK's own drivers do, because the measure
     *            is a cosine and degenerates on a problem whose residual
     *            vanishes at the minimum
     * @param maxEvaluations
     *            budget for evaluations of the residuals, {@code 1} or greater
     * @param initialDamping
     *            the damping to start from, greater than {@code 0}. It is
     *            dimensionless: it multiplies the squared column norms of the
     *            Jacobian, so it does not have to be rescaled for the problem.
     *            Raise it towards {@code 1} when the starting point is far from
     *            the minimum and the first step overshoots
     */
    public BoundedLevenbergMarquardt(double valueTolerance, double stepTolerance, double kktTolerance,
            int maxEvaluations, double initialDamping) {
        this(valueTolerance, stepTolerance, kktTolerance, maxEvaluations, initialDamping,
                DEFAULT_FUNCTION_ACCURACY);
    }

    /**
     * Creates a solver that also knows how accurately the residuals themselves
     * can be computed, which matters only on the derivative-free path.
     *
     * @param valueTolerance
     *            as above
     * @param stepTolerance
     *            as above
     * @param kktTolerance
     *            as above
     * @param maxEvaluations
     *            as above
     * @param initialDamping
     *            as above
     * @param functionAccuracy
     *            the relative error with which {@code valueAt} computes the
     *            residuals, {@code 0} or greater, {@code 0} meaning machine
     *            precision. Ignored by
     *            {@link #solve(DiffDVectorFunction, double[], int, double[], double[])},
     *            which is handed the derivative and never has to guess a step
     *            for it.
     *            {@link #solve(DVectorFunction, double[], int, double[], double[])}
     *            does have to, and it picks that step from this number: a
     *            residual that comes out of a simulation, a quadrature or a
     *            table is accurate to far less than machine precision, and
     *            differencing it over the default step measures the error
     *            rather than the derivative
     */
    public BoundedLevenbergMarquardt(double valueTolerance, double stepTolerance, double kktTolerance,
            int maxEvaluations, double initialDamping, double functionAccuracy) {
        if (!(valueTolerance >= 0.0) || Double.isInfinite(valueTolerance)) {
            throw new IllegalArgumentException("valueTolerance must be finite and non-negative : " + valueTolerance);
        }
        if (!(stepTolerance >= 0.0) || Double.isInfinite(stepTolerance)) {
            throw new IllegalArgumentException("stepTolerance must be finite and non-negative : " + stepTolerance);
        }
        if (!(kktTolerance >= 0.0) || Double.isInfinite(kktTolerance)) {
            throw new IllegalArgumentException("kktTolerance must be finite and non-negative : " + kktTolerance);
        }
        if (maxEvaluations < 1 && maxEvaluations != EVALUATIONS_FROM_PROBLEM_SIZE) {
            throw new IllegalArgumentException("maxEvaluations must be >= 1 : " + maxEvaluations);
        }
        if (!(initialDamping > 0.0) || Double.isInfinite(initialDamping)) {
            throw new IllegalArgumentException("initialDamping must be finite and positive : " + initialDamping);
        }
        if (!(functionAccuracy >= 0.0) || Double.isInfinite(functionAccuracy)) {
            throw new IllegalArgumentException(
                    "functionAccuracy must be finite and non-negative : " + functionAccuracy);
        }
        this.valueTolerance = valueTolerance;
        this.stepTolerance = stepTolerance;
        this.kktTolerance = kktTolerance;
        this.maxEvaluations = maxEvaluations;
        this.initialDamping = initialDamping;
        this.functionAccuracy = functionAccuracy;
    }

    /**
     * Minimizes the sum of the squared residuals of {@code function} over the
     * box, starting from {@code initial}.
     * <p>
     * The starting point is not modified, and it does not have to be feasible:
     * a point outside the box is moved onto it, which is what a caller who
     * wrote {@code 0} for a parameter bounded below by {@code 1e-6} meant.
     *
     * @param function
     *            the residuals and their Jacobian, evaluated only inside the
     *            box
     * @param initial
     *            the starting parameters, of length {@code n}, all finite (not
     *            modified)
     * @param residualCount
     *            the number {@code m} of residuals {@code function} produces,
     *            which must be at least {@code n}
     * @param lower
     *            lower bounds, length {@code n}, not {@code NaN}, may be
     *            {@code -inf} (not modified)
     * @param upper
     *            upper bounds, length {@code n}, not {@code NaN}, may be
     *            {@code +inf}, each {@code upper[j] >= lower[j]} (not
     *            modified). {@code lower[j] == upper[j]} pins that parameter
     * @return the parameters found, the residuals there, which constraints bind
     *         and why the search stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, not finite, if
     *             {@code residualCount} is smaller than {@code initial.length},
     *             or if a bound is {@code NaN} or the box is empty
     */
    public Result solve(DiffDVectorFunction function, double[] initial, int residualCount, double[] lower,
            double[] upper) {
        checkArguments(function, initial, residualCount, lower, upper);
        return run(new AnalyticModel(function), initial, residualCount, lower, upper);
    }

    /**
     * Minimizes the sum of the squared residuals of {@code function} over the
     * box without being given a derivative, approximating the Jacobian by
     * finite differences instead.
     * <p>
     * The difference is taken <em>inwards</em>: the step is subtracted rather
     * than added where adding it would cross the upper bound, and it shrinks
     * rather than crossing where the box is narrower than the step on both
     * sides. A bound commonly exists because the model is undefined beyond it,
     * and a difference quotient that probes outside would measure nonsense
     * without saying so.
     * <p>
     * The step itself is the magnitude of the parameter it is taken from,
     * <em>floored at one</em> -- the step both
     * {@link math.fun.NumericalDiffDMultiFunction} and
     * {@link math.fun.NumericalDiffDVectorFunction} take. The floor is what
     * keeps a parameter near zero, but not at zero, from being handed a step
     * so small that the residuals come back unchanged in every bit. That
     * column of zeros is indistinguishable from a parameter the residuals do
     * not depend on, so the parameter stays where it started and the search
     * reports success. Where the box is narrower than the step, the step
     * shrinks to fit as described above; the floor is not lowered to meet it.
     * <p>
     * If the residuals are not computed to machine precision, say so through
     * the {@code functionAccuracy} of the six-argument constructor. The step
     * for the difference is chosen from it, and left at the default it will be
     * far too small to see past the error in the values. That failure is silent
     * and it reports success: with a step below the resolution of the residuals
     * every difference comes out zero, the approximated Jacobian is the zero
     * matrix, and a zero matrix is orthogonal to every residual -- so the
     * optimality test is satisfied at the starting point, having moved nothing.
     * There is no way to tell that from a genuine minimum after the fact, which
     * is why it has to be declared beforehand.
     *
     * @param function
     *            the residuals, evaluated only inside the box
     * @param initial
     *            the starting parameters, of length {@code n}, all finite (not
     *            modified)
     * @param residualCount
     *            the number {@code m} of residuals {@code function} produces,
     *            which must be at least {@code n}
     * @param lower
     *            lower bounds, length {@code n}, as above
     * @param upper
     *            upper bounds, length {@code n}, as above
     * @return the parameters found, the residuals there, which constraints bind
     *         and why the search stopped. {@link Result#jacobianEvaluations}
     *         counts the approximations, and the {@code n} evaluations each of
     *         them costs show up in {@link Result#functionEvaluations}
     * @throws IllegalArgumentException
     *             as above
     */
    public Result solve(DVectorFunction function, double[] initial, int residualCount, double[] lower,
            double[] upper) {
        checkArguments(function, initial, residualCount, lower, upper);
        return run(new DifferenceModel(function, residualCount, initial.length, lower, upper, functionAccuracy),
                initial, residualCount, lower, upper);
    }

    private Result run(Model model, double[] initial, int m, double[] lower, double[] upper) {
        int n = initial.length;
        int rows = m + n;
        int budget = (maxEvaluations == EVALUATIONS_FROM_PROBLEM_SIZE) ? 100 * (n + 1) : maxEvaluations;

        double[] x = new double[n];
        for (int j = 0; j < n; j++) {
            x[j] = clamp(initial[j], lower[j], upper[j]);
        }
        double[] r = new double[m];
        model.valueAt(x, r);
        int functionEvaluations = 1;
        int jacobianEvaluations = 0;
        int iterations = 0;

        double[] jacobian = new double[m * n];
        double[] columnNorm = new double[n];
        double[] scaling = new double[n];
        double[] gradient = new double[n];
        double[] stepLower = new double[n];
        double[] stepUpper = new double[n];
        double[] trial = new double[n];
        double[] trialResiduals = new double[m];

        DMatrix design = new DMatrix(rows, n);
        DMatrix target = new DMatrix(rows, 1);
        double[] a = design.getArrayUnsafe();
        double[] b = target.getArrayUnsafe();

        double damping = initialDamping;
        double nu = 2.0;
        double kkt = Double.NaN;
        double scaledKkt = Double.NaN;
        Status status;

        if (!isFinite(r, m)) {
            status = Status.NOT_FINITE;
        } else {
            double sumOfSquares = sumOfSquares(r, m);
            status = Status.TOO_MANY_EVALUATIONS;

            outer:
            while (true) {
                jacobianEvaluations++;
                functionEvaluations += model.jacobianAt(x, r, jacobian);
                if (!isFinite(jacobian, m * n)) {
                    status = Status.NOT_FINITE;
                    break outer;
                }
                for (int j = 0; j < n; j++) {
                    columnNorm[j] = euclidean(jacobian, j * m, m);
                    // MINPACK's convention: the scaling never shrinks, and a
                    // column that carries no derivative is scaled by one rather
                    // than by nothing
                    scaling[j] = Math.max(scaling[j], columnNorm[j]);
                    if (scaling[j] == 0.0) {
                        scaling[j] = 1.0;
                    }
                    double sum = 0.0;
                    int column = j * m;
                    for (int i = 0; i < m; i++) {
                        sum += jacobian[column + i] * r[i];
                    }
                    gradient[j] = sum;
                }

                kkt = kktViolation(x, gradient, lower, upper);
                scaledKkt = scaledKktViolation(x, gradient, columnNorm, lower, upper, Math.sqrt(sumOfSquares));
                if (scaledKkt <= kktTolerance) {
                    status = Status.KKT_TOLERANCE_REACHED;
                    break outer;
                }

                // the part of the augmented design that does not depend on the
                // damping, and the shifted box, which always contains zero
                for (int j = 0; j < n; j++) {
                    System.arraycopy(jacobian, j * m, a, j * rows, m);
                    for (int i = m; i < rows; i++) {
                        a[j * rows + i] = 0.0;
                    }
                    stepLower[j] = lower[j] - x[j];
                    stepUpper[j] = upper[j] - x[j];
                }
                for (int i = 0; i < m; i++) {
                    b[i] = -r[i];
                }
                for (int i = m; i < rows; i++) {
                    b[i] = 0.0;
                }

                while (true) {
                    double root = Math.sqrt(damping);
                    for (int j = 0; j < n; j++) {
                        a[j * rows + m + j] = root * scaling[j];
                    }
                    BoundedLeastSquares.Result step = BoundedLeastSquares.bounded(design, target, stepLower,
                            stepUpper);
                    double[] s = step.solution;
                    // the first m residuals of the augmented system are
                    // -(J s + r), so the predicted reduction is already there
                    double predicted = sumOfSquares - sumOfSquares(step.residuals, m);

                    double stepNorm = 0.0;
                    for (int j = 0; j < n; j++) {
                        double moved;
                        if (s[j] == stepLower[j]) {
                            moved = lower[j];
                        } else if (s[j] == stepUpper[j]) {
                            moved = upper[j];
                        } else {
                            // the sum can still round outside by an ulp
                            moved = clamp(x[j] + s[j], lower[j], upper[j]);
                        }
                        trial[j] = moved;
                        double delta = moved - x[j];
                        stepNorm += delta * delta;
                    }
                    stepNorm = Math.sqrt(stepNorm);

                    model.valueAt(trial, trialResiduals);
                    functionEvaluations++;
                    boolean usable = isFinite(trialResiduals, m);
                    double trialSumOfSquares = usable ? sumOfSquares(trialResiduals, m) : Double.POSITIVE_INFINITY;
                    boolean better = usable && trialSumOfSquares < sumOfSquares;
                    double gain = (predicted > 0.0) ? (sumOfSquares - trialSumOfSquares) / predicted
                            : (better ? 1.0 : -1.0);

                    if (better && gain > 0.0) {
                        double reduction = sumOfSquares - trialSumOfSquares;
                        double previous = sumOfSquares;
                        double norm = euclidean(x, 0, n);
                        System.arraycopy(trial, 0, x, 0, n);
                        System.arraycopy(trialResiduals, 0, r, 0, m);
                        sumOfSquares = trialSumOfSquares;
                        iterations++;
                        // Nielsen's rule: the better the step predicted itself,
                        // the further the damping falls
                        double excess = 2.0 * gain - 1.0;
                        damping = damping * Math.max(1.0 / 3.0, 1.0 - excess * excess * excess);
                        nu = 2.0;
                        if (reduction <= valueTolerance * previous) {
                            status = Status.VALUE_TOLERANCE_REACHED;
                            break outer;
                        }
                        if (stepNorm <= stepTolerance * (norm + stepTolerance)) {
                            status = Status.STEP_TOLERANCE_REACHED;
                            break outer;
                        }
                        break;
                    }

                    damping = damping * nu;
                    nu = 2.0 * nu;
                    if (functionEvaluations >= budget) {
                        status = Status.TOO_MANY_EVALUATIONS;
                        break outer;
                    }
                    if (damping > MAX_DAMPING) {
                        status = Status.STEP_SOLVE_FAILED;
                        break outer;
                    }
                }

                if (functionEvaluations >= budget) {
                    status = Status.TOO_MANY_EVALUATIONS;
                    break outer;
                }
            }
        }

        BoundedLeastSquares.Bound[] activeSet = new BoundedLeastSquares.Bound[n];
        int free = 0;
        for (int j = 0; j < n; j++) {
            if (x[j] == lower[j]) {
                activeSet[j] = BoundedLeastSquares.Bound.AT_LOWER;
            } else if (x[j] == upper[j]) {
                activeSet[j] = BoundedLeastSquares.Bound.AT_UPPER;
            } else {
                activeSet[j] = BoundedLeastSquares.Bound.FREE;
                free++;
            }
        }
        double finalSumOfSquares = isFinite(r, m) ? sumOfSquares(r, m) : Double.NaN;
        return new Result(x, r, finalSumOfSquares, activeSet, kkt, scaledKkt, iterations, functionEvaluations,
                jacobianEvaluations, status, m - free,
                covarianceOf(jacobian, activeSet, m, n, free, finalSumOfSquares));
    }

    /**
     * {@code s^2 (J_F' J_F)^-1} over the free columns {@code F}, scattered back
     * into an {@code n x n} matrix with zeros where a parameter is held. Taken
     * from the singular values rather than from the normal equations, which
     * would square the condition number.
     *
     * @return the matrix, or {@code null} if it does not exist
     */
    private static double[] covarianceOf(double[] jacobian, BoundedLeastSquares.Bound[] activeSet, int m, int n,
            int free, double sumOfSquares) {
        int degreesOfFreedom = m - free;
        if (free < 1 || degreesOfFreedom < 1 || !(sumOfSquares >= 0.0)) {
            return null;
        }
        double[] block = new double[m * free];
        int[] index = new int[free];
        int k = 0;
        for (int j = 0; j < n; j++) {
            if (activeSet[j] == BoundedLeastSquares.Bound.FREE) {
                System.arraycopy(jacobian, j * m, block, k * m, m);
                index[k] = j;
                k++;
            }
        }
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(block, m, free);
        if (!svd.converged) {
            return null;
        }
        double tolerance = svd.sigma[0] * Math.max(m, free) * MathConsts.MACH_EPS_DBL;
        for (int i = 0; i < free; i++) {
            if (!(svd.sigma[i] > tolerance)) {
                return null;
            }
        }
        double sigmaSquared = sumOfSquares / degreesOfFreedom;
        double[] covariance = new double[n * n];
        for (int c = 0; c < free; c++) {
            for (int rr = 0; rr <= c; rr++) {
                double sum = 0.0;
                for (int i = 0; i < free; i++) {
                    double d = svd.sigma[i];
                    sum += svd.V[i * free + rr] * svd.V[i * free + c] / (d * d);
                }
                double value = sigmaSquared * sum;
                covariance[index[c] * n + index[rr]] = value;
                covariance[index[rr] * n + index[c]] = value;
            }
        }
        return covariance;
    }

    /**
     * How far the point is from stationary in the directions the box allows. A
     * held parameter may only want to move outward; the equality tests are
     * exact because a held parameter <em>is</em> its bound.
     */
    private static double kktViolation(double[] x, double[] gradient, double[] lower, double[] upper) {
        double worst = 0.0;
        for (int j = 0; j < x.length; j++) {
            worst = Math.max(worst, violationAt(x[j], gradient[j], lower[j], upper[j]));
        }
        return worst;
    }

    /** The same divided by {@code ||J_j|| ||r||}, which is MINPACK's measure. */
    private static double scaledKktViolation(double[] x, double[] gradient, double[] columnNorm, double[] lower,
            double[] upper, double residualNorm) {
        double worst = 0.0;
        for (int j = 0; j < x.length; j++) {
            double denominator = columnNorm[j] * residualNorm;
            if (denominator > 0.0) {
                worst = Math.max(worst, violationAt(x[j], gradient[j], lower[j], upper[j]) / denominator);
            }
        }
        return worst;
    }

    private static double violationAt(double x, double gradient, double lower, double upper) {
        if (x == lower) {
            return Math.max(0.0, -gradient);
        }
        if (x == upper) {
            return Math.max(0.0, gradient);
        }
        return Math.abs(gradient);
    }

    private static void checkArguments(Object function, double[] initial, int residualCount, double[] lower,
            double[] upper) {
        if (function == null) {
            throw new IllegalArgumentException("function is null");
        }
        if (initial == null) {
            throw new IllegalArgumentException("initial is null");
        }
        if (lower == null || upper == null) {
            throw new IllegalArgumentException("a bound array is null");
        }
        int n = initial.length;
        if (n < 1) {
            throw new IllegalArgumentException("initial is empty");
        }
        if (residualCount < n) {
            throw new IllegalArgumentException(
                    "residualCount must be >= initial.length : " + residualCount + " < " + n);
        }
        if (lower.length != n || upper.length != n) {
            throw new IllegalArgumentException("bounds must have length " + n + " : lower " + lower.length
                    + ", upper " + upper.length);
        }
        for (int j = 0; j < n; j++) {
            if (Double.isNaN(initial[j]) || Double.isInfinite(initial[j])) {
                throw new IllegalArgumentException("initial[" + j + "] is not finite : " + initial[j]);
            }
            if (Double.isNaN(lower[j]) || Double.isNaN(upper[j])) {
                throw new IllegalArgumentException("bound " + j + " is NaN");
            }
            if (lower[j] > upper[j]) {
                throw new IllegalArgumentException(
                        "lower > upper at " + j + " : " + lower[j] + " > " + upper[j]);
            }
            if (lower[j] == Double.POSITIVE_INFINITY || upper[j] == Double.NEGATIVE_INFINITY) {
                throw new IllegalArgumentException(
                        "the box is empty at " + j + " : [" + lower[j] + ", " + upper[j] + "]");
            }
        }
    }

    private static double clamp(double value, double lower, double upper) {
        if (value < lower) {
            return lower;
        }
        if (value > upper) {
            return upper;
        }
        return value;
    }

    private static boolean isFinite(double[] values, int length) {
        for (int i = 0; i < length; i++) {
            double d = values[i];
            if (Double.isNaN(d) || Double.isInfinite(d)) {
                return false;
            }
        }
        return true;
    }

    private static double sumOfSquares(double[] values, int length) {
        double sum = 0.0;
        for (int i = 0; i < length; i++) {
            sum += values[i] * values[i];
        }
        return sum;
    }

    /** Overflow-safe Euclidean norm of {@code length} entries from {@code offset}. */
    private static double euclidean(double[] values, int offset, int length) {
        double max = 0.0;
        for (int i = 0; i < length; i++) {
            double abs = Math.abs(values[offset + i]);
            if (abs > max) {
                max = abs;
            }
        }
        if (max == 0.0) {
            return 0.0;
        }
        double sum = 0.0;
        for (int i = 0; i < length; i++) {
            double scaled = values[offset + i] / max;
            sum += scaled * scaled;
        }
        return max * Math.sqrt(sum);
    }

    /** What the loop needs of a model, with or without a derivative of its own. */
    private interface Model {

        void valueAt(double[] x, double[] residuals);

        /**
         * @return how many evaluations of the residuals this cost
         */
        int jacobianAt(double[] x, double[] residuals, double[] jacobian);
    }

    private static final class AnalyticModel implements Model {

        private final DiffDVectorFunction function;

        AnalyticModel(DiffDVectorFunction function) {
            this.function = function;
        }

        @Override
        public void valueAt(double[] x, double[] residuals) {
            function.valueAt(x, residuals);
        }

        @Override
        public int jacobianAt(double[] x, double[] residuals, double[] jacobian) {
            function.jacobianAt(x, jacobian);
            return 0;
        }
    }

    /**
     * Forward differences that never leave the box: the step is taken inwards,
     * and it shrinks rather than crossing where the box is narrow on both
     * sides.
     */
    private static final class DifferenceModel implements Model {

        private final DVectorFunction function;
        private final int m;
        private final double[] lower;
        private final double[] upper;
        private final double relativeStep;
        private final double[] probe;
        private final double[] probed;

        DifferenceModel(DVectorFunction function, int m, int n, double[] lower, double[] upper,
                double functionAccuracy) {
            this.function = function;
            this.m = m;
            this.lower = lower;
            this.upper = upper;
            this.relativeStep = Math.sqrt(Math.max(functionAccuracy, MathConsts.MACH_EPS_DBL));
            this.probe = new double[n];
            this.probed = new double[m];
        }

        @Override
        public void valueAt(double[] x, double[] residuals) {
            function.valueAt(x, residuals);
        }

        @Override
        public int jacobianAt(double[] x, double[] residuals, double[] jacobian) {
            int n = x.length;
            System.arraycopy(x, 0, probe, 0, n);
            for (int j = 0; j < n; j++) {
                double xj = x[j];
                double h = relativeStep * Math.max(Math.abs(xj), 1.0);
                double room = upper[j] - xj;
                if (h <= room) {
                    // forward, the ordinary case
                } else if (h <= xj - lower[j]) {
                    h = -h;
                } else {
                    // the box is narrower than the step on both sides: take the
                    // wider of the two and stay inside it
                    h = (room >= xj - lower[j]) ? room : lower[j] - xj;
                }
                // difference by what the sum actually moved, not by what was
                // asked for: the two part company once h reaches the last bits
                // of xj, and once the box has shrunk it to fit
                double xjPlusH = xj + h;
                h = xjPlusH - xj;
                if (h == 0.0) {
                    // the parameter is pinned, so the column is never used
                    for (int i = 0; i < m; i++) {
                        jacobian[j * m + i] = 0.0;
                    }
                    continue;
                }
                probe[j] = xjPlusH;
                function.valueAt(probe, probed);
                probe[j] = xj;
                int column = j * m;
                for (int i = 0; i < m; i++) {
                    jacobian[column + i] = (probed[i] - residuals[i]) / h;
                }
            }
            return n;
        }
    }
}
