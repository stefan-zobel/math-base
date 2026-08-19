package math.optim;

import math.MathConsts;
import math.fun.DVectorFunction;
import math.fun.DiffDVectorFunction;
import math.minpack.Lmder_fcn;
import math.minpack.Lmdif_fcn;
import math.minpack.Minpack_f77;

/**
 * Nonlinear least squares by the Levenberg-Marquardt method: given {@code m}
 * residuals over {@code n} parameters, find the parameters that minimize the
 * sum of their squares.
 * <p>
 * This is a facade over {@link Minpack_f77}, the Java translation of the
 * FORTRAN MINPACK of Garbow, Hillstrom and More, which does the actual work.
 * That translation kept the calling convention of its origin -- one-based
 * arrays, nineteen positional arguments, results returned through
 * {@code int[1]} boxes -- so the algorithm was there but only reachable by a
 * caller willing to read FORTRAN. Nothing of it is changed here; the arguments
 * are named, the arrays are zero-based, and the status is an enum rather than
 * an integer whose ordering means nothing.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm">
 *      Levenberg-Marquardt algorithm</a>
 * @since 1.5.2
 */
public final class LevenbergMarquardt {

    /**
     * Why the search stopped. The integer MINPACK reports is not ordered by
     * quality -- {@link #GRADIENT_AT_MACHINE_PRECISION} is the largest of them
     * and one of the best outcomes -- which is the main reason this enum
     * exists.
     *
     * @since 1.5.2
     */
    public enum Status {

        /** The arguments were rejected before any work was done. */
        INVALID_INPUT(0, false),
        /** The relative error in the sum of squares is at most {@code ftol}. */
        VALUE_TOLERANCE_REACHED(1, true),
        /** The relative error in the parameters is at most {@code xtol}. */
        STEP_TOLERANCE_REACHED(2, true),
        /** Both of the above hold. */
        VALUE_AND_STEP_TOLERANCE_REACHED(3, true),
        /** The scaled gradient norm is at most {@code gtol}. */
        GRADIENT_TOLERANCE_REACHED(4, true),
        /**
         * The evaluation budget ran out. This says nothing about the point that
         * comes back with it: a problem whose Jacobian is rank deficient at the
         * minimum can reach it early and then spend a hundred times as many
         * evaluations failing to satisfy the step test, so the parameters are
         * worth looking at before the run is discarded.
         */
        TOO_MANY_EVALUATIONS(5, false),
        /**
         * The sum of squares cannot be reduced any further in double precision,
         * so the requested {@code ftol} can never be met. Often that means the
         * minimum was reached to full accuracy -- but the trust region can also
         * collapse short of one on a badly modelled problem, and MINPACK does
         * not distinguish the two. Not reported as a success for that reason;
         * check the gradient if it matters.
         */
        VALUE_TOLERANCE_TOO_SMALL(6, false),
        /**
         * The parameters cannot be improved any further in double precision, so
         * the requested {@code xtol} can never be met. The same caveat as
         * above applies.
         */
        STEP_TOLERANCE_TOO_SMALL(7, false),
        /**
         * The scaled gradient norm is down at machine precision. Unlike the two
         * above, this is a statement about the gradient rather than about the
         * progress being made, so it does identify a stationary point, and it
         * is a stronger one than any {@code gtol} could ask for. MINPACK files
         * it under "gtol is too small" and its own {@code lmder1} driver
         * rewrites it to {@link #GRADIENT_TOLERANCE_REACHED} before returning.
         */
        GRADIENT_AT_MACHINE_PRECISION(8, true);

        private final int code;
        private final boolean success;

        private Status(int code, boolean success) {
            this.code = code;
            this.success = success;
        }

        /**
         * The raw MINPACK {@code info} value, for cross-referencing the FORTRAN
         * documentation.
         *
         * @return the {@code info} value this status was translated from
         */
        public int code() {
            return code;
        }

        /**
         * Whether the returned parameters are a minimum, either to the accuracy
         * that was asked for or to the accuracy the arithmetic allows.
         *
         * @return {@code true} if the search succeeded
         */
        public boolean isSuccess() {
            return success;
        }
    }

    /**
     * The outcome of a fit.
     *
     * @since 1.5.2
     */
    public static final class Result {

        /** The parameters found, of length {@code n}. */
        public final double[] parameters;
        /** The residuals at {@link #parameters}, of length {@code m}. */
        public final double[] residuals;
        /** The sum of the squared {@link #residuals}, the quantity minimized. */
        public final double sumOfSquares;
        /** How often the residuals were evaluated. */
        public final int functionEvaluations;
        /** How often the Jacobian was evaluated. */
        public final int jacobianEvaluations;
        /** Why the search stopped. */
        public final Status status;
        /** Shorthand for {@code status.isSuccess()}. */
        public final boolean converged;
        /** {@code m - n}, the degrees of freedom left for estimating the noise. */
        public final int degreesOfFreedom;
        /**
         * The estimated variance-covariance matrix of {@link #parameters},
         * {@code n} by {@code n}, flat and column-major, or {@code null} where
         * it cannot be formed -- see below.
         * <p>
         * It is {@code s^2 (J'J)^-1} with {@code s^2 = sumOfSquares /
         * degreesOfFreedom}, taken from the factorization the search already
         * produced. <strong>Two things have to hold before it means
         * anything.</strong> The residuals must be independent draws of equal
         * variance, or it estimates the wrong quantity; and for a model that is
         * nonlinear in its parameters -- which is the only reason to be using
         * this class -- it is the linearization at the solution, exact only in
         * the limit of many observations. Close to the minimum of a
         * well-conditioned problem it is a good approximation and it is what
         * every curve fitting package reports; on a badly conditioned or
         * strongly curved one it can understate the true uncertainty by a lot.
         * A bootstrap over the residuals answers the same question without
         * either assumption, and {@code math.probe.Bootstrap} is in this
         * library.
         * <p>
         * {@code null} when {@code m == n}, where there is no spread left to
         * estimate a variance from, and when the Jacobian at the solution is
         * numerically rank deficient, where the inverse does not exist. Neither
         * is reported as an error, because neither says the fit failed -- the
         * parameters are as good as they ever were.
         */
        public final double[] covariance;
        /**
         * The square roots of the diagonal of {@link #covariance}, of length
         * {@code n}, or {@code null} whenever that is. The same caveats apply.
         */
        public final double[] standardErrors;

        Result(double[] parameters, double[] residuals, double sumOfSquares, int functionEvaluations,
                int jacobianEvaluations, Status status, int degreesOfFreedom, double[] covariance) {
            this.parameters = parameters;
            this.residuals = residuals;
            this.sumOfSquares = sumOfSquares;
            this.functionEvaluations = functionEvaluations;
            this.jacobianEvaluations = jacobianEvaluations;
            this.status = status;
            this.converged = status.isSuccess();
            this.degreesOfFreedom = degreesOfFreedom;
            this.covariance = covariance;
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
    }

    /**
     * Default for both step and value tolerance. MINPACK recommends nothing
     * below the square root of the machine epsilon, because the sum of squares
     * is flat at its minimum and half the digits are gone by the time it is
     * reached.
     */
    private static final double DEFAULT_TOLERANCE = Math.sqrt(MathConsts.MACH_EPS_DBL);

    /** Default bound on the first step, MINPACK's {@code factor}. */
    private static final double DEFAULT_INITIAL_STEP_BOUND = 100.0;

    /** Stands for the budget of {@code 100 * (n + 1)} that MINPACK defaults to. */
    private static final int EVALUATIONS_FROM_PROBLEM_SIZE = 0;

    /** Residuals are assumed accurate to machine precision unless said otherwise. */
    private static final double DEFAULT_FUNCTION_ACCURACY = 0.0;

    private final double valueTolerance;
    private final double stepTolerance;
    private final double gradientTolerance;
    private final int maxEvaluations;
    private final double initialStepBound;
    private final double functionAccuracy;

    /**
     * Creates a solver with the default stopping rules: both the value and the
     * step tolerance at the square root of the machine epsilon, no gradient
     * tolerance, and a budget of {@code 100 * (n + 1)} evaluations.
     */
    public LevenbergMarquardt() {
        this(DEFAULT_TOLERANCE, DEFAULT_TOLERANCE, 0.0, EVALUATIONS_FROM_PROBLEM_SIZE,
                DEFAULT_INITIAL_STEP_BOUND);
    }

    /**
     * Creates a solver with explicit stopping rules. The search stops as soon
     * as any one of the three tolerances is met, so a tolerance of {@code 0}
     * switches that test off rather than demanding the impossible.
     *
     * @param valueTolerance
     *            relative error in the sum of squares that is accepted,
     *            {@code 0} or greater. Asking for less than the square root of
     *            the machine epsilon cannot be met and ends the run with
     *            {@link Status#VALUE_TOLERANCE_TOO_SMALL}
     * @param stepTolerance
     *            relative error in the parameters that is accepted, {@code 0}
     *            or greater, with the same caveat
     * @param gradientTolerance
     *            scaled gradient norm that is accepted, {@code 0} or greater.
     *            {@code 0} switches the test off, which is what the MINPACK
     *            drivers themselves do
     * @param maxEvaluations
     *            budget for evaluations of the residuals, {@code 1} or greater
     * @param initialStepBound
     *            bounds the first step to this multiple of the scaled length of
     *            the starting point, greater than {@code 0}. MINPACK calls it
     *            {@code factor} and recommends a value in {@code [0.1, 100]};
     *            lowering it helps when the first step overshoots into a region
     *            where the residuals cannot be evaluated
     */
    public LevenbergMarquardt(double valueTolerance, double stepTolerance, double gradientTolerance,
            int maxEvaluations, double initialStepBound) {
        this(valueTolerance, stepTolerance, gradientTolerance, maxEvaluations, initialStepBound,
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
     * @param gradientTolerance
     *            as above
     * @param maxEvaluations
     *            as above
     * @param initialStepBound
     *            as above
     * @param functionAccuracy
     *            the relative error with which {@code valueAt} computes the
     *            residuals, {@code 0} or greater, {@code 0} meaning machine
     *            precision. Ignored by
     *            {@link #solve(DiffDVectorFunction, double[], int)}, which is
     *            handed the derivative and never has to guess a step for it.
     *            {@link #solve(DVectorFunction, double[], int)} does have to,
     *            and it picks that step from this number: a residual that comes
     *            out of a simulation, a quadrature or a table is accurate to
     *            far less than machine precision, and differencing it over the
     *            default step measures the error rather than the derivative.
     *            MINPACK calls it {@code epsfcn}
     */
    public LevenbergMarquardt(double valueTolerance, double stepTolerance, double gradientTolerance,
            int maxEvaluations, double initialStepBound, double functionAccuracy) {
        if (!(functionAccuracy >= 0.0) || Double.isInfinite(functionAccuracy)) {
            throw new IllegalArgumentException(
                    "functionAccuracy must be finite and non-negative : " + functionAccuracy);
        }
        if (!(valueTolerance >= 0.0) || Double.isInfinite(valueTolerance)) {
            throw new IllegalArgumentException("valueTolerance must be finite and non-negative : " + valueTolerance);
        }
        if (!(stepTolerance >= 0.0) || Double.isInfinite(stepTolerance)) {
            throw new IllegalArgumentException("stepTolerance must be finite and non-negative : " + stepTolerance);
        }
        if (!(gradientTolerance >= 0.0) || Double.isInfinite(gradientTolerance)) {
            throw new IllegalArgumentException(
                    "gradientTolerance must be finite and non-negative : " + gradientTolerance);
        }
        if (maxEvaluations < 1 && maxEvaluations != EVALUATIONS_FROM_PROBLEM_SIZE) {
            throw new IllegalArgumentException("maxEvaluations must be >= 1 : " + maxEvaluations);
        }
        if (!(initialStepBound > 0.0) || Double.isInfinite(initialStepBound)) {
            throw new IllegalArgumentException(
                    "initialStepBound must be finite and positive : " + initialStepBound);
        }
        this.valueTolerance = valueTolerance;
        this.stepTolerance = stepTolerance;
        this.gradientTolerance = gradientTolerance;
        this.maxEvaluations = maxEvaluations;
        this.initialStepBound = initialStepBound;
        this.functionAccuracy = functionAccuracy;
    }

    /**
     * Minimizes the sum of the squared residuals of {@code function}, starting
     * from {@code initial}.
     * <p>
     * The starting point is not modified. The parameters are scaled internally
     * from the columns of the Jacobian, so a problem whose parameters differ by
     * orders of magnitude does not need to be rescaled by hand.
     *
     * @param function
     *            the residuals and their Jacobian
     * @param initial
     *            the starting parameters, of length {@code n}, all finite (not
     *            modified)
     * @param residualCount
     *            the number {@code m} of residuals {@code function} produces,
     *            which must be at least {@code n}
     * @return the parameters found, the residuals there, and why the search
     *         stopped
     * @throws IllegalArgumentException
     *             if an argument is null, empty, not finite, or if
     *             {@code residualCount} is smaller than {@code initial.length}
     */
    public Result solve(DiffDVectorFunction function, double[] initial, int residualCount) {
        checkArguments(function, initial, residualCount);
        int n = initial.length;
        int m = residualCount;

        // MINPACK indexes from one, so every array is one longer than it needs
        // to be and element zero stays unused
        double[] x = new double[n + 1];
        System.arraycopy(initial, 0, x, 1, n);
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        // mode 1 is MINPACK's internal scaling, nprint 0 turns off its tracing
        Minpack_f77.lmder_f77(new JacobianAdapter(function, m, n), m, n, x, fvec, fjac, valueTolerance,
                stepTolerance, gradientTolerance, budgetFor(n), diag, 1, initialStepBound, 0, info, nfev, njev,
                ipvt, qtf);

        return unpack(x, fvec, fjac, ipvt, n, m, info[1], nfev[1], njev[1]);
    }

    /**
     * Minimizes the sum of the squared residuals of {@code function} without
     * being given a derivative, approximating the Jacobian by forward
     * differences instead.
     * <p>
     * This costs both accuracy and evaluations -- on the collection of More,
     * Garbow and Hillstrom it reaches the same minima, but never in fewer
     * evaluations and sometimes in eighty times as many. Supply a derivative
     * through {@link #solve(DiffDVectorFunction, double[], int)} where one can
     * be had. Note that Java picks between the two by the <em>static</em> type
     * of the argument, so a {@link DiffDVectorFunction} held in a variable
     * declared {@link DVectorFunction} lands here and its derivative goes
     * unused.
     * <p>
     * If the residuals are not computed to machine precision, say so through
     * the {@code functionAccuracy} of the six-argument constructor. The step
     * for the difference is chosen from it, and left at the default it will be
     * far too small to see past the error in the values. That failure is
     * silent and it reports success: with a step below the resolution of the
     * residuals every difference comes out zero, the approximated Jacobian is
     * the zero matrix, and a zero matrix is orthogonal to every residual -- so
     * MINPACK's gradient test fires and the run ends with
     * {@link Status#GRADIENT_TOLERANCE_REACHED} on the starting point, having
     * moved nothing. There is no way to tell that from a genuine minimum after
     * the fact, which is why it has to be declared beforehand.
     *
     * @param function
     *            the residuals
     * @param initial
     *            the starting parameters, of length {@code n}, all finite (not
     *            modified)
     * @param residualCount
     *            the number {@code m} of residuals {@code function} produces,
     *            which must be at least {@code n}
     * @return the parameters found, the residuals there, and why the search
     *         stopped. {@link Result#jacobianEvaluations} is {@code 0}: the
     *         Jacobian is never evaluated, only approximated, and the cost of
     *         that shows up in {@link Result#functionEvaluations}
     * @throws IllegalArgumentException
     *             if an argument is null, empty, not finite, or if
     *             {@code residualCount} is smaller than {@code initial.length}
     */
    public Result solve(DVectorFunction function, double[] initial, int residualCount) {
        checkArguments(function, initial, residualCount);
        int n = initial.length;
        int m = residualCount;

        double[] x = new double[n + 1];
        System.arraycopy(initial, 0, x, 1, n);
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];

        Minpack_f77.lmdif_f77(new ValueAdapter(function, m, n), m, n, x, fvec, valueTolerance, stepTolerance,
                gradientTolerance, budgetFor(n), functionAccuracy, diag, 1, initialStepBound, 0, info, nfev, fjac,
                ipvt, qtf);

        return unpack(x, fvec, fjac, ipvt, n, m, info[1], nfev[1], 0);
    }

    private static void checkArguments(DVectorFunction function, double[] initial, int residualCount) {
        if (function == null) {
            throw new IllegalArgumentException("function is null");
        }
        if (initial == null) {
            throw new IllegalArgumentException("initial is null");
        }
        int n = initial.length;
        if (n < 1) {
            throw new IllegalArgumentException("initial is empty");
        }
        if (residualCount < n) {
            throw new IllegalArgumentException(
                    "residualCount must be >= initial.length : " + residualCount + " < " + n);
        }
        for (int j = 0; j < n; j++) {
            if (Double.isNaN(initial[j]) || Double.isInfinite(initial[j])) {
                throw new IllegalArgumentException("initial[" + j + "] is not finite : " + initial[j]);
            }
        }
    }

    private int budgetFor(int n) {
        return (maxEvaluations == EVALUATIONS_FROM_PROBLEM_SIZE) ? 100 * (n + 1) : maxEvaluations;
    }

    private static Result unpack(double[] x, double[] fvec, double[][] fjac, int[] ipvt, int n, int m, int info,
            int nfev, int njev) {
        double[] parameters = new double[n];
        System.arraycopy(x, 1, parameters, 0, n);
        double[] residuals = new double[m];
        System.arraycopy(fvec, 1, residuals, 0, m);
        double sumOfSquares = 0.0;
        for (int i = 0; i < m; i++) {
            sumOfSquares += residuals[i] * residuals[i];
        }
        return new Result(parameters, residuals, sumOfSquares, nfev, njev, statusOf(info), m - n,
                covarianceOf(fjac, ipvt, m, n, sumOfSquares));
    }

    /**
     * The variance-covariance matrix, from what the search left behind rather
     * than from anything recomputed. Both drivers return the upper triangular
     * {@code R} of the pivoted {@code QR} of the final Jacobian in {@code fjac}
     * and the pivoting in {@code ipvt}, with
     * {@code P' (J'J) P = R'R}, so
     * {@code (J'J)^-1 = P (R'R)^-1 P' = P R^-1 R^-T P'}.
     * <p>
     * MINPACK ships a {@code covar} routine that does this in place; it was not
     * part of the translation, and inverting the triangle by plain back
     * substitution is easier to check than reproducing that packing would be.
     *
     * @return the matrix, or {@code null} if it does not exist
     */
    private static double[] covarianceOf(double[][] fjac, int[] ipvt, int m, int n, double sumOfSquares) {
        int degreesOfFreedom = m - n;
        if (degreesOfFreedom < 1) {
            return null;
        }
        // the pivoting orders the diagonal of R by nonincreasing magnitude, so
        // the first entry is the largest and sets the scale for the rank test
        double tolerance = Math.abs(fjac[1][1]) * Math.max(m, n) * MathConsts.MACH_EPS_DBL;
        for (int k = 0; k < n; k++) {
            if (!(Math.abs(fjac[k + 1][k + 1]) > tolerance)) {
                return null;
            }
        }

        // invert the upper triangle by back substitution, column by column
        double[] inverse = new double[n * n];
        for (int j = 0; j < n; j++) {
            inverse[j * n + j] = 1.0 / fjac[j + 1][j + 1];
            for (int i = j - 1; i >= 0; i--) {
                double s = 0.0;
                for (int k = i + 1; k <= j; k++) {
                    s += fjac[i + 1][k + 1] * inverse[j * n + k];
                }
                inverse[j * n + i] = -s / fjac[i + 1][i + 1];
            }
        }

        double sigmaSquared = sumOfSquares / degreesOfFreedom;
        double[] covariance = new double[n * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i <= j; i++) {
                double s = 0.0;
                // the inverse is upper triangular, so both factors are non-zero
                // only from row j onwards
                for (int k = j; k < n; k++) {
                    s += inverse[k * n + i] * inverse[k * n + j];
                }
                double value = sigmaSquared * s;
                // undo the pivoting: position i of the factorization belongs to
                // the parameter ipvt[i] the caller passed in
                int a = ipvt[i + 1] - 1;
                int b = ipvt[j + 1] - 1;
                covariance[b * n + a] = value;
                covariance[a * n + b] = value;
            }
        }
        return covariance;
    }

    private static Status statusOf(int info) {
        Status[] all = Status.values();
        for (int k = 0; k < all.length; k++) {
            if (all[k].code() == info) {
                return all[k];
            }
        }
        // Only reachable if the residual function asked MINPACK to stop by
        // setting iflag negative, which the adapter below never does
        throw new IllegalStateException("unexpected MINPACK info value : " + info);
    }

    /**
     * Turns a {@link DiffDVectorFunction} into the callback MINPACK expects:
     * one-based arrays, a row-major two-dimensional Jacobian, and the choice
     * between value and derivative made by a flag rather than by which method
     * is called.
     */
    private static final class JacobianAdapter implements Lmder_fcn {

        private final DiffDVectorFunction function;
        private final double[] x;
        private final double[] values;
        private final double[] jacobian;

        JacobianAdapter(DiffDVectorFunction function, int m, int n) {
            this.function = function;
            this.x = new double[n];
            this.values = new double[m];
            this.jacobian = new double[m * n];
        }

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, double[][] fjac, int[] iflag) {
            System.arraycopy(xOneBased, 1, x, 0, n);
            if (iflag[1] == 1) {
                function.valueAt(x, values);
                System.arraycopy(values, 0, fvec, 1, m);
            } else {
                // column-major in, row-major out; see math.fun.DJacobian
                function.jacobianAt(x, jacobian);
                for (int j = 0; j < n; j++) {
                    int column = j * m;
                    for (int i = 0; i < m; i++) {
                        fjac[i + 1][j + 1] = jacobian[column + i];
                    }
                }
            }
        }
    }

    /**
     * The same for the derivative-free path, where MINPACK asks only for
     * values and builds the Jacobian from them itself.
     */
    private static final class ValueAdapter implements Lmdif_fcn {

        private final DVectorFunction function;
        private final double[] x;
        private final double[] values;

        ValueAdapter(DVectorFunction function, int m, int n) {
            this.function = function;
            this.x = new double[n];
            this.values = new double[m];
        }

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, int[] iflag) {
            System.arraycopy(xOneBased, 1, x, 0, n);
            function.valueAt(x, values);
            System.arraycopy(values, 0, fvec, 1, m);
        }
    }
}
