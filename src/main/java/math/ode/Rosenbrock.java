package math.ode;

import math.MathConsts;
import math.fun.DVectorField;
import math.fun.DiffDVectorField;
import math.lapack.LuFactorization;

/**
 * A linearly implicit step for a stiff system, for whichever
 * {@link RosenbrockTableau} it is constructed with.
 * <p>
 * <b>What it is for.</b> An explicit method is stable only while the step size
 * times the largest eigenvalue of the Jacobian stays inside a region a few
 * units across, so on a stiff equation the step is held down by stability and
 * not by accuracy, and no tolerance changes that. A Rosenbrock method is
 * A-stable and L-stable: its step size is limited by the accuracy it is asked
 * for and by nothing else, which is why it can cross a region where the fast
 * modes have already died out in one step instead of a million.
 * <p>
 * <b>What it costs.</b> One Jacobian, one LU factorization and one back
 * substitution per stage, against an explicit method's one field evaluation
 * per stage. That is the trade: an implicit step is expensive and a stiff
 * problem needs very few of them. Where the problem is not stiff the explicit
 * family wins outright, and {@link OdeIntegrator.Result#seemsStiff} is what
 * says which case a caller is in.
 * <p>
 * <b>The Jacobian comes from the field or from differences.</b> Constructed
 * with a {@link DiffDVectorField} it asks; constructed with a plain
 * {@link DVectorField} it takes forward differences, at {@code n + 1} further
 * evaluations per step and an accuracy of about the square root of the machine
 * epsilon. Both are counted in {@link #evaluations()}, so the number stays
 * comparable with what an explicit stepper reports.
 * <p>
 * <b>A retried step costs no second Jacobian.</b> The Jacobian and the partial
 * derivative in time depend on the point and not on the step size, so a step
 * rejected and taken again from the same {@code (t, y)} reuses both. The point
 * is recognized by comparing it, exactly as {@link ExplicitRungeKutta}
 * recognizes a carried stage. Nothing is reused across a point that moved:
 * holding a Jacobian while the solution walks away from it is a heuristic, and
 * this class does not have one.
 * <p>
 * <b>A singular matrix is a step to shrink, not a failure.</b> Where
 * {@code I / (gamma h) - J} comes out exactly singular the step returns a state
 * that is not finite, which is what {@link OdeIntegrator} already does the
 * right thing with -- it halves the step on the adaptive path and says where it
 * gave up on the fixed one.
 * <p>
 * Instances are stateful and cannot be shared between threads.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Rosenbrock_methods">Wikipedia
 * Rosenbrock methods</a>.
 *
 * @since 1.5.3
 */
public final class Rosenbrock implements OdeStepper {

    /**
     * The floor of the forward difference step. See
     * {@link #differenceJacobian(double, double[])} for why it is not the one
     * {@code math.fun} uses.
     */
    private static final double DIFF_FLOOR = 1.0e-5;

    private final RosenbrockTableau tableau;
    private final DVectorField field;
    private final DiffDVectorField analytic;
    private final int n;
    private final int s;

    private final double gamma;
    private final double[] alpha;
    private final double[][] a;
    private final double[][] c;
    private final double[] d;
    private final double[] b;
    private final double[] bError;
    private final double[][] denseWeights;

    private final double[] jacobian;
    private final double[] w;
    private final int[] pivots;
    private final double[][] k;
    private final double[] stageY;
    private final double[] stageF;
    private final double[] dfdt;
    private final double[] rhs;
    private final double[] f0;
    private final double[] differenceY;

    private final double[] cacheY;
    private double cacheT;
    private boolean cacheValid;

    private final double[] startY;
    private final double[] endY;
    private final double[] endDerivative;
    private final double[][] denseTerms;
    private double lastH;
    private double lastT;
    private boolean stepped;
    private boolean endDerivativeValid;
    private boolean denseTermsValid;

    private long evaluations;
    private long jacobians;
    private long factorizations;

    /**
     * A stepper that differences the Jacobian out of the field.
     *
     * @param tableau
     *            the coefficients of the method
     * @param field
     *            the right hand side of the equation
     * @param dimension
     *            the number of components of the state
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the dimension is not
     *             positive
     */
    public Rosenbrock(RosenbrockTableau tableau, DVectorField field, int dimension) {
        this(tableau, field, null, dimension);
    }

    /**
     * A stepper that asks the field for its own Jacobian.
     *
     * @param tableau
     *            the coefficients of the method
     * @param field
     *            the right hand side of the equation together with its two
     *            first derivatives
     * @param dimension
     *            the number of components of the state
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the dimension is not
     *             positive
     */
    public Rosenbrock(RosenbrockTableau tableau, DiffDVectorField field, int dimension) {
        this(tableau, field, field, dimension);
    }

    private Rosenbrock(RosenbrockTableau tableau, DVectorField field, DiffDVectorField analytic,
            int dimension) {
        if (tableau == null) {
            throw new IllegalArgumentException("tableau must not be null");
        }
        if (field == null) {
            throw new IllegalArgumentException("field must not be null");
        }
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be at least 1, not " + dimension);
        }
        this.tableau = tableau;
        this.field = field;
        this.analytic = analytic;
        this.n = dimension;
        this.s = tableau.stages();
        this.gamma = tableau.gamma();
        this.alpha = tableau.alpha();
        this.a = tableau.a();
        this.c = tableau.c();
        this.d = tableau.d();
        this.b = tableau.b();
        this.bError = tableau.bError();
        this.denseWeights = tableau.dense();

        this.jacobian = new double[dimension * dimension];
        this.w = new double[dimension * dimension];
        this.pivots = new int[dimension];
        this.k = new double[s][dimension];
        this.stageY = new double[dimension];
        this.stageF = new double[dimension];
        this.dfdt = new double[dimension];
        this.rhs = new double[dimension];
        this.f0 = new double[dimension];
        this.differenceY = new double[dimension];
        this.cacheY = new double[dimension];
        this.startY = new double[dimension];
        this.endY = new double[dimension];
        this.endDerivative = new double[dimension];
        this.denseTerms = (denseWeights == null) ? null : new double[denseWeights.length][dimension];
    }

    /**
     * The method this stepper executes.
     *
     * @return the tableau given to the constructor
     */
    public RosenbrockTableau tableau() {
        return tableau;
    }

    /**
     * Whether the field supplies its own Jacobian rather than having it
     * differenced out of it.
     *
     * @return {@code true} if the stepper was built on a
     *         {@link DiffDVectorField}
     */
    public boolean hasAnalyticJacobian() {
        return analytic != null;
    }

    /**
     * How often the Jacobian and the derivative in time have been taken since
     * this stepper was created. It is one per step less the steps that reused
     * the point of a rejected one.
     *
     * @return the number of linearizations, never decreasing
     */
    public long jacobians() {
        return jacobians;
    }

    /**
     * How often a matrix has been factored since this stepper was created,
     * which is once per step attempted, accepted or not.
     *
     * @return the number of factorizations, never decreasing
     */
    public long factorizations() {
        return factorizations;
    }

    /**
     * The infinity norm of the Jacobian this stepper is currently holding, the
     * largest of its absolute row sums.
     * <p>
     * <b>It is an upper bound on the largest eigenvalue and not the eigenvalue
     * itself</b>, so <code>|h| jacobianNorm()</code> bounds from above the
     * product an explicit method's stability region has to contain. The bound is
     * tight for a normal Jacobian and can be many decades loose for one that is
     * not, which makes it usable to rule stiffness <em>in</em> and never to rule
     * it out. This is why {@link #stiffnessMeasure()} does not report it:
     * measured against {@link #stiffnessThreshold()}, which is calibrated
     * against the product and not against a bound on it, it would call every
     * stiff equation violently stiff and would carry that into
     * {@link OdeIntegrator.Result#stiffness}.
     * <p>
     * <b>It is as old as the linearization.</b> A step reuses the Jacobian of
     * the point it was taken at ({@link #jacobians()} counts how often that has
     * actually been formed), so between steps this reports the norm at the start
     * of the last one and not at wherever the solution now is.
     *
     * @return <code>||J||_inf</code>, or {@link Double#NaN} if nothing has been
     *         linearized since this stepper was created or last {@link #reset()}
     */
    public double jacobianNorm() {
        if (!cacheValid) {
            return Double.NaN;
        }
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double row = 0.0;
            for (int j = 0; j < n; ++j) {
                row += Math.abs(jacobian[j * n + i]);
            }
            if (row > worst) {
                worst = row;
            }
        }
        return worst;
    }

    @Override
    public int dimension() {
        return n;
    }

    @Override
    public int order() {
        return tableau.order();
    }

    @Override
    public boolean hasErrorEstimate() {
        return tableau.hasErrorEstimate();
    }

    @Override
    public boolean hasDenseOutput() {
        return tableau.hasDenseOutput();
    }

    /**
     * Always {@link Double#NaN}. The estimate exists to tell an explicit method
     * that its step size is being held down by stability rather than by
     * accuracy; this method has no such limit, so there is nothing for it to
     * report and no pair of stages at one time to form it from. The Jacobian it
     * does hold gives a bound rather than the product this is measured by --
     * {@link #jacobianNorm()} says why that is a different reading and not this
     * one.
     *
     * @return {@link Double#NaN}
     */
    @Override
    public double stiffnessMeasure() {
        return Double.NaN;
    }

    @Override
    public void derivative(double t, double[] y, double[] dydt) {
        if (y == null || y.length != n) {
            throw new IllegalArgumentException("y must be of length " + n);
        }
        if (dydt == null || dydt.length != n) {
            throw new IllegalArgumentException("dydt must be of length " + n);
        }
        field.valueAt(t, y, dydt);
        ++evaluations;
    }

    @Override
    public void step(double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (y == null || y.length != n) {
            throw new IllegalArgumentException("y must be of length " + n);
        }
        if (yOut == null || yOut.length != n) {
            throw new IllegalArgumentException("yOut must be of length " + n);
        }
        if (yOut == y) {
            throw new IllegalArgumentException("yOut must not be the array passed as y");
        }
        if (errOut != null && errOut.length != n) {
            throw new IllegalArgumentException("errOut must be of length " + n);
        }
        if (h == 0.0) {
            throw new IllegalArgumentException("h must not be zero");
        }
        stepped = false;
        endDerivativeValid = false;
        denseTermsValid = false;

        linearize(t, y);

        double diagonal = 1.0 / (gamma * h);
        for (int i = 0; i < w.length; ++i) {
            w[i] = -jacobian[i];
        }
        for (int i = 0; i < n; ++i) {
            w[i * n + i] += diagonal;
        }
        ++factorizations;
        if (!LuFactorization.factor(n, w, pivots)) {
            fill(yOut, Double.NaN);
            if (errOut != null) {
                fill(errOut, Double.NaN);
            }
            return;
        }

        for (int i = 0; i < s; ++i) {
            if (i == 0) {
                System.arraycopy(f0, 0, rhs, 0, n);
            } else {
                System.arraycopy(y, 0, stageY, 0, n);
                double[] row = a[i];
                for (int j = 0; j < i; ++j) {
                    double weight = row[j];
                    if (weight != 0.0) {
                        double[] kj = k[j];
                        for (int m = 0; m < n; ++m) {
                            stageY[m] += weight * kj[m];
                        }
                    }
                }
                field.valueAt(t + alpha[i] * h, stageY, stageF);
                ++evaluations;
                System.arraycopy(stageF, 0, rhs, 0, n);
            }
            double[] coupling = c[i];
            for (int j = 0; j < i; ++j) {
                double weight = coupling[j] / h;
                if (weight != 0.0) {
                    double[] kj = k[j];
                    for (int m = 0; m < n; ++m) {
                        rhs[m] += weight * kj[m];
                    }
                }
            }
            double timeWeight = h * d[i];
            if (timeWeight != 0.0) {
                for (int m = 0; m < n; ++m) {
                    rhs[m] += timeWeight * dfdt[m];
                }
            }
            double[] ki = k[i];
            System.arraycopy(rhs, 0, ki, 0, n);
            LuFactorization.solve(n, 1, w, pivots, ki);
        }

        System.arraycopy(y, 0, yOut, 0, n);
        for (int i = 0; i < s; ++i) {
            double weight = b[i];
            if (weight != 0.0) {
                double[] ki = k[i];
                for (int m = 0; m < n; ++m) {
                    yOut[m] += weight * ki[m];
                }
            }
        }
        if (errOut != null && bError != null) {
            fill(errOut, 0.0);
            for (int i = 0; i < s; ++i) {
                double weight = bError[i];
                if (weight != 0.0) {
                    double[] ki = k[i];
                    for (int m = 0; m < n; ++m) {
                        errOut[m] += weight * ki[m];
                    }
                }
            }
        }

        System.arraycopy(y, 0, startY, 0, n);
        System.arraycopy(yOut, 0, endY, 0, n);
        lastH = h;
        lastT = t;
        stepped = true;
    }

    @Override
    public void interpolate(double theta, double[] out) {
        if (!stepped) {
            throw new IllegalStateException("no step has been taken");
        }
        if (out == null || out.length != n) {
            throw new IllegalArgumentException("out must be of length " + n);
        }
        if (theta == 0.0) {
            System.arraycopy(startY, 0, out, 0, n);
            return;
        }
        if (theta == 1.0) {
            System.arraycopy(endY, 0, out, 0, n);
            return;
        }
        if (denseWeights == null) {
            hermite(theta, out);
            return;
        }
        if (!denseTermsValid) {
            for (int j = 0; j < denseWeights.length; ++j) {
                double[] row = denseWeights[j];
                double[] term = denseTerms[j];
                fill(term, 0.0);
                for (int i = 0; i < s; ++i) {
                    double weight = row[i];
                    if (weight != 0.0) {
                        double[] ki = k[i];
                        for (int m = 0; m < n; ++m) {
                            term[m] += weight * ki[m];
                        }
                    }
                }
            }
            denseTermsValid = true;
        }
        double rest = 1.0 - theta;
        for (int m = 0; m < n; ++m) {
            double inner = 0.0;
            for (int j = denseTerms.length - 1; j >= 0; --j) {
                inner = inner * theta + denseTerms[j][m];
            }
            out[m] = rest * startY[m] + theta * (endY[m] + rest * inner);
        }
    }

    @Override
    public long evaluations() {
        return evaluations;
    }

    @Override
    public void reset() {
        stepped = false;
        cacheValid = false;
        endDerivativeValid = false;
        denseTermsValid = false;
    }

    /**
     * The method, whether it was given a Jacobian and the size of the system.
     */
    @Override
    public String toString() {
        return "Rosenbrock[" + tableau + ", " + (hasAnalyticJacobian() ? "analytic" : "differenced")
                + " Jacobian, dimension " + n + "]";
    }

    /**
     * Fills in the field, the Jacobian and the derivative in time at
     * {@code (t, y)}, or recognizes that they are already there.
     */
    private void linearize(double t, double[] y) {
        if (cacheValid && t == cacheT && sameContent(y, cacheY)) {
            return;
        }
        field.valueAt(t, y, f0);
        ++evaluations;
        if (analytic != null) {
            analytic.jacobianAt(t, y, jacobian, dfdt);
        } else {
            differenceJacobian(t, y);
            differenceTimeDerivative(t, y);
        }
        ++jacobians;
        cacheT = t;
        System.arraycopy(y, 0, cacheY, 0, n);
        cacheValid = true;
    }

    /**
     * One column per component, at {@code n} evaluations of the field.
     * <p>
     * <b>The step is Hairer's and not the one {@code math.fun} uses.</b>
     * {@link math.fun.NumericalDiffDVectorFunction} moves a component by
     * {@code sqrt(eps) max(|y|, 1)}, a relative step floored at one, and that
     * floor assumes the arguments are scaled to be of order one. A stiff system
     * routinely is not: Robertson's second component sits at {@code 7e-08}, and
     * a step of {@code 1.5e-08} changes it by twenty percent, which is a
     * difference quotient of something other than the derivative. Hairer's
     * {@code sqrt(eps max(|y|, 1e-5))} agrees with the other rule at a scale of
     * one and follows the component down below it. Measured on Robertson: 314
     * steps against 1652, where the written Jacobian takes 305. On van der Pol,
     * whose components are of order one, the two rules give bit for bit the
     * same run.
     */
    private void differenceJacobian(double t, double[] y) {
        System.arraycopy(y, 0, differenceY, 0, n);
        for (int j = 0; j < n; ++j) {
            double base = differenceY[j];
            double moved = base
                    + Math.sqrt(MathConsts.MACH_EPS_DBL * Math.max(Math.abs(base), DIFF_FLOOR));
            double step = moved - base;
            differenceY[j] = moved;
            field.valueAt(t, differenceY, stageF);
            ++evaluations;
            int column = j * n;
            for (int i = 0; i < n; ++i) {
                jacobian[column + i] = (stageF[i] - f0[i]) / step;
            }
            differenceY[j] = base;
        }
    }

    /** One evaluation, at the same step rule applied to the time. */
    private void differenceTimeDerivative(double t, double[] y) {
        double moved = t + Math.sqrt(MathConsts.MACH_EPS_DBL * Math.max(Math.abs(t), DIFF_FLOOR));
        double step = moved - t;
        if (step == 0.0) {
            fill(dfdt, 0.0);
            return;
        }
        field.valueAt(moved, y, stageF);
        ++evaluations;
        for (int i = 0; i < n; ++i) {
            dfdt[i] = (stageF[i] - f0[i]) / step;
        }
    }

    /**
     * The cubic through the two states and the two derivatives, for a method
     * with no continuous extension of its own. The derivative at the start is
     * the field evaluation the step began with and is free; the one at the end
     * costs an evaluation and is taken only if somebody asks.
     */
    private void hermite(double theta, double[] out) {
        if (!endDerivativeValid) {
            field.valueAt(lastT + lastH, endY, endDerivative);
            ++evaluations;
            endDerivativeValid = true;
        }
        double square = theta * theta;
        double cube = square * theta;
        double h00 = 2.0 * cube - 3.0 * square + 1.0;
        double h10 = cube - 2.0 * square + theta;
        double h01 = -2.0 * cube + 3.0 * square;
        double h11 = cube - square;
        for (int m = 0; m < n; ++m) {
            out[m] = h00 * startY[m] + h01 * endY[m] + lastH * (h10 * f0[m] + h11 * endDerivative[m]);
        }
    }

    private static void fill(double[] target, double value) {
        for (int i = 0; i < target.length; ++i) {
            target[i] = value;
        }
    }

    private static boolean sameContent(double[] x, double[] y) {
        for (int i = 0; i < x.length; ++i) {
            if (x[i] != y[i]) {
                return false;
            }
        }
        return true;
    }
}
