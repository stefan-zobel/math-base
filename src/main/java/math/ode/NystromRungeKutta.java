package math.ode;

import math.fun.DSecondOrderField;

/**
 * An explicit Runge-Kutta-Nystroem method with an adaptive step size, for the
 * second order form <code>q'' = f(t, q)</code>.
 * <p>
 * <b>This is the accurate answer over a short horizon, where
 * {@link SymplecticNystrom} is the stable one over a long horizon.</b> A
 * splitting method keeps its energy error inside a band, but only at a constant
 * step size, so it cannot also be driven by an error estimate. This carries an
 * embedded solution of a lower order instead, is judged by
 * {@link StepController} like any other adaptive method, and makes no claim
 * about anything a long run conserves.
 * <p>
 * <b>The equation class is not a mechanical one.</b> A semi-discretized wave
 * equation, a Sturm-Liouville problem and a shooting reduction of a boundary
 * value problem all have this shape and no masses in them. What the shape buys
 * is stages: because the acceleration does not read the velocity, order
 * conditions drop out and a given order needs fewer evaluations than the same
 * order on the flattened system {@code y' = f(t, y)}.
 * <p>
 * <b>What that is worth, measured</b> as the evaluations needed to reach a
 * given error at the endpoint -- matched on the accuracy achieved rather than on
 * the tolerance asked for, since the three methods spend a tolerance
 * differently:
 * <table border="1">
 * <caption>this against the same equation flattened, as a fraction of what the
 * flattened run costs</caption>
 * <tr><th>problem</th><th>vs Dormand-Prince 5(4)</th><th>vs DOP853</th></tr>
 * <tr><td>harmonic oscillator</td><td>0.38 to 0.19</td><td>0.98 to 1.46</td></tr>
 * <tr><td>circular Kepler orbit, ten turns</td><td>0.26 to 0.20</td>
 * <td>1.15 to 2.81</td></tr>
 * <tr><td>semi-discretized wave, 24 nodes</td><td>0.44 to 0.26</td>
 * <td>0.49 to 1.45</td></tr>
 * </table>
 * <p>
 * <b>So it is a clear win over the fifth order method and a crossing against the
 * eighth.</b> Against Dormand-Prince it costs a quarter to a half, and the
 * advantage grows as the tolerance tightens because the orders differ. Against
 * {@link ButcherTableau#DOP853} it wins at a loose tolerance and loses below
 * about {@code 1e-8}, which is the same crossing the first order pair already
 * has, moved by two orders. For an answer wanted to twelve digits, flattening
 * and calling {@link Ode#solveAccurate} is the cheaper route.
 * <p>
 * <b>That is also a requirement and not a preference.</b> The order holds
 * <em>because</em> the acceleration ignores the velocity it is handed. A field
 * that reads it silently yields a method of lower order whose error estimate
 * understates the damage in the same proportion, so nothing in the run reveals
 * it. This stepper therefore probes for that once, on its first step, and
 * refuses to run rather than return a plausible wrong answer. See
 * {@link #step(double, double[], double, double[], double[])} for what the
 * probe can and cannot see. {@link SymplecticNystrom} accepts such a field and
 * loses only its symplecticity; {@link Ode#solve} accepts it outright.
 * <p>
 * The state is the position and the velocity stacked, position first, so a
 * system of {@code n} coordinates is advanced as a state of {@code 2n}. The
 * error estimate is filled over the whole of it, which lets the controller scale
 * the two halves by their own tolerances.
 * <p>
 * Instances are stateful and cannot be shared between threads.
 * <p>
 * <b>See</b> <a href=
 * "https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Nystr%C3%B6m_methods">
 * Wikipedia Nystroem methods</a>.
 *
 * @since 1.5.3
 */
public final class NystromRungeKutta implements OdeStepper {

    /**
     * How far two accelerations at the same point may differ before the field
     * counts as reading the velocity. Scaled by the size of the acceleration,
     * so it is a relative test; a field that genuinely ignores the velocity
     * returns bit-identical values and clears this by every decade there is.
     */
    private static final double VELOCITY_PROBE_TOLERANCE = 1.0e-12;

    private final NystromTableau tableau;
    private final DSecondOrderField field;
    private final int dof;
    private final int dimension;
    private final int stages;

    private final double[] c;
    private final double[][] a;
    private final double[] bbar;
    private final double[] b;
    private final double[] bbarStar;
    private final double[] bStar;
    private final double[][] densePosition;
    private final double[][] denseVelocity;
    private final boolean fsal;

    private final double[][] g;
    private final double[] argQ;
    private final double[] argV;
    private final double[] probe;
    private final double[] startQ;
    private final double[] startV;
    private final double[] endQ;
    private final double[] endV;

    private double lastH;
    private boolean stepped;
    private boolean velocityChecked;
    private boolean stageZeroReady;
    private boolean carryValid;
    private double carryT;
    private final double[] carryQ;
    private final double[] carryG;
    private long evaluations;

    /**
     * A stepper for the given method and equation.
     *
     * @param tableau
     *            the coefficients of the method, which must estimate the error
     *            of its own steps
     * @param field
     *            the acceleration, which must not read the velocity it is
     *            handed
     * @param degreesOfFreedom
     *            the number of coordinates of the position, half the dimension
     *            of the state this stepper advances
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, the count is not positive, or
     *             the method carries no embedded solution
     */
    public NystromRungeKutta(NystromTableau tableau, DSecondOrderField field,
            int degreesOfFreedom) {
        if (tableau == null) {
            throw new IllegalArgumentException("tableau must not be null");
        }
        if (field == null) {
            throw new IllegalArgumentException("field must not be null");
        }
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException("degreesOfFreedom must be at least 1, not "
                    + degreesOfFreedom);
        }
        if (!tableau.hasErrorEstimate()) {
            throw new IllegalArgumentException(tableau
                    + " estimates no error of its own, so its steps cannot be judged and its step"
                    + " size cannot be adapted");
        }
        this.tableau = tableau;
        this.field = field;
        this.dof = degreesOfFreedom;
        this.dimension = 2 * degreesOfFreedom;
        this.stages = tableau.stages();
        this.c = tableau.c();
        this.a = tableau.a();
        this.bbar = tableau.bbar();
        this.b = tableau.b();
        this.bbarStar = tableau.bbarStar();
        this.bStar = tableau.bStar();
        this.densePosition = tableau.densePosition();
        this.denseVelocity = tableau.denseVelocity();
        this.fsal = tableau.isFsal();
        this.g = new double[stages][degreesOfFreedom];
        this.argQ = new double[degreesOfFreedom];
        this.argV = new double[degreesOfFreedom];
        this.probe = new double[degreesOfFreedom];
        this.startQ = new double[degreesOfFreedom];
        this.startV = new double[degreesOfFreedom];
        this.endQ = new double[degreesOfFreedom];
        this.endV = new double[degreesOfFreedom];
        this.carryQ = new double[degreesOfFreedom];
        this.carryG = new double[degreesOfFreedom];
    }

    /**
     * The method this stepper executes.
     *
     * @return the tableau given to the constructor
     */
    public NystromTableau tableau() {
        return tableau;
    }

    /**
     * The number of coordinates of the position, which is half
     * {@link #dimension()}.
     *
     * @return the degrees of freedom
     */
    public int degreesOfFreedom() {
        return dof;
    }

    @Override
    public int dimension() {
        return dimension;
    }

    @Override
    public int order() {
        return tableau.order();
    }

    @Override
    public boolean hasErrorEstimate() {
        return true;
    }

    /**
     * Whether the method carries a continuous extension of its own order rather
     * than leaving the interpolation to a Hermite polynomial over the state.
     * <p>
     * With one, the state inside a step comes back at the order of the method:
     * measured on a Kepler orbit the interpolation error falls as
     * <code>h^6</code> against the Hermite polynomial's <code>h^4</code>, which
     * at a step of {@code 0.5} is seventy times more accurate. Without one the
     * Hermite fallback is still free here, since the velocity is part of the
     * state and <code>q' = v</code> holds exactly.
     */
    @Override
    public boolean hasDenseOutput() {
        return densePosition != null;
    }

    /**
     * Always {@link Double#NaN}. The estimate is a difference quotient between
     * two stages evaluated at the same time, and no two stages here are; a
     * number invented for this would go into
     * {@link OdeIntegrator.Result#stiffness} and mean nothing there.
     */
    @Override
    public double stiffnessMeasure() {
        return Double.NaN;
    }

    @Override
    public void derivative(double t, double[] y, double[] dydt) {
        if (y == null || y.length != dimension) {
            throw new IllegalArgumentException("y must be of length " + dimension);
        }
        if (dydt == null || dydt.length != dimension) {
            throw new IllegalArgumentException("dydt must be of length " + dimension);
        }
        System.arraycopy(y, 0, argQ, 0, dof);
        System.arraycopy(y, dof, argV, 0, dof);
        System.arraycopy(argV, 0, dydt, 0, dof);
        field.valueAt(t, argQ, argV, probe);
        ++evaluations;
        System.arraycopy(probe, 0, dydt, dof, dof);
    }

    /**
     * Advances the state by one step and estimates the error of doing so.
     * <p>
     * <b>The first step also asks whether the field reads the velocity</b>, by
     * evaluating it twice at the same time and position with two different
     * velocities. One of those two is the first stage of the step, so the
     * question costs a single extra evaluation, once per {@link #reset()}, and
     * turns a silent loss of order into an exception. It is not a proof: a
     * dependence that happens to vanish at the point probed is not seen. It
     * catches the case that matters, which is a field that reads the velocity
     * at all.
     */
    @Override
    public void step(double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (y == null || y.length != dimension) {
            throw new IllegalArgumentException("y must be of length " + dimension);
        }
        if (yOut == null || yOut.length != dimension) {
            throw new IllegalArgumentException("yOut must be of length " + dimension);
        }
        System.arraycopy(y, 0, startQ, 0, dof);
        System.arraycopy(y, dof, startV, 0, dof);

        if (!velocityChecked) {
            checkTheFieldIgnoresTheVelocity(t);
            velocityChecked = true;
        }

        for (int i = 0; i < stages; ++i) {
            if (i == 0 && stageZeroReady) {
                // the velocity probe already evaluated at this very point
                stageZeroReady = false;
                continue;
            }
            if (i == 0 && fsal && carryValid && carryT == t && sameContent(startQ, carryQ)) {
                // the last stage of the previous step sat at this very point
                System.arraycopy(carryG, 0, g[0], 0, dof);
                continue;
            }
            for (int m = 0; m < dof; ++m) {
                double acc = startQ[m] + c[i] * h * startV[m];
                double[] row = a[i];
                for (int j = 0; j < i; ++j) {
                    acc += h * h * row[j] * g[j][m];
                }
                argQ[m] = acc;
                argV[m] = startV[m];
            }
            field.valueAt(t + c[i] * h, argQ, argV, g[i]);
            ++evaluations;
        }

        for (int m = 0; m < dof; ++m) {
            double sq = 0.0;
            double sv = 0.0;
            double eq = 0.0;
            double ev = 0.0;
            for (int i = 0; i < stages; ++i) {
                double gi = g[i][m];
                sq += bbar[i] * gi;
                sv += b[i] * gi;
                eq += bbarStar[i] * gi;
                ev += bStar[i] * gi;
            }
            double qn = startQ[m] + h * startV[m] + h * h * sq;
            double vn = startV[m] + h * sv;
            yOut[m] = qn;
            yOut[dof + m] = vn;
            endQ[m] = qn;
            endV[m] = vn;
            if (errOut != null) {
                errOut[m] = h * h * (sq - eq);
                errOut[dof + m] = h * (sv - ev);
            }
        }

        if (fsal) {
            // the last stage sits at the end of the step and at the state the
            // step reached, so the next step can start from it
            carryT = t + h;
            System.arraycopy(endQ, 0, carryQ, 0, dof);
            System.arraycopy(g[stages - 1], 0, carryG, 0, dof);
            carryValid = true;
        }
        lastH = h;
        stepped = true;
    }

    @Override
    public void interpolate(double theta, double[] out) {
        if (!stepped) {
            throw new IllegalStateException("no step has been taken");
        }
        if (out == null || out.length != dimension) {
            throw new IllegalArgumentException("out must be of length " + dimension);
        }
        if (theta == 0.0) {
            System.arraycopy(startQ, 0, out, 0, dof);
            System.arraycopy(startV, 0, out, dof, dof);
            return;
        }
        if (theta == 1.0) {
            System.arraycopy(endQ, 0, out, 0, dof);
            System.arraycopy(endV, 0, out, dof, dof);
            return;
        }
        if (densePosition == null) {
            hermite(theta, out);
            return;
        }
        double theta2 = theta * theta;
        for (int m = 0; m < dof; ++m) {
            double sq = 0.0;
            double sv = 0.0;
            for (int i = 0; i < stages; ++i) {
                double gi = g[i][m];
                sq += horner(densePosition[i], theta) * gi;
                sv += horner(denseVelocity[i], theta) * gi;
            }
            out[m] = startQ[m] + theta * lastH * startV[m] + lastH * lastH * theta2 * sq;
            out[dof + m] = startV[m] + lastH * theta * sv;
        }
    }

    @Override
    public long evaluations() {
        return evaluations;
    }

    /**
     * Forgets the last step, the carried stage and the answer to the velocity
     * probe, so the next step asks again. The evaluation count is left alone, as
     * the interface requires.
     */
    @Override
    public void reset() {
        stepped = false;
        carryValid = false;
        velocityChecked = false;
        stageZeroReady = false;
    }

    /**
     * The method and the size of the system it advances.
     */
    @Override
    public String toString() {
        return "NystromRungeKutta[" + tableau + ", " + dof + " degrees of freedom]";
    }

    private void checkTheFieldIgnoresTheVelocity(double t) {
        // the first of the two evaluations sits at the first stage of the step,
        // so only the second one is spent on the question
        System.arraycopy(startV, 0, argV, 0, dof);
        field.valueAt(t, startQ, argV, probe);
        ++evaluations;
        for (int m = 0; m < dof; ++m) {
            // a shift large enough to matter and finite whatever the velocity is
            argV[m] = startV[m] + 1.0 + Math.abs(startV[m]);
        }
        field.valueAt(t, startQ, argV, argQ);
        ++evaluations;
        for (int m = 0; m < dof; ++m) {
            double scale = Math.max(Math.abs(probe[m]), 1.0);
            if (Math.abs(argQ[m] - probe[m]) > VELOCITY_PROBE_TOLERANCE * scale) {
                throw new IllegalArgumentException("the acceleration of component " + m
                        + " changed from " + probe[m] + " to " + g[0][m]
                        + " when only the velocity was varied, so this field is q'' = f(t, q, v)"
                        + " and not q'' = f(t, q). A Nystroem method has its order only on the"
                        + " second form, and its error estimate would understate the loss."
                        + " SymplecticNystrom accepts such a field, and so does Ode.solve on the"
                        + " flattened system");
            }
        }
        // the first of the two evaluations is the first stage of this step
        System.arraycopy(probe, 0, g[0], 0, dof);
        stageZeroReady = true;
    }

    private void hermite(double theta, double[] out) {
        double u2 = theta * theta;
        double u3 = u2 * theta;
        double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
        double h10 = u3 - 2.0 * u2 + theta;
        double h01 = -2.0 * u3 + 3.0 * u2;
        double h11 = u3 - u2;
        double d00 = 6.0 * u2 - 6.0 * theta;
        double d10 = 3.0 * u2 - 4.0 * theta + 1.0;
        double d01 = -6.0 * u2 + 6.0 * theta;
        double d11 = 3.0 * u2 - 2.0 * theta;
        for (int m = 0; m < dof; ++m) {
            out[m] = h00 * startQ[m] + h01 * endQ[m]
                    + lastH * (h10 * startV[m] + h11 * endV[m]);
            out[dof + m] = (d00 * startQ[m] + d01 * endQ[m]) / lastH + d10 * startV[m]
                    + d11 * endV[m];
        }
    }

    private static double horner(double[] coefficients, double theta) {
        double s = 0.0;
        for (int k = coefficients.length - 1; k >= 0; --k) {
            s = s * theta + coefficients[k];
        }
        return s;
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
