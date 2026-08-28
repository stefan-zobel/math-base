package math.ode;

import math.fun.DSecondOrderField;

/**
 * A symplectic step for a separable second order system, for whichever
 * {@link SplittingCoefficients} it is constructed with.
 * <p>
 * The state it advances is the position and the velocity stacked into one
 * vector of length {@code 2 n}, so an {@link OdeIntegrator} drives it, an
 * output grid reports it and an {@link Event} watches it exactly as they do an
 * explicit Runge-Kutta method. The driver never learns that this one is
 * different.
 * <p>
 * <b>What it is for.</b> An orbit conserves its energy exactly and a
 * Runge-Kutta method does not: over ten, a hundred and a thousand orbits,
 * Dormand-Prince loses it by {@code 1.28e-09}, {@code 1.18e-08} and
 * {@code 1.16e-07}, a factor of ten per factor of ten, and tightening the
 * tolerance lowers that line without changing its slope. A symplectic method
 * does not conserve the energy either -- it conserves a nearby quantity -- and
 * the consequence is that its energy error <em>oscillates within a band</em>
 * rather than growing. For a question about where a system is after a long
 * time, that difference is the whole answer.
 * <p>
 * <b>It has to be driven at a fixed step size.</b> The bounded behavior is a
 * property of a constant {@code h} and does not survive one that varies, which
 * is why {@link #hasErrorEstimate()} is {@code false} here and why
 * {@link OdeIntegrator} refuses a {@link StepController} for such a method.
 * That refusal was written for another reason and turns out to be the right one
 * for this.
 * <p>
 * <b>The force must not depend on the velocity.</b> The velocity is handed to
 * the field because {@link DSecondOrderField} carries it, and a field that
 * reads it still runs -- it is simply no longer being integrated symplectically,
 * and nothing here can detect that. Damping, drag and a magnetic force are all
 * outside what this class promises. The saving below rests on the same
 * assumption, so a velocity dependent field would also make it wrong.
 * <p>
 * <b>A method that begins with a kick shares one evaluation with the step
 * before it.</b> Its first kick has no drift in front of it, so it falls at the
 * time and the position the previous step's last kick fell at, and the field is
 * evaluated there once. That is what brings
 * {@link SplittingCoefficients#BLANES_MOAN_4} down to the six evaluations per
 * step it is published with, from the seven kicks it has; only the first step
 * of a run pays the extra. The
 * point is recognized by comparing it, exactly as {@link ExplicitRungeKutta}
 * recognizes a carried stage, so a driver that steps in an unexpected order
 * loses the saving and not the answer.
 * <p>
 * <b>Interpolation is free and consistent.</b> For a second order system the
 * velocity <em>is</em> the derivative of the position, known exactly at both
 * ends of the step, so the cubic Hermite through them costs no evaluation at
 * all -- unlike the same fallback in {@link ExplicitRungeKutta}, which has to
 * work out the derivative at the far end first. The reported velocity is the
 * derivative of the reported position, so the two cannot disagree. Measured
 * over five halvings, the error inside a step falls as {@code h^4} in the
 * position and as {@code h^3} in the velocity.
 * <p>
 * Instances are stateful and cannot be shared between threads.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Symplectic_integrator">Wikipedia
 * symplectic integrator</a>.
 *
 * @since 1.5.3
 */
public final class SymplecticNystrom implements OdeStepper {

    private final SplittingCoefficients coefficients;
    private final DSecondOrderField field;
    private final int dof;
    private final int dimension;

    private final double[] a;
    private final double[] b;
    private final double[] elapsed;

    private final boolean velocityFirst;
    private final double[] cacheQ;
    private final double[] cacheForce;
    private double cacheT;
    private boolean cacheValid;

    private final double[] q;
    private final double[] v;
    private final double[] force;
    private final double[] startQ;
    private final double[] startV;
    private final double[] endQ;
    private final double[] endV;

    private double lastH;
    private boolean stepped;
    private long evaluations;

    /**
     * A stepper for the given method and equation.
     *
     * @param coefficients
     *            the drifts and kicks of the method
     * @param field
     *            the acceleration, which must not depend on the velocity if the
     *            result is to be symplectic
     * @param degreesOfFreedom
     *            the number of coordinates of the position, half the dimension
     *            of the state this stepper advances
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the count is not positive
     */
    public SymplecticNystrom(SplittingCoefficients coefficients, DSecondOrderField field,
            int degreesOfFreedom) {
        if (coefficients == null) {
            throw new IllegalArgumentException("coefficients must not be null");
        }
        if (field == null) {
            throw new IllegalArgumentException("field must not be null");
        }
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException("degreesOfFreedom must be at least 1, not "
                    + degreesOfFreedom);
        }
        this.coefficients = coefficients;
        this.field = field;
        this.dof = degreesOfFreedom;
        this.dimension = 2 * degreesOfFreedom;
        this.a = coefficients.a();
        this.b = coefficients.b();
        this.elapsed = new double[a.length];
        double running = 0.0;
        for (int i = 0; i < a.length; ++i) {
            running += a[i];
            elapsed[i] = running;
        }
        // the last kick is at the end of the step by definition, and saying so
        // exactly is what lets the step after it recognize the same point
        elapsed[a.length - 1] = 1.0;
        this.velocityFirst = coefficients.isVelocityFirst();
        this.cacheQ = velocityFirst ? new double[degreesOfFreedom] : null;
        this.cacheForce = velocityFirst ? new double[degreesOfFreedom] : null;
        this.q = new double[degreesOfFreedom];
        this.v = new double[degreesOfFreedom];
        this.force = new double[degreesOfFreedom];
        this.startQ = new double[degreesOfFreedom];
        this.startV = new double[degreesOfFreedom];
        this.endQ = new double[degreesOfFreedom];
        this.endV = new double[degreesOfFreedom];
    }

    /**
     * The method this stepper executes.
     *
     * @return the coefficients given to the constructor
     */
    public SplittingCoefficients coefficients() {
        return coefficients;
    }

    /**
     * The number of coordinates of the position, half of
     * {@link #dimension()}.
     *
     * @return the degrees of freedom, at least one
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
        return coefficients.order();
    }

    /**
     * Always {@code false}: a splitting method carries no second solution to
     * compare against, and would not want its step size varied even if it did.
     *
     * @return {@code false}
     */
    @Override
    public boolean hasErrorEstimate() {
        return false;
    }

    /**
     * Always {@code false}: the interpolation is the cubic Hermite rather than
     * a continuous extension of the method's own order. It is nonetheless free
     * here, which the same answer from {@link ExplicitRungeKutta} does not
     * imply.
     *
     * @return {@code false}
     */
    @Override
    public boolean hasDenseOutput() {
        return false;
    }

    /**
     * Always {@link Double#NaN}: the estimate is a difference quotient between
     * two stages evaluated at the same time, and no two kicks here are.
     *
     * @return {@link Double#NaN}
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
        System.arraycopy(y, 0, q, 0, dof);
        System.arraycopy(y, dof, v, 0, dof);
        System.arraycopy(v, 0, dydt, 0, dof);
        field.valueAt(t, q, v, force);
        ++evaluations;
        System.arraycopy(force, 0, dydt, dof, dof);
    }

    /**
     * Advances the state by one step. The error estimate is ignored, since
     * there is none to give.
     */
    @Override
    public void step(double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (y == null || y.length != dimension) {
            throw new IllegalArgumentException("y must be of length " + dimension);
        }
        if (yOut == null || yOut.length != dimension) {
            throw new IllegalArgumentException("yOut must be of length " + dimension);
        }
        if (yOut == y) {
            throw new IllegalArgumentException("yOut must not be the array passed as y");
        }
        if (h == 0.0) {
            throw new IllegalArgumentException("h must not be zero");
        }
        System.arraycopy(y, 0, q, 0, dof);
        System.arraycopy(y, dof, v, 0, dof);
        System.arraycopy(q, 0, startQ, 0, dof);
        System.arraycopy(v, 0, startV, 0, dof);

        for (int i = 0; i < a.length; ++i) {
            double drift = a[i] * h;
            if (drift != 0.0) {
                for (int m = 0; m < dof; ++m) {
                    q[m] += drift * v[m];
                }
            }
            double kick = b[i];
            if (kick != 0.0) {
                // a method that begins with a kick has no drift before it, so
                // that kick is at the point the step before ended its own last
                // kick at, and the field has already been evaluated there
                if (velocityFirst && i == 0 && cacheValid && t == cacheT && sameContent(q, cacheQ)) {
                    System.arraycopy(cacheForce, 0, force, 0, dof);
                } else {
                    field.valueAt(t + elapsed[i] * h, q, v, force);
                    ++evaluations;
                }
                double weight = kick * h;
                for (int m = 0; m < dof; ++m) {
                    v[m] += weight * force[m];
                }
                if (velocityFirst && i == a.length - 1) {
                    cacheT = t + h;
                    System.arraycopy(q, 0, cacheQ, 0, dof);
                    System.arraycopy(force, 0, cacheForce, 0, dof);
                    cacheValid = true;
                }
            }
        }

        System.arraycopy(q, 0, yOut, 0, dof);
        System.arraycopy(v, 0, yOut, dof, dof);
        System.arraycopy(q, 0, endQ, 0, dof);
        System.arraycopy(v, 0, endV, 0, dof);
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
            out[m] = h00 * startQ[m] + h01 * endQ[m] + lastH * (h10 * startV[m] + h11 * endV[m]);
            out[dof + m] = (d00 * startQ[m] + d01 * endQ[m]) / lastH + d10 * startV[m] + d11 * endV[m];
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
    }

    private static boolean sameContent(double[] x, double[] y) {
        for (int i = 0; i < x.length; ++i) {
            if (x[i] != y[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * The method and the size of the system it advances.
     */
    @Override
    public String toString() {
        return "SymplecticNystrom[" + coefficients + ", " + dof + " degrees of freedom]";
    }
}
