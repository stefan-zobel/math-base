package math.ode;

import java.util.Locale;

import math.fun.DFunction;
import math.fun.DVectorField;
import math.solve.Quadrature;

/**
 * Solves an initial value problem, for a caller who does not want to choose a
 * method.
 * <p>
 * Everything here is {@link ExplicitRungeKutta} over
 * {@link ButcherTableau#DORMAND_PRINCE_45} driven by an
 * {@link OdeIntegrator}, which is the right answer for a problem that is not
 * stiff and is the only answer this package has for one that is. A caller who
 * wants a different method, a step controller of their own, several events at
 * once or the classical fixed step method builds those three objects directly;
 * this is the shortcut and not the interface.
 * <p>
 * <b>Read {@link OdeIntegrator.Result#seemsStiff}.</b> It is the one field here
 * that says the answer may be worthless: an explicit method whose step is held
 * down by stability rather than accuracy will either take a very long time or
 * run out of its budget, and no tolerance makes that better.
 * <p>
 * The fixed step form uses {@link ButcherTableau#CLASSIC_RK4} instead, because
 * that is what a fixed step size usually means and because its four stages are
 * the cheapest fourth order step there is. Dormand-Prince at the same step
 * count costs half again as much and buys an order, which is usually the better
 * trade -- {@code new OdeIntegrator(new ExplicitRungeKutta(
 * ButcherTableau.DORMAND_PRINCE_45, f, n))} is that trade.
 *
 * @since 1.5.3
 */
public final class Ode {

    /** The tolerance used where none is given, relative and absolute alike. */
    public static final double DEFAULT_TOLERANCE = 1.0e-8;

    /**
     * Integrates {@code y' = f(t, y)} from {@code t0} to {@code t1} at the
     * default tolerance, keeping the state at the end of every step taken.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @return the solution at the times the method chose
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solve(DVectorField f, double t0, double[] y0, double t1) {
        return solve(f, t0, y0, t1, DEFAULT_TOLERANCE);
    }

    /**
     * Integrates {@code y' = f(t, y)} from {@code t0} to {@code t1}, keeping
     * the state at the end of every step taken.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param tolerance
     *            the error one step may add, relative and absolute alike
     * @return the solution at the times the method chose
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solve(DVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        return integrator(f, y0, tolerance, null).solve(t0, y0, t1);
    }

    /**
     * Integrates {@code y' = f(t, y)} and reports the state at the given times
     * and at no others, which costs nothing beyond the run itself.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param times
     *            the times to report at, running from {@code t0} towards
     *            {@code t1} without turning back
     * @param tolerance
     *            the error one step may add, relative and absolute alike
     * @return the solution at the requested times
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solve(DVectorField f, double t0, double[] y0, double t1,
            double[] times, double tolerance) {
        return integrator(f, y0, tolerance, null).solve(t0, y0, t1, times);
    }

    /**
     * Integrates until the given quantity crosses zero, or until {@code t1} if
     * it never does.
     * <p>
     * The crossing is found on the interpolant of the step that brackets it, so
     * its time is as precise as the trajectory it sits on and has nothing to do
     * with the step size. Whether it happened at all is
     * {@link OdeIntegrator.Result#stoppedByEvent}, and when is
     * {@code eventTimes[0]}.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to give up at
     * @param stop
     *            the quantity whose first zero ends the run, crossed either way
     * @param tolerance
     *            the error one step may add, relative and absolute alike
     * @return the solution up to the crossing, or up to {@code t1}
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solveUntil(DVectorField f, double t0, double[] y0, double t1,
            OdeEvent stop, double tolerance) {
        Event event = Event.terminal(stop, Event.Direction.EITHER);
        return integrator(f, y0, tolerance, new Event[] { event }).solve(t0, y0, t1);
    }

    /**
     * Integrates in equally long steps with the classical fourth order method,
     * which estimates no error and so asks for no tolerance.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param steps
     *            how many equal steps to take, at least one
     * @return the solution at {@code steps + 1} equally spaced times
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the state stops being finite
     */
    public static OdeIntegrator.Result solveFixed(DVectorField f, double t0, double[] y0, double t1,
            int steps) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        return new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, f, y0.length))
                .solve(t0, y0, t1, steps);
    }

    /**
     * The state at {@code t1} and nothing else, for a caller who wants the
     * answer rather than the path to it.
     *
     * @param f
     *            the right hand side
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param tolerance
     *            the error one step may add, relative and absolute alike
     * @return a fresh <code>double[]</code> holding the state at {@code t1}
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static double[] endpoint(DVectorField f, double t0, double[] y0, double t1, double tolerance) {
        return integrator(f, y0, tolerance, null).solve(t0, y0, t1, new double[] { t1 }).finalState();
    }

    private static OdeIntegrator integrator(DVectorField f, double[] y0, double tolerance, Event[] events) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        ExplicitRungeKutta stepper = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, f, y0.length);
        return new OdeIntegrator(stepper, new StepController(tolerance, tolerance), events);
    }

    private Ode() {
        throw new AssertionError();
    }

    /**
     * A self check: five claims about this package, each measured rather than
     * asserted, and a verdict at the end.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        boolean ok = true;

        // 1. a solution that is known exactly
        DVectorField oscillator = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = y[1];
                dydt[1] = -y[0];
            }
        };
        OdeIntegrator.Result turning = solve(oscillator, 0.0, new double[] { 1.0, 0.0 }, 20.0, 1.0e-10);
        double worst = 0.0;
        for (int i = 0; i < turning.length; ++i) {
            worst = Math.max(worst, Math.hypot(turning.y[i][0] - Math.cos(turning.t[i]),
                    turning.y[i][1] + Math.sin(turning.t[i])));
        }
        ok &= report("y'' = -y over [0, 20] against the cosine", worst, 1.0e-8);

        // 2. an equation with no y in it is a quadrature, which math.solve does
        // by a route with nothing in common
        final DFunction integrand = new DFunction() {
            @Override
            public double apply(double t) {
                return Math.exp(-t * t) * Math.cos(3.0 * t) + 0.5;
            }
        };
        DVectorField asAnEquation = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = integrand.apply(t);
            }
        };
        double byOde = endpoint(asAnEquation, 0.0, new double[] { 0.0 }, 3.0, 1.0e-12)[0];
        double byQuadrature = Quadrature.integrate(integrand, 0.0, 3.0, 1.0e-13);
        ok &= report("against Gauss-Kronrod on the same integrand", Math.abs(byOde - byQuadrature), 1.0e-11);

        // 3. an event, whose time owes nothing to the step size
        OdeIntegrator.Result landing = solveUntil(oscillator, 0.0, new double[] { 1.0, 0.0 }, 20.0,
                new OdeEvent() {
                    @Override
                    public double valueAt(double t, double[] y) {
                        return y[0];
                    }
                }, 1.0e-12);
        ok &= report("the first zero of the cosine, against pi / 2",
                Math.abs(landing.eventTimes[0] - 0.5 * Math.PI), 1.0e-11);

        // 4. what this kind of method cannot do: keep an invariant
        DVectorField kepler = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double r2 = y[0] * y[0] + y[1] * y[1];
                double r3 = r2 * Math.sqrt(r2);
                dydt[0] = y[2];
                dydt[1] = y[3];
                dydt[2] = -y[0] / r3;
                dydt[3] = -y[1] / r3;
            }
        };
        double e = 0.6;
        double[] start = { 1.0 - e, 0.0, 0.0, Math.sqrt((1.0 + e) / (1.0 - e)) };
        double energy0 = 0.5 * (start[2] * start[2] + start[3] * start[3]) - 1.0 / Math.hypot(start[0],
                start[1]);
        double[] after = endpoint(kepler, 0.0, start, 100.0 * 2.0 * Math.PI, 1.0e-10);
        double energy1 = 0.5 * (after[2] * after[2] + after[3] * after[3])
                - 1.0 / Math.hypot(after[0], after[1]);
        System.out.printf(Locale.ROOT, "%-52s %12.3e   (drifts, and is meant to)%n",
                "Kepler energy over a hundred orbits", Double.valueOf(Math.abs(energy1 - energy0)));

        // 5. and the one thing it says about its own limits
        DVectorField stiff = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = y[1];
                dydt[1] = 100.0 * (1.0 - y[0] * y[0]) * y[1] - y[0];
            }
        };
        OdeIntegrator.Result vanDerPol = solve(stiff, 0.0, new double[] { 2.0, 0.0 }, 300.0, 1.0e-6);
        System.out.printf(Locale.ROOT, "%-52s %12s   (%d steps, measure %.2f)%n",
                "van der Pol at mu = 100 reports itself stiff",
                Boolean.toString(vanDerPol.seemsStiff), Long.valueOf(vanDerPol.steps),
                Double.valueOf(vanDerPol.stiffness));
        ok &= vanDerPol.seemsStiff;

        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }

    private static boolean report(String what, double measured, double allowed) {
        boolean ok = measured <= allowed;
        System.out.printf(Locale.ROOT, "%-52s %12.3e   %s%n", what, Double.valueOf(measured),
                ok ? "ok" : "FAILED, wanted at most " + allowed);
        return ok;
    }
}
