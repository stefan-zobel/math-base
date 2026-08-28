package math.ode;

import java.util.Locale;

import math.fun.DFunction;
import math.fun.DSecondOrderField;
import math.fun.DVectorField;
import math.fun.DiffDVectorField;
import math.solve.Quadrature;

/**
 * Solves an initial value problem, for a caller who does not want to choose a
 * method.
 * <p>
 * Most of what is here is {@link ExplicitRungeKutta} over
 * {@link ButcherTableau#DORMAND_PRINCE_45} driven by an {@link OdeIntegrator},
 * which is the right answer for a problem that is not stiff. A caller who wants
 * a different method, a step controller of their own, several events at once or
 * the classical fixed step method builds those three objects directly; this is
 * the shortcut and not the interface.
 * <p>
 * Four entries are different methods for different questions.
 * {@link #solveSymplectic(math.fun.DSecondOrderField, double, double[],
 * double[], double, int)} answers not "where is it now" but "where is it still,
 * after a very long time".
 * {@link #solveStiff(math.fun.DVectorField, double, double[], double, double)}
 * answers a question the explicit family cannot answer at all.
 * {@link #solveAccurate(math.fun.DVectorField, double, double[], double, double)}
 * answers the same question as {@code solve} but at a tolerance where a fifth
 * order method has become the expensive way to get there -- below about
 * {@code 1e-07}, measured. And
 * {@link #solveAuto(math.fun.DVectorField, double, double[], double, double)}
 * answers the question when the answer changes during the run.
 * <p>
 * <b>Read {@link OdeIntegrator.Result#seemsStiff}.</b> It is the one field here
 * that says the answer may be worthless: an explicit method whose step is held
 * down by stability rather than accuracy will either take a very long time or
 * run out of its budget, and no tolerance makes that better. When it comes back
 * {@code true}, hand the same problem to {@code solveStiff} -- or, if the
 * equation is only stiff in places, to {@code solveAuto}.
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
     * Integrates a mechanical system <code>q'' = f(t, q)</code> in equally long
     * steps with a symplectic method, for a question about a long time rather
     * than about accuracy now.
     * <p>
     * The method is {@link SplittingCoefficients#BLANES_MOAN_6}, order six at
     * eleven force evaluations a step. It is the cheapest of the five shipped
     * at any accuracy tighter than about {@code 1.6e-04}, which is where its
     * two extra orders stop paying for the extra evaluations, and no run worth
     * making symplectic is looser than that. The result holds the position and
     * the velocity stacked, position first.
     * <p>
     * <b>This is not the accurate answer, it is the stable one.</b> Over ten
     * orbits of a two body problem at matched cost, the adaptive fifth order
     * method above is nearly three orders of magnitude closer to the truth.
     * What this buys instead is that the energy error stays inside a band
     * instead of growing, and that the position error grows linearly with time
     * rather than as its square -- which is worth nothing over ten orbits and
     * decisive over ten thousand.
     * <p>
     * <b>The force must not read the velocity.</b> It is handed one because
     * {@link math.fun.DSecondOrderField} carries it; a field that uses it still
     * runs and is no longer symplectic, and nothing detects that.
     *
     * @param f
     *            the acceleration, which must depend on the time and the
     *            position only
     * @param t0
     *            the time the initial state belongs to
     * @param q0
     *            the initial position (not modified)
     * @param v0
     *            the initial velocity, of the same length (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param steps
     *            how many equal steps to take, at least one; the step size has
     *            to be constant, which is what the bounded behavior rests on
     * @return the position and the velocity at {@code steps + 1} equally spaced
     *         times, stacked into one state each
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the state stops being finite
     */
    public static OdeIntegrator.Result solveSymplectic(DSecondOrderField f, double t0, double[] q0,
            double[] v0, double t1, int steps) {
        if (q0 == null || v0 == null) {
            throw new IllegalArgumentException("q0 and v0 must not be null");
        }
        if (q0.length != v0.length) {
            throw new IllegalArgumentException("q0 is of length " + q0.length + " and v0 of length "
                    + v0.length + ", which must agree");
        }
        double[] y0 = new double[2 * q0.length];
        System.arraycopy(q0, 0, y0, 0, q0.length);
        System.arraycopy(v0, 0, y0, q0.length, v0.length);
        return new OdeIntegrator(new SymplecticNystrom(SplittingCoefficients.BLANES_MOAN_6, f, q0.length))
                .solve(t0, y0, t1, steps);
    }

    /**
     * Integrates {@code y' = f(t, y)} with {@link ButcherTableau#DOP853}, the
     * eighth order method, for an answer wanted to many digits.
     * <p>
     * <b>Use this when the tolerance is tight and not otherwise.</b> An eighth
     * order step costs twelve evaluations against Dormand-Prince's six, and it
     * pays for them by needing far fewer steps as the tolerance falls: the step
     * count grows like <code>rtol^(-1/9)</code> rather than
     * <code>rtol^(-1/5)</code>. Measured on the two body problem over ten
     * orbits, the crossing sits near {@code 1e-07} -- at {@code 1e-06}
     * {@link #solve(math.fun.DVectorField, double, double[], double, double)}
     * costs 831 evaluations against this one's 891, at {@code 1e-07} it is 1149
     * against 1095, at {@code 1e-10} it is 3951 against 2091, and at
     * {@code 1e-13} it is 15669 against 4035.
     * <p>
     * It also interpolates far more accurately, at {@code 2^8} per halving
     * against {@code 2^5}, and that is worth knowing for a run with an output
     * grid or an event: the three stages its continuous extension needs are
     * evaluated only inside a step somebody looked into.
     * <p>
     * <b>It is still an explicit method.</b> On a stiff equation it fails the
     * same way Dormand-Prince does, a little later; the answer to
     * {@link OdeIntegrator.Result#seemsStiff} is
     * {@link #solveStiff(math.fun.DVectorField, double, double[], double, double)}
     * and not this.
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
    public static OdeIntegrator.Result solveAccurate(DVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        return new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DOP853, f, y0.length),
                new StepController(tolerance, tolerance)).solve(t0, y0, t1);
    }

    /**
     * Integrates a stiff {@code y' = f(t, y)} with
     * {@link RosenbrockTableau#RODAS4}, differencing the Jacobian out of the
     * field.
     * <p>
     * <b>Use this when {@link OdeIntegrator.Result#seemsStiff} came back
     * {@code true}</b>, and not before. A stiff method is not a better method:
     * it costs a Jacobian, a matrix factorization and a back substitution per
     * stage, so on a problem that is not stiff it loses outright. Measured on
     * van der Pol's oscillator, where the stiffness is a dial: at
     * {@code mu = 10} the explicit method costs 789 evaluations against this
     * one's 1763, at {@code mu = 1000} they have crossed, at {@code mu = 1e4}
     * it is 74205 against 10034, and at {@code mu = 1e5} the explicit method
     * does not finish at all.
     * <p>
     * <b>What it costs to have no Jacobian.</b> Differencing one takes
     * {@code n + 1} further evaluations of the field per step, so this form
     * becomes expensive as the system grows. The overload taking a
     * {@link DiffDVectorField} takes none of them, and on a badly scaled
     * problem it also holds the step size up.
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
    public static OdeIntegrator.Result solveStiff(DVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        return new OdeIntegrator(new Rosenbrock(RosenbrockTableau.RODAS4, f, y0.length),
                new StepController(tolerance, tolerance)).solve(t0, y0, t1);
    }

    /**
     * Integrates a stiff {@code y' = f(t, y)} with
     * {@link RosenbrockTableau#RODAS4} and the Jacobian the field supplies,
     * which is the cheaper and the more accurate of the two forms.
     * <p>
     * Measured on Robertson's reaction, whose second concentration sits around
     * {@code 7e-08}: written out, the Jacobian takes the run to {@code 1e5} in
     * 305 steps and 1837 evaluations; differenced, in 314 steps and 1891. The
     * gap widens with the dimension, since a differenced Jacobian costs one
     * evaluation per component.
     *
     * @param f
     *            the right hand side together with its two first derivatives
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
    public static OdeIntegrator.Result solveStiff(DiffDVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        return new OdeIntegrator(new Rosenbrock(RosenbrockTableau.RODAS4, f, y0.length),
                new StepController(tolerance, tolerance)).solve(t0, y0, t1);
    }

    /**
     * Integrates {@code y' = f(t, y)} with both methods at once, letting each
     * step decide which of them the equation currently wants, for an equation
     * that is stiff over part of its run and not over the rest.
     * <p>
     * <b>Use this when the answer to "is it stiff" is "sometimes".</b> A
     * relaxation oscillator, a reaction with a fast transient, a system with a
     * dial in it: on those neither {@link #solve} nor {@link #solveStiff} is
     * right everywhere, and picking one of them is picking which half of the
     * run to be wrong about. Measured against the better of the two run alone,
     * and against a solver switching at exactly the right instants, which no
     * real one can beat:
     * <table border="1">
     * <caption>evaluations</caption>
     * <tr><th>problem</th><th>this</th><th>better pure</th><th>perfect</th></tr>
     * <tr><td>a dial from 1 to 1e5 and back, tolerance 1e-6</td>
     * <td>1728</td><td>2928</td><td>1468</td></tr>
     * <tr><td>van der Pol at mu = 1000, tolerance 1e-6</td>
     * <td>4698</td><td>8902</td><td>2747</td></tr>
     * <tr><td>Robertson to 1e5, tolerance 1e-8</td>
     * <td>3210</td><td>3147</td><td>3073</td></tr>
     * </table>
     * <p>
     * <b>It is not the answer to the other two cases, and it does not need to
     * be.</b> On an equation with no stiffness in it this takes no trial at all
     * and its run is {@link #solve}'s run <em>bit for bit</em>, so nothing is
     * lost by reaching for it; on one that is stiff from beginning to end,
     * Robertson's row above is the whole story -- a few percent for insurance
     * on a claim never made.
     * <p>
     * How it decides is {@link SwitchingStepper}, and a caller who wants to see
     * what it did builds that object instead and reads
     * {@link SwitchingStepper#switches()}.
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
     * @return the solution at the times the methods chose
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solveAuto(DVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        StepController controller = new StepController(tolerance, tolerance);
        return new OdeIntegrator(new SwitchingStepper(f, y0.length, controller), controller).solve(t0, y0,
                t1);
    }

    /**
     * The same, with the Jacobian the field supplies, which the implicit half
     * of the pair is the one to benefit from.
     *
     * @param f
     *            the right hand side together with its two first derivatives
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param tolerance
     *            the error one step may add, relative and absolute alike
     * @return the solution at the times the methods chose
     * @throws IllegalArgumentException
     *             if an argument is out of shape
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public static OdeIntegrator.Result solveAuto(DiffDVectorField f, double t0, double[] y0, double t1,
            double tolerance) {
        if (y0 == null) {
            throw new IllegalArgumentException("y0 must not be null");
        }
        StepController controller = new StepController(tolerance, tolerance);
        return new OdeIntegrator(new SwitchingStepper(f, y0.length, controller), controller).solve(t0, y0,
                t1);
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
     * A self check: nine claims about this package, each measured rather than
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

        // 5. and what a symplectic method does with the same orbit: the energy
        // error is no smaller, it simply stops growing
        DSecondOrderField gravity = new DSecondOrderField() {
            @Override
            public void valueAt(double t, double[] q, double[] v, double[] acceleration) {
                double r2 = q[0] * q[0] + q[1] * q[1];
                double r3 = r2 * Math.sqrt(r2);
                acceleration[0] = -q[0] / r3;
                acceleration[1] = -q[1] / r3;
            }
        };
        OdeIntegrator.Result orbit = solveSymplectic(gravity, 0.0, new double[] { start[0], start[1] },
                new double[] { start[2], start[3] }, 100.0 * 2.0 * Math.PI, 100 * 200);
        double firstBand = band(orbit, 0, 10 * 200);
        double lastBand = band(orbit, orbit.length - 10 * 200, orbit.length);
        System.out.printf(Locale.ROOT, "%-52s %12.3f   (band %.2e, first ten orbits against last ten)%n",
                "the same orbit, symplectically: the band grows by", Double.valueOf(lastBand / firstBand),
                Double.valueOf(firstBand));
        ok &= Math.abs(lastBand / firstBand - 1.0) < 0.05;

        // 6. and the one thing it says about its own limits
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

        // 7. and what to do about that: a problem the explicit family cannot
        // finish at any tolerance, with an invariant no method is told about
        DVectorField robertson = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
                dydt[2] = 3.0e7 * y[1] * y[1];
                dydt[1] = -dydt[0] - dydt[2];
            }
        };
        double[] concentrations = { 1.0, 0.0, 0.0 };
        OdeIntegrator.Result reaction = solveStiff(robertson, 0.0, concentrations, 1.0e11, 1.0e-8);
        double wandered = 0.0;
        for (int i = 0; i < reaction.length; ++i) {
            wandered = Math.max(wandered,
                    Math.abs(reaction.y[i][0] + reaction.y[i][1] + reaction.y[i][2] - 1.0));
        }
        boolean explicitGaveUp = false;
        try {
            new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, robertson, 3),
                    new StepController(1.0e-8, 1.0e-8, 20000)).solve(0.0, concentrations, 1.0e11);
        } catch (ArithmeticException expected) {
            explicitGaveUp = true;
        }
        System.out.printf(Locale.ROOT, "%-52s %12.3e   (%d evaluations, and the explicit method %s)%n",
                "Robertson to 1e11, its invariant conserved to", Double.valueOf(wandered),
                Long.valueOf(reaction.evaluations),
                explicitGaveUp ? "gives up" : "FINISHED, which it should not");
        ok &= wandered < 1.0e-14 && explicitGaveUp;

        // 8. and the third question, which is neither the horizon nor the
        // stiffness but the tolerance
        OdeIntegrator.Result five = solve(kepler, 0.0, start, 20.0, 1.0e-12);
        OdeIntegrator.Result eight = solveAccurate(kepler, 0.0, start, 20.0, 1.0e-12);
        System.out.printf(Locale.ROOT, "%-52s %12.2f   (%d evaluations against %d, at rtol 1e-12)%n",
                "ten orbits cost DOP853 this fraction of DP45", Double.valueOf(
                        eight.evaluations / (double) five.evaluations),
                Long.valueOf(eight.evaluations), Long.valueOf(five.evaluations));
        ok &= eight.evaluations * 2L < five.evaluations;

        // 9. and the fourth question, which is none of those three: an equation
        // whose own answer to "is it stiff" changes while it is being solved,
        // where picking either method is picking which half to be wrong about
        OdeIntegrator.Result switching = solveAuto(stiff, 0.0, new double[] { 2.0, 0.0 }, 300.0, 1.0e-6);
        OdeIntegrator.Result always = solveStiff(stiff, 0.0, new double[] { 2.0, 0.0 }, 300.0, 1.0e-6);
        // and where there is nothing to switch, it is the explicit run itself
        OdeIntegrator.Result smooth = solveAuto(oscillator, 0.0, new double[] { 1.0, 0.0 }, 20.0, 1.0e-10);
        boolean freeWhenIdle = smooth.evaluations == turning.evaluations;
        System.out.printf(Locale.ROOT, "%-52s %12.2f   (%d evaluations against %d, and %s)%n",
                "switching costs this fraction of staying implicit",
                Double.valueOf(switching.evaluations / (double) always.evaluations),
                Long.valueOf(switching.evaluations), Long.valueOf(always.evaluations),
                freeWhenIdle ? "free on a smooth problem" : "NOT free on a smooth problem");
        ok &= switching.evaluations * 3L < always.evaluations * 2L && freeWhenIdle;

        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }

    /**
     * The width of the energy error over a stretch of a two body run, whose
     * state is the position and then the velocity either way it was produced.
     */
    private static double band(OdeIntegrator.Result r, int from, int to) {
        double low = Double.MAX_VALUE;
        double high = -Double.MAX_VALUE;
        for (int i = from; i < to; ++i) {
            double[] y = r.y[i];
            double energy = 0.5 * (y[2] * y[2] + y[3] * y[3]) - 1.0 / Math.hypot(y[0], y[1]);
            low = Math.min(low, energy);
            high = Math.max(high, energy);
        }
        return high - low;
    }

    private static boolean report(String what, double measured, double allowed) {
        boolean ok = measured <= allowed;
        System.out.printf(Locale.ROOT, "%-52s %12.3e   %s%n", what, Double.valueOf(measured),
                ok ? "ok" : "FAILED, wanted at most " + allowed);
        return ok;
    }
}
