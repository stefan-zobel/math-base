package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorField;

/**
 * What the switching stepper has to be true of, measured rather than assumed.
 */
public class SwitchingStepperTest {

    /** How far a run on a grid may sit from an independent implicit run. */
    private static final double ACROSS_METHODS = 1.0e-4;

    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    private static final DVectorField ROBERTSON = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            dydt[2] = 3.0e7 * y[1] * y[1];
            dydt[1] = -dydt[0] - dydt[2];
        }
    };

    /** A dial that turns the stiffness from one up to 1e5 and back down. */
    private static final DVectorField MOVING_DIAL = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            double width = 0.25;
            double lambda = 1.0 + (1.0e5 - 1.0) * 0.5 * (1.0 + Math.tanh((t - 10.0) / width)) * 0.5
                    * (1.0 + Math.tanh((20.0 - t) / width));
            dydt[0] = -lambda * (y[0] - Math.cos(t)) - Math.sin(t);
        }
    };

    /**
     * The dial again, but on a chain of {@code n} cells diffusing into each
     * other, so that the dimension can be turned up while the spectrum stays
     * where it is -- the diffusion is scaled with the mesh on purpose, since
     * what is being varied here is {@code n} and nothing else.
     */
    private static DVectorField chain(final int n) {
        final double h = 1.0 / (n + 1);
        final double[] profile = new double[n];
        for (int i = 0; i < n; ++i) {
            profile[i] = Math.sin(Math.PI * (i + 1) * h);
        }
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double width = 0.25;
                double lambda = 1.0 + (1.0e5 - 1.0) * 0.5 * (1.0 + Math.tanh((t - 10.0) / width))
                        * 0.5 * (1.0 + Math.tanh((20.0 - t) / width));
                double c = Math.cos(t);
                for (int i = 0; i < n; ++i) {
                    double left = (i == 0) ? 0.0 : y[i - 1];
                    double right = (i == n - 1) ? 0.0 : y[i + 1];
                    dydt[i] = 2.5 * (left - 2.0 * y[i] + right) - lambda * (y[i] - c * profile[i]);
                }
            }
        };
    }

    private static double[] chainStart(int n) {
        double[] y = new double[n];
        double h = 1.0 / (n + 1);
        for (int i = 0; i < n; ++i) {
            y[i] = Math.sin(Math.PI * (i + 1) * h);
        }
        return y;
    }

    /**
     * Stiff modes that <em>rotate</em> as well as decay: {@code n / 2} pairs
     * whose decay rates are log-spaced from 1 to 1e5 under a dial, each turning
     * between 30 and 3 times as fast as it decays. An implicit method has to
     * resolve a rotation and cannot ramp away from it in a step or two, so the
     * step size it reaches climbs slowly and over many steps -- which is what a
     * trial of a fixed length cannot see. The exact solution is
     * {@code y_i = cos(t)}, so cost and accuracy can both be asserted.
     */
    private static DVectorField rotating(final int n) {
        final int pairs = n / 2;
        final double[] base = new double[n];
        final double[] turn = new double[pairs];
        for (int i = 0; i < n; ++i) {
            base[i] = Math.pow(1.0e5, i / (double) (n - 1));
        }
        for (int p = 0; p < pairs; ++p) {
            turn[p] = 30.0 * Math.pow(0.1, p / (double) (pairs - 1));
        }
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double w = 0.5 * (1.0 + Math.tanh((t - 8.0) / 0.25)) * 0.5
                        * (1.0 + Math.tanh((24.0 - t) / 0.25));
                double c = Math.cos(t);
                double s = Math.sin(t);
                for (int i = 0; i < n; i += 2) {
                    double m = 1.0 + (base[i] - 1.0) * w;
                    double r = m * turn[i / 2];
                    double d0 = y[i] - c;
                    double d1 = y[i + 1] - c;
                    dydt[i] = -m * d0 + r * d1 - s;
                    dydt[i + 1] = -r * d0 - m * d1 - s;
                }
            }
        };
    }

    /**
     * Decays whose rates are log-spaced from 1 to 100 under the same dial: stiff
     * enough that the explicit method's step is held down and trials are taken,
     * never stiff enough that handing over would pay. The exact solution is
     * {@code y_i = cos(t)}.
     */
    private static DVectorField spread(final int n) {
        final double[] base = new double[n];
        for (int i = 0; i < n; ++i) {
            base[i] = Math.pow(100.0, i / (double) (n - 1));
        }
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double w = 0.5 * (1.0 + Math.tanh((t - 10.0) / 0.25)) * 0.5
                        * (1.0 + Math.tanh((20.0 - t) / 0.25));
                double c = Math.cos(t);
                double s = Math.sin(t);
                for (int i = 0; i < n; ++i) {
                    dydt[i] = -(1.0 + (base[i] - 1.0) * w) * (y[i] - c) - s;
                }
            }
        };
    }

    private static double[] ones(int n) {
        double[] y = new double[n];
        for (int i = 0; i < n; ++i) {
            y[i] = 1.0;
        }
        return y;
    }

    /** A linear equation whose one eigenvalue is a number this test knows. */
    private static final DVectorField DECAY = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = -1.0e4 * y[0];
        }
    };

    private static DVectorField vanDerPol(final double mu) {
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = y[1];
                dydt[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
            }
        };
    }

    /**
     * On an equation that is never stiff the switch never happens, and a run
     * that never switches is not merely close to a pure explicit one -- it does
     * the same arithmetic in the same order, so it agrees bit for bit.
     */
    @Test
    public void aSmoothRunIsThePureExplicitRun() {
        for (int e = 6; e <= 12; e += 2) {
            double tol = Math.pow(10.0, -e);
            StepController controller = new StepController(tol, tol);
            SwitchingStepper switcher = new SwitchingStepper(OSCILLATOR, 2, controller);
            OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                    .solve(0.0, new double[] { 1.0, 0.0 }, 50.0);
            OdeIntegrator.Result pure = new OdeIntegrator(
                    new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2),
                    new StepController(tol, tol)).solve(0.0, new double[] { 1.0, 0.0 }, 50.0);

            assertEquals("switches at 1e-" + e, 0L, switcher.switches());
            assertEquals("implicit steps at 1e-" + e, 0L, switcher.stiffSteps());
            assertEquals("evaluations at 1e-" + e, pure.evaluations, mixed.evaluations);
            assertEquals("recorded points at 1e-" + e, pure.length, mixed.length);
            for (int i = 0; i < pure.length; ++i) {
                assertEquals("t[" + i + "] at 1e-" + e, pure.t[i], mixed.t[i], 0.0);
                assertEquals("y[" + i + "][0] at 1e-" + e, pure.y[i][0], mixed.y[i][0], 0.0);
                assertEquals("y[" + i + "][1] at 1e-" + e, pure.y[i][1], mixed.y[i][1], 0.0);
            }
        }
    }

    /**
     * Robertson's reaction is stiff from its first transient onwards, so the
     * switch happens once and nothing sends it back. The invariant nothing is
     * told about is what says the answer survived the handover.
     */
    @Test
    public void robertsonHandsOverOnceAndKeepsItsInvariant() {
        StepController controller = new StepController(1.0e-8, 1.0e-8);
        SwitchingStepper switcher = new SwitchingStepper(ROBERTSON, 3, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                .solve(0.0, new double[] { 1.0, 0.0, 0.0 }, 1.0e5);

        assertEquals("hands over once and is never tempted back", 1L, switcher.switches());
        assertTrue("ends on the implicit side", switcher.isStiffActive());
        assertTrue("having asked more than once", switcher.trials() > 1L);

        double wandered = 0.0;
        for (int i = 0; i < mixed.length; ++i) {
            wandered = Math.max(wandered,
                    Math.abs(mixed.y[i][0] + mixed.y[i][1] + mixed.y[i][2] - 1.0));
        }
        assertTrue("the concentrations sum to one, off by " + wandered, wandered < 1.0e-13);

        // and the answer agrees with the method that carried most of the run
        OdeIntegrator.Result pure = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, ROBERTSON, 3),
                new StepController(1.0e-8, 1.0e-8)).solve(0.0, new double[] { 1.0, 0.0, 0.0 }, 1.0e5);
        double gap = Math.abs(mixed.finalState()[0] / pure.finalState()[0] - 1.0);
        assertTrue("relative gap on the first concentration is " + gap, gap < 1.0e-6);
        assertTrue("costs no more than a fifth again as much as the pure implicit run, "
                + mixed.evaluations + " against " + pure.evaluations,
                mixed.evaluations < 1.2 * pure.evaluations);
    }

    /**
     * Van der Pol at {@code mu = 1000} is the case switching exists for: the
     * explicit method alone needs over ten million evaluations, and the
     * implicit one is beaten outright because it is made to hand the brief
     * violent transitions back.
     */
    @Test
    public void vanDerPolBeatsBothMethodsRunAlone() {
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(vanDerPol(1000.0), 2, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                .solve(0.0, new double[] { 2.0, 0.0 }, 3000.0);
        OdeIntegrator.Result pure = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, vanDerPol(1000.0), 2),
                new StepController(1.0e-6, 1.0e-6)).solve(0.0, new double[] { 2.0, 0.0 }, 3000.0);

        assertEquals("reaches the end", 3000.0, mixed.finalTime(), 0.0);
        assertTrue("goes both ways rather than once, not " + switcher.switches(),
                switcher.switches() > 1L);
        assertTrue("and steps on both sides",
                switcher.nonStiffSteps() > 0L && switcher.stiffSteps() > 0L);
        assertTrue("costs " + mixed.evaluations + " against the pure implicit run's "
                + pure.evaluations, mixed.evaluations < 0.75 * pure.evaluations);
        assertEquals("and lands where the implicit method lands", pure.finalState()[0],
                mixed.finalState()[0], 1.0e-4);
    }

    /**
     * A dial that turns the stiffness up and back down again is the controlled
     * case: the explicit method is right on the outer thirds and hopeless in
     * the middle, so a solver that never changes its mind cannot be right
     * everywhere.
     */
    @Test
    public void theDialIsCheaperThanEitherMethodAlone() {
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(MOVING_DIAL, 1, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                .solve(0.0, new double[] { 1.0 }, 30.0);
        OdeIntegrator.Result pure = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, MOVING_DIAL, 1),
                new StepController(1.0e-6, 1.0e-6)).solve(0.0, new double[] { 1.0 }, 30.0);

        assertTrue("goes both ways, not " + switcher.switches(), switcher.switches() > 1L);
        assertTrue("costs " + mixed.evaluations + " against the pure implicit run's "
                + pure.evaluations, mixed.evaluations < 0.75 * pure.evaluations);
        // the solution the dial follows is the cosine, which it is dragged onto
        // however hard the dial is turned
        assertEquals("and follows the cosine it is pulled towards", Math.cos(30.0),
                mixed.finalState()[0], 1.0e-5);
    }

    /**
     * <b>The dimension must not decide whether the two methods change places.</b>
     * <p>
     * A trial judged on the step size the other method would propose
     * <em>next</em> can only ever hand over while one implicit step costs less
     * than a couple of explicit ones, since that one-step proposal is all the
     * ramp it can see. With a differenced Jacobian that is
     * <code>n + 7 &lt; 6 * 2.7</code>, and above {@code n = 9} such a stepper
     * never hands over at all: on this chain it took 2072668 evaluations at
     * {@code n = 10} where the implicit method alone takes 6142. Letting the
     * trial settle is what removes the bound, and this is the test that says so.
     */
    @Test
    public void itHandsOverAtEveryDimension() {
        int[] sizes = { 10, 50 };
        for (int q = 0; q < sizes.length; ++q) {
            int n = sizes[q];
            DVectorField f = chain(n);
            StepController controller = new StepController(1.0e-6, 1.0e-6);
            SwitchingStepper switcher = new SwitchingStepper(f, n, controller);
            OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller).solve(0.0,
                    chainStart(n), 30.0);
            OdeIntegrator.Result pure = new OdeIntegrator(
                    new Rosenbrock(RosenbrockTableau.RODAS4, f, n),
                    new StepController(1.0e-6, 1.0e-6)).solve(0.0, chainStart(n), 30.0);

            assertTrue("it hands over at n = " + n + ", not " + switcher.switches() + " switches",
                    switcher.switches() >= 1L);
            assertTrue("and at n = " + n + " it beats the pure implicit run, " + mixed.evaluations
                    + " against " + pure.evaluations, mixed.evaluations < 0.7 * pure.evaluations);
            assertTrue("having explored to find out", switcher.explorationEvaluations() > 0L);
            for (int i = 0; i < n; ++i) {
                assertEquals("component " + i + " at n = " + n, pure.finalState()[i],
                        mixed.finalState()[i], ACROSS_METHODS);
            }
        }
    }

    /**
     * The dimension bound came back once already, and this is the shape that
     * brought it. A trial of a <em>fixed</em> length reveals a fixed multiple of
     * the step it starts from, while the ratio a switch has to clear,
     * <code>(n + 7) / 6</code>, grows with the dimension; so any fixed length
     * has an {@code n} above which the stepper stops handing over. Four steps
     * put that at about {@code n = 21}: on these rotating pairs at {@code n = 50}
     * a four-step trial reveals {@code 4.6 h} where {@code 9.5 h} is wanted, the
     * stepper never switches, and the run dies at the step budget having spent
     * 2612749 evaluations where the implicit method alone needs 71993. Settling
     * until the answer is in has no such bound, and this is the test that says
     * so.
     */
    @Test
    public void itHandsOverWhenTheStiffModesRotate() {
        int n = 50;
        DVectorField f = rotating(n);
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(f, n, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller).solve(0.0, ones(n),
                30.0);
        OdeIntegrator.Result pure = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, f, n),
                new StepController(1.0e-6, 1.0e-6)).solve(0.0, ones(n), 30.0);

        assertTrue("it hands over, not " + switcher.switches() + " switches",
                switcher.switches() >= 1L);
        assertTrue("and beats the pure implicit run, " + mixed.evaluations + " against "
                + pure.evaluations, mixed.evaluations < pure.evaluations);
        // y_i = cos(t) solves this exactly, so being cheap is not enough
        for (int i = 0; i < n; ++i) {
            assertEquals("component " + i, Math.cos(30.0), mixed.finalState()[i], 1.0e-3);
        }
    }

    /**
     * What a trial costs when it learns nothing. A trial that ends in a switch
     * is free, because the switch is the step; only one answered "stay" is paid
     * for, and it is paid for once per {@code probeEvery} steps. On a problem
     * whose answer is always "stay" the climb dies away at once, so the trial
     * should learn that in about <em>one</em> implicit step and stop. A trial of
     * four fixed steps spent 108 evaluations here where this spends 27, and at
     * {@code n = 50} it was 228 against 95.
     */
    @Test
    public void aTrialThatLearnsNothingStopsEarly() {
        int n = 20;
        DVectorField f = spread(n);
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(f, n, controller);
        new OdeIntegrator(switcher, controller).solve(0.0, ones(n), 30.0);

        assertEquals("this problem never wants the implicit method", 0L, switcher.switches());
        assertTrue("but it does ask", switcher.trials() > 0L);
        // an implicit step here costs the six stages plus the n + 1 the
        // difference quotient takes
        double oneStep = 6 + n + 1;
        double perTrial = switcher.explorationEvaluations() / (double) switcher.trials();
        assertTrue("a trial that learns nothing costs " + perTrial + " evaluations, more than the "
                + (2.0 * oneStep) + " that two implicit steps would",
                perTrial < 2.0 * oneStep);
    }

    /**
     * An interior value is asked for on whichever side took the step, so the
     * whole grid is answered even though the method changed halfway through.
     */
    @Test
    public void aGridIsAnsweredAcrossTheSwitch() {
        int points = 200;
        double[] grid = new double[points];
        for (int i = 0; i < points; ++i) {
            grid[i] = 3000.0 * (i + 1) / points;
        }
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(vanDerPol(1000.0), 2, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                .solve(0.0, new double[] { 2.0, 0.0 }, 3000.0, grid);
        OdeIntegrator.Result pure = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, vanDerPol(1000.0), 2),
                new StepController(1.0e-6, 1.0e-6)).solve(0.0, new double[] { 2.0, 0.0 }, 3000.0, grid);

        assertEquals("points", points, mixed.length);
        double worst = 0.0;
        for (int i = 0; i < points; ++i) {
            assertEquals("t[" + i + "]", grid[i], mixed.t[i], 0.0);
            worst = Math.max(worst, Math.abs(mixed.y[i][0] - pure.y[i][0]));
        }
        assertTrue("the two runs differ by " + worst + " on the grid", worst < ACROSS_METHODS);
    }

    /**
     * The endpoints of the interpolant are exact whichever side produced them,
     * which is what keeps an output time that lands on a step boundary from
     * disagreeing with the step.
     */
    @Test
    public void theEndsOfEveryStepAreExactOnBothSides() {
        SwitchingStepper switcher = new SwitchingStepper(ButcherTableau.DORMAND_PRINCE_45,
                RosenbrockTableau.RODAS4, DECAY, 1, new StepController(1.0e-6, 1.0e-6), 3);
        double[] y = { 1.0 };
        double[] next = new double[1];
        double[] error = new double[1];
        double[] at0 = new double[1];
        double[] at1 = new double[1];
        double[] middle = new double[1];
        double t = 0.0;
        double h = 1.0e-3;
        for (int i = 0; i < 20; ++i) {
            double before = y[0];
            switcher.step(t, y, h, next, error);
            switcher.interpolate(0.0, at0);
            switcher.interpolate(1.0, at1);
            switcher.interpolate(0.5, middle);
            assertEquals("start of step " + i, before, at0[0], 0.0);
            assertEquals("end of step " + i, next[0], at1[0], 0.0);
            assertTrue("the middle of step " + i + " is a number",
                    !Double.isNaN(middle[0]) && !Double.isInfinite(middle[0]));
            double[] swap = y;
            y = next;
            next = swap;
            t += h;
        }
        assertTrue("the explicit side did the stepping here", switcher.nonStiffSteps() > 0L);

        // and the same on the implicit side, reached by driving a run there
        // first, since the steps taken to find out where the implicit method
        // settles are the scout's and never the answer
        StepController controller = new StepController(1.0e-8, 1.0e-8);
        SwitchingStepper stiffly = new SwitchingStepper(ROBERTSON, 3, controller);
        OdeIntegrator.Result driven = new OdeIntegrator(stiffly, controller).solve(0.0,
                new double[] { 1.0, 0.0, 0.0 }, 1.0e3);
        assertTrue("which ends on the implicit side", stiffly.isStiffActive());

        double[] state = driven.finalState().clone();
        double[] reached = new double[3];
        double[] estimate = new double[3];
        double[] start = new double[3];
        double[] end = new double[3];
        double[] inside = new double[3];
        double time = 1.0e3;
        for (int i = 0; i < 6; ++i) {
            double[] before = state.clone();
            stiffly.step(time, state, 10.0, reached, estimate);
            stiffly.interpolate(0.0, start);
            stiffly.interpolate(1.0, end);
            stiffly.interpolate(0.5, inside);
            for (int m = 0; m < 3; ++m) {
                assertEquals("start of implicit step " + i + ", component " + m, before[m], start[m],
                        0.0);
                assertEquals("end of implicit step " + i + ", component " + m, reached[m], end[m], 0.0);
                assertTrue("the middle of implicit step " + i + " is a number",
                        !Double.isNaN(inside[m]) && !Double.isInfinite(inside[m]));
            }
            double[] swap = state;
            state = reached;
            reached = swap;
            time += 10.0;
        }
    }

    /**
     * The order the controller scales with has to be the order of the method
     * that took the step, so it changes when the method does -- which is why
     * {@link OdeIntegrator} asks for it after every step and not once.
     */
    @Test
    public void theOrderAndTheMeasureFollowTheActiveMethod() {
        SwitchingStepper byHand = new SwitchingStepper(ButcherTableau.DORMAND_PRINCE_45,
                RosenbrockTableau.RODAS4, DECAY, 1, new StepController(1.0e-6, 1.0e-6), 1000);
        assertEquals("before any step it is the explicit order", 5, byHand.order());
        assertTrue("and there is no measure yet", Double.isNaN(byHand.stiffnessMeasure()));

        double[] y = { 1.0 };
        double[] next = new double[1];
        double[] error = new double[1];
        double t = 0.0;
        double h = 1.0e-3;
        for (int i = 0; i < 3; ++i) {
            byHand.step(t, y, h, next, error);
            assertEquals("still explicit at step " + i, 5, byHand.order());
            // h |lambda| is 1e-3 times 1e4, and the difference quotient finds it
            assertEquals("the measure at step " + i, 10.0, byHand.stiffnessMeasure(), 0.05);
            double[] swap = y;
            y = next;
            next = swap;
            t += h;
        }

        // and driven properly, on a problem that does change its mind: whatever
        // side ended up stepping is the side the two answers belong to
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper driven = new SwitchingStepper(vanDerPol(1000.0), 2, controller);
        new OdeIntegrator(driven, controller).solve(0.0, new double[] { 2.0, 0.0 }, 3000.0);
        int expected = driven.isStiffActive() ? RosenbrockTableau.RODAS4.order()
                : ButcherTableau.DORMAND_PRINCE_45.order();
        assertEquals("the order is the last stepping method's", expected, driven.order());
        assertEquals("and only the explicit side has a measure to give",
                Boolean.valueOf(driven.isStiffActive()),
                Boolean.valueOf(Double.isNaN(driven.stiffnessMeasure())));
    }

    /**
     * Steps are counted as they are attempted and not as they are kept, on
     * purpose: a run whose steps are all being rejected is a run that may want
     * the other method, and one that waited for an accepted step before asking
     * could shrink towards zero without ever having asked.
     */
    @Test
    public void aRunThatOnlyEverRetriesStillAsksTheOtherMethod() {
        SwitchingStepper switcher = new SwitchingStepper(ButcherTableau.DORMAND_PRINCE_45,
                RosenbrockTableau.RODAS4, DECAY, 1, new StepController(1.0e-6, 1.0e-6), 10);
        double[] y = { 1.0 };
        double[] next = new double[1];
        double[] error = new double[1];
        for (int i = 0; i < 40; ++i) {
            switcher.step(0.0, y, 1.0e-3, next, error);
        }
        assertTrue("forty attempts asked more than once, not " + switcher.trials(),
                switcher.trials() >= 2L);
        // forty calls, forty steps: a trial on this side is the scout's, and the
        // scout advances nothing, so it appears in neither count
        assertEquals("and every attempt is counted", 40L,
                switcher.nonStiffSteps() + switcher.stiffSteps());
        assertTrue("what the asking cost is its own figure",
                switcher.explorationEvaluations() > 0L);
    }

    /**
     * At a loose tolerance the local error test will accept a step taken well
     * outside the stability region, so a run that waited out the ordinary
     * cadence before asking would be destroyed before it ever asked. Robertson
     * at {@code rtol = 1e-4} accepts its first step at a stiffness measure of
     * {@code 6.6} against a threshold of {@code 3.25}, and twenty attempts later
     * the second concentration stands at {@code -0.12}.
     */
    @Test
    public void aLooseToleranceDoesNotLoseTheRunBeforeTheFirstTrial() {
        StepController controller = new StepController(1.0e-4, 1.0e-4);
        SwitchingStepper switcher = new SwitchingStepper(ROBERTSON, 3, controller);
        OdeIntegrator.Result mixed = new OdeIntegrator(switcher, controller)
                .solve(0.0, new double[] { 1.0, 0.0, 0.0 }, 1.0e5);

        assertEquals("it reaches the end", 1.0e5, mixed.finalTime(), 0.0);
        assertEquals("having handed over once", 1L, switcher.switches());
        assertTrue("and within the first handful of steps, not after twenty",
                switcher.nonStiffSteps() < 12L);
        // the concentrations stay concentrations, to within the tolerance the
        // run was asked for -- the measured worst excursion is a tenth of it,
        // and the failure this guards against left one of them at -0.12
        for (int i = 0; i < mixed.length; ++i) {
            for (int m = 0; m < 3; ++m) {
                assertTrue("y[" + i + "][" + m + "] is a concentration, not " + mixed.y[i][m],
                        mixed.y[i][m] > -1.0e-4 && mixed.y[i][m] < 1.0 + 1.0e-4);
            }
        }
    }

    /**
     * Acting on the free estimate has to be an edge and not a level: a run that
     * stays above the threshold while the trial keeps answering "stay" would
     * otherwise ask on every single step, which was measured at 39 trials in 40.
     */
    @Test
    public void aRunHeldAboveTheThresholdDoesNotProbeEveryStep() {
        SwitchingStepper switcher = new SwitchingStepper(ButcherTableau.DORMAND_PRINCE_45,
                RosenbrockTableau.RODAS4, DECAY, 1, new StepController(1.0e-6, 1.0e-6), 10);
        double[] y = { 1.0 };
        double[] next = new double[1];
        double[] error = new double[1];
        // |h lambda| is 10 at every one of these, so the measure never drops
        // back under the threshold and the edge never comes round again
        for (int i = 0; i < 40; ++i) {
            switcher.step(0.0, y, 1.0e-3, next, error);
        }
        assertTrue("the alarm fires once and the cadence carries the rest, not " + switcher.trials(),
                switcher.trials() <= 6L);
        assertTrue("but it does fire", switcher.trials() >= 1L);
    }

    /**
     * The free estimate is kept out of the decision but not out of the way: it
     * is what stops an equation with no stiffness in it from paying for trials
     * that could only ever come back saying "stay".
     */
    @Test
    public void anEquationWithNoStiffnessTakesNoTrialAtAll() {
        StepController controller = new StepController(1.0e-8, 1.0e-8);
        SwitchingStepper switcher = new SwitchingStepper(OSCILLATOR, 2, controller);
        new OdeIntegrator(switcher, controller).solve(0.0, new double[] { 1.0, 0.0 }, 50.0);
        assertEquals("trials", 0L, switcher.trials());
        assertEquals("implicit steps", 0L, switcher.stiffSteps());
        assertEquals("and nothing spent exploring", 0L, switcher.explorationEvaluations());

        // where it is stiff the veto lifts and the trials happen
        StepController other = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper stiffly = new SwitchingStepper(MOVING_DIAL, 1, other);
        new OdeIntegrator(stiffly, other).solve(0.0, new double[] { 1.0 }, 30.0);
        assertTrue("trials on the dial, " + stiffly.trials(), stiffly.trials() > 0L);
    }

    /**
     * The counts belong to one run and the evaluations to the stepper, which is
     * the split {@link OdeStepper#reset()} describes.
     */
    @Test
    public void resetReturnsToTheExplicitSideAndKeepsTheEvaluations() {
        StepController controller = new StepController(1.0e-8, 1.0e-8);
        SwitchingStepper switcher = new SwitchingStepper(ROBERTSON, 3, controller);
        new OdeIntegrator(switcher, controller).solve(0.0, new double[] { 1.0, 0.0, 0.0 }, 1.0e5);
        assertTrue("it went stiff", switcher.isStiffActive());
        long spent = switcher.evaluations();
        assertTrue("which cost something", spent > 0L);
        assertEquals("the two sides and the exploring together",
                switcher.nonStiffStepper().evaluations() + switcher.stiffStepper().evaluations()
                        + switcher.explorationEvaluations(), spent);

        switcher.reset();
        assertFalse("back on the explicit side", switcher.isStiffActive());
        assertEquals("switches forgotten", 0L, switcher.switches());
        assertEquals("explicit steps forgotten", 0L, switcher.nonStiffSteps());
        assertEquals("implicit steps forgotten", 0L, switcher.stiffSteps());
        assertEquals("evaluations kept", spent, switcher.evaluations());
    }

    /**
     * The threshold is the explicit method's own reach along the negative real
     * axis, and it does not move when the implicit side is stepping, because
     * {@link OdeIntegrator} reads it once at the start of a run.
     */
    @Test
    public void theThresholdIsTheExplicitMethodsOwn() {
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper five = new SwitchingStepper(DECAY, 1, controller);
        assertEquals("Dormand-Prince", 3.25, five.stiffnessThreshold(), 0.0);

        SwitchingStepper eight = new SwitchingStepper(ButcherTableau.DOP853, RosenbrockTableau.RODAS4,
                DECAY, 1, controller, 3);
        assertEquals("DOP853 reaches further", ButcherTableau.DOP853.stiffnessThreshold(),
                eight.stiffnessThreshold(), 0.0);
        assertEquals("and is an eighth order method", 8, eight.order());
        assertTrue("with two estimates of its own error", eight.hasErrorEstimate());
    }

    /**
     * The trial that decides the way back is judged by tolerances, and those
     * have to be the ones the driver accepts by.
     */
    @Test
    public void theIntegratorInsistsOnTheOneController() {
        StepController mine = new StepController(1.0e-6, 1.0e-6);
        SwitchingStepper switcher = new SwitchingStepper(DECAY, 1, mine);
        assertTrue("the stepper says which one", switcher.requiredController() == mine);

        new OdeIntegrator(switcher, mine);
        // a driver that judges nothing cannot disagree with anything
        new OdeIntegrator(switcher);
        try {
            new OdeIntegrator(switcher, new StepController(1.0e-6, 1.0e-6));
            fail("a second controller with the same settings is still a second controller");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    /** No other stepper is bound to a controller, and none should be. */
    @Test
    public void theOtherSteppersRequireNoController() {
        assertTrue(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2)
                .requiredController() == null);
        assertTrue(new Rosenbrock(RosenbrockTableau.RODAS4, DECAY, 1).requiredController() == null);
    }

    /**
     * A trial has to be judged, so both methods have to estimate the error of
     * their own steps. Nothing else is asked of them -- in particular, a method
     * that can form no stiffness estimate is fine, since none is used to decide
     * anything.
     */
    @Test
    public void bothMethodsHaveToEstimateTheirOwnError() {
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        try {
            new SwitchingStepper(ButcherTableau.CLASSIC_RK4, RosenbrockTableau.RODAS4, DECAY, 1,
                    controller, 3);
            fail("the classical method estimates no error at all");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }

        // an embedded pair with no two stages at the end of the step, so no
        // stiffness estimate ever: it works, and simply takes every trial
        double[] c = { 0.0, 1.0 };
        double[][] a = { {}, { 1.0 } };
        double[] b = { 0.5, 0.5 };
        double[] bStar = { 1.0, 0.0 };
        ButcherTableau heun = new ButcherTableau("Heun 2(1)", 2, c, a, b, bStar, 1, null);
        SwitchingStepper switcher = new SwitchingStepper(heun, RosenbrockTableau.RODAS4, DECAY, 1,
                controller, 5);
        double[] y = { 1.0 };
        double[] next = new double[1];
        double[] error = new double[1];
        switcher.step(0.0, y, 1.0e-4, next, error);
        assertEquals("it is the pair it was given", 2, switcher.order());
        assertTrue("with no stiffness to report", Double.isNaN(switcher.stiffnessMeasure()));
        assertEquals("and the default threshold", 3.25, switcher.stiffnessThreshold(), 0.0);
    }

    /** The rest of what the constructor will not accept. */
    @Test
    public void theConstructorChecksItsArguments() {
        StepController controller = new StepController(1.0e-6, 1.0e-6);
        try {
            new SwitchingStepper(DECAY, 1, null);
            fail("a stepper that judges trials needs a controller");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            new SwitchingStepper(ButcherTableau.DORMAND_PRINCE_45, RosenbrockTableau.RODAS4, DECAY, 1,
                    controller, 0);
            fail("a streak of zero steps hands over before it has looked");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            new SwitchingStepper(DECAY, 0, controller);
            fail("no components to advance");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            new SwitchingStepper((DVectorField) null, 1, controller);
            fail("no equation to solve");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    /** What it says about itself, which a caller reads in an exception. */
    @Test
    public void itNamesBothMethodsAndWhichIsStepping() {
        SwitchingStepper switcher = new SwitchingStepper(DECAY, 1,
                new StepController(1.0e-6, 1.0e-6));
        String said = switcher.toString();
        assertTrue(said, said.contains("Dormand-Prince"));
        assertTrue(said, said.contains("RODAS4"));
        assertTrue(said, said.contains("non-stiff"));
        assertEquals("the dimension it advances", 1, switcher.dimension());
    }
}
