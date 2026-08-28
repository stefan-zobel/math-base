package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.Locale;

import org.junit.Test;

import math.fun.DSecondOrderField;
import math.fun.DVectorField;
import math.fun.DiffDVectorField;

/**
 * The facade, checked against the three objects it builds.
 * <p>
 * A shortcut is only worth having if it is the same thing written shorter, so
 * what is asserted here is mostly identity: {@link Ode#solve} and the
 * {@link OdeIntegrator} a caller would have assembled agree bit for bit, and
 * the form without a tolerance is the form with {@link Ode#DEFAULT_TOLERANCE}.
 * The one place the facade decides something is
 * {@link Ode#solveFixed(math.fun.DVectorField, double, double[], double, int)},
 * which takes the classical method rather than Dormand-Prince, and its four
 * evaluations per step are what says so.
 */
public final class OdeTest {

    /** y'' = -y, exact solution (cos t, -sin t) from (1, 0). */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    /** y' = 1, which crosses nothing and reaches the far end. */
    private static final DVectorField STEADY = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = 1.0;
        }
    };

    /** The same two bodies as a first order system, position then velocity. */
    private static final DVectorField KEPLER = new DVectorField() {
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

    /** Two bodies under gravity, as an acceleration on the position alone. */
    private static final DSecondOrderField GRAVITY = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acceleration) {
            double r2 = q[0] * q[0] + q[1] * q[1];
            double r3 = r2 * Math.sqrt(r2);
            acceleration[0] = -q[0] / r3;
            acceleration[1] = -q[1] / r3;
        }
    };

    @Test
    public void testTheFacadeIsTheObjectsItBuilds() {
        OdeIntegrator.Result byHand = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2),
                new StepController(1.0e-9, 1.0e-9)).solve(0.0, start(), 12.0);
        OdeIntegrator.Result byFacade = Ode.solve(OSCILLATOR, 0.0, start(), 12.0, 1.0e-9);

        assertEquals(byHand.length, byFacade.length);
        assertEquals(byHand.steps, byFacade.steps);
        assertEquals(byHand.evaluations, byFacade.evaluations);
        for (int i = 0; i < byHand.length; ++i) {
            assertEquals(byHand.t[i], byFacade.t[i], 0.0);
            assertEquals(byHand.y[i][0], byFacade.y[i][0], 0.0);
            assertEquals(byHand.y[i][1], byFacade.y[i][1], 0.0);
        }
    }

    @Test
    public void testTheFormWithoutAToleranceIsTheFormWithTheDefaultOne() {
        assertEquals(1.0e-8, Ode.DEFAULT_TOLERANCE, 0.0);
        OdeIntegrator.Result implied = Ode.solve(OSCILLATOR, 0.0, start(), 12.0);
        OdeIntegrator.Result spelled = Ode.solve(OSCILLATOR, 0.0, start(), 12.0, Ode.DEFAULT_TOLERANCE);
        assertEquals(implied.length, spelled.length);
        assertEquals(implied.evaluations, spelled.evaluations);
        assertEquals(implied.finalState()[0], spelled.finalState()[0], 0.0);
    }

    @Test
    public void testTheEndpointIsTheLastStateAndNothingElse() {
        double[] end = Ode.endpoint(OSCILLATOR, 0.0, start(), 12.0, 1.0e-10);
        assertEquals(2, end.length);
        assertEquals(Math.cos(12.0), end[0], 1.0e-8);
        assertEquals(-Math.sin(12.0), end[1], 1.0e-8);
        OdeIntegrator.Result full = Ode.solve(OSCILLATOR, 0.0, start(), 12.0, 1.0e-10);
        assertEquals(full.finalState()[0], end[0], 0.0);
        assertEquals(full.finalState()[1], end[1], 0.0);
    }

    @Test
    public void testTheFixedFormTakesTheClassicalMethod() {
        OdeIntegrator.Result r = Ode.solveFixed(OSCILLATOR, 0.0, start(), 10.0, 100);
        assertEquals(101, r.length);
        assertEquals(100L, r.steps);
        assertEquals("four stages, none of them carried over", 400L, r.evaluations);
        assertEquals(0L, r.rejected);
        assertEquals("fourth order at a step of 0.1 over ten units", Math.cos(r.finalTime()),
                r.finalState()[0], 1.0e-5);
    }

    @Test
    public void testSolvingUntilSomethingHappensStopsThere() {
        OdeIntegrator.Result r = Ode.solveUntil(OSCILLATOR, 0.0, start(), 20.0, new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0];
            }
        }, 1.0e-12);
        assertTrue(r.stoppedByEvent);
        assertEquals(1, r.eventCount());
        assertEquals(0.5 * Math.PI, r.eventTimes[0], 1.0e-10);
        assertEquals("the run ends where the event did", r.eventTimes[0], r.finalTime(), 0.0);
    }

    @Test
    public void testSolvingUntilSomethingThatNeverHappensReachesTheEnd() {
        OdeIntegrator.Result r = Ode.solveUntil(STEADY, 0.0, new double[] { 1.0 }, 5.0, new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0];
            }
        }, 1.0e-10);
        assertFalse(r.stoppedByEvent);
        assertEquals(0, r.eventCount());
        assertEquals(5.0, r.finalTime(), 0.0);
        assertEquals(6.0, r.finalState()[0], 1.0e-9);
    }

    @Test
    public void testTheGridFormReportsExactlyTheTimesItWasGiven() {
        double[] times = { 0.0, 1.0, 2.5, 3.0, 7.0 };
        OdeIntegrator.Result r = Ode.solve(OSCILLATOR, 0.0, start(), 7.0, times, 1.0e-10);
        assertEquals(times.length, r.length);
        for (int i = 0; i < times.length; ++i) {
            assertEquals(times[i], r.t[i], 0.0);
            assertEquals(Math.cos(times[i]), r.y[i][0], 1.0e-9);
        }
    }

    @Test
    public void testTheFacadeChecksItsArguments() {
        try {
            Ode.solve(OSCILLATOR, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
        try {
            Ode.solveFixed(OSCILLATOR, 0.0, null, 1.0, 10);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
        try {
            Ode.solve(OSCILLATOR, 0.0, start(), 1.0, -1.0);
            fail("expected a refusal naming the tolerance");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("rtol"));
        }
    }

    /**
     * The self check is the class's own statement about the package, and it has
     * to keep passing. Its eight lines are the same claims the test classes
     * make at length, which is the point: a reader who runs the class gets the
     * measurements without reading the tests.
     */
    @Test
    public void testTheSelfCheckStillPasses() throws UnsupportedEncodingException {
        String page = capture();
        assertTrue(page, page.contains(">>> OK"));
        assertFalse(page, page.contains("FAILED"));
        assertTrue(page, page.contains("against Gauss-Kronrod"));
        assertTrue(page, page.contains("drifts, and is meant to"));
        assertTrue("the contrast is the point of the two Kepler lines", page.contains("symplectically"));
        assertTrue("and the stiff problem answers the line above it",
                page.contains("the explicit method gives up"));
        assertTrue("and the tolerance question is there",
                page.contains("at rtol 1e-12"));
        assertTrue("and the last line is the one where the regime changes mid-run",
                page.contains("free on a smooth problem"));
        assertEquals("nine claims and a verdict", 10, page.split("\\r?\\n").length);
    }

    /**
     * The two Kepler lines are next to each other on purpose: the same orbit,
     * the same hundred revolutions, one method whose energy error grows and one
     * whose energy error does not.
     */
    @Test
    public void testTheSelfCheckPutsTheTwoKeplerRunsSideBySide() throws UnsupportedEncodingException {
        String[] lines = capture().split("\\r?\\n");
        int drifting = -1;
        int bounded = -1;
        for (int i = 0; i < lines.length; ++i) {
            if (lines[i].contains("drifts, and is meant to")) {
                drifting = i;
            }
            if (lines[i].contains("symplectically")) {
                bounded = i;
            }
        }
        assertTrue("both lines are there", drifting >= 0 && bounded >= 0);
        assertEquals("and one follows the other", drifting + 1, bounded);
    }

    @Test
    public void testTheSymplecticEntryTakesThePositionAndTheVelocitySeparately() {
        double[] q0 = { 0.4, 0.0 };
        double[] v0 = { 0.0, 2.0 };
        OdeIntegrator.Result r = Ode.solveSymplectic(GRAVITY, 0.0, q0, v0, 2.0 * Math.PI, 400);
        assertEquals(401, r.length);
        assertEquals(400L, r.steps);
        assertEquals("eleven evaluations a step, and the first one pays for the shared kick", 4401L,
                r.evaluations);
        assertEquals("the state is the position and then the velocity", 4, r.finalState().length);
        // the orbit has semi-major axis one and so a period of 2 pi exactly
        assertEquals("one whole revolution returns to the start", q0[0], r.finalState()[0], 1.0e-11);
        assertEquals(q0[1], r.finalState()[1], 1.0e-11);
        assertEquals("and the inputs are not written into", 0.4, q0[0], 0.0);
        assertEquals(2.0, v0[1], 0.0);
    }

    /**
     * {@code y' = -L (y - cos t) - sin t}, whose solution is {@code cos t} for
     * every {@code L} and which an explicit method cannot follow at all once
     * {@code L} is large.
     */
    private static final DiffDVectorField STIFF = new DiffDVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = -1.0e5 * (y[0] - Math.cos(t)) - Math.sin(t);
        }

        @Override
        public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
            dfdy[0] = -1.0e5;
            dfdt[0] = -1.0e5 * Math.sin(t) - Math.cos(t);
        }
    };

    @Test
    public void testTheStiffEntryIsRodas4AndTheObjectsAroundIt() {
        OdeIntegrator.Result byHand = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, (DVectorField) STIFF, 1),
                new StepController(1.0e-8, 1.0e-8)).solve(0.0, new double[] { 1.0 }, 10.0);
        OdeIntegrator.Result byFacade = Ode.solveStiff((DVectorField) STIFF, 0.0, new double[] { 1.0 },
                10.0, 1.0e-8);
        assertEquals(byHand.length, byFacade.length);
        assertEquals(byHand.evaluations, byFacade.evaluations);
        assertEquals(byHand.finalState()[0], byFacade.finalState()[0], 0.0);
        assertEquals("and the answer is the cosine", Math.cos(10.0), byFacade.finalState()[0], 1.0e-6);
    }

    /**
     * The same run with the Jacobian written down. The two do not take exactly
     * the same steps -- a differenced Jacobian is an approximate one and the
     * error estimates part company a little -- but the cost of each is exact:
     * six stages a step, five on a step retried from a point already
     * linearized, and two more per point where the derivatives have to be
     * differenced out of a system of one component.
     */
    @Test
    public void testTheStiffEntryTakesAJacobianWhenOneIsOffered() {
        OdeIntegrator.Result differenced = Ode.solveStiff((DVectorField) STIFF, 0.0,
                new double[] { 1.0 }, 10.0, 1.0e-8);
        OdeIntegrator.Result written = Ode.solveStiff(STIFF, 0.0, new double[] { 1.0 }, 10.0, 1.0e-8);
        assertEquals("the two agree on the answer", differenced.finalState()[0],
                written.finalState()[0], 1.0e-9);
        assertTrue("and very nearly on the number of steps: " + written.steps + " against "
                + differenced.steps, Math.abs(written.steps - differenced.steps) < 0.05 * written.steps);
        assertEquals("six a step, five on a retry, and two for the first step",
                2L + 6L * written.steps + 5L * written.rejected, written.evaluations);
        assertEquals("and two more per point where the Jacobian is differenced",
                2L + 8L * differenced.steps + 5L * differenced.rejected, differenced.evaluations);
        assertTrue("so the written form is the cheaper one: " + written.evaluations + " against "
                + differenced.evaluations, written.evaluations < differenced.evaluations);
    }

    /**
     * What the explicit entry says about the same problem, and what happens
     * when its advice is taken. This is the pair the package comment is about.
     */
    @Test
    public void testTheStiffnessReportIsTheAdviceToUseTheStiffEntry() {
        OdeIntegrator.Result explicit = Ode.solve(STIFF, 0.0, new double[] { 1.0 }, 1.0, 1.0e-8);
        assertTrue("the explicit method should notice", explicit.seemsStiff);
        OdeIntegrator.Result implicit = Ode.solveStiff(STIFF, 0.0, new double[] { 1.0 }, 1.0, 1.0e-8);
        assertFalse("and the implicit one has no such limit to report", implicit.seemsStiff);
        assertTrue("and it should be very much cheaper: " + implicit.evaluations + " against "
                + explicit.evaluations, 10L * implicit.evaluations < explicit.evaluations);
        assertEquals(Math.cos(1.0), implicit.finalState()[0], 1.0e-7);
        assertEquals(Math.cos(1.0), explicit.finalState()[0], 1.0e-7);
    }

    @Test
    public void testTheStiffEntryChecksItsArguments() {
        try {
            Ode.solveStiff((DVectorField) STIFF, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
        try {
            Ode.solveStiff(STIFF, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
    }

    /**
     * The accurate entry is DOP853 and reaches the same answer as the ordinary
     * one, for a good deal less at a tight tolerance: 3099 evaluations against
     * 9891 over ten orbits at {@code rtol = 1e-12}.
     */
    @Test
    public void testTheAccurateEntryIsDop853AndIsCheaperWhereItShouldBe() {
        double[] start = { 0.4, 0.0, 0.0, 2.0 };
        OdeIntegrator.Result five = Ode.solve(KEPLER, 0.0, start, 20.0, 1.0e-12);
        OdeIntegrator.Result eight = Ode.solveAccurate(KEPLER, 0.0, start, 20.0, 1.0e-12);
        for (int i = 0; i < start.length; ++i) {
            assertEquals("component " + i, five.finalState()[i], eight.finalState()[i], 1.0e-8);
        }
        assertTrue("at rtol 1e-12 the eighth order method should cost less than half: "
                + eight.evaluations + " against " + five.evaluations,
                eight.evaluations * 2L < five.evaluations);
        assertFalse("a smooth orbit is not stiff at any tolerance", eight.seemsStiff);
    }

    /**
     * The automatic entry is the switching stepper, and the shortcut is the
     * same thing written shorter -- including that both objects have to be
     * handed the one controller, which is what the shortcut is chiefly for.
     */
    @Test
    public void testTheAutomaticEntryIsTheSwitchingStepper() {
        StepController controller = new StepController(1.0e-8, 1.0e-8);
        OdeIntegrator.Result byHand = new OdeIntegrator(
                new SwitchingStepper((DVectorField) STIFF, 1, controller), controller)
                        .solve(0.0, new double[] { 1.0 }, 10.0);
        OdeIntegrator.Result byFacade = Ode.solveAuto((DVectorField) STIFF, 0.0, new double[] { 1.0 },
                10.0, 1.0e-8);
        assertEquals(byHand.length, byFacade.length);
        assertEquals(byHand.evaluations, byFacade.evaluations);
        assertEquals(byHand.finalState()[0], byFacade.finalState()[0], 0.0);
        assertEquals("and the answer is the cosine", Math.cos(10.0), byFacade.finalState()[0], 1.0e-6);

        // the written Jacobian reaches the same answer for less. The two are
        // not asserted against each other beyond the accuracy either of them
        // claims: a differenced Jacobian is an approximate one, the two runs
        // therefore change methods at slightly different points, and that gap
        // is a property of the answer rather than of the switching
        OdeIntegrator.Result written = Ode.solveAuto(STIFF, 0.0, new double[] { 1.0 }, 10.0, 1.0e-8);
        assertEquals("the written form is the cosine too", Math.cos(10.0), written.finalState()[0],
                1.0e-6);
        assertTrue("and it is cheaper: " + written.evaluations + " against " + byFacade.evaluations,
                written.evaluations < byFacade.evaluations);
    }

    /**
     * Where there is no stiffness there is nothing to switch, and the automatic
     * entry costs exactly what the ordinary one costs -- not nearly, exactly,
     * because it takes no trial and the arithmetic is then the same arithmetic.
     */
    @Test
    public void testTheAutomaticEntryIsFreeOnASmoothProblem() {
        double[] start = { 1.0, 0.0 };
        OdeIntegrator.Result plain = Ode.solve(OSCILLATOR, 0.0, start, 50.0, 1.0e-10);
        OdeIntegrator.Result auto = Ode.solveAuto(OSCILLATOR, 0.0, start, 50.0, 1.0e-10);
        assertEquals("evaluations", plain.evaluations, auto.evaluations);
        assertEquals("steps", plain.steps, auto.steps);
        assertEquals("recorded points", plain.length, auto.length);
        for (int i = 0; i < plain.length; ++i) {
            assertEquals("t[" + i + "]", plain.t[i], auto.t[i], 0.0);
            assertEquals("y[" + i + "][0]", plain.y[i][0], auto.y[i][0], 0.0);
            assertEquals("y[" + i + "][1]", plain.y[i][1], auto.y[i][1], 0.0);
        }
    }

    @Test
    public void testTheAutomaticEntryChecksItsArguments() {
        try {
            Ode.solveAuto((DVectorField) null, 0.0, new double[] { 1.0 }, 1.0, 1.0e-8);
            fail("expected a refusal naming the field");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().length() > 0);
        }
        try {
            Ode.solveAuto((DVectorField) STIFF, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
        try {
            Ode.solveAuto(STIFF, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
    }

    @Test
    public void testTheAccurateEntryChecksItsArguments() {
        try {
            Ode.solveAccurate(KEPLER, 0.0, null, 1.0, 1.0e-8);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
    }

    @Test
    public void testTheSymplecticEntryChecksItsArguments() {
        try {
            Ode.solveSymplectic(GRAVITY, 0.0, null, new double[2], 1.0, 10);
            fail("expected a refusal naming q0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("q0"));
        }
        try {
            Ode.solveSymplectic(GRAVITY, 0.0, new double[2], new double[3], 1.0, 10);
            fail("expected a refusal about the lengths");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("must agree"));
        }
    }

    @Test
    public void testTheSelfCheckReadsTheSameInAnyLocale() throws UnsupportedEncodingException {
        Locale before = Locale.getDefault();
        try {
            Locale.setDefault(Locale.US);
            String english = capture();
            Locale.setDefault(Locale.GERMANY);
            String german = capture();
            assertEquals(english, german);
        } finally {
            Locale.setDefault(before);
        }
    }

    private static String capture() throws UnsupportedEncodingException {
        PrintStream before = System.out;
        ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(sink, true, "UTF-8"));
            Ode.main(new String[0]);
        } finally {
            System.setOut(before);
        }
        return sink.toString("UTF-8");
    }

    private static double[] start() {
        return new double[] { 1.0, 0.0 };
    }
}
