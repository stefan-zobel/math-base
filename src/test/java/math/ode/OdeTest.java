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
     * to keep passing. Its six lines are the same claims the test classes make
     * at length, which is the point: a reader who runs the class gets the
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
        assertEquals("six claims and a verdict", 7, page.split("\\r?\\n").length);
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
