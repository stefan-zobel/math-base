package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorField;

/**
 * Events, which are the reason the continuous extension is not an afterthought.
 * <p>
 * The question an event answers is not "what is the state at the times I
 * picked" but "at what time did this happen", and the answer is not on any
 * grid. A step whose ends straddle a zero brackets it, the interpolant turns
 * the step into a function of one variable, and
 * {@link math.solve.RootFinder#brentDekker} does the rest -- so the precision
 * of the time has nothing to do with the step size.
 * <p>
 * Measured against the analytic zeros of a cosine over {@code [0, 20]}, the six
 * crossings come back off by {@code 2.4e-06}, {@code 7.1e-10} and
 * {@code 1.1e-12} at solve tolerances of {@code 1e-6}, {@code 1e-9} and
 * {@code 1e-12}. That the numbers track the <em>solve</em> tolerance and not
 * the event tolerance is the point worth keeping: the root is found exactly on
 * a trajectory that is itself only approximate, so no event tolerance can
 * rescue a loose integration.
 * <p>
 * And watching is free where the method interpolates for free: the same run
 * with and without an event costs {@code 2979} evaluations either way.
 */
public final class OdeEventTest {

    /** y'' = -y, exact solution (cos t, -sin t) from (1, 0). */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    /** Height and vertical speed under gravity alone. */
    private static final DVectorField FALLING = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -9.81;
        }
    };

    /** Two bodies under gravity in the plane, as position then velocity. */
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

    /** The first component, whose zeros along the oscillator are the zeros of a cosine. */
    private static final OdeEvent FIRST = new OdeEvent() {
        @Override
        public double valueAt(double t, double[] y) {
            return y[0];
        }
    };

    @Test
    public void testTheZerosOfACosineAreFoundWhereTheAnalysisPutsThem() {
        OdeIntegrator.Result r = watching(1.0e-12, new Event(FIRST)).solve(0.0, start(), 20.0);
        assertEquals("six half periods fit into twenty", 6, r.eventCount());
        for (int i = 0; i < r.eventCount(); ++i) {
            assertEquals(Math.PI * (0.5 + i), r.eventTimes[i], 1.0e-10);
            assertEquals("the first component is zero there", 0.0, r.eventStates[i][0], 1.0e-11);
        }
        assertFalse(r.stoppedByEvent);
        assertEquals(20.0, r.finalTime(), 0.0);
    }

    /**
     * The root is exact on a trajectory that is not, so the time is only as
     * good as the integration underneath it. A tighter event tolerance on a
     * loose solve buys nothing, and the numbers below are the evidence.
     */
    @Test
    public void testTheEventTimeIsOnlyAsGoodAsTheTrajectoryItSitsOn() {
        double looser = 0.0;
        double[] tolerances = { 1.0e-6, 1.0e-9, 1.0e-12 };
        for (int q = 0; q < tolerances.length; ++q) {
            OdeIntegrator.Result r = watching(tolerances[q], new Event(FIRST)).solve(0.0, start(), 20.0);
            double worst = 0.0;
            for (int i = 0; i < r.eventCount(); ++i) {
                worst = Math.max(worst, Math.abs(r.eventTimes[i] - Math.PI * (0.5 + i)));
            }
            assertTrue("a tighter solve must place the crossing better", looser == 0.0 || worst < looser);
            assertTrue("at solve tolerance " + tolerances[q] + " the crossing was off by " + worst,
                    worst < 1000.0 * tolerances[q]);
            looser = worst;
        }
        assertTrue("and a loose solve cannot be rescued by the event tolerance", looser > 0.0);
    }

    @Test
    public void testTheDirectionDecidesWhichCrossingsCount() {
        OdeIntegrator.Result either = watching(1.0e-10, new Event(FIRST, Event.Direction.EITHER))
                .solve(0.0, start(), 20.0);
        OdeIntegrator.Result up = watching(1.0e-10, new Event(FIRST, Event.Direction.INCREASING))
                .solve(0.0, start(), 20.0);
        OdeIntegrator.Result down = watching(1.0e-10, new Event(FIRST, Event.Direction.DECREASING))
                .solve(0.0, start(), 20.0);

        assertEquals(6, either.eventCount());
        assertEquals(3, up.eventCount());
        assertEquals(3, down.eventCount());
        // a cosine falls through zero at pi/2 and rises through it at 3 pi/2
        for (int i = 0; i < 3; ++i) {
            assertEquals(Math.PI * (0.5 + 2 * i), down.eventTimes[i], 1.0e-8);
            assertEquals(Math.PI * (1.5 + 2 * i), up.eventTimes[i], 1.0e-8);
            assertTrue("the state rises through zero", up.eventStates[i][1] > 0.0);
            assertTrue("and falls through it", down.eventStates[i][1] < 0.0);
        }
    }

    /**
     * A ball thrown up at 20 metres a second lands after {@code 2 v / g}, and
     * the run has no business carrying on to the {@code t = 100} it was given.
     */
    @Test
    public void testATerminalEventEndsTheRunWhereItHappened() {
        Event landing = Event.terminal(FIRST, Event.Direction.DECREASING);
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, FALLING, 2),
                new StepController(1.0e-12, 1.0e-12), new Event[] { landing })
                        .solve(0.0, new double[] { 0.0, 20.0 }, 100.0);

        assertTrue(r.stoppedByEvent);
        assertEquals(1, r.eventCount());
        assertEquals(2.0 * 20.0 / 9.81, r.eventTimes[0], 1.0e-10);
        assertEquals("the run ends at the event and not at t1", r.eventTimes[0], r.finalTime(), 0.0);
        assertEquals("and the last state is the one it stopped at", 0.0, r.finalState()[0], 1.0e-9);
        assertTrue("nowhere near the hundred it was asked for", r.finalTime() < 5.0);
    }

    /**
     * A zero that lands exactly on the end of a step is that end, not a root
     * found a few units in the last place away from it -- and it fires once,
     * not again on the step that follows.
     */
    @Test
    public void testAnEventOnAStepBoundaryIsThatBoundaryExactly() {
        Event atHalf = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return t - 0.5;
            }
        });
        OdeIntegrator in = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2), null,
                new Event[] { atHalf });
        OdeIntegrator.Result r = in.solve(0.0, start(), 1.0, 10);
        assertEquals("once, not once on each side of the boundary", 1, r.eventCount());
        assertEquals(0.5, r.eventTimes[0], 0.0);
        assertEquals(Math.cos(0.5), r.eventStates[0][0], 1.0e-8);
        assertEquals(-Math.sin(0.5), r.eventStates[0][1], 1.0e-8);
    }

    /**
     * Two things an event is not. A quantity that is already zero when the run
     * starts has not crossed anything, and a quantity that reaches zero and
     * turns back has not either -- {@code cos t - 1} has a double root at
     * {@code 2 pi} and no sign change there, so nothing is reported, which is
     * the right answer and not a missed one.
     */
    @Test
    public void testStartingAtZeroIsNotAnEventAndNeitherIsTouchingIt() {
        Event throughZero = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[1];
            }
        });
        // -sin t is zero at the start, falls away, and crosses at pi and 2 pi
        OdeIntegrator.Result crossings = watching(1.0e-10, throughZero).solve(0.0, start(), 7.0);
        assertEquals("the departure from zero at t = 0 is not one of them", 2, crossings.eventCount());
        assertEquals(Math.PI, crossings.eventTimes[0], 1.0e-8);
        assertEquals(2.0 * Math.PI, crossings.eventTimes[1], 1.0e-8);

        Event touching = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0] - 1.0;
            }
        });
        OdeIntegrator.Result none = watching(1.0e-10, touching).solve(0.0, start(), 7.0);
        assertEquals("a double root is not a crossing", 0, none.eventCount());
    }

    @Test
    public void testAnEventIsOnlyWatchedAsOftenAsItMayHappen() {
        Event three = new Event(FIRST, Event.Direction.EITHER, false, Event.DEFAULT_TOLERANCE, 3);
        OdeIntegrator.Result r = watching(1.0e-10, three).solve(0.0, start(), 20.0);
        assertEquals(3, r.eventCount());
        assertEquals("and the run carries on to the end", 20.0, r.finalTime(), 0.0);
        for (int i = 0; i < 3; ++i) {
            assertEquals(Math.PI * (0.5 + i), r.eventTimes[i], 1.0e-8);
        }
    }

    /**
     * Two conditions on one run, and the orbit says what the answers should be:
     * the radial velocity is zero at apoapsis and periapsis, whose distances
     * from the focus are {@code 1 + e} and {@code 1 - e} exactly.
     */
    @Test
    public void testTwoEventsAreReportedInTheOrderTheyHappenedAndSayWhichIsWhich() {
        double e = 0.6;
        Event radial = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0] * y[2] + y[1] * y[3];
            }
        });
        Event crossingTheAxis = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[1];
            }
        });
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, KEPLER, 4),
                new StepController(1.0e-11, 1.0e-11), new Event[] { radial, crossingTheAxis })
                        .solve(0.0, kepler(e), 3.0 * 2.0 * Math.PI);

        assertTrue("both conditions fire repeatedly over three orbits", r.eventCount() >= 8);
        for (int i = 1; i < r.eventCount(); ++i) {
            assertTrue("the occurrences come in the order they happened",
                    r.eventTimes[i] >= r.eventTimes[i - 1]);
        }
        int apoapsides = 0;
        for (int i = 0; i < r.eventCount(); ++i) {
            assertTrue(r.eventIndices[i] == 0 || r.eventIndices[i] == 1);
            double distance = Math.hypot(r.eventStates[i][0], r.eventStates[i][1]);
            boolean far = Math.abs(distance - (1.0 + e)) < 1.0e-7;
            boolean near = Math.abs(distance - (1.0 - e)) < 1.0e-7;
            assertTrue("every turning point is at one apsis or the other, not at " + distance, far || near);
            if (far) {
                ++apoapsides;
            }
        }
        assertTrue("and the far one is reached about as often as the near one", apoapsides > 0);
    }

    /**
     * Watching costs nothing where the method carries its own interpolant, and
     * one evaluation per bracketing step where it does not -- the classical
     * method has to work out the derivative at the end of the step before it
     * can interpolate inside it.
     */
    @Test
    public void testWatchingIsFreeWhereTheMethodInterpolatesForFree() {
        OdeIntegrator.Result blind = watching(1.0e-10).solve(0.0, start(), 20.0);
        OdeIntegrator.Result seeing = watching(1.0e-10, new Event(FIRST)).solve(0.0, start(), 20.0);
        assertEquals(blind.steps, seeing.steps);
        assertEquals("Dormand-Prince interpolates over stages it already has", blind.evaluations,
                seeing.evaluations);

        OdeIntegrator plainRk = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2));
        OdeIntegrator watchfulRk = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2), null,
                new Event[] { new Event(FIRST) });
        long without = plainRk.solve(0.0, start(), 20.0, 200).evaluations;
        OdeIntegrator.Result withEvents = watchfulRk.solve(0.0, start(), 20.0, 200);
        assertEquals(800L, without);
        assertEquals("one per step that had to be searched", without + withEvents.eventCount(),
                withEvents.evaluations);
    }

    @Test
    public void testARunWithoutEventsReportsNone() {
        OdeIntegrator.Result r = watching(1.0e-9).solve(0.0, start(), 5.0);
        assertEquals(0, r.eventCount());
        assertEquals(0, r.eventTimes.length);
        assertEquals(0, r.eventIndices.length);
        assertEquals(0, r.eventStates.length);
        assertFalse(r.stoppedByEvent);
        assertTrue(r.toString().contains("0 events"));
    }

    @Test
    public void testTheEventArrayIsCopied() {
        Event[] given = { new Event(FIRST) };
        OdeIntegrator in = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2), null, given);
        given[0] = null;
        assertEquals(1, in.events().length);
        in.events()[0] = null;
        assertEquals(1, in.events().length);
        assertTrue(in.events()[0].toString().contains("EITHER"));
    }

    @Test
    public void testTheEventsAreChecked() {
        try {
            new Event(null);
            fail("expected a refusal naming the function");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("function"));
        }
        try {
            new Event(FIRST, null);
            fail("expected a refusal naming the direction");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("direction"));
        }
        try {
            new Event(FIRST, Event.Direction.EITHER, false, 0.0, 1);
            fail("expected a refusal naming the tolerance");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("tolerance"));
        }
        try {
            new Event(FIRST, Event.Direction.EITHER, false, 1.0e-9, 0);
            fail("expected a refusal naming the count");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("maxCount"));
        }
        try {
            new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2), null,
                    new Event[] { null });
            fail("expected a refusal naming the entry");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("events[0]"));
        }
    }

    @Test
    public void testWhatAnEventSaysAboutItself() {
        Event plain = new Event(FIRST);
        assertEquals(Event.Direction.EITHER, plain.direction());
        assertFalse(plain.isTerminal());
        assertEquals(Event.DEFAULT_TOLERANCE, plain.tolerance(), 0.0);
        assertEquals(Integer.MAX_VALUE, plain.maxCount());
        assertEquals(FIRST, plain.function());

        Event stopper = Event.terminal(FIRST, Event.Direction.INCREASING);
        assertTrue(stopper.isTerminal());
        assertEquals(1, stopper.maxCount());
        assertEquals(Event.Direction.INCREASING, stopper.direction());
        assertTrue(stopper.toString().contains("terminal"));
    }

    private static double[] start() {
        return new double[] { 1.0, 0.0 };
    }

    private static double[] kepler(double eccentricity) {
        return new double[] { 1.0 - eccentricity, 0.0, 0.0,
                Math.sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) };
    }

    private static OdeIntegrator watching(double tolerance, Event... events) {
        return new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2),
                new StepController(tolerance, tolerance), events);
    }
}
