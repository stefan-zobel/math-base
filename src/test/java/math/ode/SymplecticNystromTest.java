package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DSecondOrderField;
import math.fun.DVectorField;

/**
 * {@link SymplecticNystrom} against the two body problem, whose invariants are
 * known exactly and whose behavior over a long time is the reason the class
 * exists.
 * <p>
 * <b>The result the round was built for.</b> Over two hundred orbits at two
 * hundred steps each, the energy error of every method here occupies a band
 * whose width over the last twenty orbits is the width over the first twenty to
 * three decimal places -- a ratio of {@code 1.000}. Dormand-Prince, on the same
 * problem, has no band at all: its energy error is {@code 1.28e-09},
 * {@code 1.18e-08} and {@code 1.16e-07} over ten, a hundred and a thousand
 * orbits, a factor of ten per factor of ten, and no tolerance changes that
 * slope.
 * <p>
 * <b>And the result that keeps it honest.</b> Bounded energy is not the same as
 * a better answer. At about a thousand evaluations per orbit, the position after
 * ten orbits is off by {@code 1.9e-04} for Suzuki and by {@code 2.8e-07} for
 * Dormand-Prince -- the adaptive fifth order method is nearly three orders of
 * magnitude more accurate at the same cost. What changes with the horizon is the
 * *rate*: the symplectic error grows exactly linearly, ten times over ten times
 * the orbits, and Dormand-Prince's grows as the square, eighty times over the
 * same ten. Extrapolated, they meet at some ten thousand orbits, and past that
 * the symplectic method is ahead and stays ahead.
 * <p>
 * <b>The angular momentum is not approximately conserved, it is exactly
 * conserved.</b> A kick changes the velocity along {@code q} and a drift changes
 * the position along {@code v}, so neither can move {@code q x v}. Over forty
 * thousand steps it moves by {@code 2e-14}, which is the rounding and not the
 * method.
 * <p>
 * <b>Two of the five methods begin with a kick rather than a drift</b>, and
 * that kick falls where the previous step's last one did. Sharing the
 * evaluation between them is what brings Blanes-Moan down to the six and eleven
 * evaluations a step it is published with, and it is exact rather than an
 * approximation: over fifty steps the shared answer is bit for bit the
 * unshared one, at forty-nine fewer evaluations.
 */
public final class SymplecticNystromTest {

    /** q'' = -q, exact solution (cos t, -sin t) from (1, 0). */
    private static final DSecondOrderField OSCILLATOR = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acceleration) {
            acceleration[0] = -q[0];
        }
    };

    /** A damped oscillator, which is outside what a symplectic method promises. */
    private static final DSecondOrderField DAMPED = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acceleration) {
            acceleration[0] = -q[0] - 0.1 * v[0];
        }
    };

    /** Two bodies under gravity in the plane, as an acceleration on the position. */
    private static final DSecondOrderField KEPLER = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acceleration) {
            double r2 = q[0] * q[0] + q[1] * q[1];
            double r3 = r2 * Math.sqrt(r2);
            acceleration[0] = -q[0] / r3;
            acceleration[1] = -q[1] / r3;
        }
    };

    /** The same thing flattened, for the comparison against Dormand-Prince. */
    private static final DVectorField KEPLER_FLAT = new DVectorField() {
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

    private static final double PERIOD = 2.0 * Math.PI;

    private static SplittingCoefficients[] all() {
        return new SplittingCoefficients[] { SplittingCoefficients.VERLET, SplittingCoefficients.YOSHIDA_4,
                SplittingCoefficients.SUZUKI_4, SplittingCoefficients.BLANES_MOAN_4,
                SplittingCoefficients.BLANES_MOAN_6 };
    }

    /** What a run of {@code steps} steps costs: the first one pays the extra. */
    private static long cost(SplittingCoefficients c, int steps) {
        return (long) steps * c.evaluations() + (c.isVelocityFirst() ? 1L : 0L);
    }

    /**
     * The order over a whole interval, halving the step until the error would be
     * lost in the rounding -- which happens sooner the better the method is, so
     * the range cannot be the same for all five.
     */
    @Test
    public void testTheConvergenceOrderIsWhatTheCoefficientsClaim() {
        SplittingCoefficients[] methods = all();
        for (int m = 0; m < methods.length; ++m) {
            double previous = 0.0;
            int checked = 0;
            // a method of a higher order reaches the rounding floor at a
            // coarser step, so it has to be measured at coarser steps; one of a
            // lower order needs finer ones to leave its pre-asymptotic range
            int first = (methods[m].order() >= 6) ? 10 : 100;
            for (int steps = first; steps <= 100 * first; steps *= 2) {
                OdeIntegrator.Result r = new OdeIntegrator(new SymplecticNystrom(methods[m], OSCILLATOR, 1))
                        .solve(0.0, new double[] { 1.0, 0.0 }, 10.0, steps);
                double t = r.finalTime();
                double error = Math.hypot(r.finalState()[0] - Math.cos(t), r.finalState()[1] + Math.sin(t));
                assertEquals("the run costs what the coefficients say", cost(methods[m], steps),
                        r.evaluations);
                if (error < 1.0e-11) {
                    break;
                }
                if (previous != 0.0) {
                    assertEquals(methods[m].name() + ", halving " + (checked + 1), methods[m].order(),
                            Math.log(previous / error) / Math.log(2.0), 0.05);
                    ++checked;
                }
                previous = error;
            }
            assertTrue(methods[m].name() + " gave only " + checked + " usable halvings", checked >= 2);
        }
    }

    /**
     * The saving a method that begins with a kick gets, and the reason it is a
     * saving rather than an approximation: the kick that opens a step and the
     * one that closed the step before are at the same time and the same
     * position, so the field is evaluated once. Over fifty steps the answer is
     * <b>bit for bit</b> what it is without the sharing, at forty-nine fewer
     * evaluations.
     */
    @Test
    public void testTheSharedKickSavesOneEvaluationAndChangesNothing() {
        SplittingCoefficients[] velocityFirst = { SplittingCoefficients.BLANES_MOAN_4,
                SplittingCoefficients.BLANES_MOAN_6 };
        for (int m = 0; m < velocityFirst.length; ++m) {
            SplittingCoefficients c = velocityFirst[m];
            SymplecticNystrom sharing = new SymplecticNystrom(c, KEPLER, 2);
            SymplecticNystrom fresh = new SymplecticNystrom(c, KEPLER, 2);
            double[] y1 = kepler(0.6);
            double[] y2 = kepler(0.6);
            double[] out1 = new double[4];
            double[] out2 = new double[4];
            double t = 0.0;
            for (int i = 0; i < 50; ++i) {
                sharing.step(t, y1, 0.02, out1, null);
                // a reset throws the cached evaluation away, so this one pays
                // for every kick it has
                fresh.reset();
                fresh.step(t, y2, 0.02, out2, null);
                for (int j = 0; j < 4; ++j) {
                    assertEquals(c.name() + " at step " + i + ", component " + j, out2[j], out1[j], 0.0);
                }
                System.arraycopy(out1, 0, y1, 0, 4);
                System.arraycopy(out2, 0, y2, 0, 4);
                t += 0.02;
            }
            assertEquals(cost(c, 50), sharing.evaluations());
            assertEquals("without the sharing, one more per step", 50L * (c.evaluations() + 1),
                    fresh.evaluations());
        }
    }

    /**
     * The comparison the round was asked for: at the same accuracy on a real
     * orbit, how much less does the optimized method cost than the triple jump.
     */
    @Test
    public void testBlanesMoanReachesTheSameAccuracyForAFractionOfTheCost() {
        double end = 20.0 * PERIOD;
        double[] exact = kepler(0.6);
        // step counts chosen to put the three at about the same accuracy
        long yoshida = costToReach(SplittingCoefficients.YOSHIDA_4, end, exact, 1.0e-6);
        long suzuki = costToReach(SplittingCoefficients.SUZUKI_4, end, exact, 1.0e-6);
        long blanesMoan = costToReach(SplittingCoefficients.BLANES_MOAN_4, end, exact, 1.0e-6);
        assertTrue("Suzuki should beat the triple jump: " + suzuki + " against " + yoshida,
                suzuki < yoshida);
        assertTrue("and Blanes-Moan should beat them both: " + blanesMoan, blanesMoan < suzuki);
        assertTrue("by a wide margin against Yoshida, not " + (yoshida / (double) blanesMoan),
                yoshida > 3.0 * blanesMoan);
    }

    /**
     * A palindromic method run backwards over the step it just took returns to
     * where it started, and not approximately: the sequence of drifts and kicks
     * reversed is the sequence itself.
     */
    @Test
    public void testAStepIsUndoneExactlyByTheStepBack() {
        SplittingCoefficients[] methods = all();
        for (int m = 0; m < methods.length; ++m) {
            SymplecticNystrom stepper = new SymplecticNystrom(methods[m], KEPLER, 2);
            double[] y = kepler(0.6);
            double[] there = new double[4];
            double[] back = new double[4];
            stepper.step(0.0, y, 0.05, there, null);
            stepper.step(0.05, there, -0.05, back, null);
            for (int i = 0; i < 4; ++i) {
                assertEquals(methods[m].name() + ", component " + i, y[i], back[i], 1.0e-15);
            }
        }
    }

    /**
     * The defining property, on a problem that is not linear. A drift and a kick
     * each preserve phase space volume whatever their coefficients are, so the
     * product of them does; what is measured here is limited by the central
     * differences that produce the Jacobian, not by the method.
     */
    @Test
    public void testTheOneStepMapIsSymplecticOnANonlinearProblem() {
        SplittingCoefficients[] methods = all();
        for (int m = 0; m < methods.length; ++m) {
            double[][] jacobian = jacobianOfOneStep(new SymplecticNystrom(methods[m], KEPLER, 2), kepler(0.6),
                    0.05);
            assertEquals(methods[m].name(), 0.0, symplecticDefect(jacobian), 1.0e-8);
        }
    }

    /**
     * The headline: the energy error does not grow. Its band over the last
     * twenty of two hundred orbits is the band over the first twenty, and the
     * three methods differ only in how wide that band is.
     */
    @Test
    public void testTheEnergyErrorStaysInABandThatDoesNotGrow() {
        SplittingCoefficients[] methods = all();
        double[] widths = new double[methods.length];
        for (int m = 0; m < methods.length; ++m) {
            OdeIntegrator.Result r = new OdeIntegrator(new SymplecticNystrom(methods[m], KEPLER, 2))
                    .solve(0.0, kepler(0.6), 200 * PERIOD, 200 * 200);
            double first = energyBand(r, 0, 20 * 200);
            double last = energyBand(r, r.length - 20 * 200, r.length);
            assertEquals(methods[m].name() + ": the band must not widen", 1.0, last / first, 0.05);
            widths[m] = first;
        }
        assertTrue("a higher order gives a narrower band", widths[0] > widths[1]);
        assertTrue("and a smaller error constant a narrower one still", widths[1] > widths[2]);
    }

    /**
     * What Dormand-Prince does on the same problem, for the contrast: no band,
     * a drift proportional to the time it has been running.
     */
    @Test
    public void testWhereDormandPrinceHasNoBandAtAll() {
        double[] drift = new double[2];
        for (int i = 0; i < 2; ++i) {
            int orbits = (i == 0) ? 20 : 200;
            OdeIntegrator.Result r = new OdeIntegrator(
                    new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, KEPLER_FLAT, 4),
                    new StepController(1.0e-10, 1.0e-10, 400000)).solve(0.0, kepler(0.6), orbits * PERIOD);
            double reference = energy(kepler(0.6));
            double worst = 0.0;
            for (int j = 0; j < r.length; ++j) {
                worst = Math.max(worst, Math.abs(energy(r.y[j]) - reference));
            }
            drift[i] = worst;
        }
        assertTrue("ten times as long should drift about ten times as far, not " + drift[1] / drift[0],
                drift[1] / drift[0] > 7.0 && drift[1] / drift[0] < 13.0);
    }

    /**
     * The angular momentum of a central force problem is conserved by every
     * drift and every kick separately, so it is conserved exactly rather than
     * to the order of the method. It is therefore a sharper check than the
     * energy, and the only thing that moves it is rounding.
     */
    @Test
    public void testTheAngularMomentumIsConservedExactly() {
        SplittingCoefficients[] methods = all();
        for (int m = 0; m < methods.length; ++m) {
            OdeIntegrator.Result r = new OdeIntegrator(new SymplecticNystrom(methods[m], KEPLER, 2))
                    .solve(0.0, kepler(0.6), 200 * PERIOD, 200 * 200);
            double reference = momentum(kepler(0.6));
            double worst = 0.0;
            for (int j = 0; j < r.length; ++j) {
                worst = Math.max(worst, Math.abs(momentum(r.y[j]) - reference));
            }
            assertEquals(methods[m].name() + " over forty thousand steps", 0.0, worst, 1.0e-12);
        }
    }

    /**
     * The comparison a caller actually has to make. Bounded energy does not
     * mean a better position: at matched cost the adaptive fifth order method
     * is nearly three orders of magnitude ahead after ten orbits. What the
     * symplectic method buys is the <em>rate</em> -- linear against quadratic --
     * and that only pays off over a long enough horizon.
     */
    @Test
    public void testThePositionErrorGrowsLinearlyWhereDormandPrinceGrowsAsTheSquare() {
        double[] symplectic = new double[2];
        double[] runge = new double[2];
        long[] symplecticCost = new long[2];
        long[] rungeCost = new long[2];
        int[] horizons = { 10, 100 };
        for (int i = 0; i < horizons.length; ++i) {
            int orbits = horizons[i];
            double end = orbits * PERIOD;
            OdeIntegrator.Result s = new OdeIntegrator(
                    new SymplecticNystrom(SplittingCoefficients.SUZUKI_4, KEPLER, 2))
                            .solve(0.0, kepler(0.6), end, orbits * 200, new double[] { end });
            OdeIntegrator.Result d = new OdeIntegrator(
                    new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, KEPLER_FLAT, 4),
                    new StepController(1.0e-10, 1.0e-10, 400000))
                            .solve(0.0, kepler(0.6), end, new double[] { end });
            symplectic[i] = positionError(s.finalState());
            runge[i] = positionError(d.finalState());
            symplecticCost[i] = s.evaluations;
            rungeCost[i] = d.evaluations;
        }
        for (int i = 0; i < 2; ++i) {
            assertTrue("the two must cost about the same per orbit",
                    rungeCost[i] < 1.5 * symplecticCost[i] && symplecticCost[i] < 1.5 * rungeCost[i]);
        }
        assertEquals("ten times the orbits, ten times the error", 10.0,
                symplectic[1] / symplectic[0], 0.5);
        assertTrue("and for Dormand-Prince about a hundred, not ten: " + runge[1] / runge[0],
                runge[1] / runge[0] > 40.0);
        assertTrue("at ten orbits the adaptive method is far ahead", runge[0] < 0.01 * symplectic[0]);
    }

    @Test
    public void testInterpolationCostsNothingAndTheVelocityIsTheDerivativeOfThePosition() {
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.SUZUKI_4, OSCILLATOR, 1);
        double[] y = { Math.cos(0.3), -Math.sin(0.3) };
        double[] out = new double[2];
        double[] point = new double[2];
        double h = 0.2;
        stepper.step(0.3, y, h, out, null);
        long afterStep = stepper.evaluations();
        assertEquals(5L, afterStep);
        for (int i = 1; i < 64; ++i) {
            stepper.interpolate(i / 64.0, point);
        }
        assertEquals("an interior point must cost nothing", afterStep, stepper.evaluations());

        // the reported velocity is the slope of the reported position
        double theta = 0.375;
        double eps = 1.0e-6;
        double[] up = new double[2];
        double[] down = new double[2];
        stepper.interpolate(theta + eps, up);
        stepper.interpolate(theta - eps, down);
        stepper.interpolate(theta, point);
        assertEquals("the two halves of the interpolant agree", (up[0] - down[0]) / (2.0 * eps * h),
                point[1], 1.0e-8);
    }

    @Test
    public void testTheEndsOfTheInterpolantAreTheEndsOfTheStep() {
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.YOSHIDA_4, KEPLER, 2);
        double[] y = kepler(0.6);
        double[] out = new double[4];
        double[] point = new double[4];
        stepper.step(0.0, y, 0.1, out, null);
        stepper.interpolate(0.0, point);
        for (int i = 0; i < 4; ++i) {
            assertEquals(y[i], point[i], 0.0);
        }
        stepper.interpolate(1.0, point);
        for (int i = 0; i < 4; ++i) {
            assertEquals(out[i], point[i], 0.0);
        }
    }

    @Test
    public void testTheInterpolationErrorFallsAsTheFourthPowerInPositionAndTheThirdInVelocity() {
        double[] q = new double[5];
        double[] v = new double[5];
        for (int i = 0; i < 5; ++i) {
            double h = 0.4 / Math.pow(2.0, i);
            SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.SUZUKI_4, OSCILLATOR, 1);
            double[] out = new double[2];
            double[] point = new double[2];
            stepper.step(0.3, new double[] { Math.cos(0.3), -Math.sin(0.3) }, h, out, null);
            double worstQ = 0.0;
            double worstV = 0.0;
            for (int j = 1; j < 32; ++j) {
                double theta = j / 32.0;
                stepper.interpolate(theta, point);
                double t = 0.3 + theta * h;
                worstQ = Math.max(worstQ, Math.abs(point[0] - Math.cos(t)));
                worstV = Math.max(worstV, Math.abs(point[1] + Math.sin(t)));
            }
            q[i] = worstQ;
            v[i] = worstV;
        }
        for (int i = 1; i < 5; ++i) {
            assertEquals("position, halving " + i, 4.0, Math.log(q[i - 1] / q[i]) / Math.log(2.0), 0.15);
            assertEquals("velocity, halving " + i, 3.0, Math.log(v[i - 1] / v[i]) / Math.log(2.0), 0.15);
        }
    }

    /**
     * The bounded behavior is a property of a constant step size, so a
     * controller would take it away. The refusal was written for another reason
     * -- there is no error estimate to control on -- and turns out to be the
     * right one here for a second reason.
     */
    @Test
    public void testTheDriverRefusesToChooseTheStepSize() {
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.SUZUKI_4, KEPLER, 2);
        assertFalse(stepper.hasErrorEstimate());
        try {
            new OdeIntegrator(stepper, new StepController(1.0e-9, 1.0e-9));
            fail("a symplectic method must not be driven at a varying step size");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("estimates the error"));
        }
    }

    /**
     * Events read the state through the interpolant and nothing else, so they
     * work here unchanged: the orbit passes its periapsis once per period.
     */
    @Test
    public void testEventsWorkOnASymplecticRun() {
        Event apsis = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0] * y[2] + y[1] * y[3];
            }
        });
        OdeIntegrator.Result r = new OdeIntegrator(
                new SymplecticNystrom(SplittingCoefficients.SUZUKI_4, KEPLER, 2), null,
                new Event[] { apsis }).solve(0.0, kepler(0.6), 3.0 * PERIOD, 3 * 400);
        assertTrue("apoapsis and periapsis, twice per orbit", r.eventCount() >= 5);
        for (int i = 0; i < r.eventCount(); ++i) {
            double distance = Math.hypot(r.eventStates[i][0], r.eventStates[i][1]);
            boolean far = Math.abs(distance - 1.6) < 1.0e-4;
            boolean near = Math.abs(distance - 0.4) < 1.0e-4;
            assertTrue("every turning point is at an apsis, not at " + distance, far || near);
        }
    }

    /**
     * The caveat, made measurable. A force that reads the velocity is still
     * integrated, and the one step map stops preserving area -- which is
     * correct, because a damped oscillator does not preserve it either.
     */
    @Test
    public void testAVelocityDependentForceIsIntegratedButNoLongerSymplectically() {
        SymplecticNystrom damped = new SymplecticNystrom(SplittingCoefficients.VERLET, DAMPED, 1);
        double[][] jacobian = jacobianOfOneStep(damped, new double[] { 1.0, 0.0 }, 0.05);
        double determinant = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
        assertTrue("a damped system must lose area, and this one loses " + (1.0 - determinant),
                Math.abs(determinant - 1.0) > 1.0e-4);
        assertTrue("but it must lose it, not gain it", determinant < 1.0);

        SymplecticNystrom undamped = new SymplecticNystrom(SplittingCoefficients.VERLET, OSCILLATOR, 1);
        double[][] clean = jacobianOfOneStep(undamped, new double[] { 1.0, 0.0 }, 0.05);
        assertEquals("where the same method on the same shape of problem does not", 1.0,
                clean[0][0] * clean[1][1] - clean[0][1] * clean[1][0], 1.0e-9);
    }

    @Test
    public void testTheStepperReportsWhatItIs() {
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.YOSHIDA_4, KEPLER, 2);
        assertSame(SplittingCoefficients.YOSHIDA_4, stepper.coefficients());
        assertEquals(2, stepper.degreesOfFreedom());
        assertEquals(4, stepper.dimension());
        assertEquals(4, stepper.order());
        assertFalse(stepper.hasErrorEstimate());
        assertFalse(stepper.hasDenseOutput());
        assertTrue(Double.isNaN(stepper.stiffnessMeasure()));
        assertTrue(stepper.toString().contains("degrees of freedom"));

        double[] dydt = new double[4];
        stepper.derivative(0.0, kepler(0.6), dydt);
        assertEquals("the first half of the derivative is the velocity", kepler(0.6)[2], dydt[0], 0.0);
        assertEquals(kepler(0.6)[3], dydt[1], 0.0);
        assertEquals("and the second half the acceleration", -1.0 / (0.4 * 0.4), dydt[2], 1.0e-12);
        assertEquals(1L, stepper.evaluations());
    }

    @Test
    public void testResetForgetsTheStepButNotTheCost() {
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.VERLET, OSCILLATOR, 1);
        stepper.step(0.0, new double[] { 1.0, 0.0 }, 0.1, new double[2], null);
        assertEquals(1L, stepper.evaluations());
        stepper.reset();
        assertEquals(1L, stepper.evaluations());
        try {
            stepper.interpolate(0.5, new double[2]);
            fail("there is no step to interpolate after a reset");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("no step"));
        }
    }

    @Test
    public void testTheArgumentsAreChecked() {
        try {
            new SymplecticNystrom(null, OSCILLATOR, 1);
            fail("expected a refusal naming the coefficients");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("coefficients"));
        }
        try {
            new SymplecticNystrom(SplittingCoefficients.VERLET, null, 1);
            fail("expected a refusal naming the field");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("field"));
        }
        try {
            new SymplecticNystrom(SplittingCoefficients.VERLET, OSCILLATOR, 0);
            fail("expected a refusal naming the degrees of freedom");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("degreesOfFreedom"));
        }
        SymplecticNystrom stepper = new SymplecticNystrom(SplittingCoefficients.VERLET, OSCILLATOR, 1);
        double[] y = { 1.0, 0.0 };
        double[] out = new double[2];
        refuse(stepper, new double[3], 0.1, out, "y");
        refuse(stepper, y, 0.1, new double[3], "yOut");
        refuse(stepper, y, 0.1, y, "must not be the array passed as y");
        refuse(stepper, y, 0.0, out, "h must not be zero");
        try {
            stepper.derivative(0.0, new double[3], new double[2]);
            fail("expected a refusal naming y");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("y must be of length"));
        }
    }

    // ------------------------------------------------------------------

    private static double[] kepler(double eccentricity) {
        return new double[] { 1.0 - eccentricity, 0.0, 0.0,
                Math.sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) };
    }

    private static double energy(double[] y) {
        return 0.5 * (y[2] * y[2] + y[3] * y[3]) - 1.0 / Math.hypot(y[0], y[1]);
    }

    private static double momentum(double[] y) {
        return y[0] * y[3] - y[1] * y[2];
    }

    /**
     * The evaluations the method needs to bring a twenty orbit run of the two
     * body problem back to within {@code accuracy} of where it started, doubling
     * the step count until it gets there.
     */
    private static long costToReach(SplittingCoefficients c, double end, double[] exact, double accuracy) {
        for (int steps = 20 * 25; steps <= 20 * 6400; steps *= 2) {
            OdeIntegrator.Result r = new OdeIntegrator(new SymplecticNystrom(c, KEPLER, 2)).solve(0.0,
                    kepler(0.6), end, steps, new double[] { end });
            double worst = 0.0;
            for (int i = 0; i < 2; ++i) {
                worst = Math.max(worst, Math.abs(r.finalState()[i] - exact[i]));
            }
            if (worst <= accuracy) {
                return r.evaluations;
            }
        }
        throw new AssertionError(c.name() + " never reached " + accuracy);
    }

    private static double positionError(double[] y) {
        double[] exact = kepler(0.6);
        return Math.max(Math.abs(y[0] - exact[0]), Math.abs(y[1] - exact[1]));
    }

    private static double energyBand(OdeIntegrator.Result r, int from, int to) {
        double low = Double.MAX_VALUE;
        double high = -Double.MAX_VALUE;
        for (int i = from; i < to; ++i) {
            double e = energy(r.y[i]);
            low = Math.min(low, e);
            high = Math.max(high, e);
        }
        return high - low;
    }

    private static double[][] jacobianOfOneStep(OdeStepper stepper, double[] base, double h) {
        int n = base.length;
        double[][] jacobian = new double[n][n];
        double eps = 1.0e-6;
        double[] up = new double[n];
        double[] down = new double[n];
        for (int j = 0; j < n; ++j) {
            double[] plus = base.clone();
            double[] minus = base.clone();
            plus[j] += eps;
            minus[j] -= eps;
            stepper.step(0.0, plus, h, up, null);
            stepper.step(0.0, minus, h, down, null);
            for (int i = 0; i < n; ++i) {
                jacobian[i][j] = (up[i] - down[i]) / (2.0 * eps);
            }
        }
        return jacobian;
    }

    /** The worst entry of M' J M - J, which is zero for a symplectic map. */
    private static double symplecticDefect(double[][] m) {
        int n = m.length;
        int half = n / 2;
        double[][] j = new double[n][n];
        for (int i = 0; i < half; ++i) {
            j[i][half + i] = 1.0;
            j[half + i][i] = -1.0;
        }
        double worst = 0.0;
        for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
                double sum = 0.0;
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        sum += m[k][p] * j[k][l] * m[l][q];
                    }
                }
                worst = Math.max(worst, Math.abs(sum - j[p][q]));
            }
        }
        return worst;
    }

    private static void refuse(SymplecticNystrom stepper, double[] y, double h, double[] out, String hint) {
        try {
            stepper.step(0.0, y, h, out, null);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }
}
