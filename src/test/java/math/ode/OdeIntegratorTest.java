package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DFunction;
import math.fun.DVectorField;
import math.solve.Quadrature;

/**
 * {@link OdeIntegrator} against solutions that are known in closed form, and
 * against one algorithm that shares no line of code with it.
 * <p>
 * <b>At a fixed step size</b> the claim is the order of the method over a whole
 * interval: halving the step divides the error at the end by two to the order.
 * Measured on the harmonic oscillator over {@code [0, 2]} with the step count
 * doubling from ten, the classical method gives {@code 4.00} five times over
 * and Dormand-Prince {@code 5.01, 5.00, 5.00, 4.99}, the last of which is down
 * at {@code 1.7e-13} and running out of double precision rather than out of
 * order.
 * <p>
 * <b>Adaptively</b> the claim is that the tolerance means something. Over eight
 * decades of it, the error at the end of a run to {@code t = 20} stays between
 * {@code 6.1} and {@code 7.7} times {@code rtol} -- a constant of the equation
 * and not of the tolerance, which is what makes the setting worth having -- and
 * the number of steps grows as its fifth root, {@code 2.42, 2.48, 2.51, 2.51}
 * per factor of a hundred against the {@code 2.512} the order predicts.
 * <p>
 * <b>The free oracle</b> is that an equation with no {@code y} on the right
 * hand side is a quadrature, so {@code math.solve} answers the same question by
 * a route with nothing in common: at {@code rtol = 1e-12} the two agree to
 * {@code 4.1e-14}.
 * <p>
 * <b>And Kepler says what a method of this kind cannot do.</b> The energy of a
 * two body orbit is conserved exactly and this method loses it linearly with
 * time -- {@code 1.28e-09}, {@code 1.18e-08}, {@code 1.16e-07} over ten, a
 * hundred and a thousand orbits, a factor of ten per factor of ten. That is not
 * a defect of the implementation; it is what a symplectic method is for, and
 * this is the measurement it will be held against.
 */
public final class OdeIntegratorTest {

    /** y'' = -y as a first order pair, exact solution (cos t, -sin t) at (1, 0). */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    /** y' = y (1 - y), the logistic equation, nonlinear and scalar. */
    private static final DVectorField LOGISTIC = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[0] * (1.0 - y[0]);
        }
    };

    /** y' = y squared, which reaches infinity at t = 1 from y(0) = 1. */
    private static final DVectorField BLOW_UP = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[0] * y[0];
        }
    };

    /** A right hand side that is not a number anywhere. */
    private static final DVectorField NOWHERE = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = Double.NaN;
        }
    };

    /** An integrand rather than an equation: y' = f(t) with no y in it. */
    private static final DFunction INTEGRAND = new DFunction() {
        @Override
        public double apply(double t) {
            return Math.exp(-t * t) * Math.cos(3.0 * t) + 0.5;
        }
    };

    private static final DVectorField QUADRATURE = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = INTEGRAND.apply(t);
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

    private static final double LV_A = 1.5;
    private static final double LV_B = 1.0;
    private static final double LV_C = 3.0;
    private static final double LV_D = 1.0;

    /** Predator and prey, which has a conserved quantity but no closed form. */
    private static final DVectorField LOTKA = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = LV_A * y[0] - LV_B * y[0] * y[1];
            dydt[1] = -LV_C * y[1] + LV_D * y[0] * y[1];
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

    // ------------------------------------------------------------------
    // a fixed step size
    // ------------------------------------------------------------------

    @Test
    public void testTheClassicalMethodConvergesAtFourthOrder() {
        assertSlopes(globalErrors(ButcherTableau.CLASSIC_RK4, 10, 6), 4.0, 0.05);
    }

    @Test
    public void testDormandPrinceConvergesAtFifthOrder() {
        assertSlopes(globalErrors(ButcherTableau.DORMAND_PRINCE_45, 10, 5), 5.0, 0.05);
    }

    /**
     * The output grid and the steps are two names for the same states when they
     * fall on the same times, and that has to hold exactly: a value asked for
     * at a step boundary is that step, not an interpolation that happens to
     * land a few units in the last place away from it.
     */
    @Test
    public void testAGridOnTheStepsReproducesThemBitForBit() {
        ButcherTableau[] both = { ButcherTableau.CLASSIC_RK4, ButcherTableau.DORMAND_PRINCE_45 };
        for (int q = 0; q < both.length; ++q) {
            OdeIntegrator.Result plain = fixed(both[q], OSCILLATOR, 2).solve(0.0, start(), 4.0, 40);
            OdeIntegrator.Result grid = fixed(both[q], OSCILLATOR, 2)
                    .solve(0.0, start(), 4.0, 40, plain.t.clone());

            assertEquals(plain.length, grid.length);
            for (int i = 0; i < plain.length; ++i) {
                assertEquals(both[q].name(), plain.t[i], grid.t[i], 0.0);
                assertEquals(both[q].name(), plain.y[i][0], grid.y[i][0], 0.0);
                assertEquals(both[q].name(), plain.y[i][1], grid.y[i][1], 0.0);
            }
            assertEquals("a grid on the steps must not cost an evaluation", plain.evaluations,
                    grid.evaluations);
        }
    }

    @Test
    public void testAGridOffTheStepsIsNoWorseThanTheStepsThemselves() {
        ButcherTableau[] both = { ButcherTableau.CLASSIC_RK4, ButcherTableau.DORMAND_PRINCE_45 };
        for (int q = 0; q < both.length; ++q) {
            double[] times = evenlySpaced(0.0, 4.0, 97);
            OdeIntegrator.Result grid = fixed(both[q], OSCILLATOR, 2).solve(0.0, start(), 4.0, 40, times);
            OdeIntegrator.Result plain = fixed(both[q], OSCILLATOR, 2).solve(0.0, start(), 4.0, 40);
            assertTrue(both[q].name(), worstAgainstTheCosine(grid) <= worstAgainstTheCosine(plain) * 1.05);
            for (int i = 0; i < times.length; ++i) {
                assertEquals("the times come back as they were given", times[i], grid.t[i], 0.0);
            }
        }
    }

    /**
     * What the continuous extension is worth at driver level: interpolating
     * onto a grid finer than the steps costs Dormand-Prince nothing at all,
     * while the classical method pays one evaluation per step for the
     * derivative its cubic fallback needs.
     */
    @Test
    public void testWhatInterpolatingOntoAFinerGridCosts() {
        double[] times = evenlySpaced(0.0, 4.0, 97);
        OdeIntegrator.Result dp = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2)
                .solve(0.0, start(), 4.0, 40, times);
        OdeIntegrator.Result dpPlain = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2)
                .solve(0.0, start(), 4.0, 40);
        assertEquals(dpPlain.evaluations, dp.evaluations);

        OdeIntegrator.Result rk = fixed(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2)
                .solve(0.0, start(), 4.0, 40, times);
        OdeIntegrator.Result rkPlain = fixed(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2)
                .solve(0.0, start(), 4.0, 40);
        assertEquals(160L, rkPlain.evaluations);
        assertEquals(200L, rk.evaluations);
    }

    @Test
    public void testIntegratingBackwardsUndoesIntegratingForwards() {
        ButcherTableau[] both = { ButcherTableau.CLASSIC_RK4, ButcherTableau.DORMAND_PRINCE_45 };
        for (int q = 0; q < both.length; ++q) {
            OdeIntegrator.Result there = fixed(both[q], OSCILLATOR, 2).solve(0.0, start(), 3.0, 200);
            OdeIntegrator.Result back = fixed(both[q], OSCILLATOR, 2)
                    .solve(3.0, there.finalState(), 0.0, 200);
            assertEquals(0.0, back.finalTime(), 1.0e-13);
            assertEquals(both[q].name(), 1.0, back.finalState()[0], 1.0e-9);
            assertEquals(both[q].name(), 0.0, back.finalState()[1], 1.0e-9);
        }
    }

    @Test
    public void testANonlinearScalarEquation() {
        double y0 = 0.01;
        OdeIntegrator.Result r = fixed(ButcherTableau.DORMAND_PRINCE_45, LOGISTIC, 1)
                .solve(0.0, new double[] { y0 }, 12.0, 200);
        double exact = 1.0 / (1.0 + (1.0 / y0 - 1.0) * Math.exp(-r.finalTime()));
        assertEquals(exact, r.finalState()[0], 1.0e-11);
    }

    /**
     * A fixed step size cannot react to a solution that runs off to infinity,
     * so the honest thing is to say where it happened rather than hand back a
     * state full of {@code Infinity}.
     */
    @Test
    public void testAStateThatStopsBeingFiniteIsReported() {
        try {
            fixed(ButcherTableau.CLASSIC_RK4, BLOW_UP, 1).solve(0.0, new double[] { 1.0 }, 10.0, 10);
            fail("a solution reaching infinity at t = 1 must not be integrated to t = 10");
        } catch (ArithmeticException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("not finite"));
            assertTrue(expected.getMessage(), expected.getMessage().contains("step "));
        }
    }

    @Test
    public void testTheResultReportsWhatItDid() {
        OdeIntegrator.Result r = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2)
                .solve(0.0, start(), 2.0, 20);
        assertEquals(21, r.length);
        assertEquals(21, r.t.length);
        assertEquals(21, r.y.length);
        assertEquals(0.0, r.t[0], 0.0);
        assertEquals(2.0, r.finalTime(), 1.0e-15);
        assertEquals(20L, r.steps);
        assertEquals("a fixed step size throws nothing away", 0L, r.rejected);
        assertEquals(0.0, r.rejectionRate(), 0.0);
        // seven for the first step, six for each of the nineteen that follow
        assertEquals(121L, r.evaluations);
        assertSame(r.y[r.length - 1], r.finalState());
        assertTrue(r.toString().contains("20 steps"));
    }

    @Test
    public void testTheInitialStateIsNotWrittenInto() {
        double[] y0 = start();
        fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2).solve(0.0, y0, 2.0, 20);
        assertEquals(1.0, y0[0], 0.0);
        assertEquals(0.0, y0[1], 0.0);
    }

    @Test
    public void testTheArgumentsAreChecked() {
        OdeIntegrator in = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        refuse(in, 0.0, new double[3], 1.0, 10, null, "y0");
        refuse(in, 0.0, start(), 1.0, 0, null, "steps");
        refuse(in, 0.0, start(), 0.0, 10, null, "must differ");
        refuse(in, 0.0, start(), Double.POSITIVE_INFINITY, 10, null, "finite");
        refuse(in, 0.0, new double[] { Double.NaN, 0.0 }, 1.0, 10, null, "y0 must be finite");
        refuse(in, 0.0, start(), 1.0, 10, new double[] { 0.5, 1.5 }, "outside the interval");
        refuse(in, 0.0, start(), 1.0, 10, new double[] { 0.5, 0.25 }, "turns back");
        refuse(in, 0.0, start(), 1.0, 10, new double[] { Double.NaN }, "not finite");
        try {
            in.solve(0.0, start(), 1.0, 10, null);
            fail("expected a refusal naming the times");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("times"));
        }
        try {
            new OdeIntegrator(null);
            fail("expected a refusal naming the stepper");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("stepper"));
        }
        assertNull(in.controller());
        assertEquals(2, in.stepper().dimension());
    }

    @Test
    public void testABackwardGridIsCheckedInItsOwnDirection() {
        OdeIntegrator in = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        OdeIntegrator.Result r = in.solve(1.0, start(), 0.0, 10, new double[] { 1.0, 0.5, 0.0 });
        assertEquals(3, r.length);
        assertEquals(1.0, r.t[0], 0.0);
        assertEquals(0.0, r.finalTime(), 0.0);
        refuse(in, 1.0, start(), 0.0, 10, new double[] { 0.5, 0.75 }, "turns back");
    }

    // ------------------------------------------------------------------
    // choosing the step size
    // ------------------------------------------------------------------

    /**
     * What a tolerance is worth: the error at the end is a fixed multiple of it
     * over eight decades, so halving the setting halves the answer's error.
     * That the multiple is about six rather than one is the difference between
     * bounding what one step may add and bounding what the run accumulates.
     */
    @Test
    public void testTheErrorAtTheEndIsProportionalToTheTolerance() {
        double smallest = Double.MAX_VALUE;
        double largest = 0.0;
        for (int i = 0; i < 5; ++i) {
            double rtol = Math.pow(10.0, -(4 + 2 * i));
            OdeIntegrator.Result r = adaptive(OSCILLATOR, 2, rtol, 0.04).solve(0.0, start(), 20.0);
            double ratio = worstAgainstTheCosine(r) / rtol;
            assertTrue("rtol " + rtol + " gave " + ratio + " times it", ratio >= 1.0 && ratio <= 25.0);
            smallest = Math.min(smallest, ratio);
            largest = Math.max(largest, ratio);
        }
        assertTrue("the multiple is a constant of the equation, not of the tolerance: " + smallest + " to "
                + largest, largest / smallest <= 2.0);
    }

    /**
     * The cost of a tolerance follows from the order: to divide the error by
     * {@code 2^5} the step has to halve, so the step count grows as the fifth
     * root of the tolerance and not as its reciprocal.
     */
    @Test
    public void testTheCostGrowsAsTheFifthRootOfTheTolerance() {
        long[] steps = new long[5];
        for (int i = 0; i < steps.length; ++i) {
            double rtol = Math.pow(10.0, -(4 + 2 * i));
            steps[i] = adaptive(OSCILLATOR, 2, rtol, 0.04).solve(0.0, start(), 20.0).steps;
        }
        double predicted = Math.pow(100.0, 1.0 / 5.0);
        for (int i = 1; i < steps.length; ++i) {
            double ratio = steps[i] / (double) steps[i - 1];
            assertEquals("a hundredfold tighter tolerance", predicted, ratio, 0.25);
        }
    }

    /**
     * An equation with no {@code y} on its right hand side is an integral, and
     * {@link Quadrature} answers it by adaptive Gauss-Kronrod: a different
     * algorithm, a different error estimate, no shared line of code.
     */
    @Test
    public void testAgainstGaussKronrodOnAnEquationThatIsReallyAnIntegral() {
        double byQuadrature = Quadrature.integrate(INTEGRAND, 0.0, 3.0, 1.0e-13);
        OdeIntegrator.Result loose = adaptive(QUADRATURE, 1, 1.0e-6, 0.04)
                .solve(0.0, new double[] { 0.0 }, 3.0);
        OdeIntegrator.Result tight = adaptive(QUADRATURE, 1, 1.0e-12, 0.04)
                .solve(0.0, new double[] { 0.0 }, 3.0);
        assertEquals(byQuadrature, loose.finalState()[0], 1.0e-5);
        assertEquals(byQuadrature, tight.finalState()[0], 1.0e-11);
    }

    /**
     * The energy and the angular momentum of a two body orbit are conserved
     * exactly, and an explicit Runge-Kutta method loses them steadily: the
     * drift after a thousand orbits is a hundred times the drift after ten.
     * A symplectic method is chosen precisely to break that proportionality,
     * and this is the measurement it has to beat.
     */
    @Test
    public void testKeplerLosesEnergyInProportionToTheTimeItRunsFor() {
        double period = 2.0 * Math.PI;
        double[] drift = new double[3];
        int[] orbits = { 10, 100, 1000 };
        for (int i = 0; i < orbits.length; ++i) {
            OdeIntegrator.Result r = adaptive(KEPLER, 4, 1.0e-10, 0.04)
                    .solve(0.0, kepler(0.6), orbits[i] * period);
            double energy0 = keplerEnergy(kepler(0.6));
            double momentum0 = keplerMomentum(kepler(0.6));
            double worstEnergy = 0.0;
            double worstMomentum = 0.0;
            for (int j = 0; j < r.length; ++j) {
                worstEnergy = Math.max(worstEnergy, Math.abs(keplerEnergy(r.y[j]) - energy0));
                worstMomentum = Math.max(worstMomentum, Math.abs(keplerMomentum(r.y[j]) - momentum0));
            }
            assertTrue("the angular momentum drifts no faster than the energy",
                    worstMomentum <= worstEnergy);
            drift[i] = worstEnergy;
        }
        assertTrue("ten orbits should still be accurate: " + drift[0], drift[0] < 1.0e-8);
        for (int i = 1; i < drift.length; ++i) {
            double ratio = drift[i] / drift[i - 1];
            assertTrue("ten times as long should drift about ten times as far, not " + ratio,
                    ratio > 7.0 && ratio < 13.0);
        }
    }

    /**
     * Lotka-Volterra has no closed form solution and a quantity that does not
     * change along it, which is an error measure with nothing to compare
     * against.
     */
    @Test
    public void testLotkaVolterraKeepsTheQuantityThatDoesNotChange() {
        double[] y0 = { 10.0, 5.0 };
        double before = lotkaInvariant(y0);
        double looser = 0.0;
        for (int p = 6; p <= 12; p += 3) {
            double tol = Math.pow(10.0, -p);
            OdeIntegrator.Result r = adaptive(LOTKA, 2, tol, 0.04).solve(0.0, y0.clone(), 30.0);
            double worst = 0.0;
            for (int i = 0; i < r.length; ++i) {
                worst = Math.max(worst, Math.abs(lotkaInvariant(r.y[i]) - before));
                assertTrue("both populations stay positive", r.y[i][0] > 0.0 && r.y[i][1] > 0.0);
            }
            assertTrue("the drift must fall with the tolerance", looser == 0.0 || worst < looser);
            assertTrue("at tolerance " + tol + " the invariant moved by " + worst, worst < 1000.0 * tol);
            looser = worst;
        }
    }

    @Test
    public void testTheLastStepLandsExactlyOnTheTargetInEitherDirection() {
        assertEquals(20.0, adaptive(OSCILLATOR, 2, 1.0e-9, 0.04).solve(0.0, start(), 20.0).finalTime(), 0.0);
        assertEquals(-7.5, adaptive(OSCILLATOR, 2, 1.0e-9, 0.04).solve(0.0, start(), -7.5).finalTime(), 0.0);
    }

    @Test
    public void testIntegratingBackwardsAdaptivelyUndoesIntegratingForwards() {
        OdeIntegrator.Result there = adaptive(KEPLER, 4, 1.0e-11, 0.04).solve(0.0, kepler(0.6), 12.0);
        OdeIntegrator.Result back = adaptive(KEPLER, 4, 1.0e-11, 0.04).solve(12.0, there.finalState(), 0.0);
        double[] origin = kepler(0.6);
        for (int i = 0; i < origin.length; ++i) {
            assertEquals("component " + i, origin[i], back.finalState()[i], 1.0e-7);
        }
    }

    /**
     * Where an adaptive run puts its steps is not where a caller would have
     * put them, so asking for a grid is the normal thing to do -- and it has to
     * be free, or the continuous extension was not worth carrying.
     */
    @Test
    public void testAnAdaptiveGridIsFreeAndNoWorseThanTheStepsItInterpolates() {
        double[] times = evenlySpaced(0.0, 20.0, 41);
        OdeIntegrator.Result grid = adaptive(OSCILLATOR, 2, 1.0e-10, 0.04).solve(0.0, start(), 20.0, times);
        OdeIntegrator.Result steps = adaptive(OSCILLATOR, 2, 1.0e-10, 0.04).solve(0.0, start(), 20.0);
        assertEquals(41, grid.length);
        assertEquals("the grid must not cost an evaluation", steps.evaluations, grid.evaluations);
        assertEquals("nor change what was computed", steps.steps, grid.steps);
        assertTrue(worstAgainstTheCosine(grid) <= worstAgainstTheCosine(steps) * 1.05);
        for (int i = 0; i < times.length; ++i) {
            assertEquals(times[i], grid.t[i], 0.0);
        }
    }

    /**
     * The finding the default {@code beta} rests on: damping costs a few
     * percent where the step size moves slowly and saves a factor of ten in
     * rejected steps where it has to move fast.
     */
    @Test
    public void testDampingIsWhatSavesAStiffProblemFromItsOwnRejections() {
        OdeIntegrator.Result textbook = adaptive(vanDerPol(100.0), 2, 1.0e-6, 0.0)
                .solve(0.0, new double[] { 2.0, 0.0 }, 300.0);
        OdeIntegrator.Result damped = adaptive(vanDerPol(100.0), 2, 1.0e-6, 0.04)
                .solve(0.0, new double[] { 2.0, 0.0 }, 300.0);
        assertTrue("the two must agree on where they got to",
                Math.abs(textbook.finalState()[0] - damped.finalState()[0]) < 1.0e-3);
        assertTrue("damping should cut the rejections by a large factor, not " + textbook.rejected + " to "
                + damped.rejected, damped.rejected * 5L < textbook.rejected);
        assertTrue("and should not cost evaluations here", damped.evaluations <= textbook.evaluations);
    }

    @Test
    public void testAdaptiveSolvingWithoutAControllerIsRefused() {
        OdeIntegrator in = fixed(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        try {
            in.solve(0.0, start(), 1.0);
            fail("expected a refusal naming the controller");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("step controller"));
        }
    }

    @Test
    public void testAControllerNeedsAMethodThatEstimatesItsOwnError() {
        try {
            new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2),
                    new StepController(1.0e-8, 1.0e-8));
            fail("the classical method has no embedded pair and cannot be controlled");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("estimates the error"));
        }
    }

    /**
     * Approaching a singularity, the step size falls until it is the resolution
     * of the time itself. Saying so, with the time it happened at, is more use
     * than a state full of infinities: the run below stops within
     * {@code 1.1e-09} of where the solution actually ends.
     */
    @Test
    public void testAStepSizeThatCollapsesIsReportedWithWhereItHappened() {
        try {
            adaptive(BLOW_UP, 1, 1.0e-8, 0.04).solve(0.0, new double[] { 1.0 }, 2.0);
            fail("a solution that ends at t = 1 must not be integrated to t = 2");
        } catch (ArithmeticException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("collapsed"));
            assertTrue(expected.getMessage(), expected.getMessage().contains("t = 1.0"));
        }
    }

    @Test
    public void testARightHandSideThatIsNotFiniteAtTheStartIsReported() {
        try {
            adaptive(NOWHERE, 1, 1.0e-8, 0.04).solve(0.0, new double[] { 1.0 }, 1.0);
            fail("there is no first step to be found here");
        } catch (ArithmeticException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("not finite"));
        }
    }

    /**
     * The step budget is the backstop, and what it usually means is that the
     * equation is stiff -- which the message says, because an explicit method
     * on a stiff problem does not get slow, it stops.
     */
    @Test
    public void testTheStepBudgetIsEnforcedAndSaysWhatItSuspects() {
        StepController tiny = new StepController(1.0e-6, 1.0e-6, 2000);
        OdeIntegrator in = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, vanDerPol(1000.0), 2), tiny);
        try {
            in.solve(0.0, new double[] { 2.0, 0.0 }, 3000.0);
            fail("van der Pol at mu = 1000 does not fit in two thousand steps");
        } catch (ArithmeticException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("budget of 2000"));
            assertTrue(expected.getMessage(), expected.getMessage().contains("stiff"));
        }
    }

    @Test
    public void testAnUpperBoundOnTheStepSizeIsKept() {
        StepController bounded = new StepController(1.0e-6, 1.0e-6, 0.9, 0.2, 10.0, 0.04, 0.05, 100000);
        OdeIntegrator in = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2), bounded);
        OdeIntegrator.Result r = in.solve(0.0, start(), 10.0);
        for (int i = 1; i < r.length; ++i) {
            assertTrue("step " + i + " ran from " + r.t[i - 1] + " to " + r.t[i],
                    r.t[i] - r.t[i - 1] <= 0.05 + 1.0e-12);
        }
        assertTrue("at most 0.05 per step over an interval of 10", r.steps >= 200L);
        assertSame(bounded, in.controller());
    }

    @Test
    public void testTheAdaptiveArgumentsAreChecked() {
        OdeIntegrator in = adaptive(OSCILLATOR, 2, 1.0e-8, 0.04);
        try {
            in.solve(0.0, start(), 1.0, (double[]) null);
            fail("expected a refusal naming the times");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("times"));
        }
        try {
            in.solve(0.0, start(), 1.0, new double[] { 0.5, 1.5 });
            fail("expected a refusal about the interval");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("outside the interval"));
        }
        try {
            in.solve(0.0, new double[3], 1.0);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("y0"));
        }
    }

    // ------------------------------------------------------------------
    // what the method says about its own limits
    // ------------------------------------------------------------------

    /**
     * The non-stiffness test, on a family where the stiffness is a dial. The
     * equation is {@code y' = -lambda (y - cos t) - sin t}, whose solution is
     * {@code cos t} for every {@code lambda} and whose Jacobian is
     * {@code -lambda}: the answer does not change and the difficulty of getting
     * it does. Measured, the estimate comes out at {@code 0.37}, {@code 2.77}
     * and {@code 4.27} for {@code lambda} of one, a hundred and ten thousand,
     * against a threshold of {@code 3.25}.
     */
    @Test
    public void testTheStiffnessEstimateFollowsTheEigenvalueAndNotTheAnswer() {
        double[] measures = new double[3];
        boolean[] flagged = new boolean[3];
        double[] lambdas = { 1.0, 100.0, 10000.0 };
        for (int i = 0; i < lambdas.length; ++i) {
            OdeIntegrator.Result r = adaptive(nearlyCosine(lambdas[i]), 1, 1.0e-6, 0.04)
                    .solve(0.0, new double[] { 1.0 }, 10.0);
            assertEquals("the solution is a cosine whatever lambda is", Math.cos(r.finalTime()),
                    r.finalState()[0], 1.0e-5);
            measures[i] = r.stiffness;
            flagged[i] = r.seemsStiff;
        }
        for (int i = 1; i < measures.length; ++i) {
            assertTrue("a larger eigenvalue must look stiffer: " + measures[i - 1] + " then " + measures[i],
                    measures[i] > measures[i - 1]);
        }
        assertFalse("lambda of one is not stiff", flagged[0]);
        assertTrue("lambda of ten thousand is", flagged[2]);
        assertTrue("and it is past the threshold", measures[2] > 3.25);
    }

    @Test
    public void testASmoothEquationDoesNotLookStiff() {
        OdeIntegrator.Result r = adaptive(OSCILLATOR, 2, 1.0e-6, 0.04).solve(0.0, start(), 50.0);
        assertFalse(r.seemsStiff);
        assertTrue("well inside the stability region: " + r.stiffness, r.stiffness < 1.0);
        assertFalse(r.toString().contains("stiff"));
    }

    @Test
    public void testVanDerPolSaysWhatItIs() {
        OdeIntegrator.Result r = adaptive(vanDerPol(100.0), 2, 1.0e-6, 0.04)
                .solve(0.0, new double[] { 2.0, 0.0 }, 300.0);
        assertTrue("van der Pol at mu = 100 is stiff", r.seemsStiff);
        assertTrue("and sits against the stability limit: " + r.stiffness, r.stiffness > 3.25);
        assertTrue(r.toString().contains("stiff"));
    }

    /**
     * The estimate is a difference quotient between two stages that sit at the
     * same time, and the classical method has only one such stage. Saying
     * nothing is better than saying something made up.
     */
    @Test
    public void testAMethodWithNoTwoStagesAtTheEndOfTheStepHasNoEstimate() {
        OdeIntegrator.Result r = fixed(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2)
                .solve(0.0, start(), 10.0, 100);
        assertTrue(Double.isNaN(r.stiffness));
        assertFalse(r.seemsStiff);
        assertTrue(Double.isNaN(fixed(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2).stepper()
                .stiffnessMeasure()));
    }

    // ------------------------------------------------------------------

    /** y' = -lambda (y - cos t) - sin t, whose solution is cos t for every lambda. */
    private static DVectorField nearlyCosine(final double lambda) {
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = -lambda * (y[0] - Math.cos(t)) - Math.sin(t);
            }
        };
    }

    private static double[] start() {
        return new double[] { 1.0, 0.0 };
    }

    private static double[] kepler(double eccentricity) {
        return new double[] { 1.0 - eccentricity, 0.0, 0.0,
                Math.sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) };
    }

    private static double keplerEnergy(double[] y) {
        return 0.5 * (y[2] * y[2] + y[3] * y[3]) - 1.0 / Math.hypot(y[0], y[1]);
    }

    private static double keplerMomentum(double[] y) {
        return y[0] * y[3] - y[1] * y[2];
    }

    private static double lotkaInvariant(double[] y) {
        return LV_D * y[0] - LV_C * Math.log(y[0]) + LV_B * y[1] - LV_A * Math.log(y[1]);
    }

    private static double[] evenlySpaced(double from, double to, int count) {
        double[] out = new double[count];
        for (int i = 0; i < count; ++i) {
            out[i] = from + (to - from) * i / (count - 1.0);
        }
        return out;
    }

    private static OdeIntegrator fixed(ButcherTableau tableau, DVectorField field, int dimension) {
        return new OdeIntegrator(new ExplicitRungeKutta(tableau, field, dimension));
    }

    private static OdeIntegrator adaptive(DVectorField field, int dimension, double tol, double beta) {
        StepController control = new StepController(tol, tol, 0.9, 0.2, 10.0, beta,
                Double.POSITIVE_INFINITY, 400000);
        return new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, field, dimension),
                control);
    }

    private static void refuse(OdeIntegrator in, double t0, double[] y0, double t1, int steps,
            double[] times, String hint) {
        try {
            if (times == null) {
                in.solve(t0, y0, t1, steps);
            } else {
                in.solve(t0, y0, t1, steps, times);
            }
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }

    private static double worstAgainstTheCosine(OdeIntegrator.Result r) {
        double worst = 0.0;
        for (int i = 0; i < r.length; ++i) {
            worst = Math.max(worst,
                    Math.hypot(r.y[i][0] - Math.cos(r.t[i]), r.y[i][1] + Math.sin(r.t[i])));
        }
        return worst;
    }

    private static double[] globalErrors(ButcherTableau tableau, int firstSteps, int count) {
        double[] out = new double[count];
        for (int i = 0; i < count; ++i) {
            OdeIntegrator.Result r = fixed(tableau, OSCILLATOR, 2).solve(0.0, start(), 2.0, firstSteps << i);
            double t = r.finalTime();
            out[i] = Math.hypot(r.finalState()[0] - Math.cos(t), r.finalState()[1] + Math.sin(t));
        }
        return out;
    }

    private static void assertSlopes(double[] errors, double expected, double tolerance) {
        for (int i = 1; i < errors.length; ++i) {
            double slope = Math.log(errors[i - 1] / errors[i]) / Math.log(2.0);
            assertEquals("halving " + i + " of " + (errors.length - 1), expected, slope, tolerance);
        }
    }
}
