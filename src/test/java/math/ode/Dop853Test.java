package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.fun.DVectorField;

/**
 * What the eighth order method is for, measured against the fifth order one.
 * <p>
 * {@link ButcherTableauTest} settles whether the coefficients are the ones
 * intended -- all 200 order conditions through order eight, from a generator
 * rather than a table -- and {@link ExplicitRungeKuttaTest} settles what a step
 * costs. What is left, and what this class is about, is the question a caller
 * actually has: <b>when is an eighth order method worth its twelve evaluations
 * a step?</b>
 * <p>
 * The answer is a tolerance. Dormand-Prince needs steps proportional to
 * <code>rtol^(-1/5)</code> and DOP853 to <code>rtol^(-1/9)</code>, so there is a
 * crossing, and measured on the two body problem over ten orbits it sits at
 * about {@code 1e-07}: at {@code 1e-06} Dormand-Prince costs 831 evaluations
 * against 891, at {@code 1e-07} it is 1149 against 1095, and by {@code 1e-13}
 * it is 15669 against 4035.
 */
public final class Dop853Test {

    /** The two body problem at eccentricity 0.6, in the usual four components. */
    private static final DVectorField ORBIT = new DVectorField() {
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

    /** <code>y'' = -y</code>, whose solution is a cosine and is known exactly. */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    /** <code>y' = -y^2</code>, whose solution is <code>1 / (1 + t)</code>. */
    private static final DVectorField SQUARE = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = -y[0] * y[0];
        }
    };

    private static double[] orbitStart() {
        double e = 0.6;
        return new double[] { 1.0 - e, 0.0, 0.0, Math.sqrt((1.0 + e) / (1.0 - e)) };
    }

    // ------------------------------------------------------------------
    // the order, and why the error estimate needs two halves
    // ------------------------------------------------------------------

    /**
     * The order conditions say the method is of order eight; this says the
     * arithmetic agrees. Over <code>y'' = -y</code> to {@code t = 40} from
     * twenty fixed steps, halving four times, the error falls by
     * {@code 2^8.24}, {@code 2^8.18}, {@code 2^8.10} and {@code 2^8.05} -- and
     * one halving further it stops falling, because {@code 3.0e-15} is where
     * double precision runs out. That narrow window is what an eighth order
     * method looks like: it reaches the accuracy of the arithmetic in five
     * halvings, where the fifth order method is still at {@code 1.3e-09}.
     */
    @Test
    public void testTheFixedStepErrorFallsLikeTheEighthPower() {
        double t1 = 40.0;
        double[] exact = { Math.cos(t1), -Math.sin(t1) };
        double[] errors = new double[5];
        int steps = 20;
        for (int i = 0; i < errors.length; ++i) {
            double[] end = new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DOP853, OSCILLATOR, 2))
                    .solve(0.0, new double[] { 1.0, 0.0 }, t1, steps).finalState();
            errors[i] = Math.max(Math.abs(end[0] - exact[0]), Math.abs(end[1] - exact[1]));
            steps *= 2;
        }
        for (int i = 0; i + 1 < errors.length; ++i) {
            double slope = Math.log(errors[i] / errors[i + 1]) / Math.log(2.0);
            assertEquals("halving " + i, 8.0, slope, 0.3);
        }
    }

    /**
     * <b>The reason the driver carries a second error estimate at all.</b>
     * <p>
     * DOP853's embedded solutions are of order five and of order three. Taken
     * alone, the fifth order one makes an estimate that falls like
     * <code>h^6</code>, and a step controlled on it is a step held to the
     * accuracy of a solution that is computed and then discarded. Combined as
     * {@link StepController#errorNorm(double[], double[], double[], double[])}
     * combines them, the number falls like <code>h^8</code>, which is the local
     * error of the solution the method actually advances.
     * <p>
     * Measured on <code>y' = -y^2</code> over six halvings from {@code h = 0.8}:
     * the combined estimate at {@code 7.56, 7.80, 7.90, 7.94, 7.97, 7.99} and
     * the first one alone at {@code 5.40, 5.68, 5.83, 5.91, 5.95, 5.98}. Both
     * are still climbing to their limit, so the last of each is the one to
     * read.
     */
    @Test
    public void testTheCombinedEstimateFallsTwoOrdersFasterThanTheFirstOneAlone() {
        StepController c = new StepController(1.0e-8, 1.0e-8);
        double t0 = 0.3;
        int count = 7;
        double[] combined = new double[count];
        double[] alone = new double[count];
        double h = 0.8;
        for (int i = 0; i < count; ++i) {
            ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DOP853, SQUARE, 1);
            double[] y = { 1.0 / (1.0 + t0) };
            double[] yOut = new double[1];
            double[] err = new double[1];
            rk.step(t0, y, h, yOut, err);
            combined[i] = c.errorNorm(err, rk.secondaryError(), y, yOut);
            alone[i] = c.errorNorm(err, y, yOut);
            h *= 0.5;
        }
        assertEquals("the combined estimate", 8.0, lastSlope(combined), 0.1);
        assertEquals("the first estimate alone", 6.0, lastSlope(alone), 0.1);
    }

    private static double lastSlope(double[] e) {
        return Math.log(e[e.length - 2] / e[e.length - 1]) / Math.log(2.0);
    }

    /**
     * And what that difference is worth in evaluations. The same runs with the
     * fifth order estimate alone -- which is what a tableau carrying only one
     * embedded solution would give -- cost {@code 1.42} times as much at
     * {@code rtol = 1e-06}, {@code 2.03} at {@code 1e-10} and {@code 4.44} at
     * {@code 1e-16}, because the step is being held to what a fifth order
     * solution can do while an eighth order one is what gets propagated.
     */
    @Test
    public void testControllingOnTheFirstEstimateAloneWouldCostSeveralTimesMore() {
        ButcherTableau d = ButcherTableau.DOP853;
        ButcherTableau one = new ButcherTableau("DOP853 with one estimate", 8, d.c(), d.a(), d.b(),
                d.bStar(), 5, d.dense());
        double[] tolerances = { 1.0e-6, 1.0e-10, 1.0e-14 };
        double[] least = { 1.3, 1.9, 3.3 };
        for (int i = 0; i < tolerances.length; ++i) {
            long both = orbitRun(d, tolerances[i]).evaluations;
            long single = orbitRun(one, tolerances[i]).evaluations;
            double ratio = single / (double) both;
            assertTrue("at rtol " + tolerances[i] + " one estimate cost " + ratio + " times as much",
                    ratio > least[i]);
        }
    }

    // ------------------------------------------------------------------
    // the crossing against Dormand-Prince
    // ------------------------------------------------------------------

    /**
     * Both sides of the crossing, because either alone would be a half truth.
     * At a loose tolerance the eighth order method is the more expensive way to
     * an answer nobody needed that badly.
     */
    @Test
    public void testDormandPrinceIsCheaperAtALooseToleranceAndDop853AtATightOne() {
        assertTrue("at rtol 1e-06 Dormand-Prince should still win",
                orbitRun(ButcherTableau.DORMAND_PRINCE_45, 1.0e-6).evaluations
                        < orbitRun(ButcherTableau.DOP853, 1.0e-6).evaluations);

        long dp = orbitRun(ButcherTableau.DORMAND_PRINCE_45, 1.0e-11).evaluations;
        long d8 = orbitRun(ButcherTableau.DOP853, 1.0e-11).evaluations;
        assertTrue("at rtol 1e-11 DOP853 should win, and by a lot: " + dp + " against " + d8,
                d8 * 2L < dp);
    }

    /**
     * The crossing again, as the shape rather than as two points: over seven
     * decades of tolerance the ratio of the two costs falls the whole way, from
     * {@code 0.93} in Dormand-Prince's favour to {@code 3.92} in DOP853's.
     */
    @Test
    public void testTheAdvantageGrowsTheTighterTheToleranceGets() {
        double rtol = 1.0e-6;
        double previous = 0.0;
        for (int i = 0; i < 8; ++i) {
            double ratio = orbitRun(ButcherTableau.DORMAND_PRINCE_45, rtol).evaluations
                    / (double) orbitRun(ButcherTableau.DOP853, rtol).evaluations;
            assertTrue("the ratio must grow, and at rtol " + rtol + " it went " + previous + " to " + ratio,
                    ratio > previous);
            previous = ratio;
            rtol *= 0.1;
        }
        assertTrue("and end well above one: " + previous, previous > 3.5);
    }

    /**
     * The two methods have to agree about the answer, whatever they cost to get
     * it. At {@code rtol = 1e-10} over ten orbits they agree to {@code 3.9e-08},
     * which is far below either one's distance from the truth.
     */
    @Test
    public void testTheTwoMethodsReachTheSamePlace() {
        double[] a = orbitRun(ButcherTableau.DORMAND_PRINCE_45, 1.0e-10).finalState();
        double[] b = orbitRun(ButcherTableau.DOP853, 1.0e-10).finalState();
        for (int i = 0; i < a.length; ++i) {
            assertEquals("component " + i, a[i], b[i], 1.0e-6);
        }
    }

    // ------------------------------------------------------------------
    // the interpolant
    // ------------------------------------------------------------------

    /**
     * The continuous extension is one order below the step it spans, which is
     * the same thing this package measures for Dormand-Prince and for the cubic
     * fallback: the worst error strictly inside a step falls by {@code 2^8.01},
     * {@code 2^8.00} and {@code 2^7.99} over halvings from {@code h = 0.8},
     * against {@code 2^5.0} for Dormand-Prince and {@code 2^4.0} for the cubic.
     * <p>
     * Published as an order seven interpolant, which is the same statement
     * counted the other way: an eighth order step has a local error of
     * <code>h^9</code>, and {@code h^8} in between is one order below it.
     * <p>
     * The window is narrow and sits high. One halving further the error is
     * {@code 1.2e-14} and the halving after that stops improving, so a range
     * chosen without looking would measure the arithmetic rather than the
     * method.
     */
    @Test
    public void testTheInterpolantIsOneOrderBelowTheStep() {
        double[] errors = new double[4];
        double h = 0.8;
        for (int i = 0; i < errors.length; ++i) {
            errors[i] = worstInsideAStep(h);
            h *= 0.5;
        }
        for (int i = 0; i + 1 < errors.length; ++i) {
            double slope = Math.log(errors[i] / errors[i + 1]) / Math.log(2.0);
            assertEquals("halving " + i, 8.0, slope, 0.15);
        }
    }

    /**
     * And what the three extra stages are worth: over the same steps the cubic
     * fallback, which is what a method with no continuous extension of its own
     * gets, is of order four. At {@code h = 0.2} that is {@code 4.2e-06}
     * against {@code 3.1e-12}.
     */
    @Test
    public void testTheThreeExtraStagesBuyFourOrdersOverTheCubicFallback() {
        double h = 0.2;
        double septic = worstInsideAStep(h);
        ExplicitRungeKutta cubic = new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2);
        double[] yOut = new double[2];
        double[] out = new double[2];
        cubic.step(0.0, new double[] { 1.0, 0.0 }, h, yOut, null);
        double worst = 0.0;
        for (int i = 1; i < 20; ++i) {
            double theta = i / 20.0;
            cubic.interpolate(theta, out);
            double t = theta * h;
            worst = Math.max(worst, Math.abs(out[0] - Math.cos(t)));
            worst = Math.max(worst, Math.abs(out[1] + Math.sin(t)));
        }
        assertTrue("septic " + septic + " against cubic " + worst, septic < 1.0e-6 * worst);
    }

    /**
     * However much accuracy the interpolant has in between, at the two ends it
     * is not an approximation at all: it is the state the step started from and
     * the state it reached, to the last bit.
     */
    @Test
    public void testTheEndsOfTheInterpolantAreTheEndsOfTheStepExactly() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DOP853, ORBIT, 4);
        double[] y = orbitStart();
        double[] yOut = new double[4];
        double[] err = new double[4];
        double[] out = new double[4];
        rk.step(0.0, y, 0.3, yOut, err);
        rk.interpolate(0.0, out);
        for (int i = 0; i < 4; ++i) {
            assertEquals(y[i], out[i], 0.0);
        }
        rk.interpolate(1.0, out);
        for (int i = 0; i < 4; ++i) {
            assertEquals(yOut[i], out[i], 0.0);
        }
    }

    private static double worstInsideAStep(double h) {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DOP853, OSCILLATOR, 2);
        double[] yOut = new double[2];
        double[] err = new double[2];
        double[] out = new double[2];
        rk.step(0.0, new double[] { 1.0, 0.0 }, h, yOut, err);
        double worst = 0.0;
        for (int i = 1; i < 20; ++i) {
            double theta = i / 20.0;
            rk.interpolate(theta, out);
            double t = theta * h;
            worst = Math.max(worst, Math.abs(out[0] - Math.cos(t)));
            worst = Math.max(worst, Math.abs(out[1] + Math.sin(t)));
        }
        return worst;
    }

    // ------------------------------------------------------------------
    // what it says about the equation
    // ------------------------------------------------------------------

    /**
     * The stiffness test was calibrated for a fifth order method and the number
     * it compares against is not universal: DOP853's stability region reaches
     * about twice as far, so its measure runs about twice as large on the same
     * equation and the threshold moves with it. On the family
     * <code>y' = -lambda (y - cos t) - sin t</code> the two methods measure
     * {@code 0.367} against {@code 1.142}, {@code 2.735} against {@code 5.619}
     * and {@code 4.919} against {@code 7.457} -- different numbers, and the
     * same three verdicts.
     */
    @Test
    public void testTheStiffnessVerdictSurvivesTheHigherThreshold() {
        double[] lambdas = { 1.0, 100.0, 10000.0 };
        boolean[] expected = { false, false, true };
        for (int i = 0; i < lambdas.length; ++i) {
            OdeIntegrator.Result dp = stiffnessRun(ButcherTableau.DORMAND_PRINCE_45, lambdas[i]);
            OdeIntegrator.Result d8 = stiffnessRun(ButcherTableau.DOP853, lambdas[i]);
            assertEquals("Dormand-Prince at lambda " + lambdas[i], expected[i], dp.seemsStiff);
            assertEquals("DOP853 at lambda " + lambdas[i], expected[i], d8.seemsStiff);
            assertTrue("DOP853 takes longer steps and so measures more: " + d8.stiffness + " against "
                    + dp.stiffness, d8.stiffness > dp.stiffness);
        }
    }

    @Test
    public void testTheThresholdReachesTheDriverThroughTheStepper() {
        assertEquals(6.1, new ExplicitRungeKutta(ButcherTableau.DOP853, ORBIT, 4).stiffnessThreshold(),
                0.0);
        assertEquals(3.25,
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, ORBIT, 4).stiffnessThreshold(),
                0.0);
        // and a method that cannot measure stiffness answers the default, which
        // costs nothing because its measure is never a number
        assertEquals(3.25, new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, ORBIT, 4)
                .stiffnessThreshold(), 0.0);
    }

    /**
     * The controller's constants were chosen for Dormand-Prince, and Hairer
     * runs {@code dop853} with different ones -- a floor of {@code 0.333}
     * rather than {@code 0.2}, a ceiling of {@code 6} rather than {@code 10},
     * and no proportional-integral term at all. Measured over four tolerances
     * on the orbit the two are within a tenth of each other in cost, and the
     * defaults here reject far fewer steps: 50 against 121 in total. So they
     * stay, and nothing here gains a second set of constructor arguments.
     */
    @Test
    public void testTheDefaultsSuitTheEighthOrderMethodToo() {
        double rtol = 1.0e-8;
        long ourRejections = 0L;
        long hisRejections = 0L;
        for (int i = 0; i < 4; ++i) {
            OdeIntegrator.Result ours = new OdeIntegrator(
                    new ExplicitRungeKutta(ButcherTableau.DOP853, ORBIT, 4),
                    new StepController(rtol, rtol)).solve(0.0, orbitStart(), 20.0);
            OdeIntegrator.Result his = new OdeIntegrator(
                    new ExplicitRungeKutta(ButcherTableau.DOP853, ORBIT, 4),
                    new StepController(rtol, rtol, 0.9, 0.333, 6.0, 0.0, Double.POSITIVE_INFINITY, 100000))
                            .solve(0.0, orbitStart(), 20.0);
            double ratio = ours.evaluations / (double) his.evaluations;
            assertTrue("at rtol " + rtol + " the two cost " + ours.evaluations + " and " + his.evaluations,
                    ratio > 0.85 && ratio < 1.15);
            ourRejections += ours.rejected;
            hisRejections += his.rejected;
            rtol *= 0.01;
        }
        assertTrue("the damped controller should reject fewer: " + ourRejections + " against "
                + hisRejections, ourRejections < hisRejections);
    }

    /**
     * A smooth equation at a tight tolerance is what the method is for, and it
     * must not report the equation as the problem.
     */
    @Test
    public void testATightToleranceIsNotMistakenForStiffness() {
        OdeIntegrator.Result r = orbitRun(ButcherTableau.DOP853, 1.0e-12);
        assertFalse(r.seemsStiff);
        assertTrue("well inside the stability region: " + r.stiffness, r.stiffness < 6.1);
    }

    private static OdeIntegrator.Result orbitRun(ButcherTableau tableau, double rtol) {
        return new OdeIntegrator(new ExplicitRungeKutta(tableau, ORBIT, 4),
                new StepController(rtol, rtol, 400000)).solve(0.0, orbitStart(), 20.0);
    }

    private static OdeIntegrator.Result stiffnessRun(ButcherTableau tableau, final double lambda) {
        DVectorField f = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = -lambda * (y[0] - Math.cos(t)) - Math.sin(t);
            }
        };
        return new OdeIntegrator(new ExplicitRungeKutta(tableau, f, 1),
                new StepController(1.0e-6, 1.0e-6, 400000)).solve(0.0, new double[] { 1.0 }, 10.0);
    }
}
