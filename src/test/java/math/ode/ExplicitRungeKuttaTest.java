package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorField;

/**
 * One step of {@link ExplicitRungeKutta}, against the harmonic oscillator,
 * whose solution is a cosine and needs no reference implementation.
 * <p>
 * What is asserted here is the <em>order</em> of a single step rather than any
 * particular error: halving the step size divides the local error by two to the
 * order plus one, and that ratio is the statement the coefficients of a method
 * make about themselves. Measured over five halvings from {@code h = 0.4}, the
 * classical method gives {@code 5.00, 5.00, 5.00, 5.00, 5.00} and
 * Dormand-Prince {@code 6.06, 6.01, 6.00, 6.00}; its error estimate falls one
 * order slower, at {@code 5.01, 5.00, 5.00, 5.00}, which is what it is for --
 * an estimate that shrank as fast as the error it estimates would be measuring
 * the wrong thing.
 * <p>
 * The interpolation is one order below the step it spans in both cases: the
 * continuous extension of Dormand-Prince at {@code 5.09, 5.02, 5.01, 5.00,
 * 5.00} against a method of order five, and the cubic fallback at {@code 4.09,
 * 4.02, 4.01, 4.00} against a method of order four.
 */
public final class ExplicitRungeKuttaTest {

    /** y'' = -y as a first order pair, exact solution (cos t, -sin t) at (1, 0). */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    private static final double T0 = 0.3;

    @Test
    public void testTheLocalErrorOfTheClassicalMethodIsOfFifthOrder() {
        double[] slopes = slopes(localErrors(ButcherTableau.CLASSIC_RK4, 0.4, 6, false));
        assertSlopes(slopes, 5.0, 0.1);
    }

    @Test
    public void testTheLocalErrorOfDormandPrinceIsOfSixthOrder() {
        double[] slopes = slopes(localErrors(ButcherTableau.DORMAND_PRINCE_45, 0.4, 5, false));
        assertSlopes(slopes, 6.0, 0.2);
    }

    /**
     * The embedded pair estimates the error of the <em>fourth</em> order
     * solution, so it falls one order slower than the fifth order step it is
     * attached to and is a safe overestimate of it. Measured, it exceeds the
     * true local error by a factor between 6.9 and 116 over this range, and the
     * factor grows as the step shrinks precisely because the two fall at
     * different rates.
     */
    @Test
    public void testTheErrorEstimateFallsOneOrderSlowerThanTheErrorItEstimates() {
        double[] estimates = localErrors(ButcherTableau.DORMAND_PRINCE_45, 0.4, 5, true);
        double[] errors = localErrors(ButcherTableau.DORMAND_PRINCE_45, 0.4, 5, false);
        assertSlopes(slopes(estimates), 5.0, 0.1);
        for (int i = 0; i < errors.length; ++i) {
            assertTrue("the estimate must not fall below the error it estimates",
                    estimates[i] > errors[i]);
        }
    }

    @Test
    public void testTheContinuousExtensionOfDormandPrinceIsOneOrderBelowTheStep() {
        double[] slopes = slopes(interiorErrors(ButcherTableau.DORMAND_PRINCE_45, 0.4, 6));
        assertSlopes(slopes, 5.0, 0.2);
    }

    @Test
    public void testTheCubicFallbackIsOneOrderBelowTheClassicalStep() {
        double[] slopes = slopes(interiorErrors(ButcherTableau.CLASSIC_RK4, 0.2, 5));
        assertSlopes(slopes, 4.0, 0.2);
    }

    /**
     * An output point that falls on a step boundary must not disagree with the
     * step itself, so the two ends of the interpolant are copies rather than
     * evaluations of the polynomial, which would land a few units in the last
     * place away.
     */
    @Test
    public void testTheEndsOfTheInterpolantAreTheEndsOfTheStep() {
        ButcherTableau[] both = { ButcherTableau.CLASSIC_RK4, ButcherTableau.DORMAND_PRINCE_45 };
        for (int q = 0; q < both.length; ++q) {
            ExplicitRungeKutta rk = new ExplicitRungeKutta(both[q], OSCILLATOR, 2);
            double[] y = { 1.0, 0.0 };
            double[] yOut = new double[2];
            double[] out = new double[2];
            rk.step(0.0, y, 0.25, yOut, null);
            rk.interpolate(0.0, out);
            assertEquals(both[q].name(), y[0], out[0], 0.0);
            assertEquals(both[q].name(), y[1], out[1], 0.0);
            rk.interpolate(1.0, out);
            assertEquals(both[q].name(), yOut[0], out[0], 0.0);
            assertEquals(both[q].name(), yOut[1], out[1], 0.0);
        }
    }

    @Test
    public void testAnAcceptedDormandPrinceStepCostsOneEvaluationLessThanItHasStages() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        double t = 0.0;
        for (int i = 0; i < 10; ++i) {
            rk.step(t, y, 0.1, yOut, null);
            t += 0.1;
            System.arraycopy(yOut, 0, y, 0, 2);
        }
        // seven for the first step, six for each of the nine that follow
        assertEquals(61L, rk.evaluations());
    }

    @Test
    public void testTheClassicalMethodPaysForEveryStageOfEveryStep() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        double t = 0.0;
        for (int i = 0; i < 10; ++i) {
            rk.step(t, y, 0.1, yOut, null);
            t += 0.1;
            System.arraycopy(yOut, 0, y, 0, 2);
        }
        assertEquals(40L, rk.evaluations());
    }

    /**
     * A rejected step is retried from the point it started at, where the first
     * stage is what it already was. This holds for any method and not only for
     * a first-same-as-last one, and it is recognized by comparing the point
     * rather than by the driver announcing what it did.
     */
    @Test
    public void testARetryFromTheSamePointDoesNotEvaluateTheFirstStageTwice() {
        ButcherTableau[] both = { ButcherTableau.CLASSIC_RK4, ButcherTableau.DORMAND_PRINCE_45 };
        for (int q = 0; q < both.length; ++q) {
            ExplicitRungeKutta rk = new ExplicitRungeKutta(both[q], OSCILLATOR, 2);
            double[] y = { 1.0, 0.0 };
            double[] yOut = new double[2];
            rk.step(0.0, y, 0.2, yOut, null);
            long afterFirst = rk.evaluations();
            rk.step(0.0, y, 0.1, yOut, null);
            assertEquals(both[q].name(), both[q].stages() - 1L, rk.evaluations() - afterFirst);
        }
    }

    /**
     * The derivative at the end of the step that the cubic fallback needs is
     * evaluated when an interior point is first asked for, and not at all
     * otherwise.
     */
    @Test
    public void testTheCubicFallbackCostsOneEvaluationAndOnlyWhenItIsUsed() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        double[] out = new double[2];
        rk.step(0.0, y, 0.2, yOut, null);
        assertEquals(4L, rk.evaluations());
        rk.interpolate(0.5, out);
        assertEquals(5L, rk.evaluations());
        rk.interpolate(0.25, out);
        assertEquals("the second interior point of the same step is free", 5L, rk.evaluations());
    }

    /**
     * Dormand-Prince pays nothing for an interior point: its continuous
     * extension uses the seven stages the step computed anyway.
     */
    @Test
    public void testTheContinuousExtensionCostsNothing() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        double[] out = new double[2];
        rk.step(0.0, y, 0.2, yOut, null);
        long afterStep = rk.evaluations();
        for (int i = 1; i < 20; ++i) {
            rk.interpolate(i / 20.0, out);
        }
        assertEquals(afterStep, rk.evaluations());
    }

    @Test
    public void testResetForgetsTheStepButNotTheCost() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        double[] out = new double[2];
        rk.step(0.0, y, 0.1, yOut, null);
        assertEquals(7L, rk.evaluations());
        rk.reset();
        assertEquals("the evaluation count runs across a reset", 7L, rk.evaluations());
        try {
            rk.interpolate(0.5, out);
            fail("interpolation after a reset must not answer");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("no step"));
        }
        rk.step(0.0, y, 0.1, yOut, null);
        assertEquals("the first stage is computed again after a reset", 14L, rk.evaluations());
    }

    @Test
    public void testInterpolationBeforeAnyStepIsRefused() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        try {
            rk.interpolate(0.5, new double[2]);
            fail("there is no step to interpolate");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("no step"));
        }
    }

    @Test
    public void testTheArgumentsOfAStepAreChecked() {
        ExplicitRungeKutta rk = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        double[] y = { 1.0, 0.0 };
        double[] yOut = new double[2];
        refuse(rk, 0.0, new double[3], 0.1, yOut, null, "y");
        refuse(rk, 0.0, y, 0.1, new double[3], null, "yOut");
        refuse(rk, 0.0, y, 0.1, y, null, "must not be the array passed as y");
        refuse(rk, 0.0, y, 0.1, yOut, new double[3], "errOut");
        refuse(rk, 0.0, y, 0.0, yOut, null, "h must not be zero");
    }

    @Test
    public void testTheStepperReportsWhatItsTableauCan() {
        ExplicitRungeKutta dp = new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
        assertSame(ButcherTableau.DORMAND_PRINCE_45, dp.tableau());
        assertEquals(2, dp.dimension());
        assertEquals(5, dp.order());
        assertTrue(dp.hasErrorEstimate());
        assertTrue(dp.hasDenseOutput());

        ExplicitRungeKutta rk4 = new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 4);
        assertEquals(4, rk4.dimension());
        assertEquals(4, rk4.order());
        assertFalse(rk4.hasErrorEstimate());
        assertFalse(rk4.hasDenseOutput());
    }

    @Test
    public void testTheConstructorChecksItsArguments() {
        try {
            new ExplicitRungeKutta(null, OSCILLATOR, 2);
            fail("expected a refusal naming the tableau");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("tableau"));
        }
        try {
            new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, null, 2);
            fail("expected a refusal naming the field");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("field"));
        }
        try {
            new ExplicitRungeKutta(ButcherTableau.CLASSIC_RK4, OSCILLATOR, 0);
            fail("expected a refusal naming the dimension");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("dimension"));
        }
    }

    private static void refuse(ExplicitRungeKutta rk, double t, double[] y, double h, double[] yOut,
            double[] errOut, String hint) {
        try {
            rk.step(t, y, h, yOut, errOut);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }

    /**
     * The error of one step from {@code T0}, or of the embedded estimate of it,
     * at a step size halved {@code count} times from {@code first}.
     */
    private static double[] localErrors(ButcherTableau tableau, double first, int count, boolean estimate) {
        double[] out = new double[count];
        for (int i = 0; i < count; ++i) {
            double h = first / Math.pow(2.0, i);
            ExplicitRungeKutta rk = new ExplicitRungeKutta(tableau, OSCILLATOR, 2);
            double[] y = { Math.cos(T0), -Math.sin(T0) };
            double[] yOut = new double[2];
            double[] err = new double[2];
            rk.step(T0, y, h, yOut, err);
            if (estimate) {
                out[i] = Math.hypot(err[0], err[1]);
            } else {
                out[i] = Math.hypot(yOut[0] - Math.cos(T0 + h), yOut[1] + Math.sin(T0 + h));
            }
        }
        return out;
    }

    /** The worst interpolation error strictly inside one step, over 63 points. */
    private static double[] interiorErrors(ButcherTableau tableau, double first, int count) {
        double[] out = new double[count];
        for (int i = 0; i < count; ++i) {
            double h = first / Math.pow(2.0, i);
            ExplicitRungeKutta rk = new ExplicitRungeKutta(tableau, OSCILLATOR, 2);
            double[] y = { Math.cos(T0), -Math.sin(T0) };
            double[] yOut = new double[2];
            double[] point = new double[2];
            rk.step(T0, y, h, yOut, null);
            double worst = 0.0;
            for (int j = 1; j < 64; ++j) {
                double theta = j / 64.0;
                rk.interpolate(theta, point);
                double t = T0 + theta * h;
                worst = Math.max(worst, Math.hypot(point[0] - Math.cos(t), point[1] + Math.sin(t)));
            }
            out[i] = worst;
        }
        return out;
    }

    private static double[] slopes(double[] errors) {
        double[] out = new double[errors.length - 1];
        for (int i = 1; i < errors.length; ++i) {
            out[i - 1] = Math.log(errors[i - 1] / errors[i]) / Math.log(2.0);
        }
        return out;
    }

    private static void assertSlopes(double[] slopes, double expected, double tolerance) {
        for (int i = 0; i < slopes.length; ++i) {
            assertEquals("halving " + (i + 1) + " of " + slopes.length, expected, slopes[i], tolerance);
        }
    }
}
