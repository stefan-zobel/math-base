package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorField;

/**
 * {@link StepController} on its own, without an equation being solved.
 * <p>
 * Everything here is a property rather than a number out of a run: a step at
 * exactly the tolerance is shortened rather than kept at its length, a smaller
 * error buys a longer step, the factor never leaves the band it was given, and
 * a rejected step never grows. The one number is the textbook controller's
 * answer to an error of exactly one, which is the safety factor and nothing
 * else -- that pins the formula without pinning the equation.
 */
public final class StepControllerTest {

    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    @Test
    public void testTheErrorNormIsTheRootMeanSquareOfTheScaledComponents() {
        StepController c = new StepController(1.0e-3, 1.0e-3);
        double[] error = { 3.0e-3, 4.0e-3 };
        double[] y = { 1.0, 1.0 };
        // each component is divided by atol + rtol * 1 = 2e-3, giving 1.5 and 2
        assertEquals(Math.sqrt((1.5 * 1.5 + 2.0 * 2.0) / 2.0), c.errorNorm(error, y, y), 1.0e-15);
    }

    /**
     * A second estimate of exactly zero has to leave the first one alone, and
     * it does so identically rather than nearly: the combination is
     * <code>s1 sqrt(1 / (n (s1 + 0.01 s2)))</code>, which at {@code s2 = 0}
     * is <code>sqrt(s1 / n)</code>, the plain root mean square.
     */
    @Test
    public void testTheCombinedNormFallsBackOnThePlainOneWhenTheSecondEstimateVanishes() {
        StepController c = new StepController(1.0e-3, 1.0e-3);
        double[] error = { 3.0e-3, 4.0e-3 };
        double[] zero = { 0.0, 0.0 };
        double[] y = { 1.0, 1.0 };
        assertEquals(c.errorNorm(error, y, y), c.errorNorm(error, zero, y, y), 1.0e-15);
    }

    /**
     * And the worked example in the other direction: the same two components,
     * scaled to 1.5 and 2, against a second estimate scaled to 15 and 20. Then
     * {@code s1} is {@code 6.25}, {@code s2} is {@code 625}, and the answer is
     * {@code 6.25 / sqrt(2 * 12.5)}, which is a good deal below the
     * {@code 1.7678} the first estimate alone would have given.
     */
    @Test
    public void testTheCombinedNormIsHairersCombination() {
        StepController c = new StepController(1.0e-3, 1.0e-3);
        double[] error = { 3.0e-3, 4.0e-3 };
        double[] secondary = { 3.0e-2, 4.0e-2 };
        double[] y = { 1.0, 1.0 };
        double s1 = 1.5 * 1.5 + 2.0 * 2.0;
        double s2 = 15.0 * 15.0 + 20.0 * 20.0;
        assertEquals(s1 * Math.sqrt(1.0 / (2.0 * (s1 + 0.01 * s2))), c.errorNorm(error, secondary, y, y),
                1.0e-15);
        assertEquals(1.7677669529663689, c.errorNorm(error, y, y), 1.0e-15);
    }

    /**
     * Two estimates that are both exactly zero would divide by zero, and the
     * guard Hairer writes turns that into a step accepted rather than a step
     * whose error is not a number.
     */
    @Test
    public void testTwoEstimatesOfZeroLeaveAnErrorOfZero() {
        StepController c = new StepController(1.0e-3, 1.0e-3);
        double[] zero = { 0.0, 0.0 };
        double[] y = { 1.0, 1.0 };
        assertEquals(0.0, c.errorNorm(zero, zero, y, y), 0.0);
    }

    @Test
    public void testTheLargerOfTheTwoStatesSetsTheScale() {
        StepController c = new StepController(1.0e-3, 1.0e-300);
        double[] error = { 1.0e-3 };
        assertEquals(1.0, c.errorNorm(error, new double[] { 1.0 }, new double[] { 0.0 }), 1.0e-12);
        assertEquals(1.0, c.errorNorm(error, new double[] { 0.0 }, new double[] { 1.0 }), 1.0e-12);
        assertEquals(0.5, c.errorNorm(error, new double[] { 2.0 }, new double[] { 1.0 }), 1.0e-12);
    }

    /**
     * The safety factor is what stops a step aimed exactly at the tolerance
     * from being rejected about half the time, and with the textbook controller
     * it is the whole of the answer at an error of one.
     */
    @Test
    public void testAStepAtExactlyTheToleranceIsShortenedByTheSafetyFactor() {
        StepController plain = new StepController(1.0e-6, 1.0e-6, 0.9, 0.2, 10.0, 0.0,
                Double.POSITIVE_INFINITY, 1000);
        assertEquals(0.9, plain.scale(1.0, 1.0, 5), 1.0e-15);
        assertEquals(0.9, plain.scale(1.0, 1.0e-8, 5), 1.0e-15);
    }

    @Test
    public void testASmallerErrorBuysALongerStep() {
        StepController c = new StepController(1.0e-6, 1.0e-6);
        double previous = 0.0;
        double[] errors = { 1.0, 0.5, 0.1, 0.01, 1.0e-4 };
        for (int i = 0; i < errors.length; ++i) {
            double factor = c.scale(errors[i], StepController.ERROR_FLOOR, 5);
            assertTrue("the factor must grow as the error shrinks", factor > previous);
            previous = factor;
        }
    }

    @Test
    public void testTheFactorNeverLeavesTheBandItWasGiven() {
        StepController c = new StepController(1.0e-6, 1.0e-6, 0.9, 0.25, 4.0, 0.04,
                Double.POSITIVE_INFINITY, 1000);
        assertEquals(4.0, c.scale(1.0e-30, StepController.ERROR_FLOOR, 5), 0.0);
        assertEquals(4.0, c.scale(0.0, StepController.ERROR_FLOOR, 5), 0.0);
        assertEquals(0.25, c.scale(1.0e30, StepController.ERROR_FLOOR, 5), 0.0);
        for (int i = -20; i <= 20; ++i) {
            double factor = c.scale(Math.pow(10.0, i), StepController.ERROR_FLOOR, 5);
            assertTrue(factor >= 0.25 && factor <= 4.0);
        }
    }

    @Test
    public void testARejectedStepNeverGrowsAndAnImpossibleOneShrinksAsFarAsItMay() {
        StepController c = new StepController(1.0e-6, 1.0e-6);
        for (int i = 0; i <= 20; ++i) {
            double factor = c.scaleAfterRejection(1.0 + i, 5);
            assertTrue("a rejected step must not grow", factor <= 1.0);
            assertTrue(factor >= c.minScale());
        }
        assertEquals(c.minScale(), c.scaleAfterRejection(Double.NaN, 5), 0.0);
        assertEquals(c.minScale(), c.scaleAfterRejection(Double.POSITIVE_INFINITY, 5), 0.0);
    }

    @Test
    public void testTheFirstStepPointsAtTheTargetAndStaysInsideTheInterval() {
        StepController c = new StepController(1.0e-8, 1.0e-8);
        OdeStepper forward = stepper();
        double hf = c.initialStep(forward, 0.0, new double[] { 1.0, 0.0 }, 5.0);
        assertTrue("a forward step is positive", hf > 0.0);
        assertTrue("and does not overshoot the interval", hf <= 5.0);

        OdeStepper backward = stepper();
        double hb = c.initialStep(backward, 0.0, new double[] { 1.0, 0.0 }, -5.0);
        assertTrue("a backward step is negative", hb < 0.0);
        assertTrue(hb >= -5.0);
        assertEquals("the direction is the only difference", hf, -hb, 1.0e-15);
    }

    @Test
    public void testTheFirstStepCostsTwoEvaluationsAndIsOneTheMethodWouldKeep() {
        StepController c = new StepController(1.0e-8, 1.0e-8);
        OdeStepper s = stepper();
        double h = c.initialStep(s, 0.0, new double[] { 1.0, 0.0 }, 5.0);
        assertEquals(2L, s.evaluations());

        double[] y = { 1.0, 0.0 };
        double[] yNext = new double[2];
        double[] error = new double[2];
        s.step(0.0, y, h, yNext, error);
        assertTrue("the guess should not be one the controller throws away straight away",
                c.errorNorm(error, y, yNext) <= 1.0);
    }

    @Test
    public void testTheFirstStepRespectsAnUpperBoundOnTheStepSize() {
        StepController bounded = new StepController(1.0e-8, 1.0e-8, 0.9, 0.2, 10.0, 0.04, 1.0e-3, 1000);
        double h = bounded.initialStep(stepper(), 0.0, new double[] { 1.0, 0.0 }, 5.0);
        assertTrue(h > 0.0 && h <= 1.0e-3);
    }

    @Test
    public void testTheDefaultsAreWhatTheyClaim() {
        StepController c = new StepController(1.0e-6, 1.0e-9);
        assertEquals(1.0e-6, c.relativeTolerance(), 0.0);
        assertEquals(1.0e-9, c.absoluteTolerance(), 0.0);
        assertEquals(0.9, c.safety(), 0.0);
        assertEquals(0.2, c.minScale(), 0.0);
        assertEquals(10.0, c.maxScale(), 0.0);
        assertEquals(0.04, c.beta(), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, c.maxStepSize(), 0.0);
        assertEquals(100000, c.maxSteps());
        assertEquals(500, new StepController(1.0e-6, 1.0e-9, 500).maxSteps());
        assertTrue(c.toString().contains("rtol"));
    }

    @Test
    public void testTheSettingsAreChecked() {
        refuse(0.0, 1.0e-6, "rtol");
        refuse(-1.0, 1.0e-6, "rtol");
        refuse(Double.POSITIVE_INFINITY, 1.0e-6, "rtol");
        refuse(1.0e-6, 0.0, "atol");
        refuse(1.0e-6, Double.NaN, "atol");
        refuseFull(0.0, 0.2, 10.0, 0.04, 1.0, 10, "safety");
        refuseFull(1.5, 0.2, 10.0, 0.04, 1.0, 10, "safety");
        refuseFull(0.9, 1.0, 10.0, 0.04, 1.0, 10, "minScale");
        refuseFull(0.9, 0.2, 0.5, 0.04, 1.0, 10, "maxScale");
        refuseFull(0.9, 0.2, 10.0, -0.1, 1.0, 10, "beta");
        refuseFull(0.9, 0.2, 10.0, 0.6, 1.0, 10, "beta");
        refuseFull(0.9, 0.2, 10.0, 0.04, 0.0, 10, "maxStepSize");
        refuseFull(0.9, 0.2, 10.0, 0.04, 1.0, 0, "maxSteps");
    }

    @Test
    public void testTheArgumentsOfTheFirstStepAreChecked() {
        StepController c = new StepController(1.0e-8, 1.0e-8);
        try {
            c.initialStep(null, 0.0, new double[] { 1.0, 0.0 }, 1.0);
            fail("expected a refusal naming the stepper");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("stepper"));
        }
        try {
            c.initialStep(stepper(), 0.0, new double[] { 1.0 }, 1.0);
            fail("expected a refusal naming y0");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("y0"));
        }
        try {
            c.initialStep(stepper(), 1.0, new double[] { 1.0, 0.0 }, 1.0);
            fail("expected a refusal about an empty interval");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("must differ"));
        }
    }

    private static OdeStepper stepper() {
        return new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, OSCILLATOR, 2);
    }

    private static void refuse(double rtol, double atol, String hint) {
        try {
            new StepController(rtol, atol);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }

    private static void refuseFull(double safety, double minScale, double maxScale, double beta,
            double maxStepSize, int maxSteps, String hint) {
        try {
            new StepController(1.0e-6, 1.0e-6, safety, minScale, maxScale, beta, maxStepSize, maxSteps);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }
}
