package math.solve;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import math.fun.DFunction;
import math.fun.DBiFunction;
import math.fun.DTriFunction;

/**
 * Unit tests for the core AdaptiveGaussKronrod engine. Validates error
 * estimating qualities and recursive sub-divisions.
 */
public class AdaptiveGaussKronrodTest {

    private static AdaptiveGaussKronrod.G7_K15 ruleSetup;
    private static final double STRICT_TOLERANCE = 1e-9;
    private static final double LOOSE_TOLERANCE = 1e-4;

    @BeforeClass
    public static void setUpBeforeClass() {
        ruleSetup = AdaptiveGaussKronrod.G7_K15.POINTS_15;
    }

    // =========================================================================
    // 1D ENGINE TESTING
    // =========================================================================

    @Test
    public void test1DRawGKStepOnSmoothPolynomial() {
        // f(x) = x^4 over [0, 1]. Exact analytical value = 1/5 = 0.2
        // G7-K15 is exact for polynomials up to degree 29, so the error estimate must be ~0.
        DFunction poly1D = x -> x * x * x * x;
        double exactValue = 0.2;

        AdaptiveGaussKronrod.IntegralResult result =
            AdaptiveGaussKronrod.integrate1D(ruleSetup, poly1D, 0.0, 1.0);

        assertNotNull("IntegralResult must not be null", result);
        assertEquals("1D G7-K15 single-step failed for smooth polynomial", exactValue, result.value, STRICT_TOLERANCE);
        assertTrue("Error estimate for degree-4 polynomial must be effectively zero (G7-K15 is exact up to degree 29)",
                   result.approximatedErrorEstimate < 1e-12);
    }

    @Test
    public void test1DAdaptiveSubdivisionOnNarrowGaussianSpike() {
        // Narrow Gaussian spike at x=0.37, positioned between G7-K15 quadrature nodes on [0, 1].
        // The single G7-K15 step severely underestimates the spike, so the error estimate is large
        // and multiple subdivisions must be triggered.
        // Regression for the integrate1DAdaptive bug: the old stopping criterion used
        // res.value <= epsTol instead of res.approximatedErrorEstimate <= epsTol. A function
        // whose GK approximation is accidentally small would have exited prematurely with
        // the wrong answer; this test verifies convergence to the correct value.
        DFunction spike1D = x -> Math.exp(-5000.0 * Math.pow(x - 0.37, 2));
        // True integral over [0,1] ? sqrt(PI/5000): boundary contributions are negligible
        // since exp(-5000*(0-0.37)^2) = exp(-684) ? 0
        double exactValue = Math.sqrt(Math.PI / 5000.0);

        double result = AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, spike1D, 0.0, 1.0, 1e-7, 25);

        assertEquals("Adaptive 1D subdivision failed: narrow Gaussian spike not resolved "
                     + "(regression: error-estimate vs value stopping criterion)",
                     exactValue, result, LOOSE_TOLERANCE);
    }

    // =========================================================================
    // 2D ENGINE TESTING
    // =========================================================================

    @Test
    public void test2DParallelRawStepAndErrorEstimation() {
        // f(x, y) = x * y over [0, 1] x [0, 1]. Exact analytical value = 0.25
        DBiFunction linear2D = (x, y) -> x * y;
        double exactValue = 0.25;

        AdaptiveGaussKronrod.IntegralResult result = 
            AdaptiveGaussKronrod.integrate2DParallel(ruleSetup, linear2D, 0.0, 1.0, 0.0, 1.0);

        assertNotNull("IntegralResult container must not be null", result);
        assertEquals("Raw parallel 2D numerical approximation calculation failed", exactValue, result.value, STRICT_TOLERANCE);

        // Linear functions should have virtually 0 error difference between G7 and K15
        assertTrue("Error estimate for simple linear functions should be extremely low", result.approximatedErrorEstimate < 1e-12);
    }

    @Test
    public void test2DAdaptiveSubdivisionOnSharpSpike() {
        // Highly localized Gaussian distribution peak function centered at (0.5, 0.5)
        DBiFunction spike2D = (x, y) -> Math.exp(-150.0 * (Math.pow(x - 0.5, 2) + Math.pow(y - 0.5, 2)));
        double expectedAnalyticalLimit = Math.PI / 150.0; // ~0.02094395

        // Let the adaptive engine isolate the spike natively
        double result = AdaptiveGaussKronrod.integrate2DAdaptive(ruleSetup, spike2D, 0.0, 1.0, 0.0, 1.0, 1e-6, 10);

        assertEquals("Adaptive 2D subdivision engine missed or miscalculated the localized peak", 
                     expectedAnalyticalLimit, result, LOOSE_TOLERANCE);
    }

    // =========================================================================
    // 3D ENGINE TESTING
    // =========================================================================

    @Test
    public void test3DParallelRawStepValue() {
        // f(x, y, z) = x + y + z over [0, 1]^3. Exact analytical value = 1.5
        DTriFunction plane3D = (x, y, z) -> x + y + z;
        double exactValue = 1.5;

        AdaptiveGaussKronrod.IntegralResult result = 
            AdaptiveGaussKronrod.integrate3DParallel(ruleSetup, plane3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

        assertNotNull("IntegralResult container must not be null", result);
        assertEquals("Raw parallel 3D calculation error detected", exactValue, result.value, STRICT_TOLERANCE);
    }

    @Test
    public void test3DAdaptiveSubdivisionOnNonSeparableField() {
        DTriFunction nonSeparable3D = (x, y, z) -> Math.exp(x + y) * (z * z + x * y);
        double exactValue = (1.0 / 3.0 * (Math.E * Math.E)) - (2.0 / 3.0 * Math.E) + (4.0 / 3.0); // ~1.98416

        // Force subdivision on the cubic framework using a tighter error tolerance goal
        double result = AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, nonSeparable3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-7, 8);

        assertEquals("Adaptive 3D subdivision failed on non-separable variable mappings", exactValue, result, 1e-5);
    }

    @Test
    public void test3DAdaptiveSubdivisionOnSharpPeak() {
        // Complex sharp 3D spatial peak function centered at (0.5, 0.5, 0.5)
        DTriFunction peak3D = (x, y, z) -> Math.exp(-100.0 * (Math.pow(x - 0.5, 2) + Math.pow(y - 0.5, 2) + Math.pow(z - 0.5, 2)));
        // Exact reference = integral over all of R^3 = pi^(3/2) / 100^(3/2) = pi^(3/2) / 1000
        // (valid approximation: peak is so narrow that mass outside [0,1]^3 is negligible)
        double expectedValue = Math.pow(Math.PI, 1.5) / 1000.0;

        double result = AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, peak3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-6, 8);

        assertEquals("Adaptive 3D subdivision engine failed to correctly isolate localized volume peak",
                     expectedValue, result, LOOSE_TOLERANCE);
    }

    @Test
    public void testNarrowPeakBetweenTheNodesIsNotMissed() {
        // A peak of width 1/sqrt(2*20000) ~ 0.005 centered at 0.53 falls between
        // the Kronrod nodes of the undivided domain. Both the Kronrod and the
        // Gauss rule then return ~0, agree, and |K - G| reports an error of
        // ~1e-26 for a value that is wrong by 100 percent. Without the forced
        // initial subdivision the 3D case returned 4.06e-27 instead of 1.97e-06.
        // The peak sits 66 standard deviations from either boundary, so the
        // integral over the unit domain equals the one over all of R.
        final double sharpness = 20000.0;
        final double center = 0.53;
        double oneAxis = Math.sqrt(Math.PI / sharpness);

        DFunction peak1D = x -> Math.exp(-sharpness * (x - center) * (x - center));
        double result1D = AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, peak1D, 0.0, 1.0, 1e-8, 20);
        assertEquals("1D: narrow peak between the nodes was missed", oneAxis, result1D, 1e-9);

        DBiFunction peak2D = (x, y) -> Math.exp(-sharpness
                * ((x - center) * (x - center) + (y - center) * (y - center)));
        double result2D = AdaptiveGaussKronrod.integrate2DAdaptive(ruleSetup, peak2D, 0.0, 1.0, 0.0, 1.0, 1e-8, 12);
        assertEquals("2D: narrow peak between the nodes was missed", oneAxis * oneAxis, result2D, 1e-9);

        // 3D needs a recursion budget beyond the 8 that MetaIntegrator passes:
        // with 6 of 8 levels forced, too little of it is left for the adaptive
        // part to resolve a peak this narrow
        DTriFunction peak3D = (x, y, z) -> Math.exp(-sharpness
                * ((x - center) * (x - center) + (y - center) * (y - center) + (z - center) * (z - center)));
        double expected3D = oneAxis * oneAxis * oneAxis;
        double result3D = AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, peak3D,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8, 14);
        assertEquals("3D: narrow peak between the nodes was missed", expected3D, result3D, 1e-2 * expected3D);
    }

    @Test
    public void testExplicitForcedSubdivisionsClosesTheResidual3DRisk() {
        // 0.38 is one of the peak positions that the 3D default of three forced
        // levels still misses: it returns ~0 and reports convergence. Six
        // levels, two bisections per axis, resolve it.
        final double sharpness = 20000.0;
        final double center = 0.38;
        double oneAxis = Math.sqrt(Math.PI / sharpness);
        double expected = oneAxis * oneAxis * oneAxis;
        DTriFunction peak3D = (x, y, z) -> Math.exp(-sharpness
                * ((x - center) * (x - center) + (y - center) * (y - center) + (z - center) * (z - center)));

        double byDefault = AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, peak3D,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8, 14);
        assertTrue("the default is supposed to miss this one, otherwise the test has lost its point",
                   Math.abs(byDefault - expected) > 0.5 * expected);

        double raised = AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, peak3D,
                0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-8, 14, 6);
        assertEquals("six forced levels must resolve the peak the default misses",
                     expected, raised, 1e-5 * expected);
    }

    @Test
    public void testExplicitForcedSubdivisionsAgreesWithTheDefaultAndClamps() {
        DFunction smooth = x -> Math.exp(x);
        double byDefault = AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 20);

        // passing the default explicitly must change nothing at all
        assertEquals("explicit default must reproduce the default bit for bit", byDefault,
                     AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 20, 3), 0.0);
        // a negative count is clamped to zero, i.e. to no forced subdivision
        assertEquals("negative forced count must behave like zero",
                     AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 20, 0),
                     AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 20, -5), 0.0);
        // a count above the recursion budget is clamped to the budget
        assertEquals("forced count above maxDepth must behave like maxDepth",
                     AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 4, 4),
                     AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, smooth, 0.0, 1.0, 1e-10, 4, 999), 0.0);
        // and all of them still integrate correctly
        assertEquals("forced subdivisions must not change the value on a smooth integrand",
                     Math.E - 1.0, byDefault, STRICT_TOLERANCE);
    }

    /**
     * The {@code WithError} forms are the same subdivision, reading off what
     * the plain forms drop: the values agree to the bit, because one is the
     * other's {@code value}.
     */
    @Test
    public void theErrorFormsAreTheSameSubdivision() {
        DFunction f = x -> Math.exp(-x * x) + 0.25 * x;
        assertEquals("1D", AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, f, -3.0, 3.0, 1.0e-10, 20),
                AdaptiveGaussKronrod.integrate1DAdaptiveWithError(ruleSetup, f, -3.0, 3.0, 1.0e-10, 20).value, 0.0);
        assertEquals("1D, forced levels named",
                AdaptiveGaussKronrod.integrate1DAdaptive(ruleSetup, f, -3.0, 3.0, 1.0e-10, 20, 6),
                AdaptiveGaussKronrod.integrate1DAdaptiveWithError(ruleSetup, f, -3.0, 3.0, 1.0e-10, 20, 6).value, 0.0);

        DBiFunction g = (x, y) -> Math.exp(-(x * x + y * y));
        assertEquals("2D", AdaptiveGaussKronrod.integrate2DAdaptive(ruleSetup, g, -3.0, 3.0, -3.0, 3.0, 1.0e-9, 10),
                AdaptiveGaussKronrod.integrate2DAdaptiveWithError(ruleSetup, g, -3.0, 3.0, -3.0, 3.0, 1.0e-9, 10).value,
                0.0);
        assertEquals("2D, forced levels named",
                AdaptiveGaussKronrod.integrate2DAdaptive(ruleSetup, g, -3.0, 3.0, -3.0, 3.0, 1.0e-9, 10, 6),
                AdaptiveGaussKronrod.integrate2DAdaptiveWithError(ruleSetup, g, -3.0, 3.0, -3.0, 3.0, 1.0e-9, 10, 6).value,
                0.0);

        DTriFunction h = (x, y, z) -> Math.exp(-(x * x + y * y + z * z));
        assertEquals("3D",
                AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, h, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, 1.0e-6, 14),
                AdaptiveGaussKronrod.integrate3DAdaptiveWithError(ruleSetup, h, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, 1.0e-6,
                        14).value, 0.0);
        assertEquals("3D, forced levels named",
                AdaptiveGaussKronrod.integrate3DAdaptive(ruleSetup, h, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, 1.0e-6, 14, 6),
                AdaptiveGaussKronrod.integrate3DAdaptiveWithError(ruleSetup, h, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, 1.0e-6,
                        14, 6).value, 0.0);
    }

    /**
     * A budget too small to meet the tolerance is the case the flag exists for:
     * a panel is handed out without meeting the tolerance it was given, and the
     * result says so instead of the caller having no way to find out.
     */
    @Test
    public void anExhaustedBudgetIsReported() {
        DFunction wiggly = x -> Math.sin(500.0 * x);
        AdaptiveGaussKronrod.IntegralResult starved =
                AdaptiveGaussKronrod.integrate1DAdaptiveWithError(ruleSetup, wiggly, 0.0, 1.0, 1.0e-12, 2);
        assertFalse("two levels cannot resolve 500 half waves, was " + starved, starved.converged);

        AdaptiveGaussKronrod.IntegralResult ample =
                AdaptiveGaussKronrod.integrate1DAdaptiveWithError(ruleSetup, wiggly, 0.0, 1.0, 1.0e-12, 24);
        assertTrue("with room to subdivide it converges, was " + ample, ample.converged);
        assertEquals("and gets the right answer", (1.0 - Math.cos(500.0)) / 500.0, ample.value, 1.0e-12);
    }

    /** A single application of the rule has no tolerance to miss. */
    @Test
    public void anUndividedRuleReportsConverged() {
        DFunction f = x -> x * x;
        assertTrue("1D", AdaptiveGaussKronrod.integrate1D(ruleSetup, f, 0.0, 3.0).converged);
        assertTrue("2D", AdaptiveGaussKronrod.integrate2DParallel(ruleSetup, (x, y) -> x * y, 0.0, 1.0, 0.0, 1.0)
                .converged);
        assertTrue("3D", AdaptiveGaussKronrod.integrate3DParallel(ruleSetup, (x, y, z) -> x * y * z, 0.0, 1.0, 0.0, 1.0,
                0.0, 1.0).converged);
    }
}
