package math.solve;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
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
}
