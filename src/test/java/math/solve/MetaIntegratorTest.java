package math.solve;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import math.fun.DFunction;
import math.fun.DBiFunction;
import math.fun.DTriFunction;

/**
 * Unit tests for the MetaIntegrator class covering smooth, oscillatory, 
 * and highly localized peak behaviors across 1D, 2D, and 3D domains.
 */
public class MetaIntegratorTest {

    private static AdaptiveGaussKronrod.G7_K15 ruleSetup;
    private static final double TOLERANCE = 1e-6;

    @BeforeClass
    public static void setUpBeforeClass() {
        ruleSetup = AdaptiveGaussKronrod.G7_K15.POINTS_15;
    }

    // =========================================================================
    // 1D INTEGRATION TESTS
    // =========================================================================

    @Test
    public void test1DSmoothFunction() {
        // f(x) = x^3 over [0, 2]. Exact analytical value = 4.0
        DFunction smooth1D = x -> x * x * x;
        double exactValue = 4.0;

        double result = MetaIntegrator.integrate1DSmart(ruleSetup, smooth1D, 0.0, 2.0, TOLERANCE);

        assertEquals("Smooth 1D polynomial integration failed", exactValue, result, TOLERANCE);
    }

    @Test
    public void test1DOscillatoryFunction() {
        // f(x) = sin(30*x) over [0, 1]. Exact analytical value = (1 - cos(30)) / 30
        DFunction oscillatory1D = x -> Math.sin(30.0 * x);
        double exactValue = (1.0 - Math.cos(30.0)) / 30.0;

        double result = MetaIntegrator.integrate1DSmart(ruleSetup, oscillatory1D, 0.0, 1.0, TOLERANCE);

        assertEquals("Highly oscillatory 1D integration failed", exactValue, result, TOLERANCE);
    }

    @Test
    public void test1DSubdivisionPathOnNarrowPeak() {
        // Narrow positive Gaussian spike at x=0.4: oscillation index ? 1.0 (no cancellation),
        // so the MetaIntegrator does NOT route to Clenshaw-Curtis.
        // The first G7-K15 step overestimates and has large error ? triggers Case 3
        // (integrate1DAdaptive). This is the only path not covered by the smooth/oscillatory tests.
        DFunction spike1D = x -> Math.exp(-500.0 * Math.pow(x - 0.4, 2));
        // True integral ? sqrt(PI/500); boundary mass is negligible (spike half-width ? 0.032)
        double exactValue = Math.sqrt(Math.PI / 500.0);

        double result = MetaIntegrator.integrate1DSmart(ruleSetup, spike1D, 0.0, 1.0, TOLERANCE);

        assertEquals("MetaIntegrator Case-3 (adaptive subdivision) path failed on narrow Gaussian spike",
                     exactValue, result, 1e-4);
    }

    @Test
    public void test1DSmallValueFunctionNotPrematurelyExited() {
        // f(x) = sin(30x) * exp(-200x) over [0, 1]: a damped oscillation whose net integral is
        // small (~7.3e-4) but genuine and must be computed accurately.
        // Regression for the firstStep.value < epsTol bug in integrate1DSmart: the old code
        // returned immediately whenever the GK approximation was below epsTol, bypassing
        // the error-estimate check. If the GK estimate happened to land near zero, the
        // returned value would be wrong. The fixed code checks error estimate, not value.
        // Analytical value: Im[integral of e^((-200+30i)x) from 0 to 1]
        //   = [30*(1 - e^(-200)*cos(30)) - 200*e^(-200)*sin(30)] / (200^2 + 30^2)
        //   ? 30 / 40900 ? 7.336e-4  (since e^(-200) ? 0)
        DFunction dampedOsc = x -> Math.sin(30.0 * x) * Math.exp(-200.0 * x);
        double exactValue = (30.0 * (1.0 - Math.exp(-200.0) * Math.cos(30.0))
                            - 200.0 * Math.exp(-200.0) * Math.sin(30.0))
                            / (200.0 * 200.0 + 30.0 * 30.0);

        double result = MetaIntegrator.integrate1DSmart(ruleSetup, dampedOsc, 0.0, 1.0, TOLERANCE);

        assertEquals("1D damped-oscillation integration failed (small-value non-premature-exit regression)",
                     exactValue, result, TOLERANCE);
    }

    // =========================================================================
    // 2D INTEGRATION TESTS
    // =========================================================================

    @Test
    public void test2DOscillatoryFunction() {
        // f(x, y) = sin(20*x) * cos(20*y) over [0, 1]x[0, 1]
        // Exact value = ((1 - cos(20)) / 20) * (sin(20) / 20)
        DBiFunction oscillatory2D = (x, y) -> Math.sin(20.0 * x) * Math.cos(20.0 * y);
        double exactValue = ((1.0 - Math.cos(20.0)) / 20.0) * (Math.sin(20.0) / 20.0);

        // This must route to performClenshawCurtis2D automatically
        double result = MetaIntegrator.integrate2DSmart(ruleSetup, oscillatory2D, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("Highly oscillatory 2D matrix calculation failed inside Clenshaw-Curtis router", 
                     exactValue, result, TOLERANCE);
    }

    @Test
    public void test2DSmoothFunction() {
        // f(x, y) = x^2 + y^2 over [0, 2] x [0, 2]. Exact value = 32/3 = 10.666666...
        DBiFunction smooth2D = (x, y) -> x * x + y * y;
        double exactValue = 32.0 / 3.0;

        double result = MetaIntegrator.integrate2DSmart(ruleSetup, smooth2D, 0.0, 2.0, 0.0, 2.0, TOLERANCE);

        assertEquals("Smooth 2D quadratic field integration failed", exactValue, result, TOLERANCE);
    }

    @Test
    public void test2DLocalizedPeakFunction() {
        // Steep Gaussian distribution focused tightly around center (0.5, 0.5)
        DBiFunction peak2D = (x, y) -> Math.exp(-150.0 * (Math.pow(x - 0.5, 2) + Math.pow(y - 0.5, 2)));
        // Approximate infinite plane integral baseline over full domain is PI / 150
        double referenceValue = Math.PI / 150.0;

        double result = MetaIntegrator.integrate2DSmart(ruleSetup, peak2D, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("Localized 2D Gaussian peak integration failed", referenceValue, result, 1e-4);
    }

    // =========================================================================
    // 3D INTEGRATION TESTS
    // =========================================================================

    @Test
    public void test3DSmoothFunction() {
        // f(x, y, z) = x * y * z over [0, 1] x [0, 1] x [0, 1]. Exact value = 1/8 = 0.125
        DTriFunction smooth3D = (x, y, z) -> x * y * z;
        double exactValue = 0.125;

        double result = MetaIntegrator.integrate3DSmart(ruleSetup, smooth3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("Smooth 3D multiplication linear field failed", exactValue, result, TOLERANCE);
    }

    @Test
    public void test3DNonSeparableFunction() {
        // f(x, y, z) = e^(x+y) * (z^2 + x*y) over [0, 1]^3.
        DTriFunction nonSeparable3D = (x, y, z) -> Math.exp(x + y) * (z * z + x * y);
        double exactValue = (1.0 / 3.0 * (Math.E * Math.E)) - (2.0 / 3.0 * Math.E) + (4.0 / 3.0);

        double result = MetaIntegrator.integrate3DSmart(ruleSetup, nonSeparable3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("Non-separable 3D mathematical field failed", exactValue, result, TOLERANCE);
    }

    // =========================================================================
    // MULTI-DIMENSIONAL OSCILLATORY TESTING
    // =========================================================================

    @Test
    public void test3DOscillatoryFunction() {
        // f(x, y, z) = sin(10*x) * cos(10*y) * sin(10*z) over [0, 1]^3
        // Exact value = ((1 - cos(10)) / 10) * (sin(10) / 10) * ((1 - cos(10)) / 10)
        DTriFunction oscillatory3D = (x, y, z) -> Math.sin(10.0 * x) * Math.cos(10.0 * y) * Math.sin(10.0 * z);
        double exactValue = ((1.0 - Math.cos(10.0)) / 10.0) * (Math.sin(10.0) / 10.0) * ((1.0 - Math.cos(10.0)) / 10.0);

        double result = MetaIntegrator.integrate3DSmart(ruleSetup, oscillatory3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("Highly oscillatory 3D tensor calculation failed",  exactValue, result, TOLERANCE);
    }

    @Test
    public void test3DHighFrequencyOscillatoryFunction() {
        // f(x, y, z) = sin(35*x) * cos(35*y) * sin(35*z) over [0, 1]^3
        // Exact analytical value = ((1 - cos(35)) / 35) * (sin(35) / 35) * ((1 - cos(35)) / 35)
        double freq = 35.0;
        DTriFunction highFreq3D = (x, y, z) -> Math.sin(freq * x) * Math.cos(freq * y) * Math.sin(freq * z);

        double exactValue = ((1.0 - Math.cos(freq)) / freq) * (Math.sin(freq) / freq) * ((1.0 - Math.cos(freq)) / freq);

        // This will force K15 to fail due to heavy under-sampling,
        // triggering the console log and routing to performClenshawCurtis3D
        double result = MetaIntegrator.integrate3DSmart(ruleSetup, highFreq3D, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOLERANCE);

        assertEquals("High frequency 3D Clenshaw-Curtis routing failed", exactValue, result, 1e-4);
    }
}
