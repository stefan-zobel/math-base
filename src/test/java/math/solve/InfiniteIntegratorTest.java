package math.solve;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import math.fun.DFunction;
import math.fun.DBiFunction;
import math.fun.DTriFunction;

/**
 * Unit tests for the InfiniteIntegrator class. Validates 1D, 2D, and 3D
 * coordinate transformations for improper integrals with infinite bounds.
 */
public class InfiniteIntegratorTest {

    private static AdaptiveGaussKronrod.G7_K15 ruleSetup;

    // Strict precision for well-behaved exponential decays
    private static final double INTEGRAL_TOLERANCE = 1e-6;
    // Slightly looser tolerance for highly compressed multi-dimensional infinite spaces
    private static final double SPATIAL_TOLERANCE = 1e-4;

    @BeforeClass
    public static void setUpBeforeClass() {
        ruleSetup = AdaptiveGaussKronrod.G7_K15.POINTS_15;
    }

    // =========================================================================
    // 1D INFINITE BOUNDS TESTING
    // =========================================================================

    @Test
    public void test1DFullyInfiniteBounds() {
        // Standard Gaussian curve: f(x) = e^(-x^2) over [-inf, +inf]
        // Analytical solution = sqrt(pi)
        DFunction standardGaussian = x -> Math.exp(-x * x);
        double exactValue = Math.sqrt(Math.PI);

        double result = InfiniteIntegrator.integrate1DInfinite(
            ruleSetup, standardGaussian, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, INTEGRAL_TOLERANCE
        );

        assertEquals("1D fully infinite bounds [-inf, +inf] transformation failed", exactValue, result, INTEGRAL_TOLERANCE);
    }

    @Test
    public void test1DSemiInfiniteUpperBounds() {
        // Gamma-like distribution peak: f(x) = x * e^(-x) over [0, +inf]
        // Analytical solution = 1.0
        DFunction upperInfinite = x -> x * Math.exp(-x);
        double exactValue = 1.0;

        double result = InfiniteIntegrator.integrate1DInfinite(
            ruleSetup, upperInfinite, 0.0, Double.POSITIVE_INFINITY, INTEGRAL_TOLERANCE
        );

        assertEquals("1D semi-infinite upper bounds [a, +inf] transformation failed", exactValue, result, INTEGRAL_TOLERANCE);
    }

    @Test
    public void test1DSemiInfiniteLowerBounds() {
        // Mirrored exponential decay: f(x) = e^x over [-inf, 0]
        // Analytical solution = 1.0
        DFunction lowerInfinite = x -> Math.exp(x);
        double exactValue = 1.0;

        double result = InfiniteIntegrator.integrate1DInfinite(
            ruleSetup, lowerInfinite, Double.NEGATIVE_INFINITY, 0.0, INTEGRAL_TOLERANCE
        );

        assertEquals("1D semi-infinite lower bounds [-inf, b] transformation failed", exactValue, result, INTEGRAL_TOLERANCE);
    }

    @Test
    public void test1DFiniteFallbackBounds() {
        // Plain fallback test: f(x) = x^2 over [0, 3]
        // Analytical solution = 3^3 / 3 = 9.0
        DFunction finitePolynomial = x -> x * x;
        double exactValue = 9.0;

        double result = InfiniteIntegrator.integrate1DInfinite(
            ruleSetup, finitePolynomial, 0.0, 3.0, INTEGRAL_TOLERANCE
        );

        assertEquals("1D finite fallback route [a, b] inside InfiniteIntegrator failed", exactValue, result, INTEGRAL_TOLERANCE);
    }

    // =========================================================================
    // 2D INFINITE BOUNDS TESTING
    // =========================================================================

    @Test
    public void test2DMixedInfiniteBounds() {
        // Mixed domain test: f(x, y) = e^(-x) * e^(-y^2) over [0, +inf] x [-inf, +inf]
        // Analytical solution = (Integral of e^-x from 0 to inf) * (Integral of e^-y^2 from -inf to inf)
        //                     = 1.0 * sqrt(pi) = sqrt(pi)
        DBiFunction mixed2D = (x, y) -> Math.exp(-x) * Math.exp(-y * y);
        double exactValue = Math.sqrt(Math.PI);

        double result = InfiniteIntegrator.integrate2DInfinite(
            ruleSetup, mixed2D, 0.0, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, INTEGRAL_TOLERANCE
        );

        assertEquals("2D mixed finite/infinite boundary matrix evaluation failed", exactValue, result, INTEGRAL_TOLERANCE);
    }

    @Test
    public void test2DFullyInfiniteBounds() {
        // f(x, y) = e^(-x^2 - y^2) over [-inf, +inf] x [-inf, +inf]
        // Analytical solution = (sqrt(pi))^2 = pi
        // Tests CASE A (doubly infinite) on both axes simultaneously.
        DBiFunction gaussian2D = (x, y) -> Math.exp(-x * x - y * y);
        double exactValue = Math.PI;

        double result = InfiniteIntegrator.integrate2DInfinite(
            ruleSetup, gaussian2D,
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
            INTEGRAL_TOLERANCE
        );

        assertEquals("2D fully infinite bounds [-inf,+inf]^2 (CASE A x CASE A) failed",
                     exactValue, result, SPATIAL_TOLERANCE);
    }

    @Test
    public void test2DLowerInfiniteBoundOnXAxis() {
        // f(x, y) = e^x * e^(-y^2) over [-inf, 0] x [-inf, +inf]
        // Analytical solution = integral(e^x, -inf, 0) * integral(e^(-y^2), -inf, +inf)
        //                     = 1.0 * sqrt(pi) = sqrt(pi)
        // Tests CASE C (semi-infinite downward) on x-axis and CASE A on y-axis.
        DBiFunction lowerMixed2D = (x, y) -> Math.exp(x) * Math.exp(-y * y);
        double exactValue = Math.sqrt(Math.PI);

        double result = InfiniteIntegrator.integrate2DInfinite(
            ruleSetup, lowerMixed2D,
            Double.NEGATIVE_INFINITY, 0.0,
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
            INTEGRAL_TOLERANCE
        );

        assertEquals("2D lower-infinite x-axis bounds [-inf,0] x [-inf,+inf] (CASE C x CASE A) failed",
                     exactValue, result, INTEGRAL_TOLERANCE);
    }

    @Test
    public void test2DFiniteFallbackBounds() {
        // f(x, y) = x + y over [0, 1] x [0, 1]
        // Analytical solution = integral(x+y, 0,1,0,1) = 1/2 + 1/2 = 1.0
        // Tests CASE D (finite fallback) on both axes via the 2D infinite entry point.
        DBiFunction linearSum2D = (x, y) -> x + y;
        double exactValue = 1.0;

        double result = InfiniteIntegrator.integrate2DInfinite(
            ruleSetup, linearSum2D,
            0.0, 1.0,
            0.0, 1.0,
            INTEGRAL_TOLERANCE
        );

        assertEquals("2D finite fallback bounds [0,1]^2 (CASE D x CASE D) inside InfiniteIntegrator failed",
                     exactValue, result, INTEGRAL_TOLERANCE);
    }

    // =========================================================================
    // 3D INFINITE BOUNDS TESTING
    // =========================================================================

    @Test
    public void test3DFullyInfiniteBounds() {
        // 3D Spatial Gaussian bell curve: f(x, y, z) = e^(-(x^2 + y^2 + z^2)) over [-inf, +inf]^3
        // Analytical solution = (sqrt(pi))^3 = pi^1.5 ~= 5.56832799
        DTriFunction spatialGaussian3D = (x, y, z) -> Math.exp(-(x * x + y * y + z * z));
        double exactValue = Math.pow(Math.PI, 1.5);

        double result = InfiniteIntegrator.integrate3DInfinite(
            ruleSetup,
            spatialGaussian3D, 
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 
            1e-6
        );

        // Uses spatial tolerance to account for the heavy compression of infinite R^3 domain space
        assertEquals("3D spatial infinite volume integration tensor mapping failed", exactValue, result, SPATIAL_TOLERANCE);
    }

    @Test
    public void test3DMixedBoundsSomeInfiniteAxes() {
        // f(x, y, z) = e^(-x) * e^(-y^2) * e^(-z) over [0, +inf] x [-inf, +inf] x [0, +inf]
        // Analytical solution = integral(e^-x, 0,+inf) * integral(e^-y^2,-inf,+inf) * integral(e^-z, 0,+inf)
        //                     = 1.0 * sqrt(pi) * 1.0 = sqrt(pi)
        // Tests CASE B on x-axis, CASE A on y-axis, CASE B on z-axis — three different
        // transformations active simultaneously in one call.
        DTriFunction separableMixed3D = (x, y, z) -> Math.exp(-x) * Math.exp(-y * y) * Math.exp(-z);
        double exactValue = Math.sqrt(Math.PI);

        double result = InfiniteIntegrator.integrate3DInfinite(
            ruleSetup,
            separableMixed3D,
            0.0,                       Double.POSITIVE_INFINITY,  // X: [0, +inf]  CASE B
            Double.NEGATIVE_INFINITY,  Double.POSITIVE_INFINITY,  // Y: [-inf,+inf] CASE A
            0.0,                       Double.POSITIVE_INFINITY,  // Z: [0, +inf]  CASE B
            INTEGRAL_TOLERANCE
        );

        assertEquals("3D mixed infinite/finite bounds (CASE B x CASE A x CASE B) failed",
                     exactValue, result, SPATIAL_TOLERANCE);
    }
}
