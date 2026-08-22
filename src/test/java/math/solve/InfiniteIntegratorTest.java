package math.solve;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

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

    // =========================================================================
    // THE GUARD AGAINST A SUBSTITUTION THAT MISSED THE INTEGRAND
    // =========================================================================

    /** A unit-height Gaussian blob centered at (cx, cy), of total mass one. */
    private static DBiFunction blob2D(final double cx, final double cy, final double sigma) {
        return (x, y) -> {
            double dx = (x - cx) / sigma;
            double dy = (y - cy) / sigma;
            return Math.exp(-0.5 * (dx * dx + dy * dy)) / (2.0 * Math.PI * sigma * sigma);
        };
    }

    /** The same in three dimensions, centered on the diagonal. */
    private static DTriFunction blob3D(final double c, final double sigma) {
        return (x, y, z) -> {
            double dx = (x - c) / sigma;
            double dy = (y - c) / sigma;
            double dz = (z - c) / sigma;
            double norm = Math.pow(2.0 * Math.PI, 1.5) * sigma * sigma * sigma;
            return Math.exp(-0.5 * (dx * dx + dy * dy + dz * dz)) / norm;
        };
    }

    private static double wholePlane(DBiFunction f) {
        return InfiniteIntegrator.integrate2DInfinite(ruleSetup, f, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0e-13);
    }

    private static double wholeSpace(DTriFunction f) {
        return InfiniteIntegrator.integrate3DInfinite(ruleSetup, f, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0e-13);
    }

    /** The point a refusal names, taken back out of its message. */
    private static double[] namedPoint(String message) {
        int open = message.indexOf("|f(");
        int close;
        if (open >= 0) {
            open += 3;
            close = message.indexOf(")|", open);
        } else {
            open = message.indexOf("splitting at (") + "splitting at (".length();
            close = message.indexOf(")", open);
        }
        String[] parts = message.substring(open, close).split(",");
        double[] point = new double[parts.length];
        for (int i = 0; i < parts.length; ++i) {
            point[i] = Double.parseDouble(parts[i].trim());
        }
        return point;
    }

    /**
     * Mass far from the origin arrives as a spike against the edge of the
     * transformed square and the subdivision never resolves it. The result used
     * to come back as a plausible-looking near-zero; it is now refused, and the
     * refusal names the place. Splitting the plane there has to recover the
     * integral, or the advice would be worthless.
     */
    @Test
    public void the2DFormRefusesMassItNeverSampled() {
        double[][] centres = { { 100.0, 100.0 }, { 137.0, 42.0 }, { -200.0, 60.0 }, { 300.0, 40.0 } };
        for (double[] centre : centres) {
            DBiFunction f = blob2D(centre[0], centre[1], centre[1] == 40.0 ? 10.0 : 1.0);
            String where = "blob at (" + centre[0] + ", " + centre[1] + ")";
            try {
                wholePlane(f);
                fail(where + ": the substitution missed the mass and should have refused");
            } catch (ArithmeticException expected) {
                double[] named = namedPoint(expected.getMessage());
                assertEquals(where + ": the refusal must point at the mass in x", centre[0], named[0],
                        0.05 * Math.abs(centre[0]));
                assertEquals(where + ": and in y", centre[1], named[1], 0.05 * Math.abs(centre[1]));

                double split = 0.0;
                for (int quadrant = 0; quadrant < 4; ++quadrant) {
                    double lowX = ((quadrant & 1) == 0) ? named[0] : Double.NEGATIVE_INFINITY;
                    double highX = ((quadrant & 1) == 0) ? Double.POSITIVE_INFINITY : named[0];
                    double lowY = ((quadrant & 2) == 0) ? named[1] : Double.NEGATIVE_INFINITY;
                    double highY = ((quadrant & 2) == 0) ? Double.POSITIVE_INFINITY : named[1];
                    split += InfiniteIntegrator.integrate2DInfinite(ruleSetup, f, lowX, highX, lowY, highY,
                            1.0e-13);
                }
                assertEquals(where + ": splitting where the refusal points must recover the integral", 1.0,
                        split, 1.0e-6);
            }
        }
    }

    /**
     * The harder half. At ten widths out the substitution does sample the mass,
     * so no probe can tell that anything went wrong -- and the answer is still
     * four parts in a thousand off. Splitting at the located peak is what
     * exposes it.
     */
    @Test
    public void the2DFormRefusesWhatItSampledButCouldNotResolve() {
        try {
            wholePlane(blob2D(10.0, 10.0, 1.0));
            fail("the substitution resolved the blob badly and should have refused");
        } catch (ArithmeticException expected) {
            String message = expected.getMessage();
            assertTrue("the refusal must distinguish this from a missed spike: " + message,
                    message.contains("sampled the integrand but did not resolve it"));
            double[] named = namedPoint(message);
            assertEquals("it must still name the mass", 10.0, named[0], 0.5);
            assertEquals("it must still name the mass", 10.0, named[1], 0.5);
        }
    }

    /**
     * The half that decides whether the guard is worth having: it must not cry
     * wolf. An integral that cancels to zero, one that is zero everywhere, one
     * that is merely very small, and heavy tails all have to pass untouched.
     */
    @Test
    public void the2DGuardDoesNotCryWolf() {
        assertEquals("a blob at the origin", 1.0, wholePlane(blob2D(0.0, 0.0, 1.0)), 1.0e-6);
        assertEquals("a blob three widths out", 1.0, wholePlane(blob2D(3.0, 3.0, 1.0)), 1.0e-6);
        assertEquals("an odd integrand cancels", 0.0, wholePlane((x, y) -> x * blob2D(0.0, 0.0, 1.0).apply(x, y)),
                1.0e-9);
        assertEquals("the zero function", 0.0, wholePlane((x, y) -> 0.0), 0.0);
        assertEquals("a genuinely tiny integral", 1.0e-200,
                wholePlane((x, y) -> 1.0e-200 * blob2D(0.0, 0.0, 1.0).apply(x, y)), 1.0e-205);
        assertEquals("heavy tails on both axes", Math.PI * Math.PI,
                wholePlane((x, y) -> 1.0 / ((1.0 + x * x) * (1.0 + y * y))), 1.0e-6);
    }

    /** The three dimensional form carries the same guard, and needs it sooner. */
    @Test
    public void the3DFormRefusesMassItCannotResolve() {
        double[][] cases = { { 30.0, 1.0 }, { 100.0, 10.0 } };
        for (double[] c : cases) {
            try {
                wholeSpace(blob3D(c[0], c[1]));
                fail("blob3D(" + c[0] + ", " + c[1] + "): should have refused");
            } catch (ArithmeticException expected) {
                double[] named = namedPoint(expected.getMessage());
                assertEquals("the refusal must point at the mass", c[0], named[0], 0.05 * c[0]);
                assertEquals("on every axis", c[0], named[2], 0.05 * c[0]);
            }
        }
    }

    /** And it must leave the well behaved three dimensional cases alone. */
    @Test
    public void the3DGuardDoesNotCryWolf() {
        assertEquals("a blob at the origin", 1.0, wholeSpace(blob3D(0.0, 1.0)), 1.0e-5);
        assertEquals("a blob three widths out", 1.0, wholeSpace(blob3D(3.0, 1.0)), 1.0e-5);
        assertEquals("the zero function", 0.0, wholeSpace((x, y, z) -> 0.0), 0.0);
        assertEquals("an odd integrand cancels", 0.0,
                wholeSpace((x, y, z) -> x * blob3D(0.0, 1.0).apply(x, y, z)), 1.0e-9);
    }
}
