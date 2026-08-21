package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.RootFinder;

/**
 * The evaluation machinery both interpolators share: value, derivatives, the
 * closed-form integral and the argument checks.
 */
public class CubicSplineTest {

    private static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }
    }

    /** A data set whose knots are irregularly spaced and whose values wander. */
    private static double[][] wobble(int n, long seed) {
        Lcg rnd = new Lcg(seed);
        double[] x = new double[n];
        double[] y = new double[n];
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] = acc;
            y[i] = 5.0 * rnd.next() - 2.0;
            acc += 0.2 + rnd.next();
        }
        return new double[][] { x, y };
    }

    private static double[][] smooth(int n) {
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i * (10.0 / (n - 1));
            y[i] = Math.sin(x[i]) + 0.3 * x[i];
        }
        return new double[][] { x, y };
    }

    /** Every builder over the same data, so the shared machinery is exercised by each. */
    private static CubicSpline[] both(double[][] data) {
        return new CubicSpline[] { SplineInterpolator.interpolate(data[0], data[1]),
                KrugerInterpolator.interpolate(data[0], data[1]),
                AkimaInterpolator.interpolate(data[0], data[1]),
                AkimaInterpolator.interpolate(data[0], data[1], AkimaInterpolator.Variant.MODIFIED) };
    }

    private static final String[] NAMES = { "natural spline", "Kruger spline", "Akima spline",
            "modified Akima spline" };

    private static String name(int i) {
        return NAMES[i];
    }

    // the piece that owns a segment, re-derived from the published coefficients

    private static double valueOfPiece(CubicSpline s, int piece, double x) {
        double[] c = s.coefficients(piece);
        double t = x - s.knot(piece);
        return t * (t * (t * c[3] + c[2]) + c[1]) + c[0];
    }

    private static double slopeOfPiece(CubicSpline s, int piece, double x) {
        double[] c = s.coefficients(piece);
        double t = x - s.knot(piece);
        return t * (t * 3.0 * c[3] + 2.0 * c[2]) + c[1];
    }

    @Test
    public void testTheSplinePassesThroughEveryDataPoint() {
        double[][][] sets = { smooth(11), smooth(2), wobble(12, 7L), wobble(41, 99L) };
        for (double[][] data : sets) {
            CubicSpline[] splines = both(data);
            for (int s = 0; s < splines.length; s++) {
                for (int i = 0; i < data[0].length; i++) {
                    assertEquals(name(s) + ", knot " + i, data[1][i], splines[s].value(data[0][i]), 1.0e-14);
                }
            }
        }
    }

    @Test
    public void testTheLastKnotIsEvaluatedByTheLastSegment() {
        double[][] data = wobble(9, 12345L);
        CubicSpline[] splines = both(data);
        int last = data[0].length - 1;
        for (int s = 0; s < splines.length; s++) {
            assertEquals(name(s), data[0][last], splines[s].upperBound(), 0.0);
            assertEquals(name(s), data[0][0], splines[s].lowerBound(), 0.0);
            assertEquals(name(s), data[1][last], splines[s].value(splines[s].upperBound()), 1.0e-14);
            assertEquals(name(s), data[0].length - 1, splines[s].segments());
        }
    }

    @Test
    public void testThePublishedCoefficientsDescribeTheSegmentTheyBelongTo() {
        double[][] data = wobble(10, 314159L);
        CubicSpline[] splines = both(data);
        for (int s = 0; s < splines.length; s++) {
            CubicSpline sp = splines[s];
            for (int piece = 0; piece < sp.segments(); piece++) {
                double lo = sp.knot(piece);
                double hi = sp.knot(piece + 1);
                for (int q = 0; q <= 16; q++) {
                    double x = lo + (hi - lo) * q / 16.0;
                    assertEquals(name(s) + ", piece " + piece, valueOfPiece(sp, piece, x), sp.value(x), 1.0e-14);
                    assertEquals(name(s) + ", piece " + piece, slopeOfPiece(sp, piece, x), sp.derivativeAt(x),
                            1.0e-14);
                }
            }
        }
    }

    @Test
    public void testTheClosedFormIntegralAgreesWithGaussKronrodOnEverySegment() {
        // a fifteen point Kronrod rule is exact for a cubic on one interval, so
        // segment by segment this is an independent evaluation to rounding
        double[][][] sets = { smooth(11), smooth(31), wobble(12, 5L), wobble(41, 6L) };
        for (double[][] data : sets) {
            CubicSpline[] splines = both(data);
            for (int s = 0; s < splines.length; s++) {
                final CubicSpline sp = splines[s];
                DFunction f = new DFunction() {
                    @Override
                    public double apply(double x) {
                        return sp.value(x);
                    }
                };
                double kronrod = 0.0;
                double area = 0.0;
                for (int i = 0; i < sp.segments(); i++) {
                    kronrod += AdaptiveGaussKronrod.integrate1D(G7_K15.POINTS_15, f, sp.knot(i), sp.knot(i + 1)).value;
                    area += Math.abs(sp.integrate(sp.knot(i), sp.knot(i + 1)));
                }
                double closed = sp.integrate(sp.lowerBound(), sp.upperBound());
                assertEquals(name(s), 0.0, Math.abs(closed - kronrod) / area, 1.0e-13);
            }
        }
    }

    @Test
    public void testTheIntegralIsAdditiveAndAntisymmetric() {
        double[][] data = wobble(14, 271828L);
        CubicSpline[] splines = both(data);
        for (int s = 0; s < splines.length; s++) {
            CubicSpline sp = splines[s];
            double a = sp.lowerBound();
            double b = sp.upperBound();
            double area = 0.0;
            for (int i = 0; i < sp.segments(); i++) {
                area += Math.abs(sp.integrate(sp.knot(i), sp.knot(i + 1)));
            }
            double whole = sp.integrate(a, b);
            assertEquals(name(s), 0.0, whole + sp.integrate(b, a), 1.0e-14 * area);
            assertEquals(name(s), 0.0, sp.integrate(a, a), 0.0);
            assertEquals(name(s), 0.0, sp.integrate(b, b), 0.0);
            // split at every knot and at points inside the segments
            for (int i = 0; i <= 40; i++) {
                double mid = a + (b - a) * i / 40.0;
                double split = sp.integrate(a, mid) + sp.integrate(mid, b);
                assertEquals(name(s) + ", split at " + mid, whole, split, 1.0e-14 * area);
            }
        }
    }

    @Test
    public void testTheDerivativeAgreesWithACentralDifference() {
        double[][][] sets = { smooth(11), smooth(31), wobble(12, 8L) };
        for (double[][] data : sets) {
            CubicSpline[] splines = both(data);
            for (int s = 0; s < splines.length; s++) {
                CubicSpline sp = splines[s];
                double worst = 0.0;
                double scale = 0.0;
                for (int piece = 0; piece < sp.segments(); piece++) {
                    double lo = sp.knot(piece);
                    double hi = sp.knot(piece + 1);
                    double h = (hi - lo) * 1.0e-4;
                    for (int q = 1; q <= 7; q++) {
                        double x = lo + (hi - lo) * q / 8.0;
                        double difference = (sp.value(x + h) - sp.value(x - h)) / (2.0 * h);
                        scale = Math.max(scale, Math.abs(difference));
                        worst = Math.max(worst, Math.abs(difference - sp.derivativeAt(x)));
                    }
                }
                assertEquals(name(s), 0.0, worst / scale, 1.0e-6);
            }
        }
    }

    @Test
    public void testTheSecondDerivativeIsTheDerivativeOfTheFirst() {
        double[][][] sets = { smooth(21), wobble(12, 4242L) };
        for (double[][] data : sets) {
            CubicSpline[] splines = both(data);
            for (int s = 0; s < splines.length; s++) {
                CubicSpline sp = splines[s];
                double worst = 0.0;
                double scale = 0.0;
                for (int piece = 0; piece < sp.segments(); piece++) {
                    double lo = sp.knot(piece);
                    double hi = sp.knot(piece + 1);
                    double h = (hi - lo) * 1.0e-4;
                    for (int q = 1; q <= 7; q++) {
                        double x = lo + (hi - lo) * q / 8.0;
                        double difference = (sp.derivativeAt(x + h) - sp.derivativeAt(x - h)) / (2.0 * h);
                        scale = Math.max(scale, Math.abs(difference));
                        worst = Math.max(worst, Math.abs(difference - sp.secondDerivativeAt(x)));
                    }
                }
                assertEquals(name(s), 0.0, worst / scale, 1.0e-9);
            }
        }
    }

    @Test
    public void testTheStationaryPointOfTheDerivativeIsWhereTheValueIsLargest() {
        // one clear interior maximum; brentDekker locates the root of the
        // derivative by a method that knows nothing about splines
        double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] y = { 0.0, 1.5, 3.0, 4.2, 3.1, 1.2, 0.1 };
        CubicSpline[] splines = both(new double[][] { x, y });
        for (int s = 0; s < splines.length; s++) {
            final CubicSpline sp = splines[s];
            DFunction derivative = new DFunction() {
                @Override
                public double apply(double t) {
                    return sp.derivativeAt(t);
                }
            };
            assertTrue(name(s), sp.derivativeAt(2.5) > 0.0);
            assertTrue(name(s), sp.derivativeAt(4.0) < 0.0);
            double root = RootFinder.brentDekker(2.5, 4.0, derivative, 1.0e-12);
            assertEquals(name(s), 0.0, sp.derivativeAt(root), 1.0e-9);
            double atRoot = sp.value(root);
            for (int i = 0; i <= 6000; i++) {
                double t = 6.0 * i / 6000.0;
                assertTrue(name(s) + ", larger value at " + t, sp.value(t) <= atRoot + 1.0e-12);
            }
        }
    }

    @Test
    public void testAnArgumentOutsideTheKnotsIsRejected() {
        double[][] data = smooth(7);
        CubicSpline[] splines = both(data);
        for (int s = 0; s < splines.length; s++) {
            CubicSpline sp = splines[s];
            double[] outside = { -1.0e-9, 10.0 + 1.0e-9, -5.0, 100.0 };
            for (int i = 0; i < outside.length; i++) {
                rejects(name(s) + ", value at " + outside[i], sp, 0, outside[i]);
                rejects(name(s) + ", derivative at " + outside[i], sp, 1, outside[i]);
                rejects(name(s) + ", second derivative at " + outside[i], sp, 2, outside[i]);
                rejects(name(s) + ", integrate to " + outside[i], sp, 3, outside[i]);
            }
        }
    }

    @Test
    public void testANonFiniteArgumentIsRejected() {
        // the range check alone lets a NaN through: it fails both comparisons
        double[][] data = smooth(7);
        CubicSpline[] splines = both(data);
        double[] bad = { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
        for (int s = 0; s < splines.length; s++) {
            for (int i = 0; i < bad.length; i++) {
                rejects(name(s) + ", value at " + bad[i], splines[s], 0, bad[i]);
                rejects(name(s) + ", derivative at " + bad[i], splines[s], 1, bad[i]);
                rejects(name(s) + ", second derivative at " + bad[i], splines[s], 2, bad[i]);
                rejects(name(s) + ", integrate to " + bad[i], splines[s], 3, bad[i]);
            }
        }
    }

    private static void rejects(String what, CubicSpline sp, int which, double x) {
        try {
            switch (which) {
            case 0:
                sp.value(x);
                break;
            case 1:
                sp.derivativeAt(x);
                break;
            case 2:
                sp.secondDerivativeAt(x);
                break;
            default:
                sp.integrate(sp.lowerBound(), x);
                break;
            }
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            // that is the point
        }
    }

    @Test
    public void testTheAccessorsHandOutCopies() {
        double[] x = { 0.0, 1.0, 2.0, 3.0 };
        double[] y = { 0.0, 1.0, 0.0, 1.0 };
        CubicSpline sp = SplineInterpolator.interpolate(x, y);
        double before = sp.value(1.5);
        // the arrays the caller passed in
        x[1] = 2.9;
        y[1] = 99.0;
        assertEquals("mutating the input changed the spline", before, sp.value(1.5), 0.0);
        // and the arrays the accessors return
        double[] knots = sp.knots();
        assertNotSame(knots, sp.knots());
        knots[0] = -5.0;
        assertEquals(0.0, sp.knot(0), 0.0);
        double[] coefficients = sp.coefficients(0);
        assertNotSame(coefficients, sp.coefficients(0));
        coefficients[0] = -5.0;
        assertEquals(0.0, sp.coefficients(0)[0], 0.0);
        assertEquals(before, sp.value(1.5), 0.0);
    }

    @Test
    public void testTheKnotAndSegmentIndicesAreChecked() {
        CubicSpline sp = SplineInterpolator.interpolate(new double[] { 0.0, 1.0, 2.0 },
                new double[] { 1.0, 2.0, 0.0 });
        assertEquals(2, sp.segments());
        assertEquals(3, sp.knots().length);
        int[] badKnots = { -1, 3, 4 };
        for (int i = 0; i < badKnots.length; i++) {
            try {
                sp.knot(badKnots[i]);
                fail("knot(" + badKnots[i] + ") was accepted");
            } catch (IllegalArgumentException expected) {
                // as intended
            }
        }
        int[] badSegments = { -1, 2, 3 };
        for (int i = 0; i < badSegments.length; i++) {
            try {
                sp.coefficients(badSegments[i]);
                fail("coefficients(" + badSegments[i] + ") was accepted");
            } catch (IllegalArgumentException expected) {
                // as intended
            }
        }
        assertEquals(4, sp.coefficients(0).length);
        assertEquals(4, sp.coefficients(1).length);
    }

    @Test
    public void testTheFunctionInterfaceIsTheSameFunction() {
        double[][] data = wobble(11, 2024L);
        CubicSpline[] splines = both(data);
        for (int s = 0; s < splines.length; s++) {
            CubicSpline sp = splines[s];
            for (int i = 0; i <= 500; i++) {
                double x = sp.lowerBound() + (sp.upperBound() - sp.lowerBound()) * i / 500.0;
                assertEquals(name(s), Double.doubleToRawLongBits(sp.value(x)),
                        Double.doubleToRawLongBits(sp.apply(x)));
            }
        }
    }
}
