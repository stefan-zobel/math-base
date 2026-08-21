package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fit.AkimaInterpolator.Variant;

/**
 * What Akima brings that neither of the other two schemes has: a slope rule
 * that looks two secants out, so a single bad value disturbs only its own
 * neighbourhood. And what the two weightings cost each other.
 */
public class AkimaInterpolatorTest {

    private static final Variant[] VARIANTS = { Variant.CLASSIC, Variant.MODIFIED };

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

    private static double[][] smooth(int n, double shift) {
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double t = i * (10.0 / (n - 1));
            x[i] = t + shift;
            y[i] = Math.sin(t) + 0.3 * t;
        }
        return new double[][] { x, y };
    }

    private static double[][] monotone(int n, long seed, boolean uneven) {
        Lcg rnd = new Lcg(seed);
        double[] x = new double[n];
        double[] y = new double[n];
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] = i + (uneven ? 0.4 * (i % 3) : 0.0);
            acc += Math.pow(rnd.next(), 3.0);
            y[i] = acc;
        }
        return new double[][] { x, y };
    }

    /** The one-sided limits at a knot, from the published coefficients. */
    private static double limit(CubicSpline s, int piece, double x, int order) {
        double[] c = s.coefficients(piece);
        double t = x - s.knot(piece);
        if (order == 0) {
            return t * (t * (t * c[3] + c[2]) + c[1]) + c[0];
        }
        if (order == 1) {
            return t * (t * 3.0 * c[3] + 2.0 * c[2]) + c[1];
        }
        return 6.0 * c[3] * t + 2.0 * c[2];
    }

    private static double magnitude(CubicSpline s, int order) {
        double m = 0.0;
        for (int i = 0; i < s.segments(); i++) {
            for (int q = 0; q <= 8; q++) {
                double x = s.knot(i) + (s.knot(i + 1) - s.knot(i)) * q / 8.0;
                m = Math.max(m, Math.abs(limit(s, i, x, order)));
            }
        }
        return m > 0.0 ? m : 1.0;
    }

    private static boolean sameSegment(CubicSpline a, CubicSpline b, int segment) {
        double[] ca = a.coefficients(segment);
        double[] cb = b.coefficients(segment);
        for (int k = 0; k < 4; k++) {
            if (Double.doubleToRawLongBits(ca[k]) != Double.doubleToRawLongBits(cb[k])) {
                return false;
            }
        }
        return true;
    }

    @Test
    public void testMovingOneValueChangesOnlyItsOwnNeighbourhood() {
        // the slope at a knot reads the two secants either side of it, so one
        // value reaches five knots and therefore six segments -- and not one
        // more, exactly, which is the property the natural spline lacks
        int n = 21;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = Math.sin(0.4 * i) + 0.1 * i;
        }
        int moved = 10;
        double[] disturbed = y.clone();
        disturbed[moved] += 5.0;

        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline before = AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
            CubicSpline after = AkimaInterpolator.interpolate(x, disturbed, VARIANTS[v]);
            int changed = 0;
            for (int seg = 0; seg < n - 1; seg++) {
                boolean same = sameSegment(before, after, seg);
                boolean inWindow = seg >= moved - 3 && seg <= moved + 2;
                if (inWindow) {
                    assertTrue(VARIANTS[v] + ": segment " + seg + " did not react", !same);
                } else {
                    assertTrue(VARIANTS[v] + ": segment " + seg + " reacted from " + Math.abs(seg - moved)
                            + " away", same);
                }
                if (!same) {
                    changed++;
                }
            }
            assertEquals(VARIANTS[v] + ": the support was not six segments", 6, changed);
        }

        // and the same disturbance in the natural spline, for contrast
        CubicSpline naturalBefore = SplineInterpolator.interpolate(x, y);
        CubicSpline naturalAfter = SplineInterpolator.interpolate(x, disturbed);
        for (int seg = 0; seg < n - 1; seg++) {
            assertTrue("the natural spline left segment " + seg + " alone",
                    !sameSegment(naturalBefore, naturalAfter, seg));
        }
    }

    @Test
    public void testAStraightLineIsReproducedExactly() {
        for (int v = 0; v < VARIANTS.length; v++) {
            for (int variant = 0; variant < 2; variant++) {
                int n = 7;
                double[] x = new double[n];
                double[] y = new double[n];
                for (int i = 0; i < n; i++) {
                    x[i] = (variant == 0) ? i : 0.3 * i * i + i;
                    y[i] = 2.5 * x[i] - 7.25;
                }
                CubicSpline s = AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
                double scale = Math.abs(2.5 * x[n - 1] - 7.25);
                for (int i = 0; i <= 2000; i++) {
                    double t = x[0] + (x[n - 1] - x[0]) * i / 2000.0;
                    assertEquals(VARIANTS[v] + ", variant " + variant, 2.5 * t - 7.25, s.value(t),
                            1.0e-14 * scale);
                }
            }
        }
    }

    @Test
    public void testTheClassicRuleReproducesAQuadraticOnAUniformGrid() {
        // on a uniform grid the two weights coincide and the rule collapses to
        // the average of the adjacent secants, which is the exact derivative
        // of a quadratic there; the modified weights break that and so do the
        // other two schemes in this package
        int n = 9;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i - 4;
            y[i] = 2.0 * x[i] * x[i] - 3.0 * x[i] + 1.0;
        }
        CubicSpline classic = AkimaInterpolator.interpolate(x, y);
        double scale = 0.0;
        double worstClassic = 0.0;
        double worstNatural = 0.0;
        CubicSpline natural = SplineInterpolator.interpolate(x, y);
        for (int i = 0; i <= 4000; i++) {
            double t = x[0] + (x[n - 1] - x[0]) * i / 4000.0;
            double exact = 2.0 * t * t - 3.0 * t + 1.0;
            scale = Math.max(scale, Math.abs(exact));
            worstClassic = Math.max(worstClassic, Math.abs(classic.value(t) - exact));
            worstNatural = Math.max(worstNatural, Math.abs(natural.value(t) - exact));
        }
        assertEquals("the quadratic was not reproduced", 0.0, worstClassic / scale, 1.0e-13);
        assertTrue("the natural spline reproduced it too, which makes this test pointless",
                worstNatural / scale > 1.0e-4);
    }

    @Test
    public void testTwoPointsGiveTheLineThroughThem() {
        double[] x = { 0.0, 1.0 };
        double[] y = { 2.0, 5.0 };
        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline s = AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
            assertEquals(1, s.segments());
            for (int i = 0; i <= 1000; i++) {
                double t = i / 1000.0;
                assertEquals(VARIANTS[v].toString(), 2.0 + 3.0 * t, s.value(t), 1.0e-14);
                assertEquals(VARIANTS[v].toString(), 3.0, s.derivativeAt(t), 1.0e-14);
            }
            assertEquals(3.5, s.integrate(0.0, 1.0), 1.0e-14);
        }
    }

    @Test
    public void testTheInterpolantIsOnceContinuouslyDifferentiable() {
        double[][][] sets = { smooth(11, 0.0), smooth(31, 0.0), monotone(9, 41L, true),
                { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 } } };
        for (double[][] data : sets) {
            for (int v = 0; v < VARIANTS.length; v++) {
                CubicSpline s = AkimaInterpolator.interpolate(data[0], data[1], VARIANTS[v]);
                for (int order = 0; order <= 1; order++) {
                    double scale = magnitude(s, order);
                    for (int i = 1; i < s.segments(); i++) {
                        double x = s.knot(i);
                        double jump = Math.abs(limit(s, i - 1, x, order) - limit(s, i, x, order));
                        assertEquals(VARIANTS[v] + ": derivative of order " + order + " jumps at knot " + i, 0.0,
                                jump / scale, 1.0e-12);
                    }
                }
            }
        }
    }

    @Test
    public void testTheSecondDerivativeDoesNotHaveToBeContinuous() {
        double[][] data = smooth(11, 0.0);
        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline s = AkimaInterpolator.interpolate(data[0], data[1], VARIANTS[v]);
            double scale = magnitude(s, 2);
            double largest = 0.0;
            for (int i = 1; i < s.segments(); i++) {
                double x = s.knot(i);
                largest = Math.max(largest, Math.abs(limit(s, i - 1, x, 2) - limit(s, i, x, 2)));
            }
            assertTrue(VARIANTS[v] + ": the second derivative was continuous after all: " + (largest / scale),
                    largest / scale > 0.1);
        }
    }

    @Test
    public void testTheModifiedRuleDoesNotTiltAPlateau() {
        // a ramp running into a plateau: both weights vanish at the corner, so
        // the classic rule averages the incoming slope with zero and lifts the
        // curve off the plateau. This is what the modified weighting is for.
        double[] x = { 0, 1, 2, 3, 4, 5, 6 };
        double[] y = { 0, 1, 2, 3, 3, 3, 3 };

        CubicSpline classic = AkimaInterpolator.interpolate(x, y, Variant.CLASSIC);
        CubicSpline modified = AkimaInterpolator.interpolate(x, y, Variant.MODIFIED);

        assertEquals("the classic rule takes the average at the corner", 0.5, classic.derivativeAt(3.0), 1.0e-14);
        assertEquals("the modified rule follows the plateau", 0.0, modified.derivativeAt(3.0), 1.0e-14);

        double classicTop = 0.0;
        double modifiedTop = 0.0;
        for (int i = 0; i <= 20000; i++) {
            double t = 6.0 * i / 20000.0;
            classicTop = Math.max(classicTop, classic.value(t));
            modifiedTop = Math.max(modifiedTop, modified.value(t));
        }
        assertTrue("the classic rule stayed on the plateau, which makes this test pointless",
                classicTop > 3.05);
        assertEquals("the modified rule left the plateau", 3.0, modifiedTop, 1.0e-12);
    }

    @Test
    public void testDataWithNothingToPreferGivesTheAverageAndNotANaN() {
        // every weight vanishes at every knot here
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 7, 7, 7, 7, 7 };
        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline s = AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
            for (int i = 0; i <= 4000; i++) {
                double t = 4.0 * i / 4000.0;
                assertEquals(VARIANTS[v].toString(), 7.0, s.value(t), 0.0);
                assertEquals(VARIANTS[v].toString(), 0.0, s.derivativeAt(t), 0.0);
            }
        }
    }

    @Test
    public void testMonotoneDataIsOnlyJustOvershot() {
        // Akima is not shape preserving the way KrugerInterpolator is, but it
        // is nowhere near the natural spline either; both bounds are asserted
        double worstAkima = 0.0;
        double worstNatural = 0.0;
        for (int seed = 41; seed < 47; seed++) {
            for (int uneven = 0; uneven < 2; uneven++) {
                double[][] data = monotone(9, seed, uneven == 1);
                double lo = Double.MAX_VALUE;
                double hi = -Double.MAX_VALUE;
                for (int i = 0; i < data[1].length; i++) {
                    lo = Math.min(lo, data[1][i]);
                    hi = Math.max(hi, data[1][i]);
                }
                double span = hi - lo;
                for (int v = 0; v < VARIANTS.length; v++) {
                    CubicSpline s = AkimaInterpolator.interpolate(data[0], data[1], VARIANTS[v]);
                    worstAkima = Math.max(worstAkima, excursion(s, lo, hi) / span);
                }
                worstNatural = Math.max(worstNatural,
                        excursion(SplineInterpolator.interpolate(data[0], data[1]), lo, hi) / span);
            }
        }
        assertTrue("Akima left the data by " + (100.0 * worstAkima) + "%", worstAkima < 0.02);
        assertTrue("the natural spline stayed inside the data, which makes the comparison pointless",
                worstNatural > 0.25);
    }

    private static double excursion(CubicSpline s, double lo, double hi) {
        double out = 0.0;
        for (int i = 0; i <= 20000; i++) {
            double t = s.lowerBound() + (s.upperBound() - s.lowerBound()) * i / 20000.0;
            double v = s.value(t);
            out = Math.max(out, Math.max(lo - v, v - hi));
        }
        return Math.max(out, 0.0);
    }

    @Test
    public void testTheInterpolationIsInvariantUnderAShiftOfTheAbscissa() {
        double[] shifts = { 1.0, 1.0e3, 1.0e6 };
        double[][] base = smooth(11, 0.0);
        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline reference = AkimaInterpolator.interpolate(base[0], base[1], VARIANTS[v]);
            for (int k = 0; k < shifts.length; k++) {
                double[][] moved = smooth(11, shifts[k]);
                CubicSpline s = AkimaInterpolator.interpolate(moved[0], moved[1], VARIANTS[v]);
                for (int i = 0; i <= 2000; i++) {
                    double t = 10.0 * i / 2000.0;
                    assertEquals(VARIANTS[v] + ", shift " + shifts[k], reference.value(t), s.value(t + shifts[k]),
                            1.0e-8);
                }
            }
        }
    }

    @Test
    public void testKnotsSpreadOverSixOrdersOfMagnitude() {
        double[] x = new double[10];
        double[] y = new double[10];
        double acc = 0.0;
        for (int i = 0; i < 10; i++) {
            x[i] = acc;
            y[i] = Math.cos(0.7 * i) + 0.2 * i;
            acc += Math.pow(10.0, (i % 7) - 3);
        }
        for (int v = 0; v < VARIANTS.length; v++) {
            CubicSpline s = AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
            for (int i = 0; i < x.length; i++) {
                assertEquals(VARIANTS[v] + ", knot " + i, y[i], s.value(x[i]), 1.0e-14);
            }
        }
    }

    @Test
    public void testTheConstructorAgreesWithTheFactoryBitForBit() {
        double[][] data = smooth(17, 0.0);
        AkimaInterpolator plain = new AkimaInterpolator(data[0], data[1]);
        CubicSpline classic = AkimaInterpolator.interpolate(data[0], data[1]);
        AkimaInterpolator explicit = new AkimaInterpolator(data[0], data[1], Variant.MODIFIED);
        CubicSpline modified = AkimaInterpolator.interpolate(data[0], data[1], Variant.MODIFIED);
        for (int i = 0; i <= 4000; i++) {
            double t = 10.0 * i / 4000.0;
            assertEquals(Double.doubleToRawLongBits(classic.value(t)), Double.doubleToRawLongBits(plain.value(t)));
            assertEquals(Double.doubleToRawLongBits(classic.value(t)),
                    Double.doubleToRawLongBits(plain.spline().value(t)));
            assertEquals(Double.doubleToRawLongBits(modified.value(t)),
                    Double.doubleToRawLongBits(explicit.value(t)));
        }
    }

    @Test
    public void testMalformedDataIsRejected() {
        double[] x = { 0.0, 1.0, 2.0, 3.0 };
        double[] y = { 0.0, 1.0, 0.0, 1.0 };
        rejects("null points", null, y);
        rejects("null values", x, null);
        rejects("length mismatch", x, new double[] { 0.0, 1.0, 2.0 });
        rejects("one point", new double[] { 1.0 }, new double[] { 5.0 });
        rejects("unsorted", new double[] { 0.0, 2.0, 1.0, 3.0 }, y);
        rejects("duplicate knot", new double[] { 0.0, 1.0, 1.0, 3.0 }, y);
        rejects("NaN knot", new double[] { 0.0, Double.NaN, 2.0, 3.0 }, y);
        rejects("infinite value", x, new double[] { 0.0, 1.0, Double.POSITIVE_INFINITY, 1.0 });
        try {
            AkimaInterpolator.interpolate(x, y, null);
            fail("a null variant was accepted by the factory");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
        try {
            new AkimaInterpolator(x, y, null);
            fail("a null variant was accepted by the constructor");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }

    private static void rejects(String what, double[] x, double[] y) {
        for (int v = 0; v < VARIANTS.length; v++) {
            try {
                AkimaInterpolator.interpolate(x, y, VARIANTS[v]);
                fail(what + " was accepted by the factory");
            } catch (IllegalArgumentException expected) {
                assertTrue(what + " threw without a message", expected.getMessage() != null);
            }
        }
        try {
            new AkimaInterpolator(x, y);
            fail(what + " was accepted by the constructor");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
