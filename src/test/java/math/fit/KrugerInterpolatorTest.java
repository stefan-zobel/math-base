package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * What makes the constrained cubic spline worth having beside the natural one:
 * it does not overshoot the data. And the two ways it used to get the answer
 * wrong, one of them silently.
 */
public class KrugerInterpolatorTest {

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

    @Test
    public void testTheInterpolantIsOnceContinuouslyDifferentiable() {
        double[][][] sets = { smooth(11, 0.0), smooth(31, 0.0), monotone(9, 41L, false), monotone(9, 42L, true) };
        for (double[][] data : sets) {
            CubicSpline s = KrugerInterpolator.interpolate(data[0], data[1]);
            for (int order = 0; order <= 1; order++) {
                double scale = magnitude(s, order);
                for (int i = 1; i < s.segments(); i++) {
                    double x = s.knot(i);
                    double jump = Math.abs(limit(s, i - 1, x, order) - limit(s, i, x, order));
                    assertEquals("derivative of order " + order + " jumps at knot " + i, 0.0, jump / scale,
                            1.0e-12);
                }
            }
        }
    }

    @Test
    public void testTheSecondDerivativeDoesNotHaveToBeContinuous() {
        // the price of not overshooting, and the reason both classes exist
        double[][] data = smooth(11, 0.0);
        CubicSpline s = KrugerInterpolator.interpolate(data[0], data[1]);
        double scale = magnitude(s, 2);
        double largest = 0.0;
        for (int i = 1; i < s.segments(); i++) {
            double x = s.knot(i);
            largest = Math.max(largest, Math.abs(limit(s, i - 1, x, 2) - limit(s, i, x, 2)));
        }
        assertTrue("the second derivative was continuous after all: " + (largest / scale),
                largest / scale > 0.1);
    }

    @Test
    public void testMonotoneDataIsNeverOvershot() {
        double[][][] sets = { { { 0, 1, 2, 3, 4, 5 }, { 0, 0, 0, 1, 1, 1 } },
                { { 0, 1, 2, 3, 4, 5, 6 }, { 0, 0.1, 0.15, 3.0, 3.05, 3.06, 8.0 } }, monotone(9, 51L, false),
                monotone(9, 52L, true), monotone(15, 53L, false) };
        for (double[][] data : sets) {
            CubicSpline s = KrugerInterpolator.interpolate(data[0], data[1]);
            double lo = Double.MAX_VALUE;
            double hi = -Double.MAX_VALUE;
            for (int i = 0; i < data[1].length; i++) {
                lo = Math.min(lo, data[1][i]);
                hi = Math.max(hi, data[1][i]);
            }
            double slack = 1.0e-12 * (hi - lo);
            for (int i = 0; i <= 20000; i++) {
                double t = s.lowerBound() + (s.upperBound() - s.lowerBound()) * i / 20000.0;
                double v = s.value(t);
                assertTrue("undershoot at " + t + " : " + v + " < " + lo, v >= lo - slack);
                assertTrue("overshoot at " + t + " : " + v + " > " + hi, v <= hi + slack);
            }
        }
    }

    @Test
    public void testTheSlopeIsZeroWhereTheSecantsChangeSign() {
        // the rule that keeps the interpolant inside the data
        double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0 };
        double[] y = { 0.0, 2.0, 1.0, 3.0, 2.5 };
        CubicSpline s = KrugerInterpolator.interpolate(x, y);
        assertEquals("at the peak", 0.0, s.derivativeAt(1.0), 0.0);
        assertEquals("at the trough", 0.0, s.derivativeAt(2.0), 0.0);
        assertEquals("at the second peak", 0.0, s.derivativeAt(3.0), 0.0);
    }

    @Test
    public void testAStraightLineIsReproducedExactly() {
        for (int variant = 0; variant < 2; variant++) {
            int n = 7;
            double[] x = new double[n];
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = (variant == 0) ? i : 0.3 * i * i + i;
                y[i] = 2.5 * x[i] - 7.25;
            }
            CubicSpline s = KrugerInterpolator.interpolate(x, y);
            double scale = Math.abs(2.5 * x[n - 1] - 7.25);
            for (int i = 0; i <= 2000; i++) {
                double t = x[0] + (x[n - 1] - x[0]) * i / 2000.0;
                assertEquals("variant " + variant, 2.5 * t - 7.25, s.value(t), 1.0e-14 * scale);
            }
        }
    }

    @Test
    public void testTwoPointsGiveTheLineThroughThem() {
        // the end slope formulas of the paper need an interior slope, and with
        // two points there is none: this used to return 3.78125 at the midpoint
        double[] x = { 0.0, 1.0 };
        double[] y = { 2.0, 5.0 };
        CubicSpline s = KrugerInterpolator.interpolate(x, y);
        assertEquals(1, s.segments());
        for (int i = 0; i <= 1000; i++) {
            double t = i / 1000.0;
            assertEquals(2.0 + 3.0 * t, s.value(t), 1.0e-14);
            assertEquals(3.0, s.derivativeAt(t), 1.0e-14);
        }
        assertEquals(3.5, s.integrate(0.0, 1.0), 1.0e-14);
    }

    @Test
    public void testTheInterpolationIsInvariantUnderAShiftOfTheAbscissa() {
        // the polynomials used to be written in the absolute abscissa, where
        // the constant term cancels the rest away: at a shift of 1e6 the
        // interpolant missed its own data points by 75
        double[][] base = smooth(11, 0.0);
        CubicSpline reference = KrugerInterpolator.interpolate(base[0], base[1]);
        double[] shifts = { 1.0, 1.0e3, 1.0e6 };
        for (int k = 0; k < shifts.length; k++) {
            double[][] moved = smooth(11, shifts[k]);
            CubicSpline s = KrugerInterpolator.interpolate(moved[0], moved[1]);
            for (int i = 0; i <= 2000; i++) {
                double t = 10.0 * i / 2000.0;
                assertEquals("shift " + shifts[k], reference.value(t), s.value(t + shifts[k]), 1.0e-8);
            }
            for (int i = 0; i < moved[0].length; i++) {
                assertEquals("shift " + shifts[k] + ", knot " + i, moved[1][i], s.value(moved[0][i]), 1.0e-14);
            }
        }
        // and at the scale of a millisecond timestamp
        double[][] epoch = smooth(11, 1.7e12);
        CubicSpline s = KrugerInterpolator.interpolate(epoch[0], epoch[1]);
        for (int i = 0; i < epoch[0].length; i++) {
            assertEquals("epoch knot " + i, epoch[1][i], s.value(epoch[0][i]), 1.0e-3);
        }
    }

    @Test
    public void testKnotsSpreadOverSixOrdersOfMagnitude() {
        // spacings of 1e-3 next to spacings of 1e3 in one data set; the
        // absolute-coordinate form returned -1024 for a value of 1.586512
        double[] x = new double[10];
        double[] y = new double[10];
        double acc = 0.0;
        for (int i = 0; i < 10; i++) {
            x[i] = acc;
            y[i] = Math.cos(0.7 * i) + 0.2 * i;
            acc += Math.pow(10.0, (i % 7) - 3);
        }
        CubicSpline s = KrugerInterpolator.interpolate(x, y);
        for (int i = 0; i < x.length; i++) {
            assertEquals("knot " + i, y[i], s.value(x[i]), 1.0e-14);
        }
    }

    @Test
    public void testTheConstructorAgreesWithTheFactoryBitForBit() {
        double[][] data = smooth(17, 0.0);
        KrugerInterpolator interpolator = new KrugerInterpolator(data[0], data[1]);
        CubicSpline s = KrugerInterpolator.interpolate(data[0], data[1]);
        for (int i = 0; i <= 5000; i++) {
            double t = 10.0 * i / 5000.0;
            assertEquals(Double.doubleToRawLongBits(s.value(t)),
                    Double.doubleToRawLongBits(interpolator.value(t)));
            assertEquals(Double.doubleToRawLongBits(s.value(t)),
                    Double.doubleToRawLongBits(interpolator.spline().value(t)));
        }
    }

    @Test
    public void testFlatDataIsHandledWithoutDividingByZero() {
        // the slope rule divides by the secant of the neighbouring interval
        double[] dys = { 1.0e-3, 1.0e-8, 1.0e-14, 1.0e-300, 0.0 };
        for (int k = 0; k < dys.length; k++) {
            double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
            double[] y = new double[6];
            for (int i = 0; i < 6; i++) {
                y[i] = 1.0 + i * dys[k];
            }
            CubicSpline s = KrugerInterpolator.interpolate(x, y);
            for (int i = 0; i <= 2000; i++) {
                double t = 5.0 * i / 2000.0;
                double v = s.value(t);
                assertTrue("dy = " + dys[k] + " gave " + v, !Double.isNaN(v) && !Double.isInfinite(v));
            }
            for (int i = 0; i < x.length; i++) {
                assertEquals("dy = " + dys[k] + ", knot " + i, y[i], s.value(x[i]), 0.0);
            }
        }
    }

    @Test
    public void testMalformedDataIsRejected() {
        double[] x = { 0.0, 1.0, 2.0, 3.0 };
        double[] y = { 0.0, 1.0, 0.0, 1.0 };
        rejects("null points", null, y);
        rejects("null values", x, null);
        rejects("length mismatch", x, new double[] { 0.0, 1.0, 2.0 });
        rejects("no points", new double[0], new double[0]);
        rejects("one point", new double[] { 1.0 }, new double[] { 5.0 });
        rejects("unsorted", new double[] { 0.0, 2.0, 1.0, 3.0 }, y);
        rejects("duplicate knot", new double[] { 0.0, 1.0, 1.0, 3.0 }, y);
        rejects("NaN knot", new double[] { 0.0, Double.NaN, 2.0, 3.0 }, y);
        rejects("infinite value", x, new double[] { 0.0, 1.0, Double.POSITIVE_INFINITY, 1.0 });
    }

    private static void rejects(String what, double[] x, double[] y) {
        try {
            KrugerInterpolator.interpolate(x, y);
            fail(what + " was accepted by the factory");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + " threw without a message", expected.getMessage() != null);
        }
        try {
            new KrugerInterpolator(x, y);
            fail(what + " was accepted by the constructor");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
