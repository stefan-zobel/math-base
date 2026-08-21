package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * What makes the natural cubic spline the natural cubic spline: two continuous
 * derivatives, a vanishing second derivative at both ends, and linearity in the
 * data.
 */
public class SplineInterpolatorTest {

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

    /** The one-sided limits at a knot, from the coefficients of the two pieces. */
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
    public void testTheInterpolantIsTwiceContinuouslyDifferentiable() {
        double[][][] sets = { smooth(11), smooth(31), wobble(12, 11L), wobble(41, 12L) };
        for (double[][] data : sets) {
            CubicSpline s = SplineInterpolator.interpolate(data[0], data[1]);
            for (int order = 0; order <= 2; order++) {
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
    public void testTheSecondDerivativeVanishesAtBothEnds() {
        // this is the condition the word "natural" names
        double[][][] sets = { smooth(11), smooth(4), wobble(12, 13L), wobble(41, 14L) };
        for (double[][] data : sets) {
            CubicSpline s = SplineInterpolator.interpolate(data[0], data[1]);
            double scale = magnitude(s, 2);
            assertEquals("at the left end", 0.0, s.secondDerivativeAt(s.lowerBound()) / scale, 1.0e-12);
            assertEquals("at the right end", 0.0, s.secondDerivativeAt(s.upperBound()) / scale, 1.0e-12);
        }
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
            CubicSpline s = SplineInterpolator.interpolate(x, y);
            double scale = Math.abs(2.5 * x[n - 1] - 7.25);
            for (int i = 0; i <= 2000; i++) {
                double t = x[0] + (x[n - 1] - x[0]) * i / 2000.0;
                assertEquals("variant " + variant, 2.5 * t - 7.25, s.value(t), 1.0e-14 * scale);
                assertEquals("variant " + variant, 2.5, s.derivativeAt(t), 1.0e-12);
            }
        }
    }

    @Test
    public void testTwoPointsGiveTheLineThroughThem() {
        double[] x = { 0.0, 1.0 };
        double[] y = { 2.0, 5.0 };
        CubicSpline s = SplineInterpolator.interpolate(x, y);
        assertEquals(1, s.segments());
        for (int i = 0; i <= 1000; i++) {
            double t = i / 1000.0;
            assertEquals(2.0 + 3.0 * t, s.value(t), 1.0e-14);
        }
        assertEquals(3.5, s.integrate(0.0, 1.0), 1.0e-14);
    }

    @Test
    public void testTheInterpolantIsLinearInTheData() {
        // the natural spline solves a linear system in the values, so the
        // interpolant of a sum is the sum of the interpolants -- which the
        // constrained spline of KrugerInterpolator is not
        double[][] a = wobble(13, 21L);
        double[][] b = wobble(13, 22L);
        double[] x = a[0];
        double[] sum = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            sum[i] = a[1][i] + b[1][i];
        }
        CubicSpline sa = SplineInterpolator.interpolate(x, a[1]);
        CubicSpline sb = SplineInterpolator.interpolate(x, b[1]);
        CubicSpline ss = SplineInterpolator.interpolate(x, sum);
        double scale = 0.0;
        for (int i = 0; i < sum.length; i++) {
            scale = Math.max(scale, Math.abs(sum[i]));
        }
        for (int i = 0; i <= 4000; i++) {
            double t = x[0] + (x[x.length - 1] - x[0]) * i / 4000.0;
            assertEquals(sa.value(t) + sb.value(t), ss.value(t), 1.0e-12 * scale);
        }
    }

    @Test
    public void testTheConstructorAgreesWithTheFactoryBitForBit() {
        double[][] data = wobble(15, 31L);
        SplineInterpolator interpolator = new SplineInterpolator(data[0], data[1]);
        CubicSpline s = SplineInterpolator.interpolate(data[0], data[1]);
        for (int i = 0; i <= 5000; i++) {
            double t = data[0][0] + (data[0][data[0].length - 1] - data[0][0]) * i / 5000.0;
            assertEquals(Double.doubleToRawLongBits(s.value(t)),
                    Double.doubleToRawLongBits(interpolator.value(t)));
            assertEquals(Double.doubleToRawLongBits(s.value(t)),
                    Double.doubleToRawLongBits(interpolator.spline().value(t)));
        }
    }

    @Test
    public void testTheInterpolationIsInvariantUnderAShiftOfTheAbscissa() {
        double[] shifts = { 1.0, 1.0e3, 1.0e6 };
        double[][] base = smooth(11);
        CubicSpline reference = SplineInterpolator.interpolate(base[0], base[1]);
        for (int k = 0; k < shifts.length; k++) {
            double[] moved = new double[base[0].length];
            for (int i = 0; i < moved.length; i++) {
                moved[i] = base[0][i] + shifts[k];
            }
            CubicSpline s = SplineInterpolator.interpolate(moved, base[1]);
            for (int i = 0; i <= 2000; i++) {
                double t = 10.0 * i / 2000.0;
                assertEquals("shift " + shifts[k], reference.value(t), s.value(t + shifts[k]), 1.0e-8);
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
        rejects("descending", new double[] { 3.0, 2.0, 1.0, 0.0 }, y);
        rejects("duplicate knot", new double[] { 0.0, 1.0, 1.0, 3.0 }, y);
        rejects("NaN knot", new double[] { 0.0, Double.NaN, 2.0, 3.0 }, y);
        rejects("infinite knot", new double[] { 0.0, 1.0, 2.0, Double.POSITIVE_INFINITY }, y);
        rejects("NaN value", x, new double[] { 0.0, Double.NaN, 0.0, 1.0 });
        rejects("infinite value", x, new double[] { 0.0, 1.0, Double.NEGATIVE_INFINITY, 1.0 });
    }

    private static void rejects(String what, double[] x, double[] y) {
        try {
            SplineInterpolator.interpolate(x, y);
            fail(what + " was accepted by the factory");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + " threw without a message", expected.getMessage() != null);
        }
        try {
            new SplineInterpolator(x, y);
            fail(what + " was accepted by the constructor");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
