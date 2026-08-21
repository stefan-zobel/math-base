package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.fit.SuccessiveInterpolator.Scheme;

/**
 * The tensor product of natural splines. Its correctness rests on the two
 * constructions of the same surface -- the coefficient sweep here and the value
 * sweep of {@link SuccessiveInterpolator} -- agreeing, which they can only do
 * if the linearity argument that licenses the first one holds.
 */
public class BicubicInterpolatorTest {

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

    private static double at(double lo, double hi, int p, int res) {
        return Math.min(lo + (hi - lo) * p / res, hi);
    }

    private static double[] uniform(int n, double lo, double hi) {
        double[] t = new double[n];
        for (int i = 0; i < n; i++) {
            t[i] = at(lo, hi, i, n - 1);
        }
        return t;
    }

    private static double[] uneven(int n, long seed) {
        Lcg rnd = new Lcg(seed);
        double[] t = new double[n];
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            t[i] = acc;
            acc += 0.3 + rnd.next();
        }
        return t;
    }

    private static double[][] smoothValues(double[] x, double[] y) {
        double[][] z = new double[x.length][y.length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                z[i][j] = Math.sin(x[i]) * Math.cos(y[j]) + 0.2 * x[i] * y[j];
            }
        }
        return z;
    }

    /** A grid whose values step from 0 to 1 in the far quadrant. */
    private static double[][] stepValues(int nx, int ny) {
        double[][] z = new double[nx][ny];
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                z[i][j] = (i >= nx / 2 && j >= ny / 2) ? 1.0 : 0.0;
            }
        }
        return z;
    }

    @Test
    public void testTheCoefficientSweepAndTheValueSweepAgree() {
        // BicubicInterpolator interpolates the coefficients of the first pass;
        // SuccessiveInterpolator interpolates the values. The two are the same
        // surface only because the natural spline is linear in its data, and
        // they share no line of code
        int[][] shapes = { { 7, 5 }, { 11, 9 }, { 2, 6 }, { 6, 2 }, { 2, 2 }, { 3, 4 } };
        for (int s = 0; s < shapes.length; s++) {
            for (int unevenGrid = 0; unevenGrid < 2; unevenGrid++) {
                double[] x = unevenGrid == 1 ? uneven(shapes[s][0], 11L) : uniform(shapes[s][0], 0.0, 4.0);
                double[] y = unevenGrid == 1 ? uneven(shapes[s][1], 12L) : uniform(shapes[s][1], 0.0, 3.0);
                double[][] z = smoothValues(x, y);

                PolynomialSurface tensor = BicubicInterpolator.interpolate(x, y, z);
                SuccessiveSurface sweep = SuccessiveInterpolator.interpolate(x, y, z, Scheme.NATURAL);

                double worst = 0.0;
                double scale = 0.0;
                for (int p = 0; p <= 40; p++) {
                    for (int q = 0; q <= 40; q++) {
                        double px = at(x[0], x[x.length - 1], p, 40);
                        double py = at(y[0], y[y.length - 1], q, 40);
                        double u = tensor.value(px, py);
                        scale = Math.max(scale, Math.abs(u));
                        worst = Math.max(worst, Math.abs(u - sweep.value(px, py)));
                    }
                }
                assertEquals(shapes[s][0] + "x" + shapes[s][1] + ", grid " + unevenGrid, 0.0, worst / scale,
                        1.0e-12);
            }
        }
    }

    @Test
    public void testTheResultDoesNotDependOnWhichDirectionIsSweptFirst() {
        // the transposed problem must give the transposed surface, because the
        // natural spline is linear in the data
        int[][] shapes = { { 7, 5 }, { 9, 9 }, { 4, 6 } };
        for (int s = 0; s < shapes.length; s++) {
            double[] x = uneven(shapes[s][0], 21L);
            double[] y = uneven(shapes[s][1], 22L);
            double[][] z = smoothValues(x, y);
            double[][] transposed = new double[y.length][x.length];
            for (int i = 0; i < x.length; i++) {
                for (int j = 0; j < y.length; j++) {
                    transposed[j][i] = z[i][j];
                }
            }
            PolynomialSurface a = BicubicInterpolator.interpolate(x, y, z);
            PolynomialSurface b = BicubicInterpolator.interpolate(y, x, transposed);
            double worst = 0.0;
            double scale = 0.0;
            for (int p = 0; p <= 40; p++) {
                for (int q = 0; q <= 40; q++) {
                    double px = at(x[0], x[x.length - 1], p, 40);
                    double py = at(y[0], y[y.length - 1], q, 40);
                    double u = a.value(px, py);
                    scale = Math.max(scale, Math.abs(u));
                    worst = Math.max(worst, Math.abs(u - b.value(py, px)));
                }
            }
            assertEquals(shapes[s][0] + "x" + shapes[s][1], 0.0, worst / scale, 1.0e-12);
        }
    }

    @Test
    public void testTheSurfaceIsSmoothAcrossTheGridLines() {
        // twice continuously differentiable in each variable, inherited from
        // the natural spline in both sweeps
        double[] x = uneven(8, 31L);
        double[] y = uneven(7, 32L);
        PolynomialSurface s = BicubicInterpolator.interpolate(x, y, smoothValues(x, y));
        double scale = 0.0;
        for (int p = 0; p <= 20; p++) {
            for (int q = 0; q <= 20; q++) {
                scale = Math.max(scale, Math.abs(s.dx(at(x[0], x[x.length - 1], p, 20),
                        at(y[0], y[y.length - 1], q, 20))));
            }
        }
        double eps = 1.0e-7;
        for (int i = 1; i < s.cellsX(); i++) {
            double line = s.knotX(i);
            for (int q = 1; q < 20; q++) {
                double py = at(y[0], y[y.length - 1], q, 20);
                assertEquals("value across the grid line at x = " + line, s.value(line - eps, py),
                        s.value(line + eps, py), 1.0e-6 * scale);
                assertEquals("dx across the grid line at x = " + line, s.dx(line - eps, py),
                        s.dx(line + eps, py), 1.0e-5 * scale);
            }
        }
    }

    @Test
    public void testTheSurfaceOvershootsTheDataAndTheBilinearOneDoesNot() {
        // the price of smoothness, and the reason BilinearInterpolator exists
        double worstBicubic = 0.0;
        double worstBilinear = 0.0;
        int[][] shapes = { { 7, 7 }, { 9, 5 }, { 11, 11 } };
        for (int s = 0; s < shapes.length; s++) {
            int nx = shapes[s][0];
            int ny = shapes[s][1];
            double[] x = uniform(nx, 0.0, nx - 1.0);
            double[] y = uniform(ny, 0.0, ny - 1.0);
            double[][] z = stepValues(nx, ny);
            PolynomialSurface bicubic = BicubicInterpolator.interpolate(x, y, z);
            PolynomialSurface bilinear = BilinearInterpolator.interpolate(x, y, z);
            for (int p = 0; p <= 120; p++) {
                for (int q = 0; q <= 120; q++) {
                    double px = at(x[0], x[nx - 1], p, 120);
                    double py = at(y[0], y[ny - 1], q, 120);
                    double a = bicubic.value(px, py);
                    double b = bilinear.value(px, py);
                    worstBicubic = Math.max(worstBicubic, Math.max(-a, a - 1.0));
                    worstBilinear = Math.max(worstBilinear, Math.max(-b, b - 1.0));
                }
            }
        }
        assertTrue("the bicubic surface stayed inside the data, which makes the comparison pointless",
                worstBicubic > 0.1);
        assertTrue("the bilinear surface left the data by " + worstBilinear, worstBilinear < 1.0e-12);
    }

    @Test
    public void testTheInterpolationIsInvariantUnderAShiftOfBothAxes() {
        double[] shifts = { 1.0, 1.0e3, 1.0e6 };
        double[] x = uniform(9, 0.0, 4.0);
        double[] y = uniform(7, 0.0, 3.0);
        double[][] z = smoothValues(x, y);
        PolynomialSurface reference = BicubicInterpolator.interpolate(x, y, z);
        for (int k = 0; k < shifts.length; k++) {
            double sh = shifts[k];
            double[] mx = new double[x.length];
            double[] my = new double[y.length];
            for (int i = 0; i < mx.length; i++) {
                mx[i] = x[i] + sh;
            }
            for (int j = 0; j < my.length; j++) {
                my[j] = y[j] + sh;
            }
            PolynomialSurface moved = BicubicInterpolator.interpolate(mx, my, z);
            for (int p = 0; p <= 30; p++) {
                for (int q = 0; q <= 30; q++) {
                    double px = at(x[0], x[x.length - 1], p, 30);
                    double py = at(y[0], y[y.length - 1], q, 30);
                    assertEquals("shift " + sh, reference.value(px, py), moved.value(px + sh, py + sh), 1.0e-8);
                }
            }
        }
    }

    @Test
    public void testAGridWhoseSpacingsSpanOrdersOfMagnitude() {
        double[] x = new double[9];
        double[] y = new double[7];
        double acc = 0.0;
        for (int i = 0; i < x.length; i++) {
            x[i] = acc;
            acc += Math.pow(10.0, (i % 6) - 3);
        }
        acc = 0.0;
        for (int j = 0; j < y.length; j++) {
            y[j] = acc;
            acc += Math.pow(10.0, (j % 5) - 2);
        }
        double[][] z = smoothValues(x, y);
        PolynomialSurface s = BicubicInterpolator.interpolate(x, y, z);
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                assertEquals("grid point (" + i + ", " + j + ")", z[i][j], s.value(x[i], y[j]), 1.0e-12);
            }
        }
    }

    @Test
    public void testTheConstructorAgreesWithTheFactoryBitForBit() {
        double[] x = uneven(7, 41L);
        double[] y = uneven(6, 42L);
        double[][] z = smoothValues(x, y);
        BicubicInterpolator bicubic = new BicubicInterpolator(x, y, z);
        PolynomialSurface factory = BicubicInterpolator.interpolate(x, y, z);
        BilinearInterpolator bilinear = new BilinearInterpolator(x, y, z);
        PolynomialSurface bilinearFactory = BilinearInterpolator.interpolate(x, y, z);
        for (int p = 0; p <= 30; p++) {
            for (int q = 0; q <= 30; q++) {
                double px = at(x[0], x[x.length - 1], p, 30);
                double py = at(y[0], y[y.length - 1], q, 30);
                assertEquals(Double.doubleToRawLongBits(factory.value(px, py)),
                        Double.doubleToRawLongBits(bicubic.value(px, py)));
                assertEquals(Double.doubleToRawLongBits(factory.value(px, py)),
                        Double.doubleToRawLongBits(bicubic.surface().value(px, py)));
                assertEquals(Double.doubleToRawLongBits(bilinearFactory.value(px, py)),
                        Double.doubleToRawLongBits(bilinear.value(px, py)));
                assertEquals(Double.doubleToRawLongBits(bilinearFactory.value(px, py)),
                        Double.doubleToRawLongBits(bilinear.surface().value(px, py)));
            }
        }
    }
}
