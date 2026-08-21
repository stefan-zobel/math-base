package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;

/**
 * The evaluation machinery both grid interpolators share: value, the two first
 * partial derivatives, the closed-form double integral and the argument checks.
 * Bilinear interpolation is checked here too, since it is the same surface with
 * the higher coefficients zero.
 */
public class PolynomialSurfaceTest {

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

    /** A sample point that cannot round past the far end of the axis. */
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

    /** Both closed-form builders over the same grid. */
    private static PolynomialSurface[] both(double[] x, double[] y, double[][] z) {
        return new PolynomialSurface[] { BilinearInterpolator.interpolate(x, y, z),
                BicubicInterpolator.interpolate(x, y, z) };
    }

    private static String name(int i) {
        return i == 0 ? "bilinear surface" : "bicubic surface";
    }

    @Test
    public void testEveryGridPointIsReproduced() {
        int[][] shapes = { { 7, 5 }, { 11, 9 }, { 2, 2 }, { 2, 6 }, { 6, 2 }, { 3, 3 } };
        for (int s = 0; s < shapes.length; s++) {
            double[] x = uniform(shapes[s][0], 0.0, 4.0);
            double[] y = uniform(shapes[s][1], 0.0, 3.0);
            double[][] z = smoothValues(x, y);
            PolynomialSurface[] surfaces = both(x, y, z);
            for (int k = 0; k < surfaces.length; k++) {
                for (int i = 0; i < x.length; i++) {
                    for (int j = 0; j < y.length; j++) {
                        assertEquals(name(k) + " at (" + i + ", " + j + ")", z[i][j], surfaces[k].value(x[i], y[j]),
                                1.0e-12);
                    }
                }
                assertEquals(name(k), x.length - 1, surfaces[k].cellsX());
                assertEquals(name(k), y.length - 1, surfaces[k].cellsY());
                assertEquals(name(k), x[0], surfaces[k].lowerBoundX(), 0.0);
                assertEquals(name(k), x[x.length - 1], surfaces[k].upperBoundX(), 0.0);
                assertEquals(name(k), y[0], surfaces[k].lowerBoundY(), 0.0);
                assertEquals(name(k), y[y.length - 1], surfaces[k].upperBoundY(), 0.0);
            }
        }
    }

    @Test
    public void testThePublishedCoefficientsDescribeTheCellTheyBelongTo() {
        double[] x = uneven(8, 31L);
        double[] y = uneven(6, 32L);
        double[][] z = smoothValues(x, y);
        PolynomialSurface[] surfaces = both(x, y, z);
        for (int k = 0; k < surfaces.length; k++) {
            PolynomialSurface s = surfaces[k];
            for (int i = 0; i < s.cellsX(); i++) {
                for (int j = 0; j < s.cellsY(); j++) {
                    double[] c = s.coefficients(i, j);
                    assertEquals(16, c.length);
                    for (int p = 0; p <= 4; p++) {
                        for (int q = 0; q <= 4; q++) {
                            double px = at(s.knotX(i), s.knotX(i + 1), p, 4);
                            double py = at(s.knotY(j), s.knotY(j + 1), q, 4);
                            double t = px - s.knotX(i);
                            double v = py - s.knotY(j);
                            double expected = 0.0;
                            double ta = 1.0;
                            for (int a = 0; a < 4; a++) {
                                double vb = 1.0;
                                for (int b = 0; b < 4; b++) {
                                    expected += c[a * 4 + b] * ta * vb;
                                    vb *= v;
                                }
                                ta *= t;
                            }
                            assertEquals(name(k) + " cell (" + i + ", " + j + ")", expected, s.value(px, py),
                                    1.0e-12);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testTheClosedFormIntegralAgreesWithGaussKronrodOnEveryCell() {
        // a fifteen point Kronrod rule in each direction is exact for a
        // bicubic over one cell, so this is an independent evaluation
        int[][] shapes = { { 7, 5 }, { 11, 9 } };
        for (int s = 0; s < shapes.length; s++) {
            for (int uneven = 0; uneven < 2; uneven++) {
                double[] x = uneven == 1 ? uneven(shapes[s][0], 41L) : uniform(shapes[s][0], 0.0, 4.0);
                double[] y = uneven == 1 ? uneven(shapes[s][1], 42L) : uniform(shapes[s][1], 0.0, 3.0);
                double[][] z = smoothValues(x, y);
                PolynomialSurface[] surfaces = both(x, y, z);
                for (int k = 0; k < surfaces.length; k++) {
                    PolynomialSurface surface = surfaces[k];
                    double kronrod = 0.0;
                    double area = 0.0;
                    for (int i = 0; i < surface.cellsX(); i++) {
                        for (int j = 0; j < surface.cellsY(); j++) {
                            kronrod += AdaptiveGaussKronrod.integrate2DParallel(G7_K15.POINTS_15, surface,
                                    surface.knotX(i), surface.knotX(i + 1), surface.knotY(j),
                                    surface.knotY(j + 1)).value;
                            area += Math.abs(surface.integrate(surface.knotX(i), surface.knotX(i + 1),
                                    surface.knotY(j), surface.knotY(j + 1)));
                        }
                    }
                    double closed = surface.integrate(surface.lowerBoundX(), surface.upperBoundX(),
                            surface.lowerBoundY(), surface.upperBoundY());
                    assertEquals(name(k), 0.0, Math.abs(closed - kronrod) / area, 1.0e-12);
                }
            }
        }
    }

    @Test
    public void testTheIntegralIsAdditiveAndAntisymmetric() {
        double[] x = uneven(9, 51L);
        double[] y = uneven(7, 52L);
        double[][] z = smoothValues(x, y);
        PolynomialSurface[] surfaces = both(x, y, z);
        for (int k = 0; k < surfaces.length; k++) {
            PolynomialSurface s = surfaces[k];
            double ax = s.lowerBoundX();
            double bx = s.upperBoundX();
            double ay = s.lowerBoundY();
            double by = s.upperBoundY();
            double area = 0.0;
            for (int i = 0; i < s.cellsX(); i++) {
                for (int j = 0; j < s.cellsY(); j++) {
                    area += Math.abs(s.integrate(s.knotX(i), s.knotX(i + 1), s.knotY(j), s.knotY(j + 1)));
                }
            }
            double whole = s.integrate(ax, bx, ay, by);
            assertEquals(name(k) + " antisymmetry in x", -whole, s.integrate(bx, ax, ay, by), 1.0e-14 * area);
            assertEquals(name(k) + " antisymmetry in y", -whole, s.integrate(ax, bx, by, ay), 1.0e-14 * area);
            assertEquals(name(k) + " both flipped", whole, s.integrate(bx, ax, by, ay), 1.0e-14 * area);
            assertEquals(name(k) + " degenerate in x", 0.0, s.integrate(ax, ax, ay, by), 0.0);
            assertEquals(name(k) + " degenerate in y", 0.0, s.integrate(ax, bx, ay, ay), 0.0);
            for (int p = 1; p < 8; p++) {
                double mx = at(ax, bx, p, 8);
                for (int q = 1; q < 8; q++) {
                    double my = at(ay, by, q, 8);
                    double quarters = s.integrate(ax, mx, ay, my) + s.integrate(mx, bx, ay, my)
                            + s.integrate(ax, mx, my, by) + s.integrate(mx, bx, my, by);
                    assertEquals(name(k) + " split at (" + mx + ", " + my + ")", whole, quarters, 1.0e-12 * area);
                }
            }
        }
    }

    @Test
    public void testThePartialDerivativesAgreeWithCentralDifferences() {
        double[] x = uneven(8, 61L);
        double[] y = uneven(6, 62L);
        double[][] z = smoothValues(x, y);
        PolynomialSurface[] surfaces = both(x, y, z);
        for (int k = 0; k < surfaces.length; k++) {
            PolynomialSurface s = surfaces[k];
            double worstX = 0.0;
            double worstY = 0.0;
            double scaleX = 0.0;
            double scaleY = 0.0;
            for (int i = 0; i < s.cellsX(); i++) {
                for (int j = 0; j < s.cellsY(); j++) {
                    double hx = (s.knotX(i + 1) - s.knotX(i)) * 1.0e-4;
                    double hy = (s.knotY(j + 1) - s.knotY(j)) * 1.0e-4;
                    for (int p = 1; p <= 3; p++) {
                        for (int q = 1; q <= 3; q++) {
                            double px = at(s.knotX(i), s.knotX(i + 1), p, 4);
                            double py = at(s.knotY(j), s.knotY(j + 1), q, 4);
                            double cdx = (s.value(px + hx, py) - s.value(px - hx, py)) / (2.0 * hx);
                            double cdy = (s.value(px, py + hy) - s.value(px, py - hy)) / (2.0 * hy);
                            scaleX = Math.max(scaleX, Math.abs(cdx));
                            scaleY = Math.max(scaleY, Math.abs(cdy));
                            worstX = Math.max(worstX, Math.abs(cdx - s.dx(px, py)));
                            worstY = Math.max(worstY, Math.abs(cdy - s.dy(px, py)));
                        }
                    }
                }
            }
            assertEquals(name(k) + " dx", 0.0, worstX / scaleX, 1.0e-6);
            assertEquals(name(k) + " dy", 0.0, worstY / scaleY, 1.0e-6);
        }
    }

    @Test
    public void testABilinearFunctionIsReproducedExactlyByBoth() {
        for (int unevenGrid = 0; unevenGrid < 2; unevenGrid++) {
            double[] x = unevenGrid == 1 ? uneven(6, 71L) : uniform(6, 0.0, 4.0);
            double[] y = unevenGrid == 1 ? uneven(5, 72L) : uniform(5, 0.0, 3.0);
            double[][] z = new double[x.length][y.length];
            double scale = 0.0;
            for (int i = 0; i < x.length; i++) {
                for (int j = 0; j < y.length; j++) {
                    z[i][j] = 2.0 - 0.5 * x[i] + 1.5 * y[j] + 0.75 * x[i] * y[j];
                    scale = Math.max(scale, Math.abs(z[i][j]));
                }
            }
            PolynomialSurface[] surfaces = both(x, y, z);
            for (int k = 0; k < surfaces.length; k++) {
                for (int p = 0; p <= 40; p++) {
                    for (int q = 0; q <= 40; q++) {
                        double px = at(x[0], x[x.length - 1], p, 40);
                        double py = at(y[0], y[y.length - 1], q, 40);
                        double exact = 2.0 - 0.5 * px + 1.5 * py + 0.75 * px * py;
                        assertEquals(name(k) + ", grid " + unevenGrid, exact, surfaces[k].value(px, py),
                                1.0e-13 * scale);
                    }
                }
            }
        }
    }

    @Test
    public void testTheBilinearSurfaceCannotLeaveTheRangeOfTheData() {
        // every value it produces is a convex combination of four grid values,
        // which is the one guarantee no other scheme in this package gives
        double worst = 0.0;
        for (int t = 0; t < 40; t++) {
            Lcg rnd = new Lcg(7000L + t);
            int nx = 4 + (int) (6 * rnd.next());
            int ny = 4 + (int) (6 * rnd.next());
            double[] x = new double[nx];
            double[] y = new double[ny];
            double acc = 0.0;
            for (int i = 0; i < nx; i++) {
                x[i] = acc;
                acc += 0.2 + rnd.next();
            }
            acc = 0.0;
            for (int j = 0; j < ny; j++) {
                y[j] = acc;
                acc += 0.2 + rnd.next();
            }
            double[][] z = new double[nx][ny];
            double lo = Double.MAX_VALUE;
            double hi = -Double.MAX_VALUE;
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    z[i][j] = 5.0 * rnd.next() - 2.0;
                    lo = Math.min(lo, z[i][j]);
                    hi = Math.max(hi, z[i][j]);
                }
            }
            PolynomialSurface s = BilinearInterpolator.interpolate(x, y, z);
            for (int p = 0; p <= 50; p++) {
                for (int q = 0; q <= 50; q++) {
                    double v = s.value(at(x[0], x[nx - 1], p, 50), at(y[0], y[ny - 1], q, 50));
                    worst = Math.max(worst, Math.max(lo - v, v - hi) / (hi - lo));
                }
            }
        }
        assertTrue("the bilinear surface left the data by " + worst, worst < 1.0e-12);
    }

    @Test
    public void testAnArgumentOutsideTheGridIsRejected() {
        double[] x = uniform(5, 0.0, 4.0);
        double[] y = uniform(4, 0.0, 3.0);
        PolynomialSurface[] surfaces = both(x, y, smoothValues(x, y));
        double[][] outside = { { -1.0e-9, 1.0 }, { 4.0 + 1.0e-9, 1.0 }, { 1.0, -1.0e-9 }, { 1.0, 3.0 + 1.0e-9 },
                { Double.NaN, 1.0 }, { 1.0, Double.NaN }, { Double.POSITIVE_INFINITY, 1.0 },
                { 1.0, Double.NEGATIVE_INFINITY } };
        for (int k = 0; k < surfaces.length; k++) {
            for (int i = 0; i < outside.length; i++) {
                rejects(name(k) + " value", surfaces[k], 0, outside[i][0], outside[i][1]);
                rejects(name(k) + " dx", surfaces[k], 1, outside[i][0], outside[i][1]);
                rejects(name(k) + " dy", surfaces[k], 2, outside[i][0], outside[i][1]);
                rejects(name(k) + " integrate", surfaces[k], 3, outside[i][0], outside[i][1]);
            }
        }
    }

    private static void rejects(String what, PolynomialSurface s, int which, double px, double py) {
        try {
            switch (which) {
            case 0:
                s.value(px, py);
                break;
            case 1:
                s.dx(px, py);
                break;
            case 2:
                s.dy(px, py);
                break;
            default:
                s.integrate(s.lowerBoundX(), px, s.lowerBoundY(), py);
                break;
            }
            fail(what + " accepted (" + px + ", " + py + ")");
        } catch (IllegalArgumentException expected) {
            // that is the point
        }
    }

    @Test
    public void testTheAccessorsHandOutCopiesAndCheckTheirIndices() {
        double[] x = { 0.0, 1.0, 2.0 };
        double[] y = { 0.0, 1.0 };
        double[][] z = { { 0.0, 1.0 }, { 1.0, 3.0 }, { 2.0, 4.0 } };
        PolynomialSurface s = BicubicInterpolator.interpolate(x, y, z);
        double before = s.value(1.5, 0.5);
        x[1] = 1.9;
        z[1][0] = 99.0;
        assertEquals("mutating the input changed the surface", before, s.value(1.5, 0.5), 0.0);
        double[] c = s.coefficients(0, 0);
        assertNotSame(c, s.coefficients(0, 0));
        c[0] = -5.0;
        assertEquals(0.0, s.coefficients(0, 0)[0], 0.0);

        int[] badX = { -1, 3, 4 };
        for (int i = 0; i < badX.length; i++) {
            try {
                s.knotX(badX[i]);
                fail("knotX(" + badX[i] + ") accepted");
            } catch (IllegalArgumentException expected) {
                // as intended
            }
        }
        int[] badCellX = { -1, 2, 3 };
        for (int i = 0; i < badCellX.length; i++) {
            try {
                s.coefficients(badCellX[i], 0);
                fail("coefficients(" + badCellX[i] + ", 0) accepted");
            } catch (IllegalArgumentException expected) {
                // as intended
            }
        }
        try {
            s.coefficients(0, 1);
            fail("coefficients(0, 1) accepted on a surface with one y cell");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
        try {
            s.knotY(2);
            fail("knotY(2) accepted");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }

    @Test
    public void testTheFunctionInterfaceIsTheSameFunction() {
        double[] x = uneven(7, 81L);
        double[] y = uneven(5, 82L);
        PolynomialSurface[] surfaces = both(x, y, smoothValues(x, y));
        for (int k = 0; k < surfaces.length; k++) {
            PolynomialSurface s = surfaces[k];
            for (int p = 0; p <= 30; p++) {
                for (int q = 0; q <= 30; q++) {
                    double px = at(s.lowerBoundX(), s.upperBoundX(), p, 30);
                    double py = at(s.lowerBoundY(), s.upperBoundY(), q, 30);
                    assertEquals(name(k), Double.doubleToRawLongBits(s.value(px, py)),
                            Double.doubleToRawLongBits(s.apply(px, py)));
                }
            }
        }
    }

    @Test
    public void testMalformedGridsAreRejected() {
        double[] x = { 0.0, 1.0, 2.0 };
        double[] y = { 0.0, 1.0 };
        double[][] z = { { 0.0, 1.0 }, { 1.0, 2.0 }, { 2.0, 3.0 } };
        rejectsGrid("null x", null, y, z);
        rejectsGrid("null y", x, null, z);
        rejectsGrid("null z", x, y, null);
        rejectsGrid("z too short", x, y, new double[][] { { 0.0, 1.0 }, { 1.0, 2.0 } });
        rejectsGrid("ragged z", x, y, new double[][] { { 0.0, 1.0 }, { 1.0, 2.0, 3.0 }, { 2.0, 3.0 } });
        rejectsGrid("null row", x, y, new double[][] { { 0.0, 1.0 }, null, { 2.0, 3.0 } });
        rejectsGrid("one x", new double[] { 1.0 }, y, new double[][] { { 0.0, 1.0 } });
        rejectsGrid("one y", x, new double[] { 1.0 },
                new double[][] { { 0.0 }, { 1.0 }, { 2.0 } });
        rejectsGrid("unsorted x", new double[] { 0.0, 2.0, 1.0 }, y, z);
        rejectsGrid("duplicate y", x, new double[] { 1.0, 1.0 }, z);
        rejectsGrid("NaN x", new double[] { 0.0, Double.NaN, 2.0 }, y, z);
        rejectsGrid("infinite y", x, new double[] { 0.0, Double.POSITIVE_INFINITY }, z);
        rejectsGrid("NaN in z", x, y, new double[][] { { 0.0, 1.0 }, { Double.NaN, 2.0 }, { 2.0, 3.0 } });
    }

    private static void rejectsGrid(String what, double[] x, double[] y, double[][] z) {
        try {
            BilinearInterpolator.interpolate(x, y, z);
            fail(what + " accepted by the bilinear factory");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + " threw without a message", expected.getMessage() != null);
        }
        try {
            BicubicInterpolator.interpolate(x, y, z);
            fail(what + " accepted by the bicubic factory");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
        try {
            new BicubicInterpolator(x, y, z);
            fail(what + " accepted by the constructor");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
