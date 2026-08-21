package math.fit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fit.SuccessiveInterpolator.Scheme;

/**
 * The sweep that carries the non-linear one-dimensional rules into two
 * dimensions. What it costs is the tensor-product structure -- so the sweep
 * order matters and there is nothing to integrate in closed form -- and what it
 * buys is that the Kruger rule keeps the surface inside the data where every
 * other scheme here leaves it.
 */
public class SuccessiveInterpolatorTest {

    private static final Scheme[] SCHEMES = Scheme.values();

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

    private static double[][] transpose(double[][] z, int nx, int ny) {
        double[][] t = new double[ny][nx];
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                t[j][i] = z[i][j];
            }
        }
        return t;
    }

    @Test
    public void testEveryGridPointIsReproduced() {
        int[][] shapes = { { 7, 5 }, { 11, 9 }, { 2, 2 }, { 2, 6 }, { 6, 2 }, { 3, 3 } };
        for (int s = 0; s < shapes.length; s++) {
            double[] x = uneven(shapes[s][0], 101L + s);
            double[] y = uneven(shapes[s][1], 201L + s);
            double[][] z = smoothValues(x, y);
            for (int k = 0; k < SCHEMES.length; k++) {
                SuccessiveSurface surface = SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
                assertEquals(SCHEMES[k], surface.scheme());
                assertEquals(x[0], surface.lowerBoundX(), 0.0);
                assertEquals(x[x.length - 1], surface.upperBoundX(), 0.0);
                assertEquals(y[0], surface.lowerBoundY(), 0.0);
                assertEquals(y[y.length - 1], surface.upperBoundY(), 0.0);
                for (int i = 0; i < x.length; i++) {
                    for (int j = 0; j < y.length; j++) {
                        assertEquals(SCHEMES[k] + " at (" + i + ", " + j + ")", z[i][j], surface.value(x[i], y[j]),
                                1.0e-12);
                    }
                }
            }
        }
    }

    @Test
    public void testTheNaturalSweepIsDirectionIndependentAndTheOthersAreNot() {
        // the whole reason there are two constructions in this package: the
        // natural spline is linear in the data and the other rules are not
        double[] x = uneven(8, 31L);
        double[] y = uneven(6, 32L);
        double[][] z = smoothValues(x, y);
        double[][] zt = transpose(z, x.length, y.length);
        double natural = -1.0;
        double worstOther = 0.0;
        for (int k = 0; k < SCHEMES.length; k++) {
            SuccessiveSurface a = SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
            SuccessiveSurface b = SuccessiveInterpolator.interpolate(y, x, zt, SCHEMES[k]);
            double worst = 0.0;
            double scale = 0.0;
            for (int p = 0; p <= 30; p++) {
                for (int q = 0; q <= 30; q++) {
                    double px = at(x[0], x[x.length - 1], p, 30);
                    double py = at(y[0], y[y.length - 1], q, 30);
                    double u = a.value(px, py);
                    scale = Math.max(scale, Math.abs(u));
                    worst = Math.max(worst, Math.abs(u - b.value(py, px)));
                }
            }
            if (SCHEMES[k] == Scheme.NATURAL) {
                natural = worst / scale;
            } else {
                worstOther = Math.max(worstOther, worst / scale);
            }
        }
        assertEquals("the natural sweep depended on the direction", 0.0, natural, 1.0e-12);
        assertTrue("no scheme depended on the sweep direction, which makes the documented order meaningless: "
                + worstOther, worstOther > 1.0e-3);
    }

    @Test
    public void testTheKrugerSweepNeverLeavesTheRangeOfTheData() {
        // range preservation composes: each one-dimensional Kruger interpolant
        // stays within the values it was built from, so the intermediate values
        // stay in the global range and the second pass cannot escape it either
        double worstKruger = 0.0;
        double worstNatural = 0.0;
        for (int t = 0; t < 40; t++) {
            Lcg rnd = new Lcg(1000L + t);
            int nx = 4 + (int) (8 * rnd.next());
            int ny = 4 + (int) (8 * rnd.next());
            double[] x = new double[nx];
            double[] y = new double[ny];
            double acc = 0.0;
            for (int i = 0; i < nx; i++) {
                x[i] = acc;
                acc += 0.2 + 2.0 * rnd.next();
            }
            acc = 0.0;
            for (int j = 0; j < ny; j++) {
                y[j] = acc;
                acc += 0.2 + 2.0 * rnd.next();
            }
            double[][] z = new double[nx][ny];
            double lo = Double.MAX_VALUE;
            double hi = -Double.MAX_VALUE;
            int kind = t % 4;
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    double v;
                    if (kind == 0) {
                        v = rnd.next();
                    } else if (kind == 1) {
                        v = (i >= nx / 2 && j >= ny / 2) ? 1.0 : 0.0;
                    } else if (kind == 2) {
                        v = (i + j) % 3 == 0 ? 1.0 : 0.0;
                    } else {
                        v = Math.floor(4.0 * rnd.next()) / 4.0;
                    }
                    z[i][j] = v;
                    lo = Math.min(lo, v);
                    hi = Math.max(hi, v);
                }
            }
            double span = hi - lo;
            if (!(span > 0.0)) {
                continue;
            }
            SuccessiveSurface kruger = SuccessiveInterpolator.interpolate(x, y, z, Scheme.KRUGER);
            SuccessiveSurface natural = SuccessiveInterpolator.interpolate(x, y, z, Scheme.NATURAL);
            for (int p = 0; p <= 40; p++) {
                for (int q = 0; q <= 40; q++) {
                    double px = at(x[0], x[nx - 1], p, 40);
                    double py = at(y[0], y[ny - 1], q, 40);
                    double a = kruger.value(px, py);
                    double b = natural.value(px, py);
                    worstKruger = Math.max(worstKruger, Math.max(lo - a, a - hi) / span);
                    worstNatural = Math.max(worstNatural, Math.max(lo - b, b - hi) / span);
                }
            }
        }
        assertTrue("the Kruger sweep left the data by " + (100.0 * worstKruger) + "%", worstKruger < 1.0e-12);
        assertTrue("the natural sweep stayed inside the data, which makes the comparison pointless",
                worstNatural > 0.5);
    }

    @Test
    public void testTheAkimaSweepReproducesASeparableQuadraticThatTheBicubicDoesNot() {
        // classic Akima is exact for a quadratic on a uniform grid, and the
        // sweep of an exact rule through a separable function is exact too
        double[] x = uniform(9, -4.0, 4.0);
        double[] y = uniform(9, -4.0, 4.0);
        double[][] z = new double[9][9];
        double scale = 0.0;
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 9; j++) {
                z[i][j] = x[i] * x[i] * y[j] * y[j];
                scale = Math.max(scale, Math.abs(z[i][j]));
            }
        }
        SuccessiveSurface akima = SuccessiveInterpolator.interpolate(x, y, z, Scheme.AKIMA);
        PolynomialSurface bicubic = BicubicInterpolator.interpolate(x, y, z);
        double worstAkima = 0.0;
        double worstBicubic = 0.0;
        for (int p = 0; p <= 40; p++) {
            for (int q = 0; q <= 40; q++) {
                double px = at(-4.0, 4.0, p, 40);
                double py = at(-4.0, 4.0, q, 40);
                double exact = px * px * py * py;
                worstAkima = Math.max(worstAkima, Math.abs(akima.value(px, py) - exact));
                worstBicubic = Math.max(worstBicubic, Math.abs(bicubic.value(px, py) - exact));
            }
        }
        assertEquals("the Akima sweep did not reproduce it", 0.0, worstAkima / scale, 1.0e-13);
        assertTrue("the bicubic surface reproduced it too, which makes this test pointless",
                worstBicubic / scale > 1.0e-4);
    }

    @Test
    public void testABilinearFunctionIsReproducedExactlyByEveryScheme() {
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
            for (int k = 0; k < SCHEMES.length; k++) {
                SuccessiveSurface s = SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
                for (int p = 0; p <= 30; p++) {
                    for (int q = 0; q <= 30; q++) {
                        double px = at(x[0], x[x.length - 1], p, 30);
                        double py = at(y[0], y[y.length - 1], q, 30);
                        double exact = 2.0 - 0.5 * px + 1.5 * py + 0.75 * px * py;
                        assertEquals(SCHEMES[k] + ", grid " + unevenGrid, exact, s.value(px, py), 1.0e-13 * scale);
                    }
                }
            }
        }
    }

    @Test
    public void testTheFunctionInterfaceIsTheSameFunction() {
        double[] x = uneven(7, 81L);
        double[] y = uneven(5, 82L);
        double[][] z = smoothValues(x, y);
        for (int k = 0; k < SCHEMES.length; k++) {
            SuccessiveSurface s = SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
            for (int p = 0; p <= 20; p++) {
                for (int q = 0; q <= 20; q++) {
                    double px = at(s.lowerBoundX(), s.upperBoundX(), p, 20);
                    double py = at(s.lowerBoundY(), s.upperBoundY(), q, 20);
                    assertEquals(SCHEMES[k].toString(), Double.doubleToRawLongBits(s.value(px, py)),
                            Double.doubleToRawLongBits(s.apply(px, py)));
                }
            }
        }
    }

    @Test
    public void testAnArgumentOutsideTheGridIsRejected() {
        double[] x = uniform(5, 0.0, 4.0);
        double[] y = uniform(4, 0.0, 3.0);
        double[][] z = smoothValues(x, y);
        double[][] outside = { { -1.0e-9, 1.0 }, { 4.0 + 1.0e-9, 1.0 }, { 1.0, -1.0e-9 }, { 1.0, 3.0 + 1.0e-9 },
                { Double.NaN, 1.0 }, { 1.0, Double.NaN } };
        for (int k = 0; k < SCHEMES.length; k++) {
            SuccessiveSurface s = SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
            for (int i = 0; i < outside.length; i++) {
                try {
                    s.value(outside[i][0], outside[i][1]);
                    fail(SCHEMES[k] + " accepted (" + outside[i][0] + ", " + outside[i][1] + ")");
                } catch (IllegalArgumentException expected) {
                    // as intended
                }
            }
        }
    }

    @Test
    public void testMalformedArgumentsAreRejected() {
        double[] x = { 0.0, 1.0, 2.0 };
        double[] y = { 0.0, 1.0 };
        double[][] z = { { 0.0, 1.0 }, { 1.0, 2.0 }, { 2.0, 3.0 } };
        try {
            SuccessiveInterpolator.interpolate(x, y, z, null);
            fail("a null scheme was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage() != null);
        }
        rejects("null z", x, y, null);
        rejects("z too short", x, y, new double[][] { { 0.0, 1.0 }, { 1.0, 2.0 } });
        rejects("ragged z", x, y, new double[][] { { 0.0, 1.0 }, { 1.0, 2.0, 3.0 }, { 2.0, 3.0 } });
        rejects("one x", new double[] { 1.0 }, y, new double[][] { { 0.0, 1.0 } });
        rejects("unsorted y", x, new double[] { 1.0, 0.0 }, z);
        rejects("NaN in z", x, y, new double[][] { { 0.0, 1.0 }, { Double.NaN, 2.0 }, { 2.0, 3.0 } });
    }

    private static void rejects(String what, double[] x, double[] y, double[][] z) {
        for (int k = 0; k < SCHEMES.length; k++) {
            try {
                SuccessiveInterpolator.interpolate(x, y, z, SCHEMES[k]);
                fail(what + " was accepted for " + SCHEMES[k]);
            } catch (IllegalArgumentException expected) {
                assertTrue(what + " threw without a message", expected.getMessage() != null);
            }
        }
    }
}
