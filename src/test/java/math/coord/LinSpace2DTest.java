package math.coord;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.NoSuchElementException;

import org.junit.Test;

import math.fit.BilinearInterpolator;

/**
 * {@code LinSpace2D} is the tensor product of two {@code LinSpace} axes and
 * owns no abscissa arithmetic of its own. These tests are written around that
 * claim: what they check first is that the axes really are the ones that were
 * asked for, bit for bit, and that the four properties the one-dimensional
 * class was repaired for -- exact endpoints, a degenerate range that is not
 * folded onto one point, one grid behind every accessor, and a span beyond
 * {@code Double.MAX_VALUE} that stays finite -- are inherited rather than
 * reimplemented and broken a second time.
 * <p>
 * The rest covers what is genuinely two-dimensional: the column-major layout
 * with the first variable running fastest, the block copy behind {@code slice},
 * the cuts back into the one-dimensional world, and the cell lookup, which is
 * arithmetic rather than a search and is measured from whichever endpoint is
 * the nearer one -- the forward form alone walks up to 499 positions on a span
 * wider than the double range.
 */
public class LinSpace2DTest {

    // ----- layout ----------------------------------------------------------

    @Test
    public void theLayoutIsColumnMajorWithTheFirstVariableFastest() {
        LinSpace2D g = LinSpace2D.linspace2D(0.0, 1.0, 4, 10.0, 20.0, 3).allocate();
        int nx = g.sizeX();
        assertEquals(4, nx);
        assertEquals(3, g.sizeY());
        assertEquals(12, g.size());
        assertEquals(12, g.values().length);
        double v = 0.0;
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                g.setValue(i, j, v += 1.0);
            }
        }
        double[] flat = g.values();
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("value(" + i + ", " + j + ")", g.value(i, j), flat[(j - 1) * nx + (i - 1)]);
            }
        }
        // the first variable runs fastest, so the first sizeX() entries are
        // the first cut along x
        for (int i = 1; i <= nx; ++i) {
            assertSameBits("first column", g.value(i, 1), flat[i - 1]);
        }
    }

    @Test
    public void setValueWritesExactlyOneEntry() {
        LinSpace2D g = LinSpace2D.linspace2D(0.0, 1.0, 5, 0.0, 1.0, 7).allocate();
        g.setValue(3, 4, 42.0);
        int touched = 0;
        double[] flat = g.values();
        for (int k = 0; k < flat.length; ++k) {
            if (flat[k] != 0.0) {
                ++touched;
                assertEquals(42.0, flat[k], 0.0);
                assertEquals((4 - 1) * g.sizeX() + (3 - 1), k);
            }
        }
        assertEquals(1, touched);
    }

    // ----- the axes are the axes -------------------------------------------

    @Test
    public void theAxesAreReproducedBitForBit() {
        long lcg = 1234567L;
        int checked = 0;
        for (int t = 0; t < 400; ++t) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = ((t % 5) == 0) ? a : spread(lcg);
            lcg = next(lcg);
            int nx = 1 + (int) ((lcg >>> 1) % 12);
            lcg = next(lcg);
            double c = spread(lcg);
            lcg = next(lcg);
            double d = ((t % 7) == 0) ? c : spread(lcg);
            lcg = next(lcg);
            int ny = 1 + (int) ((lcg >>> 1) % 12);

            LinSpace ax = LinSpace.linspace(a, b, nx);
            LinSpace ay = LinSpace.linspace(c, d, ny);
            LinSpace2D g = LinSpace2D.over(ax, ay);
            assertEquals(nx, g.sizeX());
            assertEquals(ny, g.sizeY());
            for (int i = 1; i <= nx; ++i) {
                assertSameBits("pointX", ax.point(i), g.pointX(i));
                assertSameBits("x().point", ax.point(i), g.x().point(i));
                ++checked;
            }
            for (int j = 1; j <= ny; ++j) {
                assertSameBits("pointY", ay.point(j), g.pointY(j));
                assertSameBits("y().point", ay.point(j), g.y().point(j));
                ++checked;
            }
        }
        assertTrue("checked " + checked + " abscissas", checked > 4000);
    }

    @Test
    public void overDoesNotTakeTheValuesOfItsArguments() {
        LinSpace ax = LinSpace.linspace(0.0, 1.0, 4).allocate();
        ax.setValue(1, 99.0);
        LinSpace ay = LinSpace.linspace(0.0, 1.0, 3);
        LinSpace2D g = LinSpace2D.over(ax, ay);
        assertFalse(g.hasValues());
        assertFalse(g.x().hasValues());
        // and the argument is untouched
        assertEquals(99.0, ax.value(1), 0.0);
    }

    @Test
    public void theWideSpanAxesSurviveTheRoundTrip() {
        double m = Double.MAX_VALUE;
        double[][] spans = { { -1e308, 1e308 }, { 1e308, -1e308 }, { -m, m }, { m, -m }, { -m / 2.0, m } };
        for (double[] s : spans) {
            for (int n : new int[] { 2, 3, 5, 17, 129 }) {
                LinSpace ax = LinSpace.linspace(s[0], s[1], n);
                LinSpace2D g = LinSpace2D.over(ax, LinSpace.linspace(0.0, 1.0, 3));
                for (int i = 1; i <= n; ++i) {
                    assertSameBits("wide axis " + ax, ax.point(i), g.pointX(i));
                    assertTrue(ax + " point " + i + " is " + g.pointX(i), isFinite(g.pointX(i)));
                }
            }
        }
    }

    // ----- what the one-dimensional class was repaired for ------------------

    @Test
    public void aDegenerateGridKeepsEveryPointItWasAskedFor() {
        LinSpace2D g = LinSpace2D.linspace2D(2.0, 2.0, 3, 5.0, 5.0, 4);
        assertEquals(3, g.sizeX());
        assertEquals(4, g.sizeY());
        assertEquals(12, g.size());
        for (int i = 1; i <= 3; ++i) {
            assertEquals(2.0, g.pointX(i), 0.0);
        }
        for (int j = 1; j <= 4; ++j) {
            assertEquals(5.0, g.pointY(j), 0.0);
        }
        LinSpace2D f = g.eval((x, y) -> x + y);
        assertEquals(12, f.values().length);
        for (int k = 0; k < 12; ++k) {
            assertEquals(7.0, f.values()[k], 0.0);
        }
        // and it can be sliced without the collapse the 1D class used to have
        LinSpace2D s = f.slice(1, 2, 2, 4);
        assertEquals(2, s.sizeX());
        assertEquals(3, s.sizeY());
        assertEquals(6, s.values().length);
    }

    @Test
    public void aSpanWiderThanTheDoubleRangeStaysFinite() {
        LinSpace2D g = LinSpace2D.linspace2D(-1e308, 1e308, 5, -1e308, 1e308, 5);
        for (int i = 1; i <= 5; ++i) {
            assertTrue("pointX " + i + " is " + g.pointX(i), isFinite(g.pointX(i)));
            assertTrue("pointY " + i + " is " + g.pointY(i), isFinite(g.pointY(i)));
        }
        assertEquals(0.0, g.pointX(3), 0.0);
        assertEquals(5.0E307, g.pointX(4), 0.0);
        LinSpace2D f = g.eval((x, y) -> 1.0);
        for (int k = 0; k < f.values().length; ++k) {
            assertEquals(1.0, f.values()[k], 0.0);
        }
    }

    @Test
    public void descendingAxesAreEvaluatedInTheirOwnDirection() {
        LinSpace2D g = LinSpace2D.compute(10.0, 0.0, 5, 4.0, -4.0, 3, (x, y) -> x - y);
        for (int j = 1; j <= 3; ++j) {
            for (int i = 1; i <= 5; ++i) {
                assertSameBits("value(" + i + ", " + j + ")", g.pointX(i) - g.pointY(j), g.value(i, j));
            }
        }
        assertEquals(10.0, g.pointX(1), 0.0);
        assertEquals(0.0, g.pointX(5), 0.0);
        assertEquals(4.0, g.pointY(1), 0.0);
        assertEquals(-4.0, g.pointY(3), 0.0);
    }

    // ----- eval, slice, cuts, transpose -------------------------------------

    @Test
    public void evalSamplesTheFunctionOnceForEveryPoint() {
        LinSpace2D g = LinSpace2D.linspace2D(-3.0, 2.0, 7, 1.0, 9.0, 5);
        final int[] calls = new int[1];
        LinSpace2D f = g.eval((x, y) -> {
            ++calls[0];
            return 3.0 * x - y;
        });
        assertEquals(35, calls[0]);
        assertEquals(35, f.size());
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("eval", 3.0 * g.pointX(i) - g.pointY(j), f.value(i, j));
            }
        }
        assertFalse("eval must not give the receiver values", g.hasValues());
    }

    @Test
    public void sliceKeepsItsBlockAndStaysOnTheGrid() {
        LinSpace2D g = LinSpace2D.compute(0.0, 6.0, 7, -2.0, 2.0, 5, (x, y) -> 100.0 * x + y);
        int slices = 0;
        for (int fx = 1; fx <= g.sizeX(); ++fx) {
            for (int tx = fx; tx <= g.sizeX(); ++tx) {
                for (int fy = 1; fy <= g.sizeY(); ++fy) {
                    for (int ty = fy; ty <= g.sizeY(); ++ty) {
                        LinSpace2D s = g.slice(fx, tx, fy, ty);
                        assertEquals(1 + tx - fx, s.sizeX());
                        assertEquals(1 + ty - fy, s.sizeY());
                        assertEquals(s.sizeX() * s.sizeY(), s.values().length);
                        assertSameBits("slice x start", g.pointX(fx), s.pointX(1));
                        assertSameBits("slice x end", g.pointX(tx), s.pointX(s.sizeX()));
                        assertSameBits("slice y start", g.pointY(fy), s.pointY(1));
                        assertSameBits("slice y end", g.pointY(ty), s.pointY(s.sizeY()));
                        for (int j = 1; j <= s.sizeY(); ++j) {
                            for (int i = 1; i <= s.sizeX(); ++i) {
                                assertSameBits("slice value", g.value(fx + i - 1, fy + j - 1), s.value(i, j));
                            }
                        }
                        ++slices;
                    }
                }
            }
        }
        assertEquals(28 * 15, slices);
    }

    @Test
    public void sliceTakesItsPositionsInEitherOrder() {
        LinSpace2D g = LinSpace2D.compute(0.0, 6.0, 7, -2.0, 2.0, 5, (x, y) -> 100.0 * x + y);
        LinSpace2D a = g.slice(2, 5, 1, 4);
        LinSpace2D b = g.slice(5, 2, 4, 1);
        assertEquals(a.sizeX(), b.sizeX());
        assertEquals(a.sizeY(), b.sizeY());
        for (int k = 0; k < a.values().length; ++k) {
            assertSameBits("either order", a.values()[k], b.values()[k]);
        }
    }

    @Test
    public void theCutsAgreeWithTheGrid() {
        LinSpace2D g = LinSpace2D.compute(-1.0, 1.0, 6, 3.0, 8.0, 4, (x, y) -> x * y);
        for (int j = 1; j <= g.sizeY(); ++j) {
            LinSpace line = g.alongX(j);
            assertEquals(g.sizeX(), line.size());
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("alongX abscissa", g.pointX(i), line.point(i));
                assertSameBits("alongX value", g.value(i, j), line.value(i));
            }
        }
        for (int i = 1; i <= g.sizeX(); ++i) {
            LinSpace line = g.alongY(i);
            assertEquals(g.sizeY(), line.size());
            for (int j = 1; j <= g.sizeY(); ++j) {
                assertSameBits("alongY abscissa", g.pointY(j), line.point(j));
                assertSameBits("alongY value", g.value(i, j), line.value(j));
            }
        }
    }

    @Test
    public void transposeIsItsOwnInverse() {
        LinSpace2D g = LinSpace2D.compute(-1.0, 1.0, 6, 3.0, 8.0, 4, (x, y) -> x - 7.0 * y);
        LinSpace2D t = g.transpose();
        assertEquals(g.sizeX(), t.sizeY());
        assertEquals(g.sizeY(), t.sizeX());
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("transpose", g.value(i, j), t.value(j, i));
                assertSameBits("transpose abscissa x", g.pointX(i), t.pointY(i));
                assertSameBits("transpose abscissa y", g.pointY(j), t.pointX(j));
            }
        }
        LinSpace2D back = t.transpose();
        for (int k = 0; k < g.values().length; ++k) {
            assertSameBits("transpose twice", g.values()[k], back.values()[k]);
        }
    }

    // ----- the cell lookup --------------------------------------------------

    @Test
    public void everyCellBracketsThePointItWasAskedFor() {
        long lcg = 987654321L;
        int lookups = 0;
        for (int t = 0; t < 600; ++t) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = spread(lcg);
            lcg = next(lcg);
            int n = 2 + (int) ((lcg >>> 1) % 30);
            LinSpace2D g = LinSpace2D.linspace2D(a, b, n, 0.0, 1.0, 2);
            for (int k = 0; k < 10; ++k) {
                lcg = next(lcg);
                double q;
                if (k < 3) {
                    q = g.pointX(1 + (int) ((lcg >>> 1) % n));
                } else if (k == 3) {
                    q = g.pointX(1);
                } else if (k == 4) {
                    q = g.pointX(n);
                } else {
                    double u = (lcg >>> 11) * 0x1.0p-53;
                    q = g.pointX(1) + u * (g.pointX(n) - g.pointX(1));
                }
                if (!isFinite(q)) {
                    continue;
                }
                assertBrackets(g, q);
                ++lookups;
            }
        }
        assertTrue("only " + lookups + " lookups", lookups > 5000);
    }

    @Test
    public void aWideSpanIsLocatedWithoutASearch() {
        double m = Double.MAX_VALUE;
        double[][] spans = { { -1e308, 1e308 }, { 1e308, -1e308 }, { -m, m }, { m, -m }, { -m / 2.0, m } };
        int lookups = 0;
        for (double[] s : spans) {
            for (int n : new int[] { 2, 3, 5, 17, 129 }) {
                LinSpace2D g = LinSpace2D.linspace2D(s[0], s[1], n, 0.0, 1.0, 2);
                for (int pos = 1; pos <= n; ++pos) {
                    assertBrackets(g, g.pointX(pos));
                    ++lookups;
                }
            }
        }
        assertEquals(5 * (2 + 3 + 5 + 17 + 129), lookups);
    }

    @Test
    public void theCellLookupAgreesWithTheSecondVariable() {
        LinSpace2D g = LinSpace2D.linspace2D(0.0, 1.0, 2, -10.0, 30.0, 9);
        assertEquals(1, g.cellY(-10.0));
        assertEquals(1, g.cellY(-9.0));
        assertEquals(8, g.cellY(30.0));
        assertEquals(8, g.cellY(29.0));
        // outside the interval answers the first or the last cell
        assertEquals(1, g.cellY(-1e300));
        assertEquals(8, g.cellY(1e300));
        for (int j = 1; j <= 9; ++j) {
            double q = g.pointY(j);
            int cell = g.cellY(q);
            assertTrue("cell " + cell, cell >= 1 && cell <= 8);
            assertTrue(q + " not in cell " + cell, g.pointY(cell) <= q && q <= g.pointY(cell + 1));
        }
    }

    @Test
    public void aDegenerateAxisHasOneCellAndASinglePointAxisHasNone() {
        LinSpace2D degenerate = LinSpace2D.linspace2D(3.0, 3.0, 5, 0.0, 1.0, 2);
        assertEquals(0.0, degenerate.x().spacing(), 0.0);
        assertEquals(1, degenerate.cellX(3.0));
        assertEquals(1, degenerate.cellX(0.0));

        LinSpace2D single = LinSpace2D.linspace2D(2.0, 9.0, 1, 0.0, 1.0, 2);
        assertEquals(2.0, single.pointX(1), 0.0);
        expect(IllegalStateException.class, () -> single.cellX(2.0));
        expect(IllegalArgumentException.class, () -> single.cellY(Double.NaN));
        expect(IllegalArgumentException.class, () -> single.cellY(Double.POSITIVE_INFINITY));
    }

    // ----- an independent oracle -------------------------------------------

    @Test
    public void theGridFeedsAnInterpolatorThatReproducesIt() {
        LinSpace2D g = LinSpace2D.compute(-2.0, 3.0, 6, 0.5, 4.5, 5, (x, y) -> 2.0 * x - 3.0 * y + 1.0);
        BilinearInterpolator bi = new BilinearInterpolator(g.x().points(), g.y().points(), g.toArrays());
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertEquals("knot (" + i + ", " + j + ")", g.value(i, j),
                        bi.value(g.pointX(i), g.pointY(j)), 1e-12);
            }
        }
        // the interpolant of a bilinear function is that function, so the
        // grid also predicts the values between its own knots
        for (int j = 1; j < g.sizeY(); ++j) {
            for (int i = 1; i < g.sizeX(); ++i) {
                double px = 0.5 * (g.pointX(i) + g.pointX(i + 1));
                double py = 0.5 * (g.pointY(j) + g.pointY(j + 1));
                assertEquals("midpoint", 2.0 * px - 3.0 * py + 1.0, bi.value(px, py), 1e-12);
                // and the cell lookup names the cell the interpolator used
                assertEquals(i, g.cellX(px));
                assertEquals(j, g.cellY(py));
            }
        }
    }

    // ----- traversals -------------------------------------------------------

    @Test
    public void theTraversalsVisitEveryPointOnceInStorageOrder() {
        LinSpace2D g = LinSpace2D.compute(0.0, 3.0, 4, 10.0, 12.0, 3, (x, y) -> x * y);
        final int[] k = new int[1];
        final double[] xs = new double[12];
        final double[] ys = new double[12];
        g.forEachPoint((x, y) -> {
            xs[k[0]] = x;
            ys[k[0]] = y;
            ++k[0];
        });
        assertEquals(12, k[0]);

        final int[] k2 = new int[1];
        final double[] vs = new double[12];
        g.forEachValue((x, y, v) -> {
            assertSameBits("forEachValue x", xs[k2[0]], x);
            assertSameBits("forEachValue y", ys[k2[0]], y);
            vs[k2[0]] = v;
            ++k2[0];
        });
        assertEquals(12, k2[0]);

        int at = 0;
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("order x", g.pointX(i), xs[at]);
                assertSameBits("order y", g.pointY(j), ys[at]);
                assertSameBits("order value", g.value(i, j), vs[at]);
                ++at;
            }
        }
        assertEquals(12, at);
    }

    // ----- validation -------------------------------------------------------

    @Test
    public void aGridTooLargeForAnIntIsRejected() {
        expect(IllegalArgumentException.class, () -> LinSpace2D.linspace2D(0.0, 1.0, 100000, 0.0, 1.0, 100000));
        expect(IllegalArgumentException.class,
                () -> LinSpace2D.over(LinSpace.linspace(0.0, 1.0, Integer.MAX_VALUE),
                        LinSpace.linspace(0.0, 1.0, 3)));
        // right at the edge it is accepted and size() tells the truth
        LinSpace2D g = LinSpace2D.linspace2D(0.0, 1.0, 46340, 0.0, 1.0, 46340);
        assertEquals(46340 * 46340, g.size());
        assertFalse(g.hasValues());
    }

    @Test
    public void validationFollowsTheOneDimensionalCase() {
        expect(IllegalArgumentException.class, () -> LinSpace2D.linspace2D(0.0, 1.0, 0, 0.0, 1.0, 3));
        expect(IllegalArgumentException.class, () -> LinSpace2D.linspace2D(0.0, 1.0, 3, 0.0, 1.0, -1));
        expect(IllegalArgumentException.class, () -> LinSpace2D.linspace2D(Double.NaN, 1.0, 3, 0.0, 1.0, 3));
        expect(IllegalArgumentException.class,
                () -> LinSpace2D.linspace2D(0.0, 1.0, 3, 0.0, Double.POSITIVE_INFINITY, 3));

        LinSpace2D g = LinSpace2D.linspace2D(0.0, 1.0, 3, 0.0, 1.0, 4);
        expect(IndexOutOfBoundsException.class, () -> g.pointX(0));
        expect(IndexOutOfBoundsException.class, () -> g.pointX(4));
        expect(IndexOutOfBoundsException.class, () -> g.pointY(5));
        expect(NoSuchElementException.class, () -> g.values());
        expect(NoSuchElementException.class, () -> g.value(1, 1));
        expect(NoSuchElementException.class, () -> g.setValue(1, 1, 0.0));
        expect(NoSuchElementException.class, () -> g.alongX(1));
        expect(NoSuchElementException.class, () -> g.alongY(1));
        expect(NoSuchElementException.class, () -> g.toArrays());
        expect(NoSuchElementException.class, () -> g.forEachValue((x, y, v) -> {
        }));
        // the position is checked before the values are
        expect(IndexOutOfBoundsException.class, () -> g.value(9, 1));

        LinSpace2D h = g.allocate();
        assertTrue(h.hasValues());
        assertEquals(12, h.values().length);
        expect(IndexOutOfBoundsException.class, () -> h.slice(0, 2, 1, 2));
        expect(IndexOutOfBoundsException.class, () -> h.slice(1, 2, 1, 9));
    }

    // ----- centered grids ---------------------------------------------------

    @Test
    public void theCenteredGridsMatchTheOneDimensionalFactories() {
        for (int nx = 1; nx <= 6; ++nx) {
            for (int ny = 1; ny <= 6; ++ny) {
                double[][] data = new double[nx][ny];
                LinSpace2D gi = LinSpace2D.centeredIntIndexed(data);
                LinSpace2D gd = LinSpace2D.centeredDoubleIndexed(data);
                assertEquals(nx, gi.sizeX());
                assertEquals(ny, gi.sizeY());
                LinSpace xi = LinSpace.centeredIntIndexed(new double[nx]);
                LinSpace xd = LinSpace.centeredDoubleIndexed(new double[nx]);
                LinSpace yi = LinSpace.centeredIntIndexed(new double[ny]);
                LinSpace yd = LinSpace.centeredDoubleIndexed(new double[ny]);
                for (int i = 1; i <= nx; ++i) {
                    assertSameBits("int x", xi.point(i), gi.pointX(i));
                    assertSameBits("double x", xd.point(i), gd.pointX(i));
                }
                for (int k = 1; k <= ny; ++k) {
                    assertSameBits("int y", yi.point(k), gi.pointY(k));
                    assertSameBits("double y", yd.point(k), gd.pointY(k));
                }
            }
        }
        // the two conventions, spelled out: the integer-indexed grid puts one
        // more point below zero than above for an even length, the other one
        // stays symmetric
        LinSpace2D even = LinSpace2D.centeredIntIndexed(new double[4][2]);
        assertEquals(-2.0, even.pointX(1), 0.0);
        assertEquals(1.0, even.pointX(4), 0.0);
        assertEquals(-1.0, even.pointY(1), 0.0);
        assertEquals(0.0, even.pointY(2), 0.0);
        LinSpace2D evenD = LinSpace2D.centeredDoubleIndexed(new double[4][2]);
        assertEquals(-1.5, evenD.pointX(1), 0.0);
        assertEquals(1.5, evenD.pointX(4), 0.0);
        assertEquals(-0.5, evenD.pointY(1), 0.0);
        assertEquals(0.5, evenD.pointY(2), 0.0);
    }

    @Test
    public void theFirstIndexOfAKernelRunsAlongTheFirstVariable() {
        // the trap the convention sets: read row by row, this literal looks
        // like a derivative in y. It is one in x
        double[][] sobel = { { 1.0, 2.0, 1.0 }, { 0.0, 0.0, 0.0 }, { -1.0, -2.0, -1.0 } };
        LinSpace2D g = LinSpace2D.centeredIntIndexed(sobel);
        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                assertSameBits("(" + i + ", " + j + ")", sobel[i - 1][j - 1], g.value(i, j));
            }
        }
        // the cut along x at y == 0 carries the whole difference ...
        double[] acrossX = g.alongX(2).values();
        assertArrayEquals(new double[] { 2.0, 0.0, -2.0 }, acrossX, 0.0);
        // ... and the cut along y at x == 0 is flat
        double[] acrossY = g.alongY(2).values();
        assertArrayEquals(new double[] { 0.0, 0.0, 0.0 }, acrossY, 0.0);
        // the values are a copy
        sobel[0][0] = 99.0;
        assertEquals(1.0, g.value(1, 1), 0.0);
    }

    @Test
    public void aCenteredGridRoundTripsThroughToArrays() {
        // this grid is already centered on integers, so the abscissas have to
        // come back as well as the values
        LinSpace2D g = LinSpace2D.compute(-2.0, 2.0, 5, -1.0, 1.0, 3, (x, y) -> 10.0 * x + y);
        LinSpace2D back = LinSpace2D.centeredIntIndexed(g.toArrays());
        assertEquals(g.sizeX(), back.sizeX());
        assertEquals(g.sizeY(), back.sizeY());
        for (int j = 1; j <= g.sizeY(); ++j) {
            for (int i = 1; i <= g.sizeX(); ++i) {
                assertSameBits("round trip value", g.value(i, j), back.value(i, j));
            }
        }
        for (int i = 1; i <= g.sizeX(); ++i) {
            assertSameBits("round trip x", g.pointX(i), back.pointX(i));
        }
        for (int j = 1; j <= g.sizeY(); ++j) {
            assertSameBits("round trip y", g.pointY(j), back.pointY(j));
        }
    }

    @Test
    public void theCenteredFactoriesRejectAShapeTheyCannotUse() {
        expect(IllegalArgumentException.class, () -> LinSpace2D.centeredIntIndexed(new double[0][0]));
        expect(IllegalArgumentException.class, () -> LinSpace2D.centeredIntIndexed(new double[3][0]));
        expect(IllegalArgumentException.class, () -> LinSpace2D.centeredDoubleIndexed(new double[0][0]));
        // ragged is the one way of getting the shape wrong that 1D cannot have
        expect(IllegalArgumentException.class,
                () -> LinSpace2D.centeredIntIndexed(new double[][] { { 1.0, 2.0 }, { 3.0 } }));
        expect(IllegalArgumentException.class,
                () -> LinSpace2D.centeredDoubleIndexed(new double[][] { { 1.0 }, { 2.0, 3.0 } }));
        expect(NullPointerException.class, () -> LinSpace2D.centeredIntIndexed(null));
        expect(NullPointerException.class, () -> LinSpace2D.centeredIntIndexed(new double[][] { { 1.0 }, null }));
    }

    // ----- helpers ----------------------------------------------------------

    private static void assertBrackets(LinSpace2D g, double q) {
        int n = g.sizeX();
        boolean up = g.pointX(n) >= g.pointX(1);
        int i = g.cellX(q);
        assertTrue(g + " cell " + i, i >= 1 && i <= n - 1);
        double lo = g.pointX(i);
        double hi = g.pointX(i + 1);
        boolean ok = up ? (lo <= q && q <= hi) : (hi <= q && q <= lo);
        if (!ok) {
            // outside the interval the first or the last cell is the answer
            boolean outside = up ? (q < g.pointX(1) || q > g.pointX(n)) : (q > g.pointX(1) || q < g.pointX(n));
            if (!(outside && (i == 1 || i == n - 1))) {
                fail(g + ": " + q + " is not in cell " + i + " [" + lo + ", " + hi + "]");
            }
        }
    }

    private static boolean isFinite(double x) {
        return !Double.isInfinite(x) && !Double.isNaN(x);
    }

    private static long next(long lcg) {
        return lcg * 6364136223846793005L + 1442695040888963407L;
    }

    private static double spread(long lcg) {
        return ((lcg >>> 11) * 0x1.0p-53) * 2000.0 - 1000.0;
    }

    private static void assertSameBits(String what, double expected, double actual) {
        if (Double.doubleToRawLongBits(expected) != Double.doubleToRawLongBits(actual)) {
            fail(what + ": expected " + expected + " but was " + actual);
        }
    }

    private static void expect(Class<? extends RuntimeException> type, Runnable call) {
        try {
            call.run();
            fail("no " + type.getSimpleName());
        } catch (RuntimeException e) {
            if (!type.isInstance(e)) {
                fail("expected " + type.getSimpleName() + " but got " + e.getClass().getSimpleName() + ": "
                        + e.getMessage());
            }
        }
    }
}
