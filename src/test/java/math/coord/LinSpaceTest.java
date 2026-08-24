package math.coord;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.NoSuchElementException;

import org.junit.Test;

import math.fun.DIndexIterator;

/**
 * The private constructor used to fold {@code start == end} onto a single
 * point, throwing away the number of points the caller had asked for. That was
 * visible directly, as {@code linspace(2, 2, 100).size() == 1}, and indirectly
 * through {@code slice}, {@code sliceTo} and {@code sliceFrom}, which built
 * their endpoints by arithmetic: two endpoints that coincide by rounding folded
 * the slice onto one point as well, and threw
 * {@code IllegalStateException: inconsistent vector dimension : 2 != 2} when
 * the slice carried data. A {@code LinSpace} of one point also kept
 * {@code end}, so {@code linspace(0, 10, 1)} sat on 10 rather than on 0.
 * <p>
 * The abscissas then came from three different routes. {@code point(int)} and
 * {@code slice} multiplied, {@code points()}, {@code iterator()},
 * {@code forEach()} and {@code forEachBi()} accumulated a running sum, and
 * {@code eval} accumulated without pinning the last point onto {@code stop} --
 * with a guard written for ascending ranges only, so a descending
 * {@code compute} returned {@code f(stop)} at every point. Measured against the
 * exact form, the accumulating route was about 42 times the less accurate of
 * the two. All six accessors go through one private {@code abscissa(pos)} now.
 * <p>
 * The step was formed as {@code stop - start} before the division, so a span
 * beyond {@code Double.MAX_VALUE} -- reachable from about {@code +/-MAX/2} --
 * made the intermediate infinite and with it every abscissa but the two
 * endpoints. Halving the bounds is only half the repair: the product
 * {@code (pos - 1) * step} leaves the range a second time for a position far
 * from {@code start}, which is why those abscissas are measured from whichever
 * endpoint is the nearer one. Ascending and descending alike: for a symmetric
 * interval the descending grid is the ascending one read backwards, bit for
 * bit.
 * <p>
 * The tests for ordinary ranges are pins: they pass before and after, and
 * deliberately use grids whose points are exactly representable.
 */
public class LinSpaceTest {

    // ----- the collapse ------------------------------------------------

    @Test
    public void aDegenerateRangeKeepsTheRequestedCount() {
        LinSpace ls = LinSpace.linspace(2.0, 2.0, 100);
        assertEquals(100, ls.size());
        assertEquals(100, ls.points().length);
        assertEquals("1x100 :  [2.0  ...  2.0]", ls.toString());
        assertEquals(100, ls.allocate().values().length);
    }

    @Test
    public void everyPointOfADegenerateRangeIsTheEndpoint() {
        LinSpace ls = LinSpace.linspace(2.0, 2.0, 4);
        assertEquals(0.0, ls.spacing(), 0.0);
        assertArrayEquals(new double[] { 2.0, 2.0, 2.0, 2.0 }, ls.points(), 0.0);
        for (int pos = 1; pos <= 4; ++pos) {
            assertEquals(2.0, ls.point(pos), 0.0);
        }

        DIndexIterator it = ls.iterator();
        int seen = 0;
        while (it.hasNext()) {
            assertEquals(++seen, it.nextIndex());
            assertEquals(2.0, it.next(), 0.0);
        }
        assertEquals(4, seen);

        final int[] count = new int[1];
        ls.forEach().forEachRemaining(x -> {
            assertEquals(2.0, x, 0.0);
            ++count[0];
        });
        assertEquals(4, count[0]);
    }

    @Test
    public void computeSamplesTheFunctionOnceForEveryPoint() {
        final int[] calls = new int[1];
        LinSpace result = LinSpace.compute(2.0, 2.0, 4, x -> {
            ++calls[0];
            return x * x;
        });
        assertArrayEquals(new double[] { 4.0, 4.0, 4.0, 4.0 }, result.values(), 0.0);
        assertEquals(4, calls[0]);
    }

    @Test
    public void signedZeroEndpointsDoNotCollapse() {
        // -0.0 == 0.0 is true, so the collapse used to fire here as well
        LinSpace ls = LinSpace.linspace(-0.0, 0.0, 5);
        assertEquals(5, ls.size());
        assertEquals(0.0, ls.spacing(), 0.0);
        double[] points = ls.points();
        assertEquals(5, points.length);
        for (int i = 0; i < points.length; ++i) {
            assertEquals(0.0, points[i], 0.0);
        }
    }

    // ----- the slice trap the collapse caused --------------------------

    @Test
    public void aSliceWhoseEndpointsCoincideKeepsItsLength() {
        // eight points across one ulp of 100: the step is smaller than half an
        // ulp, so every derived endpoint rounds back onto the one before it
        LinSpace valued = LinSpace.linspace(100.0, 100.00000000000001, 8).allocate();
        assertEquals(2, valued.slice(1, 2).size());
        assertEquals(2, valued.sliceTo(2).size());
        assertEquals(4, valued.slice(1, 4).size());
        assertEquals(7, valued.sliceFrom(2).size());

        LinSpace bare = LinSpace.linspace(100.0, 100.00000000000001, 8);
        assertEquals(2, bare.slice(1, 2).size());
        assertArrayEquals(new double[] { 100.0, 100.0 }, bare.slice(1, 2).points(), 0.0);
    }

    @Test
    public void aDegenerateSliceKeepsItsValues() {
        LinSpace ls = LinSpace.linspace(100.0, 100.00000000000001, 8).allocate();
        for (int pos = 1; pos <= ls.size(); ++pos) {
            ls.setValue(pos, 10.0 * pos);
        }
        assertArrayEquals(new double[] { 20.0, 30.0, 40.0, 50.0 }, ls.slice(2, 5).values(), 0.0);
        assertArrayEquals(new double[] { 10.0, 20.0 }, ls.sliceTo(2).values(), 0.0);
        assertArrayEquals(new double[] { 70.0, 80.0 }, ls.sliceFrom(7).values(), 0.0);
    }

    @Test
    public void noSliceIsShorterThanAsked() {
        double[] starts = { 0.0, 1.0, 100.0, 1e16, -3.5, 1e-320, 2.0 };
        int checked = 0;
        for (int si = 0; si < starts.length; ++si) {
            double from = starts[si];
            double[] ends = { from, Math.nextUp(from), from + 4e-16, from + 1e-15, from + 1.0 };
            for (int ei = 0; ei < ends.length; ++ei) {
                for (int n = 2; n <= 12; ++n) {
                    LinSpace ls = LinSpace.linspace(from, ends[ei], n).allocate();
                    assertEquals(n, ls.size());
                    for (int lo = 1; lo <= n; ++lo) {
                        for (int hi = lo; hi <= n; ++hi) {
                            assertEquals(1 + hi - lo, ls.slice(lo, hi).size());
                            ++checked;
                        }
                    }
                    assertEquals(n, ls.sliceTo(n).size());
                    assertEquals(n, ls.sliceFrom(1).size());
                }
            }
        }
        assertEquals(12705, checked);
    }

    // ----- a single point sits on start --------------------------------

    @Test
    public void aSinglePointRangeSitsOnTheStart() {
        LinSpace ls = LinSpace.linspace(0.0, 10.0, 1);
        assertEquals(1, ls.size());
        assertEquals(0.0, ls.start(), 0.0);
        assertEquals(0.0, ls.end(), 0.0);
        assertEquals(0.0, ls.point(1), 0.0);
        assertEquals(0.0, ls.spacing(), 0.0);
        assertArrayEquals(new double[] { 0.0 }, ls.points(), 0.0);
        assertArrayEquals(new double[] { 0.0 }, LinSpace.compute(0.0, 10.0, 1, x -> x * x).values(), 0.0);

        // and the same the other way round
        assertEquals(10.0, LinSpace.linspace(10.0, 0.0, 1).point(1), 0.0);
        assertEquals(10.0, LinSpace.linspace(10.0, 0.0, 1).end(), 0.0);
    }

    @Test
    public void aSinglePointRangeAnswersTheSameThroughEveryAccessor() {
        LinSpace ls = LinSpace.linspace(-3.0, 7.0, 1);
        double expected = ls.point(1);
        assertEquals(expected, ls.start(), 0.0);
        assertEquals(expected, ls.end(), 0.0);
        assertEquals(expected, ls.points()[0], 0.0);
        assertEquals(expected, ls.iterator().next(), 0.0);
        final double[] first = new double[1];
        ls.forEach().tryAdvance(x -> first[0] = x);
        assertEquals(expected, first[0], 0.0);
    }

    // ----- pins: none of this may move ---------------------------------

    @Test
    public void ordinaryRangesAreUnchanged() {
        // exactly representable grids, so that this pin survives a later
        // change of the abscissa route
        LinSpace up = LinSpace.linspace(0.0, 4.0, 5);
        assertArrayEquals(new double[] { 0.0, 1.0, 2.0, 3.0, 4.0 }, up.points(), 0.0);
        assertEquals(1.0, up.spacing(), 0.0);

        LinSpace down = LinSpace.linspace(4.0, 0.0, 5);
        assertArrayEquals(new double[] { 4.0, 3.0, 2.0, 1.0, 0.0 }, down.points(), 0.0);
        assertEquals(1.0, down.spacing(), 0.0);

        for (int pos = 1; pos <= 5; ++pos) {
            assertEquals(pos - 1.0, up.point(pos), 0.0);
            assertEquals(5.0 - pos, down.point(pos), 0.0);
        }
        assertArrayEquals(new double[] { 0.0, 1.0, 4.0, 9.0, 16.0 },
                LinSpace.compute(0.0, 4.0, 5, x -> x * x).values(), 0.0);
    }

    @Test
    public void ordinarySlicesAreUnchanged() {
        LinSpace ls = LinSpace.linspace(0.0, 8.0, 5);
        assertArrayEquals(new double[] { 0.0, 2.0, 4.0, 6.0, 8.0 }, ls.points(), 0.0);
        assertArrayEquals(new double[] { 2.0, 4.0, 6.0 }, ls.slice(2, 4).points(), 0.0);
        assertArrayEquals(new double[] { 2.0, 4.0, 6.0 }, ls.slice(4, 2).points(), 0.0);
        assertArrayEquals(new double[] { 0.0, 2.0, 4.0 }, ls.sliceTo(3).points(), 0.0);
        assertArrayEquals(new double[] { 4.0, 6.0, 8.0 }, ls.sliceFrom(3).points(), 0.0);

        assertEquals(1, ls.slice(3, 3).size());
        assertEquals(4.0, ls.slice(3, 3).point(1), 0.0);
        assertEquals(1, ls.sliceTo(1).size());
        assertEquals(0.0, ls.sliceTo(1).point(1), 0.0);
        assertEquals(1, ls.sliceFrom(5).size());
        assertEquals(8.0, ls.sliceFrom(5).point(1), 0.0);
        assertEquals(5, ls.slice(1, 5).size());
    }

    @Test
    public void theCenteredFactoriesAreUnchanged() {
        assertArrayEquals(new double[] { -1.0, 0.0 }, LinSpace.centeredIntIndexed(new double[2]).points(), 0.0);
        assertArrayEquals(new double[] { -1.0, 0.0, 1.0 }, LinSpace.centeredIntIndexed(new double[3]).points(), 0.0);
        assertArrayEquals(new double[] { -2.0, -1.0, 0.0, 1.0 },
                LinSpace.centeredIntIndexed(new double[4]).points(), 0.0);
        assertArrayEquals(new double[] { -0.5, 0.5 }, LinSpace.centeredDoubleIndexed(new double[2]).points(), 0.0);
        assertArrayEquals(new double[] { -1.5, -0.5, 0.5, 1.5 },
                LinSpace.centeredDoubleIndexed(new double[4]).points(), 0.0);

        LinSpace one = LinSpace.centeredIntIndexed(new double[] { 7.0 });
        assertEquals(1, one.size());
        assertEquals(0.0, one.point(1), 0.0);
        assertEquals(7.0, one.value(1), 0.0);
    }

    @Test
    public void validationIsUnchanged() {
        expect(IllegalArgumentException.class, () -> LinSpace.linspace(0.0, 1.0, 0));
        expect(IllegalArgumentException.class, () -> LinSpace.linspace(0.0, 1.0, -3));
        expect(IllegalArgumentException.class, () -> LinSpace.linspace(Double.NaN, 1.0, 5));
        expect(IllegalArgumentException.class, () -> LinSpace.linspace(0.0, Double.POSITIVE_INFINITY, 5));
        expect(IllegalArgumentException.class, () -> LinSpace.centeredIntIndexed(new double[0]));
        expect(IllegalArgumentException.class, () -> LinSpace.centeredDoubleIndexed(new double[0]));
        expect(IndexOutOfBoundsException.class, () -> LinSpace.linspace(0.0, 1.0, 5).point(0));
        expect(IndexOutOfBoundsException.class, () -> LinSpace.linspace(0.0, 1.0, 5).point(6));
        expect(IndexOutOfBoundsException.class, () -> LinSpace.linspace(0.0, 1.0, 5).slice(0, 3));
        expect(NoSuchElementException.class, () -> LinSpace.linspace(0.0, 1.0, 5).value(1));
        expect(NoSuchElementException.class, () -> LinSpace.linspace(0.0, 1.0, 5).values());
    }

    // ----- one grid for every accessor ---------------------------------

    @Test
    public void aDescendingRangeIsEvaluatedInItsOwnDirection() {
        assertArrayEquals(new double[] { 16.0, 9.0, 4.0, 1.0, 0.0 },
                LinSpace.compute(4.0, 0.0, 5, x -> x * x).values(), 0.0);
        assertArrayEquals(new double[] { 1.0, 0.5, 0.0, -0.5, -1.0 },
                LinSpace.compute(1.0, -1.0, 5, x -> x).values(), 0.0);
        assertArrayEquals(new double[] { 8.0, 6.0, 4.0, 2.0, 0.0 },
                LinSpace.linspace(8.0, 0.0, 5).eval(x -> x).values(), 0.0);
        // through a slice of a descending range as well
        assertArrayEquals(new double[] { 6.0, 4.0, 2.0 },
                LinSpace.linspace(8.0, 0.0, 5).slice(2, 4).eval(x -> x).values(), 0.0);
    }

    @Test
    public void evalLandsOnBothEndpoints() {
        double[][] ranges = { { 0.0, 1.0 }, { 1.0, 0.0 }, { -3.5, 7.25 }, { 7.25, -3.5 } };
        for (int r = 0; r < ranges.length; ++r) {
            double a = ranges[r][0];
            double b = ranges[r][1];
            for (int n = 2; n <= 300; ++n) {
                double[] y = LinSpace.linspace(a, b, n).eval(x -> x).values();
                assertEquals("first point of linspace(" + a + ", " + b + ", " + n + ")", a, y[0], 0.0);
                assertEquals("last point of linspace(" + a + ", " + b + ", " + n + ")", b, y[n - 1], 0.0);
            }
        }
        assertEquals(1.0, LinSpace.compute(0.0, 1.0, 11, Math::sqrt).values()[10], 0.0);
    }

    @Test
    public void everyAccessorReturnsTheSameAbscissa() {
        long lcg = 987654321L;
        int checked = 0;
        for (int trial = 0; trial < 500; ++trial) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = spread(lcg);
            lcg = next(lcg);
            int n = 2 + (int) ((lcg >>> 11) % 64);
            LinSpace ls = LinSpace.linspace(a, b, n);

            double[] byArray = ls.points();
            double[] byEval = ls.eval(x -> x).values();
            final double[] byForEach = new double[n];
            final int[] k = new int[1];
            ls.forEach().forEachRemaining(x -> byForEach[k[0]++] = x);
            assertEquals(n, k[0]);

            DIndexIterator it = ls.iterator();
            for (int i = 0; i < n; ++i) {
                assertEquals(i + 1, it.nextIndex());
                double byIterator = it.next();
                double expected = ls.point(i + 1);
                assertSameBits("points()", expected, byArray[i]);
                assertSameBits("eval()", expected, byEval[i]);
                assertSameBits("forEach()", expected, byForEach[i]);
                assertSameBits("iterator()", expected, byIterator);
                ++checked;
            }
            assertFalse(it.hasNext());
        }
        assertTrue("expected a few thousand abscissas, got " + checked, checked > 10000);
    }

    @Test
    public void theGridIsMonotoneAndInsideTheInterval() {
        long lcg = 24680L;
        for (int trial = 0; trial < 500; ++trial) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = spread(lcg);
            lcg = next(lcg);
            int n = 2 + (int) ((lcg >>> 11) % 128);
            double lo = Math.min(a, b);
            double hi = Math.max(a, b);
            double[] p = LinSpace.linspace(a, b, n).points();
            assertEquals(a, p[0], 0.0);
            assertEquals(b, p[n - 1], 0.0);
            for (int i = 0; i < n; ++i) {
                assertTrue("point " + i + " of linspace(" + a + ", " + b + ", " + n + ") = " + p[i],
                        p[i] >= lo && p[i] <= hi);
                if (i > 0) {
                    assertTrue("not monotone at " + i, (b >= a) ? p[i] >= p[i - 1] : p[i] <= p[i - 1]);
                }
            }
        }
    }

    @Test
    public void sliceEndpointsLieOnTheParentGrid() {
        long lcg = 13579L;
        int checked = 0;
        for (int trial = 0; trial < 120; ++trial) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = spread(lcg);
            lcg = next(lcg);
            int n = 3 + (int) ((lcg >>> 11) % 20);
            LinSpace ls = LinSpace.linspace(a, b, n);
            for (int from = 1; from <= n; ++from) {
                for (int to = from; to <= n; ++to) {
                    LinSpace sub = ls.slice(from, to);
                    assertEquals(1 + to - from, sub.size());
                    assertSameBits("slice start", ls.point(from), sub.point(1));
                    assertSameBits("slice end", ls.point(to), sub.point(sub.size()));
                    ++checked;
                }
            }
        }
        assertTrue(checked > 10000);
    }

    @Test
    public void pointAndSpacingDidNotMove() {
        // the multiplying route is the one that stayed; this pins it against
        // the arithmetic as it was written before the accessors were unified
        long lcg = 112233L;
        for (int trial = 0; trial < 2000; ++trial) {
            lcg = next(lcg);
            double a = spread(lcg);
            lcg = next(lcg);
            double b = spread(lcg);
            lcg = next(lcg);
            int n = 1 + (int) ((lcg >>> 11) % 96);
            LinSpace ls = LinSpace.linspace(a, b, n);
            double aa = ls.start();
            double bb = ls.end();
            for (int pos = 1; pos <= ls.size(); ++pos) {
                assertSameBits("point(" + pos + ") of linspace(" + a + ", " + b + ", " + n + ")",
                        formerPoint(aa, bb, ls.size(), pos), ls.point(pos));
            }
            // by value, not bit for bit: the former spacing() could answer
            // -0.0 for linspace(0.0, -0.0, n), Math.abs(step) answers 0.0
            assertEquals(formerSpacing(aa, bb, ls.size()), ls.spacing(), 0.0);
        }
    }

    // ----- a span wider than the double range --------------------------

    @Test
    public void aSpanWiderThanTheDoubleRangeStaysFinite() {
        LinSpace up = LinSpace.linspace(-1e308, 1e308, 5);
        LinSpace down = LinSpace.linspace(1e308, -1e308, 5);
        assertArrayEquals(new double[] { -1e308, -1e308 / 2.0, 0.0, 1e308 / 2.0, 1e308 }, up.points(), 0.0);
        assertArrayEquals(new double[] { 1e308, 1e308 / 2.0, 0.0, -1e308 / 2.0, -1e308 }, down.points(), 0.0);
        assertEquals(1e308 / 2.0, up.spacing(), 0.0);
        assertEquals(1e308 / 2.0, down.spacing(), 0.0);
        assertArrayEquals(new double[] { -1e308, -1e308 / 2.0, 0.0, 1e308 / 2.0, 1e308 },
                up.eval(x -> x).values(), 0.0);
        assertArrayEquals(new double[] { 1e308, 1e308 / 2.0, 0.0, -1e308 / 2.0, -1e308 },
                down.eval(x -> x).values(), 0.0);

        // the widest range there is, either way round
        double max = Double.MAX_VALUE;
        assertArrayEquals(new double[] { -max, -max / 2.0, 0.0, max / 2.0, max },
                LinSpace.linspace(-max, max, 5).points(), 0.0);
        assertArrayEquals(new double[] { max, max / 2.0, 0.0, -max / 2.0, -max },
                LinSpace.linspace(max, -max, 5).points(), 0.0);

        // a count at which repairing the step alone would still overflow for
        // every position past the middle
        assertAllFinite("ascending, 1000 points", LinSpace.linspace(-1e308, 1e308, 1000).points());
        assertAllFinite("descending, 1000 points", LinSpace.linspace(1e308, -1e308, 1000).points());

        // and off-center, where the two halves are not the same width
        assertAllFinite("off-center ascending", LinSpace.linspace(-max, 1.0, 7).points());
        assertAllFinite("off-center descending", LinSpace.linspace(1.0, -max, 7).points());
        assertAllFinite("off-center ascending", LinSpace.linspace(-1e308, 1.7e308, 7).points());
        assertAllFinite("off-center descending", LinSpace.linspace(1.7e308, -1e308, 7).points());
    }

    @Test
    public void theWideGridIsSymmetricAndStrictlyOrdered() {
        int[] sizes = { 2, 3, 5, 17, 129, 1000, 1001 };
        for (int s = 0; s < sizes.length; ++s) {
            int n = sizes[s];
            double[] up = LinSpace.linspace(-1e308, 1e308, n).points();
            double[] down = LinSpace.linspace(1e308, -1e308, n).points();
            for (int i = 0; i < n; ++i) {
                assertTrue("n = " + n + ", ascending point " + i + " is " + up[i], isFinite(up[i]));
                assertTrue("n = " + n + ", descending point " + i + " is " + down[i], isFinite(down[i]));
                // the interval is symmetric about zero, so the grid must be
                assertEquals("n = " + n + ", symmetry at " + i, -up[n - 1 - i], up[i], 0.0);
                // and for a symmetric interval the descending grid is the
                // ascending one read backwards, exactly
                assertSameBits("n = " + n + ", mirror at " + i, up[n - 1 - i], down[i]);
                if (i > 0) {
                    assertTrue("n = " + n + ", not increasing at " + i, up[i] > up[i - 1]);
                    assertTrue("n = " + n + ", not decreasing at " + i, down[i] < down[i - 1]);
                }
            }
            if (n % 2 == 1) {
                assertEquals("n = " + n + ", ascending midpoint", 0.0, up[n / 2], 0.0);
                assertEquals("n = " + n + ", descending midpoint", 0.0, down[n / 2], 0.0);
            }
        }
    }

    @Test
    public void everyAccessorAgreesOnAWideSpan() {
        assertAccessorsAgree(LinSpace.linspace(-1e308, 1e308, 17));
        assertAccessorsAgree(LinSpace.linspace(1e308, -1e308, 17));
        assertAccessorsAgree(LinSpace.linspace(-Double.MAX_VALUE, 1.0, 17));
        assertAccessorsAgree(LinSpace.linspace(1.0, -Double.MAX_VALUE, 17));
    }

    @Test
    public void aWideSpanCanBeSliced() {
        assertSlicesStayOnTheGrid(LinSpace.linspace(-1e308, 1e308, 9));
        assertSlicesStayOnTheGrid(LinSpace.linspace(1e308, -1e308, 9));
        assertSlicesStayOnTheGrid(LinSpace.linspace(-Double.MAX_VALUE, 1.0, 9));
        assertSlicesStayOnTheGrid(LinSpace.linspace(1.0, -Double.MAX_VALUE, 9));
    }

    @Test
    public void aSubnormalSpacingIsNotHalvedAway() {
        // halving the bounds unconditionally would round this step down to
        // zero, which is why it only happens when the span is infinite
        LinSpace tiny = LinSpace.linspace(-Double.MIN_VALUE, Double.MIN_VALUE, 3);
        assertEquals(Double.MIN_VALUE, tiny.spacing(), 0.0);
        assertArrayEquals(new double[] { -Double.MIN_VALUE, 0.0, Double.MIN_VALUE }, tiny.points(), 0.0);

        LinSpace tinyDown = LinSpace.linspace(Double.MIN_VALUE, -Double.MIN_VALUE, 3);
        assertEquals(Double.MIN_VALUE, tinyDown.spacing(), 0.0);
        assertArrayEquals(new double[] { Double.MIN_VALUE, 0.0, -Double.MIN_VALUE }, tinyDown.points(), 0.0);

        LinSpace small = LinSpace.linspace(0.0, 8.0 * Double.MIN_VALUE, 5);
        assertEquals(2.0 * Double.MIN_VALUE, small.spacing(), 0.0);
        LinSpace smallDown = LinSpace.linspace(8.0 * Double.MIN_VALUE, 0.0, 5);
        assertEquals(2.0 * Double.MIN_VALUE, smallDown.spacing(), 0.0);
    }

    // ----- helpers -----------------------------------------------------

    private static boolean isFinite(double x) {
        return !Double.isInfinite(x) && !Double.isNaN(x);
    }

    private static void assertAllFinite(String what, double[] points) {
        for (int i = 0; i < points.length; ++i) {
            assertTrue(what + ", point " + i + " is " + points[i], isFinite(points[i]));
        }
    }

    private static void assertAccessorsAgree(LinSpace lsp) {
        int n = lsp.size();
        double[] byArray = lsp.points();
        double[] byEval = lsp.eval(x -> x).values();
        final double[] byForEach = new double[n];
        final int[] k = new int[1];
        lsp.forEach().forEachRemaining(x -> byForEach[k[0]++] = x);
        assertEquals(n, k[0]);

        DIndexIterator it = lsp.iterator();
        for (int i = 0; i < n; ++i) {
            assertEquals(i + 1, it.nextIndex());
            double byIterator = it.next();
            double expected = lsp.point(i + 1);
            assertTrue(lsp + ", point " + i + " is " + expected, isFinite(expected));
            assertSameBits(lsp + " points()", expected, byArray[i]);
            assertSameBits(lsp + " eval()", expected, byEval[i]);
            assertSameBits(lsp + " forEach()", expected, byForEach[i]);
            assertSameBits(lsp + " iterator()", expected, byIterator);
        }
        assertFalse(it.hasNext());
    }

    private static void assertSlicesStayOnTheGrid(LinSpace lsp) {
        int n = lsp.size();
        for (int from = 1; from <= n; ++from) {
            for (int to = from; to <= n; ++to) {
                LinSpace sub = lsp.slice(from, to);
                assertEquals(1 + to - from, sub.size());
                assertSameBits(lsp + " slice start", lsp.point(from), sub.point(1));
                assertSameBits(lsp + " slice end", lsp.point(to), sub.point(sub.size()));
                for (int pos = 1; pos <= sub.size(); ++pos) {
                    assertTrue(lsp + " slice(" + from + ", " + to + ").point(" + pos + ") is " + sub.point(pos),
                            isFinite(sub.point(pos)));
                }
            }
        }
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

    private static double formerSpacing(double a, double b, int n) {
        if (n == 1) {
            return 0.0;
        }
        if (b < a) {
            return (a - b) / (n - 1);
        }
        return (b - a) / (n - 1);
    }

    private static double formerStep(double a, double b, int n) {
        double spacing = formerSpacing(a, b, n);
        return (a > b) ? -spacing : spacing;
    }

    private static double formerPoint(double a, double b, int n, int pos) {
        if (pos == 1) {
            return a;
        }
        if (pos == n) {
            return b;
        }
        return a + ((pos - 1) * formerStep(a, b, n));
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
