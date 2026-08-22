package math.demo.hubble;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

/**
 * The slopes are properties of 24 fixed numbers and are asserted to the digits
 * the demo prints. The two identities in section 3 are asserted as identities,
 * because that is what they are: the geometric mean regression is the ratio of
 * the two spreads, and total least squares on standardized columns is the
 * geometric mean regression.
 */
public class HubbleDemoTest {

    @Test
    public void dataIsIntact() {
        assertEquals(24, Datasets.size());
        assertEquals(21.873, Datasets.distanceChecksum(), 1.0e-9);
        assertEquals(8955.0, Datasets.velocityChecksum(), 1.0e-9);

        double[] r = Datasets.distance();
        double[] v = Datasets.velocity();
        double[] ms = Datasets.brightestStarMagnitude();
        double[] mt = Datasets.visualMagnitude();
        assertEquals(24, r.length);
        assertEquals(24, v.length);
        assertEquals(24, ms.length);
        assertEquals(24, mt.length);

        int withMagnitude = 0;
        for (int i = 0; i < r.length; ++i) {
            assertTrue("distance " + i, r[i] > 0.0 && r[i] <= 2.0);
            assertTrue("velocity " + i, v[i] >= -220.0 && v[i] <= 1090.0);
            assertTrue("visual magnitude " + i, mt[i] > 0.0 && mt[i] < 13.0);
            if (!Double.isNaN(ms[i])) {
                withMagnitude++;
                assertTrue("brightest star magnitude " + i, ms[i] > 15.0 && ms[i] < 21.0);
            }
            if (i > 0) {
                assertTrue("the paper orders the table by distance", r[i] >= r[i - 1]);
            }
            assertTrue("every object is named", Datasets.name(i).length() > 0);
        }
        assertEquals("the paper prints two dots for ten objects", 14, withMagnitude);
    }

    @Test
    public void accessorsHandOutCopies() {
        double[] first = Datasets.distance();
        assertNotSame(first, Datasets.distance());
        first[0] = -1.0;
        assertEquals(0.032, Datasets.distance()[0], 0.0);
    }

    /** Four estimators, four numbers, all of them properties of the same 24 points. */
    @Test
    public void theFourSlopes() {
        double[] s = HubbleDemo.slopes(HubbleDemo.everything());
        assertEquals("ordinary", 454.158, s[0], 1.0e-3);
        assertEquals("reverse", 728.366, s[1], 1.0e-3);
        assertEquals("total least squares", 728.365, s[2], 1.0e-3);
        assertEquals("geometric mean", 575.147, s[3], 1.0e-3);

        assertTrue("the ordinary fit is the shallowest", s[0] < s[3]);
        assertTrue("the geometric mean lies between", s[3] < s[1]);
        assertEquals("total least squares sits on the reverse regression here", s[1], s[2], 1.0e-2);
    }

    /**
     * The finding of section 3, as an assertion: the same estimator returns two
     * different answers depending on the units of the columns it is handed.
     */
    @Test
    public void totalLeastSquaresIsNotScaleInvariant() {
        double[] r = Datasets.distance();
        double[] v = Datasets.velocity();
        double sdR = HubbleDemo.standardDeviation(r);
        double sdV = HubbleDemo.standardDeviation(v);

        double raw = HubbleDemo.totalLeastSquaresBySvd(r, v);
        double scaled = HubbleDemo.totalLeastSquaresBySvd(HubbleDemo.standardized(r), HubbleDemo.standardized(v))
                * sdV / sdR;
        double[] s = HubbleDemo.slopes(HubbleDemo.everything());

        assertEquals("on the raw columns it is the reverse regression", s[1], raw, 1.0e-2);
        assertEquals("on scaled columns it is the geometric mean", s[3], scaled, 1.0e-6);
        assertTrue("and the two are far apart: " + raw + " against " + scaled, raw - scaled > 100.0);
    }

    /** The scale-invariant estimator is the ratio of the two spreads, exactly. */
    @Test
    public void theGeometricMeanIsTheRatioOfTheSpreads() {
        double ratio = HubbleDemo.standardDeviation(Datasets.velocity())
                / HubbleDemo.standardDeviation(Datasets.distance());
        assertEquals(ratio, HubbleDemo.slopes(HubbleDemo.everything())[3], 1.0e-9);
    }

    /** What the ordinary fit reports about itself, and the route it takes to get there. */
    @Test
    public void theOrdinaryFitReportsItself() {
        HubbleDemo.Ordinary o = HubbleDemo.ordinary();
        assertEquals(454.158, o.slope, 1.0e-3);
        assertEquals(-40.784, o.intercept, 1.0e-3);
        assertEquals(75.237, o.slopeStandardError, 1.0e-3);
        assertEquals(0.6235, o.rSquared, 1.0e-4);

        assertEquals("t is the slope over its standard error", o.slope / o.slopeStandardError, o.tValue, 1.0e-9);
        assertTrue("p = " + o.pValue, o.pValue < 1.0e-5);
        assertTrue("the slope must lie inside its own interval",
                o.slopeInterval[0] < o.slope && o.slope < o.slopeInterval[1]);
        assertEquals("the interval must be centred on the estimate", o.slope,
                0.5 * (o.slopeInterval[0] + o.slopeInterval[1]), 1.0e-9);
        assertEquals("SvdLeastSquares is the route OLS takes, so it must agree exactly", o.slope, o.slopeBySvd, 0.0);
    }

    /** A seeded bootstrap has to reproduce exactly, or the demo cannot print it. */
    @Test
    public void theBootstrapReproduces() {
        HubbleDemo.Interval first = HubbleDemo.interval(0);
        HubbleDemo.Interval second = HubbleDemo.interval(0);
        assertEquals(first.percentile[0], second.percentile[0], 0.0);
        assertEquals(first.percentile[1], second.percentile[1], 0.0);
        assertEquals(first.bca[0], second.bca[0], 0.0);
        assertEquals(first.bca[1], second.bca[1], 0.0);
        assertTrue("the interval must contain its own point estimate",
                first.percentile[0] < first.point && first.point < first.percentile[1]);
    }

    /** The claim the demo exists to make. */
    @Test
    public void noIntervalReachesTheModernValue() {
        HubbleDemo.Interval[] intervals = new HubbleDemo.Interval[4];
        for (int k = 0; k < intervals.length; ++k) {
            intervals[k] = HubbleDemo.interval(k);
            assertFalse(intervals[k].name + " covers " + Datasets.MODERN, intervals[k].covers(Datasets.MODERN));
            assertTrue(intervals[k].name + " must lie above it", intervals[k].percentile[0] > Datasets.MODERN);
            assertTrue(intervals[k].name + " BCa must lie above it too", intervals[k].bca[0] > Datasets.MODERN);
        }
        double nearest = HubbleDemo.nearest(intervals);
        assertTrue("the nearest bound is " + nearest, nearest > 4.0 * Datasets.MODERN);
    }

    /** Hubble's own value divides the four estimators, which is the difficulty stated as data. */
    @Test
    public void hubblesOwnValueIsCoveredByTwoOfTheFour() {
        int covering = 0;
        for (int k = 0; k < 4; ++k) {
            covering += HubbleDemo.interval(k).covers(Datasets.HUBBLE_K) ? 1 : 0;
        }
        assertEquals("the table does not decide between the estimators", 2, covering);
    }

    /** Everything the demo prints has to come out the same way twice. */
    @Test
    public void theDemoPrintsTheSameThingTwice() {
        String first = run();
        String second = run();
        assertEquals(first, second);
        assertTrue("the demo printed almost nothing", first.length() > 2000);
        assertFalse("a locale slipped into the output", first.contains("454,"));
        assertTrue("the demo must say what it established", first.contains("what this run established"));
        assertTrue("and it must name the four estimators", first.contains("geometric mean"));
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            HubbleDemo.main(new String[0]);
        } catch (UnsupportedEncodingException e) {
            throw new AssertionError(e);
        } finally {
            System.setOut(out);
        }
        try {
            return buffer.toString("UTF-8");
        } catch (UnsupportedEncodingException e) {
            throw new AssertionError(e);
        }
    }
}
