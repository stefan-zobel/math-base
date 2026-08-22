package math.demo.sthelens;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

import math.fun.DBiFunction;

/**
 * The demo's claims, asserted against the 3721 fixed elevations they are
 * properties of.
 * <p>
 * Two of them are the point of the whole exercise and are asserted as
 * inequalities rather than to a tolerance, because it is their direction that
 * matters: the ground leaves the box of its four coarse samples more often than
 * the bicubic surface does, and where the bicubic surface leaves it, it is
 * nearer the truth than the bilinear one. If either ever flips, the demo's
 * conclusion has flipped with it and the text has to change.
 */
public class SurfaceDemoTest {

    @Test
    public void dataIsIntact() {
        assertEquals(61, Datasets.size());
        double[][] z = Datasets.elevation();
        assertEquals(61, z.length);
        double sum = 0.0;
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;
        for (int i = 0; i < z.length; ++i) {
            assertEquals("row " + i, 61, z[i].length);
            for (int j = 0; j < z[i].length; ++j) {
                sum += z[i][j];
                lo = Math.min(lo, z[i][j]);
                hi = Math.max(hi, z[i][j]);
            }
        }
        assertEquals("the checksum in the Datasets javadoc", 7509002.138, sum, 5.0e-3);
        assertEquals("lowest sample", 1387.69, lo, 5.0e-3);
        assertEquals("highest sample", 2530.95, hi, 5.0e-3);
    }

    @Test
    public void theGeometryIsConsistent() {
        double[] x = Datasets.eastings();
        double[] y = Datasets.northings();
        assertEquals(61, x.length);
        assertEquals(61, y.length);
        for (int i = 1; i < x.length; ++i) {
            assertTrue("eastings must increase", x[i] > x[i - 1]);
            assertTrue("northings must increase", y[i] > y[i - 1]);
        }
        // the window is square to within one per cent
        double east = x[x.length - 1];
        double north = y[y.length - 1];
        assertEquals("the window is not square", 1.0, east / north, 0.013);
        assertEquals(Datasets.LAT_SOUTH, Datasets.latitude(0), 0.0);
        assertEquals(Datasets.LAT_NORTH, Datasets.latitude(60), 1.0e-12);
        assertEquals(Datasets.LON_WEST, Datasets.longitude(0), 0.0);
        assertEquals(Datasets.LON_EAST, Datasets.longitude(60), 1.0e-12);
    }

    @Test
    public void accessorsHandOutCopies() {
        assertNotSame(Datasets.elevation(), Datasets.elevation());
        assertNotSame(Datasets.eastings(), Datasets.eastings());
        double[][] z = Datasets.elevation();
        double before = z[0][0];
        z[0][0] = -12345.0;
        assertEquals("the stored grid was reachable", before, Datasets.elevation()[0][0], 0.0);
    }

    /** The coarse survey is a subset of the lidar, not a resampling of it. */
    @Test
    public void theCoarseSurveyIsEveryFourthSample() {
        double[][] fine = Datasets.elevation();
        double[][] coarse = SurfaceDemo.coarse();
        assertEquals(16, SurfaceDemo.coarseSize());
        assertEquals(16, coarse.length);
        for (int i = 0; i < coarse.length; ++i) {
            for (int j = 0; j < coarse[i].length; ++j) {
                assertEquals("(" + i + ", " + j + ")", fine[i * SurfaceDemo.STRIDE][j * SurfaceDemo.STRIDE],
                        coarse[i][j], 0.0);
            }
        }
        assertEquals(61 * 61 - 16 * 16, SurfaceDemo.withheldPoints());
        assertEquals(3465, SurfaceDemo.withheldPoints());
    }

    /** Every one of them interpolates: at a sample it returns that sample. */
    @Test
    public void everySurfaceReproducesTheSamplesItWasGiven() {
        double[][] coarse = SurfaceDemo.coarse();
        double[] x = SurfaceDemo.coarseEastings();
        double[] y = SurfaceDemo.coarseNorthings();
        String[] names = SurfaceDemo.names();
        DBiFunction[] surfaces = SurfaceDemo.surfaces();
        for (int k = 0; k < surfaces.length; ++k) {
            for (int i = 0; i < x.length; ++i) {
                for (int j = 0; j < y.length; ++j) {
                    assertEquals(names[k] + " at (" + i + ", " + j + ")", coarse[i][j],
                            surfaces[k].apply(x[i], y[j]), 1.0e-6);
                }
            }
        }
    }

    /**
     * Bilinear is a convex combination of the four samples around a point, so
     * it cannot leave their range. This is the guarantee the demo is about.
     */
    @Test
    public void bilinearNeverLeavesTheBoxOfItsFourCorners() {
        SurfaceDemo.Score bilinear = SurfaceDemo.score(SurfaceDemo.surfaces()[0]);
        assertEquals("bilinear left the box", 0, bilinear.outside);
        assertEquals(0.0, bilinear.worstExcursion, 0.0);
    }

    /** And the mountain leaves it more than twice as often as bicubic does. */
    @Test
    public void theGroundLeavesThatBoxMoreOftenThanBicubicDoes() {
        int ground = SurfaceDemo.groundOutsideItsCorners();
        SurfaceDemo.Score bicubic = SurfaceDemo.score(SurfaceDemo.surfaces()[1]);
        assertTrue("bicubic never left the box, so the demo has no subject", bicubic.outside > 0);
        assertTrue("the ground stayed inside its own samples: " + ground + " against " + bicubic.outside,
                ground > bicubic.outside);
        assertEquals(238, ground);
        assertEquals(99, bicubic.outside);
    }

    /** Where it does leave the box, it is the nearer of the two, not the wilder. */
    @Test
    public void whereBicubicBreaksTheBoundItIsCloserThanBilinear() {
        DBiFunction[] s = SurfaceDemo.surfaces();
        double bicubic = SurfaceDemo.score(s[1]).rmsWhereOutside;
        double bilinear = SurfaceDemo.rmsWhere(s[0], s[1]);
        assertTrue("bicubic " + bicubic + " is not closer than bilinear " + bilinear,
                bicubic < bilinear);
        assertEquals(23.18, bicubic, 0.01);
        assertEquals(35.51, bilinear, 0.01);
    }

    /**
     * The correction the demo records about itself: the advantage is
     * concentrated at those points but not confined to them.
     */
    @Test
    public void theGainIsConcentratedThereButNotConfinedToIt() {
        double share = SurfaceDemo.shareOfGainAtExcursions();
        assertTrue("a fifth of the gain is not at the excursions: " + share, share > 0.10);
        assertTrue("nearly all of the gain is at the excursions: " + share, share < 0.35);
        assertEquals(0.196, share, 0.01);
    }

    /** Two entries of the menu, one method. */
    @Test
    public void bicubicAndTheNaturalSweepAreOneMethod() {
        double gap = SurfaceDemo.bicubicAgainstNatural();
        assertTrue("they are no longer the same surface, gap " + gap + " m", gap < 1.0e-9);
        SurfaceDemo.Score bicubic = SurfaceDemo.score(SurfaceDemo.surfaces()[1]);
        SurfaceDemo.Score natural = SurfaceDemo.score(SurfaceDemo.surfaces()[2]);
        assertEquals(bicubic.rms, natural.rms, 1.0e-9);
        assertEquals(bicubic.outside, natural.outside);
    }

    /** Kruger is the one that is bounded and still better than bilinear. */
    @Test
    public void krugerIsBoundedAndStillBeatsBilinear() {
        DBiFunction[] s = SurfaceDemo.surfaces();
        SurfaceDemo.Score kruger = SurfaceDemo.score(s[3]);
        SurfaceDemo.Score bilinear = SurfaceDemo.score(s[0]);
        assertEquals("Kruger left the box", 0, kruger.outside);
        assertTrue("Kruger " + kruger.rms + " no longer beats bilinear " + bilinear.rms,
                kruger.rms < bilinear.rms);
        // but it does not fix the worst case, which is the honest half of it
        assertTrue("Kruger now fixes the extreme too, which the demo denies",
                kruger.maxError > 0.9 * bilinear.maxError);
    }

    /** The ordering of the six, which is what the table in section 3 says. */
    @Test
    public void theTableIsInTheOrderTheDemoClaims() {
        DBiFunction[] s = SurfaceDemo.surfaces();
        double bilinear = SurfaceDemo.score(s[0]).rms;
        double bicubic = SurfaceDemo.score(s[1]).rms;
        double kruger = SurfaceDemo.score(s[3]).rms;
        double akima = SurfaceDemo.score(s[4]).rms;
        assertTrue("bicubic is not the most accurate", bicubic < kruger && bicubic < akima);
        assertTrue("bilinear is not the least accurate", bilinear > kruger && bilinear > akima);
        assertEquals(21.30, bilinear, 0.01);
        assertEquals(18.66, bicubic, 0.01);
    }

    @Test
    public void theDemoPrintsTheSameThingTwice() {
        String first = run();
        String second = run();
        assertEquals(first, second);
        assertTrue("the demo printed almost nothing", first.length() > 2000);
        assertFalse("a locale slipped into the output", first.contains("18,66"));
        assertTrue("the demo must say what it established", first.contains("what this run established"));
        assertTrue("and it must name the guarantee it is about", first.contains("convex combination"));
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            SurfaceDemo.main(new String[0]);
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
