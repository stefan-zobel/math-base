package math.demo.michelson;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

import math.distribution.StudentT;

/**
 * The bounds here were measured before they were asserted. The ones that are
 * properties of the table -- the mean, the counts, the checksum -- are exact;
 * the ones that are round-off are asserted as bounds, because their digits
 * differ between the two source trees and the demo does not print them.
 */
public class MichelsonDemoTest {

    @Test
    public void dataIsIntact() {
        assertEquals(100, Datasets.size());
        assertEquals(29985.24, Datasets.checksum(), 1.0e-3);
        assertEquals(24, Datasets.SETS);
        assertEquals(299.792458, Datasets.ACCEPTED, 0.0);

        double[] speed = Datasets.speed();
        double[] temperature = Datasets.temperature();
        double[] day = Datasets.day();
        double[] afternoon = Datasets.afternoon();
        int[] set = Datasets.set();
        assertEquals(100, speed.length);
        assertEquals(100, temperature.length);
        assertEquals(100, day.length);
        assertEquals(100, afternoon.length);
        assertEquals(100, set.length);

        int pm = 0;
        for (int i = 0; i < speed.length; ++i) {
            assertTrue("speed " + i + " is out of range: " + speed[i], speed[i] > 299.0 && speed[i] < 301.0);
            assertTrue("temperature " + i, temperature[i] >= 58.0 && temperature[i] <= 90.0);
            assertTrue("day " + i, day[i] >= 1.0 && day[i] <= 28.0);
            assertTrue("afternoon " + i, afternoon[i] == 0.0 || afternoon[i] == 1.0);
            assertTrue("set " + i, set[i] >= 1 && set[i] <= Datasets.SETS);
            pm += (int) afternoon[i];
        }
        assertEquals("afternoon runs", 63, pm);
        assertEquals("morning runs", 37, speed.length - pm);
    }

    @Test
    public void accessorsHandOutCopies() {
        double[] first = Datasets.speed();
        assertNotSame(first, Datasets.speed());
        double original = first[0];
        first[0] = -1.0;
        assertEquals(original, Datasets.speed()[0], 0.0);
    }

    /** The sample is a fact about 1879 and does not move. */
    @Test
    public void theSampleIsWhatItHasAlwaysBeen() {
        MichelsonDemo.Description d = MichelsonDemo.describe();
        assertEquals(100, d.runs);
        assertEquals(37, d.morning);
        assertEquals(63, d.afternoon);
        assertEquals(299.8524, d.mean, 1.0e-9);
        assertEquals(0.0790105478190557, d.sd, 1.0e-12);
        assertEquals(d.sd / 10.0, d.standardError, 1.0e-15);
        assertEquals(59.942, d.biasKmPerSecond(), 1.0e-3);
    }

    /**
     * The claim the demo exists to make: three intervals computed three
     * different ways agree with each other and none of them covers a value
     * that is known exactly.
     */
    @Test
    public void everyIntervalMissesTheAcceptedValue() {
        MichelsonDemo.Intervals in = MichelsonDemo.intervals();
        double[][] all = { in.classical, in.percentile, in.bca };
        String[] names = { "classical", "percentile", "BCa" };

        for (int k = 0; k < all.length; ++k) {
            assertFalse(names[k] + " covers the accepted value", MichelsonDemo.Intervals.contains(all[k]));
            assertTrue(names[k] + " must lie above it", all[k][0] > Datasets.ACCEPTED);
            assertTrue(names[k] + " misses by " + MichelsonDemo.Intervals.halfWidthsAway(all[k]) + " half widths",
                    MichelsonDemo.Intervals.halfWidthsAway(all[k]) > 2.0);
            assertTrue(names[k] + " must contain the sample mean",
                    all[k][0] <= in.mean && in.mean <= all[k][1]);
        }

        // and they agree with each other far more closely than they miss
        for (int k = 1; k < all.length; ++k) {
            assertEquals("lower bounds", all[0][0], all[k][0], 0.001);
            assertEquals("upper bounds", all[0][1], all[k][1], 0.001);
        }
    }

    /** A seeded bootstrap has to reproduce exactly, or the demo cannot print it. */
    @Test
    public void theBootstrapReproduces() {
        MichelsonDemo.Intervals first = MichelsonDemo.intervals();
        MichelsonDemo.Intervals second = MichelsonDemo.intervals();
        assertEquals(first.percentile[0], second.percentile[0], 0.0);
        assertEquals(first.percentile[1], second.percentile[1], 0.0);
        assertEquals(first.bca[0], second.bca[0], 0.0);
        assertEquals(first.bca[1], second.bca[1], 0.0);
        assertEquals(first.bootstrapMean, second.bootstrapMean, 0.0);
        assertEquals("the bootstrap mean must sit on the sample mean", first.mean, first.bootstrapMean, 0.001);
    }

    /** The critical value the interval is built from is the one every table prints. */
    @Test
    public void theCriticalValueIsTheTabulatedOne() {
        MichelsonDemo.Intervals in = MichelsonDemo.intervals();
        assertEquals(1.98422, in.criticalT, 1.0e-5);
        assertEquals(in.criticalT, new StudentT(99.0).inverseCdf(0.975), 0.0);
    }

    /**
     * Section 4: the density, the distribution function and the quantile are
     * separate pieces of code and have to agree. This is the assertion behind
     * the bound the demo prints instead of round-off digits.
     */
    @Test
    public void theLibraryAgreesWithItself() {
        MichelsonDemo.CrossCheck c = MichelsonDemo.crossCheck();
        assertBounded("cdf against the integral of the pdf", c.cdfAgainstIntegral);
        assertBounded("Normal.inverseCdf against an independent root", c.normalQuantileAgainstRoot);
        assertBounded("StudentT.inverseCdf against an independent root", c.studentQuantileAgainstRoot);
        assertBounded("the mean as an integral", c.meanAgainstIntegral);
        assertBounded("the variance as an integral", c.varianceAgainstIntegral);
        assertEquals(1.98422, c.tabulatedQuantile, 1.0e-5);
    }

    private static void assertBounded(String what, double deviation) {
        assertTrue(what + " is not a number", !Double.isNaN(deviation));
        assertTrue(what + " deviates by " + deviation, deviation <= MichelsonDemo.AGREEMENT);
    }

    /** On a hundred values a sketch stores everything, and still answers by its own rule. */
    @Test
    public void theSketchesAgreeWithTheOrderStatistics() {
        MichelsonDemo.SketchCheck s = MichelsonDemo.sketches();
        assertEquals("a compression of 100 cannot compress 100 points", 100, s.centroids);
        assertEquals("the median is the one quantile both conventions share", s.exact[2], s.digest[2], 1.0e-12);
        assertTrue("t-digest deviates by " + s.digestWorst, s.digestWorst < 1.0e-2);
        assertTrue("P-squared deviates by " + s.medianDeviation, s.medianDeviation < 1.0e-2);
        for (int k = 1; k < s.exact.length; ++k) {
            assertTrue("the exact quantiles must increase", s.exact[k] >= s.exact[k - 1]);
            assertTrue("so must the sketch's", s.digest[k] >= s.digest[k - 1]);
        }
    }

    /** The difference the eye sees does not survive a model that knows when the runs happened. */
    @Test
    public void theAfternoonEffectDoesNotSurviveTheModel() {
        assertTrue("the raw difference is real enough", MichelsonDemo.rawAfternoonEffect() > 40.0);

        MichelsonDemo.Model m = MichelsonDemo.model();
        int afternoon = m.indexOf("afternoon");
        assertTrue("p = " + m.pValues[afternoon], m.pValues[afternoon] > 0.1);
        assertTrue("the interval must cover zero",
                m.intervals[afternoon][0] < 0.0 && m.intervals[afternoon][1] > 0.0);
    }

    /** What does survive is the drift over the month. */
    @Test
    public void theDriftSurvivesTheModel() {
        MichelsonDemo.Model m = MichelsonDemo.model();
        int day = m.indexOf("day");
        assertTrue("the drift is downward: " + m.beta[day], m.beta[day] < 0.0);
        assertTrue("p = " + m.pValues[day], m.pValues[day] < 0.001);
        assertTrue("the interval must exclude zero", m.intervals[day][1] < 0.0);
        assertEquals("about four km/s per day", -4.4, 1000.0 * m.beta[day], 0.5);
    }

    /** The model explains some of the scatter and none of the bias. */
    @Test
    public void theModelExplainsScatterButNotTheBias() {
        MichelsonDemo.Model m = MichelsonDemo.model();
        assertTrue("R^2 = " + m.rSquared, m.rSquared > 0.15 && m.rSquared < 0.30);
        assertTrue("the residual spread must be smaller than the raw one", m.residualSd < m.rawSd);
        assertEquals(96, m.degreesOfFreedom);

        // the fit passes through the mean of the sample, so it cannot move the bias
        MichelsonDemo.Description d = MichelsonDemo.describe();
        assertTrue("the bias is what it was", Math.abs(d.biasKmPerSecond()) > 50.0);
    }

    /** Everything the demo prints has to come out the same way twice. */
    @Test
    public void theDemoPrintsTheSameThingTwice() {
        String first = run();
        String second = run();
        assertEquals(first, second);
        assertTrue("the demo printed almost nothing", first.length() > 2000);
        assertFalse("a locale slipped into the output", first.contains("299,8"));
        assertTrue("the demo must say what it established", first.contains("what this run established"));
        assertTrue("and it must name the value it is judged against", first.contains("299.792458"));
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            MichelsonDemo.main(new String[0]);
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
