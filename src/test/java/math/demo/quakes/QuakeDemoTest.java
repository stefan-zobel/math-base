package math.demo.quakes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.Map;

import org.junit.Test;

/**
 * The b values and the counts are properties of 2115 fixed numbers and are
 * asserted to the digits the demo prints. The agreements between two routes to
 * the same quantity are asserted as bounds, because their digits are round-off
 * and differ between the two source trees.
 */
public class QuakeDemoTest {

    @Test
    public void dataIsIntact() {
        assertEquals(2115, Datasets.size());
        assertEquals(11232.5, Datasets.magnitudeChecksum(), 1.0e-6);
        assertEquals(95169.969, Datasets.depthChecksum(), 1.0e-6);
        assertEquals(148, Datasets.regionNames().length);
        assertEquals(4, Datasets.scaleNames().length);
        assertEquals("mww", Datasets.MOMENT_MAGNITUDE);

        double[] magnitude = Datasets.magnitude();
        double[] depth = Datasets.depth();
        assertEquals(2115, magnitude.length);
        assertEquals(2115, depth.length);
        for (int i = 0; i < magnitude.length; ++i) {
            assertTrue("magnitude " + i + " is " + magnitude[i], magnitude[i] >= 5.0 && magnitude[i] <= 8.8);
            assertTrue("depth " + i + " is " + depth[i], depth[i] > 0.0 && depth[i] < 700.0);
            assertTrue("scale " + i, Datasets.scale(i).length() >= 2);
            assertTrue("region " + i, Datasets.region(i).length() > 0);
        }
        assertEquals("moment magnitudes", 1156, Datasets.magnitudesOn("mww").length);
        assertEquals("body wave magnitudes", 949, Datasets.magnitudesOn("mb").length);
    }

    @Test
    public void accessorsHandOutCopies() {
        double[] first = Datasets.magnitude();
        assertNotSame(first, Datasets.magnitude());
        first[0] = -1.0;
        assertEquals(5.5, Datasets.magnitude()[0], 0.0);
        String[] names = Datasets.regionNames();
        names[0] = "nowhere";
        assertEquals("Russia", Datasets.regionNames()[0]);
    }

    /** The finding of section 2: one column, two scales, three answers. */
    @Test
    public void theMixtureIsNotTheMomentMagnitude() {
        double[] mixed = QuakeDemo.above(5.0, null);
        double[] moment = Datasets.magnitudesOn("mww");
        double[] body = Datasets.magnitudesOn("mb");

        assertEquals("read as one quantity", 1.2034, QuakeDemo.akiB(mixed, 5.0), 1.0e-4);
        assertEquals("moment magnitude alone", 0.8684, QuakeDemo.akiB(moment, 5.0), 1.0e-4);
        assertEquals("body wave magnitude alone", 2.2967, QuakeDemo.akiB(body, 5.0), 1.0e-4);

        assertEquals("the body wave scale saturates", 6.3, QuakeDemo.max(body), 1.0e-9);
        assertEquals("the moment scale does not", 8.8, QuakeDemo.max(moment), 1.0e-9);
        assertTrue("the mixture looks plausible, which is the trap",
                Math.abs(QuakeDemo.akiB(mixed, 5.0) - 1.0) < 0.25);
    }

    /** The bias falls away as the threshold rises, and the standard error grows to meet it. */
    @Test
    public void bSettlesOnceTheCatalogueIsComplete() {
        double previousError = 0.0;
        for (int k = 0; k < QuakeDemo.THRESHOLDS.length; ++k) {
            double mc = QuakeDemo.THRESHOLDS[k];
            double[] m = QuakeDemo.above(mc, Datasets.MOMENT_MAGNITUDE);
            double b = QuakeDemo.akiB(m, mc);
            double error = b / Math.sqrt(m.length);

            if (mc <= 5.05) {
                assertTrue("at Mc " + mc + " incompleteness must still bias b low: " + b, b < 0.9);
            }
            if (mc >= 5.4) {
                assertTrue("at Mc " + mc + " b must have settled near one: " + b, b > 0.92 && b < 1.05);
            }
            assertTrue("the standard error must grow with the threshold", error > previousError);
            previousError = error;
        }
    }

    /** Two independent routes to the same number. */
    @Test
    public void bothRoutesAgreeOnB() {
        QuakeDemo.Fit fit = QuakeDemo.fit();
        assertEquals(437, fit.events);
        assertEquals(0.9862, fit.closedForm, 1.0e-4);
        assertTrue("the optimizer deviates by " + fit.agreement, fit.agreement <= QuakeDemo.AGREEMENT);
        assertTrue("it should not need many iterations: " + fit.iterations, fit.iterations < 50);
        assertTrue("the search must report a convergence, whichever rule fired: " + fit.termination,
                fit.termination.equals("VALUE_TOLERANCE") || fit.termination.equals("GRADIENT_TOLERANCE")
                        || fit.termination.equals("ZERO_GRADIENT"));
        assertTrue("the gradient must be small at the exit: " + fit.gradientNorm, fit.gradientNorm < 1.0e-6);
    }

    /** The correct comparison for quantized data reverses the naive one. */
    @Test
    public void theTailIsExponentialOnTheBinsTheDataHave() {
        QuakeDemo.Models models = QuakeDemo.models();
        assertEquals("exponential", models.names[models.best]);
        for (int k = 0; k < models.names.length; ++k) {
            if (k != models.best) {
                assertTrue(models.names[k] + " must not beat the exponential on AIC",
                        models.aic[k] > models.aic[models.best]);
            }
        }

        int bestContinuous = 0;
        for (int k = 1; k < models.continuous.length; ++k) {
            if (models.continuous[k] > models.continuous[bestContinuous]) {
                bestContinuous = k;
            }
        }
        assertEquals("scored as points, the lognormal wins on an artefact of the rounding", "lognormal",
                models.names[bestContinuous]);
    }

    /** What the two quantile sketches do with a column that is pinned to one value. */
    @Test
    public void theQuantileSketchesOnTheDepths() {
        QuakeDemo.SketchCheck s = QuakeDemo.sketches(100.0);

        assertTrue("t-digest deviates by " + s.digestWorst + " km", s.digestWorst < 1.0);
        assertEquals("and it gets the median exactly", s.exact[2], s.digest[2], 0.0);
        assertEquals("which is the depth the USGS assigns when it cannot resolve one",
                Datasets.UNRESOLVED_DEPTH, s.exact[2], 0.0);
        assertTrue("a large part of the catalogue carries that depth: " + s.unresolved,
                s.unresolved > Datasets.size() / 3);

        assertTrue("P-squared misses the median by " + s.squaredMiss + " km, which is the point",
                s.squaredMiss > 1.0);
        assertTrue("but not by more than that", s.squaredMiss < 5.0);
        assertTrue("the depths are not continuous either: " + s.distinctDepths,
                s.distinctDepths < Datasets.size() / 2);
    }

    /** The only guarantee a Count-Min sketch makes is that it never undercounts. */
    @Test
    public void theRegionSketchNeverUndercounts() {
        QuakeDemo.RegionSketch sketch = QuakeDemo.regionSketch(4, 128);
        assertEquals(148, sketch.keys);
        assertEquals(0, sketch.undercounts);
        assertTrue("worst overestimate " + sketch.worstOverestimate, sketch.worstOverestimate > 0);
        assertTrue("but a small one", sketch.worstOverestimate < 50);
        assertTrue("the mean overestimate is " + sketch.meanOverestimate, sketch.meanOverestimate < 5.0);

        assertEquals("the heaviest hitter is found", "Russia", sketch.topTen.get(0));
        for (int i = 0; i < sketch.estimated.length; ++i) {
            assertTrue(sketch.topTen.get(i) + " was undercounted", sketch.estimated[i] >= sketch.exact[i]);
        }

        Map<String, Integer> exact = QuakeDemo.exactRegionCounts();
        assertEquals(148, exact.size());
        assertEquals(508, exact.get("Russia").intValue());

        // a wider sketch is exact on a key space this small
        QuakeDemo.RegionSketch wide = QuakeDemo.regionSketch(4, 512);
        assertEquals("at four times the width nothing collides", 0L, wide.worstOverestimate);
    }

    /** A seeded bootstrap has to reproduce exactly, or the demo cannot print it. */
    @Test
    public void theBootstrapReproduces() {
        QuakeDemo.Fit first = QuakeDemo.fit();
        QuakeDemo.Fit second = QuakeDemo.fit();
        assertEquals(first.percentile[0], second.percentile[0], 0.0);
        assertEquals(first.percentile[1], second.percentile[1], 0.0);
        assertEquals(first.bca[0], second.bca[0], 0.0);
        assertEquals(first.bca[1], second.bca[1], 0.0);
        assertTrue("the interval must contain the estimate",
                first.percentile[0] < first.closedForm && first.closedForm < first.percentile[1]);
        assertTrue("and it must be narrow enough to be worth printing",
                first.percentile[1] - first.percentile[0] < 0.3);
    }

    /** The extrapolation is arithmetically sound and physically impossible. */
    @Test
    public void theReturnLevelIsImpossible() {
        QuakeDemo.Fit fit = QuakeDemo.fit();
        QuakeDemo.ReturnLevel level = QuakeDemo.returnLevel(fit);

        assertTrue("the two routes deviate by " + level.agreement, level.agreement <= QuakeDemo.AGREEMENT);
        assertEquals(437.0, level.rate, 0.0);
        assertEquals(8.8, level.largestObserved, 1.0e-9);

        assertTrue("the hundred year magnitude is " + level.byClosedForm,
                level.byClosedForm > QuakeDemo.STRONGEST_EVER);
        assertTrue("even the lower end of the interval is impossible: " + level.lower,
                level.lower > QuakeDemo.STRONGEST_EVER);
        assertTrue("the interval has width", level.upper > level.lower);
    }

    /** Everything the demo prints has to come out the same way twice. */
    @Test
    public void theDemoPrintsTheSameThingTwice() {
        String first = run();
        String second = run();
        assertEquals(first, second);
        assertTrue("the demo printed almost nothing", first.length() > 3000);
        assertFalse("a locale slipped into the output", first.contains("0,9"));
        assertTrue("the demo must say what it established", first.contains("what this run established"));
        assertTrue("and it must name the scale that saturates", first.contains("mb"));
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            QuakeDemo.main(new String[0]);
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
