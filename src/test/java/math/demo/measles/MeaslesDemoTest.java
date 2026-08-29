package math.demo.measles;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

import math.probe.ACF;

/**
 * What the demo claims is either a property of four fixed tables or a property
 * of two systems of differential equations, and both kinds are asserted here as
 * what they are rather than as digits copied off a page.
 * <p>
 * The sharpest of them is a negative: the mean interval between epidemic peaks
 * is the same for the annual attractor and the biennial one, so the statistic
 * that looks like it measures the period cannot. That is asserted directly,
 * because it is the reason the demo classifies on peak heights.
 */
public class MeaslesDemoTest {

    /** How close two doubles have to be to count as the same peak height. */
    private static final double SAME = 1.0e-3;

    @Test
    public void dataIsIntact() {
        double[] ma = Datasets.massachusetts();
        double[] ny = Datasets.newYork();
        double[] ct = Datasets.connecticut();
        double[] mt = Datasets.montana();

        assertEquals(2184, ma.length);
        assertEquals(1820, ny.length);
        assertEquals(1820, ct.length);
        assertEquals(1820, mt.length);
        assertEquals(Datasets.COMPARISON_WEEKS, ny.length);

        assertEquals(17663.05, Datasets.checksum(ma), 1.0e-9);
        assertEquals(12018.60, Datasets.checksum(ny), 1.0e-9);
        assertEquals(17517.17, Datasets.checksum(ct), 1.0e-9);
        assertEquals(21503.96, Datasets.checksum(mt), 1.0e-9);

        assertEquals(75, Datasets.weeksWithNoReport(ma));
        assertEquals(62, Datasets.weeksWithNoReport(ny));
        assertEquals(75, Datasets.weeksWithNoReport(ct));
        assertEquals(130, Datasets.weeksWithNoReport(mt));

        // an incidence is not negative and not a percentage of the population
        for (double[] series : new double[][] { ma, ny, ct, mt }) {
            for (int i = 0; i < series.length; ++i) {
                if (!Double.isNaN(series[i])) {
                    assertTrue("week " + i, series[i] >= 0.0 && series[i] < 1000.0);
                }
            }
        }

        assertEquals(1928, Datasets.yearOf(0));
        assertEquals(1928, Datasets.yearOf(51));
        assertEquals(1929, Datasets.yearOf(52));
        assertEquals(Datasets.MASSACHUSETTS_LAST_YEAR, Datasets.yearOf(ma.length - 1));
        assertEquals(Datasets.COMPARISON_LAST_YEAR, Datasets.yearOf(ny.length - 1));

        assertNotSame("the accessor hands out a copy", Datasets.massachusetts(),
                Datasets.massachusetts());
    }

    /**
     * A week with no report is not a week with no cases. The source prints a
     * dash and the array carries {@code NaN}, and nothing in the demo may turn
     * one into a zero -- which would be a claim rather than a gap.
     */
    @Test
    public void aMissingWeekIsNotAZero() {
        double[] ma = Datasets.massachusetts();
        int zeros = 0;
        for (int i = 0; i < ma.length; ++i) {
            if (ma[i] == 0.0) {
                ++zeros;
            }
        }
        assertTrue("the record does contain reported zeros", zeros > 0);
        assertTrue("and separately, weeks with no report", Datasets.weeksWithNoReport(ma) > 0);
    }

    // ------------------------------------------------------------ the models

    /**
     * The finding the demo was built to test, and it fails: the model without a
     * latent class period-doubles too, so a two year cycle is evidence of
     * nothing.
     */
    @Test
    public void bothModelsPeriodDouble() {
        assertEquals("SIR is annual at a mild forcing", 1, periodOfSir(0.02));
        assertEquals("SEIR is annual at a mild forcing", 1, periodOfSeir(0.02));
        assertEquals("SIR doubles too", 2, periodOfSir(0.08));
        assertEquals("SEIR doubles", 2, periodOfSeir(0.08));
    }

    /**
     * Why that took measuring. A biennial attractor keeps its peaks locked to
     * the forcing, so they stay a year apart and only their size alternates.
     */
    @Test
    public void theMeanIntervalCannotTellTheTwoApart() {
        MeaslesDemo.Attractor annual = sir(0.02);
        MeaslesDemo.Attractor biennial = sir(0.08);
        assertEquals(1, annual.period);
        assertEquals(2, biennial.period);

        assertEquals("the annual attractor steps one year", 1.0, annual.meanGap, 1.0e-3);
        assertEquals("and so does the biennial one", 1.0, biennial.meanGap, 1.0e-2);

        assertTrue("while the peak heights are what differ: " + biennial.largestPeak + " against "
                + biennial.smallestPeak, biennial.largestPeak > 4.0 * biennial.smallestPeak);
        assertEquals("and do not, in the annual one", annual.largestPeak, annual.smallestPeak,
                SAME * annual.largestPeak);
    }

    /**
     * What the latent period does decide: the forcing at which the model stops
     * having a periodic answer at all. The SEIR gives up first.
     */
    @Test
    public void theLatentPeriodMovesWhereTheModelBreaks() {
        assertTrue("the SEIR is irregular at 0.20", periodOfSeir(0.20) == 0);
        assertTrue("where the SIR is still periodic", periodOfSir(0.20) != 0);
    }

    /**
     * A trough is a property of the attractor while the attractor is a cycle and
     * a property of the trajectory once it is not. Two solver tolerances one
     * decade apart say which.
     */
    @Test
    public void theTroughIsTheModelsOnlyWhilePeriodic() {
        MeaslesDemo.Attractor periodic = seir(0.15);
        double tight = MeaslesDemo.troughAtTolerance(0.15, 1.0e-13, 1.0e-20);
        assertEquals(2, periodic.period);
        assertEquals("a periodic trough does not move with the tolerance", periodic.trough, tight,
                SAME * periodic.trough);

        MeaslesDemo.Attractor irregular = seir(0.20);
        double tightChaos = MeaslesDemo.troughAtTolerance(0.20, 1.0e-13, 1.0e-20);
        assertEquals(0, irregular.period);
        double ratio = Math.max(irregular.trough, tightChaos) / Math.min(irregular.trough,
                tightChaos);
        assertTrue("a chaotic one moves by orders of magnitude, not " + ratio, ratio > 100.0);
    }

    /** And where it is not the model's, it is under one infective either way. */
    @Test
    public void theIrregularRegimeIsBelowOneInfective() {
        double people = seir(0.20).trough * MeaslesDemo.MASSACHUSETTS_PEOPLE;
        assertTrue("the trough is " + people + " people", people < 1.0);
        double tightened = MeaslesDemo.troughAtTolerance(0.20, 1.0e-13, 1.0e-20)
                * MeaslesDemo.MASSACHUSETTS_PEOPLE;
        assertTrue("under one is the one thing both tolerances agree on: " + tightened,
                tightened < MeaslesDemo.MASSACHUSETTS_PEOPLE * 1.0e-3);
    }

    // -------------------------------------------------------------- the data

    /**
     * The cycle is strongest in the largest population and collapses in the
     * smallest, which is the opposite of what aggregation would do and is what
     * the model's own trough predicts. Monotone at all three lags.
     */
    @Test
    public void theCycleGrowsWithThePopulation() {
        double[][] series = Datasets.byDecreasingPopulation();
        String[] names = Datasets.byDecreasingPopulationNames();
        assertEquals(4, series.length);
        assertEquals(4, names.length);

        for (int lagYears = 1; lagYears <= 3; ++lagYears) {
            double previous = Double.MAX_VALUE;
            for (int i = 0; i < series.length; ++i) {
                double[] logged = MeaslesDemo
                        .logIncidence(MeaslesDemo.fillSingleGaps(series[i]));
                double[] a = ACF.acf(logged, 3 * Datasets.WEEKS_PER_YEAR + 1);
                double value = a[lagYears * Datasets.WEEKS_PER_YEAR];
                assertTrue(names[i] + " at " + lagYears + " years is " + value + ", above "
                        + previous, value < previous);
                previous = value;
            }
        }
    }

    /** Filling a single week cannot move a value that was never missing. */
    @Test
    public void fillingAGapTouchesNothingElse() {
        double[] raw = Datasets.massachusetts();
        double[] filled = MeaslesDemo.fillSingleGaps(raw);
        assertEquals(raw.length, filled.length);
        int changed = 0;
        for (int i = 0; i < raw.length; ++i) {
            if (Double.isNaN(raw[i])) {
                assertFalse("week " + i + " is still missing", Double.isNaN(filled[i]));
                ++changed;
            } else {
                assertEquals("week " + i, raw[i], filled[i], 0.0);
            }
        }
        assertEquals(Datasets.weeksWithNoReport(raw), changed);
    }

    // ------------------------------------------------------- the inverse step

    /**
     * The ringing period and the reproduction number are one function and its
     * inverse, so composing them has to give back what went in.
     */
    @Test
    public void theRecurrenceInvertsToAReproductionNumber() {
        // asserted in R0 and not in the period, because that is the quantity the
        // root finder is given a tolerance in. Going the other way would be
        // asserting that tolerance times a slope, plus the noise floor, which is
        // three quantities where one will do.
        for (double r0 : new double[] { 8.0, 12.0, 17.0, 20.0 }) {
            double period = MeaslesDemo.ringingPeriod(r0);
            assertEquals("R0 = " + r0 + " through the recurrence and back", r0,
                    MeaslesDemo.reproductionNumberFor(period), 0.05);
        }
        for (double period : new double[] { 2.0, 2.5, 3.0 }) {
            double r0 = MeaslesDemo.reproductionNumberFor(period);
            assertTrue("R0 = " + r0 + " is outside anything measles could be",
                    r0 > 5.0 && r0 < 25.0);
        }
        assertTrue("a shorter recurrence means a larger R0",
                MeaslesDemo.reproductionNumberFor(2.0) > MeaslesDemo.reproductionNumberFor(3.0));
    }

    /**
     * The instrument's own limit, measured. The ringing period is an event time
     * and not an integral, so it carries a noise floor that does not shrink with
     * the step: a millionth in {@code R0} moves it as much as a thousandth does.
     * That is why the demo asks Brent-Dekker for a hundredth and prints one
     * decimal place.
     */
    @Test
    public void theRingingPeriodHasANoiseFloorThatDoesNotShrink() {
        double small = spread(18.62, 1.0e-6);
        double large = spread(18.62, 1.0e-3);
        assertTrue("a 1e-6 step already moves the period by " + small,
                small > 0.2 * MeaslesDemo.PERIOD_NOISE);
        assertTrue("and a 1e-3 step moves it no further: " + small + " against " + large,
                small > 0.2 * large);
    }

    /** How far the ringing period moves when R0 is stepped by {@code h}. */
    private static double spread(double centre, double h) {
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;
        for (int i = -3; i <= 3; ++i) {
            double p = MeaslesDemo.ringingPeriod(centre + i * h);
            lo = Math.min(lo, p);
            hi = Math.max(hi, p);
        }
        return hi - lo;
    }

    /**
     * The two year reading puts the herd immunity threshold in the nineties,
     * which is the number the last section is about.
     */
    @Test
    public void theThresholdIsWhereTheLiteraturePutsIt() {
        double r0 = MeaslesDemo.reproductionNumberFor(2.0);
        double threshold = 1.0 - 1.0 / r0;
        assertTrue("R0 from a two year recurrence is " + r0, r0 > 17.0 && r0 < 20.0);
        assertTrue("threshold " + threshold, threshold > 0.94 && threshold < 0.95);
    }

    // ------------------------------------------------------------- the page

    /**
     * The page has to be the same page twice. It is the one property the whole
     * of section 5 exists to protect, and the empty cells are how.
     */
    @Test
    public void thePageIsTheSamePageTwice() {
        assertEquals(run(), run());
    }

    /**
     * An iteration count, an evaluation count and a wall time are properties of
     * the route rather than of the answer, and the tenth demo printed one and
     * had to take it back. None of the three is here.
     */
    @Test
    public void noRouteDependentQuantityReachesThePage() {
        String page = run().toLowerCase(java.util.Locale.ROOT);
        assertFalse("iterations", page.contains("iteration"));
        assertFalse("evaluations", page.contains("evaluation"));
        assertFalse("milliseconds", page.contains(" ms"));
        assertFalse("seconds", page.contains(" seconds"));
    }

    /** Everything the prose asserts about the tables is on the page. */
    @Test
    public void thePageSaysWhatTheMeasurementsSay() {
        String page = run();
        assertTrue(page.contains("1. The record"));
        assertTrue(page.contains("3. What the model settles into, and how that is measured"));
        assertTrue(page.contains("7. The number that decides the campaign"));
        assertTrue(page.contains(Integer.toString(Datasets.VACCINE_YEAR)));
        assertTrue("the calibrated row is marked", page.contains("<- the record"));
        assertTrue("and the caveat that makes the check worth running",
                page.contains("for running the check, not for knowing when to"));
        assertTrue("a cell whose two tolerances disagree carries no number",
                page.contains("--          no"));
    }

    /**
     * The calibration is the demo's spine: the record's heavy-to-light ratio is
     * a monotone function of the one free parameter, so it inverts. Asserted as
     * the identity it is.
     */
    @Test
    public void theRecordSetsTheOneFreeParameter() {
        for (double amplitude : new double[] { 0.07, 0.10, 0.13, 0.16 }) {
            double ratio = MeaslesDemo.seirRatio(amplitude);
            assertEquals("swing " + amplitude + " through the ratio and back", amplitude,
                    MeaslesDemo.forcingFor(ratio), 1.0e-3);
        }
        double previous = 0.0;
        for (double amplitude = 0.06; amplitude <= 0.181; amplitude += 0.02) {
            double ratio = MeaslesDemo.seirRatio(amplitude);
            assertTrue("the ratio has to grow with the swing to be invertible at all: " + ratio
                    + " after " + previous, ratio > previous);
            previous = ratio;
        }
    }

    /**
     * And it does not rest on the latent class, which is what lets section 4
     * blame the latent class for something else.
     */
    @Test
    public void theCalibrationDoesNotNeedTheLatentClass() {
        for (double amplitude : new double[] { 0.09, 0.12, 0.15, 0.18 }) {
            double with = MeaslesDemo.seirRatio(amplitude);
            double without = MeaslesDemo.sirRatio(amplitude);
            assertEquals("at swing " + amplitude, with, without, 0.15 * with);
        }
    }

    /**
     * The claim section 3 makes, asserted as what it is. Two starting states a
     * factor of three apart in the infectious count have to end up doing the
     * same thing, or the word attractor is not earned and none of the numbers
     * below it are properties of the equation.
     */
    @Test
    public void theModelForgetsWhereItStarted() {
        double amplitude = MeaslesDemo.forcingFor(MeaslesDemo.loudRatio());
        math.fun.DVectorField f = MeaslesDemo.seirField(amplitude);
        double[] a = MeaslesDemo.seirStart();
        double[] b = { a[0] * 1.10, a[1] * 0.5, a[2] * 3.0 };

        double atStart = gapBetween(f, a, b, 0.0);
        assertTrue("the two starts have to differ to begin with: " + atStart, atStart > 0.1);

        double previous = atStart;
        for (double years : new double[] { 25.0, 50.0, 100.0 }) {
            double gap = gapBetween(f, a, b, years);
            assertTrue("the gap has to keep shrinking, " + gap + " after " + previous,
                    gap < previous);
            previous = gap;
        }
        double settled = gapBetween(f, a, b, MeaslesDemo.SETTLE_YEARS);
        assertTrue("after " + MeaslesDemo.SETTLE_YEARS + " years the two agree to " + settled,
                settled < 1.0e-9);
    }

    /** How far apart the next epidemic of two runs is, after settling for a while. */
    private static double gapBetween(math.fun.DVectorField f, double[] a, double[] b,
            double years) {
        double pa = MeaslesDemo.nextPeak(f, MeaslesDemo.advance(f, a, years), years);
        double pb = MeaslesDemo.nextPeak(f, MeaslesDemo.advance(f, b, years), years);
        return Math.abs(pa - pb) / Math.max(pa, pb);
    }

    /** The record moves, which is the reason section 4 exists. */
    @Test
    public void theRecordIsNotInOneState() {
        double quiet = MeaslesDemo.quietRatio();
        double loud = MeaslesDemo.loudRatio();
        assertTrue("the loud stretch is " + loud + " against " + quiet, loud > 4.0 * quiet);
        assertTrue("and it calibrates higher",
                MeaslesDemo.forcingFor(loud) > 2.0 * MeaslesDemo.forcingFor(quiet));
        assertTrue("but still inside the window where the model answers",
                MeaslesDemo.agrees(MeaslesDemo.forcingFor(loud), true));
    }

    // ------------------------------------------------------------- machinery

    private static MeaslesDemo.Attractor sir(double amplitude) {
        return MeaslesDemo.attractorOf(MeaslesDemo.sirField(amplitude), 2, MeaslesDemo.sirStart(),
                amplitude);
    }

    private static MeaslesDemo.Attractor seir(double amplitude) {
        return MeaslesDemo.attractorOf(MeaslesDemo.seirField(amplitude), 3, MeaslesDemo.seirStart(),
                amplitude);
    }

    private static int periodOfSir(double amplitude) {
        return sir(amplitude).period;
    }

    private static int periodOfSeir(double amplitude) {
        return seir(amplitude).period;
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            MeaslesDemo.main(new String[0]);
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
