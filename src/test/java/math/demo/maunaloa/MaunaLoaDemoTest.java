package math.demo.maunaloa;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

/**
 * Nails down what {@link MaunaLoaDemo} prints, so that the demo cannot rot
 * silently as the classes underneath it change.
 * <p>
 * The bounds are measurements, not guesses: each was read off a run over the
 * embedded record and then loosened enough that ordinary numerical drift
 * cannot trip it, while a real change in behaviour still can.
 */
public final class MaunaLoaDemoTest {

    // ------------------------------------------------------------ the data

    @Test
    public void testTheEmbeddedRecordIsInternallyConsistent() {
        int n = Datasets.size();
        assertEquals(821, n);
        assertEquals(n, Datasets.co2().length);
        assertEquals(n, Datasets.uncertainty().length);
        assertEquals(n, Datasets.deseasonalized().length);
        assertEquals(n, Datasets.decimalDate().length);
        assertEquals("1958-03", Datasets.label(0));
        assertEquals("2026-07", Datasets.label(n - 1));

        double[] t = Datasets.decimalDate();
        for (int i = 1; i < n; ++i) {
            assertTrue("the abscissa must be strictly increasing at " + i, t[i] > t[i - 1]);
        }
        // the checksum recorded in the javadoc of Datasets: it says whether a
        // refresh of the NOAA file changed history or only appended to it
        double total = 0.0;
        double[] y = Datasets.co2();
        for (int i = 0; i < n; ++i) {
            total += y[i];
        }
        assertEquals(296610.72, total, 0.005);
    }

    @Test
    public void testTheUncertaintyIsUsableWhereItIsReported() {
        int usable = 0;
        double[] unc = Datasets.uncertainty();
        for (int i = 0; i < Datasets.size(); ++i) {
            if (Datasets.hasUncertainty(i)) {
                assertTrue("a usable uncertainty must be positive at " + i, unc[i] > 0.0);
                ++usable;
            }
        }
        assertEquals(625, usable);
        // the two months NOAA interpolated carry an uncertainty field of
        // exactly zero, which would become an infinite regression weight
        int[] interpolated = Datasets.noaaInterpolated();
        assertArrayEquals(new int[] { 213, 313 }, interpolated);
        assertEquals("1975-12", Datasets.label(interpolated[0]));
        assertEquals("1984-04", Datasets.label(interpolated[1]));
        for (int i = 0; i < interpolated.length; ++i) {
            assertEquals(0.0, unc[interpolated[i]], 0.0);
            assertTrue("an interpolated month must not count as weighted",
                    !Datasets.hasUncertainty(interpolated[i]));
        }
    }

    // -------------------------------------------------- 2. interpolation

    @Test
    public void testEveryInterpolatorRebuildsASingleMissingMonthWell() {
        // measured: 0.32 to 0.39 ppm RMS, against a series whose own seasonal
        // swing is 6 ppm
        for (int s = 0; s < MaunaLoaDemo.SCHEMES.length; ++s) {
            MaunaLoaDemo.HoldOut h = MaunaLoaDemo.holdOut(s, 1);
            assertEquals(819, h.cases);
            assertTrue(h.scheme + " RMS was " + h.rms, h.rms < 0.5);
            assertTrue(h.scheme + " worst was " + h.worst, h.worst < 4.0);
        }
    }

    @Test
    public void testTheErrorGrowsWithTheGap() {
        int[] gaps = { 1, 3, 6, 12 };
        for (int s = 0; s < MaunaLoaDemo.SCHEMES.length; ++s) {
            double previous = 0.0;
            for (int g = 0; g < gaps.length; ++g) {
                MaunaLoaDemo.HoldOut h = MaunaLoaDemo.holdOut(s, gaps[g]);
                assertTrue(h.scheme + " at gap " + gaps[g] + " did not get worse", h.rms > previous);
                previous = h.rms;
            }
        }
    }

    @Test
    public void testTheSchemesTradePlacesAsTheGapGrows() {
        // the finding worth keeping from step 2: over one month the natural
        // spline is the most accurate and Kruger the least, and over a year it
        // is the other way round. A shape preserving rule gives up a little on
        // short gaps and degrades far more gracefully on long ones
        double[] shortGap = new double[MaunaLoaDemo.SCHEMES.length];
        double[] longGap = new double[MaunaLoaDemo.SCHEMES.length];
        for (int s = 0; s < MaunaLoaDemo.SCHEMES.length; ++s) {
            shortGap[s] = MaunaLoaDemo.holdOut(s, 1).rms;
            longGap[s] = MaunaLoaDemo.holdOut(s, 12).rms;
        }
        int spline = 0;
        int kruger = MaunaLoaDemo.SCHEMES.length - 1;
        assertEquals("natural spline", MaunaLoaDemo.SCHEMES[spline]);
        assertEquals("Kruger", MaunaLoaDemo.SCHEMES[kruger]);
        assertTrue("over one month the spline should win", shortGap[spline] < shortGap[kruger]);
        assertTrue("over twelve months Kruger should win", longGap[kruger] < longGap[spline]);
    }

    @Test
    public void testOurFillingOfNoaasOwnTwoGapsIsCloseToTheirs() {
        for (int s = 0; s < MaunaLoaDemo.SCHEMES.length; ++s) {
            double[] errors = MaunaLoaDemo.versusNoaa(s);
            assertEquals(2, errors.length);
            for (int m = 0; m < errors.length; ++m) {
                assertTrue(MaunaLoaDemo.SCHEMES[s] + " differed by " + errors[m], Math.abs(errors[m]) < 1.5);
            }
        }
    }

    // -------------------------------------------------- 3. and 4. the fit

    @Test
    public void testTheTrendAndSeasonFitTheRecord() {
        MaunaLoaDemo.Fit fit = MaunaLoaDemo.fitOrdinary();
        assertEquals(821, fit.observations);
        assertTrue("R^2 was " + fit.rSquared, fit.rSquared > 0.999);
        // the seasonal swing at Mauna Loa is about 6 ppm peak to peak
        assertEquals(6.0, fit.seasonalPeakToPeak, 1.0);
        assertEquals(2.8, fit.annualAmplitude, 0.5);
        assertTrue("the semi annual term should be the smaller one",
                fit.semiAnnualAmplitude < fit.annualAmplitude);
        // the rise is accelerating, which is the whole reason for the squared
        // term, and the present rate is about 2.4 to 2.6 ppm per year
        double[] tc = MaunaLoaDemo.centredTime();
        double now = fit.growth(tc[tc.length - 1]);
        assertTrue("growth in 1958 was " + fit.growth(0.0), fit.growth(0.0) < now);
        assertEquals(2.5, now, 0.4);
    }

    @Test
    public void testOurSeasonalAdjustmentAgreesWithNoaasOwn() {
        // an oracle from outside this repository: NOAA ships its own
        // deseasonalized column, computed by a different method
        double rms = MaunaLoaDemo.versusNoaaDeseasonalized(MaunaLoaDemo.fitOrdinary());
        assertTrue("RMS gap to NOAA was " + rms + " ppm", rms < 0.5);
    }

    @Test
    public void testWeightingMovesTheFitButNotFarOnThisData() {
        MaunaLoaDemo.Fit ordinary = MaunaLoaDemo.fitOnUncertainMonths(false);
        MaunaLoaDemo.Fit weighted = MaunaLoaDemo.fitOnUncertainMonths(true);
        assertEquals(625, ordinary.observations);
        assertEquals(625, weighted.observations);
        for (int j = 0; j < ordinary.beta.length; ++j) {
            double shift = Math.abs(weighted.beta[j] - ordinary.beta[j]) / ordinary.standardErrors[j];
            assertTrue("coefficient " + j + " moved " + shift + " standard errors", shift < 3.0);
        }
        // the weights are the reported uncertainties, which span a factor of
        // about ten, so the two fits must not come out identical either
        double biggest = 0.0;
        for (int j = 0; j < ordinary.beta.length; ++j) {
            biggest = Math.max(biggest, Math.abs(weighted.beta[j] - ordinary.beta[j]));
        }
        assertTrue("weighting changed nothing at all", biggest > 0.0);
    }

    // -------------------------------------------------- 5. and 6.

    @Test
    public void testTheResidualsAreStronglyAutocorrelated() {
        // this is what makes the standard errors of step 3 too small, and it
        // is the reason step 6 does not resample them
        MaunaLoaDemo.Correlation c = MaunaLoaDemo.correlate(MaunaLoaDemo.fitOrdinary());
        assertEquals(1.0, c.acf[0], 1.0e-12);
        assertTrue("lag 1 was " + c.acf[1], c.acf[1] > 0.8);
        for (int k = 1; k <= 24; ++k) {
            assertTrue("lag " + k + " should be significant, was " + c.acf[k], Math.abs(c.acf[k]) > c.band);
        }
    }

    @Test
    public void testTheAnnualIncrementsAreWeaklyEnoughRelatedToResample() {
        // the justification for the iid bootstrap in step 6. Over the whole
        // record the same number is large, but that is the long run rise of
        // the increments rather than dependence
        MaunaLoaDemo.Correlation c = MaunaLoaDemo.correlate(MaunaLoaDemo.fitOrdinary());
        assertTrue("recent increments lag 1 was " + c.incrementLag1,
                Math.abs(c.incrementLag1) < c.incrementBand);
        assertTrue("the whole record should look worse, was " + c.allIncrementsLag1,
                c.allIncrementsLag1 > c.incrementLag1);
    }

    @Test
    public void testTheGrowthIntervalIsReproducibleAndAgreesWithTheClassicalOne() {
        MaunaLoaDemo.Growth first = MaunaLoaDemo.growth(MaunaLoaDemo.GROWTH_YEARS);
        MaunaLoaDemo.Growth second = MaunaLoaDemo.growth(MaunaLoaDemo.GROWTH_YEARS);
        // a printed confidence interval that cannot be reproduced is not a
        // result; tolerance 0.0 on purpose
        assertEquals(first.percentile[0], second.percentile[0], 0.0);
        assertEquals(first.percentile[1], second.percentile[1], 0.0);
        assertEquals(first.bca[0], second.bca[0], 0.0);
        assertEquals(first.bca[1], second.bca[1], 0.0);

        assertEquals(30, first.years);
        assertTrue("percentile must bracket the observed mean",
                first.percentile[0] < first.observed && first.observed < first.percentile[1]);
        assertTrue("BCa must bracket the observed mean",
                first.bca[0] < first.observed && first.observed < first.bca[1]);
        // measured at about 2.2 ppm per year over the last thirty years
        assertEquals(2.2, first.observed, 0.4);
        // once the dependence has been checked, the bootstrap and the
        // classical formula have to say the same thing
        double bootstrapHalfWidth = 0.5 * (first.percentile[1] - first.percentile[0]);
        assertEquals(1.0, bootstrapHalfWidth / first.classicalHalfWidth, 0.2);
    }

    // -------------------------------------------------- 7. and 8.

    @Test
    public void testTheTwoIntegrationRoutesAgree() {
        // the spline integrates itself in closed form; Gauss-Kronrod K15 is
        // exact for a cubic on one interval. Neither knows about the other,
        // so agreement to round-off is a real check on both
        MaunaLoaDemo.Excess e = MaunaLoaDemo.excess(MaunaLoaDemo.fitOrdinary());
        assertTrue("the excess should be positive, was " + e.closedForm, e.closedForm > 0.0);
        assertTrue("relative difference was " + e.relativeDifference, e.relativeDifference < 1.0e-12);
    }

    @Test
    public void testTheCrossingIsARootOfTheTrend() {
        MaunaLoaDemo.Fit fit = MaunaLoaDemo.fitOrdinary();
        MaunaLoaDemo.Crossing c = MaunaLoaDemo.crossing(fit, 450.0);
        assertTrue("the root should be within the next few decades, was " + c.year,
                c.year > 2027.0 && c.year < 2060.0);
        assertEquals("Brent should land on the level", 0.0, c.residual, 1.0e-8);
        double t0 = Datasets.decimalDate()[0];
        assertEquals(450.0, fit.trend(c.year - t0), 1.0e-8);
    }

    @Test
    public void testTheTransformFindsTheSeasonNobodyToldItAbout() {
        // the length is trimmed to whole years so that a twelve month cycle
        // lands on a bin exactly; 816 is not a power of two, so this is the
        // Bluestein path
        MaunaLoaDemo.Spectrum s = MaunaLoaDemo.spectrum(MaunaLoaDemo.fitOrdinary());
        assertEquals(816, s.months);
        assertEquals(0, s.months % 12);
        assertEquals("the strongest line must be the annual cycle", 12.0,
                s.periods[MaunaLoaDemo.strongest(s, 0)], 1.0e-9);
        assertEquals("the second must be its first harmonic", 6.0,
                s.periods[MaunaLoaDemo.strongest(s, 1)], 1.0e-9);
    }

    // -------------------------------------------------- the demo as a whole

    @Test
    public void testTheDemoRunsAndSaysTheSameThingEveryTime() throws UnsupportedEncodingException {
        String first = runMain();
        String second = runMain();
        assertEquals("the demo must be reproducible from end to end", first, second);
        assertTrue("the demo printed nothing", first.length() > 2000);
        String[] required = { "1. the record", "2. take months out", "3. a trend and a seasonal cycle",
                "math.linalg.Wls", "5. the residuals are not independent", "math.probe.Bootstrap",
                "7. cumulative excess", "8. when does the trend reach 450 ppm" };
        for (int i = 0; i < required.length; ++i) {
            assertTrue("the demo no longer prints \"" + required[i] + "\"", first.contains(required[i]));
        }
    }

    private static String runMain() throws UnsupportedEncodingException {
        PrintStream saved = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, "US-ASCII"));
            MaunaLoaDemo.main(new String[0]);
        } finally {
            System.setOut(saved);
        }
        return captured.toString("US-ASCII");
    }
}
