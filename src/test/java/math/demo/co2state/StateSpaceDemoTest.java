package math.demo.co2state;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.Locale;

import org.junit.Test;

import math.demo.maunaloa.Datasets;
import math.linalg.DMatrix;
import math.ts.KalmanFilter;

/**
 * Guards what {@link StateSpaceDemo} prints.
 * <p>
 * The fit is the expensive part -- Nelder-Mead over four variances, each
 * evaluation a pass of the filter over 821 months -- so it is done once and
 * shared. Everything downstream of it is milliseconds.
 * <p>
 * The assertions are of two kinds. Some pin a number the demo prints, so that a
 * change in the library that moves it has to be noticed; the tolerances there
 * are loose enough to survive a different rounding and tight enough to fail a
 * different answer. The rest are the structural claims the demo makes in
 * words -- that the innovations are white where the deterministic residuals
 * were not, that a repeated schedule and a missing month are different things,
 * that the interval brackets the point estimate the other demo gives -- and
 * those are what would actually be wrong if the model were wrong.
 */
public final class StateSpaceDemoTest {

    /** The one expensive thing, fitted once for the whole class. */
    private static StateSpaceDemo.Fit fitted;

    private static synchronized StateSpaceDemo.Fit fit() {
        if (fitted == null) {
            fitted = StateSpaceDemo.fit(StateSpaceDemo.HARMONICS);
        }
        return fitted;
    }

    // -------------------------------------------------------- the model itself

    @Test
    public void testTheTransitionIsAnIdentityAtAZeroGap() {
        // F(0) has to be the identity: no time passes, nothing moves. That is
        // what makes the dt scaling of step 7 a generalization of step 1
        // rather than a second model
        DMatrix f = StateSpaceDemo.transition(StateSpaceDemo.HARMONICS, 0.0);
        int d = StateSpaceDemo.dimension(StateSpaceDemo.HARMONICS);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                assertEquals("at (" + i + ", " + j + ")", i == j ? 1.0 : 0.0, f.get(i, j), 1.0e-15);
            }
        }
    }

    @Test
    public void testASeasonalPairIsARotationAndKeepsItsLength() {
        DMatrix f = StateSpaceDemo.transition(StateSpaceDemo.HARMONICS, 1.0);
        for (int j = 0; j < StateSpaceDemo.HARMONICS; ++j) {
            int at = 2 + 2 * j;
            double c = f.get(at, at);
            double s = f.get(at, at + 1);
            assertEquals("harmonic " + (j + 1), 1.0, c * c + s * s, 1.0e-14);
            assertEquals("the block must be a rotation", c, f.get(at + 1, at + 1), 1.0e-15);
            assertEquals("the block must be a rotation", -s, f.get(at + 1, at), 1.0e-15);
        }
    }

    @Test
    public void testTwelveMonthsOfRotationComeBackToWhereTheyStarted() {
        // a seasonal component that did not return after a year would not be
        // seasonal, and it is the one property of the model that can be
        // checked without any data at all
        DMatrix year = StateSpaceDemo.transition(StateSpaceDemo.HARMONICS, 12.0);
        for (int j = 0; j < StateSpaceDemo.HARMONICS; ++j) {
            int at = 2 + 2 * j;
            assertEquals("harmonic " + (j + 1), 1.0, year.get(at, at), 1.0e-14);
            assertEquals("harmonic " + (j + 1), 0.0, year.get(at, at + 1), 1.0e-14);
        }
    }

    @Test
    public void testTheProcessNoiseScalesWithTheGap() {
        double[] v = fit().logVariance;
        DMatrix one = StateSpaceDemo.processNoise(StateSpaceDemo.HARMONICS, v, 1.0);
        DMatrix five = StateSpaceDemo.processNoise(StateSpaceDemo.HARMONICS, v, 5.0);
        int d = StateSpaceDemo.dimension(StateSpaceDemo.HARMONICS);
        for (int i = 0; i < d; ++i) {
            assertEquals("at " + i, 5.0 * one.get(i, i), five.get(i, i), 1.0e-18);
        }
    }

    // ------------------------------------------------------------ 1. the fit

    @Test
    public void testTheFitConverges() {
        StateSpaceDemo.Fit fit = fit();
        assertTrue("Nelder-Mead should converge, not run out of budget", fit.converged);
        assertEquals(StateSpaceDemo.HARMONICS, fit.harmonics);
        assertEquals(4, fit.logVariance.length);
        // measured -255.687 with three harmonics
        assertEquals(-255.687, fit.logLikelihood, 0.05);
        for (int i = 0; i < 4; ++i) {
            assertTrue("every variance must be positive and finite at " + i,
                    fit.sigma(i) > 0.0 && fit.sigma(i) < 10.0);
        }
    }

    @Test
    public void testTheLevelIsAllowedToMoveAboutAsMuchAsTheMeasurementIsWrong() {
        // the four variances only matter through their ratio, and this is the
        // ratio: measured 0.685, so the level disturbance is a little smaller
        // than the observation error rather than larger
        StateSpaceDemo.Fit fit = fit();
        assertEquals("sigma level", 0.1668, fit.sigma(0), 0.01);
        assertEquals("sigma observation", 0.2015, fit.sigma(3), 0.01);
        assertTrue("the level disturbance is the smaller of the two, "
                + fit.sigma(0) + " against " + fit.sigma(3), fit.sigma(0) < fit.sigma(3));
        assertEquals("signal to noise", 0.685, fit.signalToNoise(), 0.1);
        assertTrue("the slope should wander far more slowly than the level",
                fit.sigma(1) < 0.05 * fit.sigma(0));
    }

    @Test
    public void testTheFilterSplitsEachSurpriseAboutInHalf() {
        // what the ratio buys, and the number the demo prints: a gain of zero
        // is a rigid curve and a gain of one is no smoothing at all
        double gain = StateSpaceDemo.steadyStateLevelGain(fit());
        // measured 0.534
        assertEquals(0.534, gain, 0.05);
        assertTrue("a gain has to lie between the two corners", gain > 0.0 && gain < 1.0);
    }

    @Test
    public void testTheThirdHarmonicIsWorthKeeping() {
        // step 1 is a model choice made by the likelihood rather than by hand,
        // so what it chose has to stay chosen
        StateSpaceDemo.Fit two = StateSpaceDemo.fit(StateSpaceDemo.DETERMINISTIC_HARMONICS);
        assertTrue("a third harmonic must improve the likelihood, was " + two.logLikelihood
                + " against " + fit().logLikelihood, fit().logLikelihood > two.logLikelihood);
        // measured 7.64 nats
        assertEquals(7.64, fit().logLikelihood - two.logLikelihood, 0.5);
    }

    // ------------------------------------------------------ 2. the decomposition

    @Test
    public void testTheLevelAgreesWithNoaasOwnDeseasonalizedSeries() {
        // the sharpest check in the demo: NOAA removes the seasonal cycle by a
        // method that has nothing in common with a Kalman smoother, and
        // publishes the result beside the data
        StateSpaceDemo.Decomposition parts = StateSpaceDemo.decompose(fit());
        assertEquals(Datasets.size(), parts.level.length);
        // measured 0.1658 ppm on a series that runs from 315 to 429 ppm
        assertTrue("rms against NOAA was " + parts.rmsAgainstNoaa, parts.rmsAgainstNoaa < 0.25);
        assertTrue("worst against NOAA was " + parts.worstAgainstNoaa, parts.worstAgainstNoaa < 1.0);
    }

    @Test
    public void testTheDecompositionAddsUp() {
        // level plus season is what the model claims the observation is, so
        // the two have to reconstruct the data to within the measurement error
        StateSpaceDemo.Decomposition parts = StateSpaceDemo.decompose(fit());
        double[] co2 = Datasets.co2();
        double sum = 0.0;
        for (int t = 0; t < co2.length; ++t) {
            double gap = co2[t] - (parts.level[t] + parts.season[t]);
            sum += gap * gap;
        }
        double rms = Math.sqrt(sum / co2.length);
        assertTrue("the parts should reconstruct the series to about the observation error, "
                + "rms was " + rms, rms < 3.0 * fit().sigma(3));
    }

    @Test
    public void testTheRiseHasAccelerated() {
        StateSpaceDemo.Decomposition parts = StateSpaceDemo.decompose(fit());
        int n = Datasets.size();
        double recent = 12.0 * parts.slope[n - 1];
        double early = 12.0 * parts.slope[24];
        assertTrue("the recent slope should be the larger, " + recent + " against " + early,
                recent > 2.0 * early);
        // measured 2.407 ppm per year now against 0.896 sixty-six years ago
        assertEquals(2.407, recent, 0.3);
        assertEquals(0.896, early, 0.3);
        // and the seasonal cycle is a few ppm peak to trough, measured 6.611
        assertEquals(6.611, parts.seasonalAmplitude, 0.5);
    }

    // -------------------------------------------------------- 3. the whiteness

    @Test
    public void testTheInnovationsAreWhiteWhereTheResidualsWereNot() {
        StateSpaceDemo.Whiteness white = StateSpaceDemo.whiteness(fit());
        // MaunaLoaDemo step 5 reports this same 0.906 about the same design
        assertEquals("the deterministic fit's acf(1)", 0.906, white.deterministicAcf1, 0.01);
        assertEquals("every lag of it leaves the band", 24, white.deterministicOutside);
        // and the structural model's innovations, measured 0.0304
        assertTrue("the innovations should be far less correlated, acf(1) was "
                + white.structuralAcf1, Math.abs(white.structuralAcf1) < 0.1);
        assertTrue("acf(1) should sit inside the band", Math.abs(white.structuralAcf1) < white.bound);
        assertTrue("far fewer lags should leave the band, " + white.structuralOutside + " did",
                white.structuralOutside < 10);
    }

    @Test
    public void testTheInnovationsAreStandardizedAsClaimed() {
        // e / sqrt(S) is standard normal if the model is right, which is what
        // makes it an anomaly score rather than a residual
        StateSpaceDemo.Whiteness white = StateSpaceDemo.whiteness(fit());
        assertEquals("mean", 0.0, white.innovationMean, 0.15);
        assertEquals("standard deviation", 1.0, white.innovationSd, 0.1);
    }

    // ---------------------------------------------- 4. the interpolated months

    @Test
    public void testTheSmootherPutsBackWhatNoaaInterpolated() {
        StateSpaceDemo.Filled filled = StateSpaceDemo.fillTheInterpolatedMonths(fit());
        assertEquals(2, filled.index.length);
        assertEquals("1975-12", Datasets.label(filled.index[0]));
        assertEquals("1984-04", Datasets.label(filled.index[1]));
        for (int i = 0; i < filled.index.length; ++i) {
            double difference = Math.abs(filled.smoothed[i] - filled.noaa[i]);
            assertTrue("month " + Datasets.label(filled.index[i]) + " was off by " + difference,
                    difference < 0.6);
            assertTrue("the smoother must report its own uncertainty", filled.sd[i] > 0.0);
            // and its answer has to be consistent with that uncertainty
            assertTrue("the difference should sit inside two standard deviations",
                    difference < 2.5 * filled.sd[i]);
        }
    }

    @Test
    public void testABlankedMonthIsAStepWithNoCorrection() {
        // what the phrase means arithmetically. Note what is *not* asserted:
        // that blanking a month lowers the likelihood. A log density is not a
        // log probability and is positive whenever the density exceeds one,
        // which at an innovation standard deviation of 0.2 ppm most of these
        // terms do -- so dropping a month can move the total either way
        double[] co2 = Datasets.co2();
        double[] blanked = co2.clone();
        blanked[400] = Double.NaN;
        KalmanFilter.Result holed = KalmanFilter.filter(fit().model(),
                StateSpaceDemo.column(blanked));

        assertEquals("nothing was observed there", 0, holed.observedComponents[400]);
        assertEquals("so there is no innovation", 0, holed.innovation[400].length);
        assertTrue("and no anomaly score", Double.isNaN(holed.squaredMahalanobis[400]));
        // the exact statement: with nothing to correct against, the filtered
        // state is the predicted one, bit for bit
        for (int i = 0; i < StateSpaceDemo.dimension(StateSpaceDemo.HARMONICS); ++i) {
            assertEquals("coordinate " + i, holed.predictedMean[400][i], holed.filteredMean[400][i],
                    0.0);
            assertEquals("variance " + i, holed.predictedCovariance[400].get(i, i),
                    holed.filteredCovariance[400].get(i, i), 0.0);
        }
        // and every other month still contributes
        assertEquals(1, holed.observedComponents[399]);
        assertEquals(1, holed.observedComponents[401]);
    }

    // ------------------------------------------------------- 5. the anomalies

    @Test
    public void testTheMoveToMaunakeaLeftNoOffset() {
        StateSpaceDemo.Anomalies found = StateSpaceDemo.anomalies(fit());
        assertEquals("the observatory moved for eight months", 8, found.moveMonths);
        assertTrue("a mean of " + found.meanDuringMove + " should be inside one standard error of "
                + found.moveStandardError(), Math.abs(found.meanDuringMove) < found.moveStandardError());
    }

    @Test
    public void testTheLargestAnomalyIsTheElNino() {
        StateSpaceDemo.Anomalies found = StateSpaceDemo.anomalies(fit());
        // measured 2016-04 at 4.05 standard deviations
        assertEquals("2016-04", Datasets.label(found.largestAt));
        assertTrue("it should stand out clearly, was " + found.largest, found.largest > 3.5);
        assertTrue("and it should be the largest, including during the move",
                found.largest > found.largestDuringMove);
    }

    @Test
    public void testTheMonthsOfTheMoveAreTheOnesTheDatasetNames() {
        int during = 0;
        for (int t = 0; t < Datasets.size(); ++t) {
            if (StateSpaceDemo.measuredAtMaunakea(t)) {
                ++during;
            }
        }
        assertEquals(8, during);
        assertTrue(StateSpaceDemo.measuredAtMaunakea(indexOf("2022-12")));
        assertTrue(StateSpaceDemo.measuredAtMaunakea(indexOf("2023-07")));
        assertTrue("the month before must not count",
                !StateSpaceDemo.measuredAtMaunakea(indexOf("2022-11")));
        assertTrue("nor the month after",
                !StateSpaceDemo.measuredAtMaunakea(indexOf("2023-08")));
    }

    private static int indexOf(String label) {
        for (int t = 0; t < Datasets.size(); ++t) {
            if (Datasets.label(t).equals(label)) {
                return t;
            }
        }
        throw new IllegalArgumentException("no such month : " + label);
    }

    // --------------------------------------------- 6. the published uncertainty

    @Test
    public void testThePublishedUncertaintyEarnsItsPlace() {
        StateSpaceDemo.Published published = StateSpaceDemo.publishedUncertainty(fit());
        assertEquals("the months that carry one", 625, published.months);
        // the published spread brackets the single constant the fit chose,
        // which is what makes the comparison fair rather than a rescaling
        assertTrue("published " + published.minSd + " to " + published.maxSd + " should bracket "
                + published.fittedSd,
                published.minSd < published.fittedSd && published.fittedSd < published.maxSd);
        // measured 23.25 nats, for no extra fitted parameter
        assertTrue("the published uncertainties should raise the likelihood, gain was "
                + published.gain(), published.gain() > 10.0);
        // but they barely move the level, which is the point the demo makes
        assertTrue("the level should hardly move, moved " + published.worstLevelShift,
                published.worstLevelShift < 1.0);
    }

    // ------------------------------------------------- 7. the irregular schedule

    @Test
    public void testAnIrregularScheduleRecoversTheLevel() {
        StateSpaceDemo.Irregular thin = StateSpaceDemo.thinToAnIrregularSchedule(fit());
        assertEquals(Datasets.size(), thin.total);
        assertTrue("it should keep roughly three in ten, kept " + thin.kept,
                thin.kept > 200 && thin.kept < 290);
        assertTrue("the schedule should be genuinely irregular, largest gap was " + thin.largestGap,
                thin.largestGap > 6);
        // measured 0.0932 ppm rms against the answer from the whole record
        assertTrue("rms against the full record was " + thin.rmsAgainstFull,
                thin.rmsAgainstFull < 0.2);
        // and, the claim the demo actually makes, barely worse against NOAA
        StateSpaceDemo.Decomposition parts = StateSpaceDemo.decompose(fit());
        assertTrue("throwing the months away should cost little against NOAA: "
                + thin.rmsAgainstNoaa + " against " + parts.rmsAgainstNoaa,
                thin.rmsAgainstNoaa - parts.rmsAgainstNoaa < 0.05);
    }

    // ------------------------------------------------------- 8. the crossing

    @Test
    public void testTheCrossingBracketsTheDeterministicAnswer() {
        StateSpaceDemo.Crossing crossing = StateSpaceDemo.crossing(fit(), 450.0);
        assertTrue("the three horizons must be ordered",
                crossing.earliest < crossing.central && crossing.central < crossing.latest);
        // MaunaLoaDemo finds 2034.42 by rooting its fitted trend; the interval
        // here has to contain it, or the two models disagree about the future
        assertEquals("2034-04", crossing.label(crossing.central));
        assertTrue("the band should contain the deterministic answer",
                crossing.label(crossing.earliest).compareTo("2034-06") < 0
                        && crossing.label(crossing.latest).compareTo("2034-06") > 0);
        // measured 5.0 years wide
        assertTrue("the window was " + crossing.widthInYears() + " years",
                crossing.widthInYears() > 2.0 && crossing.widthInYears() < 12.0);
    }

    @Test
    public void testForecastingWidensTheIntervalAndNeverNarrowsIt() {
        // what a forecast is: predict() with nothing to correct against, so
        // the covariance can only grow
        StateSpaceDemo.Fit fit = fit();
        KalmanFilter kf = new KalmanFilter(fit.model());
        double[] co2 = Datasets.co2();
        double[] row = new double[1];
        for (int t = 0; t < co2.length; ++t) {
            if (t > 0) {
                kf.predict();
            }
            row[0] = co2[t];
            kf.update(row);
        }
        double previous = kf.covariance().get(0, 0);
        for (int step = 0; step < 120; ++step) {
            kf.predict();
            double now = kf.covariance().get(0, 0);
            assertTrue("the level variance narrowed at step " + step, now >= previous);
            previous = now;
        }
    }

    // -------------------------------------------------------------- the output

    /**
     * The demo runs twice in this class and no more: once here under
     * {@link Locale#ROOT} and cached, and once under a comma locale. Two runs
     * is the floor, because the second is what the locale comparison needs, and
     * each is two maximum likelihood fits over 821 months.
     */
    private static String rootOutput;

    private static synchronized String output() throws UnsupportedEncodingException {
        if (rootOutput == null) {
            Locale saved = Locale.getDefault();
            try {
                Locale.setDefault(Locale.ROOT);
                rootOutput = runMain();
            } finally {
                Locale.setDefault(saved);
            }
        }
        return rootOutput;
    }

    @Test
    public void testTheDemoPrintsTheSamePageUnderAnyLocale() throws UnsupportedEncodingException {
        // both house requirements for a demo at once: the page has to come out
        // identical twice, and a locale whose decimal separator is a comma must
        // not change it. That is the invariant math.OutputLocaleTest asserts
        // for the library's own printers, applied to two dozen numbers here
        String root = output();
        Locale saved = Locale.getDefault();
        String comma;
        try {
            Locale.setDefault(Locale.GERMANY);
            comma = runMain();
        } finally {
            Locale.setDefault(saved);
        }
        assertEquals("the demo must print the same page twice, in any locale", root, comma);
        assertTrue("the demo printed nothing", root.length() > 2000);
        assertTrue("the page should be full of decimal points", root.contains("450 ppm"));
    }

    @Test
    public void testTheDemoStillPrintsItsEightSteps() throws UnsupportedEncodingException {
        String output = output();
        String[] required = { "1. the residuals the other demo could not get rid of",
                "2. what the four variances are", "3. level, slope and season",
                "4. the two months NOAA interpolated", "5. the eruption",
                "6. the uncertainty NOAA publishes", "7. seven months in ten thrown away",
                "8. when does it reach 450 ppm" };
        for (int i = 0; i < required.length; ++i) {
            assertTrue("the demo no longer prints \"" + required[i] + "\"",
                    output.contains(required[i]));
        }
    }

    /**
     * The page has to state the question it answers, and answer it first. The
     * demo exists because {@code MaunaLoaDemo} diagnoses a defect about itself,
     * and until 2026-08-30 that sentence lived only in the class comment while
     * the page opened on a procedure and buried the answer at position three.
     *
     * @throws UnsupportedEncodingException
     *             if UTF-8 is not available, which it is
     */
    @Test
    public void testThePageAsksItsQuestionBeforeAnsweringIt()
            throws UnsupportedEncodingException {
        String output = output();
        int question = output.indexOf("acf(1) = 0.906");
        int answer = output.indexOf("1. the residuals the other demo could not get rid of");
        int machinery = output.indexOf("2. what the four variances are");
        assertTrue("the page never states the defect it exists to repair", question >= 0);
        assertTrue("and it has to come before any section", question < answer);
        assertTrue("the answer comes before the machinery", answer < machinery);
    }

    /**
     * Section 2's sweep is the demo explaining its own vocabulary by turning the
     * one knob from end to end, so the shape of it is asserted rather than the
     * digits: the first autocorrelation has to fall monotonically as the level
     * is let loose, pass through white where the likelihood put it, and come
     * back with the other sign once the level is chasing the noise.
     */
    @Test
    public void testTheSweepGoesFromRigidThroughWhiteToOvershooting() {
        StateSpaceDemo.Knob[] sweep = StateSpaceDemo.knobSweep(fit());
        assertTrue("the sweep needs both ends and the middle", sweep.length >= 5);

        double previous = Double.MAX_VALUE;
        int fitted = -1;
        for (int i = 0; i < sweep.length; ++i) {
            assertTrue(sweep[i].label + " does not continue the sweep: " + sweep[i].acf1
                    + " after " + previous, sweep[i].acf1 < previous);
            previous = sweep[i].acf1;
            if (sweep[i].fitted) {
                fitted = i;
            }
        }
        assertTrue("exactly one row is the fitted one", fitted > 0 && fitted < sweep.length - 1);

        assertTrue("a frozen level leaves the innovations correlated: " + sweep[0].acf1,
                sweep[0].acf1 > 0.9);
        assertEquals("and every lag outside the band", 24, sweep[0].outside);
        assertTrue("the fitted row is white: " + sweep[fitted].acf1,
                Math.abs(sweep[fitted].acf1) < 0.05);
        assertTrue("and the far end overshoots: " + sweep[sweep.length - 1].acf1,
                sweep[sweep.length - 1].acf1 < -0.1);

        assertTrue("the gain has to rise with the freedom",
                sweep[0].gain < sweep[fitted].gain
                        && sweep[fitted].gain < sweep[sweep.length - 1].gain);
        assertTrue("a frozen level is a coefficient, so its gain is nearly nothing",
                sweep[0].gain < 0.01);
    }

    /**
     * And the point of the sweep: two criteria that share no arithmetic agree.
     * The likelihood never sees an autocorrelation, and it lands where the
     * autocorrelation is smallest.
     */
    @Test
    public void testTheLikelihoodAndTheWhitenessAgreeWithoutBeingTold() {
        StateSpaceDemo.Knob[] sweep = StateSpaceDemo.knobSweep(fit());
        int whitest = 0;
        for (int i = 1; i < sweep.length; ++i) {
            if (Math.abs(sweep[i].acf1) < Math.abs(sweep[whitest].acf1)) {
                whitest = i;
            }
        }
        assertTrue("the whitest row is not the one the likelihood chose",
                sweep[whitest].fitted);
    }

    /**
     * The four capability sections have one claim in common, and saying it is
     * what keeps them from reading as a feature tour.
     *
     * @throws UnsupportedEncodingException
     *             if UTF-8 is not available, which it is
     */
    @Test
    public void testTheCapabilitySectionsShareAStatedReason()
            throws UnsupportedEncodingException {
        String output = output();
        int reason = output.indexOf("a fitted curve can do");
        int missing = output.indexOf("4. the two months NOAA interpolated");
        assertTrue("the four sections carry no shared claim", reason >= 0);
        assertTrue("which has to be made before them", reason < missing);
    }

    @Test
    public void testThePageCarriesNothingThatDependsOnTheSearchPath()
            throws UnsupportedEncodingException {
        // the regression guard for a measured defect: Nelder-Mead takes 439
        // iterations on JDK 8 and 432 on JDK 25, because Math.exp and Math.cos
        // are allowed an ulp and take different intrinsics, which flips one
        // simplex comparison. Both land within 6e-12 of the same optimum, so
        // every figure on the page agrees -- but a count of iterations is a
        // property of the route, and printing it would break the rule that a
        // demo prints the same page on both runtimes
        assertTrue("the iteration count must stay off the page",
                !output().contains("iteration"));
        // it is still available to a caller that wants it
        assertTrue("but the fit still reports it", fit().iterations > 0);
    }

    private static String runMain() throws UnsupportedEncodingException {
        PrintStream saved = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, "UTF-8"));
            StateSpaceDemo.main(new String[0]);
        } finally {
            System.setOut(saved);
        }
        String text = captured.toString("UTF-8");
        assertNotNull(text);
        return text;
    }
}
