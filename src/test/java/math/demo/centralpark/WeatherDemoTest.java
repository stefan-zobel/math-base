package math.demo.centralpark;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.Locale;

import org.junit.Test;

import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.ts.KalmanFilter;
import math.ts.LinearGaussianModel;
import math.ts.RtsSmoother;

/**
 * Guards what {@link WeatherDemo} prints.
 * <p>
 * The four fits are the expensive part -- Nelder-Mead over up to five
 * parameters, each evaluation a pass of the filter over 365 days of a state of
 * eleven -- so the ladder is fitted once and shared, and so is the coarsened
 * refit of step 4. Everything else is milliseconds.
 * <p>
 * The assertions are of two kinds. Some pin a number the demo prints, so that a
 * change in the library that moves it has to be noticed; the tolerances there
 * are loose enough to survive a different rounding and tight enough to fail a
 * different answer. The rest are the claims the demo makes in words -- that a
 * filter with no process noise is least squares, that the off-diagonal of
 * {@code R} means one thing at M1 and another at M3, that a trace day is a dry
 * day, that blanking one reading leaves the days before it untouched -- and
 * those are what would actually be wrong if the model were wrong.
 */
public final class WeatherDemoTest {

    /** The expensive thing, fitted once for the whole class. */
    private static WeatherDemo.Fit[] rungs;

    /** The second expensive thing: step 4 refits on coarsened readings. */
    private static WeatherDemo.Grid resolution;

    /** The page as it comes out under the root locale. */
    private static String rootOutput;

    private static synchronized WeatherDemo.Fit[] ladder() {
        if (rungs == null) {
            rungs = WeatherDemo.ladder();
        }
        return rungs;
    }

    private static WeatherDemo.Fit working() {
        return ladder()[WeatherDemo.WORKING];
    }

    private static synchronized WeatherDemo.Grid resolution() {
        if (resolution == null) {
            resolution = WeatherDemo.grid(working());
        }
        return resolution;
    }

    // -------------------------------------------------------- the model itself

    @Test
    public void testTheSeasonalBlocksAreRotations() {
        DMatrix f = WeatherDemo.transition(WeatherDemo.WORKING, new double[] { 0.0, 0.0, 0.0, 0.0 });
        int[] blocks = { WeatherDemo.ANNUAL, WeatherDemo.SEMIANNUAL, WeatherDemo.HALF_ANNUAL,
                WeatherDemo.HALF_SEMIANNUAL };
        for (int i = 0; i < blocks.length; ++i) {
            int at = blocks[i];
            double c = f.get(at, at);
            double s = f.get(at, at + 1);
            assertEquals("block at " + at + " must preserve length", 1.0, c * c + s * s, 1.0e-14);
            assertEquals("block at " + at + " must be a rotation", c, f.get(at + 1, at + 1), 1.0e-15);
            assertEquals("block at " + at + " must be a rotation", -s, f.get(at + 1, at), 1.0e-15);
        }
    }

    @Test
    public void testTheAnnualBlockTurnsOnceAYear() {
        // the one property of the seasonal part that needs no data at all: the
        // angle has to close after exactly one period, or the component is not
        // seasonal and the whole of step 2 means nothing
        DMatrix f = WeatherDemo.transition(0, new double[] { 0.0, 0.0 });
        double angle = Math.atan2(f.get(WeatherDemo.ANNUAL, WeatherDemo.ANNUAL + 1),
                f.get(WeatherDemo.ANNUAL, WeatherDemo.ANNUAL));
        assertEquals("one turn per period", 2.0 * Math.PI, angle * WeatherDemo.PERIOD, 1.0e-12);
        double half = Math.atan2(f.get(WeatherDemo.SEMIANNUAL, WeatherDemo.SEMIANNUAL + 1),
                f.get(WeatherDemo.SEMIANNUAL, WeatherDemo.SEMIANNUAL));
        assertEquals("two turns per period", 4.0 * Math.PI, half * WeatherDemo.PERIOD, 1.0e-12);
    }

    @Test
    public void testTheTwoRowsOfHDifferOnlyWhereTheSpreadDoes() {
        // the structural claim of the whole demo: the two readings agree on the
        // temperature they share and disagree in the sign of the spread
        DMatrix h = WeatherDemo.observation(WeatherDemo.WORKING);
        int[] shared = { WeatherDemo.LEVEL, WeatherDemo.ANNUAL, WeatherDemo.SEMIANNUAL,
                WeatherDemo.ANOMALY };
        for (int i = 0; i < shared.length; ++i) {
            assertEquals("shared at " + shared[i], h.get(0, shared[i]), h.get(1, shared[i]), 0.0);
        }
        int[] opposed = { WeatherDemo.HALF_RANGE, WeatherDemo.HALF_ANNUAL,
                WeatherDemo.HALF_SEMIANNUAL };
        for (int i = 0; i < opposed.length; ++i) {
            assertEquals("opposed at " + opposed[i], h.get(0, opposed[i]), -h.get(1, opposed[i]), 0.0);
        }
    }

    @Test
    public void testTheProcessNoiseIsSingularAtEveryRung() {
        // one year is one seasonal cycle, so the seasonal part must not be
        // allowed to wander. That makes Q rank one at best, which
        // LinearGaussianModel permits and this demo depends on
        for (int rung = 0; rung < WeatherDemo.RUNGS; ++rung) {
            DMatrix q = WeatherDemo.processNoise(rung, new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 });
            int nonZero = 0;
            for (int i = 0; i < q.numRows(); ++i) {
                for (int j = 0; j < q.numColumns(); ++j) {
                    if (q.get(i, j) != 0.0) {
                        ++nonZero;
                    }
                }
            }
            assertTrue("rung " + rung + " has " + nonZero + " non-zero entries in Q",
                    nonZero <= 1);
            assertTrue("Q must be smaller than full rank", 1 < WeatherDemo.dimension(rung));
        }
    }

    @Test
    public void testTheObservationNoiseIsPositiveDefiniteWhereverTheSearchGoes() {
        // R is parameterized as L L', which is the reason the fit needs no
        // constraint on it. If that ever stopped holding, the fit would wander
        // into covariances that are not ones and the box would not catch it
        double[] scale = { -30.0, -5.0, 0.0, 5.0, 30.0 };
        double[] offDiagonal = { -1.0e4, -7.0, 0.0, 7.0, 1.0e4 };
        for (int a = 0; a < scale.length; ++a) {
            for (int b = 0; b < scale.length; ++b) {
                for (int c = 0; c < offDiagonal.length; ++c) {
                    double[] p = { 0.0, 0.0, scale[a], scale[b], offDiagonal[c] };
                    DMatrix r = WeatherDemo.observationNoise(WeatherDemo.UNBOUNDED, p);
                    double determinant = r.get(0, 0) * r.get(1, 1) - r.get(0, 1) * r.get(1, 0);
                    assertTrue("diagonal must be positive at " + scale[a] + ", " + scale[b],
                            r.get(0, 0) > 0.0 && r.get(1, 1) > 0.0);
                    assertTrue("determinant must not be negative", determinant >= 0.0);
                }
            }
        }
    }

    @Test
    public void testThePersistenceStaysInsideTheUnitInterval() {
        double[] u = { -1.0e6, -40.0, -1.0, 0.0, 1.0, 40.0, 1.0e6 };
        for (int i = 0; i < u.length; ++i) {
            double phi = WeatherDemo.persistence(u[i]);
            assertTrue("at u = " + u[i] + " the coefficient is " + phi, phi >= 0.0 && phi <= 1.0);
        }
        assertEquals("the middle of the range is a half", 0.5, WeatherDemo.persistence(0.0), 1.0e-15);
    }

    @Test
    public void testTheAnomalyExistsOnlyFromTheThirdRung() {
        assertEquals(10, WeatherDemo.dimension(0));
        assertEquals(10, WeatherDemo.dimension(1));
        assertEquals(11, WeatherDemo.dimension(WeatherDemo.WORKING));
        assertEquals(11, WeatherDemo.dimension(WeatherDemo.UNBOUNDED));
        assertEquals("a rung with no anomaly has no persistence", 0.0, ladder()[0].phi(), 0.0);
        assertEquals("nor a shock", 0.0, ladder()[1].sigmaAnomaly(), 0.0);
    }

    // ----------------------------------------------- the least squares oracle

    /** Rung 0 with equal noise on the two readings and a prior of the given width. */
    private static LinearGaussianModel flat(double priorVariance) {
        DMatrix r = new DMatrix(2, 2);
        r.set(0, 0, 9.0);
        r.set(1, 1, 9.0);
        int d = WeatherDemo.dimension(0);
        DMatrix p0 = new DMatrix(d, d);
        for (int i = 0; i < d; ++i) {
            p0.set(i, i, priorVariance);
        }
        double[] zero = { 0.0, 0.0 };
        return new LinearGaussianModel(WeatherDemo.transition(0, zero),
                WeatherDemo.processNoise(0, zero), WeatherDemo.observation(0), r, new double[d], p0);
    }

    private static double worstAgainstLeastSquares(double priorVariance, LSSummary ols) {
        RtsSmoother.Result s = RtsSmoother.smooth(flat(priorVariance), WeatherDemo.observations());
        DMatrix h = WeatherDemo.observation(0);
        int n = s.length;
        double worst = 0.0;
        for (int t = 0; t < n; ++t) {
            for (int c = 0; c < 2; ++c) {
                double implied = 0.0;
                for (int j = 0; j < WeatherDemo.dimension(0); ++j) {
                    implied += h.get(c, j) * s.mean[t][j];
                }
                worst = Math.max(worst, Math.abs(implied - ols.getYHat().get(c == 0 ? t : n + t)));
            }
        }
        return worst;
    }

    @Test
    public void testAFilterWithNoProcessNoiseIsLeastSquares() {
        // the demo's claim in step 1 and the reason rung 0 is a model rather
        // than a straw man. Compare the fitted values and not the state: the
        // seasonal components rotate, so the coefficients live in a frame that
        // turns with the year and the last state is not the OLS coefficient
        LSSummary ols = WeatherDemo.ordinaryLeastSquares();
        assertEquals("the design has ten terms", WeatherDemo.TERMS, ols.getCoefficientsCount());
        assertEquals("measured 5.2e-07", 0.0, worstAgainstLeastSquares(1.0e6, ols), 1.0e-5);
    }

    @Test
    public void testWhatIsLeftBetweenThemIsTheFinitePrior() {
        // and the proof that it is the prior and nothing else: the gap is the
        // ridge a proper prior imposes, so it has to fall like 1 / P0. Measured
        // 5.135e-04, 5.135e-05, 5.135e-06 at 1e3, 1e4, 1e5
        LSSummary ols = WeatherDemo.ordinaryLeastSquares();
        double coarse = worstAgainstLeastSquares(1.0e3, ols);
        double fine = worstAgainstLeastSquares(1.0e4, ols);
        double finer = worstAgainstLeastSquares(1.0e5, ols);
        assertEquals("ten times the prior, a tenth of the gap", 10.0, coarse / fine, 0.5);
        assertEquals("and again", 10.0, fine / finer, 0.5);
    }

    // --------------------------------------------------------------- the ladder

    @Test
    public void testEveryRungImprovesTheLikelihood() {
        WeatherDemo.Fit[] fits = ladder();
        for (int i = 1; i < fits.length; ++i) {
            assertTrue("rung " + i + " must not be worse than rung " + (i - 1),
                    fits[i].logLikelihood > fits[i - 1].logLikelihood);
            assertTrue("rung " + i + " converged", fits[i].converged);
        }
        assertEquals("M0, measured", -2032.17, fits[0].logLikelihood, 0.05);
        assertEquals("M1, measured", -1878.90, fits[1].logLikelihood, 0.05);
        assertEquals("M2, measured", -1788.90, fits[WeatherDemo.WORKING].logLikelihood, 0.05);
        assertEquals("M3, measured", -1787.18, fits[WeatherDemo.UNBOUNDED].logLikelihood, 0.05);
    }

    @Test
    public void testTheOffDiagonalMeansOneThingAtM1AndAnotherAtM3() {
        // the centre of the demo. The same parameter, added to two different
        // models, and worth two orders of magnitude apart
        WeatherDemo.Fit[] fits = ladder();
        double atM1 = fits[1].logLikelihood - fits[0].logLikelihood;
        double atM3 = fits[WeatherDemo.UNBOUNDED].logLikelihood - fits[WeatherDemo.WORKING].logLikelihood;
        assertEquals("measured 153.27", 153.27, atM1, 0.1);
        assertEquals("measured 1.72", 1.72, atM3, 0.1);
        assertTrue("it has to be worth far less once a state carries the signal",
                atM1 > 50.0 * atM3);
    }

    @Test
    public void testTheOffDiagonalChangesSignOnTheWayUp() {
        WeatherDemo.Fit[] fits = ladder();
        assertEquals("measured +0.7571", 0.7571, fits[1].noiseCorrelation(), 0.001);
        assertTrue("at M1 the two readings look alike", fits[1].noiseCorrelation() > 0.5);
        assertTrue("at M3 they look opposed",
                fits[WeatherDemo.UNBOUNDED].noiseCorrelation() < -0.5);
    }

    @Test
    public void testTheTopRungSitsOnTheBoundaryOfItsParameterSpace() {
        // the reason the demo reports M3 and works with M2. A correlation of
        // exactly -1 is not an estimate, it is the edge of the set
        WeatherDemo.Fit top = ladder()[WeatherDemo.UNBOUNDED];
        assertEquals("the correlation is the boundary", -1.0, top.noiseCorrelation(), 1.0e-5);
        DMatrix r = WeatherDemo.observationNoise(top.rung, top.parameters);
        double determinant = r.get(0, 0) * r.get(1, 1) - r.get(0, 1) * r.get(1, 0);
        assertEquals("so R is singular", 0.0, determinant / (r.get(0, 0) * r.get(1, 1)), 1.0e-5);
    }

    @Test
    public void testAkaikePrefersTheModelOnTheBoundary() {
        // stated so that it stays uncomfortable: the criterion picks the model
        // the demo declines to use, and by a margin small enough to see
        WeatherDemo.Fit[] fits = ladder();
        assertTrue("AIC really does prefer M3",
                fits[WeatherDemo.UNBOUNDED].aic() < fits[WeatherDemo.WORKING].aic());
        assertEquals("measured 1.44", 1.44,
                fits[WeatherDemo.WORKING].aic() - fits[WeatherDemo.UNBOUNDED].aic(), 0.1);
    }

    @Test
    public void testTheMinimumIsMeasuredFarMoreSharplyThanTheMaximum() {
        WeatherDemo.Fit fit = working();
        assertEquals("measured 2.7392", 2.7392, fit.sigmaMaximum(), 0.005);
        assertEquals("measured 0.4255", 0.4255, fit.sigmaMinimum(), 0.005);
        assertTrue("the asymmetry is what steps 3 and 6 rest on",
                fit.sigmaMaximum() > 5.0 * fit.sigmaMinimum());
    }

    // -------------------------------------------------------- the decomposition

    @Test
    public void testTheDecompositionReproducesTheSpreadExactly() {
        // the two fitted readings are built from the same three parts with the
        // half range entering twice with opposite signs, so their difference is
        // that component doubled and nothing else. Bit for bit
        WeatherDemo.Decomposition parts = WeatherDemo.decompose(working());
        for (int t = 0; t < parts.maximum.length; ++t) {
            assertEquals("at day " + t, 2.0 * parts.halfRange[t], parts.maximum[t] - parts.minimum[t],
                    1.0e-12);
        }
    }

    @Test
    public void testTheLargestAnomalyIsTheColdestNightOfTheYear() {
        WeatherDemo.Decomposition parts = WeatherDemo.decompose(working());
        int[] tmin = Datasets.tminTenths();
        int coldest = 0;
        for (int t = 1; t < tmin.length; ++t) {
            if (tmin[t] < tmin[coldest]) {
                coldest = t;
            }
        }
        assertEquals("the demo says so in step 2", coldest, parts.largestAnomalyDay);
        assertEquals("2025-01-22", Datasets.label(parts.largestAnomalyDay));
        assertTrue("and it is a cold anomaly", parts.largestAnomaly < 0.0);
    }

    @Test
    public void testTheFittedPersistenceExceedsWhatTheResidualsShow() {
        // an autocorrelation measured through observation noise is diluted by
        // it. If this ever came out the other way round the model would be
        // taking weather for noise rather than the reverse
        double[] residuals = WeatherDemo.deterministicResiduals();
        double acf1 = math.probe.ACF.acf(residuals, 1)[1];
        assertEquals("measured 0.5799", 0.5799, acf1, 0.002);
        assertTrue("the model must find more memory than the residuals show",
                working().phi() > acf1);
        assertEquals("measured 0.6398", 0.6398, working().phi(), 0.002);
    }

    // ----------------------------------------------------------- the innovations

    @Test
    public void testTheInnovationsAreWhiterThanTheResiduals() {
        WeatherDemo.Whiteness white = WeatherDemo.whiteness(working());
        for (int c = 0; c < 2; ++c) {
            assertTrue("reading " + c + " must lose memory, not gain it",
                    Math.abs(white.innovationAcf1[c]) < Math.abs(white.deterministicAcf1[c]));
            assertTrue("and lose lags outside the band",
                    white.innovationOutside[c] <= white.deterministicOutside[c]);
        }
        assertEquals("the minimum comes out white, measured 0.0246", 0.0246, white.innovationAcf1[1],
                0.005);
        assertTrue("the maximum does not quite, measured 0.2157", white.innovationAcf1[0] > 0.15);
    }

    @Test
    public void testTheInnovationsAreStandardizedAsClaimed() {
        WeatherDemo.Whiteness white = WeatherDemo.whiteness(working());
        double error = 1.0 / Math.sqrt(Datasets.size());
        for (int c = 0; c < 2; ++c) {
            assertEquals("mean of reading " + c, 0.0, white.innovationMean[c], 4.0 * error);
            assertEquals("standard deviation of reading " + c, 1.0, white.innovationSd[c], 0.05);
        }
    }

    // ---------------------------------------------------------------- the grid

    @Test
    public void testRoundingIsASeventhOfTheMinimumAndNothingOfTheMaximum() {
        WeatherDemo.Grid grid = resolution();
        assertEquals("one degree Fahrenheit over root twelve", 0.1604, grid.quantizationSd, 1.0e-4);
        assertEquals("measured 14.2 percent", 0.1420, grid.shareOfMinimum, 0.005);
        assertTrue("and next to nothing of the maximum", grid.shareOfMaximum < 0.01);
    }

    @Test
    public void testTheAdditiveQuantizationRuleMisses() {
        // the finding of step 4, kept as a test so that it stays a measurement
        // and not an anecdote: d^2/12 assumes a rounding error independent of
        // what it rounds, and these columns were already on a lattice
        WeatherDemo.Grid grid = resolution();
        assertTrue("the coarsening is not what d^2/12 says it is",
                Math.abs(grid.measuredCoarseSd - grid.coarseSd) > 0.02);
        assertTrue("and the refit lands outside both predictions",
                grid.refittedMinimum > grid.predictedMinimum
                        && grid.refittedMinimum > grid.predictedFromMeasured);
        assertEquals("measured 0.5608", 0.5608, grid.refittedMinimum, 0.01);
        assertEquals("the maximum barely moves", working().sigmaMaximum(), grid.refittedMaximum, 0.01);
    }

    // ---------------------------------------------------------- the trace days

    @Test
    public void testATraceDayIsADryDay() {
        WeatherDemo.Trace kinds = WeatherDemo.trace(WeatherDemo.decompose(working()));
        assertTrue("a trace is indistinguishable from a dry day",
                Math.abs(kinds.traceAgainstDry) < 1.0);
        assertTrue("and does differ from a wet one, though not by much",
                Math.abs(kinds.traceAgainstWet) > 1.5);
        assertTrue("a wet day has the narrower spread", kinds.meanSpread[0] < kinds.meanSpread[2]);
    }

    @Test
    public void testTheThreeKindsOfDayCoverTheYear() {
        WeatherDemo.Trace kinds = WeatherDemo.trace(WeatherDemo.decompose(working()));
        assertEquals("measurable rain", 121, kinds.count[0]);
        assertEquals("a trace", 35, kinds.count[1]);
        assertEquals("dry", 209, kinds.count[2]);
        assertEquals("and nothing else", Datasets.size(),
                kinds.count[0] + kinds.count[1] + kinds.count[2]);
    }

    // ------------------------------------------------------------ the readings

    @Test
    public void testTheMinimumAloneIsWorthNearlyBothReadings() {
        WeatherDemo.Channels pair = WeatherDemo.channels(working());
        assertTrue("both readings must beat either alone", pair.fromBoth <= pair.fromMinimum);
        assertTrue("the minimum alone is nearly as good", pair.fromMinimum < 1.05 * pair.fromBoth);
        assertTrue("the maximum alone is not", pair.fromMaximum > 2.0 * pair.fromBoth);
        assertEquals("measured 0.9121", 0.9121, pair.fromBoth, 0.005);
    }

    @Test
    public void testADiagonalNoiseLeavesNothingForTheMissingTerm() {
        // the exact form of the demo's point in step 7: the term is
        // proportional to R10, so a diagonal R makes it identically zero and
        // the reading the filter implies is already the conditional mean
        WeatherDemo.Blanked gap = WeatherDemo.blank(working(), WeatherDemo.BLANK_FROM,
                WeatherDemo.BLANK_LENGTH);
        assertEquals("exactly zero, not nearly", 0.0, gap.largestTerm, 0.0);
        assertEquals("so the two answers coincide", gap.plainError, gap.correctedError, 0.0);
        WeatherDemo.Blanked correlated = WeatherDemo.blank(ladder()[WeatherDemo.UNBOUNDED],
                WeatherDemo.BLANK_FROM, WeatherDemo.BLANK_LENGTH);
        assertTrue("and a correlated R makes it something", correlated.largestTerm > 0.05);
    }

    @Test
    public void testBlankingOneReadingLeavesTheDaysBeforeItUntouched() {
        // a missing value must not reach backwards. The filter is a forward
        // recursion, so this has to hold bit for bit, and if it ever did not
        // every blanking experiment in the demo would be meaningless
        LinearGaussianModel m = working().model();
        DMatrix full = WeatherDemo.observations();
        DMatrix holed = WeatherDemo.observations();
        for (int t = WeatherDemo.BLANK_FROM; t < WeatherDemo.BLANK_FROM + WeatherDemo.BLANK_LENGTH; ++t) {
            holed.set(t, 1, Double.NaN);
        }
        KalmanFilter.Result a = KalmanFilter.filter(m, full);
        KalmanFilter.Result b = KalmanFilter.filter(m, holed);
        for (int t = 0; t < WeatherDemo.BLANK_FROM; ++t) {
            for (int j = 0; j < WeatherDemo.dimension(working().rung); ++j) {
                assertEquals("day " + t + ", component " + j, a.filteredMean[t][j],
                        b.filteredMean[t][j], 0.0);
            }
        }
        assertEquals("and the day itself is where they part", 2, a.observedComponents[WeatherDemo.BLANK_FROM]);
        assertEquals("with one component gone", 1, b.observedComponents[WeatherDemo.BLANK_FROM]);
    }

    @Test
    public void testTheChangeOfSourceLeavesNoOffset() {
        WeatherDemo.Break change = WeatherDemo.sourceBreak(working());
        for (int c = 0; c < 2; ++c) {
            assertTrue("reading " + c + " shifted by " + change.ratio[c] + " standard errors",
                    Math.abs(change.ratio[c]) < 2.0);
        }
    }

    // --------------------------------------------------------------- the data

    @Test
    public void testTheChecksumsAreWhatTheClassCommentSays() {
        assertEquals(365, Datasets.size());
        assertEquals(6229.8, Datasets.checksumTmax(), 1.0e-9);
        assertEquals(3397.3, Datasets.checksumTmin(), 1.0e-9);
        assertEquals(1006.1, Datasets.checksumPrecipitation(), 1.0e-9);
    }

    @Test
    public void testEveryTraceDayIsAZeroInTheNumberColumn() {
        // the point of keeping the flag at all: these days are not visible in
        // the number, so a reader who never looks at the flag cannot find them
        int[] trace = Datasets.trace();
        int[] prcp = Datasets.precipitationTenths();
        assertEquals(35, trace.length);
        for (int i = 0; i < trace.length; ++i) {
            assertEquals("day " + trace[i] + " should read zero", 0, prcp[trace[i]]);
            assertTrue(Datasets.isTrace(trace[i]));
            assertFalse(Datasets.isWet(trace[i]));
            if (i > 0) {
                assertTrue("the indices are ascending", trace[i] > trace[i - 1]);
            }
        }
    }

    @Test
    public void testTheReadingsAfterTheSourceBreakLieExactlyOnTheFahrenheitGrid() {
        // the claim in the Datasets comment, as a measurement. The CF6 form is
        // printed in whole degrees Fahrenheit, and 254 of 254 values say so
        int[] tmax = Datasets.tmaxTenths();
        int[] tmin = Datasets.tminTenths();
        int on = 0;
        int total = 0;
        for (int t = Datasets.SOURCE_BREAK; t < Datasets.size(); ++t) {
            for (int c = 0; c < 2; ++c) {
                int v = c == 0 ? tmax[t] : tmin[t];
                ++total;
                for (int f = -40; f <= 120; ++f) {
                    double tenths = (f - 32) * Datasets.FAHRENHEIT_STEP;
                    if (Math.round(tenths) == v || (long) tenths == v) {
                        ++on;
                        break;
                    }
                }
            }
        }
        assertEquals("the stretch is 127 days of two readings", 254, total);
        assertEquals("and every one of them is a whole degree Fahrenheit", 254, on);
    }

    @Test
    public void testTheMinimumNeverExceedsTheMaximum() {
        int[] tmax = Datasets.tmaxTenths();
        int[] tmin = Datasets.tminTenths();
        for (int t = 0; t < tmax.length; ++t) {
            assertTrue("at " + Datasets.label(t), tmin[t] < tmax[t]);
        }
    }

    // --------------------------------------------------------------- the page

    private static synchronized String root() throws UnsupportedEncodingException {
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
        // every printf carries Locale.ROOT, and this is what says so. A German
        // default would otherwise turn every decimal point into a comma
        String comma;
        Locale saved = Locale.getDefault();
        try {
            Locale.setDefault(Locale.GERMANY);
            comma = runMain();
        } finally {
            Locale.setDefault(saved);
        }
        assertEquals("the page must not depend on the default locale", root(), comma);
    }

    @Test
    public void testTheDemoStillPrintsItsEightSteps() throws UnsupportedEncodingException {
        String text = root();
        for (int step = 1; step <= 8; ++step) {
            assertTrue("step " + step + " is missing", text.contains("=== " + step + "."));
        }
        assertTrue("the station has to be named", text.contains(Datasets.STATION));
        assertTrue("and the ladder printed", text.contains("+ linked errors again"));
        assertTrue("the page has to lead with the two readings",
                text.indexOf("one temperature, read twice") < text.indexOf("how the model was chosen"));
    }

    @Test
    public void testThePageCarriesNoIterationCount() throws UnsupportedEncodingException {
        // an iteration count is a property of the route and not of the answer.
        // The other state space demo printed one and the page then differed
        // between JDK 8 and JDK 25, because Math.exp is allowed an ulp and one
        // ulp flips a simplex comparison
        String text = root().toLowerCase(Locale.ROOT);
        assertFalse("an iteration count would make the page runtime dependent",
                text.contains("iteration"));
    }

    private static String runMain() throws UnsupportedEncodingException {
        PrintStream saved = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, "UTF-8"));
            WeatherDemo.main(new String[0]);
        } finally {
            System.setOut(saved);
        }
        String text = captured.toString("UTF-8");
        assertNotNull(text);
        return text;
    }
}
