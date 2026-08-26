package math.stats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.Locale;

import org.junit.Test;

import math.cern.Arithmetic;
import math.distribution.ContinuousDistribution;
import math.distribution.Exponential;
import math.distribution.Hypergeometric;
import math.distribution.Normal;
import math.distribution.StudentT;
import math.distribution.Uniform;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.Wls;
import math.list.DoubleArrayList;
import math.stats.gof.AndersonDarling;
import math.stats.gof.CramerVonMises;
import math.stats.gof.DurbinWatson;
import math.stats.gof.KolmogorovSmirnov;
import math.stats.gof.KolmogorovSmirnovTwoSample;
import math.stats.gof.Lilliefors;
import math.stats.gof.WeightedChiSquare;
import math.stats.mle.MLE;
import math.stats.mle.ParNormal;

/**
 * A p-value is not checked by remembering one. It is checked by asking whether
 * it is uniform when the null hypothesis holds, whether its tails add up, and
 * whether the confidence interval it comes with excludes exactly what it says
 * it excludes.
 */
public final class HypothesisTestsTest {

    /** Deterministic uniforms in {@code (0,1)}. */
    private static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) + 0.5) * 0x1.0p-53;
        }
    }

    private static double[] sample(ContinuousDistribution d, int n, long seed) {
        Lcg lcg = new Lcg(seed);
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = d.inverseCdf(lcg.next());
        }
        return x;
    }

    // ------------------------------------------------------------ the p-value --

    @Test
    public void testThePValueIsUniformUnderTheNullHypothesis() {
        // this is what actually establishes a p-value: tested against the mean
        // the data really has, it must be uniform, so the proportion below any
        // threshold is that threshold
        for (int n : new int[] { 3, 10, 50 }) {
            int reps = 8000;
            int below01 = 0;
            int below05 = 0;
            int below50 = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double p = HypothesisTests.tOneSample(x, 3.0, Alternative.TWO_SIDED, 0.95).test.pValue;
                assertTrue("p = " + p, p >= 0.0 && p <= 1.0);
                if (p < 0.01) {
                    below01++;
                }
                if (p < 0.05) {
                    below05++;
                }
                if (p < 0.50) {
                    below50++;
                }
            }
            // measured over 20000 replications: 0.0101 .. 0.0102, 0.0497 ..
            // 0.0525 and 0.5009 .. 0.5054
            assertRate("n=" + n + " P(p<0.01)", below01 / (double) reps, 0.01, 0.006);
            assertRate("n=" + n + " P(p<0.05)", below05 / (double) reps, 0.05, 0.015);
            assertRate("n=" + n + " P(p<0.50)", below50 / (double) reps, 0.50, 0.040);
        }
    }

    private static void assertRate(String what, double got, double want, double tolerance) {
        assertTrue(what + " = " + got + ", wanted " + want + " within " + tolerance,
                Math.abs(got - want) <= tolerance);
    }

    @Test
    public void testTheOneSidedTailsSplitTheTwoSidedOne() {
        // the two one-sided p-values partition the probability, and twice the
        // smaller of them is the two-sided one
        for (int n : new int[] { 3, 10, 50 }) {
            for (long seed = 1; seed <= 200; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double less = HypothesisTests.tOneSample(x, 3.0, Alternative.LESS, 0.95).test.pValue;
                double greater = HypothesisTests.tOneSample(x, 3.0, Alternative.GREATER, 0.95).test.pValue;
                double twoSided = HypothesisTests.tOneSample(x, 3.0, Alternative.TWO_SIDED, 0.95).test.pValue;
                String at = "n=" + n + " seed=" + seed;
                // measured as exactly zero, but Math.exp differs by one unit in
                // the last place between JDK 8 and JDK 25, so this is not the
                // place for a bit-for-bit assertion
                assertEquals(at + ": the tails do not partition", 1.0, less + greater, 1.0e-15);
                assertEquals(at + ": twice the smaller tail is not the two-sided value", twoSided,
                        2.0 * Math.min(less, greater), 1.0e-15);
            }
        }
    }

    @Test
    public void testTheStatisticAndItsPValueAreTheDefinition() {
        // the same arithmetic, so nothing may differ
        for (int n : new int[] { 3, 10, 200 }) {
            for (long seed = 1; seed <= 60; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                TTestResult r = HypothesisTests.tOneSample(x, 2.5, Alternative.TWO_SIDED, 0.95);
                double[] moments = HypothesisTests.meanAndDeviation(x);
                double standardError = moments[1] / Math.sqrt(n);
                double t = (moments[0] - 2.5) / standardError;
                String at = "n=" + n + " seed=" + seed;
                assertEquals(at + ": estimate", moments[0], r.estimate, 0.0);
                assertEquals(at + ": standard error", standardError, r.standardError, 0.0);
                assertEquals(at + ": statistic", t, r.test.statistic, 0.0);
                assertEquals(at + ": degrees of freedom", n - 1.0, r.test.degreesOfFreedom, 0.0);
                assertEquals(at + ": p-value", 2.0 * new StudentT(n - 1.0).cdf(-Math.abs(t)), r.test.pValue, 0.0);
            }
        }
    }

    // ---------------------------------------------------------- the interval --

    @Test
    public void testTheConfidenceIntervalAgreesWithTheDecision() {
        // the interval and the decision are two views of one computation: the
        // interval excludes the null value exactly when the test rejects. It
        // holds for the one-sided alternatives too, which is why the interval
        // is built with one infinite end for those
        for (double level : new double[] { 0.90, 0.95, 0.99 }) {
            double alpha = 1.0 - level;
            for (Alternative alternative : Alternative.values()) {
                for (int n : new int[] { 3, 10, 50 }) {
                    for (long seed = 1; seed <= 200; seed++) {
                        // walk the null value across the sample so the boundary
                        // is approached from both sides
                        double mu = 3.0 + 0.6 * (seed % 11 - 5);
                        double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                        TTestResult r = HypothesisTests.tOneSample(x, mu, alternative, level);
                        boolean rejects = r.test.rejectsAt(alpha);
                        boolean outside = mu < r.lower || mu > r.upper;
                        assertEquals(alternative + " level=" + level + " n=" + n + " seed=" + seed + " mu=" + mu
                                + " p=" + r.test.pValue + " interval=[" + r.lower + ", " + r.upper + "]", rejects,
                                outside);
                    }
                }
            }
        }
    }

    @Test
    public void testAOneSidedIntervalIsOpenAtOneEnd() {
        double[] x = sample(new Normal(3.0, 2.0), 20, 11L);
        TTestResult less = HypothesisTests.tOneSample(x, 3.0, Alternative.LESS, 0.95);
        assertEquals("lower", Double.NEGATIVE_INFINITY, less.lower, 0.0);
        assertTrue("upper", less.upper > less.estimate);
        TTestResult greater = HypothesisTests.tOneSample(x, 3.0, Alternative.GREATER, 0.95);
        assertTrue("lower", greater.lower < greater.estimate);
        assertEquals("upper", Double.POSITIVE_INFINITY, greater.upper, 0.0);
        TTestResult both = HypothesisTests.tOneSample(x, 3.0, Alternative.TWO_SIDED, 0.95);
        assertEquals("the two-sided interval is not centred on the estimate",
                both.upper - both.estimate, both.estimate - both.lower, 1.0e-12);
    }

    // ----------------------------------------------------------- the moments --

    @Test
    public void testTheDeviationAgreesWithTheListImplementation() {
        // math.list.DoubleArrayList.stddev is the same corrected two-pass
        // without Bessel's correction, and it is tested in its own package, so
        // it is the independent algorithm to check against
        double worst = 0.0;
        String at = "";
        for (double centre : new double[] { 0.0, 1.0, 1.0e8, 1.0e14 }) {
            double[] x = new double[16];
            for (int i = 0; i < x.length; i++) {
                x[i] = centre + (i - 7.5);
            }
            double rel = deviationDeviation(x);
            if (rel > worst) {
                worst = rel;
                at = "centre " + centre;
            }
        }
        // a variance would have overflowed at 1e300 and underflowed at 1e-320;
        // a deviation, and the scale invariant statistic built from it, do not
        for (double centre : new double[] { 1.0e300, 1.0e-320 }) {
            double[] x = new double[16];
            for (int i = 0; i < x.length; i++) {
                x[i] = centre * (1.0 + 0.01 * (i - 7.5));
            }
            double[] moments = HypothesisTests.meanAndDeviation(x);
            assertTrue("centre " + centre + " gave a deviation of " + moments[1],
                    moments[1] > 0.0 && moments[1] < Double.POSITIVE_INFINITY);
            double rel = deviationDeviation(x);
            if (rel > worst) {
                worst = rel;
                at = "centre " + centre;
            }
        }
        double rel = deviationDeviation(sample(new Normal(3.0, 2.0), 200, 99L));
        if (rel > worst) {
            worst = rel;
            at = "normal(3, 2), n = 200";
        }
        // measured worst: 1.9e-16
        assertTrue("worst relative deviation " + worst + " at " + at, worst < 1.0e-14);
    }

    /** Relative distance between our deviation and the one the list computes. */
    private static double deviationDeviation(double[] x) {
        double ours = HypothesisTests.meanAndDeviation(x)[1];
        double population = new DoubleArrayList(x).stddev();
        double theirs = population * Math.sqrt(x.length / (x.length - 1.0));
        return Math.abs(ours - theirs) / theirs;
    }

    @Test
    public void testTheMeanSurvivesDataThatWouldOverflowItsSum() {
        // sixteen values at 1e300 sum to 1.6e301, which is fine; sixteen at
        // 1e307 do not, and the mean is still an ordinary number
        double[] x = new double[16];
        for (int i = 0; i < x.length; i++) {
            x[i] = 1.0e307 * (1.0 + 0.01 * (i - 7.5));
        }
        double[] moments = HypothesisTests.meanAndDeviation(x);
        assertEquals("mean", 1.0e307, moments[0], 1.0e295);
        assertTrue("deviation " + moments[1], moments[1] > 0.0 && moments[1] < Double.POSITIVE_INFINITY);
    }

    @Test
    public void testTheStatisticIsScaleInvariant() {
        // t is a ratio of two things that scale alike, so scaling the sample
        // and the null value together may not move it
        double[] x = sample(new Normal(3.0, 2.0), 40, 5L);
        double reference = HypothesisTests.tOneSample(x, 2.5, Alternative.TWO_SIDED, 0.95).test.statistic;
        for (double factor : new double[] { 0x1.0p-900, 0x1.0p-40, 0x1.0p40, 0x1.0p900 }) {
            double[] scaled = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                scaled[i] = x[i] * factor;
            }
            double got = HypothesisTests.tOneSample(scaled, 2.5 * factor, Alternative.TWO_SIDED, 0.95).test.statistic;
            assertEquals("factor " + factor, reference, got, 1.0e-12 * Math.abs(reference));
        }
    }

    // ------------------------------------------------------------ rejectsAt --

    @Test
    public void testRejectsAtIsConsistentWithThePValue() {
        double[] x = sample(new Normal(3.0, 2.0), 25, 17L);
        TestResult r = HypothesisTests.tOneSample(x, 2.0, Alternative.TWO_SIDED, 0.95).test;
        for (double alpha : new double[] { 0.001, 0.01, 0.05, 0.10, 0.5, 0.9 }) {
            assertEquals("alpha=" + alpha, r.pValue < alpha, r.rejectsAt(alpha));
        }
    }

    @Test
    public void testRejectsAtRejectsASignificanceLevelOutsideTheUnitInterval() {
        TestResult r = HypothesisTests.tOneSample(sample(new Normal(3.0, 2.0), 25, 17L), 2.0, Alternative.TWO_SIDED,
                0.95).test;
        for (double alpha : new double[] { 0.0, 1.0, -0.1, 1.5, Double.NaN }) {
            try {
                r.rejectsAt(alpha);
                fail("alpha = " + alpha + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage(), expected.getMessage().contains("alpha"));
            }
        }
    }

    // -------------------------------------------------------------- printing --

    @Test
    public void testTheOutputDoesNotDependOnTheLocale() {
        // math.OutputLocaleTest states the rule; this is the same invariant for
        // the two new value types
        Locale before = Locale.getDefault();
        try {
            double[] x = sample(new Normal(3.4, 1.5), 12, 7L);
            TTestResult r = HypothesisTests.tOneSample(x, 3.0, Alternative.TWO_SIDED, 0.95);
            Locale.setDefault(Locale.ROOT);
            String root = r.toString();
            String rootTest = r.test.toString();
            Locale.setDefault(Locale.GERMANY);
            // the invariant itself, not a pattern match on the text, which is
            // how math.OutputLocaleTest puts it: the commas in this line are
            // separators between fields, and a decimal comma would show up as a
            // difference between the two strings
            assertEquals("TTestResult.toString", root, r.toString());
            assertEquals("TestResult.toString", rootTest, r.test.toString());
            Locale.setDefault(Locale.ROOT);
            SimulatedTestResult fitted = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL, 200, 3L);
            String rootFitted = fitted.toString();
            Locale.setDefault(Locale.GERMANY);
            assertEquals("SimulatedTestResult.toString", rootFitted, fitted.toString());
        } finally {
            Locale.setDefault(before);
        }
    }

    @Test
    public void testAWholeNumberOfDegreesOfFreedomPrintsAsOne() {
        String line = HypothesisTests.tOneSample(sample(new Normal(3.0, 2.0), 12, 7L), 3.0, Alternative.TWO_SIDED,
                0.95).test.toString();
        assertTrue("expected 'df = 11' in : " + line, line.contains("df = 11,") || line.contains("df = 11 "));
        assertTrue("the alternative is not named in : " + line, line.contains("two-sided"));
    }

    // --------------------------------------------------------- the guard rail --

    @Test
    public void testTheSampleIsChecked() {
        rejects("null", null, 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("an empty sample", new double[0], 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("a single observation", new double[] { 1.0 }, 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("a NaN", new double[] { 1.0, Double.NaN, 2.0 }, 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("an infinity", new double[] { 1.0, Double.POSITIVE_INFINITY, 2.0 }, 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("a constant sample", new double[] { 2.0, 2.0, 2.0, 2.0 }, 0.0, Alternative.TWO_SIDED, 0.95);
        rejects("a non-finite mu", new double[] { 1.0, 2.0, 3.0 }, Double.NaN, Alternative.TWO_SIDED, 0.95);
        rejects("a null alternative", new double[] { 1.0, 2.0, 3.0 }, 0.0, null, 0.95);
        rejects("a level of zero", new double[] { 1.0, 2.0, 3.0 }, 0.0, Alternative.TWO_SIDED, 0.0);
        rejects("a level of one", new double[] { 1.0, 2.0, 3.0 }, 0.0, Alternative.TWO_SIDED, 1.0);
    }

    private static void rejects(String what, double[] x, double mu, Alternative alternative, double level) {
        try {
            HypothesisTests.tOneSample(x, mu, alternative, level);
            fail("tOneSample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------- two samples --

    @Test
    public void testTheTwoSamplePValuesAreUniformUnderTheNullHypothesis() {
        for (int nx : new int[] { 4, 12, 40 }) {
            int ny = Math.max(3, nx / 2);
            int reps = 5000;
            int welch = 0;
            int pooled = 0;
            int paired = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] a = sample(new Normal(3.0, 2.0), nx, seed * 7919L + 1);
                double[] b = sample(new Normal(3.0, 2.0), ny, seed * 104729L + 3);
                double[] c = sample(new Normal(3.0, 2.0), nx, seed * 15485863L + 7);
                if (HypothesisTests.tTwoSample(a, b, Alternative.TWO_SIDED, 0.95).test.rejectsAt(0.05)) {
                    welch++;
                }
                if (HypothesisTests.tTwoSamplePooled(a, b, Alternative.TWO_SIDED, 0.95).test.rejectsAt(0.05)) {
                    pooled++;
                }
                if (HypothesisTests.tPaired(a, c, Alternative.TWO_SIDED, 0.95).test.rejectsAt(0.05)) {
                    paired++;
                }
            }
            // measured over 12000 replications with the variances equal:
            // Welch 0.0418 .. 0.0518, pooled 0.0483 .. 0.0515, paired
            // 0.0484 .. 0.0522. Welch is a little conservative at the
            // smallest sizes, which is its known behaviour
            assertRate("nx=" + nx + " Welch", welch / (double) reps, 0.05, 0.018);
            assertRate("nx=" + nx + " pooled", pooled / (double) reps, 0.05, 0.015);
            assertRate("nx=" + nx + " paired", paired / (double) reps, 0.05, 0.015);
        }
    }

    @Test
    public void testPooledLosesItsCalibrationWhereWelchKeepsIt() {
        // why Welch is the default, and why the pooled test carries its
        // assumption in its name. With the variances four to one and the
        // larger of them on the smaller sample, the pooled test rejects a
        // true null far more often than the level it was asked for
        int reps = 5000;
        int welch = 0;
        int pooled = 0;
        for (long seed = 1; seed <= reps; seed++) {
            double[] a = sample(new Normal(3.0, 1.0), 12, seed * 7919L + 1);
            double[] b = sample(new Normal(3.0, 4.0), 6, seed * 104729L + 3);
            if (HypothesisTests.tTwoSample(a, b, Alternative.TWO_SIDED, 0.95).test.rejectsAt(0.05)) {
                welch++;
            }
            if (HypothesisTests.tTwoSamplePooled(a, b, Alternative.TWO_SIDED, 0.95).test.rejectsAt(0.05)) {
                pooled++;
            }
        }
        double welchRate = welch / (double) reps;
        double pooledRate = pooled / (double) reps;
        // measured over 12000 replications: 0.0540 and 0.1729
        assertTrue("Welch is no longer calibrated either : " + welchRate, welchRate < 0.08);
        assertTrue("the pooled test did not lose its calibration : " + pooledRate, pooledRate > 0.10);
    }

    @Test
    public void testThePairedTestIsTheOneSampleTestOnTheDifferences() {
        // the same arithmetic on the same numbers, so every field agrees to
        // the last bit
        for (int n : new int[] { 3, 10, 100 }) {
            for (long seed = 1; seed <= 40; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double[] y = sample(new Normal(2.5, 1.5), n, seed * 104729L + 3);
                double[] differences = new double[n];
                for (int i = 0; i < n; i++) {
                    differences[i] = x[i] - y[i];
                }
                TTestResult paired = HypothesisTests.tPaired(x, y, Alternative.TWO_SIDED, 0.95);
                TTestResult single = HypothesisTests.tOneSample(differences, 0.0, Alternative.TWO_SIDED, 0.95);
                String at = "n=" + n + " seed=" + seed;
                assertEquals(at + ": statistic", single.test.statistic, paired.test.statistic, 0.0);
                assertEquals(at + ": p-value", single.test.pValue, paired.test.pValue, 0.0);
                assertEquals(at + ": df", single.test.degreesOfFreedom, paired.test.degreesOfFreedom, 0.0);
                assertEquals(at + ": estimate", single.estimate, paired.estimate, 0.0);
                assertEquals(at + ": standard error", single.standardError, paired.standardError, 0.0);
                assertEquals(at + ": lower", single.lower, paired.lower, 0.0);
                assertEquals(at + ": upper", single.upper, paired.upper, 0.0);
                assertEquals(at + ": name", "paired t", paired.test.test);
            }
        }
    }

    @Test
    public void testWelchAndPooledAgreeWhenTheSamplesAreIdentical() {
        // equal sizes and equal spreads is the case where the two formulas
        // reduce to each other, so only the rounding may differ
        double worstError = 0.0;
        double worstDf = 0.0;
        for (int n : new int[] { 3, 10, 100 }) {
            for (long seed = 1; seed <= 25; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double[] y = x.clone();
                TTestResult welch = HypothesisTests.tTwoSample(x, y, Alternative.TWO_SIDED, 0.95);
                TTestResult pooled = HypothesisTests.tTwoSamplePooled(x, y, Alternative.TWO_SIDED, 0.95);
                worstError = Math.max(worstError,
                        Math.abs(welch.standardError - pooled.standardError) / pooled.standardError);
                worstDf = Math.max(worstDf, Math.abs(welch.test.degreesOfFreedom - pooled.test.degreesOfFreedom)
                        / pooled.test.degreesOfFreedom);
                assertEquals("df should be 2n - 2", 2.0 * n - 2.0, pooled.test.degreesOfFreedom, 0.0);
            }
        }
        // measured: 2.7e-16 and 1.4e-16
        assertTrue("the standard errors differ by " + worstError, worstError < 1.0e-14);
        assertTrue("the degrees of freedom differ by " + worstDf, worstDf < 1.0e-14);
    }

    @Test
    public void testTheWelchDegreesOfFreedomLieBetweenTheirBounds() {
        // Welch-Satterthwaite never falls below the smaller sample and never
        // exceeds what the pooled test would claim
        for (int nx : new int[] { 3, 8, 40 }) {
            for (int ny : new int[] { 3, 8, 40 }) {
                for (double sy : new double[] { 0.1, 1.0, 30.0 }) {
                    for (long seed = 1; seed <= 25; seed++) {
                        double[] x = sample(new Normal(3.0, 2.0), nx, seed * 7919L + 1);
                        double[] y = sample(new Normal(3.0, sy), ny, seed * 104729L + 3);
                        double df = HypothesisTests.tTwoSample(x, y, Alternative.TWO_SIDED, 0.95)
                                .test.degreesOfFreedom;
                        String at = "nx=" + nx + " ny=" + ny + " sy=" + sy + " seed=" + seed + " df=" + df;
                        assertTrue(at, df >= Math.min(nx, ny) - 1.0 - 1.0e-9);
                        assertTrue(at, df <= nx + ny - 2.0 + 1.0e-9);
                    }
                }
            }
        }
    }

    @Test
    public void testTheTwoSampleStatisticIsScaleInvariant() {
        // scaling by a power of two is exact, and the whole computation is
        // built from quantities that scale alike, so not one bit may move --
        // even where the scaling sends the sum of squares out of range and
        // the deviation has to be taken the long way round
        double[] x = sample(new Normal(3.0, 2.0), 30, 5L);
        double[] y = sample(new Normal(2.5, 3.0), 25, 7L);
        double welch = HypothesisTests.tTwoSample(x, y, Alternative.TWO_SIDED, 0.95).test.statistic;
        double pooled = HypothesisTests.tTwoSamplePooled(x, y, Alternative.TWO_SIDED, 0.95).test.statistic;
        for (double factor : new double[] { 0x1.0p-900, 0x1.0p-40, 0x1.0p40, 0x1.0p900 }) {
            double[] xs = new double[x.length];
            double[] ys = new double[y.length];
            for (int i = 0; i < x.length; i++) {
                xs[i] = x[i] * factor;
            }
            for (int i = 0; i < y.length; i++) {
                ys[i] = y[i] * factor;
            }
            assertEquals("Welch at factor " + factor, welch,
                    HypothesisTests.tTwoSample(xs, ys, Alternative.TWO_SIDED, 0.95).test.statistic, 0.0);
            assertEquals("pooled at factor " + factor, pooled,
                    HypothesisTests.tTwoSamplePooled(xs, ys, Alternative.TWO_SIDED, 0.95).test.statistic, 0.0);
        }
    }

    @Test
    public void testTheTwoSampleStatisticIsUnchangedByAShift() {
        // the limit here is the input, not the arithmetic: at a shift of 1e12
        // data of spread two is quantized to 1.2e-4 by the double grid it
        // lands on, and four digits of the answer are gone before the test
        // sees the numbers. 1e6 still leaves ten
        double[] x = sample(new Normal(3.0, 2.0), 30, 5L);
        double[] y = sample(new Normal(2.5, 3.0), 25, 7L);
        double reference = HypothesisTests.tTwoSample(x, y, Alternative.TWO_SIDED, 0.95).test.statistic;
        double[] shifts = { 1.0, 1.0e6 };
        double[] tolerances = { 1.0e-11, 1.0e-6 };
        for (int k = 0; k < shifts.length; k++) {
            double[] xs = new double[x.length];
            double[] ys = new double[y.length];
            for (int i = 0; i < x.length; i++) {
                xs[i] = x[i] + shifts[k];
            }
            for (int i = 0; i < y.length; i++) {
                ys[i] = y[i] + shifts[k];
            }
            double got = HypothesisTests.tTwoSample(xs, ys, Alternative.TWO_SIDED, 0.95).test.statistic;
            assertEquals("shift " + shifts[k], reference, got, tolerances[k] * Math.abs(reference));
        }
    }

    @Test
    public void testTheEstimateIsTheFirstMeanMinusTheSecond() {
        // the sign convention the alternative is read against
        double[] x = { 10.0, 11.0, 12.0, 13.0 };
        double[] y = { 1.0, 2.0, 3.0, 4.0 };
        TTestResult r = HypothesisTests.tTwoSample(x, y, Alternative.GREATER, 0.95);
        assertEquals("estimate", 9.0, r.estimate, 1.0e-12);
        assertTrue("x is clearly the larger, so GREATER must not reject nothing", r.test.pValue < 0.01);
        TTestResult reversed = HypothesisTests.tTwoSample(y, x, Alternative.GREATER, 0.95);
        assertEquals("estimate reversed", -9.0, reversed.estimate, 1.0e-12);
        assertTrue("and the other way round it must not reject", reversed.test.pValue > 0.99);
    }

    @Test
    public void testTheTwoSampleGuardRail() {
        double[] ok = { 1.0, 2.0, 3.0 };
        double[] constant = { 2.0, 2.0, 2.0 };
        rejectsTwo("null x", null, ok);
        rejectsTwo("null y", ok, null);
        rejectsTwo("an empty x", new double[0], ok);
        rejectsTwo("a single observation in y", ok, new double[] { 1.0 });
        rejectsTwo("a NaN", new double[] { 1.0, Double.NaN, 3.0 }, ok);
        rejectsTwo("two constant samples", constant, constant);

        // one constant sample is not a problem: its variance is zero and
        // Welch simply gives all the degrees of freedom to the other one
        TTestResult oneFlat = HypothesisTests.tTwoSample(ok, constant, Alternative.TWO_SIDED, 0.95);
        assertEquals("df", 2.0, oneFlat.test.degreesOfFreedom, 1.0e-12);

        try {
            HypothesisTests.tPaired(ok, new double[] { 1.0, 2.0 }, Alternative.TWO_SIDED, 0.95);
            fail("tPaired accepted samples of different lengths");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("paired"));
        }
    }

    private static void rejectsTwo(String what, double[] x, double[] y) {
        try {
            HypothesisTests.tTwoSample(x, y, Alternative.TWO_SIDED, 0.95);
            fail("tTwoSample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            HypothesisTests.tTwoSamplePooled(x, y, Alternative.TWO_SIDED, 0.95);
            fail("tTwoSamplePooled accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------- the variance ratio --

    @Test
    public void testTheVarianceRatioPValuesAreUniformUnderTheNullHypothesis() {
        for (int nx : new int[] { 5, 15, 60 }) {
            int ny = Math.max(3, nx / 2);
            int reps = 5000;
            int below = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] a = sample(new Normal(3.0, 2.0), nx, seed * 7919L + 1);
                double[] b = sample(new Normal(-1.0, 2.0), ny, seed * 104729L + 3);
                if (HypothesisTests.fVarianceRatio(a, b, Alternative.TWO_SIDED).test.rejectsAt(0.05)) {
                    below++;
                }
            }
            // measured over 12000 replications: 0.0461 .. 0.0505
            assertRate("nx=" + nx, below / (double) reps, 0.05, 0.015);
        }
    }

    @Test
    public void testSwappingTheSamplesInvertsTheVarianceRatio() {
        // F(x, y) = 1 / F(y, x) with the degrees of freedom swapped, so the
        // two-sided answer may not depend on which sample came first and the
        // one-sided ones have to change places. This is what pins the
        // two-sided p-value as twice the smaller tail rather than twice the
        // lower one, which is the classic way to get it wrong
        double worstTwoSided = 0.0;
        double worstOneSided = 0.0;
        double worstProduct = 0.0;
        for (int nx : new int[] { 3, 8, 40 }) {
            for (int ny : new int[] { 3, 8, 40 }) {
                for (double spread : new double[] { 0.2, 1.0, 5.0 }) {
                    for (long seed = 1; seed <= 25; seed++) {
                        double[] a = sample(new Normal(3.0, 2.0), nx, seed * 7919L + 1);
                        double[] b = sample(new Normal(3.0, 2.0 * spread), ny, seed * 104729L + 3);
                        FTestResult ab = HypothesisTests.fVarianceRatio(a, b, Alternative.TWO_SIDED);
                        FTestResult ba = HypothesisTests.fVarianceRatio(b, a, Alternative.TWO_SIDED);
                        assertEquals("the degrees of freedom did not swap", ab.numeratorDf, ba.denominatorDf, 0.0);
                        assertEquals("the degrees of freedom did not swap", ab.denominatorDf, ba.numeratorDf, 0.0);
                        worstTwoSided = Math.max(worstTwoSided,
                                Math.abs(ab.test.pValue - ba.test.pValue) / ab.test.pValue);
                        worstProduct = Math.max(worstProduct,
                                Math.abs(ab.test.statistic * ba.test.statistic - 1.0));
                        double greater = HypothesisTests.fVarianceRatio(a, b, Alternative.GREATER).test.pValue;
                        double less = HypothesisTests.fVarianceRatio(b, a, Alternative.LESS).test.pValue;
                        worstOneSided = Math.max(worstOneSided, Math.abs(greater - less) / greater);
                    }
                }
            }
        }
        // measured: 1.1e-14, 1.0e-14 and 4.4e-16
        assertTrue("the two-sided p-value moved by " + worstTwoSided, worstTwoSided < 1.0e-12);
        assertTrue("GREATER and LESS disagree by " + worstOneSided, worstOneSided < 1.0e-12);
        assertTrue("F(x,y) * F(y,x) is off by " + worstProduct, worstProduct < 1.0e-13);
    }

    @Test
    public void testTheVarianceRatioTailsAddUp() {
        double worst = 0.0;
        for (int n : new int[] { 5, 30 }) {
            for (long seed = 1; seed <= 100; seed++) {
                double[] a = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double[] b = sample(new Normal(3.0, 2.0), n, seed * 104729L + 3);
                double less = HypothesisTests.fVarianceRatio(a, b, Alternative.LESS).test.pValue;
                double greater = HypothesisTests.fVarianceRatio(a, b, Alternative.GREATER).test.pValue;
                worst = Math.max(worst, Math.abs(less + greater - 1.0));
            }
        }
        // measured: 1.6e-15
        assertTrue("the tails do not partition, off by " + worst, worst < 1.0e-13);
    }

    @Test
    public void testTheUpperVarianceRatioTailSurvivesFarOut() {
        // the upper tail is taken through the reciprocal distribution rather
        // than as 1 - cdf, which would have returned a flat zero here
        double[] wide = sample(new Normal(0.0, 100.0), 40, 5L);
        double[] narrow = sample(new Normal(0.0, 1.0), 40, 7L);
        TestResult r = HypothesisTests.fVarianceRatio(wide, narrow, Alternative.GREATER).test;
        assertTrue("the ratio should be enormous : " + r.statistic, r.statistic > 1.0e3);
        // measured: 1.6e-71
        assertTrue("the p-value collapsed to " + r.pValue, r.pValue > 0.0 && r.pValue < 1.0e-40);
    }

    @Test
    public void testTheVarianceRatioIsScaleInvariant() {
        double[] a = sample(new Normal(3.0, 2.0), 20, 5L);
        double[] b = sample(new Normal(1.0, 3.0), 15, 7L);
        double reference = HypothesisTests.fVarianceRatio(a, b, Alternative.TWO_SIDED).test.statistic;
        for (double factor : new double[] { 0x1.0p-900, 0x1.0p900 }) {
            double[] as = new double[a.length];
            double[] bs = new double[b.length];
            for (int i = 0; i < a.length; i++) {
                as[i] = a[i] * factor;
            }
            for (int i = 0; i < b.length; i++) {
                bs[i] = b[i] * factor;
            }
            assertEquals("factor " + factor, reference,
                    HypothesisTests.fVarianceRatio(as, bs, Alternative.TWO_SIDED).test.statistic, 0.0);
        }
    }

    @Test
    public void testThePooledTStatisticSquaredIsTheVarianceRatioOfTheAnova() {
        // t^2 = F is the identity that ties the two-sample t to the F family,
        // and it would break on a wrong pooled variance or a wrong denominator
        double worst = 0.0;
        for (int n : new int[] { 3, 10, 60 }) {
            for (long seed = 1; seed <= 40; seed++) {
                double[] p = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                double[] q = sample(new Normal(2.0, 2.0), n, seed * 104729L + 3);
                double t = HypothesisTests.tTwoSamplePooled(p, q, Alternative.TWO_SIDED, 0.95).test.statistic;
                worst = Math.max(worst, Math.abs(t * t - anovaF(p, q)) / anovaF(p, q));
            }
        }
        // measured: 1.3e-15
        assertTrue("t^2 and F differ by " + worst, worst < 1.0e-12);
    }

    /** The one-way ANOVA F for two groups, from its definition. */
    private static double anovaF(double[] p, double[] q) {
        int n = p.length + q.length;
        double mp = mean(p);
        double mq = mean(q);
        double grand = (mp * p.length + mq * q.length) / n;
        double between = p.length * (mp - grand) * (mp - grand) + q.length * (mq - grand) * (mq - grand);
        double within = 0.0;
        for (int i = 0; i < p.length; i++) {
            within += (p[i] - mp) * (p[i] - mp);
        }
        for (int i = 0; i < q.length; i++) {
            within += (q[i] - mq) * (q[i] - mq);
        }
        return between / (within / (n - 2.0));
    }

    private static double mean(double[] x) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        return sum / x.length;
    }

    // --------------------------------------------------------- chi-squared --

    @Test
    public void testTheChiSquaredPValuesAreUniformUnderTheNullHypothesis() {
        for (int categories : new int[] { 2, 5, 20 }) {
            for (int total : new int[] { 60, 400 }) {
                int reps = 4000;
                int below = 0;
                double[] expected = new double[categories];
                for (int i = 0; i < categories; i++) {
                    expected[i] = total / (double) categories;
                }
                for (long seed = 1; seed <= reps; seed++) {
                    long[] observed = multinomial(categories, total, seed * 7919L + 1);
                    if (HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 0).rejectsAt(0.05)) {
                        below++;
                    }
                }
                // measured over 8000 replications: 0.0506 .. 0.0578, the high
                // end at two categories where the statistic is coarsest
                assertRate("k=" + categories + " n=" + total, below / (double) reps, 0.05, 0.020);
            }
        }
    }

    @Test
    public void testTheIndependencePValuesAreUniformUnderTheNullHypothesis() {
        for (int[] shape : new int[][] { { 2, 2 }, { 3, 4 }, { 5, 5 } }) {
            int reps = 4000;
            int below = 0;
            for (long seed = 1; seed <= reps; seed++) {
                long[][] table = independentTable(shape[0], shape[1], 400, seed * 7919L + 1);
                if (HypothesisTests.chiSquaredIndependence(table).rejectsAt(0.05)) {
                    below++;
                }
            }
            // measured over 8000 replications: 0.0476 .. 0.0520
            assertRate(shape[0] + "x" + shape[1], below / (double) reps, 0.05, 0.015);
        }
    }

    @Test
    public void testTwoByTwoIndependenceAgreesWithFishersExactTest() {
        // an entirely different algorithm for the same question: Fisher sums
        // the exact null distribution of the table from Hypergeometric, where
        // the chi-squared test approximates it. They must reach the same
        // verdict on all but the borderline tables
        int disagree = 0;
        int total = 0;
        double worstGap = 0.0;
        for (long seed = 1; seed <= 1500; seed++) {
            long[][] table = independentTable(2, 2, 200, seed * 7919L + 1);
            double approximate = HypothesisTests.chiSquaredIndependence(table).pValue;
            double exact = fisherExact(table);
            total++;
            if ((approximate < 0.05) != (exact < 0.05)) {
                disagree++;
            }
            worstGap = Math.max(worstGap, Math.abs(approximate - exact));
        }
        // measured over 4000 tables: 24 disagreements and a worst gap of 0.114
        assertTrue("the two disagree on " + disagree + " of " + total, disagree < total / 20);
        assertTrue("the p-values are " + worstGap + " apart", worstGap < 0.20);

        // and on a table that is plainly dependent, both must be decisive
        long[][] dependent = { { 90L, 10L }, { 10L, 90L } };
        assertTrue("chi-squared", HypothesisTests.chiSquaredIndependence(dependent).pValue < 1.0e-20);
        assertTrue("Fisher", fisherExact(dependent) < 1.0e-20);
    }

    /** The two-sided Fisher exact p-value of a 2x2 table, summed from the null. */
    private static double fisherExact(long[][] t) {
        int a = (int) t[0][0];
        int b = (int) t[0][1];
        int c = (int) t[1][0];
        int d = (int) t[1][1];
        int firstRow = a + b;
        int firstColumn = a + c;
        int n = a + b + c + d;
        Hypergeometric nullDistribution = new Hypergeometric(n, firstColumn, firstRow);
        double observed = nullDistribution.pmf(a);
        double p = 0.0;
        int lo = Math.max(0, firstRow + firstColumn - n);
        int hi = Math.min(firstRow, firstColumn);
        for (int k = lo; k <= hi; k++) {
            double pk = nullDistribution.pmf(k);
            if (pk <= observed * (1.0 + 1.0e-9)) {
                p += pk;
            }
        }
        return Math.min(1.0, p);
    }

    @Test
    public void testTheChiSquaredStatisticIsTheDefinition() {
        long[] observed = { 18L, 22L, 30L, 30L };
        double[] expected = { 25.0, 25.0, 25.0, 25.0 };
        double byHand = 0.0;
        for (int i = 0; i < observed.length; i++) {
            double residual = observed[i] - expected[i];
            byHand += residual * residual / expected[i];
        }
        TestResult r = HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 0);
        assertEquals("statistic", byHand, r.statistic, 0.0);
        assertEquals("degrees of freedom", 3.0, r.degreesOfFreedom, 0.0);
        assertEquals("the tail is the upper one", Alternative.GREATER, r.alternative);

        // a perfect fit has a statistic of zero and a p-value of one
        TestResult perfect = HypothesisTests.chiSquaredGoodnessOfFit(new long[] { 25L, 25L, 25L, 25L }, expected,
                0);
        assertEquals("statistic", 0.0, perfect.statistic, 0.0);
        assertEquals("p-value", 1.0, perfect.pValue, 0.0);
    }

    @Test
    public void testEstimatedParametersCostDegreesOfFreedom() {
        // fitting the expected counts to the same sample uses up degrees of
        // freedom, and a test that does not subtract them is too forgiving
        long[] observed = { 18L, 22L, 30L, 30L };
        double[] expected = { 25.0, 25.0, 25.0, 25.0 };
        TestResult none = HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 0);
        TestResult one = HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 1);
        assertEquals("df without", 3.0, none.degreesOfFreedom, 0.0);
        assertEquals("df with one", 2.0, one.degreesOfFreedom, 0.0);
        assertEquals("the statistic must not change", none.statistic, one.statistic, 0.0);
        assertTrue("subtracting a degree of freedom has to lower the p-value", one.pValue < none.pValue);

        try {
            HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 3);
            fail("three estimated parameters from four categories was accepted");
        } catch (IllegalArgumentException expectedFailure) {
            assertTrue(expectedFailure.getMessage(), expectedFailure.getMessage().contains("degrees of freedom"));
        }
    }

    @Test
    public void testTheIndependenceStatisticIsTheDefinition() {
        long[][] table = { { 30L, 20L }, { 15L, 35L } };
        double total = 100.0;
        double[] rows = { 50.0, 50.0 };
        double[] columns = { 45.0, 55.0 };
        double byHand = 0.0;
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                double expected = rows[i] * columns[k] / total;
                double residual = table[i][k] - expected;
                byHand += residual * residual / expected;
            }
        }
        TestResult r = HypothesisTests.chiSquaredIndependence(table);
        assertEquals("statistic", byHand, r.statistic, 1.0e-12);
        assertEquals("degrees of freedom", 1.0, r.degreesOfFreedom, 0.0);

        // transposing a table cannot change whether its two directions are
        // independent
        long[][] transposed = { { 30L, 15L }, { 20L, 35L } };
        assertEquals("transposing changed the statistic", r.statistic,
                HypothesisTests.chiSquaredIndependence(transposed).statistic, 1.0e-12);
    }

    // --------------------------------------------------- the guard rail, 3 --

    @Test
    public void testTheChiSquaredGuardRail() {
        long[] observed = { 10L, 20L, 30L };
        double[] expected = { 20.0, 20.0, 20.0 };
        rejectsChiSquared("null observed", null, expected, 0);
        rejectsChiSquared("null expected", observed, null, 0);
        rejectsChiSquared("mismatched lengths", observed, new double[] { 30.0, 30.0 }, 0);
        rejectsChiSquared("one category", new long[] { 10L }, new double[] { 10.0 }, 0);
        rejectsChiSquared("a negative count", new long[] { -1L, 20L, 41L }, expected, 0);
        rejectsChiSquared("a zero expected count", observed, new double[] { 0.0, 30.0, 30.0 }, 0);
        rejectsChiSquared("probabilities instead of counts", observed,
                new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 }, 0);
        rejectsChiSquared("a negative parameter count", observed, expected, -1);

        rejectsIndependence("null", null);
        rejectsIndependence("one row", new long[][] { { 1L, 2L } });
        rejectsIndependence("one column", new long[][] { { 1L }, { 2L } });
        rejectsIndependence("a null row", new long[][] { { 1L, 2L }, null });
        rejectsIndependence("a ragged table", new long[][] { { 1L, 2L }, { 3L } });
        rejectsIndependence("a negative count", new long[][] { { 1L, 2L }, { -3L, 4L } });
        rejectsIndependence("an empty row", new long[][] { { 1L, 2L }, { 0L, 0L } });
        rejectsIndependence("an empty column", new long[][] { { 0L, 2L }, { 0L, 4L } });
    }

    private static void rejectsChiSquared(String what, long[] observed, double[] expected, int parameters) {
        try {
            HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, parameters);
            fail("chiSquaredGoodnessOfFit accepted " + what);
        } catch (IllegalArgumentException expectedFailure) {
            assertTrue("threw without saying what was wrong", expectedFailure.getMessage() != null);
        }
    }

    private static void rejectsIndependence(String what, long[][] table) {
        try {
            HypothesisTests.chiSquaredIndependence(table);
            fail("chiSquaredIndependence accepted " + what);
        } catch (IllegalArgumentException expectedFailure) {
            assertTrue("threw without saying what was wrong", expectedFailure.getMessage() != null);
        }
    }

    @Test
    public void testTheVarianceRatioGuardRail() {
        double[] ok = { 1.0, 2.0, 3.0 };
        double[] constant = { 2.0, 2.0, 2.0 };
        rejectsF("null x", null, ok);
        rejectsF("null y", ok, null);
        rejectsF("a single observation", new double[] { 1.0 }, ok);
        rejectsF("a NaN", new double[] { 1.0, Double.NaN, 3.0 }, ok);
        rejectsF("a constant denominator", ok, constant);

        // a constant numerator is fine: the ratio is zero, and zero variance
        // is as far from the null as the lower tail reaches
        FTestResult flat = HypothesisTests.fVarianceRatio(constant, ok, Alternative.LESS);
        assertEquals("statistic", 0.0, flat.test.statistic, 0.0);
        assertEquals("p-value", 0.0, flat.test.pValue, 0.0);
    }

    private static void rejectsF(String what, double[] x, double[] y) {
        try {
            HypothesisTests.fVarianceRatio(x, y, Alternative.TWO_SIDED);
            fail("fVarianceRatio accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------- the fixtures --

    /** Counts from a uniform multinomial, so the null hypothesis really holds. */
    private static long[] multinomial(int categories, int total, long seed) {
        Lcg lcg = new Lcg(seed);
        long[] counts = new long[categories];
        for (int i = 0; i < total; i++) {
            counts[(int) (lcg.next() * categories)]++;
        }
        return counts;
    }

    /** A table whose row and column really are independent. */
    private static long[][] independentTable(int rows, int columns, int total, long seed) {
        Lcg lcg = new Lcg(seed);
        long[][] table = new long[rows][columns];
        for (int i = 0; i < total; i++) {
            table[(int) (lcg.next() * rows)][(int) (lcg.next() * columns)]++;
        }
        // an empty margin is not a table, so nudge a count into one
        for (int i = 0; i < rows; i++) {
            long row = 0L;
            for (int k = 0; k < columns; k++) {
                row += table[i][k];
            }
            if (row == 0L) {
                table[i][0] = 1L;
            }
        }
        for (int k = 0; k < columns; k++) {
            long column = 0L;
            for (int i = 0; i < rows; i++) {
                column += table[i][k];
            }
            if (column == 0L) {
                table[0][k] = 1L;
            }
        }
        return table;
    }

    // -------------------------------------------------- Kolmogorov-Smirnov --

    @Test
    public void testTheKolmogorovSmirnovPValuesAreUniformUnderTheNullHypothesis() {
        Uniform uniform = new Uniform(0.0, 1.0);
        for (int n : new int[] { 1, 5, 30, 200 }) {
            int reps = 4000;
            for (Alternative alternative : Alternative.values()) {
                int below = 0;
                for (long seed = 1; seed <= reps; seed++) {
                    double[] x = sample(uniform, n, seed * 7919L + 1);
                    if (HypothesisTests.kolmogorovSmirnov(x, uniform, alternative).rejectsAt(0.05)) {
                        below++;
                    }
                }
                // measured over 12000 replications: 0.0487 .. 0.0570
                assertRate("n=" + n + " " + alternative, below / (double) reps, 0.05, 0.018);
            }
        }
    }

    @Test
    public void testTheKolmogorovSmirnovStatisticIsTheDefinition() {
        for (int n : new int[] { 1, 7, 60 }) {
            for (long seed = 1; seed <= 60; seed++) {
                ContinuousDistribution d = new Normal(3.0, 2.0);
                double[] x = sample(d, n, seed * 7919L + 1);
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = d.cdf(x[i]);
                }
                Arrays.sort(u);
                double dPlus = 0.0;
                double dMinus = 0.0;
                for (int i = 0; i < n; i++) {
                    dPlus = Math.max(dPlus, (i + 1.0) / n - u[i]);
                    dMinus = Math.max(dMinus, u[i] - i / (double) n);
                }
                String at = "n=" + n + " seed=" + seed;
                assertEquals(at + ": LESS is D+", dPlus,
                        HypothesisTests.kolmogorovSmirnov(x, d, Alternative.LESS).statistic, 0.0);
                assertEquals(at + ": GREATER is D-", dMinus,
                        HypothesisTests.kolmogorovSmirnov(x, d, Alternative.GREATER).statistic, 0.0);
                assertEquals(at + ": two-sided is the larger of them", Math.max(dPlus, dMinus),
                        HypothesisTests.kolmogorovSmirnov(x, d, Alternative.TWO_SIDED).statistic, 0.0);
            }
        }
    }

    @Test
    public void testTheKolmogorovSmirnovAgreesWithThePublishedCriticalValues() {
        // the two-sided critical values of D_n, so barF at each has to give
        // back the level it was tabulated for
        double[][] critical = { { 5, 0.56328, 0.66853 }, { 10, 0.40925, 0.48893 }, { 20, 0.29408, 0.35241 },
                { 50, 0.18841, 0.22585 }, { 100, 0.13403, 0.16073 } };
        double worst = 0.0;
        for (int i = 0; i < critical.length; i++) {
            int n = (int) critical[i][0];
            worst = Math.max(worst, Math.abs(KolmogorovSmirnov.barF(n, critical[i][1]) - 0.05));
            worst = Math.max(worst, Math.abs(KolmogorovSmirnov.barF(n, critical[i][2]) - 0.01));
        }
        // measured: 8.8e-5, and that is the table rounding, not the algorithm
        assertTrue("the worst deviation from the table is " + worst, worst < 2.0e-4);
    }

    @Test
    public void testTheOneSidedAlternativesPointTheWayTheyClaim() {
        // GREATER asks whether the sample is stochastically larger, which is
        // its distribution function lying below the hypothesized one. Getting
        // this backwards is a silent error, so it is pinned here
        ContinuousDistribution standard = new Normal(0.0, 1.0);
        double[] high = sample(new Normal(0.6, 1.0), 200, 11L);
        assertTrue("a sample shifted up must look larger",
                HypothesisTests.kolmogorovSmirnov(high, standard, Alternative.GREATER).pValue < 1.0e-6);
        assertTrue("and must not look smaller",
                HypothesisTests.kolmogorovSmirnov(high, standard, Alternative.LESS).pValue > 0.9);

        double[] low = sample(new Normal(-0.6, 1.0), 200, 11L);
        assertTrue("a sample shifted down must look smaller",
                HypothesisTests.kolmogorovSmirnov(low, standard, Alternative.LESS).pValue < 1.0e-6);
        assertTrue("and must not look larger",
                HypothesisTests.kolmogorovSmirnov(low, standard, Alternative.GREATER).pValue > 0.9);
    }

    @Test
    public void testASingleObservationGetsTheExactPValue() {
        // the superseded facade in math.stats.gof stored the statistic where
        // the p-value belongs for a sample of one, which inverted the
        // decision. For n = 1 the null distributions are closed forms:
        // P[D_1 >= d] = 2 - 2d and P[D_1^+ >= d] = 1 - d
        Uniform uniform = new Uniform(0.0, 1.0);
        for (double u : new double[] { 0.01, 0.25, 0.5, 0.75, 0.99 }) {
            double[] one = { u };
            TestResult two = HypothesisTests.kolmogorovSmirnov(one, uniform, Alternative.TWO_SIDED);
            TestResult less = HypothesisTests.kolmogorovSmirnov(one, uniform, Alternative.LESS);
            TestResult greater = HypothesisTests.kolmogorovSmirnov(one, uniform, Alternative.GREATER);
            String at = "u=" + u;
            assertEquals(at + ": D+", 1.0 - u, less.statistic, 1.0e-15);
            assertEquals(at + ": D-", u, greater.statistic, 1.0e-15);
            assertEquals(at + ": two-sided p", 2.0 - 2.0 * Math.max(u, 1.0 - u), two.pValue, 1.0e-15);
            assertEquals(at + ": p(LESS)", u, less.pValue, 1.0e-15);
            assertEquals(at + ": p(GREATER)", 1.0 - u, greater.pValue, 1.0e-15);
        }
    }

    @Test
    public void testTheKolmogorovSmirnovTestIsInvariantUnderAMonotoneTransform() {
        // the transform by the distribution function is what the statistic is
        // built on, so testing a shifted and scaled sample against the shifted
        // and scaled distribution has to give the same number
        double worst = 0.0;
        for (long seed = 1; seed <= 40; seed++) {
            double[] raw = sample(new Normal(0.0, 1.0), 50, seed * 7919L + 1);
            double[] moved = new double[raw.length];
            for (int i = 0; i < raw.length; i++) {
                moved[i] = 4.0 + 3.0 * raw[i];
            }
            double before = HypothesisTests.kolmogorovSmirnov(raw, new Normal(0.0, 1.0), Alternative.TWO_SIDED)
                    .statistic;
            double after = HypothesisTests.kolmogorovSmirnov(moved, new Normal(4.0, 3.0), Alternative.TWO_SIDED)
                    .statistic;
            worst = Math.max(worst, Math.abs(before - after));
        }
        // measured: 1.1e-16
        assertTrue("the statistic moved by " + worst, worst < 1.0e-14);
    }

    @Test
    public void testTheKolmogorovSmirnovTestRejectsTheWrongDistribution() {
        double[] x = sample(new Normal(0.0, 1.0), 100, 13L);
        TestResult right = HypothesisTests.kolmogorovSmirnov(x, new Normal(0.0, 1.0), Alternative.TWO_SIDED);
        TestResult wrong = HypothesisTests.kolmogorovSmirnov(x, new Exponential(1.0), Alternative.TWO_SIDED);
        assertTrue("the truth was rejected : " + right.pValue, right.pValue > 0.05);
        assertTrue("an exponential passed for a normal : " + wrong.pValue, wrong.pValue < 1.0e-15);
    }

    @Test
    public void testTheKolmogorovSmirnovGuardRail() {
        Uniform uniform = new Uniform(0.0, 1.0);
        double[] ok = { 0.2, 0.4, 0.6 };
        rejectsKs("null", null, uniform, Alternative.TWO_SIDED);
        rejectsKs("an empty sample", new double[0], uniform, Alternative.TWO_SIDED);
        rejectsKs("a NaN", new double[] { 0.2, Double.NaN }, uniform, Alternative.TWO_SIDED);
        rejectsKs("an infinity", new double[] { 0.2, Double.POSITIVE_INFINITY }, uniform, Alternative.TWO_SIDED);
        rejectsKs("a null distribution", ok, null, Alternative.TWO_SIDED);
        rejectsKs("a null alternative", ok, uniform, null);
    }

    private static void rejectsKs(String what, double[] x, ContinuousDistribution d, Alternative alternative) {
        try {
            HypothesisTests.kolmogorovSmirnov(x, d, alternative);
            fail("kolmogorovSmirnov accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ----------------------------------------------------- Anderson-Darling --

    @Test
    public void testTheAndersonDarlingStatisticIsTheDefinition() {
        Normal standard = new Normal(0.0, 1.0);
        int saturated = 0;
        int compared = 0;
        for (int n : new int[] { 1, 2, 3, 9, 60, 400 }) {
            for (long seed = 1; seed <= 40; seed++) {
                // a wide sample read through a narrow distribution, so that
                // some of the transformed values land on 0 or 1 exactly
                double[] x = sample(new Normal(2.0, 3.0), n, seed * 7919L + 1);
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = standard.cdf(x[i]);
                }
                Arrays.sort(u);
                double sum = 0.0;
                for (int i = 0; i < n; i++) {
                    sum += (2 * i + 1) * (Math.log(u[i]) + Math.log(1.0 - u[n - 1 - i]));
                }
                double definition = -n - sum / n;
                String at = "n=" + n + " seed=" + seed;
                double got = HypothesisTests.andersonDarling(x, standard).statistic;

                if (Double.isInfinite(definition)) {
                    // the definition itself has no value here; what the
                    // method must not do is pass the infinity on
                    saturated++;
                    assertTrue(at + ": a saturated sample gave " + got,
                            got > 0.0 && !Double.isInfinite(got) && !Double.isNaN(got));
                } else {
                    compared++;
                    assertEquals(at, definition, got, 1.0e-12);
                }
            }
        }
        // measured: 163 of the 240 samples have a value the definition can
        // state, and 77 do not
        assertTrue("nothing was compared against the definition", compared > 150);
        assertTrue("no sample saturated the distribution function", saturated > 0);
    }

    @Test
    public void testTheAndersonDarlingNullAgreesWithAnInversionOfItsOwnLimit() {
        // A_n^2 has the same shape as W_n^2 in the limit, a weighted sum of
        // chi-squares, and its weights are 1 / (j (j + 1)). That gives a second
        // route to the same distribution which shares nothing with the first:
        // one is an interpolation table with an empirical correction in 1/n,
        // the other an inversion of a characteristic function
        int terms = 200;
        double[] lambda = new double[terms];
        double kept = 0.0;
        for (int j = 1; j <= terms; j++) {
            lambda[j - 1] = 1.0 / (j * (j + 1.0));
            kept += lambda[j - 1];
        }
        // the whole series sums to one, and what is left out is subtracted
        // from the argument rather than ignored
        double dropped = 1.0 - kept;

        double worst = 0.0;
        double[] points = { 1.933, 2.492, 3.070, 3.857 };
        for (int i = 0; i < points.length; i++) {
            double viaTable = AndersonDarling.barF(100000, points[i]);
            double viaInversion = WeightedChiSquare.barF(lambda, points[i] - dropped);
            worst = Math.max(worst, Math.abs(viaTable - viaInversion));
        }
        // measured: 2.1e-6, which is well inside the three decimal digits the
        // table claims for itself
        assertTrue("the two routes to the limit differ by " + worst, worst < 1.0e-5);
    }

    @Test
    public void testTheAndersonDarlingNullDistributionAgreesWithASimulatedOne() {
        // the tabulated null distribution is an interpolation with an
        // empirical correction in 1/n, so it is checked against the only
        // thing that is not an approximation: the statistic itself, drawn
        Uniform uniform = new Uniform(0.0, 1.0);
        int reps = 100000;
        double worst = 0.0;
        for (int n : new int[] { 2, 5, 20, 100 }) {
            double[] drawn = new double[reps];
            for (int b = 0; b < reps; b++) {
                drawn[b] = HypothesisTests.andersonDarling(sample(uniform, n, (b + 1L) * 7919L + 1),
                        uniform).statistic;
            }
            Arrays.sort(drawn);
            for (double tail : new double[] { 0.50, 0.25, 0.10, 0.05, 0.01 }) {
                double quantile = drawn[(int) ((1.0 - tail) * reps)];
                worst = Math.max(worst, Math.abs(AndersonDarling.barF(n, quantile) - tail));
            }
        }
        // measured over three independent runs: 0.0045, and the largest
        // deviations are at n = 2 where the correction has least to work with
        assertTrue("the tabulated tail is off by " + worst, worst < 0.008);
    }

    @Test
    public void testTheAndersonDarlingRejectsAtItsLevelUnderTheNullHypothesis() {
        Uniform uniform = new Uniform(0.0, 1.0);
        for (int n : new int[] { 1, 5, 30, 200 }) {
            int reps = 4000;
            int atFive = 0;
            int atOne = 0;
            for (long seed = 1; seed <= reps; seed++) {
                TestResult result = HypothesisTests.andersonDarling(sample(uniform, n, seed * 7919L + 1),
                        uniform);
                if (result.rejectsAt(0.05)) {
                    atFive++;
                }
                if (result.rejectsAt(0.01)) {
                    atOne++;
                }
            }
            // measured: 0.0460 .. 0.0525 and 0.0095 .. 0.0105
            assertRate("n=" + n + " at 0.05", atFive / (double) reps, 0.05, 0.018);
            assertRate("n=" + n + " at 0.01", atOne / (double) reps, 0.01, 0.008);
        }
    }

    @Test
    public void testTheAndersonDarlingSeesATailDepartureTheKolmogorovSmirnovMisses() {
        // this is the reason the test exists. The contamination lengthens
        // five percent of the observations and leaves the rest alone, so it
        // shows up where the empirical distribution function has almost no
        // room to move and the supremum of an unweighted difference is blind
        Exponential unit = new Exponential(1.0);
        for (int n : new int[] { 50, 200 }) {
            int reps = 2000;
            int andersonDarling = 0;
            int kolmogorovSmirnov = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = contaminatedExponential(n, seed * 7919L + 1);
                if (HypothesisTests.andersonDarling(x, unit).rejectsAt(0.05)) {
                    andersonDarling++;
                }
                if (HypothesisTests.kolmogorovSmirnov(x, unit, Alternative.TWO_SIDED).rejectsAt(0.05)) {
                    kolmogorovSmirnov++;
                }
            }
            // measured: 0.1270 against 0.0670 at n = 50, 0.2415 against
            // 0.1080 at n = 200
            assertTrue("n=" + n + ": AD rejected " + andersonDarling + " of " + reps + ", KS "
                    + kolmogorovSmirnov, andersonDarling > kolmogorovSmirnov * 1.5);
        }
    }

    /** A unit exponential with five percent of a five times longer one mixed in. */
    private static double[] contaminatedExponential(int n, long seed) {
        Lcg lcg = new Lcg(seed);
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            double u = lcg.next();
            double scale = lcg.next() < 0.05 ? 5.0 : 1.0;
            x[i] = -scale * Math.log1p(-u);
        }
        return x;
    }

    @Test
    public void testTheAndersonDarlingStatisticDoesNotDependOnTheOrderOfTheSample() {
        Normal standard = new Normal(0.0, 1.0);
        for (long seed = 1; seed <= 40; seed++) {
            double[] x = sample(new Normal(1.0, 2.0), 40, seed * 7919L + 1);
            double reference = HypothesisTests.andersonDarling(x, standard).statistic;

            double[] reversed = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                reversed[i] = x[x.length - 1 - i];
            }
            double[] sorted = x.clone();
            Arrays.sort(sorted);

            // the statistic reads the sample only through its order
            // statistics, so these do identical arithmetic
            assertEquals("reversed, seed=" + seed, reference,
                    HypothesisTests.andersonDarling(reversed, standard).statistic, 0.0);
            assertEquals("sorted, seed=" + seed, reference,
                    HypothesisTests.andersonDarling(sorted, standard).statistic, 0.0);
        }
    }

    @Test
    public void testTheAndersonDarlingSurvivesADistributionFunctionThatSaturates() {
        // a double-valued distribution function reaches 0 and 1 long before
        // the mathematics does, and the logarithm of either is infinite. An
        // infinite statistic would be a p-value of exactly zero, which says
        // more than the arithmetic knows
        Normal standard = new Normal(0.0, 1.0);
        double[] extreme = { -50.0, -40.0, 0.0, 40.0, 50.0 };
        TestResult result = HypothesisTests.andersonDarling(extreme, standard);
        assertTrue("the statistic is " + result.statistic, result.statistic > 0.0
                && !Double.isInfinite(result.statistic) && !Double.isNaN(result.statistic));
        // measured: 55.1652 and 6.4e-25
        assertTrue("the p-value is " + result.pValue, result.pValue > 0.0 && result.pValue < 1.0e-20);

        // and the sample really does saturate the distribution function
        assertEquals("the far left of a standard normal", 0.0, standard.cdf(-50.0), 0.0);
        assertEquals("the far right of a standard normal", 1.0, standard.cdf(50.0), 0.0);
    }

    @Test
    public void testTheAndersonDarlingTailFallsExceptAtOneSeam() {
        // the tabulated null distribution changes method at A2 = 5, and the
        // two pieces do not meet: the tail steps back up there. It is pinned
        // rather than fixed because the file is a port that is kept as it was,
        // and the step is far smaller than the accuracy the port claims
        for (int n : new int[] { 2, 6, 20, 100 }) {
            double previous = Double.MAX_VALUE;
            double worstRise = 0.0;
            double at = Double.NaN;
            for (int step = 1; step < 30000; step++) {
                double a2 = step * 0.001;
                double tail = AndersonDarling.barF(n, a2);
                if (tail > previous && tail - previous > worstRise) {
                    worstRise = tail - previous;
                    at = a2;
                }
                previous = tail;
            }
            // measured: 1.7e-4 at n = 2 falling to 1.1e-8 at n = 100, always
            // at A2 = 5 exactly
            assertTrue("n=" + n + ": the tail rose by " + worstRise + " near A2=" + at,
                    worstRise < 2.0e-4);
            assertTrue("n=" + n + ": the rise is not at the seam but at " + at,
                    Double.isNaN(at) || Math.abs(at - 5.0) < 0.002);
        }
    }

    @Test
    public void testAndersonDarlingRejectsWhatItCannotTest() {
        Normal standard = new Normal(0.0, 1.0);
        rejectsAndersonDarling("a null sample", null, standard);
        rejectsAndersonDarling("an empty sample", new double[0], standard);
        rejectsAndersonDarling("a NaN observation", new double[] { 1.0, Double.NaN }, standard);
        rejectsAndersonDarling("an infinite observation", new double[] { 1.0, Double.POSITIVE_INFINITY },
                standard);
        rejectsAndersonDarling("a null distribution", new double[] { 1.0, 2.0 }, null);

        // and the null distribution on its own
        try {
            AndersonDarling.barF(0, 1.0);
            fail("barF accepted a sample size of zero");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsAndersonDarling(String what, double[] x, ContinuousDistribution d) {
        try {
            HypothesisTests.andersonDarling(x, d);
            fail("andersonDarling accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // --------------------------------------------------------- Fisher exact --

    @Test
    public void testTheExactTestKeepsItsLevelWhereTheApproximationDoesNot() {
        // the whole reason the exact test exists. A test that holds its level
        // never rejects more often than the level allows when the null
        // hypothesis is true; the chi-squared approximation does, and the
        // sparser the table the worse it gets
        int[] totals = { 12, 20, 30, 50, 100, 400 };
        double[] levels = { 0.01, 0.05, 0.10 };
        int fisherAbove = 0;
        int chiSquaredAbove = 0;
        double leastPower = 1.0;
        for (int i = 0; i < totals.length; i++) {
            int usable = 0;
            int[] fisher = new int[levels.length];
            int[] chiSquared = new int[levels.length];
            for (long seed = 1; seed <= 3000; seed++) {
                long[][] table = independentTable(2, 2, totals[i], seed * 7919L + 1);
                if (hasAnEmptyMargin(table)) {
                    continue;
                }
                usable++;
                double exact = HypothesisTests.fisherExact(table, Alternative.TWO_SIDED).pValue;
                double approximate = HypothesisTests.chiSquaredIndependence(table).pValue;
                for (int k = 0; k < levels.length; k++) {
                    if (exact <= levels[k]) {
                        fisher[k]++;
                    }
                    if (approximate <= levels[k]) {
                        chiSquared[k]++;
                    }
                }
            }
            for (int k = 0; k < levels.length; k++) {
                if (fisher[k] / (double) usable > levels[k]) {
                    fisherAbove++;
                }
                if (chiSquared[k] / (double) usable > levels[k]) {
                    chiSquaredAbove++;
                }
            }
            leastPower = Math.min(leastPower, fisher[2] / (double) usable);
        }
        // measured over the 18 combinations: the exact test is above its own
        // level in none of them, the approximation in nine, worst 0.145 at a
        // nominal 0.10
        assertEquals("the exact test rejected too often", 0, fisherAbove);
        assertTrue("the approximation held its level in " + (18 - chiSquaredAbove) + " of 18 cases",
                chiSquaredAbove >= 6);
        // and it is conservative, not inert: the smallest rate at a nominal 0.10
        // was 0.0454
        assertTrue("the exact test rejects almost nothing : " + leastPower, leastPower > 0.03);
    }

    @Test
    public void testTheTeaTastingTableGivesThePublishedPValues() {
        // Fisher's own example: four cups of each kind, three of each guessed
        // correctly. Of the 70 ways to pick four cups, 16 give three hits and
        // one gives four, so the upper tail is 17/70 = 0.2428571..., and the
        // two-sided value is twice that, the table being its own mirror image
        long[][] tea = { { 3L, 1L }, { 1L, 3L } };
        TestResult greater = HypothesisTests.fisherExact(tea, Alternative.GREATER);
        assertEquals("the published one-sided p-value", 0.24285714285714285, greater.pValue, 1.0e-12);
        assertEquals("the published two-sided p-value", 0.4857142857142857,
                HypothesisTests.fisherExact(tea, Alternative.TWO_SIDED).pValue, 1.0e-12);
        assertEquals("the other tail", 0.9857142857142857,
                HypothesisTests.fisherExact(tea, Alternative.LESS).pValue, 1.0e-12);

        assertEquals("the statistic is the observed count", 3.0, greater.statistic, 0.0);
        assertTrue("the null distribution of the count has no degrees of freedom",
                Double.isNaN(greater.degreesOfFreedom));
        assertEquals("the alternative travels with the result", Alternative.GREATER, greater.alternative);
        assertEquals("the name", "Fisher exact", greater.test);
    }

    @Test
    public void testTheTwoOneSidedTailsOverlapInExactlyTheObservedTable() {
        // both tails contain the observed table, so together they cover the
        // whole null once and that table twice
        double worst = 0.0;
        for (int total : new int[] { 15, 60, 300 }) {
            for (long seed = 1; seed <= 500; seed++) {
                long[][] table = independentTable(2, 2, total, seed * 7919L + 1);
                if (hasAnEmptyMargin(table)) {
                    continue;
                }
                double less = HypothesisTests.fisherExact(table, Alternative.LESS).pValue;
                double greater = HypothesisTests.fisherExact(table, Alternative.GREATER).pValue;
                double observed = tableProbability(table);
                worst = Math.max(worst, Math.abs(less + greater - 1.0 - observed));
            }
        }
        // measured worst case 2.8e-13, the drift of a few hundred log-factorials
        assertTrue("the two tails do not cover the null : " + worst, worst < 1.0e-12);
    }

    @Test
    public void testTheAnswerDoesNotDependOnHowTheTableIsWritten() {
        double worstTranspose = 0.0;
        double worstSwap = 0.0;
        for (int total : new int[] { 15, 60, 300, 1500 }) {
            for (long seed = 1; seed <= 400; seed++) {
                long[][] table = independentTable(2, 2, total, seed * 7919L + 1);
                if (hasAnEmptyMargin(table)) {
                    continue;
                }
                long[][] transposed = { { table[0][0], table[1][0] }, { table[0][1], table[1][1] } };
                long[][] swapped = { { table[0][1], table[0][0] }, { table[1][1], table[1][0] } };
                worstTranspose = Math.max(worstTranspose,
                        Math.abs(HypothesisTests.fisherExact(table, Alternative.TWO_SIDED).pValue
                                - HypothesisTests.fisherExact(transposed, Alternative.TWO_SIDED).pValue));
                // exchanging the columns turns a large first cell into a small
                // one, so the two one-sided tests change places
                worstSwap = Math.max(worstSwap,
                        Math.abs(HypothesisTests.fisherExact(table, Alternative.LESS).pValue
                                - HypothesisTests.fisherExact(swapped, Alternative.GREATER).pValue));
            }
        }
        // transposing reaches the same probability by a different sum of
        // logarithms, so it is 5.1e-13 rather than exact; exchanging the columns
        // reuses the same terms and came out at 8.9e-16
        assertTrue("transposing changed the answer by " + worstTranspose, worstTranspose < 1.0e-12);
        assertTrue("exchanging the columns changed the answer by " + worstSwap, worstSwap < 1.0e-14);
    }

    @Test
    public void testItAgreesWithTheFactorialFormula() {
        // an independent route to the same number: the probability of a table
        // with fixed margins is r1! r2! c1! c2! / (N! a! b! c! d!), which shares
        // no line of code with the hypergeometric mass function summed inside
        Alternative[] alternatives = { Alternative.TWO_SIDED, Alternative.LESS, Alternative.GREATER };
        double worst = 0.0;
        for (int total : new int[] { 12, 40, 150, 600 }) {
            for (long seed = 1; seed <= 300; seed++) {
                long[][] table = independentTable(2, 2, total, seed * 7919L + 1);
                if (hasAnEmptyMargin(table)) {
                    continue;
                }
                for (int i = 0; i < alternatives.length; i++) {
                    double shipped = HypothesisTests.fisherExact(table, alternatives[i]).pValue;
                    double formula = byTheFactorialFormula(table, alternatives[i]);
                    worst = Math.max(worst, Math.abs(shipped - formula) / formula);
                }
            }
        }
        // measured worst relative difference 1.2e-12
        assertTrue("the two routes disagree by " + worst, worst < 5.0e-12);
    }

    @Test
    public void testTheFarTailIsNotRoundedToZero() {
        // the p-value of a perfectly separated table is 1 / C(2n, n) doubled,
        // which leaves the range a p-value has if it is ever formed as one
        // minus a distribution function
        long[][] separated = { { 500L, 0L }, { 0L, 500L } };
        double two = HypothesisTests.fisherExact(separated, Alternative.TWO_SIDED).pValue;
        double greater = HypothesisTests.fisherExact(separated, Alternative.GREATER).pValue;
        assertTrue("the tail collapsed to " + two, two > 0.0);
        assertTrue("the tail is not small : " + two, two < 1.0e-290);
        // the table is its own mirror image, so the two-sided sum is the
        // one-sided one twice -- the same two terms added in the same order
        assertEquals("the two-sided value is not the one-sided one doubled", two / 2.0, greater, 0.0);

        // and one that is merely lopsided still lands far below any level
        long[][] lopsided = { { 90L, 10L }, { 10L, 90L } };
        assertTrue("a plainly dependent table was not rejected",
                HypothesisTests.fisherExact(lopsided, Alternative.TWO_SIDED).pValue < 1.0e-30);
    }

    @Test
    public void testFisherRejectsWhatItCannotTest() {
        long[][] ok = { { 3L, 1L }, { 1L, 3L } };
        rejectsFisher("a null table", null, Alternative.TWO_SIDED);
        rejectsFisher("a null alternative", ok, null);
        rejectsFisher("a table with three rows", new long[][] { { 1L, 1L }, { 1L, 1L }, { 1L, 1L } },
                Alternative.TWO_SIDED);
        rejectsFisher("a table with three columns", new long[][] { { 1L, 1L, 1L }, { 1L, 1L, 1L } },
                Alternative.TWO_SIDED);
        rejectsFisher("a null row", new long[][] { { 1L, 1L }, null }, Alternative.TWO_SIDED);
        rejectsFisher("a negative count", new long[][] { { -1L, 2L }, { 3L, 4L } }, Alternative.TWO_SIDED);
        rejectsFisher("an empty row", new long[][] { { 0L, 0L }, { 3L, 4L } }, Alternative.TWO_SIDED);
        rejectsFisher("an empty column", new long[][] { { 0L, 2L }, { 0L, 4L } }, Alternative.TWO_SIDED);
        rejectsFisher("a count past the int range",
                new long[][] { { 3_000_000_000L, 1L }, { 1L, 1L } }, Alternative.TWO_SIDED);
        rejectsFisher("a total past the int range",
                new long[][] { { 2_000_000_000L, 1_000_000_000L }, { 1L, 1L } }, Alternative.TWO_SIDED);
    }

    private static void rejectsFisher(String what, long[][] table, Alternative alternative) {
        try {
            HypothesisTests.fisherExact(table, alternative);
            fail("fisherExact accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    /** A margin of zero leaves nothing to test, and both tests reject it. */
    private static boolean hasAnEmptyMargin(long[][] t) {
        return t[0][0] + t[0][1] == 0L || t[1][0] + t[1][1] == 0L || t[0][0] + t[1][0] == 0L
                || t[0][1] + t[1][1] == 0L;
    }

    /** {@code r1! r2! c1! c2! / (N! a! b! c! d!)}, the textbook formula. */
    private static double tableProbability(long[][] t) {
        int a = (int) t[0][0];
        int b = (int) t[0][1];
        int c = (int) t[1][0];
        int d = (int) t[1][1];
        return Math.exp(Arithmetic.logFactorial(a + b) + Arithmetic.logFactorial(c + d)
                + Arithmetic.logFactorial(a + c) + Arithmetic.logFactorial(b + d)
                - Arithmetic.logFactorial(a + b + c + d) - Arithmetic.logFactorial(a)
                - Arithmetic.logFactorial(b) - Arithmetic.logFactorial(c) - Arithmetic.logFactorial(d));
    }

    private static double byTheFactorialFormula(long[][] t, Alternative alternative) {
        int a = (int) t[0][0];
        int rowTotal = (int) (t[0][0] + t[0][1]);
        int columnTotal = (int) (t[0][0] + t[1][0]);
        int total = (int) (t[0][0] + t[0][1] + t[1][0] + t[1][1]);
        double observed = tableProbability(t);
        double p = 0.0;
        for (int k = Math.max(0, rowTotal + columnTotal - total); k <= Math.min(rowTotal, columnTotal); k++) {
            long[][] candidate = { { k, rowTotal - k },
                    { columnTotal - k, total - rowTotal - columnTotal + k } };
            double mass = tableProbability(candidate);
            boolean take;
            if (alternative == Alternative.LESS) {
                take = k <= a;
            } else if (alternative == Alternative.GREATER) {
                take = k >= a;
            } else {
                take = mass <= observed * (1.0 + 1.0e-9);
            }
            if (take) {
                p += mass;
            }
        }
        return Math.min(1.0, p);
    }

    // ----------------------------------------- two-sample Kolmogorov-Smirnov --

    @Test
    public void testTheTwoSamplePValueIsCalibratedUnderTheNullHypothesis() {
        // both samples from the same distribution, so the p-value may not
        // fall below its level more often than the level allows
        double[] levels = { 0.01, 0.05, 0.50 };
        int[][] shapes = { { 5, 5 }, { 8, 12 }, { 30, 30 }, { 50, 70 } };
        Uniform uniform = new Uniform(0.0, 1.0);
        for (int s = 0; s < shapes.length; s++) {
            double[] rates = rejectionRates(shapes[s][0], shapes[s][1], uniform, uniform, 4000, levels);
            for (int k = 0; k < rates.length; k++) {
                double level = levels[k % levels.length];
                double noise = 2.5 * Math.sqrt(level * (1.0 - level) / 4000.0);
                assertTrue("m=" + shapes[s][0] + " n=" + shapes[s][1] + " rejected " + rates[k]
                        + " of the time at level " + level, rates[k] <= level + noise);
            }
            if (shapes[s][0] >= 50) {
                // and at the largest exact shape it is not merely inert:
                // measured 0.4835, 0.4928 and 0.4935 at a nominal 0.5
                for (int k = 2; k < rates.length; k += levels.length) {
                    assertTrue("the median p-value is off at " + rates[k],
                            rates[k] > 0.40 && rates[k] < 0.60);
                }
            }
        }
        // the asymptotic branch, above the exact limit: 0.4333, 0.4533, 0.5367
        double[] big = rejectionRates(1200, 1200, uniform, uniform, 300, levels);
        for (int k = 2; k < big.length; k += levels.length) {
            assertTrue("the approximation is off at " + big[k], big[k] > 0.40 && big[k] < 0.60);
        }
    }

    @Test
    public void testTheTwoSampleStatisticIsTheDefinition() {
        // the merge walk against the largest gap read off every pooled point
        // by counting, which is the definition and quadratic
        for (int[] size : new int[][] { { 6, 9 }, { 40, 25 }, { 200, 200 } }) {
            for (long seed = 1; seed <= 100; seed++) {
                double[] x = sample(new Normal(0.0, 1.0), size[0], seed * 7919L + 1);
                double[] y = sample(new Normal(0.3, 2.0), size[1], seed * 104729L + 7);
                double[] byHand = largestGaps(x, y);
                // the same integers divided by the same integers, so bit for bit
                assertEquals("D+", byHand[0],
                        HypothesisTests.kolmogorovSmirnovTwoSample(x, y, Alternative.LESS).statistic, 0.0);
                assertEquals("D-", byHand[1],
                        HypothesisTests.kolmogorovSmirnovTwoSample(x, y, Alternative.GREATER).statistic,
                        0.0);
                assertEquals("D", Math.max(byHand[0], byHand[1]),
                        HypothesisTests.kolmogorovSmirnovTwoSample(x, y, Alternative.TWO_SIDED).statistic,
                        0.0);
            }
        }
    }

    @Test
    public void testTheTwoSampleTestOnlySeesTheOrderOfTheObservations() {
        // shifting both samples by the same amount and scaling both by the
        // same power of two leaves every interleaving where it was
        for (int[] size : new int[][] { { 7, 9 }, { 30, 40 }, { 300, 300 } }) {
            for (long seed = 1; seed <= 100; seed++) {
                double[] x = sample(new Normal(0.0, 1.0), size[0], seed * 7919L + 1);
                double[] y = sample(new Normal(0.4, 1.0), size[1], seed * 104729L + 7);
                double[] movedX = new double[x.length];
                double[] movedY = new double[y.length];
                for (int i = 0; i < x.length; i++) {
                    movedX[i] = (x[i] + 1000.0) * 0x1.0p900;
                }
                for (int i = 0; i < y.length; i++) {
                    movedY[i] = (y[i] + 1000.0) * 0x1.0p900;
                }
                for (Alternative alternative : Alternative.values()) {
                    TestResult before = HypothesisTests.kolmogorovSmirnovTwoSample(x, y, alternative);
                    TestResult after = HypothesisTests.kolmogorovSmirnovTwoSample(movedX, movedY,
                            alternative);
                    assertEquals("statistic", before.statistic, after.statistic, 0.0);
                    assertEquals("p-value", before.pValue, after.pValue, 0.0);
                }
            }
        }
    }

    @Test
    public void testTheExactTailAgreesWithEveryInterleaving() {
        // the independent algorithm: walk all C(m+n, m) interleavings and
        // count the ones whose statistic reaches d
        double worst = 0.0;
        for (int[] size : new int[][] { { 3, 5 }, { 4, 6 }, { 5, 8 }, { 7, 7 }, { 6, 10 } }) {
            int m = size[0];
            int n = size[1];
            for (int k = 1; k <= m * n; k++) {
                double d = k / (double) (m * n);
                double counted = countInterleavings(m, n, d, false);
                if (counted > 0.0) {
                    worst = Math.max(worst,
                            Math.abs(KolmogorovSmirnovTwoSample.barFExact(m, n, d) - counted) / counted);
                }
                double countedOneSided = countInterleavings(m, n, d, true);
                if (countedOneSided > 0.0) {
                    worst = Math.max(worst,
                            Math.abs(KolmogorovSmirnovTwoSample.barFExactOneSided(m, n, d)
                                    - countedOneSided) / countedOneSided);
                }
            }
        }
        // measured worst case 3.3e-15 over the whole range of every shape
        assertTrue("the recursion and the enumeration differ by " + worst, worst < 1.0e-13);
    }

    @Test
    public void testASeparatedPairGetsTheCombinatorialPValue() {
        // when every x lies below every y, exactly two of the C(m+n, m)
        // interleavings are as extreme, so the p-value is 2 / C(m+n, m)
        for (int[] size : new int[][] { { 8, 8 }, { 10, 10 }, { 15, 15 }, { 20, 30 } }) {
            int m = size[0];
            int n = size[1];
            double[] x = new double[m];
            double[] y = new double[n];
            for (int i = 0; i < m; i++) {
                x[i] = i;
            }
            for (int i = 0; i < n; i++) {
                y[i] = m + i;
            }
            TestResult r = HypothesisTests.kolmogorovSmirnovTwoSample(x, y, Alternative.TWO_SIDED);
            assertEquals("the statistic is one", 1.0, r.statistic, 0.0);
            double reference = 2.0 / Arithmetic.binomial((long) (m + n), (long) m);
            // measured worst relative deviation 5.4e-15
            assertTrue("p = " + r.pValue + " against 2/C(m+n,m) = " + reference,
                    Math.abs(r.pValue - reference) / reference < 5.0e-14);
        }
    }

    @Test
    public void testTiesMoveTheTestOntoTheApproximationAndSaySo() {
        double[] other = { 1.5, 2.5, 3.5, 4.5, 5.5 };
        double[] withATie = { 1.0, 2.0, 2.0, 3.0, 4.0 };
        double[] withoutOne = { 1.0, 2.0, 2.25, 3.0, 4.0 };
        assertTrue("a repeated value must not be tested exactly",
                HypothesisTests.kolmogorovSmirnovTwoSample(withATie, other, Alternative.TWO_SIDED).test
                        .endsWith("asymptotic"));
        assertTrue("distinct values must be tested exactly",
                HypothesisTests.kolmogorovSmirnovTwoSample(withoutOne, other, Alternative.TWO_SIDED).test
                        .endsWith("exact"));
        // a value shared between the samples is a tie as much as one repeated
        // inside a sample
        assertTrue("a shared value must not be tested exactly", HypothesisTests
                .kolmogorovSmirnovTwoSample(new double[] { 1.0, 2.0, 3.0 }, new double[] { 2.0, 4.0 },
                        Alternative.TWO_SIDED).test.endsWith("asymptotic"));
        // and one that repeats past the point where the shorter sample runs
        // out, which the merge walk never reaches
        assertTrue("a repetition in the tail must not be tested exactly", HypothesisTests
                .kolmogorovSmirnovTwoSample(new double[] { 1.0, 2.0, 3.0, 9.0, 9.0 },
                        new double[] { 1.5, 2.5 }, Alternative.TWO_SIDED).test.endsWith("asymptotic"));
    }

    @Test
    public void testTheTwoSampleAlternativesPointAtTheSample() {
        // GREATER asks whether x is stochastically larger, so a sample shifted
        // upwards has to be found by GREATER and not by LESS
        double[] y = sample(new Normal(0.0, 1.0), 60, 999L);
        double[] low = sample(new Normal(-0.8, 1.0), 60, 11L);
        double[] high = sample(new Normal(0.8, 1.0), 60, 11L);
        assertTrue("a sample shifted down is not stochastically smaller",
                HypothesisTests.kolmogorovSmirnovTwoSample(low, y, Alternative.LESS).pValue < 0.001);
        assertTrue("and it is certainly not larger",
                HypothesisTests.kolmogorovSmirnovTwoSample(low, y, Alternative.GREATER).pValue > 0.5);
        assertTrue("a sample shifted up is not stochastically larger",
                HypothesisTests.kolmogorovSmirnovTwoSample(high, y, Alternative.GREATER).pValue < 0.001);
        assertTrue("two samples from the same distribution look alike",
                HypothesisTests.kolmogorovSmirnovTwoSample(sample(new Normal(0.0, 1.0), 60, 11L), y,
                        Alternative.TWO_SIDED).pValue > 0.1);
    }

    @Test
    public void testTheApproximationDoesNotUnderstateTheExactTail() {
        // it may be too large -- that only costs power -- but a p-value that
        // is too small would turn an unremarkable result into a significant one
        double worst = Double.MAX_VALUE;
        for (int[] size : new int[][] { { 40, 40 }, { 100, 100 }, { 30, 300 }, { 60, 200 } }) {
            int m = size[0];
            int n = size[1];
            for (int k = 1; k <= m * n; k += Math.max(1, m * n / 60)) {
                double d = k / (double) (m * n);
                double exact = KolmogorovSmirnovTwoSample.barFExact(m, n, d);
                if (exact > 1.0e-250) {
                    worst = Math.min(worst, KolmogorovSmirnovTwoSample.barFAsymptotic(m, n, d) / exact);
                }
            }
        }
        // measured smallest ratio 0.999958, at m = n = 100
        assertTrue("the approximation fell to " + worst + " of the exact tail", worst > 0.999);
    }

    @Test
    public void testTheTwoSidedTailIsNoLargerThanBothOneSidedOnes() {
        double worstExcess = 0.0;
        for (int[] size : new int[][] { { 5, 5 }, { 12, 20 }, { 60, 60 }, { 100, 100 } }) {
            int m = size[0];
            int n = size[1];
            for (int k = 1; k <= m * n; k += Math.max(1, m * n / 40)) {
                double d = k / (double) (m * n);
                double bound = 2.0 * KolmogorovSmirnovTwoSample.barFExactOneSided(m, n, d);
                if (bound > 1.0e-250) {
                    worstExcess = Math.max(worstExcess,
                            (KolmogorovSmirnovTwoSample.barFExact(m, n, d) - bound) / bound);
                }
            }
        }
        // the two events overlap only where neither can happen, so the sum is
        // an upper bound; measured excess 2.6e-14, which is round-off
        assertTrue("the two-sided tail exceeds both one-sided ones by " + worstExcess,
                worstExcess < 1.0e-12);
    }

    @Test
    public void testTheTwoSampleTestRejectsWhatItCannotTest() {
        double[] ok = { 1.0, 2.0, 3.0 };
        rejectsKsTwo("a null first sample", null, ok, Alternative.TWO_SIDED);
        rejectsKsTwo("a null second sample", ok, null, Alternative.TWO_SIDED);
        rejectsKsTwo("an empty sample", new double[0], ok, Alternative.TWO_SIDED);
        rejectsKsTwo("an infinite observation", new double[] { 1.0, Double.POSITIVE_INFINITY }, ok,
                Alternative.TWO_SIDED);
        rejectsKsTwo("a NaN observation", ok, new double[] { Double.NaN }, Alternative.TWO_SIDED);
        rejectsKsTwo("a null alternative", ok, ok, null);

        // and the distribution itself
        rejectsExactTail("a sample of no observations", 0, 5, 0.5);
        rejectsExactTail("a second sample of no observations", 5, 0, 0.5);
        rejectsExactTail("a statistic that is not a number", 5, 5, Double.NaN);
        rejectsExactTail("a lattice past the exact limit", 2000, 2000, 0.5);
    }

    private static void rejectsKsTwo(String what, double[] x, double[] y, Alternative alternative) {
        try {
            HypothesisTests.kolmogorovSmirnovTwoSample(x, y, alternative);
            fail("kolmogorovSmirnovTwoSample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsExactTail(String what, int m, int n, double d) {
        try {
            KolmogorovSmirnovTwoSample.barFExact(m, n, d);
            fail("barFExact accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    /** How often each alternative falls below each level, in that order. */
    private static double[] rejectionRates(int m, int n, ContinuousDistribution first,
            ContinuousDistribution second, int replications, double[] levels) {
        Alternative[] alternatives = Alternative.values();
        double[] rates = new double[alternatives.length * levels.length];
        for (long seed = 1; seed <= replications; seed++) {
            double[] x = sample(first, m, seed * 7919L + 1);
            double[] y = sample(second, n, seed * 104729L + 7);
            for (int a = 0; a < alternatives.length; a++) {
                double p = HypothesisTests.kolmogorovSmirnovTwoSample(x, y, alternatives[a]).pValue;
                for (int k = 0; k < levels.length; k++) {
                    if (p <= levels[k]) {
                        rates[a * levels.length + k] += 1.0;
                    }
                }
            }
        }
        for (int i = 0; i < rates.length; i++) {
            rates[i] /= replications;
        }
        return rates;
    }

    /** {@code {D+, D-}} read off every pooled observation by counting. */
    private static double[] largestGaps(double[] x, double[] y) {
        double[] pooled = new double[x.length + y.length];
        System.arraycopy(x, 0, pooled, 0, x.length);
        System.arraycopy(y, 0, pooled, x.length, y.length);
        double dPlus = 0.0;
        double dMinus = 0.0;
        for (int i = 0; i < pooled.length; i++) {
            double below = 0.0;
            for (int k = 0; k < x.length; k++) {
                if (x[k] <= pooled[i]) {
                    below++;
                }
            }
            double belowY = 0.0;
            for (int k = 0; k < y.length; k++) {
                if (y[k] <= pooled[i]) {
                    belowY++;
                }
            }
            double gap = below / x.length - belowY / y.length;
            dPlus = Math.max(dPlus, gap);
            dMinus = Math.max(dMinus, -gap);
        }
        return new double[] { dPlus, dMinus };
    }

    private static double interleavingsSeen;
    private static double interleavingsExtreme;

    /** The share of the interleavings whose statistic reaches {@code d}. */
    private static double countInterleavings(int m, int n, double d, boolean oneSided) {
        interleavingsSeen = 0.0;
        interleavingsExtreme = 0.0;
        walkInterleavings(0, 0, m, n, 0.0, 0.0, d, oneSided);
        return interleavingsExtreme / interleavingsSeen;
    }

    private static void walkInterleavings(int i, int k, int m, int n, double dPlus, double dMinus,
            double d, boolean oneSided) {
        if (i == m && k == n) {
            interleavingsSeen++;
            if ((oneSided ? dPlus : Math.max(dPlus, dMinus)) >= d - 1.0e-12) {
                interleavingsExtreme++;
            }
            return;
        }
        if (i < m) {
            double gap = (i + 1.0) / m - k / (double) n;
            walkInterleavings(i + 1, k, m, n, Math.max(dPlus, gap), Math.max(dMinus, -gap), d, oneSided);
        }
        if (k < n) {
            double gap = i / (double) m - (k + 1.0) / n;
            walkInterleavings(i, k + 1, m, n, Math.max(dPlus, gap), Math.max(dMinus, -gap), d, oneSided);
        }
    }

    // ------------------------------------------------------------ Lilliefors --

    @Test
    public void testLillieforsHoldsItsLevelWhereThePlainTestCollapses() {
        // the whole reason this test exists. Fit a normal distribution to a
        // normal sample and test the sample against the fit: the Lilliefors
        // p-value is calibrated, the ordinary Kolmogorov-Smirnov p-value is so
        // far off that it never rejects at all
        for (int n : new int[] { 5, 10, 20, 50 }) {
            int replications = 600;
            int fittedRejections = 0;
            int plainRejections = 0;
            for (long seed = 1; seed <= replications; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                if (HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL, 2000, 7L).test.pValue <= 0.05) {
                    fittedRejections++;
                }
                ParNormal fit = MLE.getNormalMLE(x);
                if (HypothesisTests.kolmogorovSmirnov(x, new Normal(fit.mean, fit.stdDev),
                        Alternative.TWO_SIDED).pValue <= 0.05) {
                    plainRejections++;
                }
            }
            double fittedRate = fittedRejections / (double) replications;
            // measured 0.0467, 0.0450, 0.0483 and 0.0600 at a nominal 0.05
            assertTrue("n=" + n + " the fitted test rejected " + fittedRate + " of the time",
                    fittedRate > 0.02 && fittedRate < 0.09);
            // and the plain test rejected none of 2400 samples
            assertEquals("n=" + n + " the plain test rejected something", 0, plainRejections);
        }
    }

    @Test
    public void testTheSimulatedCriticalValuesAgreeWithTheLiterature() {
        // Stephens (1974) gives the upper tail of D (sqrt(n) - 0.01 +
        // 0.85/sqrt(n)) as 0.895 at five per cent and 1.035 at one, for the
        // normal family fitted in both parameters. Those points are for the
        // n-1 standard deviation and this package fits by maximum likelihood,
        // so they have to be carried across by the Bessel factor first
        for (int n : new int[] { 50, 100, 200 }) {
            double modifier = Math.sqrt(n) - 0.01 + 0.85 / Math.sqrt(n);
            double bessel = Math.sqrt(n / (n - 1.0));
            double at05 = Lilliefors.barF(Lilliefors.Family.NORMAL, n, 0.895 / modifier * bessel, 20000,
                    5L);
            double at01 = Lilliefors.barF(Lilliefors.Family.NORMAL, n, 1.035 / modifier * bessel, 20000,
                    5L);
            // measured 0.0469, 0.0498, 0.0529 and 0.0092, 0.0125, 0.0119;
            // the modification is a fit to moderate n and drifts at small n,
            // which is why this starts at 50
            assertTrue("n=" + n + " the five per cent point came out at " + at05,
                    at05 > 0.042 && at05 < 0.059);
            assertTrue("n=" + n + " the one per cent point came out at " + at01,
                    at01 > 0.007 && at01 < 0.015);
        }
    }

    @Test
    public void testTheSimulatedPValueIsReproducibleAndSaysHowUncertainItIs() {
        double[] x = sample(new Normal(0.0, 1.0), 40, 11L);
        SimulatedTestResult first = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL);
        SimulatedTestResult second = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL);
        assertEquals("the short form is a function of the sample", first.test.pValue,
                second.test.pValue, 0.0);
        assertEquals("and it says how many replications went into it", Lilliefors.DEFAULT_REPLICATIONS,
                first.replications);
        assertEquals("the standard error is the binomial one", Math.sqrt(
                first.test.pValue * (1.0 - first.test.pValue) / first.replications),
                first.monteCarloStandardError, 0.0);

        // another seed moves the answer, and more replications move it less
        double[] spread = new double[3];
        int[] sizes = { 1000, 4000, 16000 };
        for (int i = 0; i < sizes.length; i++) {
            double lo = 1.0;
            double hi = 0.0;
            for (long seed = 1; seed <= 12; seed++) {
                double p = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL, sizes[i], seed)
                        .test.pValue;
                lo = Math.min(lo, p);
                hi = Math.max(hi, p);
            }
            spread[i] = hi - lo;
        }
        // measured 0.0340, 0.0235 and 0.0119 over twelve seeds
        assertTrue("twelve seeds gave the same answer", spread[0] > 0.0);
        assertTrue("sixteen times the work did not halve the spread : " + spread[2] + " against "
                + spread[0], spread[2] < 0.5 * spread[0]);
    }

    @Test
    public void testTheFittedStatisticIsTheDistanceToTheFit() {
        // the statistic the simulation measures on its own samples has to be
        // the one the test measures on the real one, or the comparison means
        // nothing
        for (int n : new int[] { 5, 20, 100 }) {
            for (long seed = 1; seed <= 200; seed++) {
                double[] x = sample(new Normal(1.0, 3.0), n, seed * 7919L + 1);
                ParNormal fit = MLE.getNormalMLE(x);
                Normal fitted = new Normal(fit.mean, fit.stdDev);
                assertEquals("the fitted distance",
                        HypothesisTests.kolmogorovSmirnov(x, fitted, Alternative.TWO_SIDED).statistic,
                        Lilliefors.statistic(Lilliefors.Family.NORMAL, x), 0.0);
                // the same for the second statistic, which is written out
                // twice: once in the facade and once in the simulation
                assertEquals("the fitted tail-weighted distance",
                        HypothesisTests.andersonDarling(x, fitted).statistic,
                        Lilliefors.statistic(Lilliefors.Statistic.ANDERSON_DARLING,
                                Lilliefors.Family.NORMAL, x), 0.0);
            }
        }
    }

    @Test
    public void testTheLogNormalFamilyIsTheNormalOneOnTheLogarithms() {
        // exponentiating a sample and fitting a log-normal is fitting a normal
        // to the logarithms, and the statistic is invariant under a monotone
        // transformation, so both p-values have to come out identical
        for (int n : new int[] { 6, 30 }) {
            for (long seed = 1; seed <= 100; seed++) {
                double[] normal = sample(new Normal(0.5, 1.5), n, seed * 7919L + 1);
                double[] exponentiated = new double[n];
                for (int i = 0; i < n; i++) {
                    exponentiated[i] = Math.exp(normal[i]);
                }
                assertEquals("the p-value",
                        HypothesisTests.lilliefors(normal, Lilliefors.Family.NORMAL, 500, 5L).test.pValue,
                        HypothesisTests.lilliefors(exponentiated, Lilliefors.Family.LOGNORMAL, 500, 5L)
                                .test.pValue,
                        0.0);
                assertEquals("the tail-weighted p-value",
                        HypothesisTests.andersonDarling(normal, Lilliefors.Family.NORMAL, 500, 5L)
                                .test.pValue,
                        HypothesisTests.andersonDarling(exponentiated, Lilliefors.Family.LOGNORMAL, 500, 5L)
                                .test.pValue,
                        0.0);
            }
        }
    }

    @Test
    public void testLillieforsFindsADistributionThatIsNotTheOneItFits() {
        int heavyTails = 0;
        int skewed = 0;
        for (long seed = 1; seed <= 200; seed++) {
            if (HypothesisTests.lilliefors(sample(new StudentT(3.0), 100, seed * 7919L + 1),
                    Lilliefors.Family.NORMAL, 1000, 5L).test.pValue <= 0.05) {
                heavyTails++;
            }
            if (HypothesisTests.lilliefors(sample(new Exponential(1.0), 100, seed * 104729L + 7),
                    Lilliefors.Family.NORMAL, 1000, 5L).test.pValue <= 0.05) {
                skewed++;
            }
        }
        // measured 0.620 against t(3) and 1.000 against a skewed sample
        assertTrue("t(3) was found only " + heavyTails + " times in 200", heavyTails > 80);
        assertTrue("an exponential sample was found only " + skewed + " times in 200", skewed > 180);
    }

    @Test
    public void testLillieforsRejectsWhatItCannotTest() {
        double[] ok = { 1.0, 2.0, 3.0, 4.0 };
        rejectsLilliefors("a null sample", null, Lilliefors.Family.NORMAL, 100, 1L);
        rejectsLilliefors("a null family", ok, null, 100, 1L);
        rejectsLilliefors("two observations and two parameters", new double[] { 1.0, 2.0 },
                Lilliefors.Family.NORMAL, 100, 1L);
        rejectsLilliefors("a NaN observation", new double[] { 1.0, 2.0, Double.NaN },
                Lilliefors.Family.NORMAL, 100, 1L);
        rejectsLilliefors("no replications", ok, Lilliefors.Family.NORMAL, 0, 1L);
        rejectsLilliefors("a negative observation in a positive family",
                new double[] { 1.0, 2.0, -3.0 }, Lilliefors.Family.EXPONENTIAL, 100, 1L);
        rejectsLilliefors("a zero in a log-normal sample", new double[] { 1.0, 2.0, 0.0 },
                Lilliefors.Family.LOGNORMAL, 100, 1L);

        // and the null distribution on its own
        try {
            Lilliefors.barF(Lilliefors.Family.NORMAL, 2, 0.3, 100, 1L);
            fail("barF accepted a sample of two");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            Lilliefors.barF(Lilliefors.Family.NORMAL, 10, Double.NaN, 100, 1L);
            fail("barF accepted a statistic that is not a number");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsLilliefors(String what, double[] x, Lilliefors.Family family,
            int replications, long seed) {
        try {
            HypothesisTests.lilliefors(x, family, replications, seed);
            fail("lilliefors accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // --------------------------------------------- Anderson-Darling, fitted --

    @Test
    public void testTheFittedAndersonDarlingHoldsItsLevelWhereTheFullySpecifiedOneCollapses() {
        for (int n : new int[] { 5, 10, 30, 100 }) {
            int reps = 1500;
            int fitted = 0;
            int specified = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                if (HypothesisTests.andersonDarling(x, Lilliefors.Family.NORMAL, 2000, 7L).test
                        .rejectsAt(0.05)) {
                    fitted++;
                }
                // the same statistic, read against the same fit, but through
                // the null distribution that assumes the fit was not looked at
                ParNormal fit = MLE.getNormalMLE(x);
                if (HypothesisTests.andersonDarling(x, new Normal(fit.mean, fit.stdDev)).rejectsAt(0.05)) {
                    specified++;
                }
            }
            // measured: 0.0440 .. 0.0547 fitted, and 0 of 6000 for the other
            assertRate("n=" + n + " fitted", fitted / (double) reps, 0.05, 0.018);
            assertTrue("n=" + n + ": the fully specified test rejected " + specified + " of " + reps
                    + " true normals", specified == 0);
        }
    }

    @Test
    public void testTheFittedAndersonDarlingReproducesStephensPointsUnderTheirOwnConvention() {
        // Stephens (1974) tabulates A2 (1 + 0.75/n + 2.25/n^2) for a normal
        // fitted in both parameters. The table is built on the n - 1 standard
        // deviation and math.stats.mle fits by maximum likelihood, so the
        // construction is checked by simulating it the way the table was
        // built, and the difference the library convention makes is measured
        // rather than assumed
        double[] levels = { 0.10, 0.05, 0.025, 0.01 };
        double[] points = { 0.631, 0.752, 0.873, 1.035 };
        Normal standard = new Normal(0.0, 1.0);
        int reps = 40000;
        for (int n : new int[] { 20, 100 }) {
            double modifier = 1.0 + 0.75 / n + 2.25 / (n * (double) n);
            double[] drawn = new double[reps];
            for (int b = 0; b < reps; b++) {
                drawn[b] = besselFittedAndersonDarling(sample(standard, n, (b + 1L) * 7919L + 1)) * modifier;
            }
            Arrays.sort(drawn);
            for (int i = 0; i < levels.length; i++) {
                double simulated = upperTail(drawn, points[i]);
                // measured over three independent runs: 0.0021 at worst
                assertEquals("n=" + n + " at the " + levels[i] + " point", levels[i], simulated, 0.005);
            }
        }

        // and now the convention. At n = 20 the library null, fitted by
        // maximum likelihood, puts more mass beyond the tabulated points than
        // the table says: reading a maximum likelihood fit off that table
        // rejects about eleven percent of the time at the ten percent point
        double modifier = 1.0 + 0.75 / 20.0 + 2.25 / 400.0;
        double atTen = Lilliefors.barF(Lilliefors.Statistic.ANDERSON_DARLING, Lilliefors.Family.NORMAL, 20,
                0.631 / modifier, 40000, 3L);
        double atFive = Lilliefors.barF(Lilliefors.Statistic.ANDERSON_DARLING, Lilliefors.Family.NORMAL, 20,
                0.752 / modifier, 40000, 3L);
        // measured: 0.1106 and 0.0553
        assertTrue("the maximum likelihood null gives " + atTen + " at the 0.10 point",
                atTen > 0.10 && atTen < 0.13);
        assertTrue("the maximum likelihood null gives " + atFive + " at the 0.05 point",
                atFive > 0.05 && atFive < 0.065);

        // by n = 100 the two conventions have almost met
        double wide = 1.0 + 0.75 / 100.0 + 2.25 / 10000.0;
        double atFiveWide = Lilliefors.barF(Lilliefors.Statistic.ANDERSON_DARLING, Lilliefors.Family.NORMAL,
                100, 0.752 / wide, 40000, 3L);
        // measured: 0.0504
        assertEquals("at n = 100 the convention hardly matters", 0.05, atFiveWide, 0.006);
    }

    /** {@code A_n^2} against a normal fitted with the {@code n - 1} standard deviation. */
    private static double besselFittedAndersonDarling(double[] x) {
        int n = x.length;
        double mean = 0.0;
        for (int i = 0; i < n; i++) {
            mean += x[i];
        }
        mean /= n;
        double sumSquares = 0.0;
        for (int i = 0; i < n; i++) {
            double deviation = x[i] - mean;
            sumSquares += deviation * deviation;
        }
        return HypothesisTests.andersonDarling(x, new Normal(mean, Math.sqrt(sumSquares / (n - 1.0))))
                .statistic;
    }

    /** The share of a sorted array at or beyond {@code at}. */
    private static double upperTail(double[] sorted, double at) {
        int count = 0;
        for (int i = sorted.length - 1; i >= 0 && sorted[i] >= at; i--) {
            count++;
        }
        return count / (double) sorted.length;
    }

    @Test
    public void testTheFittedAndersonDarlingSeesAHeavyTailSoonerThanLilliefors() {
        for (int n : new int[] { 20, 50, 200 }) {
            int reps = 800;
            int andersonDarling = 0;
            int lilliefors = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] heavy = sample(new StudentT(3.0), n, seed * 7919L + 1);
                if (HypothesisTests.andersonDarling(heavy, Lilliefors.Family.NORMAL, 1000, 5L).test
                        .rejectsAt(0.05)) {
                    andersonDarling++;
                }
                if (HypothesisTests.lilliefors(heavy, Lilliefors.Family.NORMAL, 1000, 5L).test
                        .rejectsAt(0.05)) {
                    lilliefors++;
                }
            }
            // measured: 0.2975 against 0.2300, 0.5638 against 0.4675, 0.9875
            // against 0.9413
            assertTrue("n=" + n + ": the tail-weighted test found " + andersonDarling + " of " + reps
                    + " and the other " + lilliefors, andersonDarling > lilliefors);
        }
    }

    @Test
    public void testTheFittedAndersonDarlingIsReproducibleAndSaysHowUncertainItIs() {
        double[] x = sample(new Normal(0.0, 1.0), 40, 11L);
        SimulatedTestResult first = HypothesisTests.andersonDarling(x, Lilliefors.Family.NORMAL);
        SimulatedTestResult second = HypothesisTests.andersonDarling(x, Lilliefors.Family.NORMAL);
        assertEquals("the short form is a function of the sample", first.test.pValue, second.test.pValue,
                0.0);
        assertEquals("and it says how many replications went into it", Lilliefors.DEFAULT_REPLICATIONS,
                first.replications);
        assertEquals("the standard error is the binomial one",
                Math.sqrt(first.test.pValue * (1.0 - first.test.pValue) / first.replications),
                first.monteCarloStandardError, 0.0);
        assertTrue("the result names the statistic and the family: " + first.test.test,
                first.test.test.contains("Anderson-Darling") && first.test.test.contains("normal"));

        // the two fitted tests measure different things, so they must not
        // return the same statistic
        assertTrue("the two fitted statistics coincide",
                first.test.statistic != HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL)
                        .test.statistic);
    }

    @Test
    public void testTheFittedAndersonDarlingRejectsWhatItCannotTest() {
        double[] ok = { 1.0, 2.0, 3.0, 4.0 };
        rejectsFittedAndersonDarling("a null sample", null, Lilliefors.Family.NORMAL, 100, 1L);
        rejectsFittedAndersonDarling("a null family", ok, null, 100, 1L);
        rejectsFittedAndersonDarling("two observations and two parameters", new double[] { 1.0, 2.0 },
                Lilliefors.Family.NORMAL, 100, 1L);
        rejectsFittedAndersonDarling("a NaN observation", new double[] { 1.0, 2.0, Double.NaN },
                Lilliefors.Family.NORMAL, 100, 1L);
        rejectsFittedAndersonDarling("no replications", ok, Lilliefors.Family.NORMAL, 0, 1L);
        rejectsFittedAndersonDarling("a negative observation in a positive family",
                new double[] { 1.0, 2.0, -3.0 }, Lilliefors.Family.EXPONENTIAL, 100, 1L);

        // and the null distribution on its own
        try {
            Lilliefors.barF(null, Lilliefors.Family.NORMAL, 10, 0.3, 100, 1L);
            fail("barF accepted a null statistic");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            Lilliefors.statistic(null, Lilliefors.Family.NORMAL, ok);
            fail("statistic accepted a null statistic");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsFittedAndersonDarling(String what, double[] x, Lilliefors.Family family,
            int replications, long seed) {
        try {
            HypothesisTests.andersonDarling(x, family, replications, seed);
            fail("andersonDarling accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ----------------------------------------------------- Cramer-von Mises --

    /**
     * {@code W_n^2} of a sorted sample of transformed values, written out. The
     * facade computes a p-value with it, and that costs twelve milliseconds,
     * so anything drawing a null distribution measures the statistic directly.
     */
    private static double cramerVonMisesStatistic(double[] sortedUniform) {
        int n = sortedUniform.length;
        double sum = 1.0 / (12.0 * n);
        for (int i = 0; i < n; i++) {
            double gap = sortedUniform[i] - (2 * i + 1) / (2.0 * n);
            sum += gap * gap;
        }
        return sum;
    }

    /** The value of the statistic whose upper tail is {@code alpha}. */
    private static double cramerVonMisesCritical(int n, double alpha) {
        double lo = 0.0;
        double hi = 10.0;
        for (int i = 0; i < 40; i++) {
            double mid = 0.5 * (lo + hi);
            if (CramerVonMises.barF(n, mid) > alpha) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    @Test
    public void testTheCramerVonMisesStatisticIsTheDefinition() {
        Normal standard = new Normal(0.0, 1.0);
        for (int n : new int[] { 1, 2, 7, 40 }) {
            for (long seed = 1; seed <= 25; seed++) {
                double[] x = sample(new Normal(0.5, 2.0), n, seed * 7919L + 1);
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = standard.cdf(x[i]);
                }
                Arrays.sort(u);
                TestResult result = HypothesisTests.cramerVonMises(x, standard);
                assertEquals("n=" + n + " seed=" + seed, cramerVonMisesStatistic(u), result.statistic,
                        1.0e-14);
                // and the facade reads its p-value off the null distribution
                assertEquals("n=" + n + " seed=" + seed + ": the p-value",
                        CramerVonMises.barF(n, result.statistic), result.pValue, 0.0);
            }
        }
    }

    @Test
    public void testTheCramerVonMisesNullAgreesWithThePublishedPoints() {
        // the asymptotic percentage points of W^2, which is the check that the
        // eigenvalues 1 / (pi^2 j^2) and the dropped mean of the truncated
        // series are both right
        double[] levels = { 0.10, 0.05, 0.025, 0.01 };
        double[] points = { 0.347, 0.461, 0.581, 0.743 };
        double worst = 0.0;
        for (int i = 0; i < levels.length; i++) {
            worst = Math.max(worst, Math.abs(CramerVonMises.barF(100000, points[i]) - levels[i]));
        }
        // measured: 2.0e-4, and that is the table rounded to three decimals
        assertTrue("the published points come back as " + worst + " off", worst < 5.0e-4);
    }

    @Test
    public void testTheCramerVonMisesTailIsRightWhereItDecides() {
        // Stephens' modification is a tail device: from the ten percent point
        // down it is right to a fraction of a percent at every sample size, and
        // above it it carries a small sample too far. Both halves are pinned,
        // the second one loosely, so that neither can change unnoticed
        Uniform uniform = new Uniform(0.0, 1.0);
        int reps = 100000;
        for (int n : new int[] { 5, 20, 100 }) {
            double[] drawn = new double[reps];
            for (int b = 0; b < reps; b++) {
                double[] u = sample(uniform, n, (b + 1L) * 7919L + 1);
                Arrays.sort(u);
                drawn[b] = cramerVonMisesStatistic(u);
            }
            Arrays.sort(drawn);
            for (double tail : new double[] { 0.10, 0.05, 0.025, 0.01 }) {
                double quantile = drawn[(int) ((1.0 - tail) * reps)];
                // measured over ten sample sizes: 0.0094 at worst, and that at
                // n = 3; from n = 5 up it is 0.0025
                assertEquals("n=" + n + " at the " + tail + " point", tail,
                        CramerVonMises.barF(n, quantile), 0.006);
            }
            // and the body, where the modification is only a rough guide
            double median = drawn[reps / 2];
            double atMedian = CramerVonMises.barF(n, median);
            // measured: 0.693 at n = 5, 0.554 at n = 20, 0.510 at n = 100
            assertTrue("n=" + n + ": the median came back at " + atMedian,
                    atMedian > 0.45 && atMedian < 0.75);
        }
    }

    @Test
    public void testTheCramerVonMisesRejectsAtItsLevel() {
        Uniform uniform = new Uniform(0.0, 1.0);
        for (int n : new int[] { 1, 5, 30 }) {
            double atFive = cramerVonMisesCritical(n, 0.05);
            double atOne = cramerVonMisesCritical(n, 0.01);
            int reps = 4000;
            int five = 0;
            int one = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] u = sample(uniform, n, seed * 7919L + 1);
                Arrays.sort(u);
                double statistic = cramerVonMisesStatistic(u);
                if (statistic >= atFive) {
                    five++;
                }
                if (statistic >= atOne) {
                    one++;
                }
            }
            // measured: 0.0498 .. 0.0580 and 0.0088 .. 0.0113
            assertRate("n=" + n + " at 0.05", five / (double) reps, 0.05, 0.018);
            assertRate("n=" + n + " at 0.01", one / (double) reps, 0.01, 0.008);
        }
    }

    @Test
    public void testTheCramerVonMisesSeesMoreThanTheSupremumDoes() {
        // W^2 adds up the whole discrepancy where D_n keeps only its widest
        // point, so it finds a departure that is spread out more often. It does
        // not beat A_n^2 anywhere that was tried, including the shift this test
        // uses -- that was measured rather than assumed, and the assertion says
        // only what held
        Normal standard = new Normal(0.0, 1.0);
        for (int n : new int[] { 20, 50 }) {
            double critical = cramerVonMisesCritical(n, 0.05);
            int reps = 2000;
            int kolmogorovSmirnov = 0;
            int cramerVonMises = 0;
            int andersonDarling = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(0.35, 1.0), n, seed * 7919L + 1);
                if (HypothesisTests.kolmogorovSmirnov(x, standard, Alternative.TWO_SIDED).rejectsAt(0.05)) {
                    kolmogorovSmirnov++;
                }
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = standard.cdf(x[i]);
                }
                Arrays.sort(u);
                if (cramerVonMisesStatistic(u) >= critical) {
                    cramerVonMises++;
                }
                if (HypothesisTests.andersonDarling(x, standard).rejectsAt(0.05)) {
                    andersonDarling++;
                }
            }
            // measured: 0.2675 / 0.3105 / 0.3220 at n = 20 and
            // 0.5675 / 0.6410 / 0.6630 at n = 50
            assertTrue("n=" + n + ": CvM found " + cramerVonMises + " of " + reps + " and KS "
                    + kolmogorovSmirnov, cramerVonMises > kolmogorovSmirnov);
            assertTrue("n=" + n + ": AD found " + andersonDarling + " and CvM " + cramerVonMises,
                    andersonDarling >= cramerVonMises);
        }
    }

    @Test
    public void testTheFittedCramerVonMisesHoldsItsLevelWhereTheFullySpecifiedOneCollapses() {
        for (int n : new int[] { 5, 30 }) {
            double critical = cramerVonMisesCritical(n, 0.05);
            int reps = 1200;
            int fitted = 0;
            int specified = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(3.0, 2.0), n, seed * 7919L + 1);
                if (HypothesisTests.cramerVonMises(x, Lilliefors.Family.NORMAL, 2000, 7L).test
                        .rejectsAt(0.05)) {
                    fitted++;
                }
                ParNormal fit = MLE.getNormalMLE(x);
                Normal fitDistribution = new Normal(fit.mean, fit.stdDev);
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = fitDistribution.cdf(x[i]);
                }
                Arrays.sort(u);
                if (cramerVonMisesStatistic(u) >= critical) {
                    specified++;
                }
            }
            // measured: 0.0460 .. 0.0573 fitted, and 0 of 6000 for the other
            assertRate("n=" + n + " fitted", fitted / (double) reps, 0.05, 0.018);
            assertTrue("n=" + n + ": the fully specified test rejected " + specified + " of " + reps,
                    specified == 0);
        }
    }

    @Test
    public void testTheCramerVonMisesAtOneObservationIsExact() {
        // W_1^2 = 1/12 + (u - 1/2)^2 lives on [1/12, 1/3] and its distribution
        // function is 2 sqrt(x - 1/12) exactly, which is where the modification
        // used above has nothing to say: read through the limit instead, a
        // statistic of 0.5 would earn a p-value for a value it cannot reach
        double worst = 0.0;
        for (int step = 0; step <= 500; step++) {
            double x = 1.0 / 12.0 + step * (0.25 / 500.0);
            double exact = Math.min(1.0, 2.0 * Math.sqrt(x - 1.0 / 12.0));
            worst = Math.max(worst, Math.abs(CramerVonMises.cdf(1, x) - exact));
            worst = Math.max(worst, Math.abs(CramerVonMises.barF(1, x) - (1.0 - exact)));
        }
        // measured: 1.1e-16
        assertTrue("the closed form is off by " + worst, worst < 1.0e-14);
        assertEquals("below the support", 1.0, CramerVonMises.barF(1, 0.05), 0.0);
        assertEquals("above the support", 0.0, CramerVonMises.barF(1, 0.4), 0.0);

        // and the deep tail does not cancel: near 1/3 the complement is
        // 4 (1/3 - x) / (1 + 2 sqrt(x - 1/12)), which stays positive
        double tail = CramerVonMises.barF(1, 1.0 / 3.0 - 1.0e-12);
        // measured: 2.0e-12
        assertEquals("the deep tail at n = 1", 2.0e-12, tail, 1.0e-14);
    }

    @Test
    public void testTheFittedCramerVonMisesMeasuresWhatTheFacadeMeasures() {
        for (int n : new int[] { 5, 20, 100 }) {
            for (long seed = 1; seed <= 60; seed++) {
                double[] x = sample(new Normal(1.0, 3.0), n, seed * 7919L + 1);
                ParNormal fit = MLE.getNormalMLE(x);
                assertEquals("the fitted unweighted distance",
                        HypothesisTests.cramerVonMises(x, new Normal(fit.mean, fit.stdDev)).statistic,
                        Lilliefors.statistic(Lilliefors.Statistic.CRAMER_VON_MISES,
                                Lilliefors.Family.NORMAL, x),
                        0.0);
            }
        }

        // exponentiating a sample and fitting a log-normal is fitting a normal
        // to the logarithms, and this statistic is invariant under a monotone
        // transformation too
        for (long seed = 1; seed <= 60; seed++) {
            double[] normal = sample(new Normal(0.5, 1.5), 20, seed * 7919L + 1);
            double[] exponentiated = new double[normal.length];
            for (int i = 0; i < normal.length; i++) {
                exponentiated[i] = Math.exp(normal[i]);
            }
            assertEquals("the p-value",
                    HypothesisTests.cramerVonMises(normal, Lilliefors.Family.NORMAL, 500, 5L).test.pValue,
                    HypothesisTests.cramerVonMises(exponentiated, Lilliefors.Family.LOGNORMAL, 500, 5L)
                            .test.pValue,
                    0.0);
        }
    }

    @Test
    public void testCramerVonMisesRejectsWhatItCannotTest() {
        Normal standard = new Normal(0.0, 1.0);
        try {
            HypothesisTests.cramerVonMises(null, standard);
            fail("cramerVonMises accepted a null sample");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            HypothesisTests.cramerVonMises(new double[] { 1.0, Double.NaN }, standard);
            fail("cramerVonMises accepted a NaN observation");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            HypothesisTests.cramerVonMises(new double[] { 1.0, 2.0 }, (ContinuousDistribution) null);
            fail("cramerVonMises accepted a null distribution");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            HypothesisTests.cramerVonMises(new double[] { 1.0, 2.0 }, (Lilliefors.Family) null);
            fail("cramerVonMises accepted a null family");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            CramerVonMises.barF(0, 0.5);
            fail("barF accepted a sample size of zero");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // ----------------------------------------------------- WeightedChiSquare --

    @Test
    public void testTheWeightedChiSquareAgreesWithASingleChiSquare() {
        // one weight leaves an ordinary chi-square with one degree of freedom,
        // read through a scale that may be negative, and both tails of that are
        // error functions
        for (double weight : new double[] { 1.0, 3.5, 0.25, -2.0 }) {
            for (double x : new double[] { -3.0, -0.5, 0.0, 0.5, 2.0, 9.0 }) {
                double q = x / weight;
                double lower;
                if (q < 0.0) {
                    lower = weight > 0.0 ? 0.0 : 1.0;
                } else {
                    // P[z^2 <= q] = 2 Phi(sqrt(q)) - 1, an independent route
                    double chi = 2.0 * new Normal(0.0, 1.0).cdf(Math.sqrt(q)) - 1.0;
                    lower = weight > 0.0 ? chi : 1.0 - chi;
                }
                String at = "weight=" + weight + " x=" + x;
                assertEquals(at, lower, WeightedChiSquare.cdf(new double[] { weight }, x), 1.0e-12);
                assertEquals(at + ": the tails partition", 1.0,
                        WeightedChiSquare.cdf(new double[] { weight }, x)
                                + WeightedChiSquare.barF(new double[] { weight }, x),
                        1.0e-12);
            }
        }
    }

    @Test
    public void testTheWeightedChiSquareIgnoresZeroWeights() {
        double[] plain = { 2.0, 1.5, 1.0, 0.75, 0.5 };
        double[] padded = { 2.0, 0.0, 1.5, 1.0, 0.0, 0.0, 0.75, 0.5 };
        for (double x : new double[] { 1.5, 6.0 }) {
            assertEquals("a zero weight contributes nothing at x = " + x,
                    WeightedChiSquare.cdf(plain, x), WeightedChiSquare.cdf(padded, x), 0.0);
        }
        // and a sum of nothing is zero, which every tail has to answer
        double[] allZero = { 0.0, 0.0 };
        assertEquals("P[0 <= 1]", 1.0, WeightedChiSquare.cdf(allZero, 1.0), 0.0);
        assertEquals("P[0 <= -1]", 0.0, WeightedChiSquare.cdf(allZero, -1.0), 0.0);
        assertEquals("P[0 >= -1]", 1.0, WeightedChiSquare.barF(allZero, -1.0), 0.0);
        assertEquals("P[0 >= 1]", 0.0, WeightedChiSquare.barF(allZero, 1.0), 0.0);
    }

    @Test
    public void testTheWeightedChiSquareShiftMovesTheDistribution() {
        // the shift is the one thing the Durbin-Watson form never exercises,
        // since a ratio is the x = 0 case. Every call here is a quadrature, and
        // with a shift and few weights each one is expensive -- see the class
        // javadoc -- so this walks a handful of points rather than a grid
        double[] lambda = { 2.0, 1.5, 1.25, 1.0, 0.75, 0.5, 0.25, 0.1 };
        double mean = 0.0;
        for (int i = 0; i < lambda.length; i++) {
            mean += lambda[i];
        }
        double previous = -1.0;
        for (double x : new double[] { -1.0, 0.0, 2.0, 5.0, 7.35, 10.0, 16.0, 25.0 }) {
            double lower = WeightedChiSquare.cdf(lambda, x);
            assertTrue("the distribution function fell at x = " + x, lower >= previous - 1.0e-12);
            assertEquals("the tails partition at x = " + x, 1.0,
                    lower + WeightedChiSquare.barF(lambda, x), 1.0e-12);
            previous = lower;
        }
        // every weight is positive, so the sum is too. The quadrature reaches
        // that as the difference of two numbers near one half, so it is nothing
        // to round-off rather than nothing exactly
        assertEquals("nothing is below zero", 0.0, WeightedChiSquare.cdf(lambda, 0.0), 1.0e-15);
        assertEquals("and everything is above it", 1.0, WeightedChiSquare.barF(lambda, 0.0), 1.0e-15);

        // the mean of the sum is sum lambda, and a positive weighted sum of
        // chi-squares is right-skewed, so its mean sits above its median
        double atMean = WeightedChiSquare.cdf(lambda, mean);
        assertTrue("the mean of " + mean + " came back at " + atMean, atMean > 0.5 && atMean < 0.8);
    }

    @Test
    public void testTheWeightedChiSquareRejectsWhatItCannotCompute() {
        double[] ok = { 1.0, 2.0 };
        try {
            WeightedChiSquare.cdf(null, 1.0);
            fail("cdf accepted a null weight vector");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            WeightedChiSquare.cdf(new double[0], 1.0);
            fail("cdf accepted an empty weight vector");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            WeightedChiSquare.cdf(new double[] { 1.0, Double.POSITIVE_INFINITY }, 1.0);
            fail("cdf accepted a weight that is not finite");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            WeightedChiSquare.cdf(ok, Double.NaN);
            fail("cdf accepted a threshold that is not a number");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            WeightedChiSquare.cdf(ok, 1.0, 0.0);
            fail("cdf accepted a tolerance of zero");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    // --------------------------------------------------------- Durbin-Watson --

    /** An intercept, a trend and however many smooth columns are wanted. */
    private static DMatrix durbinWatsonDesign(int n, int k, long seed) {
        Lcg lcg = new Lcg(seed);
        Normal standard = new Normal(0.0, 1.0);
        DMatrix X = new DMatrix(n, k);
        for (int i = 0; i < n; i++) {
            X.set(i, 0, 1.0);
        }
        for (int column = 1; column < k; column++) {
            for (int i = 0; i < n; i++) {
                X.set(i, column, column == 1 ? i / (double) n : standard.inverseCdf(lcg.next()));
            }
        }
        return X;
    }

    private static double[] olsResiduals(DMatrix X, double[] y) {
        DMatrix regressand = new DMatrix(y.length, 1);
        for (int i = 0; i < y.length; i++) {
            regressand.set(i, 0, y[i]);
        }
        LSSummary fit = OLS.estimate(0.05, X, regressand);
        double[] e = new double[y.length];
        for (int i = 0; i < e.length; i++) {
            e[i] = fit.getResiduals().get(i);
        }
        return e;
    }

    /** Errors from a first-order autoregression with the given correlation. */
    private static double[] autocorrelated(int n, double rho, long seed) {
        double[] innovation = sample(new Normal(0.0, 1.0), n, seed);
        double[] e = new double[n];
        e[0] = innovation[0] / Math.sqrt(1.0 - rho * rho);
        for (int i = 1; i < n; i++) {
            e[i] = rho * e[i - 1] + innovation[i];
        }
        return e;
    }

    @Test
    public void testTheDurbinWatsonStatisticIsTheDefinition() {
        for (int n : new int[] { 8, 25, 60 }) {
            DMatrix X = durbinWatsonDesign(n, 2, 5L);
            for (long seed = 1; seed <= 40; seed++) {
                double[] e = sample(new Normal(0.0, 3.0), n, seed * 7919L + 1);
                double numerator = 0.0;
                double denominator = e[0] * e[0];
                for (int i = 1; i < n; i++) {
                    double step = e[i] - e[i - 1];
                    numerator += step * step;
                    denominator += e[i] * e[i];
                }
                assertEquals("n=" + n + " seed=" + seed, numerator / denominator,
                        HypothesisTests.durbinWatson(e, X, Alternative.TWO_SIDED).statistic, 0.0);
            }
        }
    }

    @Test
    public void testTheDurbinWatsonSpectrumOfAnInterceptOnlyDesignHasAClosedForm() {
        // the column space of a design that is only an intercept is spanned by
        // the eigenvector of the zero eigenvalue of A, so what is left is the
        // rest of the spectrum of A itself, which is known in closed form. This
        // checks the decomposition, the assembly and the ordering at once
        double worst = 0.0;
        for (int n : new int[] { 6, 20, 60 }) {
            double[] intercept = new double[n];
            Arrays.fill(intercept, 1.0);
            double[] nu = DurbinWatson.nullEigenvalues(intercept, n, 1);
            assertEquals("n=" + n + ": there are n - 1 of them", n - 1, nu.length);
            for (int i = 0; i < nu.length; i++) {
                // they come back descending, so nu[i] belongs to j = n - 1 - i
                double closedForm = 2.0 * (1.0 - Math.cos(Math.PI * (n - 1 - i) / (double) n));
                worst = Math.max(worst, Math.abs(nu[i] - closedForm));
            }
            // and their sum is the trace of A, which the intercept takes nothing
            // from because its eigenvalue is zero
            double sum = 0.0;
            for (int i = 0; i < nu.length; i++) {
                sum += nu[i];
            }
            assertEquals("n=" + n + ": the trace", 2.0 * n - 2.0, sum, 1.0e-12);
        }
        // measured: 4.9e-15
        assertTrue("the spectrum is off by " + worst, worst < 1.0e-13);
    }

    @Test
    public void testTheDurbinWatsonTailAgreesWithTheArcsineLaw() {
        // with two eigenvalues d is nu2 + (nu1 - nu2) w with w ~ Beta(1/2, 1/2),
        // so the distribution function is an arcsine and Imhof can be checked
        // against algebra rather than against another approximation
        double[] nu = { 3.0, 1.0 };
        for (double p : new double[] { 0.5, 0.1, 1.0e-2, 1.0e-3, 1.0e-4 }) {
            double share = Math.sin(0.5 * Math.PI * p);
            share = share * share;
            double exact = 2.0 / Math.PI * Math.asin(Math.sqrt(share));
            double got = DurbinWatson.cdf(nu, 1.0 + 2.0 * share);
            // measured: 4.9e-14 at 1e-2, 7.0e-10 at 1e-4
            assertEquals("at p = " + p, exact, got, 1.0e-8 * exact);
        }

        // and this is where it stops: the tails are 1/2 minus the same integral,
        // so a p-value far below the level anyone decides at is only an order of
        // magnitude. The floor is pinned so that it cannot quietly rise
        double share = Math.sin(0.5 * Math.PI * 1.0e-6);
        share = share * share;
        double deep = DurbinWatson.cdf(nu, 1.0 + 2.0 * share);
        double exactDeep = 2.0 / Math.PI * Math.asin(Math.sqrt(share));
        // measured: 8.4e-6 relative
        assertEquals("at p = 1e-6", exactDeep, deep, 1.0e-4 * exactDeep);
    }

    @Test
    public void testTheDurbinWatsonTailAgreesWithASimulationOfTheRatio() {
        // Imhof against the definition it inverts: draw z, form the ratio
        Normal standard = new Normal(0.0, 1.0);
        int reps = 100000;
        double worst = 0.0;
        for (int[] shape : new int[][] { { 15, 1 }, { 20, 3 } }) {
            DMatrix X = durbinWatsonDesign(shape[0], shape[1], 5L);
            double[] nu = DurbinWatson.nullEigenvalues(X.getArrayUnsafe(), shape[0], shape[1]);
            double[] drawn = new double[reps];
            for (int b = 0; b < reps; b++) {
                double[] z = sample(standard, nu.length, (b + 1L) * 7919L + 1);
                double numerator = 0.0;
                double denominator = 0.0;
                for (int i = 0; i < nu.length; i++) {
                    double square = z[i] * z[i];
                    numerator += nu[i] * square;
                    denominator += square;
                }
                drawn[b] = numerator / denominator;
            }
            Arrays.sort(drawn);
            for (double tail : new double[] { 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99 }) {
                worst = Math.max(worst,
                        Math.abs(DurbinWatson.cdf(nu, drawn[(int) (tail * reps)]) - tail));
            }
        }
        // measured: 0.00154
        assertTrue("the drawn distribution differs by " + worst, worst < 0.006);
    }

    @Test
    public void testTheDurbinWatsonPValuesAreUniformUnderTheNullHypothesis() {
        // the whole chain at once: a fit, its residuals, the decomposition of
        // its design and the inversion. A p-value that is uniform is right for
        // reasons no single value could show
        Normal standard = new Normal(0.0, 1.0);
        Uniform uniform = new Uniform(0.0, 1.0);
        int reps = 2000;
        for (int[] shape : new int[][] { { 15, 1 }, { 30, 3 } }) {
            int n = shape[0];
            DMatrix X = durbinWatsonDesign(n, shape[1], 5L);
            double[] nu = DurbinWatson.nullEigenvalues(X.getArrayUnsafe(), n, shape[1]);
            double[] p = new double[reps];
            int atFive = 0;
            int atOne = 0;
            for (int b = 0; b < reps; b++) {
                double[] e = olsResiduals(X, sample(standard, n, (b + 1L) * 7919L + 1));
                double numerator = 0.0;
                double denominator = e[0] * e[0];
                for (int i = 1; i < n; i++) {
                    double step = e[i] - e[i - 1];
                    numerator += step * step;
                    denominator += e[i] * e[i];
                }
                p[b] = DurbinWatson.cdf(nu, numerator / denominator);
                if (p[b] <= 0.05) {
                    atFive++;
                }
                if (p[b] <= 0.01) {
                    atOne++;
                }
            }
            String at = "n=" + n + " k=" + shape[1];
            // measured over four designs: 0.0470 .. 0.0545 and 0.0085 .. 0.0140
            assertRate(at + " at 0.05", atFive / (double) reps, 0.05, 0.018);
            assertRate(at + " at 0.01", atOne / (double) reps, 0.01, 0.008);
            // measured: 0.0347 at worst, against a 0.001 critical value of 0.043
            assertTrue(at + ": the p-values are not uniform",
                    HypothesisTests.kolmogorovSmirnov(p, uniform, Alternative.TWO_SIDED).statistic < 0.06);
        }
    }

    @Test
    public void testTheDurbinWatsonAlternativesAreAboutTheAutocorrelation() {
        // GREATER asks whether the errors are positively autocorrelated, which
        // is a SMALL statistic. Getting this backwards is a silent error, so it
        // is pinned here against errors that really are autocorrelated
        DMatrix X = durbinWatsonDesign(24, 3, 5L);
        int reps = 400;
        int[] greater = new int[3];
        int[] less = new int[3];
        double[] meanStatistic = new double[3];
        double[] rhos = { 0.0, 0.5, -0.5 };
        for (int r = 0; r < rhos.length; r++) {
            for (long seed = 1; seed <= reps; seed++) {
                double[] e = autocorrelated(24, rhos[r], seed * 7919L + 1);
                double[] y = new double[24];
                for (int i = 0; i < 24; i++) {
                    y[i] = 1.0 + 2.0 * X.get(i, 1) + e[i];
                }
                double[] residuals = olsResiduals(X, y);
                TestResult positive = HypothesisTests.durbinWatson(residuals, X, Alternative.GREATER);
                TestResult negative = HypothesisTests.durbinWatson(residuals, X, Alternative.LESS);
                meanStatistic[r] += positive.statistic;
                if (positive.rejectsAt(0.05)) {
                    greater[r]++;
                }
                if (negative.rejectsAt(0.05)) {
                    less[r]++;
                }
            }
            meanStatistic[r] /= reps;
        }
        // measured: mean d 2.0932 / 1.3207 / 2.8542, GREATER rejects
        // 0.0625 / 0.6375 / 0.0000 and LESS 0.0325 / 0.0000 / 0.6600
        assertRate("with no autocorrelation", greater[0] / (double) reps, 0.05, 0.035);
        assertTrue("positive autocorrelation must push the statistic down, mean was " + meanStatistic[1],
                meanStatistic[1] < 1.6);
        assertTrue("and must be found by GREATER", greater[1] > 0.4 * reps);
        assertTrue("and never by LESS", less[1] == 0);
        assertTrue("negative autocorrelation must push it up, mean was " + meanStatistic[2],
                meanStatistic[2] > 2.5);
        assertTrue("and must be found by LESS", less[2] > 0.4 * reps);
        assertTrue("and never by GREATER", greater[2] == 0);
    }

    @Test
    public void testTheDurbinWatsonTailsPartitionTheDistribution() {
        DMatrix X = durbinWatsonDesign(24, 3, 5L);
        double[] nu = DurbinWatson.nullEigenvalues(X.getArrayUnsafe(), 24, 3);
        double[] e = sample(new Normal(0.0, 1.0), 24, 17L);
        double previous = -1.0;
        for (double d : new double[] { 0.4, 0.8, 1.4, 2.0, 2.6, 3.2, 3.6 }) {
            double lower = DurbinWatson.cdf(nu, d);
            double upper = DurbinWatson.barF(nu, d);
            assertEquals("the two tails at d = " + d, 1.0, lower + upper, 1.0e-12);
            assertTrue("the lower tail is not increasing at d = " + d, lower > previous);
            previous = lower;
        }

        // and the facade divides them the way it says it does
        TestResult positive = HypothesisTests.durbinWatson(e, X, Alternative.GREATER);
        TestResult negative = HypothesisTests.durbinWatson(e, X, Alternative.LESS);
        TestResult both = HypothesisTests.durbinWatson(e, X, Alternative.TWO_SIDED);
        assertEquals("the one-sided p-values partition", 1.0, positive.pValue + negative.pValue, 1.0e-12);
        assertEquals("the two-sided one is twice the smaller",
                Math.min(1.0, 2.0 * Math.min(positive.pValue, negative.pValue)), both.pValue, 0.0);
        assertEquals("and they all report the same statistic", positive.statistic, both.statistic, 0.0);
    }

    @Test
    public void testTheDurbinWatsonStatisticIsInvariantUnderScaleAndSign() {
        DMatrix X = durbinWatsonDesign(30, 2, 5L);
        for (long seed = 1; seed <= 30; seed++) {
            double[] e = sample(new Normal(0.0, 1.0), 30, seed * 7919L + 1);
            double reference = HypothesisTests.durbinWatson(e, X, Alternative.GREATER).statistic;
            double[] scaled = new double[e.length];
            double[] negated = new double[e.length];
            for (int i = 0; i < e.length; i++) {
                // a power of two changes the exponent and nothing else, so the
                // ratio has to come back bit for bit
                scaled[i] = e[i] * 1024.0;
                negated[i] = -e[i];
            }
            assertEquals("scaled, seed=" + seed, reference,
                    HypothesisTests.durbinWatson(scaled, X, Alternative.GREATER).statistic, 0.0);
            assertEquals("negated, seed=" + seed, reference,
                    HypothesisTests.durbinWatson(negated, X, Alternative.GREATER).statistic, 0.0);
        }
    }

    @Test
    public void testTheExactDurbinWatsonTailLiesBetweenTheInterlacingBounds() {
        // Poincare separation: whatever the k columns of X are, the eigenvalues
        // of the restriction interlace those of A itself, so the bottom m and
        // the top m of 2 (1 - cos(pi j / n)) bound the distribution. That is the
        // classical pair of bounds, and it holds without knowing the design
        for (int n : new int[] { 20, 40 }) {
            for (int k : new int[] { 1, 2, 4 }) {
                DMatrix X = durbinWatsonDesign(n, k, 5L);
                double[] exact = DurbinWatson.nullEigenvalues(X.getArrayUnsafe(), n, k);
                int m = n - k;
                double[] all = new double[n];
                for (int index = 0; index < n; index++) {
                    all[index] = 2.0 * (1.0 - Math.cos(Math.PI * index / (double) n));
                }
                Arrays.sort(all);
                double[] bottom = new double[m];
                double[] top = new double[m];
                for (int i = 0; i < m; i++) {
                    bottom[i] = all[i];
                    top[i] = all[i + k];
                }
                // the eigenvalues themselves interlace, ascending
                for (int i = 0; i < m; i++) {
                    double ascending = exact[m - 1 - i];
                    assertTrue("n=" + n + " k=" + k + " i=" + i + ": " + ascending + " left ["
                            + bottom[i] + ", " + top[i] + "]",
                            ascending >= bottom[i] - 1.0e-9 && ascending <= top[i] + 1.0e-9);
                }
                for (double d : new double[] { 1.0, 1.5, 2.0, 2.5 }) {
                    double atMost = DurbinWatson.cdf(top, d);
                    double atLeast = DurbinWatson.cdf(bottom, d);
                    double got = DurbinWatson.cdf(exact, d);
                    assertTrue("n=" + n + " k=" + k + " d=" + d + ": " + got + " left [" + atMost + ", "
                            + atLeast + "]", got >= atMost - 1.0e-12 && got <= atLeast + 1.0e-12);
                }
            }
        }
    }

    @Test
    public void testTheTwoDurbinWatsonEntryPointsAgree() {
        int n = 24;
        DMatrix X = durbinWatsonDesign(n, 3, 5L);
        DMatrix y = new DMatrix(n, 1);
        double[] noise = sample(new Normal(0.0, 1.0), n, 13L);
        for (int i = 0; i < n; i++) {
            y.set(i, 0, 1.0 + 2.0 * X.get(i, 1) + noise[i]);
        }
        LSSummary fit = OLS.estimate(0.05, X, y);
        double[] residuals = new double[n];
        for (int i = 0; i < n; i++) {
            residuals[i] = fit.getResiduals().get(i);
        }
        for (Alternative alternative : Alternative.values()) {
            TestResult viaFit = HypothesisTests.durbinWatson(fit, alternative);
            TestResult viaArrays = HypothesisTests.durbinWatson(residuals, X, alternative);
            assertEquals(alternative + ": statistic", viaArrays.statistic, viaFit.statistic, 0.0);
            assertEquals(alternative + ": p-value", viaArrays.pValue, viaFit.pValue, 0.0);
            assertTrue(alternative + ": the result names the test", viaFit.test.contains("Durbin-Watson"));
        }
    }

    @Test
    public void testDurbinWatsonRejectsWhatItCannotTest() {
        DMatrix X = durbinWatsonDesign(10, 2, 5L);
        double[] ok = sample(new Normal(0.0, 1.0), 10, 3L);
        rejectsDurbinWatson("a null residual vector", null, X, Alternative.GREATER);
        rejectsDurbinWatson("one residual", new double[] { 1.0 }, durbinWatsonDesign(1, 1, 5L),
                Alternative.GREATER);
        rejectsDurbinWatson("a NaN residual", withNaN(ok), X, Alternative.GREATER);
        rejectsDurbinWatson("a null design matrix", ok, null, Alternative.GREATER);
        rejectsDurbinWatson("a null alternative", ok, X, null);
        rejectsDurbinWatson("a design matrix of the wrong height", ok, durbinWatsonDesign(11, 2, 5L),
                Alternative.GREATER);
        rejectsDurbinWatson("residuals that are all zero", new double[10], X, Alternative.GREATER);

        // a design matrix whose columns repeat has no residual space of the
        // size the test assumes, and there is no honest answer to give
        DMatrix repeated = new DMatrix(10, 2);
        for (int i = 0; i < 10; i++) {
            repeated.set(i, 0, 1.0);
            repeated.set(i, 1, 2.0);
        }
        rejectsDurbinWatson("a rank deficient design", ok, repeated, Alternative.GREATER);

        // the summary form
        try {
            HypothesisTests.durbinWatson((LSSummary) null, Alternative.GREATER);
            fail("durbinWatson accepted a null fit");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }

        DMatrix y = new DMatrix(10, 1);
        for (int i = 0; i < 10; i++) {
            y.set(i, 0, ok[i]);
        }
        double[] weights = new double[10];
        Arrays.fill(weights, 1.0);
        weights[3] = 4.0;
        try {
            HypothesisTests.durbinWatson(Wls.estimate(0.05, X, y, weights), Alternative.GREATER);
            fail("durbinWatson accepted a weighted fit");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message does not say it was the weights: " + expected.getMessage(),
                    expected.getMessage() != null && expected.getMessage().contains("weighted"));
        }

        LSSummary cleared = OLS.estimate(0.05, X, y);
        cleared.clearTemporaries();
        try {
            HypothesisTests.durbinWatson(cleared, Alternative.GREATER);
            fail("durbinWatson accepted a fit that had released its residuals");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message does not name clearTemporaries: " + expected.getMessage(),
                    expected.getMessage() != null && expected.getMessage().contains("clearTemporaries"));
        }

        // and the null distribution on its own
        try {
            DurbinWatson.cdf(new double[0], 2.0);
            fail("cdf accepted an empty spectrum");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            DurbinWatson.nullEigenvalues(new double[20], 10, 10);
            fail("nullEigenvalues accepted as many columns as rows");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static double[] withNaN(double[] x) {
        double[] broken = x.clone();
        broken[2] = Double.NaN;
        return broken;
    }

    private static void rejectsDurbinWatson(String what, double[] residuals, DMatrix X,
            Alternative alternative) {
        try {
            HypothesisTests.durbinWatson(residuals, X, alternative);
            fail("durbinWatson accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
