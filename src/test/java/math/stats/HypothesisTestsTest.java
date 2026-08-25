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
import math.list.DoubleArrayList;
import math.stats.gof.KolmogorovSmirnov;
import math.stats.gof.KolmogorovSmirnovTwoSample;
import math.stats.gof.Lilliefors;
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
            LillieforsResult fitted = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL, 200, 3L);
            String rootFitted = fitted.toString();
            Locale.setDefault(Locale.GERMANY);
            assertEquals("LillieforsResult.toString", rootFitted, fitted.toString());
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
                        assertEquals("the degrees of freedom did not swap", ab.numeratorDf, ba.denominatorDf);
                        assertEquals("the degrees of freedom did not swap", ab.denominatorDf, ba.numeratorDf);
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
        LillieforsResult first = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL);
        LillieforsResult second = HypothesisTests.lilliefors(x, Lilliefors.Family.NORMAL);
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
                assertEquals("the fitted distance", HypothesisTests.kolmogorovSmirnov(x,
                        new Normal(fit.mean, fit.stdDev), Alternative.TWO_SIDED).statistic,
                        Lilliefors.statistic(Lilliefors.Family.NORMAL, x), 0.0);
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
}
