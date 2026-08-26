package math.stats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.distribution.Normal;
import math.distribution.StudentT;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;

/**
 * Each of the four k-sample tests is checked against something that computes
 * the same number by another route, and then against the assumption it rests
 * on -- which is where the choice between them is actually made, and which is
 * measured here rather than asserted.
 */
public final class KSampleTestsTest {

    /**
     * A SplitMix64 finalizer over the in-test generator's seed. Streams started
     * at neighbouring seeds are not independent enough to measure a rejection
     * rate with, which the Kruskal-Wallis calibration found out by measuring
     * one.
     */
    private static long mix(long z) {
        z = z + 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    private static final Normal NORMAL = new Normal(0.0, 1.0);
    private static final StudentT HEAVY = new StudentT(3.0);

    /**
     * {@code k} deterministic groups, normal or {@code t(3)}, each with its own
     * mean and spread. One generator feeds all of them, so the groups are
     * independent of one another rather than merely differently seeded.
     */
    private static double[][] groups(long seed, int[] sizes, double[] means, double[] spreads, boolean heavy) {
        long state = mix(seed);
        double[][] g = new double[sizes.length][];
        for (int i = 0; i < sizes.length; i++) {
            g[i] = new double[sizes[i]];
            for (int j = 0; j < sizes[i]; j++) {
                state = state * 6364136223846793005L + 1442695040888963407L;
                double u = ((state >>> 11) + 0.5) * 0x1.0p-53;
                g[i][j] = means[i] + spreads[i] * (heavy ? HEAVY.inverseCdf(u) : NORMAL.inverseCdf(u));
            }
        }
        return g;
    }

    /** {@code k} equal group sizes. */
    private static int[] sizes(int k, int n) {
        int[] a = new int[k];
        Arrays.fill(a, n);
        return a;
    }

    /** {@code k} copies of one number, for the means or the spreads. */
    private static double[] each(int k, double value) {
        double[] a = new double[k];
        Arrays.fill(a, value);
        return a;
    }

    /** Standard normal groups of equal size, all with the same mean. */
    private static double[][] nullGroups(long seed, int k, int n, boolean heavy) {
        return groups(seed, sizes(k, n), each(k, 0.0), each(k, 1.0), heavy);
    }

    // -------------------------------------------------- against another route --

    @Test
    public void testTwoGroupsMakeTheAnovaThePooledTTestSquared() {
        // the analysis of variance of two groups is the pooled t-test, and the
        // two live in different parts of this class: one sums squares about a
        // grand mean, the other forms a difference over a standard error. The
        // relation is exact in arithmetic, so what is measured here is only how
        // much the two spellings differ
        double worstF = 0.0;
        double worstP = 0.0;
        for (long seed = 1; seed <= 800; seed++) {
            int nx = 2 + (int) (seed % 30);
            int ny = 2 + (int) ((seed / 7) % 30);
            double[][] two = groups(seed, new int[] { nx, ny }, new double[] { 0.0, 0.3 }, each(2, 1.0), false);
            FTestResult f = HypothesisTests.oneWayAnova(two);
            TTestResult t = HypothesisTests.tTwoSamplePooled(two[0], two[1], Alternative.TWO_SIDED, 0.95);

            String at = "seed=" + seed + " sizes=" + nx + "," + ny;
            double squared = t.test.statistic * t.test.statistic;
            worstF = Math.max(worstF, Math.abs(f.test.statistic - squared) / squared);
            worstP = Math.max(worstP, Math.abs(f.test.pValue - t.test.pValue));
            assertEquals(at + ": numerator df", 1.0, f.numeratorDf, 0.0);
            assertEquals(at + ": denominator df", t.test.degreesOfFreedom, f.denominatorDf, 0.0);
        }
        // measured over 3000 pairs: 9.8e-16 relative on the statistic and
        // 6.5e-12 absolute on the p-value, the latter being two entirely
        // different tail integrals -- a Student t against an incomplete beta
        assertTrue("the statistics differ by " + worstF + " relative", worstF < 1.0e-13);
        assertTrue("the p-values differ by " + worstP, worstP < 1.0e-10);
    }

    @Test
    public void testTwoGroupsMakeWelchsAnovaTheWelchTTestSquared() {
        // and the same for the unequal-variance pair, which also pins the
        // fractional degrees of freedom: welchAnova reaches them through
        // 1 / lambda and tTwoSample through the Welch-Satterthwaite expression,
        // and the two share no code
        double worstF = 0.0;
        double worstDf = 0.0;
        double worstP = 0.0;
        for (long seed = 1; seed <= 800; seed++) {
            int nx = 2 + (int) (seed % 30);
            int ny = 2 + (int) ((seed / 7) % 30);
            double sx = Math.pow(10.0, (mix(seed) % 5) - 2);
            double sy = Math.pow(10.0, (mix(seed + 1) % 5) - 2);
            double[][] two = groups(seed, new int[] { nx, ny }, new double[] { 0.0, 0.3 },
                    new double[] { sx, sy }, false);
            FTestResult f = HypothesisTests.welchAnova(two);
            TTestResult t = HypothesisTests.tTwoSample(two[0], two[1], Alternative.TWO_SIDED, 0.95);

            double squared = t.test.statistic * t.test.statistic;
            worstF = Math.max(worstF, Math.abs(f.test.statistic - squared) / squared);
            worstDf = Math.max(worstDf,
                    Math.abs(f.denominatorDf - t.test.degreesOfFreedom) / t.test.degreesOfFreedom);
            worstP = Math.max(worstP, Math.abs(f.test.pValue - t.test.pValue));
            assertEquals("seed=" + seed + ": numerator df", 1.0, f.numeratorDf, 0.0);
        }
        // measured over 3000 pairs: 7.1e-16, 7.0e-16 and 9.1e-13
        assertTrue("the statistics differ by " + worstF + " relative", worstF < 1.0e-13);
        assertTrue("the degrees of freedom differ by " + worstDf + " relative", worstDf < 1.0e-13);
        assertTrue("the p-values differ by " + worstP, worstP < 1.0e-10);
    }

    @Test
    public void testTheAnovaIsTheRegressionOnGroupIndicators() {
        // the k-group cross-check, and the one that owes nothing to the
        // two-sample tests: fitting the response on an intercept and k-1
        // indicator columns gives an R^2 whose F is the analysis of variance,
        // by way of a least squares solve in math.linalg that shares no line
        // with a sum-of-squares decomposition
        double worstF = 0.0;
        for (long seed = 1; seed <= 400; seed++) {
            int k = 2 + (int) (seed % 4);
            int[] sizes = new int[k];
            double[] means = new double[k];
            int total = 0;
            for (int i = 0; i < k; i++) {
                sizes[i] = 2 + (int) ((mix(seed * 31L + i) >>> 3) % 12);
                means[i] = 0.2 * i;
                total += sizes[i];
            }
            double[][] g = groups(seed, sizes, means, each(k, 1.0), false);
            FTestResult f = HypothesisTests.oneWayAnova(g);

            DMatrix design = new DMatrix(total, k);
            DMatrix response = new DMatrix(total, 1);
            int at = 0;
            for (int i = 0; i < k; i++) {
                for (int r = 0; r < sizes[i]; r++) {
                    design.set(at, 0, 1.0);
                    for (int c = 1; c < k; c++) {
                        design.set(at, c, (c == i) ? 1.0 : 0.0);
                    }
                    response.set(at, 0, g[i][r]);
                    at++;
                }
            }
            LSSummary fit = OLS.estimate(0.05, design, response);
            double r2 = fit.getRSquared();
            double regression = (r2 / (k - 1.0)) / ((1.0 - r2) / (total - (double) k));

            assertEquals("seed=" + seed + ": denominator df", total - (double) k, f.denominatorDf, 0.0);
            assertEquals("seed=" + seed + ": residual df", total - k, fit.getDegreesOfFreedom());
            worstF = Math.max(worstF, Math.abs(f.test.statistic - regression) / regression);
        }
        // measured over 1500 fits: 1.1e-13 relative, which is a QR solve
        // against a pair of sums and not a rounding difference
        assertTrue("the two differ by " + worstF + " relative", worstF < 1.0e-11);
    }

    @Test
    public void testLeveneIsTheAnovaOfTheDistancesFromTheGroupCenter() {
        // both tests of spread are defined as an analysis of variance of the
        // distances, and are implemented as one. What this pins is that the
        // centers are the ones claimed -- the mean for Levene and the median
        // for Brown-Forsythe -- since the distances are formed here and the
        // analysis is then asked for through the public method
        for (long seed = 1; seed <= 600; seed++) {
            int k = 2 + (int) (seed % 3);
            int[] sizes = new int[k];
            for (int i = 0; i < k; i++) {
                sizes[i] = 2 + (int) ((mix(seed * 31L + i) >>> 3) % 15);
            }
            double[][] g = groups(seed, sizes, new double[k], each(k, 1.0), false);
            String at = "seed=" + seed;
            // the same arithmetic on the same numbers, so bit for bit
            assertEquals(at + ": Levene", HypothesisTests.oneWayAnova(distances(g, false)).test.statistic,
                    HypothesisTests.levene(g).test.statistic, 0.0);
            assertEquals(at + ": Brown-Forsythe", HypothesisTests.oneWayAnova(distances(g, true)).test.statistic,
                    HypothesisTests.brownForsythe(g).test.statistic, 0.0);
        }
    }

    /** {@code |x - center|} per group, about the median or about the mean. */
    private static double[][] distances(double[][] g, boolean aboutTheMedian) {
        double[][] z = new double[g.length][];
        for (int i = 0; i < g.length; i++) {
            double center = aboutTheMedian ? median(g[i]) : mean(g[i]);
            z[i] = new double[g[i].length];
            for (int j = 0; j < g[i].length; j++) {
                z[i][j] = Math.abs(g[i][j] - center);
            }
        }
        return z;
    }

    private static double mean(double[] x) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        return sum / x.length;
    }

    /**
     * The median, spelled as the implementation spells it. The midpoint form
     * and {@code (a + b) / 2} differ by one unit in the last place, and that is
     * enough to move a statistic whose denominator is nearly zero by 1e32 --
     * measured, once in two thousand samples.
     */
    private static double median(double[] x) {
        double[] sorted = x.clone();
        Arrays.sort(sorted);
        int n = sorted.length;
        if (n % 2 == 1) {
            return sorted[n / 2];
        }
        return sorted[n / 2 - 1] + 0.5 * (sorted[n / 2] - sorted[n / 2 - 1]);
    }

    // ------------------------------------------------------- what cannot move --

    @Test
    public void testAPowerOfTwoOnTheDataChangesNothing() {
        // every statistic here is a ratio in which a common factor cancels, and
        // the implementations scale by a power of two before squaring anything
        // for exactly that reason. A power of two is exact, so this is bit for
        // bit and not nearly so
        for (long seed = 1; seed <= 200; seed++) {
            int k = 2 + (int) (seed % 3);
            double[][] g = groups(seed, sizes(k, 4 + (int) (seed % 9)), new double[k], each(k, 1.0), false);
            for (int exponent : new int[] { -900, -400, -60, -1, 1, 60, 400, 900 }) {
                double factor = Math.scalb(1.0, exponent);
                double[][] scaled = new double[k][];
                for (int i = 0; i < k; i++) {
                    scaled[i] = new double[g[i].length];
                    for (int j = 0; j < g[i].length; j++) {
                        scaled[i][j] = g[i][j] * factor;
                    }
                }
                String at = "seed=" + seed + " 2^" + exponent;
                assertEquals(at + ": one-way ANOVA", HypothesisTests.oneWayAnova(g).test.statistic,
                        HypothesisTests.oneWayAnova(scaled).test.statistic, 0.0);
                assertEquals(at + ": Welch ANOVA", HypothesisTests.welchAnova(g).test.statistic,
                        HypothesisTests.welchAnova(scaled).test.statistic, 0.0);
                assertEquals(at + ": Bartlett", HypothesisTests.bartlett(g).statistic,
                        HypothesisTests.bartlett(scaled).statistic, 0.0);
                assertEquals(at + ": Levene", HypothesisTests.levene(g).test.statistic,
                        HypothesisTests.levene(scaled).test.statistic, 0.0);
                assertEquals(at + ": Brown-Forsythe", HypothesisTests.brownForsythe(g).test.statistic,
                        HypothesisTests.brownForsythe(scaled).test.statistic, 0.0);
            }
        }
    }

    @Test
    public void testTheStatisticsSurviveDataFarFromOne() {
        // the sums square what meanAndDeviation goes to some trouble not to
        // square. Without the rescale the direct route returns NaN from 1e155
        // upwards and again below 1e-200, which is what this covers
        double[][] base = nullGroups(7L, 3, 12, false);
        double anova = HypothesisTests.oneWayAnova(base).test.statistic;
        double welch = HypothesisTests.welchAnova(base).test.statistic;
        double bartlett = HypothesisTests.bartlett(base).statistic;
        double forsythe = HypothesisTests.brownForsythe(base).test.statistic;
        for (double factor : new double[] { 1.0e150, 1.0e155, 1.0e200, 1.0e300, 1.0e-150, 1.0e-200, 1.0e-300 }) {
            double[][] g = new double[3][];
            for (int i = 0; i < 3; i++) {
                g[i] = new double[base[i].length];
                for (int j = 0; j < base[i].length; j++) {
                    g[i][j] = base[i][j] * factor;
                }
            }
            String at = "x " + factor;
            // not a power of two, so the data itself is rounded and only the
            // scaling inside is exact
            assertEquals(at + ": one-way ANOVA", anova, HypothesisTests.oneWayAnova(g).test.statistic, 1.0e-12);
            assertEquals(at + ": Welch ANOVA", welch, HypothesisTests.welchAnova(g).test.statistic, 1.0e-12);
            assertEquals(at + ": Bartlett", bartlett, HypothesisTests.bartlett(g).statistic, 1.0e-12);
            assertEquals(at + ": Brown-Forsythe", forsythe, HypothesisTests.brownForsythe(g).test.statistic,
                    1.0e-12);
        }
    }

    @Test
    public void testAShiftLeavesTheTestsOfSpreadWhereTheyWere() {
        // the three tests of spread cannot see a location shift in arithmetic.
        // In floating point they can, because adding an offset and taking it
        // away again is not exact, and this is what that costs
        double worst = 0.0;
        for (long seed = 1; seed <= 300; seed++) {
            int k = 2 + (int) (seed % 3);
            double[][] g = groups(seed, sizes(k, 4 + (int) (seed % 9)), new double[k], each(k, 1.0), false);
            double bartlett = HypothesisTests.bartlett(g).statistic;
            double forsythe = HypothesisTests.brownForsythe(g).test.statistic;
            for (double by : new double[] { 1.0, 1000.0, 1.0e6 }) {
                double[][] shifted = new double[k][];
                for (int i = 0; i < k; i++) {
                    shifted[i] = new double[g[i].length];
                    for (int j = 0; j < g[i].length; j++) {
                        shifted[i][j] = g[i][j] + by;
                    }
                }
                worst = Math.max(worst, relative(bartlett, HypothesisTests.bartlett(shifted).statistic));
                worst = Math.max(worst,
                        relative(forsythe, HypothesisTests.brownForsythe(shifted).test.statistic));
            }
        }
        // measured: 2.8e-10 relative at an offset of 1e6 on unit data, which is
        // the six digits the offset costs and not a defect of the test
        assertTrue("a shift moved a test of spread by " + worst + " relative", worst < 1.0e-8);
    }

    private static double relative(double a, double b) {
        return Math.abs(a - b) / Math.max(1.0, Math.abs(a));
    }

    // ------------------------------------------------------------- the levels --

    private static final int ANOVA = 0;
    private static final int WELCH = 1;
    private static final int BARTLETT = 2;
    private static final int LEVENE = 3;
    private static final int FORSYTHE = 4;
    private static final int RANK = 5;

    /**
     * How often each of the six tests rejects at five percent, over
     * {@code reps} samples of one shape.
     * <p>
     * All six read the same sample, which is what keeps this file quick: the
     * generation costs more than the tests do, and asking one question per
     * pass paid for it six times over. It also makes the comparisons below
     * paired rather than independent, so a difference between two rates is
     * measured more sharply than either rate is.
     */
    private static double[] rates(int reps, int[] sizes, double[] means, double[] spreads, boolean heavy) {
        int[] hit = new int[6];
        for (long seed = 1; seed <= reps; seed++) {
            double[][] g = groups(seed * 131L, sizes, means, spreads, heavy);
            double[] p = { HypothesisTests.oneWayAnova(g).test.pValue,
                    HypothesisTests.welchAnova(g).test.pValue, HypothesisTests.bartlett(g).pValue,
                    HypothesisTests.levene(g).test.pValue, HypothesisTests.brownForsythe(g).test.pValue,
                    HypothesisTests.kruskalWallis(g).pValue };
            for (int i = 0; i < 6; i++) {
                if (p[i] <= 0.05) {
                    hit[i]++;
                }
            }
        }
        double[] rate = new double[6];
        for (int i = 0; i < 6; i++) {
            rate[i] = hit[i] / (double) reps;
        }
        return rate;
    }

    @Test
    public void testTheTestsHoldTheirLevelWhereTheirAssumptionsHold() {
        // normal groups, equal means, equal spreads: every null hypothesis here
        // is true and every assumption is met, so all five should spend about
        // the five percent they are asked for. Measured over 20000 samples;
        // three standard errors of a rate of 0.05 at 10000 is 0.0065
        double[][] measured = { rates(10000, sizes(3, 10), each(3, 0.0), each(3, 1.0), false),
                rates(10000, sizes(4, 25), each(4, 0.0), each(4, 1.0), false) };
        for (int i = 0; i < 2; i++) {
            String at = (i == 0) ? "3 x 10" : "4 x 25";
            double[] r = measured[i];
            assertTrue(at + ": one-way ANOVA spent " + r[ANOVA], r[ANOVA] > 0.04 && r[ANOVA] < 0.06);
            assertTrue(at + ": Welch ANOVA spent " + r[WELCH], r[WELCH] > 0.04 && r[WELCH] < 0.06);
            assertTrue(at + ": Bartlett spent " + r[BARTLETT], r[BARTLETT] > 0.04 && r[BARTLETT] < 0.06);
            assertTrue(at + ": Kruskal-Wallis spent " + r[RANK], r[RANK] > 0.035 && r[RANK] < 0.06);
        }
    }

    @Test
    public void testBartlettIsUnusableOnHeavyTails() {
        // the measurement that decides which test of equal variance to reach
        // for, recorded rather than forbidden. The groups below are equally
        // variable, so every rejection is a wrong one -- and Bartlett gets
        // worse as the sample grows, because what it is really detecting is the
        // tails and a larger sample shows it more of them
        double[] small = rates(5000, sizes(3, 10), each(3, 0.0), each(3, 1.0), true);
        double[] large = rates(5000, sizes(4, 25), each(4, 0.0), each(4, 1.0), true);
        // measured over 20000 samples: 0.3021 and 0.5406
        assertTrue("Bartlett spent " + small[BARTLETT] + " at three t(3) groups of ten",
                small[BARTLETT] > 0.25 && small[BARTLETT] < 0.36);
        assertTrue("Bartlett spent " + large[BARTLETT] + " at four t(3) groups of 25",
                large[BARTLETT] > 0.48 && large[BARTLETT] < 0.60);
        assertTrue("it did not get worse with the sample", large[BARTLETT] > small[BARTLETT]);

        // and the robust form on the very same samples stays under the level
        // -- measured: 0.0307 and 0.0349
        assertTrue("Brown-Forsythe spent " + small[FORSYTHE], small[FORSYTHE] < 0.05);
        assertTrue("Brown-Forsythe spent " + large[FORSYTHE], large[FORSYTHE] < 0.05);

        // and so do the two tests of the means, which is why heavy tails are
        // not a reason to reach past them -- measured: 0.0420 and 0.0419
        assertTrue("one-way ANOVA spent " + small[ANOVA], small[ANOVA] < 0.06);
        assertTrue("one-way ANOVA spent " + large[ANOVA], large[ANOVA] < 0.06);
    }

    @Test
    public void testLeveneOverspendsItsLevelAndBrownForsytheUnderspendsIt() {
        // Levene's distances are measured from a mean the same data provided,
        // so they are not independent of it and the analysis of variance run on
        // them is only approximate. It errs towards rejecting, which is the
        // direction that matters, and the median form errs the other way
        double[] ten = rates(10000, sizes(3, 10), each(3, 0.0), each(3, 1.0), false);
        double[] five = rates(10000, sizes(3, 5), each(3, 0.0), each(3, 1.0), false);
        // measured over 20000 samples: 0.0655, 0.0848, 0.0325 and 0.0046
        assertTrue("Levene spent " + ten[LEVENE] + " at three normal groups of ten",
                ten[LEVENE] > 0.05 && ten[LEVENE] < 0.08);
        assertTrue("Levene spent " + five[LEVENE] + " at three normal groups of five",
                five[LEVENE] > 0.05 && five[LEVENE] < 0.11);
        assertTrue("Brown-Forsythe spent " + ten[FORSYTHE], ten[FORSYTHE] > 0.02 && ten[FORSYTHE] < 0.05);
        assertTrue("Brown-Forsythe spent " + five[FORSYTHE] + " at groups of five", five[FORSYTHE] < 0.02);
        assertTrue("the median form was not the more conservative of the two",
                five[FORSYTHE] < five[LEVENE] && ten[FORSYTHE] < ten[LEVENE]);
    }

    @Test
    public void testTheClassicalAnovaFailsBothWaysOnUnequalVariances() {
        // the reason Welch's test is here. The means below are all equal, so
        // every rejection is wrong, and which way the classical test fails
        // depends on which group carries the spread: a small variable group
        // makes it reject far too often, a large one makes it stop rejecting
        // at all
        double[] wide = { 1.0, 1.0, 4.0 };
        double[] smallGroupIsWide = rates(5000, new int[] { 30, 10, 5 }, each(3, 0.0), wide, false);
        double[] largeGroupIsWide = rates(5000, new int[] { 5, 10, 30 }, each(3, 0.0), wide, false);
        // measured over 20000 samples: 0.3357 and 0.0014
        assertTrue("the classical test spent " + smallGroupIsWide[ANOVA] + " where the small group was wide",
                smallGroupIsWide[ANOVA] > 0.25);
        assertTrue("the classical test spent " + largeGroupIsWide[ANOVA] + " where the large group was wide",
                largeGroupIsWide[ANOVA] < 0.01);

        // measured: 0.0590 and 0.0519 -- from the same samples, so the two
        // tests are being compared on identical data and not merely on data
        // of the same shape
        assertTrue("Welch spent " + smallGroupIsWide[WELCH],
                smallGroupIsWide[WELCH] > 0.04 && smallGroupIsWide[WELCH] < 0.07);
        assertTrue("Welch spent " + largeGroupIsWide[WELCH],
                largeGroupIsWide[WELCH] > 0.04 && largeGroupIsWide[WELCH] < 0.07);

        // and the rank test is no answer to this one: it is a test of
        // distributions, and groups of different spread have different
        // distributions whatever their means do -- measured 0.1177 and 0.0103
        assertTrue("Kruskal-Wallis was unaffected by the unequal spreads",
                smallGroupIsWide[RANK] > 0.07 && largeGroupIsWide[RANK] < 0.03);
    }

    @Test
    public void testTheRankTestSeesTheShiftTheParametricOnesLose() {
        // the parametric and non-parametric twins, side by side. Under
        // normality the analysis of variance is the better of the two and the
        // rank test gives up a few points of power for assuming nothing; under
        // a t(3) that trade reverses, and by more
        int[] three = sizes(3, 10);
        double[] moved = { 0.0, 0.0, 1.0 };
        double[] normal = rates(5000, three, moved, each(3, 1.0), false);
        double[] heavy = rates(5000, three, moved, each(3, 1.0), true);
        double anovaNormal = normal[ANOVA];
        double rankNormal = normal[RANK];
        double anovaHeavy = heavy[ANOVA];
        double rankHeavy = heavy[RANK];

        assertTrue("normal: ANOVA found " + anovaNormal + ", the rank test " + rankNormal,
                anovaNormal > rankNormal);
        assertTrue("t(3): the rank test found " + rankHeavy + ", ANOVA " + anovaHeavy, rankHeavy > anovaHeavy);
        // measured over 20000 samples at three groups of ten and a shift of
        // one: 0.5764 against 0.5361 under normality, and 0.3004 against 0.3641
        // under a t(3). The gap the rank test wins by is the larger of the two
        assertTrue("the rank test gave up more than it gained", rankHeavy - anovaHeavy > anovaNormal - rankNormal);
    }

    // ---------------------------------------------------------- the whole thing --

    @Test
    public void testTheWorkflowFindsTheGroupThatMoved() {
        // what the slice is for, run once end to end: check the assumption,
        // ask the omnibus question, and follow a rejection with pairwise
        // comparisons corrected for having made three of them. Only the third
        // group is moved, so two of the three pairs are real and one is not
        int reps = 1000;
        int bothFound = 0;
        int spurious = 0;
        int spreadsCalledUnequal = 0;
        for (long seed = 1; seed <= reps; seed++) {
            double[][] g = groups(seed * 131L, sizes(3, 12), new double[] { 0.0, 0.0, 2.0 }, each(3, 1.0), false);
            // the variances are equal here, so the step that checks that should
            // let almost every sample through -- counted rather than asserted
            // per sample, since a five percent test is wrong one time in twenty
            // by construction and there are a thousand of them
            if (HypothesisTests.brownForsythe(g).test.rejectsAt(0.05)) {
                spreadsCalledUnequal++;
            }
            if (!HypothesisTests.welchAnova(g).test.rejectsAt(0.05)) {
                continue;
            }
            double[] pairwise = { HypothesisTests.tTwoSample(g[0], g[1], Alternative.TWO_SIDED, 0.95).test.pValue,
                    HypothesisTests.tTwoSample(g[0], g[2], Alternative.TWO_SIDED, 0.95).test.pValue,
                    HypothesisTests.tTwoSample(g[1], g[2], Alternative.TWO_SIDED, 0.95).test.pValue };
            double[] adjusted = MultipleTesting.holmBonferroni(pairwise);
            for (int i = 0; i < 3; i++) {
                assertTrue("seed=" + seed + ": the correction made a p-value smaller",
                        adjusted[i] >= pairwise[i]);
            }
            if (adjusted[1] <= 0.05 && adjusted[2] <= 0.05) {
                bothFound++;
            }
            if (adjusted[0] <= 0.05) {
                spurious++;
            }
        }
        // measured over 2000 runs: both real pairs found 98.2 percent of the
        // time, and the pair that is not real called different 4.6 percent --
        // which is the five percent the correction promises, over the one true
        // null hypothesis in the family
        assertTrue("only " + bothFound + " of " + reps + " found both real differences", bothFound > 0.95 * reps);
        assertTrue("the equal pair was called different " + spurious + " times", spurious < 0.07 * reps);
        assertTrue("the spreads were called unequal " + spreadsCalledUnequal + " times, and they are equal",
                spreadsCalledUnequal < 0.05 * reps);
    }

    // ------------------------------------------------------------ the corners --

    @Test
    public void testGroupsThatAreConstantAndDifferAreCertain() {
        // no spread inside any group and the groups do not agree: the within
        // sum is zero, the statistic is infinite and the p-value is exactly
        // zero. That is a real answer, not a failure, and it is the same
        // decision pearsonCorrelation takes for an infinite t
        FTestResult f = HypothesisTests.oneWayAnova(new double[][] { { 1.0, 1.0, 1.0 }, { 2.0, 2.0, 2.0 } });
        assertEquals("the statistic", Double.POSITIVE_INFINITY, f.test.statistic, 0.0);
        assertEquals("the p-value", 0.0, f.test.pValue, 0.0);
        assertTrue("it did not reject", f.test.rejectsAt(1.0e-300));
    }

    @Test
    public void testBartlettIsExactlyZeroWhenTheSpreadsAgree() {
        // three groups that are shifts of one another have identical variances,
        // so the statistic is zero by Jensen's inequality with equality. What
        // rounding leaves a little below zero there would read as a p-value of
        // zero, since chiSquareComplemented takes a negative argument to be off
        // the scale
        double[][] g = { { 1.0, 2.0, 3.0, 4.0 }, { 11.0, 12.0, 13.0, 14.0 }, { -6.0, -5.0, -4.0, -3.0 } };
        TestResult b = HypothesisTests.bartlett(g);
        assertEquals("the statistic", 0.0, b.statistic, 0.0);
        assertEquals("the p-value", 1.0, b.pValue, 0.0);
        assertEquals("the degrees of freedom", 2.0, b.degreesOfFreedom, 0.0);
    }

    @Test
    public void testTheResultsSayWhichTestRanAndOnWhatDegreesOfFreedom() {
        double[][] g = nullGroups(3L, 4, 7, false);
        FTestResult anova = HypothesisTests.oneWayAnova(g);
        assertEquals("the name", "one-way ANOVA", anova.test.test);
        assertEquals("numerator df", 3.0, anova.numeratorDf, 0.0);
        assertEquals("denominator df", 24.0, anova.denominatorDf, 0.0);
        assertEquals("the alternative", Alternative.GREATER, anova.test.alternative);
        assertTrue("the degrees of freedom are not in the TestResult",
                Double.isNaN(anova.test.degreesOfFreedom));
        assertTrue("the two degrees of freedom are not printed", anova.toString().contains("df = (3, 24)"));

        FTestResult welch = HypothesisTests.welchAnova(g);
        assertEquals("the name", "Welch ANOVA", welch.test.test);
        assertEquals("numerator df", 3.0, welch.numeratorDf, 0.0);
        // Welch's is fractional, which is what FTestResult carries doubles for,
        // and it never exceeds what the classical test would have used
        assertTrue("Welch's denominator df was " + welch.denominatorDf,
                welch.denominatorDf > 3.0 && welch.denominatorDf < 24.0);
        assertTrue("a whole number came out of Welch's denominator",
                welch.denominatorDf != Math.rint(welch.denominatorDf));

        assertEquals("the name", "Bartlett", HypothesisTests.bartlett(g).test);
        assertEquals("the name", "Levene", HypothesisTests.levene(g).test.test);
        assertEquals("the name", "Brown-Forsythe", HypothesisTests.brownForsythe(g).test.test);
    }

    // ---------------------------------------------------------- the guard rail --

    @Test
    public void testTheKSampleTestsRejectWhatTheyCannotCompute() {
        refusedByAll("a null argument", null);
        refusedByAll("one group", new double[][] { { 1.0, 2.0, 3.0 } });
        refusedByAll("a null group", new double[][] { { 1.0, 2.0 }, null });
        refusedByAll("an empty group", new double[][] { { 1.0, 2.0 }, {} });
        refusedByAll("an infinite observation",
                new double[][] { { 1.0, 2.0 }, { 3.0, Double.POSITIVE_INFINITY } });
        refusedByAll("a NaN observation", new double[][] { { 1.0, 2.0 }, { 3.0, Double.NaN } });

        // one value everywhere leaves nothing to analyse for any of them
        refusedByAll("one value in all the groups", new double[][] { { 5.0, 5.0 }, { 5.0, 5.0 } });

        // the two that need a variance per group need two observations in each
        // and some spread in each; the other three are content without
        double[][] singleton = { { 1.0, 2.0, 3.0 }, { 4.0 } };
        refuses("welchAnova on a group of one", singleton, 1);
        refuses("bartlett on a group of one", singleton, 2);
        assertTrue("oneWayAnova refused a group of one",
                HypothesisTests.oneWayAnova(singleton).test.pValue >= 0.0);
        assertTrue("brownForsythe refused a group of one",
                HypothesisTests.brownForsythe(singleton).test.pValue >= 0.0);

        double[][] flat = { { 1.0, 2.0, 3.0 }, { 4.0, 4.0, 4.0 } };
        refuses("welchAnova on a group with no spread", flat, 1);
        refuses("bartlett on a group with no spread", flat, 2);
        assertTrue("oneWayAnova refused a group with no spread",
                HypothesisTests.oneWayAnova(flat).test.pValue >= 0.0);

        // and an F test needs something left over once the k means are fitted
        refuses("oneWayAnova with one observation per group", new double[][] { { 1.0 }, { 2.0 } }, 0);
        refuses("levene with one observation per group", new double[][] { { 1.0 }, { 2.0 } }, 3);
    }

    private static void refusedByAll(String what, double[][] groups) {
        for (int which = 0; which <= 4; which++) {
            refuses(what, groups, which);
        }
    }

    private static void refuses(String what, double[][] groups, int which) {
        String name;
        try {
            switch (which) {
            case 0:
                name = "oneWayAnova";
                HypothesisTests.oneWayAnova(groups);
                break;
            case 1:
                name = "welchAnova";
                HypothesisTests.welchAnova(groups);
                break;
            case 2:
                name = "bartlett";
                HypothesisTests.bartlett(groups);
                break;
            case 3:
                name = "levene";
                HypothesisTests.levene(groups);
                break;
            default:
                name = "brownForsythe";
                HypothesisTests.brownForsythe(groups);
                break;
            }
            fail(name + " accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        } catch (NullPointerException fromInside) {
            fail("a null slipped through the guard rail for " + what);
        }
    }

    @Test
    public void testTheInputIsNotModified() {
        double[][] g = nullGroups(21L, 3, 6, false);
        double[][] before = new double[3][];
        for (int i = 0; i < 3; i++) {
            before[i] = g[i].clone();
        }
        HypothesisTests.oneWayAnova(g);
        HypothesisTests.welchAnova(g);
        HypothesisTests.bartlett(g);
        HypothesisTests.levene(g);
        HypothesisTests.brownForsythe(g);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < g[i].length; j++) {
                assertEquals("group " + i + " observation " + j, before[i][j], g[i][j], 0.0);
            }
        }
    }
}
