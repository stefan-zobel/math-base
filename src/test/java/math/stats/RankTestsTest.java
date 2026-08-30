package math.stats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.ContinuousDistribution;
import math.distribution.Exponential;
import math.distribution.Normal;
import math.stats.rank.MannWhitneyU;
import math.stats.rank.Ranks;
import math.stats.rank.WilcoxonSignedRank;

/**
 * A rank test keeps the order of the observations and throws the values away.
 * That is what makes it checkable without a table: the null distribution is a
 * count of equally likely orderings, so it can be counted a second way, and
 * the test itself cannot notice any strictly increasing transform of its own
 * data.
 */
public final class RankTestsTest {

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

    // --------------------------------------------------------- the ranking --

    @Test
    public void testTiedValuesShareTheAverageOfTheRanksTheySpan() {
        double[] v = { 3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0 };
        double[] want = { 4.5, 1.5, 6.0, 1.5, 8.0, 11.0, 3.0, 10.0, 8.0, 4.5, 8.0 };
        Ranks.Result r = Ranks.of(v);
        for (int i = 0; i < want.length; i++) {
            assertEquals("rank of " + v[i], want[i], r.ranks[i], 0.0);
        }
        // two pairs and one triple: 6 + 6 + 24
        assertEquals("tieSum", 36.0, r.tieSum, 0.0);
        assertTrue("hasTies", r.hasTies);
    }

    @Test
    public void testTheRanksSumToTheSameTotalWhetherOrNotThereAreTies() {
        // the midrank of a run is the average of the ranks it consumes, so no
        // amount of tying can move the total
        for (int n : new int[] { 1, 2, 7, 50, 301 }) {
            for (long seed = 1; seed <= 40; seed++) {
                double[] x = sample(new Normal(0.0, 1.0), n, seed * 7919L + 1);
                Ranks.Result plain = Ranks.of(x);
                double[] coarse = x.clone();
                for (int i = 0; i < coarse.length; i++) {
                    coarse[i] = Math.rint(coarse[i]);
                }
                Ranks.Result tied = Ranks.of(coarse);
                double want = n * (n + 1.0) / 2.0;
                assertEquals("n=" + n + " untied", want, sum(plain.ranks), 1.0e-9);
                assertEquals("n=" + n + " tied", want, sum(tied.ranks), 1.0e-9);
                assertEquals("n=" + n + ": an untied sample has no tie sum", 0.0, plain.tieSum, 0.0);
                assertFalse("n=" + n + ": an untied sample has no ties", plain.hasTies);
            }
        }
    }

    @Test
    public void testTheRanksDoNotChangeUnderAnIncreasingTransform() {
        // bit for bit: a rank cannot see anything but the order
        for (long seed = 1; seed <= 200; seed++) {
            double[] x = sample(new Normal(0.0, 2.0), 23, seed * 7919L + 1);
            double[] moved = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                moved[i] = Math.exp(x[i]);
            }
            double[] a = Ranks.of(x).ranks;
            double[] b = Ranks.of(moved).ranks;
            for (int i = 0; i < a.length; i++) {
                assertEquals("seed=" + seed + " i=" + i, a[i], b[i], 0.0);
            }
        }
    }

    @Test
    public void testTheRankingRejectsWhatHasNoOrder() {
        rejects("null", null);
        rejects("empty", new double[0]);
        rejects("NaN", new double[] { 1.0, Double.NaN, 3.0 });
        // an infinity does have a place in an ordering, and keeps it
        Ranks.Result r = Ranks.of(new double[] { 1.0, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY });
        assertEquals("the finite value is in the middle", 2.0, r.ranks[0], 0.0);
        assertEquals("+Infinity is the largest", 3.0, r.ranks[1], 0.0);
        assertEquals("-Infinity is the smallest", 1.0, r.ranks[2], 0.0);
    }

    private static void rejects(String what, double[] values) {
        try {
            Ranks.of(values);
            fail(what + ": accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
    }

    private static double sum(double[] x) {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += x[i];
        }
        return s;
    }

    // ------------------------------------------- the exact Mann-Whitney null --

    /** {@code P[U >= u]} for every {@code u}, counted over every subset. */
    private static double[] bruteTails(int m, int n) {
        int total = m + n;
        double[] counts = new double[m * n + 1];
        double all = 0.0;
        for (int mask = 0; mask < (1 << total); mask++) {
            if (Integer.bitCount(mask) != m) {
                continue;
            }
            int s = 0;
            for (int b = 0; b < total; b++) {
                if ((mask & (1 << b)) != 0) {
                    s += b + 1;
                }
            }
            counts[s - m * (m + 1) / 2] += 1.0;
            all += 1.0;
        }
        double[] tails = new double[m * n + 1];
        double reached = 0.0;
        for (int u = counts.length - 1; u >= 0; u--) {
            reached += counts[u];
            tails[u] = reached / all;
        }
        return tails;
    }

    @Test
    public void testTheExactNullIsWhatCountingTheOrderingsGives() {
        // the independent algorithm: enumerate every way of dealing the ranks
        // out to the two samples. Measured worst departure 4.4e-16
        for (int m = 1; m <= 13; m++) {
            for (int n = 1; m + n <= 14; n++) {
                double[] want = bruteTails(m, n);
                for (int u = 0; u <= m * n; u++) {
                    String at = "m=" + m + " n=" + n + " u=" + u;
                    assertEquals(at, want[u], MannWhitneyU.barFExact(m, n, u), 1.0e-15);
                }
            }
        }
    }

    @Test
    public void testTheExactNullDoesNotCareWhichSampleIsRanked() {
        // U depends on the two sizes only as an unordered pair, which is what
        // lets the recursion always carry the smaller table -- bit for bit
        int[][] shapes = { { 2, 9 }, { 3, 20 }, { 5, 40 }, { 7, 60 }, { 10, 100 }, { 1, 300 }, { 12, 13 } };
        for (int i = 0; i < shapes.length; i++) {
            int m = shapes[i][0];
            int n = shapes[i][1];
            for (int u = 0; u <= m * n; u += 7) {
                assertEquals("m=" + m + " n=" + n + " u=" + u, MannWhitneyU.barFExact(m, n, u),
                        MannWhitneyU.barFExact(n, m, u), 0.0);
            }
        }
    }

    @Test
    public void testTheTwoTailsOfTheExactNullPartitionIt() {
        // a discrete statistic counts its own value in both tails, so the two
        // of them come to one plus that value and not to one. Measured worst
        // departure 4.4e-16
        for (int m = 1; m <= 11; m++) {
            for (int n = 1; m + n <= 18; n++) {
                double mn = (double) m * n;
                for (int u = 0; u <= m * n; u++) {
                    double upper = MannWhitneyU.barFExact(m, n, u);
                    double lower = MannWhitneyU.barFExact(m, n, mn - u);
                    double atU = upper - MannWhitneyU.barFExact(m, n, u + 1.0);
                    assertEquals("m=" + m + " n=" + n + " u=" + u, 1.0 + atU, lower + upper, 1.0e-15);
                }
            }
        }
    }

    @Test
    public void testTheExactNullEndsWhereItShould() {
        for (int m = 1; m <= 6; m++) {
            for (int n = 1; n <= 6; n++) {
                String at = "m=" + m + " n=" + n;
                assertEquals(at + ": nothing is below zero", 1.0, MannWhitneyU.barFExact(m, n, 0.0), 0.0);
                assertEquals(at + ": nothing is above mn", 0.0, MannWhitneyU.barFExact(m, n, m * n + 1.0), 0.0);
                assertTrue(at + ": the largest value is attainable", MannWhitneyU.barFExact(m, n, m * n) > 0.0);
            }
        }
    }

    @Test
    public void testTheExactNullRefusesWhatItCannotAfford() {
        // 2500 is the measured limit: 11 ms at its worst shape, m = n = 50
        assertTrue("50 x 50 is inside", MannWhitneyU.barFExact(50, 50, 1250.0) > 0.0);
        try {
            MannWhitneyU.barFExact(51, 50, 1275.0);
            fail("a product above the limit was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message names the limit", expected.getMessage().contains("2500"));
        }
        try {
            MannWhitneyU.barFExact(0, 5, 0.0);
            fail("an empty sample was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        try {
            MannWhitneyU.barFAsymptotic(5, 5, 12.0, -1.0);
            fail("a negative tie sum was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
    }

    @Test
    public void testTheNormalApproximationClosesOnTheExactNull() {
        // over the p-values a reader acts on, and it is the shape rather than
        // the product that governs: measured 4.2e-3 at 10 x 10 and 8.2e-4 at
        // 50 x 50
        int[][] shapes = { { 10, 10 }, { 20, 20 }, { 30, 30 }, { 50, 50 }, { 10, 250 }, { 25, 100 } };
        double[] allowed = { 5.0e-3, 2.5e-3, 1.6e-3, 1.0e-3, 3.3e-3, 1.4e-3 };
        for (int i = 0; i < shapes.length; i++) {
            int m = shapes[i][0];
            int n = shapes[i][1];
            double worst = 0.0;
            // the window sits above the middle of a null that is symmetric
            // about it, and the exact tail costs a whole table per call, so it
            // is walked at a stride rather than one value at a time
            int stride = Math.max(1, m * n / 80);
            for (int u = m * n / 2; u <= m * n; u += stride) {
                double exact = MannWhitneyU.barFExact(m, n, u);
                if (exact < 0.001) {
                    break;
                }
                if (exact > 0.2) {
                    continue;
                }
                worst = Math.max(worst, Math.abs(MannWhitneyU.barFAsymptotic(m, n, u, 0.0) - exact));
            }
            assertTrue("m=" + m + " n=" + n + ": worst difference " + worst, worst > 0.0 && worst <= allowed[i]);
        }
    }

    @Test
    public void testEverythingTiedLeavesTheApproximationWithNothingToSay() {
        // the tie correction takes the whole variance away, and a tail of one
        // is the only honest answer left
        int m = 6;
        int n = 9;
        double total = m + n;
        double everything = total * total * total - total;
        assertEquals(1.0, MannWhitneyU.barFAsymptotic(m, n, 27.0, everything), 0.0);
    }

    // ------------------------------------------------- the two-sample test --

    /** {@code U} counted pair by pair, a tie counting half. */
    private static double dominatingPairs(double[] x, double[] y) {
        double u = 0.0;
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                if (x[i] > y[j]) {
                    u += 1.0;
                } else if (x[i] == y[j]) {
                    u += 0.5;
                }
            }
        }
        return u;
    }

    @Test
    public void testTheStatisticCountsThePairsItSaysItCounts() {
        // the rank sum and the pair count are the same number by an identity,
        // and the second way round is the independent algorithm
        for (long seed = 1; seed <= 300; seed++) {
            double[] x = sample(new Normal(0.0, 1.0), 9, seed * 7919L + 1);
            double[] y = sample(new Exponential(1.0), 12, seed * 104729L + 7);
            if (seed % 2 == 0) {
                for (int i = 0; i < x.length; i++) {
                    x[i] = Math.rint(x[i]);
                }
                for (int i = 0; i < y.length; i++) {
                    y[i] = Math.rint(y[i]);
                }
            }
            TestResult r = HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED);
            assertEquals("seed=" + seed, dominatingPairs(x, y), r.statistic, 1.0e-9);
        }
    }

    @Test
    public void testTheTestIsBlindToAnIncreasingTransform() {
        // bit for bit, statistic and p-value alike, on both paths
        for (long seed = 1; seed <= 200; seed++) {
            double[] x = sample(new Normal(2.0, 1.0), 11, seed * 7919L + 1);
            double[] y = sample(new Normal(2.3, 1.0), 13, seed * 104729L + 7);
            double[] ex = new double[x.length];
            double[] ey = new double[y.length];
            for (int i = 0; i < x.length; i++) {
                ex[i] = Math.exp(x[i]);
            }
            for (int i = 0; i < y.length; i++) {
                ey[i] = Math.exp(y[i]);
            }
            for (Alternative alternative : Alternative.values()) {
                TestResult a = HypothesisTests.mannWhitneyU(x, y, alternative);
                TestResult b = HypothesisTests.mannWhitneyU(ex, ey, alternative);
                String at = "seed=" + seed + " " + alternative;
                assertEquals(at + ": statistic", a.statistic, b.statistic, 0.0);
                assertEquals(at + ": p-value", a.pValue, b.pValue, 0.0);
                assertEquals(at + ": path", a.test, b.test);
            }
        }
    }

    @Test
    public void testSwappingTheSamplesSwapsTheOneSidedTails() {
        // U_y is mn - U_x, so the statement "x is the larger" run one way is
        // the statement "y is the larger" run the other -- the same number
        for (long seed = 1; seed <= 200; seed++) {
            double[] x = sample(new Normal(0.0, 1.0), 8, seed * 7919L + 1);
            double[] y = sample(new Normal(0.4, 1.5), 11, seed * 104729L + 7);
            String at = "seed=" + seed;
            assertEquals(at, HypothesisTests.mannWhitneyU(x, y, Alternative.GREATER).pValue,
                    HypothesisTests.mannWhitneyU(y, x, Alternative.LESS).pValue, 0.0);
            assertEquals(at, HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED).pValue,
                    HypothesisTests.mannWhitneyU(y, x, Alternative.TWO_SIDED).pValue, 0.0);
            assertEquals(at + ": the statistics add to mn", 8 * 11,
                    HypothesisTests.mannWhitneyU(x, y, Alternative.LESS).statistic
                            + HypothesisTests.mannWhitneyU(y, x, Alternative.LESS).statistic,
                    1.0e-9);
        }
    }

    @Test
    public void testTheTwoSidedPValueIsTwiceTheSmallerTail() {
        for (long seed = 1; seed <= 400; seed++) {
            double[] x = sample(new Normal(0.0, 1.0), 7, seed * 7919L + 1);
            double[] y = sample(new Exponential(1.0), 9, seed * 104729L + 7);
            double less = HypothesisTests.mannWhitneyU(x, y, Alternative.LESS).pValue;
            double greater = HypothesisTests.mannWhitneyU(x, y, Alternative.GREATER).pValue;
            double two = HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED).pValue;
            assertEquals("seed=" + seed, Math.min(1.0, 2.0 * Math.min(less, greater)), two, 0.0);
        }
    }

    @Test
    public void testATieOrASampleTooLargeSendsTheTestDownTheOtherPath() {
        double[] x = { 1.0, 2.0, 3.0, 4.0 };
        double[] y = { 5.0, 6.0, 7.0 };
        assertEquals("no ties, small samples", "Mann-Whitney U, exact",
                HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED).test);
        double[] tied = { 1.0, 2.0, 3.0, 5.0 };
        assertEquals("one tie is enough", "Mann-Whitney U, asymptotic",
                HypothesisTests.mannWhitneyU(tied, y, Alternative.TWO_SIDED).test);
        double[] big = sample(new Normal(0.0, 1.0), 51, 12345L);
        double[] other = sample(new Normal(0.0, 1.0), 50, 999L);
        assertEquals("a product of 2550 is past the limit", "Mann-Whitney U, asymptotic",
                HypothesisTests.mannWhitneyU(big, other, Alternative.TWO_SIDED).test);
        double[] just = sample(new Normal(0.0, 1.0), 50, 12345L);
        assertEquals("a product of 2500 is not", "Mann-Whitney U, exact",
                HypothesisTests.mannWhitneyU(just, other, Alternative.TWO_SIDED).test);
    }

    @Test
    public void testThePValueIsValidUnderTheNullHypothesis() {
        // a rank test on a discrete null cannot spend its whole significance
        // level, so the statement is not that the p-value is uniform but that
        // it does not overspend. Measured over 20000 replications, at the five
        // percent level: 0.0288 at 4 against 4, 0.0429 at 10 against 10, 0.0494
        // at 50 against 50 and 0.0505 at 200 against 200
        int[][] shapes = { { 4, 4 }, { 5, 7 }, { 10, 10 }, { 60, 60 }, { 200, 200 } };
        boolean[] rounded = { false, false, false, false, true };
        int reps = 8000;
        for (int i = 0; i < shapes.length; i++) {
            int m = shapes[i][0];
            int n = shapes[i][1];
            int below01 = 0;
            int below05 = 0;
            int below10 = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(0.0, 1.0), m, seed * 7919L + 1);
                double[] y = sample(new Normal(0.0, 1.0), n, seed * 104729L + 7);
                if (rounded[i]) {
                    for (int k = 0; k < x.length; k++) {
                        x[k] = Math.rint(2.0 * x[k]) / 2.0;
                    }
                    for (int k = 0; k < y.length; k++) {
                        y[k] = Math.rint(2.0 * y[k]) / 2.0;
                    }
                }
                double p = HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED).pValue;
                assertTrue("p = " + p, p >= 0.0 && p <= 1.0);
                if (p < 0.01) {
                    below01++;
                }
                if (p < 0.05) {
                    below05++;
                }
                if (p < 0.10) {
                    below10++;
                }
            }
            String at = "m=" + m + " n=" + n;
            // the Monte Carlo standard error at 8000 replications is 0.0011 at
            // the one percent level and 0.0024 at the five percent level
            assertTrue(at + " P(p<0.01) = " + below01 / (double) reps, below01 / (double) reps <= 0.01 + 0.004);
            assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps <= 0.05 + 0.008);
            assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps <= 0.10 + 0.012);
            if (m >= 10) {
                // and a null fine enough to spend it does spend it
                assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps >= 0.05 - 0.012);
                assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps >= 0.10 - 0.018);
            }
        }
    }

    @Test
    public void testItSeesASeparationThatIsThere() {
        // not a calibration, a sanity check: two samples that barely overlap
        double[] low = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] high = { 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0 };
        TestResult apart = HypothesisTests.mannWhitneyU(low, high, Alternative.LESS);
        assertTrue("a complete separation is not significant: " + apart.pValue, apart.pValue < 0.001);
        assertEquals("U is zero when every pair goes the same way", 0.0, apart.statistic, 0.0);
        TestResult wrongWay = HypothesisTests.mannWhitneyU(low, high, Alternative.GREATER);
        assertTrue("the other tail is nearly all of it: " + wrongWay.pValue, wrongWay.pValue > 0.999);
        double[] mixed = { 1.0, 4.0, 5.0, 8.0, 9.0, 12.0, 13.0 };
        double[] between = { 2.0, 3.0, 6.0, 7.0, 10.0, 11.0, 14.0 };
        TestResult same = HypothesisTests.mannWhitneyU(mixed, between, Alternative.TWO_SIDED);
        assertTrue("interleaved samples are not evidence: " + same.pValue, same.pValue > 0.5);
    }

    @Test
    public void testTheTestRejectsWhatItCannotAnswer() {
        double[] ok = { 1.0, 2.0, 3.0 };
        refuses("x null", null, ok, Alternative.TWO_SIDED);
        refuses("y null", ok, null, Alternative.TWO_SIDED);
        refuses("x empty", new double[0], ok, Alternative.TWO_SIDED);
        refuses("y empty", ok, new double[0], Alternative.TWO_SIDED);
        refuses("x not finite", new double[] { 1.0, Double.NaN }, ok, Alternative.TWO_SIDED);
        refuses("y infinite", ok, new double[] { Double.POSITIVE_INFINITY }, Alternative.TWO_SIDED);
        refuses("no alternative", ok, ok, null);
        refuses("no spread at all", new double[] { 2.0, 2.0 }, new double[] { 2.0 }, Alternative.TWO_SIDED);
        // one observation each is a sample this test can still answer, which
        // is not true of the t-test beside it
        TestResult r = HypothesisTests.mannWhitneyU(new double[] { 1.0 }, new double[] { 2.0 },
                Alternative.TWO_SIDED);
        assertEquals("U", 0.0, r.statistic, 0.0);
        assertEquals("nothing can be significant with one observation each", 1.0, r.pValue, 0.0);
    }

    private static void refuses(String what, double[] x, double[] y, Alternative alternative) {
        try {
            HypothesisTests.mannWhitneyU(x, y, alternative);
            fail(what + ": accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
    }

    // ------------------------------------------- the exact signed rank null --

    /** {@code P[W+ >= w]} for every {@code w}, counted over every sign pattern. */
    private static double[] bruteSignedTails(int n) {
        int most = n * (n + 1) / 2;
        double[] counts = new double[most + 1];
        double all = 0.0;
        for (int mask = 0; mask < (1 << n); mask++) {
            int s = 0;
            for (int b = 0; b < n; b++) {
                if ((mask & (1 << b)) != 0) {
                    s += b + 1;
                }
            }
            counts[s] += 1.0;
            all += 1.0;
        }
        double[] tails = new double[most + 1];
        double reached = 0.0;
        for (int w = most; w >= 0; w--) {
            reached += counts[w];
            tails[w] = reached / all;
        }
        return tails;
    }

    @Test
    public void testTheExactSignedRankNullIsWhatCountingTheSignPatternsGives() {
        // bit for bit, and it can be: every probability here is a count over
        // 2^n with both of them well inside the exact integers of a double, so
        // neither the halving recursion nor the counting rounds anything at
        // all. Measured over n = 1 .. 20
        for (int n = 1; n <= 18; n++) {
            double[] want = bruteSignedTails(n);
            for (int w = 0; w < want.length; w++) {
                assertEquals("n=" + n + " w=" + w, want[w], WilcoxonSignedRank.barFExact(n, w), 0.0);
            }
        }
    }

    @Test
    public void testTheTwoTailsOfTheSignedRankNullPartitionIt() {
        // measured worst departure 1.6e-15 over n = 1 .. 60
        for (int n = 1; n <= 40; n++) {
            double most = n * (n + 1.0) / 2.0;
            for (int w = 0; w <= most; w++) {
                double upper = WilcoxonSignedRank.barFExact(n, w);
                double lower = WilcoxonSignedRank.barFExact(n, most - w);
                double atW = upper - WilcoxonSignedRank.barFExact(n, w + 1.0);
                assertEquals("n=" + n + " w=" + w, 1.0 + atW, lower + upper, 5.0e-15);
            }
        }
    }

    @Test
    public void testTheSignedRankApproximationClosesOnTheExactNull() {
        // measured 8.0e-3 at n = 10 and 8.0e-4 at n = 100, over the p-values a
        // reader acts on. The exact tail costs a whole table per call, so the
        // window is walked at a stride
        int[] sizes = { 10, 15, 20, 30, 50, 100 };
        double[] allowed = { 9.0e-3, 6.0e-3, 4.5e-3, 3.0e-3, 1.8e-3, 1.0e-3 };
        for (int i = 0; i < sizes.length; i++) {
            int n = sizes[i];
            int most = n * (n + 1) / 2;
            int stride = Math.max(1, most / 60);
            double worst = 0.0;
            for (int w = most / 2; w <= most; w += stride) {
                double exact = WilcoxonSignedRank.barFExact(n, w);
                if (exact < 0.001) {
                    break;
                }
                if (exact > 0.2) {
                    continue;
                }
                worst = Math.max(worst, Math.abs(WilcoxonSignedRank.barFAsymptotic(n, w, 0.0) - exact));
            }
            assertTrue("n=" + n + ": worst difference " + worst, worst > 0.0 && worst <= allowed[i]);
        }
    }

    @Test
    public void testTheSignedRankNullRefusesWhatItCannotAfford() {
        assertTrue("300 is inside", WilcoxonSignedRank.barFExact(300, 22500.0) > 0.0);
        try {
            WilcoxonSignedRank.barFExact(301, 22500.0);
            fail("a sample above the limit was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message names the limit", expected.getMessage().contains("300"));
        }
        try {
            WilcoxonSignedRank.barFExact(0, 0.0);
            fail("an empty sample was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        try {
            WilcoxonSignedRank.barFAsymptotic(5, 7.0, 1.0e9);
            fail("an impossible tie sum was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        // the ties can take the variance down but never all the way: the signs
        // are still free, whatever the magnitudes do
        double everything = 5.0 * 5.0 * 5.0 - 5.0;
        double tail = WilcoxonSignedRank.barFAsymptotic(5, 15.0, everything);
        assertTrue("a tail of " + tail, tail > 0.0 && tail < 1.0);
    }

    // ---------------------------------------------- the one-sample test --

    /**
     * {@code W+} counted as the Walsh averages that come out positive: for a
     * pair {@code i <= j} the sign of {@code d_i + d_j} is the sign of
     * whichever of the two is the larger in magnitude, so the count of them is
     * the rank sum of the positive differences without a ranking anywhere.
     */
    private static double positiveWalshAverages(double[] d) {
        double w = 0.0;
        for (int i = 0; i < d.length; i++) {
            for (int k = i; k < d.length; k++) {
                if (d[i] + d[k] > 0.0) {
                    w += 1.0;
                }
            }
        }
        return w;
    }

    @Test
    public void testTheStatisticIsTheRankSumItSaysItIs() {
        for (long seed = 1; seed <= 300; seed++) {
            double[] x = sample(new Normal(0.3, 1.0), 17, seed * 7919L + 1);
            TestResult r = HypothesisTests.wilcoxonSignedRank(x, 0.0, Alternative.TWO_SIDED);
            assertEquals("seed=" + seed, positiveWalshAverages(x), r.statistic, 0.0);
        }
    }

    @Test
    public void testDifferencesOfExactlyZeroAreDroppedBeforeTheRanking() {
        double[] withZeros = { 0.0, 1.5, 0.0, -2.5, 3.5, 0.0, 4.5 };
        double[] without = { 1.5, -2.5, 3.5, 4.5 };
        for (Alternative alternative : Alternative.values()) {
            TestResult a = HypothesisTests.wilcoxonSignedRank(withZeros, 0.0, alternative);
            TestResult b = HypothesisTests.wilcoxonSignedRank(without, 0.0, alternative);
            assertEquals(alternative + ": statistic", b.statistic, a.statistic, 0.0);
            assertEquals(alternative + ": p-value", b.pValue, a.pValue, 0.0);
            assertEquals(alternative + ": path", b.test, a.test);
        }
        // three of the seven observations bought nothing: the answer is the
        // one for a sample of four
        TestResult r = HypothesisTests.wilcoxonSignedRank(withZeros, 0.0, Alternative.TWO_SIDED);
        assertEquals("W+", 8.0, r.statistic, 0.0);
        assertEquals("p", 0.375, r.pValue, 1.0e-15);
    }

    @Test
    public void testTheSignedRankTestIsBlindToAnOddIncreasingTransform() {
        // the cube keeps every sign and every ordering of the magnitudes, so
        // it keeps the statistic and the p-value bit for bit. The exponential
        // would not: it does not fix the zero the signs are taken about
        for (long seed = 1; seed <= 300; seed++) {
            double[] x = sample(new Normal(0.4, 1.0), 17, seed * 7919L + 1);
            double[] cubed = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                cubed[i] = x[i] * x[i] * x[i];
            }
            for (Alternative alternative : Alternative.values()) {
                TestResult a = HypothesisTests.wilcoxonSignedRank(x, 0.0, alternative);
                TestResult b = HypothesisTests.wilcoxonSignedRank(cubed, 0.0, alternative);
                String at = "seed=" + seed + " " + alternative;
                assertEquals(at + ": statistic", a.statistic, b.statistic, 0.0);
                assertEquals(at + ": p-value", a.pValue, b.pValue, 0.0);
            }
        }
    }

    @Test
    public void testThePairedFormIsTheOneSampleFormOfTheDifferences() {
        for (long seed = 1; seed <= 200; seed++) {
            double[] x = sample(new Normal(1.0, 1.0), 14, seed * 7919L + 1);
            double[] y = sample(new Normal(0.7, 1.3), 14, seed * 104729L + 7);
            double[] differences = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                differences[i] = x[i] - y[i];
            }
            for (Alternative alternative : Alternative.values()) {
                TestResult a = HypothesisTests.wilcoxonSignedRankPaired(x, y, alternative);
                TestResult b = HypothesisTests.wilcoxonSignedRank(differences, 0.0, alternative);
                String at = "seed=" + seed + " " + alternative;
                assertEquals(at + ": statistic", b.statistic, a.statistic, 0.0);
                assertEquals(at + ": p-value", b.pValue, a.pValue, 0.0);
            }
        }
    }

    @Test
    public void testTheSignedRankPValueIsValidUnderTheNullHypothesis() {
        // measured over 20000 replications at the five percent level: 0.0000
        // at n = 5, where the null is too coarse to spend it, 0.0398 at n = 8,
        // 0.0506 at n = 15 and 0.0510 at n = 400
        int[] sizes = { 5, 8, 15, 40, 400 };
        int reps = 8000;
        for (int i = 0; i < sizes.length; i++) {
            int n = sizes[i];
            int below01 = 0;
            int below05 = 0;
            int below10 = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(new Normal(0.0, 1.0), n, seed * 7919L + 1);
                double p = HypothesisTests.wilcoxonSignedRank(x, 0.0, Alternative.TWO_SIDED).pValue;
                assertTrue("p = " + p, p >= 0.0 && p <= 1.0);
                if (p < 0.01) {
                    below01++;
                }
                if (p < 0.05) {
                    below05++;
                }
                if (p < 0.10) {
                    below10++;
                }
            }
            String at = "n=" + n;
            assertTrue(at + " P(p<0.01) = " + below01 / (double) reps, below01 / (double) reps <= 0.01 + 0.005);
            assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps <= 0.05 + 0.008);
            assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps <= 0.10 + 0.015);
            if (n >= 15) {
                assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps >= 0.05 - 0.012);
                assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps >= 0.10 - 0.018);
            }
        }
    }

    @Test
    public void testTheSignedRankTestRejectsWhatItCannotAnswer() {
        double[] ok = { 1.0, 2.0, 3.0 };
        refusesSignedRank("null", null, 0.0, Alternative.TWO_SIDED);
        refusesSignedRank("empty", new double[0], 0.0, Alternative.TWO_SIDED);
        refusesSignedRank("not finite", new double[] { 1.0, Double.NaN }, 0.0, Alternative.TWO_SIDED);
        refusesSignedRank("mu not finite", ok, Double.NaN, Alternative.TWO_SIDED);
        refusesSignedRank("mu infinite", ok, Double.POSITIVE_INFINITY, Alternative.TWO_SIDED);
        refusesSignedRank("no alternative", ok, 0.0, null);
        refusesSignedRank("nothing but zeros", new double[] { 2.0, 2.0, 2.0 }, 2.0, Alternative.TWO_SIDED);
        try {
            HypothesisTests.wilcoxonSignedRankPaired(ok, new double[] { 1.0, 2.0 }, Alternative.TWO_SIDED);
            fail("unpaired samples were accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        // one observation is a sample this test can answer, and one that says
        // nothing: a single difference is as likely to be positive as not
        TestResult r = HypothesisTests.wilcoxonSignedRank(new double[] { 5.0 }, 0.0, Alternative.TWO_SIDED);
        assertEquals("W+", 1.0, r.statistic, 0.0);
        assertEquals("p", 1.0, r.pValue, 0.0);
    }

    private static void refusesSignedRank(String what, double[] x, double mu, Alternative alternative) {
        try {
            HypothesisTests.wilcoxonSignedRank(x, mu, alternative);
            fail(what + ": accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
    }

    // ------------------------------------------------- the k-sample test --

    /**
     * A SplitMix64 finalizer. Group seeds that differ by one make streams of
     * the in-test generator that differ by one too, and those are not
     * independent enough to measure a p-value with -- which is a thing this
     * test found out by measuring one.
     */
    private static long mix(long z) {
        z = z + 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /** {@code H} through the squared rank sums, the same number by algebra. */
    private static double kruskalWallisAgain(double[][] groups) {
        int total = 0;
        for (int g = 0; g < groups.length; g++) {
            total += groups[g].length;
        }
        double[] pooled = new double[total];
        int at = 0;
        for (int g = 0; g < groups.length; g++) {
            System.arraycopy(groups[g], 0, pooled, at, groups[g].length);
            at += groups[g].length;
        }
        Ranks.Result ranked = Ranks.of(pooled);
        double n = total;
        double sum = 0.0;
        at = 0;
        for (int g = 0; g < groups.length; g++) {
            double rankSum = 0.0;
            for (int i = 0; i < groups[g].length; i++) {
                rankSum += ranked.ranks[at++];
            }
            sum += rankSum * rankSum / groups[g].length;
        }
        double h = 12.0 / (n * (n + 1.0)) * sum - 3.0 * (n + 1.0);
        return h / (1.0 - ranked.tieSum / (n * n * n - n));
    }

    @Test
    public void testTheKruskalWallisStatisticIsWhatTheOtherFormGives() {
        // the mean ranks and the squared rank sums are the same statistic by
        // algebra and different arithmetic. Measured worst 5.0e-14 relative
        for (int t = 0; t < 2; t++) {
            for (long seed = 1; seed <= 200; seed++) {
                double[][] groups = { rounded(sample(new Normal(0.0, 1.0), 6, mix(seed * 31L)), t == 1),
                        rounded(sample(new Normal(0.0, 1.0), 9, mix(seed * 31L + 1)), t == 1),
                        rounded(sample(new Normal(0.0, 1.0), 14, mix(seed * 31L + 2)), t == 1),
                        rounded(sample(new Normal(0.0, 1.0), 40, mix(seed * 31L + 3)), t == 1) };
                double want = kruskalWallisAgain(groups);
                double got = HypothesisTests.kruskalWallis(groups).statistic;
                assertEquals("seed=" + seed, want, got, 1.0e-12 * Math.max(1.0, want));
            }
        }
    }

    @Test
    public void testKruskalWallisOnTwoGroupsIsTheSquareOfTheStandardizedU() {
        // the two tests are one test where they overlap, ties and all, which
        // is what neither of them can drift away from. Measured 1.1e-15
        // relative
        for (int t = 0; t < 2; t++) {
            for (long seed = 1; seed <= 300; seed++) {
                double[] x = rounded(sample(new Normal(0.0, 1.0), 11, mix(seed * 31L)), t == 1);
                double[] y = rounded(sample(new Normal(0.2, 1.0), 14, mix(seed * 31L + 1)), t == 1);
                double h = HypothesisTests.kruskalWallis(new double[][] { x, y }).statistic;
                double u = HypothesisTests.mannWhitneyU(x, y, Alternative.TWO_SIDED).statistic;
                double[] pooled = new double[x.length + y.length];
                System.arraycopy(x, 0, pooled, 0, x.length);
                System.arraycopy(y, 0, pooled, x.length, y.length);
                double tieSum = Ranks.of(pooled).tieSum;
                double m = x.length;
                double n = y.length;
                double total = m + n;
                double variance = m * n / 12.0 * ((total + 1.0) - tieSum / (total * (total - 1.0)));
                double z = (u - m * n / 2.0) / Math.sqrt(variance);
                assertEquals("seed=" + seed, z * z, h, 1.0e-12 * Math.max(1.0, z * z));
            }
        }
    }

    @Test
    public void testKruskalWallisIsBlindToTheTransformAndToTheOrderOfTheGroups() {
        for (long seed = 1; seed <= 300; seed++) {
            double[] a = sample(new Normal(0.0, 1.0), 7, mix(seed * 31L));
            double[] b = sample(new Normal(0.5, 1.0), 9, mix(seed * 31L + 1));
            double[] c = sample(new Normal(1.0, 2.0), 11, mix(seed * 31L + 2));
            TestResult plain = HypothesisTests.kruskalWallis(new double[][] { a, b, c });
            double[][] moved = { exp(a), exp(b), exp(c) };
            TestResult after = HypothesisTests.kruskalWallis(moved);
            String at = "seed=" + seed;
            assertEquals(at + ": the transform moved the statistic", plain.statistic, after.statistic, 0.0);
            assertEquals(at + ": the transform moved the p-value", plain.pValue, after.pValue, 0.0);
            TestResult swapped = HypothesisTests.kruskalWallis(new double[][] { c, a, b });
            // only the order the group sums are added up in changes, which is
            // the one thing here that is not bit for bit
            assertEquals(at + ": the order moved the statistic", plain.statistic, swapped.statistic,
                    1.0e-12 * Math.max(1.0, plain.statistic));
            assertEquals(at + ": degrees of freedom", 2.0, plain.degreesOfFreedom, 0.0);
            assertEquals(at + ": the alternative is about the statistic", Alternative.GREATER, plain.alternative);
        }
    }

    @Test
    public void testTheKruskalWallisPValueIsValidUnderTheNullHypothesis() {
        // measured over 20000 replications at the five percent level: 0.0412
        // at three groups of four, 0.0463 at three of ten and 0.0495 at three
        // of a hundred. The chi-squared null is the large-sample one, so what
        // it does at the small end is fall short rather than overspend
        int[][] shapes = { { 4, 4, 4 }, { 5, 7, 9 }, { 10, 10, 10 }, { 20, 20, 20, 20 } };
        int reps = 8000;
        for (int i = 0; i < shapes.length; i++) {
            int[] shape = shapes[i];
            int below01 = 0;
            int below05 = 0;
            int below10 = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[][] groups = new double[shape.length][];
                for (int g = 0; g < shape.length; g++) {
                    groups[g] = sample(new Normal(0.0, 1.0), shape[g], mix(seed * 31L + g));
                }
                double p = HypothesisTests.kruskalWallis(groups).pValue;
                assertTrue("p = " + p, p >= 0.0 && p <= 1.0);
                if (p < 0.01) {
                    below01++;
                }
                if (p < 0.05) {
                    below05++;
                }
                if (p < 0.10) {
                    below10++;
                }
            }
            String at = "groups of " + shape.length;
            assertTrue(at + " P(p<0.01) = " + below01 / (double) reps, below01 / (double) reps <= 0.01 + 0.005);
            assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps <= 0.05 + 0.010);
            assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps <= 0.10 + 0.018);
            if (shape[0] >= 10) {
                assertTrue(at + " P(p<0.05) = " + below05 / (double) reps, below05 / (double) reps >= 0.05 - 0.012);
                assertTrue(at + " P(p<0.10) = " + below10 / (double) reps, below10 / (double) reps >= 0.10 - 0.020);
            }
        }
    }

    @Test
    public void testKruskalWallisSeesAGroupOutOfPlace() {
        double[] low = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        double[] middle = { 6.0, 7.0, 8.0, 9.0, 10.0 };
        double[] high = { 11.0, 12.0, 13.0, 14.0, 15.0 };
        TestResult apart = HypothesisTests.kruskalWallis(new double[][] { low, middle, high });
        assertTrue("three separated groups are not significant: " + apart.pValue, apart.pValue < 0.005);
        double[] a = { 1.0, 4.0, 7.0, 10.0, 13.0 };
        double[] b = { 2.0, 5.0, 8.0, 11.0, 14.0 };
        double[] c = { 3.0, 6.0, 9.0, 12.0, 15.0 };
        TestResult mixed = HypothesisTests.kruskalWallis(new double[][] { a, b, c });
        assertTrue("interleaved groups are evidence: " + mixed.pValue, mixed.pValue > 0.5);
        // and a group of identical values is a legitimate group
        double[] flat = { 3.0, 3.0, 3.0, 3.0 };
        TestResult withFlat = HypothesisTests.kruskalWallis(new double[][] { low, flat, high });
        assertTrue("a flat group broke it", withFlat.pValue >= 0.0 && withFlat.pValue <= 1.0);
    }

    @Test
    public void testKruskalWallisRejectsWhatItCannotAnswer() {
        double[] ok = { 1.0, 2.0, 3.0 };
        refusesKruskalWallis("null", null);
        refusesKruskalWallis("one group", new double[][] { ok });
        refusesKruskalWallis("a null group", new double[][] { ok, null });
        refusesKruskalWallis("an empty group", new double[][] { ok, new double[0] });
        refusesKruskalWallis("a value that is not finite", new double[][] { ok, { 1.0, Double.NaN } });
        refusesKruskalWallis("no spread at all", new double[][] { { 2.0, 2.0 }, { 2.0 } });
    }

    private static void refusesKruskalWallis(String what, double[][] groups) {
        try {
            HypothesisTests.kruskalWallis(groups);
            fail(what + ": accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
    }

    private static double[] exp(double[] x) {
        double[] y = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            y[i] = Math.exp(x[i]);
        }
        return y;
    }

    private static double[] rounded(double[] x, boolean coarse) {
        if (!coarse) {
            return x;
        }
        double[] y = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            y[i] = Math.rint(2.0 * x[i]) / 2.0;
        }
        return y;
    }
}
