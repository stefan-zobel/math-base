package math.stats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.probe.CovarianceAccumulator;
import math.stats.rank.KendallTau;
import math.stats.rank.Ranks;
import math.stats.rank.SpearmanRho;

/**
 * A correlation is checked against a second way of computing it, and a
 * correlation test against the regression it is secretly the same as. What is
 * left over -- that the rank version cannot see a transform the other one can
 * -- is the whole reason for having both.
 */
public final class RankCorrelationTest {

    /**
     * A SplitMix64 finalizer over the in-test generator's seed. Streams started
     * at neighbouring seeds are not independent enough to measure a p-value
     * with, which the Kruskal-Wallis calibration found out by measuring one.
     */
    private static long mix(long z) {
        z = z + 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /** Deterministic standard normals, coarsened into ties on request. */
    private static double[] sample(int n, long seed, boolean rounded) {
        math.distribution.Normal d = new math.distribution.Normal(0.0, 1.0);
        long state = seed;
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            x[i] = d.inverseCdf(((state >>> 11) + 0.5) * 0x1.0p-53);
            if (rounded) {
                x[i] = Math.rint(2.0 * x[i]) / 2.0;
            }
        }
        return x;
    }

    /** Pearson through the accumulator in {@code math.probe}. */
    private static double pearson(double[] a, double[] b) {
        CovarianceAccumulator moments = new CovarianceAccumulator(2);
        double[] pair = new double[2];
        for (int i = 0; i < a.length; i++) {
            pair[0] = a[i];
            pair[1] = b[i];
            moments.add(pair);
        }
        return moments.correlation()[1];
    }

    // ------------------------------------------------------ the coefficient --

    @Test
    public void testRhoIsPearsonOnTheMidranks() {
        // the definition, checked against an independent Pearson that lives in
        // another package and shares no code with it. Measured worst 3.3e-16
        for (int t = 0; t < 2; t++) {
            for (int n : new int[] { 4, 9, 17, 60, 300 }) {
                for (long seed = 1; seed <= 40; seed++) {
                    double[] x = sample(n, mix(seed * 31L), t == 1);
                    double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                    Ranks.Result rx = Ranks.of(x);
                    Ranks.Result ry = Ranks.of(y);
                    if (allTied(rx) || allTied(ry)) {
                        continue;
                    }
                    String at = "n=" + n + " seed=" + seed + (t == 1 ? " rounded" : "");
                    assertEquals(at, pearson(rx.ranks, ry.ranks), SpearmanRho.coefficient(rx, ry), 1.0e-15);
                }
            }
        }
    }

    @Test
    public void testRhoIsTheFamiliarFormulaOnlyWhenNothingIsTied() {
        // 1 - 6 sum d^2 / (n^3 - n) is the definition's equal on untied data
        // and not otherwise, which is why the class computes the definition.
        // Measured worst 1.1e-16 where they do agree
        double worstTied = 0.0;
        for (int t = 0; t < 2; t++) {
            for (int n : new int[] { 5, 17, 60 }) {
                for (long seed = 1; seed <= 40; seed++) {
                    double[] x = sample(n, mix(seed * 31L), t == 1);
                    double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                    Ranks.Result rx = Ranks.of(x);
                    Ranks.Result ry = Ranks.of(y);
                    if (allTied(rx) || allTied(ry)) {
                        continue;
                    }
                    double s = 0.0;
                    for (int i = 0; i < n; i++) {
                        double d = rx.ranks[i] - ry.ranks[i];
                        s += d * d;
                    }
                    double simple = 1.0 - 6.0 * s / (n * (double) n * n - n);
                    double rho = SpearmanRho.coefficient(rx, ry);
                    if (rx.hasTies || ry.hasTies) {
                        worstTied = Math.max(worstTied, Math.abs(rho - simple));
                    } else {
                        assertEquals("n=" + n + " seed=" + seed, simple, rho, 1.0e-15);
                    }
                }
            }
        }
        assertTrue("with ties the two forms should part company, and they did not", worstTied > 1.0e-3);
    }

    private static boolean allTied(Ranks.Result r) {
        return r.ranks[0] == r.ranks[r.ranks.length - 1];
    }

    @Test
    public void testRhoIsOneAtAgreementAndMinusOneAtReversal() {
        double[] up = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] steeper = { 1.0, 8.0, 27.0, 64.0, 125.0, 216.0, 343.0 };
        double[] down = { 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0 };
        assertEquals("the same ordering", 1.0, SpearmanRho.coefficient(up, steeper), 0.0);
        assertEquals("the opposite ordering", -1.0, SpearmanRho.coefficient(up, down), 0.0);
        // Pearson agrees only where the relation is a straight line
        assertEquals("Pearson on the straight reversal", -1.0, pearson(up, down), 1.0e-15);
        // measured at 0.9344: the cube is monotone, so the ranks are in step,
        // but it is not a straight line and Pearson charges for the bend
        assertTrue("Pearson on the curve is short of one", pearson(up, steeper) < 0.94);
    }

    // -------------------------------------------------------- the exact null --

    @Test
    public void testTheExactNullAtThreeIsTheOneThatCanBeCountedByHand() {
        // six pairings: one gives rho = 1, two give 0.5, two give -0.5, one
        // gives -1, and rho = 0 cannot happen at all
        assertEquals("P[rho >= 1]", 1.0 / 6.0, SpearmanRho.barFExact(3, 1.0), 1.0e-15);
        assertEquals("P[rho >= 0.5]", 0.5, SpearmanRho.barFExact(3, 0.5), 1.0e-15);
        assertEquals("P[rho >= 0] is P[rho >= 0.5], since nothing lies between",
                SpearmanRho.barFExact(3, 0.5), SpearmanRho.barFExact(3, 0.0), 0.0);
        assertEquals("P[rho >= -1]", 1.0, SpearmanRho.barFExact(3, -1.0), 1.0e-15);
    }

    @Test
    public void testTheTwoTailsOfTheExactNullPartitionIt() {
        // a discrete statistic counts its own value in both tails. Measured
        // worst departure 4.4e-16
        for (int n = 3; n <= 9; n++) {
            double span = n * (double) n * n - n;
            int atoms = (n * n * n - n) / 6 + 1;
            for (int k = 0; k < atoms; k++) {
                double rho = 1.0 - 12.0 * k / span;
                double upper = SpearmanRho.barFExact(n, rho);
                double lower = SpearmanRho.barFExact(n, -rho);
                double previous = (k == 0) ? 0.0 : SpearmanRho.barFExact(n, 1.0 - 12.0 * (k - 1) / span);
                assertEquals("n=" + n + " k=" + k, 1.0 + (upper - previous), lower + upper, 1.0e-15);
            }
        }
    }

    @Test
    public void testTheExactNullRefusesWhatItCannotAfford() {
        // nine is the measured limit: 7.2 ms against 67 ms at ten
        assertTrue("nine is inside", SpearmanRho.barFExact(9, 0.7) > 0.0);
        try {
            SpearmanRho.barFExact(10, 0.7);
            fail("a sample above the limit was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message names the limit", expected.getMessage().contains("9"));
        }
        try {
            SpearmanRho.barFExact(1, 0.0);
            fail("a single observation was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        try {
            SpearmanRho.barFAsymptotic(2, 0.0);
            fail("two observations leave no degrees of freedom");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        // a perfect agreement has no tail beyond it and a perfect reversal has
        // all of it, and no finite t says so
        assertEquals(0.0, SpearmanRho.barFAsymptotic(9, 1.0), 0.0);
        assertEquals(1.0, SpearmanRho.barFAsymptotic(9, -1.0), 0.0);
    }

    @Test
    public void testTheApproximationClosesOnTheExactNull() {
        // measured 3.3e-2 at n = 5 and 4.6e-3 at n = 9, over the p-values a
        // reader acts on. It is a poor approximation down here, which is what
        // the exact path is for
        int[] sizes = { 5, 6, 7, 8, 9 };
        double[] allowed = { 3.5e-2, 1.8e-2, 1.1e-2, 7.5e-3, 5.0e-3 };
        for (int i = 0; i < sizes.length; i++) {
            int n = sizes[i];
            double span = n * (double) n * n - n;
            int atoms = (n * n * n - n) / 6 + 1;
            double worst = 0.0;
            for (int k = 0; k < atoms; k++) {
                double rho = 1.0 - 12.0 * k / span;
                double exact = SpearmanRho.barFExact(n, rho);
                if (exact < 0.001 || exact > 0.2) {
                    continue;
                }
                worst = Math.max(worst, Math.abs(SpearmanRho.barFAsymptotic(n, rho) - exact));
            }
            assertTrue("n=" + n + ": worst difference " + worst, worst > 0.0 && worst <= allowed[i]);
        }
    }

    // --------------------------------------------------------- the two tests --

    @Test
    public void testThePearsonTestIsTheSlopeOfASimpleRegression() {
        // the same number by a least squares solve that shares no code with
        // it: R^2 is r^2, and the t of the slope is the t of r. Measured 8.9e-16
        // on the first, 4.0e-15 relative on the second and 6.7e-16 on the
        // p-value
        for (long seed = 1; seed <= 200; seed++) {
            int n = 5 + (int) (seed % 40);
            double[] x = sample(n, mix(seed * 31L), false);
            double[] noise = sample(n, mix(seed * 31L + 1), false);
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = 0.5 + 0.8 * x[i] + noise[i];
            }
            DMatrix design = new DMatrix(n, 2);
            DMatrix response = new DMatrix(n, 1);
            for (int i = 0; i < n; i++) {
                design.set(i, 0, 1.0);
                design.set(i, 1, x[i]);
                response.set(i, 0, y[i]);
            }
            LSSummary fit = OLS.estimate(0.05, design, response);
            TestResult r = HypothesisTests.pearsonCorrelation(x, y, Alternative.TWO_SIDED);
            String at = "seed=" + seed + " n=" + n;
            assertEquals(at + ": R^2", fit.getRSquared(), r.statistic * r.statistic, 1.0e-14);
            assertEquals(at + ": p-value", fit.getPValues().get(1), r.pValue, 1.0e-14);
            assertEquals(at + ": degrees of freedom", n - 2.0, r.degreesOfFreedom, 0.0);
            double t = r.statistic * Math.sqrt((n - 2.0) / (1.0 - r.statistic * r.statistic));
            assertEquals(at + ": t", fit.getTValues().get(1), t, 1.0e-13 * Math.max(1.0, Math.abs(t)));
        }
    }

    @Test
    public void testSpearmanIsBlindToAnIncreasingTransformAndPearsonIsNot() {
        // bit for bit on one side, and demonstrably not on the other -- which
        // is the entire reason both tests are here
        double moved = 0.0;
        for (long seed = 1; seed <= 200; seed++) {
            double[] x = sample(13, mix(seed * 31L), false);
            double[] y = sample(13, mix(seed * 31L + 1), false);
            double[] ex = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                ex[i] = Math.exp(x[i]);
            }
            for (Alternative alternative : Alternative.values()) {
                TestResult before = HypothesisTests.spearmanCorrelation(x, y, alternative);
                TestResult after = HypothesisTests.spearmanCorrelation(ex, y, alternative);
                String at = "seed=" + seed + " " + alternative;
                assertEquals(at + ": statistic", before.statistic, after.statistic, 0.0);
                assertEquals(at + ": p-value", before.pValue, after.pValue, 0.0);
                assertEquals(at + ": path", before.test, after.test);
            }
            moved = Math.max(moved, Math.abs(HypothesisTests.pearsonCorrelation(x, y, Alternative.TWO_SIDED).statistic
                    - HypothesisTests.pearsonCorrelation(ex, y, Alternative.TWO_SIDED).statistic));
        }
        assertTrue("Pearson should have moved and did not: " + moved, moved > 0.1);
    }

    @Test
    public void testTheTwoSidedValueIsTwiceTheSmallerTail() {
        for (int t = 0; t < 2; t++) {
            for (long seed = 1; seed <= 300; seed++) {
                int n = (t == 0) ? 8 : 25;
                double[] x = sample(n, mix(seed * 31L), t == 1);
                double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                double less = HypothesisTests.spearmanCorrelation(x, y, Alternative.LESS).pValue;
                double greater = HypothesisTests.spearmanCorrelation(x, y, Alternative.GREATER).pValue;
                double two = HypothesisTests.spearmanCorrelation(x, y, Alternative.TWO_SIDED).pValue;
                assertEquals("seed=" + seed, Math.min(1.0, 2.0 * Math.min(less, greater)), two, 0.0);
                double pLess = HypothesisTests.pearsonCorrelation(x, y, Alternative.LESS).pValue;
                double pGreater = HypothesisTests.pearsonCorrelation(x, y, Alternative.GREATER).pValue;
                assertEquals("seed=" + seed + ": the Pearson tails partition", 1.0, pLess + pGreater, 1.0e-14);
            }
        }
    }

    @Test
    public void testATieOrASampleTooLargeSendsSpearmanDownTheOtherPath() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        double[] y = { 2.0, 1.0, 5.0, 3.0, 4.0 };
        assertEquals("no ties, small sample", "Spearman rank correlation, exact",
                HypothesisTests.spearmanCorrelation(x, y, Alternative.TWO_SIDED).test);
        double[] tied = { 2.0, 1.0, 5.0, 3.0, 3.0 };
        assertEquals("one tie is enough", "Spearman rank correlation, asymptotic",
                HypothesisTests.spearmanCorrelation(x, tied, Alternative.TWO_SIDED).test);
        double[] nine = sample(9, mix(11L), false);
        double[] alsoNine = sample(9, mix(12L), false);
        assertEquals("nine is inside the limit", "Spearman rank correlation, exact",
                HypothesisTests.spearmanCorrelation(nine, alsoNine, Alternative.TWO_SIDED).test);
        double[] ten = sample(10, mix(11L), false);
        double[] alsoTen = sample(10, mix(12L), false);
        TestResult past = HypothesisTests.spearmanCorrelation(ten, alsoTen, Alternative.TWO_SIDED);
        assertEquals("ten is past it", "Spearman rank correlation, asymptotic", past.test);
        assertEquals("and then the degrees of freedom are real", 8.0, past.degreesOfFreedom, 0.0);
    }

    @Test
    public void testThePValuesAreValidUnderTheNullHypothesis() {
        // measured over 20000 replications at the five percent level: Spearman
        // spends 0.0161 at n = 5, where the null has too few atoms to spend it,
        // 0.0475 at n = 7 and 0.0510 at n = 200; Pearson stays between 0.0470
        // and 0.0514 throughout
        int[] sizes = { 5, 7, 12, 30, 200 };
        int reps = 8000;
        for (int i = 0; i < sizes.length; i++) {
            int n = sizes[i];
            int[] rho = new int[3];
            int[] r = new int[3];
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(n, mix(seed * 31L), false);
                double[] y = sample(n, mix(seed * 31L + 1), false);
                count(rho, HypothesisTests.spearmanCorrelation(x, y, Alternative.TWO_SIDED).pValue);
                count(r, HypothesisTests.pearsonCorrelation(x, y, Alternative.TWO_SIDED).pValue);
            }
            String at = "n=" + n;
            // the Monte Carlo standard error at 8000 replications is 0.0011 at
            // the one percent level and 0.0024 at the five
            assertRate(at + " Spearman P(p<0.01)", rho[0] / (double) reps, 0.01 + 0.005);
            assertRate(at + " Spearman P(p<0.05)", rho[1] / (double) reps, 0.05 + 0.009);
            assertRate(at + " Spearman P(p<0.10)", rho[2] / (double) reps, 0.10 + 0.015);
            assertRate(at + " Pearson P(p<0.01)", r[0] / (double) reps, 0.01 + 0.005);
            assertRate(at + " Pearson P(p<0.05)", r[1] / (double) reps, 0.05 + 0.009);
            assertRate(at + " Pearson P(p<0.10)", r[2] / (double) reps, 0.10 + 0.015);
            if (n >= 12) {
                // and a null fine enough to spend the level does spend it
                assertTrue(at + " Spearman P(p<0.05) = " + rho[1] / (double) reps,
                        rho[1] / (double) reps >= 0.05 - 0.012);
                assertTrue(at + " Pearson P(p<0.05) = " + r[1] / (double) reps,
                        r[1] / (double) reps >= 0.05 - 0.012);
            }
        }
    }

    @Test
    public void testTiesAtASmallSampleMakeSpearmanOverspend() {
        // recorded rather than forbidden: with ties there is no exact path, and
        // the approximation errs towards rejecting. Measured over 20000
        // replications of a true null at the five percent level: 0.0575 at
        // n = 5, 0.0553 at n = 7, and back to 0.0505 by n = 12. The bound here
        // is on how bad it may get, not on its being bad
        int reps = 8000;
        for (int n : new int[] { 5, 7, 12 }) {
            int[] rho = new int[3];
            int dropped = 0;
            for (long seed = 1; seed <= reps; seed++) {
                double[] x = sample(n, mix(seed * 31L), true);
                double[] y = sample(n, mix(seed * 31L + 1), true);
                try {
                    count(rho, HypothesisTests.spearmanCorrelation(x, y, Alternative.TWO_SIDED).pValue);
                } catch (IllegalArgumentException everythingTied) {
                    dropped++;
                }
            }
            double kept = reps - dropped;
            assertTrue("n=" + n + ": too many replications were all tied", kept > 0.99 * reps);
            assertRate("n=" + n + " tied P(p<0.05)", rho[1] / kept, 0.070);
            assertRate("n=" + n + " tied P(p<0.01)", rho[0] / kept, 0.022);
        }
    }

    private static void assertRate(String what, double got, double atMost) {
        assertTrue(what + " = " + got + ", wanted at most " + atMost, got <= atMost);
    }

    private static void count(int[] below, double p) {
        assertTrue("p = " + p, p >= 0.0 && p <= 1.0);
        if (p < 0.01) {
            below[0]++;
        }
        if (p < 0.05) {
            below[1]++;
        }
        if (p < 0.10) {
            below[2]++;
        }
    }

    @Test
    public void testABendInTheRelationIsWhereTheTwoTestsPartCompany() {
        // a monotone curve steep enough that Pearson loses most of it while the
        // ranks do not notice the bend at all
        int n = 30;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i / (double) n;
            y[i] = Math.exp(9.0 * x[i]);
        }
        TestResult rho = HypothesisTests.spearmanCorrelation(x, y, Alternative.GREATER);
        TestResult tau = HypothesisTests.kendallTau(x, y, Alternative.GREATER);
        TestResult r = HypothesisTests.pearsonCorrelation(x, y, Alternative.GREATER);
        assertEquals("the ranks are in step, so tau is one too", 1.0, tau.statistic, 0.0);
        assertEquals("and nothing beats a perfect ordering there either", 1.0 / factorial(30), tau.pValue,
                1.0e-45);
        assertEquals("the ranks are in step, so rho is one", 1.0, rho.statistic, 0.0);
        // a coefficient of one leaves no tail beyond it at all, and no finite t
        // would have said so
        assertEquals("and nothing can be more extreme than a perfect ordering", 0.0, rho.pValue, 0.0);
        // measured: r = 0.7225 with a p-value of 3.3e-6, against a rho of one.
        // Both still reject; the difference is how much of the association
        // each of them is able to see
        assertTrue("Pearson should be well short of one, and is " + r.statistic, r.statistic < 0.8);
        assertTrue("Pearson should still see it: " + r.pValue, r.pValue < 1.0e-5);
    }

    @Test
    public void testTheCorrelationTestsRejectWhatTheyCannotAnswer() {
        double[] ok = { 1.0, 2.0, 3.0, 4.0 };
        refuses("x null", null, ok);
        refuses("y null", ok, null);
        refuses("different lengths", ok, new double[] { 1.0, 2.0, 3.0 });
        refuses("too short", new double[] { 1.0, 2.0 }, new double[] { 3.0, 4.0 });
        refuses("not finite", new double[] { 1.0, 2.0, Double.NaN, 4.0 }, ok);
        refuses("x has no spread", new double[] { 2.0, 2.0, 2.0, 2.0 }, ok);
        refuses("y has no spread", ok, new double[] { 7.0, 7.0, 7.0, 7.0 });
        for (int which = 0; which < 2; which++) {
            try {
                if (which == 0) {
                    HypothesisTests.spearmanCorrelation(ok, ok, null);
                } else {
                    HypothesisTests.pearsonCorrelation(ok, ok, null);
                }
                fail("a null alternative was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue("empty message", expected.getMessage().length() > 0);
            }
        }
    }

    private static void refuses(String what, double[] x, double[] y) {
        try {
            HypothesisTests.spearmanCorrelation(x, y, Alternative.TWO_SIDED);
            fail(what + ": Spearman accepted it");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
        try {
            HypothesisTests.pearsonCorrelation(x, y, Alternative.TWO_SIDED);
            fail(what + ": Pearson accepted it");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
        try {
            HypothesisTests.kendallTau(x, y, Alternative.TWO_SIDED);
            fail(what + ": Kendall accepted it");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": empty message", expected.getMessage().length() > 0);
        }
    }

    // --------------------------------------------------------- Kendall tau --

    /** tau-b the obvious way: every pair of observations, one at a time. */
    private static double[] naiveTau(double[] x, double[] y) {
        int n = x.length;
        long concordant = 0L;
        long discordant = 0L;
        long tiedInX = 0L;
        long tiedInY = 0L;
        for (int i = 0; i < n; i++) {
            for (int k = i + 1; k < n; k++) {
                double dx = Math.signum(x[i] - x[k]);
                double dy = Math.signum(y[i] - y[k]);
                if (dx == 0.0) {
                    tiedInX++;
                }
                if (dy == 0.0) {
                    tiedInY++;
                }
                if (dx != 0.0 && dy != 0.0) {
                    if (dx * dy > 0.0) {
                        concordant++;
                    } else {
                        discordant++;
                    }
                }
            }
        }
        double total = n * (n - 1) / 2.0;
        double s = concordant - discordant;
        return new double[] { s / Math.sqrt((total - tiedInX) * (total - tiedInY)), s };
    }

    @Test
    public void testKnightsAlgorithmIsTheQuadraticCount() {
        // the O(n log n) sort-and-count against the O(n^2) loop over pairs it
        // exists to avoid -- bit for bit, with ties and without
        for (int t = 0; t < 2; t++) {
            for (int n : new int[] { 3, 9, 17, 60, 300 }) {
                for (long seed = 1; seed <= 40; seed++) {
                    double[] x = sample(n, mix(seed * 31L), t == 1);
                    double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                    KendallTau.Result got;
                    try {
                        got = KendallTau.of(x, y);
                    } catch (IllegalArgumentException everythingTied) {
                        continue;
                    }
                    double[] want = naiveTau(x, y);
                    String at = "n=" + n + " seed=" + seed + (t == 1 ? " rounded" : "");
                    assertEquals(at + ": tau", want[0], got.tau, 0.0);
                    assertEquals(at + ": S", want[1], got.s, 0.0);
                }
            }
        }
    }

    /** P[S >= s] counted over every pairing, for a small n. */
    private static double[] bruteKendallTails(int n) {
        int maxK = n * (n - 1) / 2;
        double[] counts = new double[maxK + 1];
        int[] p = new int[n];
        for (int i = 0; i < n; i++) {
            p[i] = i;
        }
        permute(p, 0, counts);
        double total = 0.0;
        for (int k = 0; k <= maxK; k++) {
            total += counts[k];
        }
        double[] tails = new double[maxK + 1];
        double reached = 0.0;
        for (int k = 0; k <= maxK; k++) {
            reached += counts[k];
            tails[k] = reached / total;
        }
        return tails;
    }

    private static void permute(int[] p, int at, double[] counts) {
        if (at == p.length) {
            int inversions = 0;
            for (int i = 0; i < p.length; i++) {
                for (int k = i + 1; k < p.length; k++) {
                    if (p[i] > p[k]) {
                        inversions++;
                    }
                }
            }
            counts[inversions] += 1.0;
            return;
        }
        for (int i = at; i < p.length; i++) {
            int t = p[at];
            p[at] = p[i];
            p[i] = t;
            permute(p, at + 1, counts);
            t = p[at];
            p[at] = p[i];
            p[i] = t;
        }
    }

    @Test
    public void testTheExactKendallNullIsWhatCountingThePairingsGives() {
        // the independent algorithm: every pairing, its inversions counted the
        // slow way. Measured worst departure 1.1e-16
        for (int n = 2; n <= 8; n++) {
            double[] want = bruteKendallTails(n);
            double total = n * (n - 1) / 2.0;
            for (int k = 0; k < want.length; k++) {
                double s = total - 2.0 * k;
                assertEquals("n=" + n + " s=" + s, want[k], KendallTau.barFExact(n, s), 1.0e-15);
            }
        }
    }

    @Test
    public void testTheExtremeAtomOfTheKendallNullIsOneOverNFactorial() {
        // the sharp end of the recursion, and the one place a carried window
        // sum could have destroyed it: only one pairing of n observations
        // agrees perfectly, so its probability is 1 / n!. It is summed from
        // this end, where nothing has been subtracted yet
        for (int n : new int[] { 5, 10, 20, 50, 100 }) {
            double total = n * (n - 1) / 2.0;
            double want = 1.0 / factorial(n);
            double got = KendallTau.barFExact(n, total);
            assertEquals("n=" + n, want, got, 1.0e-13 * want);
        }
    }

    private static double factorial(int n) {
        double f = 1.0;
        for (int i = 2; i <= n; i++) {
            f *= i;
        }
        return f;
    }

    @Test
    public void testTheTwoTailsOfTheKendallNullPartitionIt() {
        // measured worst departure 1.1e-14
        for (int n = 2; n <= 20; n++) {
            double total = n * (n - 1) / 2.0;
            for (double s = -total; s <= total; s += 2.0) {
                double upper = KendallTau.barFExact(n, s);
                double lower = KendallTau.barFExact(n, -s);
                double atom = upper - KendallTau.barFExact(n, s + 2.0);
                assertEquals("n=" + n + " s=" + s, 1.0 + atom, lower + upper, 5.0e-14);
            }
        }
    }

    @Test
    public void testTheKendallApproximationClosesOnTheExactNull() {
        // measured 6.3e-3 at n = 5 and 4.9e-4 at n = 100, over the p-values a
        // reader acts on
        int[] sizes = { 5, 10, 20, 50, 100 };
        double[] allowed = { 7.0e-3, 5.5e-3, 2.8e-3, 1.1e-3, 5.5e-4 };
        for (int i = 0; i < sizes.length; i++) {
            int n = sizes[i];
            double total = n * (n - 1) / 2.0;
            double[] inOrder = new double[n];
            for (int k = 0; k < n; k++) {
                inOrder[k] = k;
            }
            KendallTau.Result untied = KendallTau.of(inOrder, inOrder);
            double worst = 0.0;
            for (double s = total; s >= 0.0; s -= 2.0) {
                double exact = KendallTau.barFExact(n, s);
                if (exact < 0.001) {
                    continue;
                }
                if (exact > 0.2) {
                    break;
                }
                worst = Math.max(worst, Math.abs(KendallTau.barFAsymptotic(untied, s) - exact));
            }
            assertTrue("n=" + n + ": worst difference " + worst, worst > 0.0 && worst <= allowed[i]);
        }
    }

    @Test
    public void testTheKendallNullRefusesWhatItCannotAfford() {
        // two hundred is the measured limit: 11 ms for the recursion there
        assertTrue("two hundred is inside", KendallTau.barFExact(200, 2000.0) > 0.0);
        try {
            KendallTau.barFExact(201, 2000.0);
            fail("a sample above the limit was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message names the limit", expected.getMessage().contains("200"));
        }
        try {
            KendallTau.of(new double[] { 1.0, 2.0 }, new double[] { 1.0, Double.NaN });
            fail("a NaN was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
        try {
            KendallTau.of(new double[] { 1.0, 1.0, 1.0 }, new double[] { 1.0, 2.0, 3.0 });
            fail("a sample with no untied pair was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("empty message", expected.getMessage().length() > 0);
        }
    }

    @Test
    public void testTauIsOneAtAgreementAndMinusOneAtReversal() {
        double[] up = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
        double[] cubed = { 1.0, 8.0, 27.0, 64.0, 125.0, 216.0, 343.0 };
        double[] down = { 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0 };
        assertEquals("the same ordering", 1.0, KendallTau.of(up, cubed).tau, 0.0);
        assertEquals("the opposite ordering", -1.0, KendallTau.of(up, down).tau, 0.0);
        // and a tie anywhere puts it out of reach, which is what the b in
        // tau-b is about
        double[] oneTie = { 1.0, 2.0, 2.0, 4.0, 5.0, 6.0, 7.0 };
        assertTrue("a single tie should cost it", KendallTau.of(oneTie, up).tau < 1.0);
    }

    @Test
    public void testKendallIsBlindToAnIncreasingTransform() {
        for (long seed = 1; seed <= 300; seed++) {
            for (int n : new int[] { 9, 25 }) {
                double[] x = sample(n, mix(seed * 31L), false);
                double[] y = sample(n, mix(seed * 31L + 1), false);
                double[] ex = new double[n];
                for (int i = 0; i < n; i++) {
                    ex[i] = Math.exp(x[i]);
                }
                for (Alternative alternative : Alternative.values()) {
                    TestResult before = HypothesisTests.kendallTau(x, y, alternative);
                    TestResult after = HypothesisTests.kendallTau(ex, y, alternative);
                    String at = "n=" + n + " seed=" + seed + " " + alternative;
                    assertEquals(at + ": statistic", before.statistic, after.statistic, 0.0);
                    assertEquals(at + ": p-value", before.pValue, after.pValue, 0.0);
                }
                double less = HypothesisTests.kendallTau(x, y, Alternative.LESS).pValue;
                double greater = HypothesisTests.kendallTau(x, y, Alternative.GREATER).pValue;
                double two = HypothesisTests.kendallTau(x, y, Alternative.TWO_SIDED).pValue;
                assertEquals("n=" + n + " seed=" + seed, Math.min(1.0, 2.0 * Math.min(less, greater)), two, 0.0);
            }
        }
    }

    @Test
    public void testTheTwoCoefficientsAreTiedTogetherByDanielsInequality() {
        // |3 tau - 2 rho| <= 1 always, and the bound is reached -- measured at
        // exactly 1.0. That tau is the smaller of the two is the usual case and
        // not a rule: over 3997 samples |tau| exceeded |rho| in 341 of them,
        // most starkly where rho came out at exactly zero and tau did not
        double worst = 0.0;
        int exceeded = 0;
        int checked = 0;
        for (int t = 0; t < 2; t++) {
            for (int n : new int[] { 4, 5, 8, 20, 100 }) {
                for (long seed = 1; seed <= 100; seed++) {
                    double[] x = sample(n, mix(seed * 31L), t == 1);
                    double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                    double tau;
                    double rho;
                    try {
                        tau = KendallTau.of(x, y).tau;
                        rho = SpearmanRho.coefficient(x, y);
                    } catch (IllegalArgumentException everythingTied) {
                        continue;
                    }
                    checked++;
                    worst = Math.max(worst, Math.abs(3.0 * tau - 2.0 * rho));
                    if (Math.abs(tau) > Math.abs(rho) + 1.0e-12) {
                        exceeded++;
                    }
                }
            }
        }
        assertTrue("Daniels' inequality was broken: " + worst, worst <= 1.0 + 1.0e-12);
        assertTrue("the bound should be reached and was not: " + worst, worst > 0.99);
        assertTrue("checked too few pairs", checked > 500);
        assertTrue("tau should be the smaller one most of the time", exceeded < checked / 4);
    }

    @Test
    public void testATieOrASampleTooLargeSendsKendallDownTheOtherPath() {
        double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        double[] y = { 2.0, 1.0, 5.0, 3.0, 4.0 };
        assertEquals("no ties, small sample", "Kendall tau-b, exact",
                HypothesisTests.kendallTau(x, y, Alternative.TWO_SIDED).test);
        double[] tied = { 2.0, 1.0, 5.0, 3.0, 3.0 };
        assertEquals("one tie is enough", "Kendall tau-b, asymptotic",
                HypothesisTests.kendallTau(x, tied, Alternative.TWO_SIDED).test);
        // unlike Spearman, which stops at nine, this one reaches two hundred
        double[] wide = sample(200, mix(11L), false);
        double[] alsoWide = sample(200, mix(12L), false);
        assertEquals("two hundred is inside the limit", "Kendall tau-b, exact",
                HypothesisTests.kendallTau(wide, alsoWide, Alternative.TWO_SIDED).test);
        double[] wider = sample(201, mix(11L), false);
        double[] alsoWider = sample(201, mix(12L), false);
        assertEquals("and two hundred and one is not", "Kendall tau-b, asymptotic",
                HypothesisTests.kendallTau(wider, alsoWider, Alternative.TWO_SIDED).test);
    }

    @Test
    public void testTheKendallPValueIsValidUnderTheNullHypothesis() {
        // measured over 20000 replications at the five percent level: 0.0161 at
        // n = 5, where the null has too few atoms to spend it, 0.0461 at n = 9
        // and 0.0453 at n = 20. With ties -- where Spearman overspends -- this
        // one stays at or under the level, because tau-b corrects its variance
        // for them instead of ignoring them: 0.0138, 0.0409 and 0.0461
        int[] sizes = { 5, 9, 20, 400 };
        int reps = 8000;
        for (int t = 0; t < 2; t++) {
            for (int i = 0; i < sizes.length; i++) {
                int n = sizes[i];
                int[] below = new int[3];
                int dropped = 0;
                for (long seed = 1; seed <= reps; seed++) {
                    double[] x = sample(n, mix(seed * 31L), t == 1);
                    double[] y = sample(n, mix(seed * 31L + 1), t == 1);
                    try {
                        count(below, HypothesisTests.kendallTau(x, y, Alternative.TWO_SIDED).pValue);
                    } catch (IllegalArgumentException everythingTied) {
                        dropped++;
                    }
                }
                double kept = reps - dropped;
                String at = "n=" + n + (t == 1 ? " rounded" : "");
                assertTrue(at + ": too many replications were all tied", kept > 0.99 * reps);
                assertRate(at + " P(p<0.01)", below[0] / kept, 0.01 + 0.005);
                assertRate(at + " P(p<0.05)", below[1] / kept, 0.05 + 0.009);
                assertRate(at + " P(p<0.10)", below[2] / kept, 0.10 + 0.015);
                if (n >= 20) {
                    assertTrue(at + " P(p<0.05) = " + below[1] / kept, below[1] / kept >= 0.05 - 0.012);
                }
            }
        }
    }
}
