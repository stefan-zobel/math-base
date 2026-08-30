package math.stats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.distribution.Normal;

/**
 * A correction for multiple testing is not checked by remembering what it
 * returns. It is checked by asking whether it rejects the same hypotheses the
 * procedure it implements would reject, and whether the error rate it claims to
 * control is the one it controls.
 */
public final class MultipleTestingTest {

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

    private static final Normal STANDARD = new Normal(0.0, 1.0);

    /**
     * {@code 1 + 1/2 + ... + 1/m}, written out here so that the test does not
     * borrow the sum it is checking.
     */
    private static double harmonic(int m) {
        double sum = 0.0;
        for (int j = m; j >= 1; j--) {
            sum += 1.0 / j;
        }
        return sum;
    }

    /**
     * A family of p-values with ties, zeros and ones in it, since those are
     * where an ordering argument is easiest to get wrong.
     */
    private static double[] awkwardFamily(int m, Lcg lcg) {
        double[] p = new double[m];
        for (int i = 0; i < m; i++) {
            double u = lcg.next();
            if (u < 0.08) {
                p[i] = 0.0;
            } else if (u < 0.14) {
                p[i] = 1.0;
            } else if (u < 0.40) {
                // rounded to twentieths, so ties are common
                p[i] = Math.floor(u * 20.0) / 20.0;
            } else {
                p[i] = u;
            }
        }
        return p;
    }

    // ----------------------------------------------- the two that matter --

    @Test
    public void testTheAdjustedPValuesRejectExactlyWhatTheStepUpProcedureDoes() {
        // Benjamini and Hochberg defined a decision, not an adjusted p-value:
        // the largest k with p_(k) <= k alpha / m, then the k smallest. That
        // loop is written out here and shares nothing with the implementation
        // but its input. If the two ever disagree, the promise in the javadoc
        // -- that thresholding the adjusted values is the same thing -- is
        // false, and a caller who believed it has the wrong answer
        Lcg lcg = new Lcg(12345L);
        int checked = 0;
        for (int round = 0; round < 1500; round++) {
            int m = 1 + (int) (lcg.next() * 30);
            double[] p = awkwardFamily(m, lcg);
            double[] adjusted = MultipleTesting.benjaminiHochberg(p);
            for (double alpha : new double[] { 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.9 }) {
                int byDecision = stepUpRejectionCount(p, alpha);
                int byAdjusted = 0;
                for (int i = 0; i < m; i++) {
                    if (adjusted[i] <= alpha) {
                        byAdjusted++;
                    }
                }
                assertEquals("m=" + m + " round=" + round + " alpha=" + alpha, byDecision, byAdjusted);
                checked++;
            }
        }
        // measured over 28000 pairs in the probe: none disagreed
        assertTrue("nothing was checked", checked > 10000);
    }

    /** The original decision: how many of the ordered p-values are rejected. */
    private static int stepUpRejectionCount(double[] p, double alpha) {
        double[] ascending = p.clone();
        Arrays.sort(ascending);
        int largest = 0;
        for (int k = 1; k <= ascending.length; k++) {
            if (ascending[k - 1] <= k * alpha / ascending.length) {
                largest = k;
            }
        }
        return largest;
    }

    @Test
    public void testTheFalseDiscoveryRateIsTheOneThatIsControlled() {
        // the reason the procedure exists, measured rather than assumed. Of m
        // hypotheses m0 are true nulls; the false discovery proportion is the
        // share of the rejections that hit one of those. Its average has to sit
        // at or below alpha m0 / m
        int m = 200;
        double alpha = 0.05;
        int reps = 2000;
        for (int m0 : new int[] { 200, 190, 160, 100 }) {
            double bound = alpha * m0 / (double) m;
            double benjaminiHochberg = 0.0;
            double benjaminiYekutieli = 0.0;
            double uncorrected = 0.0;
            double found = 0.0;
            Lcg lcg = new Lcg(999L);
            for (int r = 0; r < reps; r++) {
                double[] p = mixedFamily(m, m0, 3.0, lcg);
                double[] bh = MultipleTesting.benjaminiHochberg(p);
                benjaminiHochberg += falseDiscoveryProportion(bh, m0, alpha);
                benjaminiYekutieli += falseDiscoveryProportion(
                        MultipleTesting.benjaminiYekutieli(p), m0, alpha);
                uncorrected += falseDiscoveryProportion(p, m0, alpha);
                found += trueDiscoveries(bh, m0, alpha);
            }
            String at = "m0=" + m0;
            // measured: 0.0530 / 0.0476 / 0.0394 / 0.0252 against bounds of
            // 0.0500 / 0.0475 / 0.0400 / 0.0250. The first sits 0.6 standard
            // errors above its bound, which is what 2000 replications buy
            assertTrue(at + ": Benjamini-Hochberg gave " + benjaminiHochberg / reps + " against " + bound,
                    benjaminiHochberg / reps <= bound + 0.012);
            // measured: 0.0042 .. 0.0078, far under, which is what the
            // harmonic factor is paid for
            assertTrue(at + ": Benjamini-Yekutieli gave " + benjaminiYekutieli / reps,
                    benjaminiYekutieli / reps <= bound);
            // and the number that says why any of this exists. Measured:
            // 1.0000 when every null is true, then 0.4967, 0.1754, 0.0511
            assertTrue(at + ": uncorrected testing gave " + uncorrected / reps,
                    uncorrected / reps > bound);

            if (m0 < m) {
                // a procedure that rejects nothing controls every error rate
                // there is, so the bound above only means something beside
                // this. Measured: 0.509, 0.713, 0.836 of the real effects
                double share = found / reps / (m - m0);
                assertTrue(at + ": only " + share + " of the real effects were found", share > 0.4);
            }
        }
    }

    /**
     * {@code m} one-sided p-values, the first {@code m0} from a true null and
     * the rest from a real effect.
     */
    private static double[] mixedFamily(int m, int m0, double shift, Lcg lcg) {
        double[] p = new double[m];
        for (int i = 0; i < m; i++) {
            double z = STANDARD.inverseCdf(lcg.next()) + (i < m0 ? 0.0 : shift);
            // the upper tail, so a real effect pushes the p-value down
            p[i] = STANDARD.cdf(-z);
        }
        return p;
    }

    private static double falseDiscoveryProportion(double[] adjusted, int m0, double alpha) {
        int rejected = 0;
        int falsely = 0;
        for (int i = 0; i < adjusted.length; i++) {
            if (adjusted[i] <= alpha) {
                rejected++;
                if (i < m0) {
                    falsely++;
                }
            }
        }
        return rejected == 0 ? 0.0 : falsely / (double) rejected;
    }

    private static double trueDiscoveries(double[] adjusted, int m0, double alpha) {
        int found = 0;
        for (int i = m0; i < adjusted.length; i++) {
            if (adjusted[i] <= alpha) {
                found++;
            }
        }
        return found;
    }

    // ------------------------------------------------------ the invariants --

    @Test
    public void testTheAdjustedPValuesAreMonotoneAndTiesAgree() {
        Lcg lcg = new Lcg(7L);
        for (int round = 0; round < 800; round++) {
            int m = 2 + (int) (lcg.next() * 40);
            double[] p = awkwardFamily(m, lcg);
            double[] bh = MultipleTesting.benjaminiHochberg(p);
            double[] by = MultipleTesting.benjaminiYekutieli(p);
            for (int i = 0; i < m; i++) {
                for (int k = 0; k < m; k++) {
                    String at = "round=" + round + " i=" + i + " k=" + k;
                    if (p[i] < p[k]) {
                        assertTrue(at + ": " + bh[i] + " > " + bh[k], bh[i] <= bh[k]);
                        assertTrue(at + ": " + by[i] + " > " + by[k], by[i] <= by[k]);
                    }
                    if (p[i] == p[k]) {
                        // the running minimum gives this without a special
                        // case, which is why there is no branch for it
                        assertEquals(at + ": tied p-values", bh[i], bh[k], 0.0);
                        assertEquals(at + ": tied p-values", by[i], by[k], 0.0);
                    }
                }
            }
        }
    }

    @Test
    public void testTheAdjustedPValuesStayBetweenTheirOwnAndBonferroni() {
        Lcg lcg = new Lcg(11L);
        for (int round = 0; round < 800; round++) {
            int m = 1 + (int) (lcg.next() * 40);
            double[] p = awkwardFamily(m, lcg);
            double[] bh = MultipleTesting.benjaminiHochberg(p);
            double[] by = MultipleTesting.benjaminiYekutieli(p);
            for (int i = 0; i < m; i++) {
                String at = "round=" + round + " i=" + i + " p=" + p[i];
                // never smaller than the p-value it came from. That holds by
                // construction and, without the repair in the sweep, would
                // fail by one unit in the last place
                assertTrue(at + ": " + bh[i] + " is below it", bh[i] >= p[i]);
                assertTrue(at + ": " + by[i] + " is below it", by[i] >= p[i]);
                // and never larger than Bonferroni, which corrects the same
                // family for the stricter guarantee
                assertTrue(at + ": " + bh[i] + " is above Bonferroni", bh[i] <= Math.min(1.0, m * p[i]));
                assertTrue(at + ": " + bh[i] + " left [0, 1]", bh[i] >= 0.0 && bh[i] <= 1.0);
                assertTrue(at + ": " + by[i] + " left [0, 1]", by[i] >= 0.0 && by[i] <= 1.0);
                // the price of dependence is paid in one direction only
                assertTrue(at + ": Benjamini-Yekutieli undercut Benjamini-Hochberg", by[i] >= bh[i]);
            }
        }
    }

    @Test
    public void testBenjaminiYekutieliIsBenjaminiHochbergScaledByTheHarmonicNumber() {
        // that is the whole difference between the two, so it is worth pinning.
        // It is not bit for bit: the two sweeps round p * factor * m / j
        // differently, and the plan for this slice claimed an exact equality
        // that the measurement refused
        Lcg lcg = new Lcg(3L);
        double worst = 0.0;
        for (int round = 0; round < 800; round++) {
            int m = 1 + (int) (lcg.next() * 40);
            double[] p = awkwardFamily(m, lcg);
            double[] bh = MultipleTesting.benjaminiHochberg(p);
            double[] by = MultipleTesting.benjaminiYekutieli(p);
            double h = harmonic(m);
            for (int i = 0; i < m; i++) {
                double wanted = Math.min(1.0, bh[i] * h);
                if (wanted > 0.0) {
                    worst = Math.max(worst, Math.abs(by[i] - wanted) / wanted);
                } else {
                    assertEquals("a zero stays zero", 0.0, by[i], 0.0);
                }
            }
        }
        // measured over 86465 adjusted values: 2.4e-16, one unit in the last
        // place
        assertTrue("the two differ by " + worst + " relative", worst < 1.0e-14);
    }

    @Test
    public void testTheOrderOfTheInputIsTheOrderOfTheOutput() {
        Lcg lcg = new Lcg(5L);
        for (int round = 0; round < 400; round++) {
            int m = 2 + (int) (lcg.next() * 20);
            double[] p = awkwardFamily(m, lcg);
            double[] adjusted = MultipleTesting.benjaminiHochberg(p);

            // a Fisher-Yates shuffle of the family, carrying the answer with it
            int[] where = new int[m];
            for (int i = 0; i < m; i++) {
                where[i] = i;
            }
            for (int i = m - 1; i > 0; i--) {
                int j = (int) (lcg.next() * (i + 1));
                int swap = where[i];
                where[i] = where[j];
                where[j] = swap;
            }
            double[] permuted = new double[m];
            for (int i = 0; i < m; i++) {
                permuted[i] = p[where[i]];
            }
            double[] adjustedPermuted = MultipleTesting.benjaminiHochberg(permuted);
            for (int i = 0; i < m; i++) {
                // the same arithmetic on the same multiset, so bit for bit
                assertEquals("round=" + round + " i=" + i, adjusted[where[i]], adjustedPermuted[i], 0.0);
            }
        }
    }

    @Test
    public void testASingleTestIsNotCorrectedAtAll() {
        for (double p : new double[] { 0.0, 1.0e-12, 0.03, 0.5, 1.0 }) {
            double[] one = { p };
            assertEquals("Benjamini-Hochberg of one", p, MultipleTesting.benjaminiHochberg(one)[0], 0.0);
            assertEquals("Benjamini-Yekutieli of one", p, MultipleTesting.benjaminiYekutieli(one)[0], 0.0);
        }
        // and nothing can be made of a family that is all ones
        double[] ones = { 1.0, 1.0, 1.0, 1.0 };
        for (int i = 0; i < ones.length; i++) {
            assertEquals("a family of ones", 1.0, MultipleTesting.benjaminiHochberg(ones)[i], 0.0);
            assertEquals("a family of ones", 1.0, MultipleTesting.benjaminiYekutieli(ones)[i], 0.0);
        }
    }

    @Test
    public void testASmallFamilyWorkedOutByHand() {
        // six p-values. Walking the ranks downwards, with m = 6:
        //   k=6  0.600 * 6/6 = 0.6000  -> 0.6000
        //   k=5  0.042 * 6/5 = 0.0504  -> 0.0504
        //   k=4  0.041 * 6/4 = 0.0615  -> 0.0504   the running minimum holds
        //   k=3  0.039 * 6/3 = 0.0780  -> 0.0504   and holds again
        //   k=2  0.008 * 6/2 = 0.0240  -> 0.0240
        //   k=1  0.001 * 6/1 = 0.0060  -> 0.0060
        // Rank four is the case the whole step-up shape exists for: on its own
        // it would earn 0.0615 and be kept at the five percent level, and the
        // rank above it pulls it down to 0.0504 and rejects it
        double[] p = { 0.001, 0.008, 0.039, 0.041, 0.042, 0.60 };
        double[] wanted = { 0.006, 0.024, 0.0504, 0.0504, 0.0504, 0.60 };
        double[] got = MultipleTesting.benjaminiHochberg(p);
        for (int i = 0; i < p.length; i++) {
            assertEquals("i=" + i, wanted[i], got[i], 1.0e-15);
        }

        // five of the six were under 0.05 to begin with; two survive the
        // correction
        int uncorrected = 0;
        int corrected = 0;
        for (int i = 0; i < p.length; i++) {
            if (p[i] <= 0.05) {
                uncorrected++;
            }
            if (got[i] <= 0.05) {
                corrected++;
            }
        }
        assertEquals("uncorrected", 5, uncorrected);
        assertEquals("corrected", 2, corrected);

        // and Benjamini-Yekutieli scales all of it by H_6 = 2.45
        double[] yekutieli = MultipleTesting.benjaminiYekutieli(p);
        assertEquals("H_6", 2.45, harmonic(6), 1.0e-15);
        for (int i = 0; i < p.length; i++) {
            assertEquals("i=" + i, Math.min(1.0, wanted[i] * 2.45), yekutieli[i], 1.0e-15);
        }
    }

    // -------------------------------------------------------------- Holm --

    @Test
    public void testHolmRejectsExactlyWhatTheStepDownProcedureDoes() {
        // Holm defined a decision, not an adjusted p-value: judge the smallest
        // p-value against alpha / m, the next against alpha / (m - 1), and stop
        // at the first one that survives. That loop is written out here and
        // shares nothing with the implementation but its input
        Lcg lcg = new Lcg(4242L);
        int checked = 0;
        for (int round = 0; round < 1500; round++) {
            int m = 1 + (int) (lcg.next() * 30);
            double[] p = awkwardFamily(m, lcg);
            double[] adjusted = MultipleTesting.holmBonferroni(p);
            for (double alpha : new double[] { 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.9 }) {
                int byDecision = stepDownRejectionCount(p, alpha);
                int byAdjusted = 0;
                for (int i = 0; i < m; i++) {
                    if (adjusted[i] <= alpha) {
                        byAdjusted++;
                    }
                }
                assertEquals("m=" + m + " round=" + round + " alpha=" + alpha, byDecision, byAdjusted);
                checked++;
            }
        }
        assertTrue("nothing was checked", checked > 10000);
    }

    /** The original decision: how many of the ordered p-values are rejected. */
    private static int stepDownRejectionCount(double[] p, double alpha) {
        double[] ascending = p.clone();
        Arrays.sort(ascending);
        int m = ascending.length;
        for (int k = 1; k <= m; k++) {
            if (ascending[k - 1] > alpha / (m - k + 1)) {
                return k - 1;
            }
        }
        return m;
    }

    @Test
    public void testHolmControlsTheFamilyWiseErrorRate() {
        // the promise the other two do not make, measured rather than assumed:
        // with every null hypothesis true, the probability of rejecting any of
        // them at all stays at or below alpha. Uncorrected it is 1 - (1 - a)^m,
        // which is 40 percent for ten tests at five
        double alpha = 0.05;
        int reps = 4000;
        for (int m : new int[] { 3, 10, 40 }) {
            int holmWrong = 0;
            int uncorrectedWrong = 0;
            Lcg lcg = new Lcg(20260826L + m);
            for (int r = 0; r < reps; r++) {
                // every p-value is a draw from the null, which is uniform
                double[] p = new double[m];
                for (int i = 0; i < m; i++) {
                    p[i] = lcg.next();
                }
                double[] adjusted = MultipleTesting.holmBonferroni(p);
                boolean anyHolm = false;
                boolean anyRaw = false;
                for (int i = 0; i < m; i++) {
                    anyHolm |= adjusted[i] <= alpha;
                    anyRaw |= p[i] <= alpha;
                }
                if (anyHolm) {
                    holmWrong++;
                }
                if (anyRaw) {
                    uncorrectedWrong++;
                }
            }
            double holmRate = holmWrong / (double) reps;
            double rawRate = uncorrectedWrong / (double) reps;
            // measured at m = 3, 10 and 40: 0.048, 0.046 and 0.052 against an
            // uncorrected 0.14, 0.40 and 0.87. Three standard errors of a rate
            // of 0.05 over 4000 draws is 0.010
            assertTrue("m=" + m + ": Holm was wrong " + holmRate + " of the time", holmRate < alpha + 0.012);
            assertTrue("m=" + m + ": uncorrected was wrong " + rawRate + " of the time",
                    rawRate > 1.0 - Math.pow(1.0 - alpha, m) - 0.03);
        }
    }

    @Test
    public void testHolmIsNeverWorseThanBonferroniAndNeverBetterThanBenjaminiHochberg() {
        // the two procedures Holm sits between. It dominates Bonferroni
        // uniformly, which is why there is no Bonferroni method here, and it is
        // dominated by the false discovery rate procedures, which is the price
        // of the stronger guarantee
        Lcg lcg = new Lcg(13L);
        for (int round = 0; round < 800; round++) {
            int m = 1 + (int) (lcg.next() * 40);
            double[] p = awkwardFamily(m, lcg);
            double[] holm = MultipleTesting.holmBonferroni(p);
            double[] bh = MultipleTesting.benjaminiHochberg(p);
            for (int i = 0; i < m; i++) {
                String at = "round=" + round + " i=" + i + " p=" + p[i];
                assertTrue(at + ": " + holm[i] + " is below its own p-value", holm[i] >= p[i]);
                assertTrue(at + ": " + holm[i] + " is above Bonferroni", holm[i] <= Math.min(1.0, m * p[i]));
                assertTrue(at + ": " + holm[i] + " left [0, 1]", holm[i] >= 0.0 && holm[i] <= 1.0);
                // measured over 4000 families: Holm fell below Benjamini-
                // Hochberg by at most 1.1e-16, one unit in the last place at a
                // rank where the two agree exactly
                assertTrue(at + ": Holm undercut Benjamini-Hochberg", holm[i] >= bh[i] - 1.0e-15);
            }
        }
    }

    @Test
    public void testHolmIsMonotoneAndTiesAgree() {
        Lcg lcg = new Lcg(17L);
        for (int round = 0; round < 800; round++) {
            int m = 2 + (int) (lcg.next() * 40);
            double[] p = awkwardFamily(m, lcg);
            double[] holm = MultipleTesting.holmBonferroni(p);
            for (int i = 0; i < m; i++) {
                for (int k = 0; k < m; k++) {
                    String at = "round=" + round + " i=" + i + " k=" + k;
                    if (p[i] < p[k]) {
                        assertTrue(at + ": " + holm[i] + " > " + holm[k], holm[i] <= holm[k]);
                    }
                    if (p[i] == p[k]) {
                        // the running maximum gives this without a special
                        // case: the later of two equal p-values carries the
                        // smaller multiplier, so the maximum does not move
                        assertEquals(at + ": tied p-values", holm[i], holm[k], 0.0);
                    }
                }
            }
        }
    }

    @Test
    public void testHolmOnASmallFamilyWorkedOutByHand() {
        // the same six p-values the step-up example uses, with m = 6:
        //   k=1  0.001 * 6 = 0.006   -> 0.006
        //   k=2  0.008 * 5 = 0.040   -> 0.040
        //   k=3  0.039 * 4 = 0.156   -> 0.156
        //   k=4  0.041 * 3 = 0.123   -> 0.156   the running maximum holds
        //   k=5  0.042 * 2 = 0.084   -> 0.156   and holds again
        //   k=6  0.600 * 1 = 0.600   -> 0.600
        // Ranks four and five are what the step-down shape exists for: once a
        // p-value survives, nothing after it can be rejected however small its
        // own product is
        double[] p = { 0.001, 0.008, 0.039, 0.041, 0.042, 0.60 };
        double[] wanted = { 0.006, 0.040, 0.156, 0.156, 0.156, 0.60 };
        double[] got = MultipleTesting.holmBonferroni(p);
        for (int i = 0; i < p.length; i++) {
            assertEquals("i=" + i, wanted[i], got[i], 1.0e-15);
        }

        // five of the six were under 0.05 to begin with; two survive Holm, the
        // same two Benjamini-Hochberg keeps -- on a family this small the two
        // guarantees cost the same, and the adjusted values still differ
        int corrected = 0;
        for (int i = 0; i < p.length; i++) {
            if (got[i] <= 0.05) {
                corrected++;
            }
        }
        assertEquals("corrected", 2, corrected);
        double[] bh = MultipleTesting.benjaminiHochberg(p);
        assertTrue("the third rank agreed with Benjamini-Hochberg", got[2] > bh[2]);
    }

    @Test
    public void testHolmOfASingleTestIsTheTestItself() {
        for (double p : new double[] { 0.0, 1.0e-12, 0.03, 0.5, 1.0 }) {
            assertEquals("Holm of one", p, MultipleTesting.holmBonferroni(new double[] { p })[0], 0.0);
        }
    }

    @Test
    public void testTheOrderOfTheInputIsTheOrderOfHolmsOutput() {
        Lcg lcg = new Lcg(23L);
        for (int round = 0; round < 400; round++) {
            int m = 2 + (int) (lcg.next() * 20);
            double[] p = awkwardFamily(m, lcg);
            double[] adjusted = MultipleTesting.holmBonferroni(p);

            int[] where = new int[m];
            for (int i = 0; i < m; i++) {
                where[i] = i;
            }
            for (int i = m - 1; i > 0; i--) {
                int at = (int) (lcg.next() * (i + 1));
                int swap = where[i];
                where[i] = where[at];
                where[at] = swap;
            }
            double[] permuted = new double[m];
            for (int i = 0; i < m; i++) {
                permuted[i] = p[where[i]];
            }
            double[] adjustedPermuted = MultipleTesting.holmBonferroni(permuted);
            for (int i = 0; i < m; i++) {
                // the same arithmetic on the same multiset, so bit for bit
                assertEquals("round=" + round + " i=" + i, adjusted[where[i]], adjustedPermuted[i], 0.0);
            }
        }
    }

    // ------------------------------------------------------ the guard rail --

    @Test
    public void testMultipleTestingRejectsWhatItCannotCorrect() {
        rejects("a null family", null);
        rejects("an empty family", new double[0]);
        rejects("a NaN p-value", new double[] { 0.1, Double.NaN });
        rejects("a negative p-value", new double[] { 0.1, -0.01 });
        rejects("a p-value above one", new double[] { 0.1, 1.0 + 1.0e-9 });
        rejects("an infinite p-value", new double[] { 0.1, Double.POSITIVE_INFINITY });
    }

    private static void rejects(String what, double[] p) {
        try {
            MultipleTesting.benjaminiHochberg(p);
            fail("benjaminiHochberg accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            MultipleTesting.benjaminiYekutieli(p);
            fail("benjaminiYekutieli accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            MultipleTesting.holmBonferroni(p);
            fail("holmBonferroni accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    @Test
    public void testTheInputIsNotModified() {
        double[] p = { 0.6, 0.041, 0.001, 0.042, 0.008, 0.039 };
        double[] before = p.clone();
        MultipleTesting.benjaminiHochberg(p);
        MultipleTesting.benjaminiYekutieli(p);
        MultipleTesting.holmBonferroni(p);
        for (int i = 0; i < p.length; i++) {
            assertEquals("i=" + i, before[i], p[i], 0.0);
        }
    }
}
