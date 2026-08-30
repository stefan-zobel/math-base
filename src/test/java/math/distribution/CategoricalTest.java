package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.rng.AliasTable;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;

/**
 * {@link Categorical}, the one law in this package whose masses are given
 * rather than computed.
 * <p>
 * That makes most of the usual invariants exact statements about bookkeeping
 * rather than about a formula, and they are asserted as such: against
 * {@link DiscreteUniform} and against {@link Multinomial} at one draw the
 * agreement is bit for bit. Where a tolerance appears it was measured over
 * several thousand weight vectors spanning ten decades and is roughly two
 * orders of magnitude above the worst figure seen.
 */
public final class CategoricalTest {

    private static final double[][] SETTINGS = { { 1.0, 1.0, 1.0 }, { 0.5, 2.0, 7.0 }, { 100.0, 1.0, 3.0 },
            { 3.0, 0.0, 5.0, 0.0, 1.0 }, { 1.0e-6, 1.0, 1.0e6 }, { 7.0 },
            { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 } };

    // ------------------------------------------------------------------------
    // mass, log mass, distribution function
    // ------------------------------------------------------------------------

    @Test
    public void theMassSumsToOne() {
        // over the random weights below the worst residual was 1.0e-15
        for (double[] w : SETTINGS) {
            Categorical c = new Categorical(w);
            double sum = 0.0;
            for (int k = 0; k < c.categories(); k++) {
                sum += c.pmf(k);
            }
            assertEquals("the mass of " + Arrays.toString(w) + " does not sum to one", 1.0, sum, 1.0e-13);
        }
        forEachRandomLaw((w, c) -> {
            double sum = 0.0;
            for (int k = 0; k < c.categories(); k++) {
                sum += c.pmf(k);
            }
            assertEquals("the mass of " + Arrays.toString(w) + " does not sum to one", 1.0, sum, 1.0e-13);
        });
    }

    @Test
    public void theCdfIsTheRunningSumOfTheMass() {
        forEachRandomLaw((w, c) -> {
            double running = 0.0;
            for (int k = 0; k < c.categories(); k++) {
                running += c.pmf(k);
                assertEquals("the cdf is not the running mass at " + k, running, c.cdf(k), 1.0e-13);
            }
            // and it closes exactly, which is what the last entry is assigned
            // rather than accumulated for
            assertEquals("the cdf does not close at the last category", 1.0, c.cdf(c.categories() - 1), 0.0);
        });
    }

    @Test
    public void theCdfIsMonotoneAndClampedOutside() {
        forEachRandomLaw((w, c) -> {
            int m = c.categories();
            assertEquals("mass below the first category", 0.0, c.cdf(-1), 0.0);
            assertEquals("mass below the first category", 0.0, c.cdf(Integer.MIN_VALUE), 0.0);
            assertEquals("the cdf is not one above the last category", 1.0, c.cdf(m), 0.0);
            assertEquals("the cdf is not one above the last category", 1.0, c.cdf(Integer.MAX_VALUE), 0.0);
            double previous = 0.0;
            for (int k = 0; k < m; k++) {
                assertTrue("the cdf fell at " + k + " for " + Arrays.toString(w), c.cdf(k) >= previous);
                previous = c.cdf(k);
            }
        });
    }

    @Test
    public void outsideTheCategoriesThereIsNoMass() {
        Categorical c = new Categorical(new double[] { 2.0, 3.0, 5.0 });
        for (int k : new int[] { -1, 3, 4, Integer.MIN_VALUE, Integer.MAX_VALUE }) {
            assertEquals("mass at " + k, 0.0, c.pmf(k), 0.0);
            assertEquals("log mass at " + k, Double.NEGATIVE_INFINITY, c.logPmf(k), 0.0);
        }
    }

    @Test
    public void theWeightsNeedNotSumToOne() {
        // counts and the probabilities they normalize to are the same law; the
        // worst gap measured was 8.9e-16
        forEachRandomLaw((w, c) -> {
            double sum = 0.0;
            for (int i = 0; i < w.length; i++) {
                sum += w[i];
            }
            double[] normalized = new double[w.length];
            for (int i = 0; i < w.length; i++) {
                normalized[i] = w[i] / sum;
            }
            Categorical n = new Categorical(normalized);
            for (int k = 0; k < w.length; k++) {
                assertEquals("the mass differs at " + k, n.pmf(k), c.pmf(k), 1.0e-14);
                assertEquals("the cdf differs at " + k, n.cdf(k), c.cdf(k), 1.0e-14);
            }
        });
    }

    @Test
    public void aCategoryOfWeightZeroIsImpossible() {
        double[] w = { 4.0, 0.0, 1.0, 0.0, 5.0 };
        Categorical c = new Categorical(w);
        assertEquals("a category of weight zero has mass", 0.0, c.pmf(1), 0.0);
        assertEquals("a category of weight zero has mass", 0.0, c.pmf(3), 0.0);
        assertEquals("a category of weight zero has a log mass", Double.NEGATIVE_INFINITY, c.logPmf(1), 0.0);
        assertEquals("a category of weight zero has a log mass", Double.NEGATIVE_INFINITY, c.logPmf(3), 0.0);

        // and it takes nothing from the others: the law on the surviving
        // categories is the same one, bit for bit
        Categorical without = new Categorical(new double[] { 4.0, 1.0, 5.0 });
        int[] survivors = { 0, 2, 4 };
        for (int i = 0; i < survivors.length; i++) {
            assertEquals("the mass moved at " + survivors[i], without.pmf(i), c.pmf(survivors[i]), 0.0);
            assertEquals("the log mass moved at " + survivors[i], without.logPmf(i), c.logPmf(survivors[i]), 0.0);
        }
    }

    @Test
    public void theLogMassAnswersWhereTheQuotientHasUnderflowed() {
        // the reason the constructor guards the quotient. Ten decades of weight
        // is what importance weights look like just before a resampling step,
        // and Math.log(1.0e-320 / 1.0e300) is -Infinity where the answer is
        // finite and perfectly ordinary
        Categorical c = new Categorical(new double[] { 1.0e-320, 1.0, 1.0e300 });

        // the smallest share is below the smallest double and is gone; the
        // middle one is 1.0e-300 and survives, which is the line the guard has
        // to find
        assertEquals("the share should have underflowed", 0.0, c.pmf(0), 0.0);
        assertEquals("the middle share should have survived", 1.0e-300, c.pmf(1), 1.0e-315);
        assertTrue("the log mass underflowed with the quotient : " + c.logPmf(0),
                c.logPmf(0) > Double.NEGATIVE_INFINITY);

        // log(1e-320) - log(1e300) and log(1) - log(1e300)
        assertEquals("the smallest weight", Math.log(1.0e-320) - Math.log(1.0e300), c.logPmf(0), 1.0e-9);
        assertEquals("the middle weight", -Math.log(1.0e300), c.logPmf(1), 1.0e-9);
        assertEquals("the dominant weight takes everything", 0.0, c.logPmf(2), 1.0e-12);
    }

    @Test
    public void permutingTheCategoriesPermutesTheMass() {
        // the mass is a lookup, so this is nearly exact -- but the sum of the
        // weights is accumulated in category order, so reversing them moves its
        // last digit. Worst gap measured: 6.7e-16 on the mass, 1.0e-15 relative
        // on the log mass
        forEachRandomLaw((w, c) -> {
            int m = w.length;
            double[] reversed = new double[m];
            for (int i = 0; i < m; i++) {
                reversed[i] = w[m - 1 - i];
            }
            Categorical r = new Categorical(reversed);
            for (int k = 0; k < m; k++) {
                assertEquals("the mass is not mirrored at " + k, r.pmf(m - 1 - k), c.pmf(k), 1.0e-14);
                double a = c.logPmf(k);
                double b = r.logPmf(m - 1 - k);
                if (Double.isInfinite(a) && Double.isInfinite(b)) {
                    continue;
                }
                assertEquals("the log mass is not mirrored at " + k, 0.0,
                        Math.abs(a - b) / Math.max(1.0, Math.abs(a)), 1.0e-13);
            }
        });
    }

    @Test
    public void theInputIsCopiedRatherThanKept() {
        double[] w = { 1.0, 2.0, 1.0 };
        Categorical c = new Categorical(w);
        double before = c.pmf(1);
        double cdfBefore = c.cdf(1);
        double meanBefore = c.mean();
        Arrays.fill(w, 99.0);
        assertEquals("the mass followed the caller's array", before, c.pmf(1), 0.0);
        assertEquals("the cdf followed the caller's array", cdfBefore, c.cdf(1), 0.0);
        assertEquals("the mean followed the caller's array", meanBefore, c.mean(), 0.0);
    }

    @Test
    public void theCategoricalRejectsWhatItCannotBe() {
        refuses("a null weight vector", null);
        refuses("an empty weight vector", new double[0]);
        refuses("a negative weight", new double[] { 1.0, -1.0 });
        refuses("a NaN weight", new double[] { 1.0, Double.NaN });
        refuses("an infinite weight", new double[] { 1.0, Double.POSITIVE_INFINITY });
        refuses("weights that are all zero", new double[] { 0.0, 0.0 });
        refuses("a single zero weight", new double[] { 0.0 });
        refuses("weights that overflow", new double[] { Double.MAX_VALUE, Double.MAX_VALUE });
    }

    // ------------------------------------------------------------------------
    // moments, support, quantile
    // ------------------------------------------------------------------------

    @Test
    public void atUniformWeightsItIsTheDiscreteUniform() {
        // two classes that share no line of code, so this is a real
        // cross-check -- and it holds bit for bit over every size up to 400,
        // mass, distribution function, both moments and every quantile
        for (int m = 1; m <= 400; m++) {
            double[] w = new double[m];
            Arrays.fill(w, 1.0);
            Categorical c = new Categorical(w);
            DiscreteUniform u = new DiscreteUniform(0, m - 1);
            for (int k = -1; k <= m; k++) {
                assertEquals("m=" + m + ": the mass differs at " + k, u.pmf(k), c.pmf(k), 0.0);
                assertEquals("m=" + m + ": the cdf differs at " + k, u.cdf(k), c.cdf(k), 0.0);
            }
            assertEquals("m=" + m + ": the mean differs", u.mean(), c.mean(), 0.0);
            assertEquals("m=" + m + ": the variance differs", u.variance(), c.variance(), 0.0);
            assertEquals("m=" + m + ": the support starts elsewhere", u.supportLowerBound(), c.supportLowerBound());
            assertEquals("m=" + m + ": the support ends elsewhere", u.supportUpperBound(), c.supportUpperBound());
            for (int q = 0; q <= 500; q++) {
                double p = q / 500.0;
                assertEquals("m=" + m + ": the quantile differs at " + p, u.quantile(p), c.quantile(p));
            }
        }
    }

    @Test
    public void atTwoCategoriesItIsABernoulli() {
        for (double q : new double[] { 0.01, 0.25, 0.5, 0.9, 0.999 }) {
            Categorical c = new Categorical(new double[] { 1.0 - q, q });
            Binomial b = new Binomial(1, q);
            for (int k = -1; k <= 2; k++) {
                assertEquals("q=" + q + ": the mass differs at " + k, b.pmf(k), c.pmf(k), 1.0e-15);
                assertEquals("q=" + q + ": the cdf differs at " + k, b.cdf(k), c.cdf(k), 1.0e-15);
            }
            assertEquals("q=" + q + ": the mean differs", b.mean(), c.mean(), 1.0e-15);
            assertEquals("q=" + q + ": the variance differs", b.variance(), c.variance(), 1.0e-15);
        }
    }

    @Test
    public void theOneDrawMultinomialIsThisLaw() {
        // Multinomial(1, w) is this distribution written as a vector of
        // indicators, and it computes its log mass through a multinomial
        // coefficient that is zero here. The two agree bit for bit
        for (double[] w : SETTINGS) {
            Categorical c = new Categorical(w);
            Multinomial mn = new Multinomial(1, w);
            for (int k = 0; k < w.length; k++) {
                int[] indicator = new int[w.length];
                indicator[k] = 1;
                assertEquals("the log mass differs at " + k + " of " + Arrays.toString(w), mn.logPmf(indicator),
                        c.logPmf(k), 0.0);
                assertEquals("the mass differs at " + k + " of " + Arrays.toString(w), mn.pmf(indicator), c.pmf(k),
                        1.0e-15);
            }
            // and the moments say the same thing one category at a time
            double[] mean = new double[w.length];
            mn.mean(mean);
            for (int k = 0; k < w.length; k++) {
                assertEquals("the marginal mean differs at " + k, mean[k], c.pmf(k), 1.0e-15);
            }
        }
    }

    @Test
    public void theQuantileInvertsTheCdf() {
        // cdf(q - 1) < p <= cdf(q), the bracket the interface promises. Checked
        // over 800000 pairs in a probe without a single violation
        forEachRandomLaw((w, c) -> {
            for (int i = 0; i <= 200; i++) {
                double p = i / 200.0;
                int k = c.quantile(p);
                assertTrue("quantile(" + p + ") = " + k + " undershoots the cdf for " + Arrays.toString(w),
                        c.cdf(k) >= p);
                if (p > 0.0 && k > 0) {
                    assertTrue("quantile(" + p + ") = " + k + " overshoots the cdf for " + Arrays.toString(w),
                            c.cdf(k - 1) < p);
                }
            }
            for (int k = 0; k < c.categories(); k++) {
                assertTrue("quantile(cdf(" + k + ")) ran past " + k, c.quantile(c.cdf(k)) <= k);
            }
        });
    }

    @Test
    public void theQuantileNeverLandsWhereThereIsNoMass() {
        double[][] holes = { { 0.0, 0.0, 3.0, 1.0 }, { 2.0, 0.0, 0.0 }, { 0.0, 5.0, 0.0 },
                { 1.0, 0.0, 2.0, 0.0, 3.0, 0.0 } };
        for (double[] w : holes) {
            Categorical c = new Categorical(w);
            for (int i = 0; i <= 1000; i++) {
                double p = i / 1000.0;
                int k = c.quantile(p);
                assertTrue("quantile(" + p + ") = " + k + " has no weight in " + Arrays.toString(w),
                        c.logPmf(k) > Double.NEGATIVE_INFINITY);
            }
        }
    }

    @Test
    public void theSupportIsWhereTheMassIs() {
        // a weight of zero is not a small probability but no support at all, so
        // the bounds skip the leading and trailing zeros -- which is what makes
        // quantile(0) and quantile(1) land on categories that can occur
        Categorical c = new Categorical(new double[] { 0.0, 0.0, 3.0, 0.0, 1.0, 0.0 });
        assertEquals("the support starts at a weight of zero", 2, c.supportLowerBound());
        assertEquals("the support ends at a weight of zero", 4, c.supportUpperBound());
        assertEquals("quantile(0) is outside the support", 2, c.quantile(0.0));
        assertEquals("quantile(1) is outside the support", 4, c.quantile(1.0));

        Categorical full = new Categorical(new double[] { 1.0, 2.0, 3.0 });
        assertEquals(0, full.supportLowerBound());
        assertEquals(2, full.supportUpperBound());
    }

    @Test
    public void theMomentsAreTheSumsTheyClaimToBe() {
        forEachRandomLaw((w, c) -> {
            double sum = 0.0;
            double first = 0.0;
            for (int i = 0; i < w.length; i++) {
                sum += w[i];
                first += i * w[i];
            }
            double mean = first / sum;
            double second = 0.0;
            for (int i = 0; i < w.length; i++) {
                double d = i - mean;
                second += w[i] * d * d;
            }
            assertEquals("the mean is not the weighted index", mean, c.mean(), 1.0e-12);
            assertEquals("the variance is not the weighted spread", second / sum, c.variance(), 1.0e-12);
            assertTrue("a negative variance : " + c.variance(), c.variance() >= 0.0);
        });
    }

    @Test
    public void aSingleCategoryTakesEverything() {
        Categorical c = new Categorical(new double[] { 7.0 });
        assertEquals("the only category does not take everything", 1.0, c.pmf(0), 0.0);
        assertEquals("the only category does not take everything", 0.0, c.logPmf(0), 0.0);
        assertEquals(1.0, c.cdf(0), 0.0);
        assertEquals("a degenerate law has a mean elsewhere", 0.0, c.mean(), 0.0);
        assertEquals("a degenerate law has a spread", 0.0, c.variance(), 0.0);
        assertEquals(0, c.supportLowerBound());
        assertEquals(0, c.supportUpperBound());
        assertEquals(0, c.quantile(0.0));
        assertEquals(0, c.quantile(0.5));
        assertEquals(0, c.quantile(1.0));
    }

    @Test
    public void theQuantileRejectsWhatItCannotAnswer() {
        Categorical c = new Categorical(new double[] { 1.0, 2.0, 3.0 });
        for (double p : new double[] { -1.0e-15, -1.0, 1.0 + 1.0e-15, 2.0, Double.NaN,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY }) {
            try {
                c.quantile(p);
                fail("the quantile accepted p = " + p);
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------------
    // the draw
    // ------------------------------------------------------------------------

    @Test
    public void theDrawIsTheTablesDraw() {
        // the distribution holds a table so that a caller need not build one;
        // from the same seed it must be the identical draw, which makes every
        // statement AliasTableTest establishes carry over whole
        double[] w = { 2.0, 5.0, 1.0, 0.0, 3.0 };
        for (long seed : new long[] { 1L, 42L, -7L, 987654321L }) {
            Categorical c = new Categorical(w);
            PseudoRandom one = DefaultRng.newPseudoRandom(seed);
            PseudoRandom two = DefaultRng.newPseudoRandom(seed);
            AliasTable table = AliasTable.of(w);
            for (int i = 0; i < 500; i++) {
                assertEquals("seed " + seed + ", draw " + i, table.sample(two), c.sample(one));
            }
        }
    }

    @Test
    public void theTableItHandsOutIsTheOneItDrawsFrom() {
        double[] w = { 4.0, 1.0, 0.0, 5.0 };
        Categorical c = new Categorical(w);
        assertEquals("the table covers other categories", c.categories(), c.aliasTable().outcomes());
        PseudoRandom one = DefaultRng.newPseudoRandom(3L);
        PseudoRandom two = DefaultRng.newPseudoRandom(3L);
        for (int i = 0; i < 200; i++) {
            assertEquals("draw " + i, c.aliasTable().sample(two), c.sample(one));
        }
    }

    @Test
    public void everyDrawIsACategoryTheMassAllows() {
        double[] w = { 3.0, 0.0, 1.0, 0.0, 6.0 };
        Categorical c = new Categorical(w);
        PseudoRandom prng = DefaultRng.newPseudoRandom(20260827L);
        int[] seen = new int[w.length];
        for (int i = 0; i < 20000; i++) {
            int k = c.sample(prng);
            assertTrue("a draw left the categories : " + k, k >= 0 && k < w.length);
            assertTrue("a category of weight zero was drawn : " + k, c.logPmf(k) > Double.NEGATIVE_INFINITY);
            seen[k]++;
        }
        for (int k = 0; k < w.length; k++) {
            if (w[k] > 0.0) {
                assertTrue("category " + k + " was never drawn", seen[k] > 0);
            }
        }
    }

    @Test
    public void theDrawRejectsWhatItCannotDrawFrom() {
        Categorical c = new Categorical(new double[] { 1.0, 1.0 });
        try {
            c.sample(null);
            fail("the draw accepted a null generator");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without a message", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------------

    /** What a body handed a random weight vector and the law over it does. */
    private interface Check {
        void run(double[] weights, Categorical law);
    }

    /**
     * Two hundred weight vectors of up to forty categories, spread over ten
     * decades, every twentieth of them holding a weight of exactly zero.
     * Deterministic: the generator is the LCG the tests in this package use.
     */
    private static void forEachRandomLaw(Check check) {
        long lcg = 20260827L;
        for (int trial = 0; trial < 200; trial++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            int m = 1 + (int) ((lcg >>> 33) % 40L);
            double[] w = new double[m];
            boolean positive = false;
            for (int i = 0; i < m; i++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double u = (lcg >>> 11) * 0x1.0p-53;
                w[i] = (trial % 20 == 0 && i % 7 == 3) ? 0.0 : Math.pow(10.0, -5.0 + 10.0 * u);
                positive |= w[i] > 0.0;
            }
            if (!positive) {
                continue;
            }
            check.run(w, new Categorical(w));
        }
    }

    private static void refuses(String what, double[] weights) {
        try {
            new Categorical(weights);
            fail("the constructor accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
