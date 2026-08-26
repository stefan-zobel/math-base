package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.distribution.Binomial;
import math.distribution.DiscreteDistribution;
import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;

/**
 * The multinomial sampler, checked against its marginals -- which are binomial
 * and have a mass function in {@code math.distribution} -- and then against the
 * two things marginals cannot see: that the counts sum to the total exactly,
 * and that the categories are negatively correlated by the right amount.
 */
public final class MultinomialSamplerTest {

    /**
     * Cells of a discrete law each carrying an expected count of at least five,
     * which is what the chi-squared approximation asks for. Six standard
     * deviations is not far enough for a pooled tail: at forty thousand draws
     * its expected count is near {@code 1e-7} there, and one draw landing in it
     * sends the statistic to a million.
     */
    static final class Bins {
        private final int[] firstOfCell;
        private final double[] share;

        private Bins(int[] firstOfCell, double[] share) {
            this.firstOfCell = firstOfCell;
            this.share = share;
        }

        static Bins of(DiscreteDistribution law, int draws) {
            double sd = Math.sqrt(law.variance());
            long span = (long) Math.ceil(12.0 * sd) + 2L;
            int lo = (int) Math.max(law.supportLowerBound(), law.mean() - span);
            int hi = (int) Math.min(law.supportUpperBound(), (long) (law.mean() + span));
            int[] edges = new int[hi - lo + 2];
            double[] shares = new double[hi - lo + 2];
            int cells = 0;
            double accumulated = law.cdf(lo - 1);
            edges[0] = Integer.MIN_VALUE;
            for (int k = lo; k <= hi; k++) {
                accumulated += law.pmf(k);
                if (accumulated * draws >= 5.0) {
                    shares[cells++] = accumulated;
                    edges[cells] = k + 1;
                    accumulated = 0.0;
                }
            }
            shares[cells++] = accumulated + 1.0 - law.cdf(hi);
            if (cells > 1 && shares[cells - 1] * draws < 5.0) {
                shares[cells - 2] += shares[cells - 1];
                cells--;
            }
            double[] share = new double[cells];
            int[] first = new int[cells];
            double total = 0.0;
            for (int i = 0; i < cells; i++) {
                share[i] = shares[i];
                first[i] = edges[i];
                total += share[i];
            }
            for (int i = 0; i < cells; i++) {
                share[i] /= total;
            }
            return new Bins(first, share);
        }

        int cells() {
            return share.length;
        }

        private int cellOf(int k) {
            int lo = 0;
            int hi = firstOfCell.length - 1;
            while (lo < hi) {
                int mid = (lo + hi + 1) >>> 1;
                if (firstOfCell[mid] <= k) {
                    lo = mid;
                } else {
                    hi = mid - 1;
                }
            }
            return lo;
        }

        double fit(int[] drawn) {
            long[] observed = new long[share.length];
            for (int i = 0; i < drawn.length; i++) {
                observed[cellOf(drawn[i])]++;
            }
            double[] expected = new double[share.length];
            for (int i = 0; i < share.length; i++) {
                expected[i] = drawn.length * share[i];
            }
            return HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 0).pValue;
        }
    }

    // ------------------------------------------------------- the marginals --

    @Test
    public void testEveryMarginalIsTheBinomialItShouldBe() {
        Uniform uniform = new Uniform(0.0, 1.0);
        double[][] cases = { { 0.25, 0.25, 0.25, 0.25 }, { 0.7, 0.2, 0.09, 0.01 }, { 0.5, 0.5 },
                { 0.4, 0.0, 0.6 } };
        int n = 40;
        int reps = 200;
        int draws = 20000;
        for (double[] probabilities : cases) {
            MultinomialSampler sampler = MultinomialSampler.of(probabilities);
            for (int category = 0; category < probabilities.length; category++) {
                if (probabilities[category] == 0.0) {
                    continue;
                }
                Bins bins = Bins.of(new Binomial(n, probabilities[category]), draws);
                assertTrue("the binning collapsed", bins.cells() >= 3);
                double[] pValues = new double[reps];
                for (int r = 0; r < reps; r++) {
                    PseudoRandom prng = DefaultRng.newPseudoRandom(r * 7919L + 1);
                    int[] counts = new int[probabilities.length];
                    int[] marginal = new int[draws];
                    for (int d = 0; d < draws; d++) {
                        sampler.sample(prng, n, counts);
                        marginal[d] = counts[category];
                    }
                    pValues[r] = bins.fit(marginal);
                }
                double uniformity = HypothesisTests
                        .kolmogorovSmirnov(pValues, uniform, Alternative.TWO_SIDED).pValue;
                // measured over the four vectors: 0.0706 .. 0.5382
                assertTrue(Arrays.toString(probabilities) + " category " + category
                        + ": uniformity p = " + uniformity, uniformity > 0.001);
            }
        }
    }

    // ------------------------------------------- what marginals cannot see --

    @Test
    public void testTheCountsSumToTheTotalOnEveryDraw() {
        // this is the invariant the conditional construction exists to give,
        // and it is exact rather than nearly: the last category takes what is
        // left rather than being drawn
        PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
        double[][] cases = { { 0.5, 0.3, 0.2 }, { 0.25, 0.25, 0.25, 0.25 }, { 0.999, 0.001 },
                { 1.0, 0.0, 0.0, 0.0 } };
        for (double[] probabilities : cases) {
            MultinomialSampler sampler = MultinomialSampler.of(probabilities);
            int[] counts = new int[probabilities.length];
            for (int n : new int[] { 0, 1, 7, 30, 5000 }) {
                for (int d = 0; d < 20000; d++) {
                    sampler.sample(prng, n, counts);
                    int total = 0;
                    for (int i = 0; i < counts.length; i++) {
                        assertTrue("a negative count", counts[i] >= 0);
                        total += counts[i];
                    }
                    assertEquals(Arrays.toString(probabilities) + " n=" + n, n, total);
                }
            }
        }
    }

    @Test
    public void testTheCategoriesAreCorrelatedByTheRightAmount() {
        // a sampler that got every marginal right and the dependence wrong
        // would pass everything above. Cov(k_i, k_j) = -n p_i p_j is the
        // statement about the joint law
        int n = 30;
        int draws = 500000;
        double[][] cases = { { 0.5, 0.3, 0.2 }, { 0.25, 0.25, 0.25, 0.25 } };
        for (double[] probabilities : cases) {
            int m = probabilities.length;
            MultinomialSampler sampler = MultinomialSampler.of(probabilities);
            PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
            double[] sum = new double[m];
            double[][] cross = new double[m][m];
            int[] counts = new int[m];
            for (int d = 0; d < draws; d++) {
                sampler.sample(prng, n, counts);
                for (int i = 0; i < m; i++) {
                    sum[i] += counts[i];
                    for (int j = i; j < m; j++) {
                        cross[i][j] += (double) counts[i] * counts[j];
                    }
                }
            }
            for (int i = 0; i < m; i++) {
                // the diagonal is the marginal variance, n p (1 - p)
                double variance = cross[i][i] / draws - (sum[i] / draws) * (sum[i] / draws);
                assertEquals("variance of category " + i, n * probabilities[i] * (1.0 - probabilities[i]),
                        variance, 0.05);
                for (int j = i + 1; j < m; j++) {
                    double covariance = cross[i][j] / draws - (sum[i] / draws) * (sum[j] / draws);
                    // measured over two million draws: within 0.008 of the
                    // truth, so half a million leaves room at 0.05
                    assertEquals("covariance of " + i + " and " + j,
                            -n * probabilities[i] * probabilities[j], covariance, 0.05);
                }
            }
        }
    }

    @Test
    public void testTwoCategoriesAreABinomial() {
        // the chain of conditionals degenerates to a single binomial draw here,
        // so the two must agree in distribution
        int n = 25;
        double p = 0.3;
        MultinomialSampler sampler = MultinomialSampler.of(new double[] { p, 1.0 - p });
        Bins bins = Bins.of(new Binomial(n, p), 200000);
        PseudoRandom prng = DefaultRng.newPseudoRandom(11L);
        int[] counts = new int[2];
        int[] first = new int[200000];
        for (int d = 0; d < first.length; d++) {
            sampler.sample(prng, n, counts);
            first[d] = counts[0];
        }
        double pValue = bins.fit(first);
        // measured: the mean came out 7.4975 against 7.5
        assertTrue("the first category does not fit Binomial(25, 0.3), p = " + pValue, pValue > 0.001);
    }

    // ------------------------------------------------------------ the edges --

    @Test
    public void testACategoryOfWeightZeroIsNeverCounted() {
        // the suffix-sum form of the conditional probabilities is what makes
        // this exact: where every later weight is zero the quotient is exactly
        // one, so the leftover cannot reach a category that should not have it
        MultinomialSampler sampler = MultinomialSampler.of(new double[] { 0.0, 0.5, 0.0, 0.5, 0.0 });
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        int[] counts = new int[5];
        long[] seen = new long[5];
        for (int d = 0; d < 200000; d++) {
            sampler.sample(prng, 10, counts);
            for (int i = 0; i < 5; i++) {
                seen[i] += counts[i];
            }
        }
        assertEquals("category 0", 0L, seen[0]);
        assertEquals("category 2", 0L, seen[2]);
        assertEquals("category 4", 0L, seen[4]);
        assertEquals("the two that carry the weight", 2000000L, seen[1] + seen[3]);

        // and the same where the trailing zeros are the last categories, which
        // is where a prefix form would leak one
        MultinomialSampler trailing = MultinomialSampler.of(new double[] { 1.0, 2.0, 0.0, 0.0 });
        assertEquals("the conditional at the last weighted category", 1.0,
                trailing.conditionalProbability(1), 0.0);
        int[] four = new int[4];
        for (int d = 0; d < 200000; d++) {
            trailing.sample(prng, 12, four);
            assertEquals("a trailing zero was counted", 0, four[2]);
            assertEquals("a trailing zero was counted", 0, four[3]);
        }
    }

    @Test
    public void testTheSmallAndTheLargeSampler() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(7L);
        int[] one = new int[1];
        MultinomialSampler single = MultinomialSampler.of(new double[] { 4.0 });
        assertEquals("categories", 1, single.categories());
        single.sample(prng, 17, one);
        assertEquals("one category takes everything", 17, one[0]);

        double[] flat = new double[1000];
        Arrays.fill(flat, 1.0);
        int[] many = new int[1000];
        MultinomialSampler.of(flat).sample(prng, 100000, many);
        int total = 0;
        int occupied = 0;
        for (int i = 0; i < many.length; i++) {
            total += many[i];
            if (many[i] > 0) {
                occupied++;
            }
        }
        assertEquals("the total", 100000, total);
        // a hundred thousand over a thousand equal categories leaves none empty
        assertEquals("every category was reached", 1000, occupied);
    }

    @Test
    public void testTheWeightsNeedNotSumToOne() {
        int[] asProbabilities = new int[3];
        int[] asCounts = new int[3];
        MultinomialSampler.of(new double[] { 0.5, 0.3, 0.2 })
                .sample(DefaultRng.newPseudoRandom(31L), 50, asProbabilities);
        MultinomialSampler.of(new double[] { 500.0, 300.0, 200.0 })
                .sample(DefaultRng.newPseudoRandom(31L), 50, asCounts);
        for (int i = 0; i < 3; i++) {
            assertEquals("at " + i, asProbabilities[i], asCounts[i]);
        }
    }

    @Test
    public void testTheSamplerIsReproducibleFromTheSeed() {
        MultinomialSampler sampler = MultinomialSampler.of(new double[] { 0.4, 0.35, 0.25 });
        int[] first = new int[3];
        int[] again = new int[3];
        PseudoRandom one = DefaultRng.newPseudoRandom(23L);
        PseudoRandom two = DefaultRng.newPseudoRandom(23L);
        for (int d = 0; d < 500; d++) {
            sampler.sample(one, 40, first);
            sampler.sample(two, 40, again);
            for (int i = 0; i < 3; i++) {
                assertEquals("draw " + d + " at " + i, first[i], again[i]);
            }
        }
    }

    // --------------------------------------------------------- the guard rail --

    @Test
    public void testTheSamplerRejectsWhatItCannotDraw() {
        rejectsBuild("null weights", null);
        rejectsBuild("no weights", new double[0]);
        rejectsBuild("a negative weight", new double[] { 1.0, -0.5 });
        rejectsBuild("a NaN weight", new double[] { 1.0, Double.NaN });
        rejectsBuild("an infinite weight", new double[] { 1.0, Double.POSITIVE_INFINITY });
        rejectsBuild("weights that are all zero", new double[] { 0.0, 0.0 });

        MultinomialSampler sampler = MultinomialSampler.of(new double[] { 1.0, 1.0, 1.0 });
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        rejectsDraw("a null generator", sampler, null, 10, new int[3]);
        rejectsDraw("a null output array", sampler, prng, 10, null);
        rejectsDraw("a negative total", sampler, prng, -1, new int[3]);
        rejectsDraw("an output array of the wrong length", sampler, prng, 10, new int[2]);
    }

    private static void rejectsBuild(String what, double[] probabilities) {
        try {
            MultinomialSampler.of(probabilities);
            fail("of accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsDraw(String what, MultinomialSampler sampler, PseudoRandom prng, int n,
            int[] counts) {
        try {
            sampler.sample(prng, n, counts);
            fail("sample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
