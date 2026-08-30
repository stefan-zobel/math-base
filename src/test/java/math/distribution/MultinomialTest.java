package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.linalg.DMatrix;
import math.rng.DefaultRng;
import math.rng.MultinomialSampler;
import math.rng.PseudoRandom;

/**
 * The multinomial against classes that know nothing about it: {@link Binomial},
 * which it reduces to at two categories and which every marginal is, and
 * {@code MultinomialSampler}, which draws from it.
 * <p>
 * Most of what is asserted here is structural -- the mass over all outcomes,
 * the identity a merged pair of categories has to satisfy, a covariance whose
 * rows sum to zero -- rather than a golden value.
 */
public class MultinomialTest {

    /**
     * Every count vector of {@code n} draws over three categories, with the
     * last taken as the remainder.
     */
    private static double totalMassOverThree(Multinomial law, int n) {
        double sum = 0.0;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                sum += law.pmf(new int[] { i, j, n - i - j });
            }
        }
        return sum;
    }

    @Test
    public void theMassSumsToOneOverEveryCountVector() {
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 0.5, 2.0, 7.0 }, { 100.0, 1.0, 3.0 }, { 1.0e-8, 1.0, 1.0 },
                { 0.0, 3.0, 1.0 } };
        int[] totals = { 5, 12, 40 };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            for (int t = 0; t < totals.length; t++) {
                int n = totals[t];
                double sum = totalMassOverThree(new Multinomial(n, settings[s]), n);
                if (Math.abs(sum - 1.0) > worst) {
                    worst = Math.abs(sum - 1.0);
                    at = Arrays.toString(settings[s]) + " n=" + n + " sums to " + sum;
                }
            }
        }
        // and once with four categories, where there are 816 of them at n = 15
        double[][] four = { { 1.5, 0.3, 4.0, 2.0 }, { 1.0, 1.0, 1.0, 1.0 } };
        for (int s = 0; s < four.length; s++) {
            int[] ns = { 6, 15 };
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                Multinomial law = new Multinomial(n, four[s]);
                double sum = 0.0;
                for (int i = 0; i <= n; i++) {
                    for (int j = 0; j <= n - i; j++) {
                        for (int k = 0; k <= n - i - j; k++) {
                            sum += law.pmf(new int[] { i, j, k, n - i - j - k });
                        }
                    }
                }
                if (Math.abs(sum - 1.0) > worst) {
                    worst = Math.abs(sum - 1.0);
                    at = Arrays.toString(four[s]) + " n=" + n + " sums to " + sum;
                }
            }
        }
        // and over weights drawn at random, so the settings above are not the
        // only ones the sum has been seen to close over
        long lcg = 20260827L;
        for (int trial = 0; trial < 200; trial++) {
            double[] w = new double[3];
            for (int i = 0; i < w.length; i++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                w[i] = (lcg >>> 11) * 0x1.0p-53 + 1.0e-3;
            }
            int n = 4 + (trial % 20);
            double sum = totalMassOverThree(new Multinomial(n, w), n);
            if (Math.abs(sum - 1.0) > worst) {
                worst = Math.abs(sum - 1.0);
                at = Arrays.toString(w) + " n=" + n + " sums to " + sum;
            }
        }
        // measured over the settings above: 2.6e-14, at the uniform law with
        // n = 40, where 861 masses are added up
        assertTrue("the mass sums to " + worst + " off one, " + at, worst < 1.0e-12);
    }

    @Test
    public void atTwoCategoriesItIsABinomial() {
        double[] ps = { 0.5, 0.1, 0.9, 0.001, 0.999, 1.0e-7 };
        int[] ns = { 1, 2, 10, 200, 5000 };
        double worst = 0.0;
        String at = "";
        for (int a = 0; a < ps.length; a++) {
            double p = ps[a];
            for (int b = 0; b < ns.length; b++) {
                int n = ns[b];
                Multinomial law = new Multinomial(n, new double[] { p, 1.0 - p });
                Binomial ref = new Binomial(n, p);
                for (int k = 0; k <= n; k++) {
                    double mine = law.logPmf(new int[] { k, n - k });
                    double theirs = ref.logPmf(k);
                    if (Double.isInfinite(mine) && Double.isInfinite(theirs)) {
                        continue;
                    }
                    double relative = Math.abs(mine - theirs) / Math.max(1.0, Math.abs(theirs));
                    if (relative > worst) {
                        worst = relative;
                        at = "p=" + p + " n=" + n + " k=" + k + " : " + mine + " vs " + theirs;
                    }
                }
            }
        }
        // not exact, and cannot be: Binomial.logPmf reads the second category
        // as log1p(-p) where this class reads log(p_1 / p_0), and the two
        // differ in the last bit. At n = 5000 that bit is multiplied by the
        // count, so the worst agreement measured is 2.2e-12 relative -- at
        // p = 0.1, k = 497, where the second term carries a factor of 4503
        assertTrue("the two categories disagree with the binomial by " + worst + ", " + at, worst < 1.0e-11);
    }

    @Test
    public void mergingCategoriesGivesAMultinomial() {
        // the counts of two categories, added, are the counts of a single
        // category of their combined weight. Nothing in the class knows this
        double[][] settings = { { 0.5, 2.0, 7.0 }, { 1.0, 1.0, 1.0 }, { 1.0e-6, 5.0, 2.0 } };
        int[] ns = { 4, 11, 30 };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            double[] w = settings[s];
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                Multinomial three = new Multinomial(n, w);
                Multinomial two = new Multinomial(n, new double[] { w[0] + w[1], w[2] });
                for (int a = 0; a <= n; a++) {
                    double split = 0.0;
                    for (int b = 0; b <= a; b++) {
                        split += three.pmf(new int[] { b, a - b, n - a });
                    }
                    double merged = two.pmf(new int[] { a, n - a });
                    if (Math.abs(split - merged) > worst) {
                        worst = Math.abs(split - merged);
                        at = Arrays.toString(w) + " n=" + n + " a=" + a + " : " + split + " vs " + merged;
                    }
                }
            }
        }
        // measured: 7.8e-16
        assertTrue("the merged category is off by " + worst + ", " + at, worst < 1.0e-14);
    }

    @Test
    public void theWeightsNeedNotSumToOne() {
        double[] w = { 0.5, 2.0, 7.0 };
        double sum = 9.5;
        double[] normalized = new double[w.length];
        for (int i = 0; i < w.length; i++) {
            normalized[i] = w[i] / sum;
        }
        int n = 9;
        Multinomial raw = new Multinomial(n, w);
        Multinomial unit = new Multinomial(n, normalized);
        for (int i = 0; i < w.length; i++) {
            assertEquals("probability " + i + " is not normalized", normalized[i], raw.probability(i), 1.0e-15);
        }
        double worst = 0.0;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                int[] c = { i, j, n - i - j };
                worst = Math.max(worst, Math.abs(raw.pmf(c) - unit.pmf(c)));
            }
        }
        assertTrue("weights and probabilities disagree by " + worst, worst < 1.0e-15);
    }

    @Test
    public void aCategoryOfWeightZeroIsImpossible() {
        // the term such a category contributes is skipped rather than added as
        // a zero, so the law on the survivors comes back bit for bit -- and if
        // the skipping were dropped the product would read 0.0 * -Infinity and
        // every one of these would be NaN
        double[][] survivors = { { 1.0, 3.0 }, { 0.5, 2.0, 7.0 }, { 1.0, 1.0, 1.0, 1.0 } };
        int[] ns = { 0, 1, 5, 17 };
        for (int s = 0; s < survivors.length; s++) {
            double[] w = survivors[s];
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                double[] padded = new double[w.length + 2];
                System.arraycopy(w, 0, padded, 1, w.length);
                Multinomial full = new Multinomial(n, padded);
                Multinomial rest = new Multinomial(n, w);
                for (int a = 0; a <= n; a++) {
                    int[] c = new int[w.length];
                    c[0] = a;
                    c[c.length - 1] = n - a;
                    int[] p = new int[padded.length];
                    System.arraycopy(c, 0, p, 1, c.length);
                    assertEquals("the padded law differs at " + Arrays.toString(p), rest.logPmf(c), full.logPmf(p),
                            0.0);
                }
            }
        }
        // and counting one is not unlikely but impossible
        Multinomial zero = new Multinomial(5, new double[] { 0.0, 1.0, 3.0 });
        assertEquals("a count in a category of weight zero has mass", 0.0, zero.pmf(new int[] { 1, 1, 3 }), 0.0);
        assertTrue("a count in a category of weight zero has finite log mass",
                zero.logPmf(new int[] { 1, 1, 3 }) == Double.NEGATIVE_INFINITY);
    }

    @Test
    public void theLogMassAnswersWhereTheMassHasUnderflowed() {
        int n = 4000;
        Multinomial law = new Multinomial(n, new double[] { 0.25, 0.25, 0.25, 0.25 });
        int[] allInOne = { n, 0, 0, 0 };
        assertEquals("the mass has not underflowed", 0.0, law.pmf(allInOne), 0.0);
        // with every draw in one category the coefficient is one and the log
        // mass is n log(1/4) exactly, which is what the logarithm still says
        double expected = n * Math.log(0.25);
        assertEquals("the log mass is not n log(1/4)", expected, law.logPmf(allInOne),
                1.0e-11 * Math.abs(expected));
        assertTrue("the log mass is not finite", !Double.isInfinite(law.logPmf(allInOne)));
        // the mode of the same law does not underflow, so the two routes can be
        // compared where both are defined
        int[] mode = { 1000, 1000, 1000, 1000 };
        assertEquals("the mass and its logarithm disagree at the mode", Math.log(law.pmf(mode)), law.logPmf(mode),
                1.0e-12);
    }

    @Test
    public void theLogMassSurvivesAWeightWhoseShareDoesNot() {
        // the share of the first category is 1.0e-620 and is below the smallest
        // double, so log(p_i / p0) is -Infinity where the answer is finite and
        // ordinary. Ten decades of weight is what importance weights look like
        // just before a resampling step, and this is the law they describe
        double[] weights = { 1.0e-320, 1.0, 1.0e300 };
        Multinomial law = new Multinomial(3, weights);
        int[] counts = { 1, 1, 1 };
        double logMass = law.logPmf(counts);
        assertTrue("the log mass underflowed with the quotient : " + logMass, !Double.isInfinite(logMass));

        // log(3!) + log(1.0e-320/p0) + log(1.0/p0) + log(1.0e300/p0), with
        // p0 == 1.0e300 to the last bit
        double expected = Math.log(6.0) + (Math.log(1.0e-320) - Math.log(1.0e300)) + (-Math.log(1.0e300));
        assertEquals("the log mass is not the sum of the shares", expected, logMass, 1.0e-9);

        // and the one draw case is exactly what Categorical says it is
        Categorical single = new Categorical(weights);
        for (int k = 0; k < weights.length; k++) {
            int[] indicator = new int[weights.length];
            indicator[k] = 1;
            assertEquals("the one draw law differs at " + k, single.logPmf(k),
                    new Multinomial(1, weights).logPmf(indicator), 0.0);
        }
    }

    @Test
    public void aVeryLargeNumberOfTrialsStaysFinite() {
        // logFactorial has no upper bound and switches to a Stirling series
        // above 29, so a count vector no double could hold the mass of still
        // has a log mass
        int[] ns = { 100000, 1000000, 10000000 };
        for (int t = 0; t < ns.length; t++) {
            int n = ns[t];
            Multinomial law = new Multinomial(n, new double[] { 0.2, 0.3, 0.5 });
            int[] mode = { n / 5, (3 * n) / 10, n - n / 5 - (3 * n) / 10 };
            double logMass = law.logPmf(mode);
            assertTrue("the log mass at n=" + n + " is " + logMass, logMass < 0.0 && logMass > -50.0);
        }
    }

    @Test
    public void permutingTheProbabilitiesPermutesEverything() {
        double[] w = { 0.5, 2.0, 7.0 };
        double[] reversed = { 7.0, 2.0, 0.5 };
        int n = 8;
        Multinomial law = new Multinomial(n, w);
        Multinomial mirror = new Multinomial(n, reversed);
        double worst = 0.0;
        String at = "";
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                int[] c = { i, j, n - i - j };
                int[] flipped = { c[2], c[1], c[0] };
                double mine = law.logPmf(c);
                double theirs = mirror.logPmf(flipped);
                double relative = Math.abs(mine - theirs) / Math.max(1.0, Math.abs(theirs));
                if (relative > worst) {
                    worst = relative;
                    at = Arrays.toString(c) + " : " + mine + " vs " + theirs;
                }
            }
        }
        // equivariant, but not bit for bit, and it cannot be: logPmf adds the
        // terms in category order and reversing the categories reverses the
        // summation, which floating point addition does not survive unchanged.
        // Measured over five weight vectors and n up to 300: 4.9e-14 relative
        assertTrue("the mirrored law is off by " + worst + ", " + at, worst < 1.0e-12);
        for (int i = 0; i < w.length; i++) {
            assertEquals("probability " + i + " is not mirrored", law.probability(i),
                    mirror.probability(w.length - 1 - i), 1.0e-15);
        }
    }

    @Test
    public void theInputIsCopiedRatherThanKept() {
        double[] w = { 1.0, 2.0, 3.0 };
        Multinomial law = new Multinomial(6, w);
        double before = law.logPmf(new int[] { 1, 2, 3 });
        Arrays.fill(w, 1.0);
        assertEquals("writing to the array the constructor was given changed the law", before,
                law.logPmf(new int[] { 1, 2, 3 }), 0.0);
    }

    @Test
    public void theCountsMustSumToTheTrials() {
        // a vector that does not sum to n is not an outcome of small
        // probability but no outcome at all, and saying so is more use than
        // returning a zero the caller has to interpret
        Multinomial law = new Multinomial(7, new double[] { 1.0, 2.0, 4.0 });
        refusesCounts("counts summing to less than n", law, new int[] { 1, 2, 3 });
        refusesCounts("counts summing to more than n", law, new int[] { 3, 3, 3 });
        refusesCounts("an empty count vector", law, new int[] { 0, 0, 0 });
        // the exact sum is enough, wherever it is spread
        law.logPmf(new int[] { 7, 0, 0 });
        law.logPmf(new int[] { 0, 0, 7 });
        law.logPmf(new int[] { 2, 2, 3 });
    }

    @Test
    public void theMultinomialRejectsWhatItCannotBe() {
        refuses("a negative number of trials", -1, new double[] { 1.0, 1.0 });
        refuses("null weights", 3, null);
        refuses("empty weights", 3, new double[0]);
        refuses("a negative weight", 3, new double[] { 1.0, -1.0e-30 });
        refuses("a NaN weight", 3, new double[] { 1.0, Double.NaN });
        refuses("an infinite weight", 3, new double[] { 1.0, Double.POSITIVE_INFINITY });
        refuses("weights that are all zero", 3, new double[] { 0.0, 0.0 });

        Multinomial law = new Multinomial(4, new double[] { 1.0, 1.0 });
        refusesCounts("null counts", law, null);
        refusesCounts("too few counts", law, new int[] { 4 });
        refusesCounts("too many counts", law, new int[] { 2, 1, 1 });
        refusesCounts("a negative count", law, new int[] { 5, -1 });

        try {
            law.probability(2);
            fail("probability accepted a category that does not exist");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            law.probability(-1);
            fail("probability accepted a negative category");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    @Test
    public void theMarginalAgreesWithTheMomentsAndTheRowsSumToZero() {
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 0.5, 2.0, 7.0 }, { 100.0, 1.0, 3.0 }, { 1.0e-8, 1.0, 1.0 },
                { 0.0, 3.0, 1.0 }, { 1.5, 0.3, 4.0, 2.0 } };
        int[] ns = { 1, 10, 1000, 1000000 };
        double worstVariance = 0.0;
        double worstRow = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            double[] w = settings[s];
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                Multinomial law = new Multinomial(n, w);
                DMatrix c = law.covariance();
                double[] mean = new double[w.length];
                law.mean(mean);
                for (int i = 0; i < w.length; i++) {
                    Binomial marginal = law.marginal(i);
                    // both are n times the same quotient, so this one is exact
                    assertEquals("the marginal mean is not the mean", marginal.mean(), mean[i], 0.0);
                    double variance = Math.abs(marginal.variance() - c.getUnsafe(i, i))
                            / Math.max(1.0, Math.abs(marginal.variance()));
                    if (variance > worstVariance) {
                        worstVariance = variance;
                        at = Arrays.toString(w) + " n=" + n + " i=" + i;
                    }
                    double row = 0.0;
                    for (int j = 0; j < w.length; j++) {
                        row += c.getUnsafe(i, j);
                    }
                    // the row sum is a difference of numbers the size of the
                    // diagonal, so it is that scale it has to be small against
                    worstRow = Math.max(worstRow, Math.abs(row) / Math.max(1.0e-300, c.getUnsafe(i, i)));
                }
            }
        }
        // measured: 5.9e-16. Not exact, and again for a stated reason --
        // Binomial.variance() is n p (1 - p) where the diagonal here is
        // n p_i (sum of the others), which is the same number by a different
        // route. Where the two genuinely part company is the test below
        assertTrue("the marginal variance is off the diagonal by " + worstVariance + ", " + at,
                worstVariance < 1.0e-14);
        // measured: 2.8e-16 of the diagonal, at four categories and n = 1e6
        assertTrue("a covariance row sums to " + worstRow + " of its diagonal", worstRow < 1.0e-14);
    }

    @Test
    public void theMarginalIsTheAggregatedMass() {
        // the marginal is a Binomial by an argument about single draws; here it
        // is checked against the multinomial mass itself, added over every
        // outcome that puts k draws in the category
        double[][] settings = { { 0.5, 2.0, 7.0 }, { 1.0, 1.0, 1.0 }, { 1.0e-6, 5.0, 2.0 }, { 0.0, 3.0, 1.0 } };
        int[] ns = { 3, 9, 25 };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                Multinomial law = new Multinomial(n, settings[s]);
                for (int i = 0; i < 3; i++) {
                    Binomial marginal = law.marginal(i);
                    for (int k = 0; k <= n; k++) {
                        double aggregated = 0.0;
                        for (int a = 0; a <= n; a++) {
                            for (int b = 0; b <= n - a; b++) {
                                int[] counts = { a, b, n - a - b };
                                if (counts[i] == k) {
                                    aggregated += law.pmf(counts);
                                }
                            }
                        }
                        if (Math.abs(aggregated - marginal.pmf(k)) > worst) {
                            worst = Math.abs(aggregated - marginal.pmf(k));
                            at = Arrays.toString(settings[s]) + " n=" + n + " i=" + i + " k=" + k;
                        }
                    }
                }
            }
        }
        // measured: 5.4e-15
        assertTrue("the marginal is off the aggregated mass by " + worst + ", " + at, worst < 1.0e-13);
    }

    @Test
    public void theCovarianceSurvivesAWeightTheDifferenceWouldNot() {
        // one category within an ulp of the whole weight. Its probability
        // rounds to exactly one, so 1 - p_i is exactly zero and the textbook
        // n p (1 - p) reports a variance of zero where there is one -- as does
        // Binomial.variance(), which computes it that way. Summing the other
        // weights instead keeps the number
        int n = 1000;
        double tiny = 1.0e-18;
        Multinomial law = new Multinomial(n, new double[] { 1.0, tiny, tiny });
        assertEquals("the leading probability does not round to one", 1.0, law.probability(0), 0.0);
        assertEquals("there is a complement left to cancel", 0.0, 1.0 - law.probability(0), 0.0);

        DMatrix c = law.covariance();
        double expected = n * (2.0 * tiny);
        assertEquals("the covariance diagonal has cancelled away", expected, c.getUnsafe(0, 0),
                1.0e-12 * expected);
        assertEquals("the naive product would not have cancelled", 0.0,
                n * law.probability(0) * (1.0 - law.probability(0)), 0.0);
        assertEquals("the binomial marginal would not have cancelled", 0.0, law.marginal(0).variance(), 0.0);
        // and the row still sums to zero, since the off-diagonals cancelled the
        // same way
        assertEquals("the row does not sum to zero", 0.0,
                c.getUnsafe(0, 0) + c.getUnsafe(0, 1) + c.getUnsafe(0, 2), 1.0e-12 * expected);
    }

    @Test
    public void noTrialsPutNothingAnywhere() {
        Multinomial law = new Multinomial(0, new double[] { 1.0, 2.0, 3.0 });
        assertEquals("the empty outcome is not certain", 1.0, law.pmf(new int[] { 0, 0, 0 }), 0.0);
        assertEquals("the empty outcome has a log mass", 0.0, law.logPmf(new int[] { 0, 0, 0 }), 0.0);
        double[] mean = new double[3];
        law.mean(mean);
        DMatrix c = law.covariance();
        for (int i = 0; i < 3; i++) {
            assertEquals("no trials but a mean at " + i, 0.0, mean[i], 0.0);
            assertEquals("no trials but a marginal mean at " + i, 0.0, law.marginal(i).mean(), 0.0);
            for (int j = 0; j < 3; j++) {
                assertEquals("no trials but a covariance at " + i + "," + j, 0.0, c.getUnsafe(i, j), 0.0);
            }
        }
    }

    @Test
    public void aSingleCategoryTakesEverything() {
        for (int n = 0; n <= 5; n++) {
            Multinomial law = new Multinomial(n, new double[] { 2.5 });
            assertEquals("the single category is not certain", 1.0, law.pmf(new int[] { n }), 0.0);
            assertEquals("the single category has a log mass", 0.0, law.logPmf(new int[] { n }), 0.0);
            assertEquals("the single probability is not one", 1.0, law.probability(0), 0.0);
            // where Dirichlet.marginal has to refuse -- a Beta has no shape of
            // zero -- a degenerate Binomial is an ordinary object
            Binomial marginal = law.marginal(0);
            assertEquals("the degenerate marginal is not certain", 1.0, marginal.pmf(n), 1.0e-15);
            assertEquals("the degenerate marginal has a mean other than n", (double) n, marginal.mean(), 0.0);
            assertEquals("the degenerate marginal has a variance", 0.0, marginal.variance(), 0.0);
            assertEquals("the one by one covariance is not zero", 0.0, law.covariance().getUnsafe(0, 0), 0.0);
        }
    }

    @Test
    public void theMomentsRejectWhatTheyCannotWrite() {
        Multinomial law = new Multinomial(4, new double[] { 1.0, 1.0, 2.0 });
        refusesOut("null", law, null);
        refusesOut("too short an array", law, new double[2]);
        refusesOut("too long an array", law, new double[4]);
        try {
            law.marginal(3);
            fail("marginal accepted a category that does not exist");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            law.marginal(-1);
            fail("marginal accepted a negative category");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    @Test
    public void theDrawIsTheSamplersDraw() {
        // this class holds a sampler so that a caller need not build one, and
        // that is all it does with it
        double[] w = { 0.7, 1.3, 4.0, 0.0, 2.5 };
        int n = 60;
        Multinomial law = new Multinomial(n, w);
        MultinomialSampler sampler = MultinomialSampler.of(w);
        for (long seed = 1L; seed <= 20L; seed++) {
            int[] viaLaw = new int[w.length];
            law.sample(DefaultRng.newPseudoRandom(seed), viaLaw);
            int[] viaSampler = new int[w.length];
            sampler.sample(DefaultRng.newPseudoRandom(seed), n, viaSampler);
            assertTrue("the two routes differ at seed " + seed + " : " + Arrays.toString(viaLaw) + " vs "
                    + Arrays.toString(viaSampler), Arrays.equals(viaLaw, viaSampler));
        }
    }

    @Test
    public void everyDrawIsAnOutcomeTheMassAllows() {
        // the sampler places conditional binomials and knows nothing of the
        // mass function; logPmf refuses anything that does not sum to n and
        // answers minus infinity for a category of weight zero. So a draw that
        // logPmf accepts and gives finite log mass to is the two halves of the
        // class agreeing on what an outcome is
        double[][] settings = { { 0.7, 1.3, 4.0, 2.5 }, { 0.0, 3.0, 1.0, 0.0 }, { 1.0e-9, 1.0, 1.0, 1.0 } };
        int[] ns = { 0, 1, 25, 4000 };
        long lcg = 20260827L;
        for (int s = 0; s < settings.length; s++) {
            double[] w = settings[s];
            for (int t = 0; t < ns.length; t++) {
                int n = ns[t];
                Multinomial law = new Multinomial(n, w);
                for (int trial = 0; trial < 50; trial++) {
                    lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                    int[] counts = new int[w.length];
                    law.sample(DefaultRng.newPseudoRandom(lcg), counts);
                    long total = 0L;
                    for (int i = 0; i < counts.length; i++) {
                        assertTrue("a negative count " + Arrays.toString(counts), counts[i] >= 0);
                        total += counts[i];
                    }
                    assertEquals("the counts do not sum to n : " + Arrays.toString(counts), n, total);
                    // would throw if the counts missed n, and would be minus
                    // infinity if a category of weight zero had been counted
                    double logMass = law.logPmf(counts);
                    assertTrue("the sampler produced an impossible outcome " + Arrays.toString(counts) + " of "
                            + Arrays.toString(w), logMass > Double.NEGATIVE_INFINITY);
                }
            }
        }
    }

    @Test
    public void theInterfaceDefaultDrawsWhatTheSamplerDoes() {
        // PseudoRandom.nextMultinomial builds a sampler per call, so it must
        // land on the same counts as one built once
        double[] w = { 2.0, 5.0, 1.0 };
        int n = 40;
        MultinomialSampler sampler = MultinomialSampler.of(w);
        for (long seed = 1L; seed <= 20L; seed++) {
            PseudoRandom p = DefaultRng.newPseudoRandom(seed);
            int[] viaDefault = new int[w.length];
            p.nextMultinomial(n, w, viaDefault);
            int[] viaSampler = new int[w.length];
            sampler.sample(DefaultRng.newPseudoRandom(seed), n, viaSampler);
            assertTrue("the interface default differs at seed " + seed + " : " + Arrays.toString(viaDefault)
                    + " vs " + Arrays.toString(viaSampler), Arrays.equals(viaDefault, viaSampler));
        }
    }

    @Test
    public void theDrawRejectsWhatItCannotWrite() {
        Multinomial law = new Multinomial(5, new double[] { 1.0, 1.0, 1.0 });
        PseudoRandom p = DefaultRng.newPseudoRandom(1L);
        refusesDraw("a null generator", law, null, new int[3]);
        refusesDraw("null counts", law, p, null);
        refusesDraw("too short an array", law, p, new int[2]);
        refusesDraw("too long an array", law, p, new int[4]);
    }

    private static void refusesDraw(String what, Multinomial law, PseudoRandom prng, int[] counts) {
        try {
            law.sample(prng, counts);
            fail("sample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refusesOut(String what, Multinomial law, double[] out) {
        try {
            law.mean(out);
            fail("mean accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refuses(String what, int n, double[] probabilities) {
        try {
            new Multinomial(n, probabilities);
            fail("Multinomial accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refusesCounts(String what, Multinomial law, int[] counts) {
        try {
            law.logPmf(counts);
            fail("logPmf accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            law.pmf(counts);
            fail("pmf accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
