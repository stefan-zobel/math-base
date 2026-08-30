package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.DiscreteDistribution;
import math.distribution.Geometric;
import math.distribution.Hypergeometric;
import math.distribution.NegativeBinomial;
import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;

/**
 * The three samplers that needed no algorithm of their own: geometric by
 * inversion, negative binomial as a gamma-Poisson mixture, hypergeometric by a
 * walk over its own recurrence.
 */
public final class WaitingTimeSamplerTest {

    private interface Draw {
        int next(PseudoRandom prng);
    }

    /**
     * The p-value of the claim that many chi-squared fits against the law's own
     * mass function have uniform p-values, which is the statement that a sampler
     * is right rather than that one run of it was not obviously wrong.
     */
    private static double uniformityOfFits(DiscreteDistribution law, Draw draw, int reps, int draws) {
        Bins bins = Bins.of(law, draws);
        assertTrue("the binning collapsed to " + bins.cells() + " cells", bins.cells() >= 3);
        double[] pValues = new double[reps];
        for (int i = 0; i < reps; i++) {
            PseudoRandom prng = DefaultRng.newPseudoRandom(i * 7919L + 1);
            int[] drawn = new int[draws];
            for (int j = 0; j < draws; j++) {
                drawn[j] = draw.next(prng);
            }
            pValues[i] = bins.fit(drawn);
        }
        return HypothesisTests.kolmogorovSmirnov(pValues, new Uniform(0.0, 1.0), Alternative.TWO_SIDED)
                .pValue;
    }

    /**
     * Cells of a discrete law, each carrying an expected count of at least
     * five. Six standard deviations is not far enough for a pooled tail: at
     * that distance its expected count is near {@code 1e-7}, and one draw
     * landing there sends the chi-squared statistic to a million.
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
            long hiWanted = Math.min(law.supportUpperBound(), (long) (law.mean() + span));
            int hi = (int) Math.min(hiWanted, lo + 400000L);

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

    // ---------------------------------------------------------- geometric --

    @Test
    public void testTheGeometricFitsItsMassFunction() {
        for (final double p : new double[] { 0.9, 0.5, 0.1, 0.01, 0.001 }) {
            Geometric law = new Geometric(p);
            double uniformity = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    return GeometricSpliterator.sampleFor(prng, p);
                }
            }, 120, 40000);
            // measured: 0.0765 .. 0.7553 over the five probabilities
            assertTrue("p=" + p + ": uniformity p = " + uniformity, uniformity > 0.001);
        }
    }

    @Test
    public void testTheGeometricEdges() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        // a certain success never fails first
        int[] certain = prng.geometric(500, 1.0).toArray();
        for (int i = 0; i < certain.length; i++) {
            assertEquals("p = 1", 0, certain[i]);
        }
        // and the smallest probability offered still fits an int
        int[] rare = prng.geometric(500, GeometricSpliterator.MIN_PROBABILITY).toArray();
        for (int i = 0; i < rare.length; i++) {
            assertTrue("a count of " + rare[i] + " is negative", rare[i] >= 0);
        }
        // the mean of a geometric is (1-p)/p, which for p = 0.25 is 3
        double mean = mean(prng.geometric(400000, 0.25).toArray());
        assertEquals("the mean at p = 0.25", 3.0, mean, 0.05);
    }

    // -------------------------------------------------- negative binomial --

    @Test
    public void testTheNegativeBinomialFitsItsMassFunction() {
        int[][] shapes = { { 1, 5 }, { 3, 5 }, { 10, 5 }, { 2, 1 }, { 50, 8 } };
        for (int[] shape : shapes) {
            final int r = shape[0];
            final double p = shape[1] / 10.0;
            NegativeBinomial law = new NegativeBinomial(r, p);
            double uniformity = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    return NegativeBinomialSpliterator.sampleFor(prng, r, p);
                }
            }, 120, 40000);
            // measured: 0.5536 .. 0.8969 over the five shapes
            assertTrue("r=" + r + " p=" + p + ": uniformity p = " + uniformity, uniformity > 0.001);
        }
    }

    @Test
    public void testTheMixtureAgreesWithTheDefinitionItStandsFor() {
        // the sampler is a gamma-Poisson mixture, which is an identity rather
        // than the definition. The definition is counting failures against
        // Bernoulli draws, and it is written out here: an independent
        // algorithm, not a rearrangement of the same one
        for (int[] shape : new int[][] { { 3, 4 }, { 8, 6 } }) {
            final int r = shape[0];
            final double p = shape[1] / 10.0;
            NegativeBinomial law = new NegativeBinomial(r, p);
            double mixture = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    return NegativeBinomialSpliterator.sampleFor(prng, r, p);
                }
            }, 120, 40000);
            double bernoulli = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    int successes = 0;
                    int failures = 0;
                    while (successes < r) {
                        if (prng.nextDouble() < p) {
                            successes++;
                        } else {
                            failures++;
                        }
                    }
                    return failures;
                }
            }, 120, 40000);
            String at = "r=" + r + " p=" + p;
            // measured: mixture 0.1700 and 0.2781, Bernoulli 0.4183 and 0.2173
            assertTrue(at + ": the mixture, uniformity p = " + mixture, mixture > 0.001);
            assertTrue(at + ": the definition, uniformity p = " + bernoulli, bernoulli > 0.001);
        }
    }

    @Test
    public void testTheNegativeBinomialEdges() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        int[] certain = prng.negativeBinomial(200, 4, 1.0).toArray();
        for (int i = 0; i < certain.length; i++) {
            assertEquals("p = 1", 0, certain[i]);
        }
        // one success to wait for is the geometric distribution
        NegativeBinomial law = new NegativeBinomial(1, 0.3);
        assertEquals("the mean", law.mean(), mean(prng.negativeBinomial(400000, 1, 0.3).toArray()), 0.05);
    }

    // ----------------------------------------------------- hypergeometric --

    @Test
    public void testTheHypergeometricFitsItsMassFunction() {
        int[][] shapes = { { 50, 20, 10 }, { 1000, 100, 50 }, { 100, 50, 50 }, { 500, 5, 400 },
                { 60, 30, 3 } };
        for (int[] shape : shapes) {
            final int population = shape[0];
            final int successes = shape[1];
            final int draws = shape[2];
            Hypergeometric law = new Hypergeometric(population, successes, draws);
            double uniformity = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    return HypergeometricSpliterator.sample(prng, population, successes, draws);
                }
            }, 120, 40000);
            // measured: 0.4612 .. 0.8111 over the five shapes
            assertTrue("N=" + population + " K=" + successes + " n=" + draws + ": uniformity p = "
                    + uniformity, uniformity > 0.001);
        }
    }

    @Test
    public void testTheHypergeometricStaysInsideItsSupport() {
        // drawing everything, or drawing from a population that is all
        // successes, leaves no freedom at all
        PseudoRandom prng = DefaultRng.newPseudoRandom(9L);
        int[] all = prng.hypergeometric(200, 40, 15, 40).toArray();
        for (int i = 0; i < all.length; i++) {
            assertEquals("drawing the whole population", 15, all[i]);
        }
        int[] none = prng.hypergeometric(200, 40, 0, 10).toArray();
        for (int i = 0; i < none.length; i++) {
            assertEquals("no successes to find", 0, none[i]);
        }
        int[] forced = prng.hypergeometric(200, 40, 35, 30).toArray();
        for (int i = 0; i < forced.length; i++) {
            // 30 drawn from 40 of which 35 succeed: at least 25 must succeed
            assertTrue("a count of " + forced[i] + " is outside [25, 30]", forced[i] >= 25 && forced[i] <= 30);
        }
        int[] empty = prng.hypergeometric(50, 0, 0, 0).toArray();
        for (int i = 0; i < empty.length; i++) {
            assertEquals("an empty draw", 0, empty[i]);
        }
    }

    // ------------------------------------------------------------ the wiring --

    @Test
    public void testTheStreamsAreReproducibleAndSplit() {
        assertReproducible(DefaultRng.newPseudoRandom(7L).geometric(200, 0.2).toArray(),
                DefaultRng.newPseudoRandom(7L).geometric(200, 0.2).toArray(), "geometric");
        assertReproducible(DefaultRng.newPseudoRandom(7L).negativeBinomial(200, 4, 0.3).toArray(),
                DefaultRng.newPseudoRandom(7L).negativeBinomial(200, 4, 0.3).toArray(), "negative binomial");
        assertReproducible(DefaultRng.newPseudoRandom(7L).hypergeometric(200, 100, 30, 20).toArray(),
                DefaultRng.newPseudoRandom(7L).hypergeometric(200, 100, 30, 20).toArray(), "hypergeometric");

        assertTrue("the geometric spliterator refused to split",
                new GeometricSpliterator(DefaultRng.newPseudoRandom(1L), 0L, 1024L, 0.2).trySplit() != null);
        assertTrue("the negative binomial spliterator refused to split",
                new NegativeBinomialSpliterator(DefaultRng.newPseudoRandom(1L), 0L, 1024L, 3, 0.4)
                        .trySplit() != null);
        assertTrue("the hypergeometric spliterator refused to split",
                new HypergeometricSpliterator(DefaultRng.newPseudoRandom(1L), 0L, 1024L, 100, 30, 20)
                        .trySplit() != null);
    }

    private static void assertReproducible(int[] first, int[] again, String what) {
        assertEquals(what + ": length", first.length, again.length);
        for (int i = 0; i < first.length; i++) {
            assertEquals(what + " at " + i, first[i], again[i]);
        }
    }

    @Test
    public void testEveryGeneratorCanDrawAllThree() {
        PseudoRandom[] generators = { new Interleaved4Stc64(1L), new Lcg64Xor1024Mix(2L),
                new MersenneTwister64(3L), new Sfc64(4L), new SplitMix64(5L), new Stc64(6L),
                new XorShift1024Star(7L), new XorShift128Plus(8L), new XorShiftRot256StarStar(9L) };
        for (int g = 0; g < generators.length; g++) {
            String name = generators[g].getAlgorithm();
            assertEquals(name + ": the geometric mean", new Geometric(0.2).mean(),
                    mean(generators[g].geometric(200000, 0.2).toArray()), 0.1);
            assertEquals(name + ": the negative binomial mean", new NegativeBinomial(4, 0.3).mean(),
                    mean(generators[g].negativeBinomial(200000, 4, 0.3).toArray()), 0.15);
            assertEquals(name + ": the hypergeometric mean", new Hypergeometric(200, 100, 30).mean(),
                    mean(generators[g].hypergeometric(200000, 200, 100, 30).toArray()), 0.05);
        }
    }

    private static double mean(int[] drawn) {
        double sum = 0.0;
        for (int i = 0; i < drawn.length; i++) {
            sum += drawn[i];
        }
        return sum / drawn.length;
    }

    @Test
    public void testTheWaitingTimeSamplersRejectWhatTheyCannotDraw() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        rejects("a probability above one", prng, 0);
        rejects("a probability below the bound", prng, 1);
        rejects("a probability that is not a number", prng, 2);
        rejects("no successes to wait for", prng, 3);
        rejects("a probability of zero", prng, 4);
        rejects("a mean beyond an int", prng, 5);
        rejects("a negative population", prng, 6);
        rejects("more successes than population", prng, 7);
        rejects("more draws than population", prng, 8);
        rejects("a negative stream size", prng, 9);
    }

    private static void rejects(String what, PseudoRandom prng, int which) {
        try {
            switch (which) {
            case 0:
                prng.geometric(1.5);
                break;
            case 1:
                prng.geometric(1.0e-12);
                break;
            case 2:
                prng.geometric(Double.NaN);
                break;
            case 3:
                prng.negativeBinomial(0, 0.5);
                break;
            case 4:
                prng.negativeBinomial(3, 0.0);
                break;
            case 5:
                prng.negativeBinomial(1000, 1.0e-9);
                break;
            case 6:
                prng.hypergeometric(-1, 0, 0);
                break;
            case 7:
                prng.hypergeometric(10, 11, 5);
                break;
            case 8:
                prng.hypergeometric(10, 5, 11);
                break;
            default:
                prng.hypergeometric(-1L, 10, 5, 5);
                break;
            }
            fail("accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
