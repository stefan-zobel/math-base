package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import math.distribution.Binomial;
import math.distribution.DiscreteDistribution;
import math.distribution.Poisson;
import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;
import math.stats.TestResult;

/**
 * The Poisson and binomial samplers, checked against the mass functions in
 * {@code math.distribution} that were written independently of them.
 * <p>
 * Both have two branches that are exact over different ranges, so both are
 * exercised on both sides of the crossover: what a crossover has to be safe to
 * move, and it can only be that if either branch is right where the other one
 * runs.
 */
public final class CountSamplerTest {

    /** What a sampler produces, so that a branch can be pinned without a stream. */
    private interface Draw {
        int next(PseudoRandom prng);
    }

    /**
     * A binning of a discrete law in which every cell carries an expected count
     * of at least five, which is what the chi-squared approximation asks for.
     * <p>
     * Six standard deviations is not far enough out here: the pooled tail of a
     * Poisson at that distance has an expected count near {@code 1e-7}, and one
     * draw landing in it sends the statistic to a million. Cells are therefore
     * merged from the bottom until each is heavy enough, and the last short one
     * is folded into its predecessor.
     */
    private static final class Binning {
        private final int[] firstOfCell;
        private final double[] share;

        private Binning(int[] firstOfCell, double[] share) {
            this.firstOfCell = firstOfCell;
            this.share = share;
        }

        static Binning of(DiscreteDistribution law, int draws) {
            double sd = Math.sqrt(law.variance());
            int lo = (int) Math.max(law.supportLowerBound(), Math.floor(law.mean() - 12.0 * sd) - 1);
            int hi = (int) Math.min(law.supportUpperBound(), Math.ceil(law.mean() + 12.0 * sd) + 1);
            List<Integer> edges = new ArrayList<Integer>();
            List<Double> shares = new ArrayList<Double>();
            double accumulated = law.cdf(lo - 1);
            edges.add(Integer.MIN_VALUE);
            for (int k = lo; k <= hi; k++) {
                accumulated += law.pmf(k);
                if (accumulated * draws >= 5.0) {
                    shares.add(accumulated);
                    edges.add(k + 1);
                    accumulated = 0.0;
                }
            }
            shares.add(accumulated + 1.0 - law.cdf(hi));

            int cells = shares.size();
            double[] share = new double[cells];
            int[] first = new int[cells];
            for (int i = 0; i < cells; i++) {
                share[i] = shares.get(i).doubleValue();
                first[i] = edges.get(i).intValue();
            }
            if (cells > 1 && share[cells - 1] * draws < 5.0) {
                share[cells - 2] += share[cells - 1];
                double[] s = new double[cells - 1];
                int[] f = new int[cells - 1];
                System.arraycopy(share, 0, s, 0, cells - 1);
                System.arraycopy(first, 0, f, 0, cells - 1);
                share = s;
                first = f;
            }
            double total = 0.0;
            for (int i = 0; i < share.length; i++) {
                total += share[i];
            }
            for (int i = 0; i < share.length; i++) {
                share[i] /= total;
            }
            return new Binning(first, share);
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

    /**
     * The p-value of the claim that the p-values of many chi-squared fits are
     * uniform, which is the statement that a sampler is right rather than that
     * one run of it was not obviously wrong.
     */
    private static double uniformityOfFits(DiscreteDistribution law, Draw draw, int reps, int draws) {
        Binning binning = Binning.of(law, draws);
        assertTrue("the binning collapsed to " + binning.cells() + " cells", binning.cells() >= 3);
        double[] pValues = new double[reps];
        for (int i = 0; i < reps; i++) {
            PseudoRandom prng = DefaultRng.newPseudoRandom(i * 7919L + 1);
            int[] drawn = new int[draws];
            for (int j = 0; j < draws; j++) {
                drawn[j] = draw.next(prng);
            }
            pValues[i] = binning.fit(drawn);
        }
        TestResult result = HypothesisTests.kolmogorovSmirnov(pValues, new Uniform(0.0, 1.0),
                Alternative.TWO_SIDED);
        return result.pValue;
    }

    // ------------------------------------------------------------- Poisson --

    @Test
    public void testBothPoissonBranchesFitTheMassFunction() {
        // the multiplication is exact everywhere; the rejection's constants
        // were derived for a mean of ten or more, and it is only asked about
        // above that
        for (final double lambda : new double[] { 0.1, 1.0, 5.0, 29.999, 30.0, 200.0 }) {
            Poisson law = new Poisson(lambda);
            double byMultiplication = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    return PoissonSpliterator.byMultiplication(prng, lambda);
                }
            }, 120, 40000);
            // measured: 0.0419 .. 0.6490 over six means
            assertTrue("lambda=" + lambda + ": multiplication, uniformity p = " + byMultiplication,
                    byMultiplication > 0.001);

            if (lambda >= 10.0) {
                double byRejection = uniformityOfFits(law, new Draw() {
                    @Override
                    public int next(PseudoRandom prng) {
                        return PoissonSpliterator.byTransformedRejection(prng, lambda);
                    }
                }, 120, 40000);
                // measured: 0.2604 .. 0.8927
                assertTrue("lambda=" + lambda + ": rejection, uniformity p = " + byRejection,
                        byRejection > 0.001);
            }
        }
    }

    @Test
    public void testThePoissonBranchesAgreeAtTheCrossover() {
        // a crossover is only safe to move if both branches are right where the
        // other one runs, so they are compared just below and just above it
        for (double lambda : new double[] { PoissonSpliterator.REJECTION_THRESHOLD - 0.001,
                PoissonSpliterator.REJECTION_THRESHOLD + 0.001 }) {
            Poisson law = new Poisson(lambda);
            double[] moments = momentsOfPoisson(lambda, true, 400000);
            double[] other = momentsOfPoisson(lambda, false, 400000);
            // measured at the crossover: means within 0.02 of each other and of
            // the truth, variances within 0.05
            assertEquals("lambda=" + lambda + ": the multiplied mean", law.mean(), moments[0], 0.05);
            assertEquals("lambda=" + lambda + ": the rejected mean", law.mean(), other[0], 0.05);
            assertEquals("lambda=" + lambda + ": the multiplied variance", law.variance(), moments[1], 0.25);
            assertEquals("lambda=" + lambda + ": the rejected variance", law.variance(), other[1], 0.25);
        }
    }

    private static double[] momentsOfPoisson(double lambda, boolean multiply, int draws) {
        PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
        PoissonSpliterator.Rejection rejection = multiply ? null : new PoissonSpliterator.Rejection(lambda);
        double sum = 0.0;
        double sumOfSquares = 0.0;
        for (int i = 0; i < draws; i++) {
            double k = multiply ? PoissonSpliterator.byMultiplication(prng, lambda) : rejection.sample(prng);
            sum += k;
            sumOfSquares += k * k;
        }
        double mean = sum / draws;
        return new double[] { mean, sumOfSquares / draws - mean * mean };
    }

    @Test
    public void testThePoissonEdges() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        for (int i = 0; i < 100; i++) {
            assertEquals("a mean of zero", 0, prng.poisson(1, 0.0).toArray()[0]);
        }
        // a mean of a millionth over a million draws: measured, none was not
        // zero, and the chance of one is about one in three
        int nonZero = 0;
        int[] tiny = prng.poisson(1000000, 1.0e-6).toArray();
        for (int i = 0; i < tiny.length; i++) {
            if (tiny[i] != 0) {
                nonZero++;
            }
        }
        assertTrue("a mean of 1e-6 produced " + nonZero + " non-zero counts", nonZero <= 20);

        // the largest mean offered still fits an int with room to spare
        int[] huge = prng.poisson(200, PoissonSpliterator.MAX_MEAN).toArray();
        for (int i = 0; i < huge.length; i++) {
            assertTrue("a count of " + huge[i] + " is not near a billion",
                    Math.abs(huge[i] - 1.0e9) < 12.0 * Math.sqrt(1.0e9));
        }
    }

    // ------------------------------------------------------------ binomial --

    @Test
    public void testBothBinomialBranchesFitTheMassFunction() {
        double[][] shapes = { { 10, 0.5 }, { 100, 0.05 }, { 1000, 0.01 }, { 50, 0.3 }, { 1000, 0.5 },
                { 100000, 0.001 }, { 20, 0.9 } };
        for (double[] shape : shapes) {
            final int n = (int) shape[0];
            final double p = shape[1];
            final boolean reflect = p > 0.5;
            final double q = reflect ? 1.0 - p : p;
            Binomial law = new Binomial(n, p);

            double byInversion = uniformityOfFits(law, new Draw() {
                @Override
                public int next(PseudoRandom prng) {
                    int k = BinomialSpliterator.byInversion(prng, n, q);
                    return reflect ? n - k : k;
                }
            }, 120, 40000);
            // measured: 0.0728 .. 0.8578 over seven shapes
            assertTrue("n=" + n + " p=" + p + ": inversion, uniformity p = " + byInversion,
                    byInversion > 0.001);

            if (n * q >= 10.0) {
                double byRejection = uniformityOfFits(law, new Draw() {
                    @Override
                    public int next(PseudoRandom prng) {
                        int k = BinomialSpliterator.byTransformedRejection(prng, n, q);
                        return reflect ? n - k : k;
                    }
                }, 120, 40000);
                // measured: 0.1176 .. 0.6521, and the one case that came back
                // at 0.005 stayed above 0.45 over three further seed offsets
                assertTrue("n=" + n + " p=" + p + ": rejection, uniformity p = " + byRejection,
                        byRejection > 0.001);
            }
        }
    }

    @Test
    public void testTheBinomialBranchesAgreeAtTheCrossover() {
        // n p sits just below and just above the crossover
        double p = 0.25;
        for (double target : new double[] { BinomialSpliterator.REJECTION_THRESHOLD - 1.0,
                BinomialSpliterator.REJECTION_THRESHOLD + 1.0 }) {
            int n = (int) Math.round(target / p);
            Binomial law = new Binomial(n, p);
            double[] inversion = momentsOfBinomial(n, p, true, 400000);
            double[] rejection = momentsOfBinomial(n, p, false, 400000);
            String at = "n=" + n + " p=" + p;
            assertEquals(at + ": the inverted mean", law.mean(), inversion[0], 0.05);
            assertEquals(at + ": the rejected mean", law.mean(), rejection[0], 0.05);
            assertEquals(at + ": the inverted variance", law.variance(), inversion[1], 0.3);
            assertEquals(at + ": the rejected variance", law.variance(), rejection[1], 0.3);
        }
    }

    private static double[] momentsOfBinomial(int n, double p, boolean inversion, int draws) {
        PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
        BinomialSpliterator.Rejection rejection = inversion ? null
                : new BinomialSpliterator.Rejection(n, p);
        double sum = 0.0;
        double sumOfSquares = 0.0;
        for (int i = 0; i < draws; i++) {
            double k = inversion ? BinomialSpliterator.byInversion(prng, n, p) : rejection.sample(prng);
            sum += k;
            sumOfSquares += k * k;
        }
        double mean = sum / draws;
        return new double[] { mean, sumOfSquares / draws - mean * mean };
    }

    @Test
    public void testTheBinomialIsSymmetricUnderReflection() {
        // counting the failures of p is counting the successes of 1 - p, which
        // is how only half the parameter range is implemented
        for (int n : new int[] { 7, 60, 400 }) {
            for (double p : new double[] { 0.05, 0.25, 0.5 }) {
                Binomial law = new Binomial(n, 1.0 - p);
                // the sized overload: binomial(n, p) would be an endless
                // stream, which toArray then tries to hold all of
                double[] moments = momentsOfSample(
                        DefaultRng.newPseudoRandom(17L).binomial(200000L, n, 1.0 - p).toArray());
                String at = "n=" + n + " p=" + (1.0 - p);
                assertEquals(at + ": the mean", law.mean(), moments[0], 0.02 * n + 0.05);
                assertEquals(at + ": the variance", law.variance(), moments[1], 0.05 * n + 0.5);
            }
        }
    }

    private static double[] momentsOfSample(int[] drawn) {
        double sum = 0.0;
        double sumOfSquares = 0.0;
        for (int i = 0; i < drawn.length; i++) {
            sum += drawn[i];
            sumOfSquares += (double) drawn[i] * drawn[i];
        }
        double mean = sum / drawn.length;
        return new double[] { mean, sumOfSquares / drawn.length - mean * mean };
    }

    @Test
    public void testTheBinomialEdges() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        assertEquals("no trials", 0, prng.binomial(1, 0, 0.5).toArray()[0]);
        int[] never = prng.binomial(100, 10, 0.0).toArray();
        int[] always = prng.binomial(100, 10, 1.0).toArray();
        for (int i = 0; i < never.length; i++) {
            assertEquals("p = 0", 0, never[i]);
            assertEquals("p = 1", 10, always[i]);
        }
        // one trial is a coin, and the two outcomes are the only ones
        int ones = 0;
        int[] coin = prng.binomial(400000, 1, 0.5).toArray();
        for (int i = 0; i < coin.length; i++) {
            assertTrue("a single trial gave " + coin[i], coin[i] == 0 || coin[i] == 1);
            ones += coin[i];
        }
        // measured over a million draws: 0.49933
        assertEquals("a fair coin", 0.5, ones / (double) coin.length, 0.01);
    }

    // ------------------------------------------------------------ the wiring --

    @Test
    public void testTheStreamsAreReproducibleAndSplit() {
        for (int form = 0; form < 2; form++) {
            final boolean poisson = form == 0;
            int[] first = draw(DefaultRng.newPseudoRandom(7L), poisson);
            int[] again = draw(DefaultRng.newPseudoRandom(7L), poisson);
            assertEquals("length", first.length, again.length);
            for (int i = 0; i < first.length; i++) {
                assertEquals((poisson ? "poisson" : "binomial") + " at " + i, first[i], again[i]);
            }
        }
        assertTrue("the Poisson spliterator refused to split",
                new PoissonSpliterator(DefaultRng.newPseudoRandom(1L), 0L, 1024L, 12.0).trySplit() != null);
        assertTrue("the binomial spliterator refused to split",
                new BinomialSpliterator(DefaultRng.newPseudoRandom(1L), 0L, 1024L, 40, 0.3)
                        .trySplit() != null);
    }

    private static int[] draw(PseudoRandom prng, boolean poisson) {
        return poisson ? prng.poisson(200, 4.5).toArray() : prng.binomial(200, 30, 0.4).toArray();
    }

    @Test
    public void testEveryGeneratorCanDrawBothCounts() {
        PseudoRandom[] generators = { new Interleaved4Stc64(1L), new Lcg64Xor1024Mix(2L),
                new MersenneTwister64(3L), new Sfc64(4L), new SplitMix64(5L), new Stc64(6L),
                new XorShift1024Star(7L), new XorShift128Plus(8L), new XorShiftRot256StarStar(9L) };
        Poisson poissonLaw = new Poisson(8.0);
        Binomial binomialLaw = new Binomial(60, 0.3);
        for (int g = 0; g < generators.length; g++) {
            String name = generators[g].getAlgorithm();
            double[] poisson = momentsOfSample(generators[g].poisson(200000, 8.0).toArray());
            assertEquals(name + ": the Poisson mean", poissonLaw.mean(), poisson[0], 0.1);
            assertEquals(name + ": the Poisson variance", poissonLaw.variance(), poisson[1], 0.3);
            double[] binomial = momentsOfSample(generators[g].binomial(200000, 60, 0.3).toArray());
            assertEquals(name + ": the binomial mean", binomialLaw.mean(), binomial[0], 0.1);
            assertEquals(name + ": the binomial variance", binomialLaw.variance(), binomial[1], 0.3);
        }
    }

    @Test
    public void testTheCountSamplersRejectWhatTheyCannotDraw() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        rejectsPoisson("a negative mean", prng, -1.0);
        rejectsPoisson("a mean that is not a number", prng, Double.NaN);
        rejectsPoisson("a mean beyond an int", prng, 2.0 * PoissonSpliterator.MAX_MEAN);
        rejectsPoisson("an infinite mean", prng, Double.POSITIVE_INFINITY);
        rejectsBinomial("a negative number of trials", prng, -1, 0.5);
        rejectsBinomial("a probability below zero", prng, 10, -0.01);
        rejectsBinomial("a probability above one", prng, 10, 1.01);
        rejectsBinomial("a probability that is not a number", prng, 10, Double.NaN);
        try {
            prng.poisson(-1L, 1.0);
            fail("poisson accepted a negative stream size");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            prng.binomial(-1L, 10, 0.5);
            fail("binomial accepted a negative stream size");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsPoisson(String what, PseudoRandom prng, double lambda) {
        try {
            prng.poisson(lambda);
            fail("poisson accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsBinomial(String what, PseudoRandom prng, int n, double p) {
        try {
            prng.binomial(n, p);
            fail("binomial accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
