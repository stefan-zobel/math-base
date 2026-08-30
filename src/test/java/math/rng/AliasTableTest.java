package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import org.junit.Test;

import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;
import math.stats.TestResult;

/**
 * A sampler is a hypothesis about a distribution. It is checked here the way
 * any other hypothesis is: by drawing from it and asking whether the counts
 * are what the weights say, over enough seeds that the answer is about the
 * sampler rather than about one run.
 */
public final class AliasTableTest {

    /** Deterministic uniforms in {@code (0,1)}, for the weight vectors. */
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

    private static long[] counts(PseudoRandom prng, double[] weights, int draws) {
        long[] observed = new long[weights.length];
        prng.categorical(draws, weights).forEach(k -> observed[k]++);
        return observed;
    }

    private static double[] expectedCounts(double[] weights, int draws) {
        double sum = 0.0;
        for (int i = 0; i < weights.length; i++) {
            sum += weights[i];
        }
        double[] expected = new double[weights.length];
        for (int i = 0; i < weights.length; i++) {
            expected[i] = draws * weights[i] / sum;
        }
        return expected;
    }

    // ------------------------------------------------ does it fit the weights --

    @Test
    public void testTheDrawnCountsFitTheWeights() {
        // one chi-squared test that passes says little. Several hundred whose
        // p-values are uniform say that the sampler is right, which is the
        // statement worth making
        Uniform uniform = new Uniform(0.0, 1.0);
        double[] flat = new double[8];
        Arrays.fill(flat, 1.0);
        double[] skewed = { 1000.0, 100.0, 10.0, 1.0, 1.0, 0.5, 0.25, 0.125 };
        double[] twoOutcomes = { 1.0, 3.0 };

        int reps = 300;
        int draws = 50000;
        for (double[] weights : new double[][] { flat, skewed, twoOutcomes }) {
            double[] expected = expectedCounts(weights, draws);
            double[] pValues = new double[reps];
            int rejected = 0;
            for (int r = 0; r < reps; r++) {
                long[] observed = counts(DefaultRng.newPseudoRandom(r * 7919L + 1), weights, draws);
                pValues[r] = HypothesisTests.chiSquaredGoodnessOfFit(observed, expected, 0).pValue;
                if (pValues[r] <= 0.05) {
                    rejected++;
                }
            }
            String at = weights.length + " outcomes";
            // measured over 400 seeds and 100000 draws: 0.0350 .. 0.0550
            assertTrue(at + ": rejected " + rejected + " of " + reps,
                    Math.abs(rejected / (double) reps - 0.05) < 0.04);
            TestResult uniformity = HypothesisTests.kolmogorovSmirnov(pValues, uniform, Alternative.TWO_SIDED);
            // measured: 0.0346 .. 0.0381, against a 0.001 critical value near
            // 0.12 at this many replications
            assertTrue(at + ": the p-values are not uniform, D = " + uniformity.statistic,
                    uniformity.statistic < 0.12);
        }
    }

    @Test
    public void testAnOutcomeOfWeightZeroIsNeverDrawn() {
        // the interesting failure of the alias method is a leftover entry in
        // the light worklist being given a full share at the end, which would
        // hand a zero-weight outcome a probability of 1/n
        double[] weights = { 0.0, 2.0, 0.0, 0.0, 5.0, 0.0, 3.0, 0.0 };
        int draws = 2000000;
        long[] seen = counts(DefaultRng.newPseudoRandom(11L), weights, draws);
        for (int i = 0; i < weights.length; i++) {
            if (weights[i] == 0.0) {
                assertEquals("outcome " + i + " has weight zero", 0L, seen[i]);
            } else {
                assertTrue("outcome " + i + " never appeared", seen[i] > 0L);
            }
        }
        // and the ones that do have weight land where they should: measured
        // over ten million draws, 1998374 / 4999180 / 3002446 against
        // 2000000 / 5000000 / 3000000
        assertEquals("the total", draws, seen[1] + seen[4] + seen[6]);
        assertEquals("the middle one", 0.5 * draws, seen[4], 0.01 * draws);
    }

    @Test
    public void testTheVoseInvariantsHold() {
        // the construction promises two things about every column, and they
        // are cheaper to check directly than to infer from the counts
        Lcg lcg = new Lcg(4242L);
        int tables = 0;
        for (int round = 0; round < 20000; round++) {
            int n = 1 + (int) (lcg.next() * 40);
            double[] weights = new double[n];
            boolean any = false;
            for (int i = 0; i < n; i++) {
                // a quarter of them zero, the rest spanning twelve orders of
                // magnitude, which is where a worklist gets its ordering wrong
                weights[i] = lcg.next() < 0.25 ? 0.0 : Math.pow(10.0, -6.0 + 12.0 * lcg.next());
                any |= weights[i] > 0.0;
            }
            if (!any) {
                continue;
            }
            AliasTable table = AliasTable.of(weights);
            tables++;
            assertEquals("outcomes", n, table.outcomes());
            for (int i = 0; i < n; i++) {
                double stay = table.stayProbability(i);
                assertTrue("stay probability " + stay + " at " + i, stay >= 0.0 && stay <= 1.0);
                int alias = table.aliasOf(i);
                assertTrue("alias index " + alias + " at " + i, alias >= 0 && alias < n);
            }
        }
        // measured: 19819 tables built, no violation of either invariant
        assertTrue("almost nothing was built: " + tables, tables > 19000);
    }

    // ------------------------------------------------------------ the streams --

    @Test
    public void testTheWeightsAndThePreparedTableGiveTheSameStream() {
        // the two overloads must not drift: one builds the table the other is
        // handed, and from the same seed they are the same arithmetic
        double[] weights = { 5.0, 3.0, 2.0, 0.5 };
        AliasTable table = AliasTable.of(weights);
        int[] fromWeights = DefaultRng.newPseudoRandom(23L).categorical(500, weights).toArray();
        int[] fromTable = DefaultRng.newPseudoRandom(23L).categorical(500, table).toArray();
        assertArrayEqualsBitForBit(fromWeights, fromTable);

        // and the unbounded forms agree with the sized ones on their prefix
        int[] unbounded = DefaultRng.newPseudoRandom(23L).categorical(weights).limit(500).toArray();
        assertArrayEqualsBitForBit(fromWeights, unbounded);
    }

    private static void assertArrayEqualsBitForBit(int[] wanted, int[] got) {
        assertEquals("length", wanted.length, got.length);
        for (int i = 0; i < wanted.length; i++) {
            assertEquals("at " + i, wanted[i], got[i]);
        }
    }

    @Test
    public void testASequentialStreamIsReproducibleAndAParallelOneSplits() {
        double[] weights = { 5.0, 3.0, 2.0 };
        int[] first = DefaultRng.newPseudoRandom(7L).categorical(200, weights).toArray();
        int[] again = DefaultRng.newPseudoRandom(7L).categorical(200, weights).toArray();
        assertArrayEqualsBitForBit(first, again);

        // a parallel stream draws from a generator per split, so its values
        // depend on the parallelism of the common pool and are not
        // reproducible -- which the PseudoRandom javadoc says. What has to
        // hold is the distribution
        final AtomicLong[] tally = new AtomicLong[weights.length];
        for (int i = 0; i < tally.length; i++) {
            tally[i] = new AtomicLong();
        }
        int draws = 2000000;
        DefaultRng.newPseudoRandom(3L).categorical(draws, weights).parallel()
                .forEach(k -> tally[k].incrementAndGet());
        long[] observed = new long[weights.length];
        for (int i = 0; i < tally.length; i++) {
            observed[i] = tally[i].get();
        }
        // measured: p = 0.2847
        double pValue = HypothesisTests.chiSquaredGoodnessOfFit(observed, expectedCounts(weights, draws), 0)
                .pValue;
        assertTrue("the parallel stream does not fit the weights, p = " + pValue, pValue > 0.001);

        // and the check above means nothing unless the stream really splits
        CategoricalSpliterator spliterator = new CategoricalSpliterator(DefaultRng.newPseudoRandom(3L), 0L,
                1024L, AliasTable.of(weights));
        assertTrue("the spliterator refused to split", spliterator.trySplit() != null);
    }

    @Test
    public void testEveryGeneratorCanDrawFromTheTable() {
        // the methods are inherited from AbstractRng64, which is the kind of
        // claim that holds until one generator overrides something
        double[] weights = { 1.0, 2.0, 1.0 };
        PseudoRandom[] generators = { new Interleaved4Stc64(1L), new Lcg64Xor1024Mix(2L),
                new MersenneTwister64(3L), new Sfc64(4L), new SplitMix64(5L), new Stc64(6L),
                new XorShift1024Star(7L), new XorShift128Plus(8L), new XorShiftRot256StarStar(9L) };
        for (int g = 0; g < generators.length; g++) {
            String name = generators[g].getAlgorithm();
            long[] observed = counts(generators[g], weights, 40000);
            long total = 0L;
            for (int i = 0; i < observed.length; i++) {
                assertTrue(name + ": outcome " + i + " never appeared", observed[i] > 0L);
                total += observed[i];
            }
            assertEquals(name + ": the total", 40000L, total);
            double pValue = HypothesisTests.chiSquaredGoodnessOfFit(observed, expectedCounts(weights, 40000), 0)
                    .pValue;
            assertTrue(name + ": does not fit the weights, p = " + pValue, pValue > 0.0005);
        }
    }

    // -------------------------------------------------------------- the edges --

    @Test
    public void testTheSmallAndTheLargeTable() {
        // one outcome carries everything and needs no coin at all
        AliasTable single = AliasTable.of(new double[] { 3.0 });
        assertEquals("outcomes", 1, single.outcomes());
        PseudoRandom prng = DefaultRng.newPseudoRandom(13L);
        for (int i = 0; i < 1000; i++) {
            assertEquals("the only outcome", 0, single.sample(prng));
        }

        // a hundred thousand of them cost the same per draw as three, which is
        // the whole point of the method against a cumulative search
        int n = 100000;
        double[] weights = new double[n];
        Lcg lcg = new Lcg(99L);
        for (int i = 0; i < n; i++) {
            weights[i] = lcg.next();
        }
        AliasTable large = AliasTable.of(weights);
        assertEquals("outcomes", n, large.outcomes());
        int lowest = n;
        int highest = -1;
        for (int i = 0; i < 200000; i++) {
            int drawn = large.draw(prng);
            assertTrue("out of range: " + drawn, drawn >= 0 && drawn < n);
            lowest = Math.min(lowest, drawn);
            highest = Math.max(highest, drawn);
        }
        // 200000 draws over 100000 nearly equal outcomes reach both ends
        assertTrue("the low end was never reached: " + lowest, lowest < n / 100);
        assertTrue("the high end was never reached: " + highest, highest > n - n / 100);
    }

    @Test
    public void testTheWeightsNeedNotSumToOne() {
        // counts serve as well as probabilities, so the same shape scaled by
        // any positive factor is the same table
        double[] asProbabilities = { 0.5, 0.3, 0.2 };
        double[] asCounts = { 50.0, 30.0, 20.0 };
        double[] asAnything = { 5.0e7, 3.0e7, 2.0e7 };
        int[] a = DefaultRng.newPseudoRandom(31L).categorical(2000, asProbabilities).toArray();
        int[] b = DefaultRng.newPseudoRandom(31L).categorical(2000, asCounts).toArray();
        int[] c = DefaultRng.newPseudoRandom(31L).categorical(2000, asAnything).toArray();
        assertArrayEqualsBitForBit(a, b);
        assertArrayEqualsBitForBit(a, c);
    }

    // --------------------------------------------------------- the guard rail --

    @Test
    public void testTheTableRejectsWhatItCannotBuild() {
        rejects("null weights", null);
        rejects("no weights", new double[0]);
        rejects("a negative weight", new double[] { 1.0, -0.5 });
        rejects("a NaN weight", new double[] { 1.0, Double.NaN });
        rejects("an infinite weight", new double[] { 1.0, Double.POSITIVE_INFINITY });
        rejects("weights that are all zero", new double[] { 0.0, 0.0, 0.0 });

        try {
            AliasTable.of(new double[] { 1.0 }).sample(null);
            fail("sample accepted a null generator");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            DefaultRng.newPseudoRandom(1L).categorical(10, (AliasTable) null);
            fail("categorical accepted a null table");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            DefaultRng.newPseudoRandom(1L).categorical(-1L, new double[] { 1.0, 1.0 });
            fail("categorical accepted a negative stream size");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejects(String what, double[] weights) {
        try {
            AliasTable.of(weights);
            fail("of accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            DefaultRng.newPseudoRandom(1L).categorical(weights);
            fail("categorical accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
