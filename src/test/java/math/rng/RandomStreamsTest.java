package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.Test;

import math.cern.ProbabilityFuncs;

/**
 * The spliterator layer: the bounds the streams promise, the sizes they
 * deliver, and the splitting that {@code .parallel()} triggers.
 * <p>
 * Splitting used to hand the same generator to both halves. Since a generator
 * of this package is not thread-safe, a parallel stream was a data race with
 * observable duplicates -- 102723 of 200000 values on
 * {@code longs(200000).parallel()}. {@code trySplit} asks the source for an
 * independent generator now and gives up rather than share one.
 */
public class RandomStreamsTest {

    /**
     * The count has to be exact and the values distinct. A parallel result is
     * not the sequential one and is not portable across machines -- the split
     * tree depends on {@code ForkJoinPool.common.parallelism} -- which is why
     * this asserts structure rather than a value.
     */
    @Test
    public void testParallelStreamsDoNotShareAGenerator() {
        long[] v = new Stc64(42L).longs(200000L).parallel().toArray();
        assertEquals("stream size", 200000, v.length);
        long[] sorted = v.clone();
        Arrays.sort(sorted);
        int duplicates = 0;
        for (int i = 1; i < sorted.length; ++i) {
            if (sorted[i] == sorted[i - 1]) {
                ++duplicates;
            }
        }
        assertEquals("duplicate values in a parallel stream", 0, duplicates);

        double[] d = new Stc64(7L).doubles(200000L).parallel().toArray();
        assertEquals("stream size", 200000, d.length);
        double[] ds = d.clone();
        Arrays.sort(ds);
        int dd = 0;
        for (int i = 1; i < ds.length; ++i) {
            if (ds[i] == ds[i - 1]) {
                ++dd;
            }
        }
        assertEquals("duplicate values in a parallel double stream", 0, dd);
    }

    /** A sequential stream is reproducible from its seed, and stays so. */
    @Test
    public void testSequentialStreamsAreReproducible() {
        long[] a = new Stc64(42L).longs(50000L).toArray();
        long[] b = new Stc64(42L).longs(50000L).toArray();
        assertTrue("the same seed gave a different stream", Arrays.equals(a, b));
        double[] c = new Stc64(42L).normal(50000L, 1.0, 2.0).toArray();
        double[] e = new Stc64(42L).normal(50000L, 1.0, 2.0).toArray();
        assertTrue("the same seed gave a different stream", Arrays.equals(c, e));
    }

    /** Every sized stream method has to deliver exactly what was asked for. */
    @Test
    public void testStreamSizesAreExact() {
        PseudoRandom prng = new Stc64(42L);
        assertEquals("ints", 1000, prng.ints(1000L).toArray().length);
        assertEquals("ints bounded", 1000, prng.ints(1000L, 0, 9).toArray().length);
        assertEquals("longs", 1000, prng.longs(1000L).toArray().length);
        assertEquals("longs bounded", 1000, prng.longs(1000L, -5L, 5L).toArray().length);
        assertEquals("doubles", 1000, prng.doubles(1000L).toArray().length);
        assertEquals("doubles bounded", 1000, prng.doubles(1000L, -1.0, 1.0).toArray().length);
        assertEquals("normal", 1000, prng.normal(1000L, 0.0, 1.0).toArray().length);
        assertEquals("gamma", 1000, prng.gamma(1000L, 2.0, 1.0).toArray().length);
        assertEquals("beta", 1000, prng.beta(1000L, 2.0, 3.0).toArray().length);
        assertEquals("chiSquare", 1000, prng.chiSquare(1000L, 3.0).toArray().length);
        assertEquals("studentT", 1000, prng.studentT(1000L, 5.0).toArray().length);
        assertEquals("fisherF", 1000, prng.fisherF(1000L, 3, 5).toArray().length);
        assertEquals("cauchy", 1000, prng.cauchy(1000L, 0.0, 1.0).toArray().length);
        assertEquals("exponential", 1000, prng.exponential(1000L, 1.5).toArray().length);
        assertEquals("logNormal", 1000, prng.logNormal(1000L, 0.0, 1.0).toArray().length);
        assertEquals("an empty stream", 0, prng.longs(0L).toArray().length);
    }

    /**
     * The integral streams draw through {@code nextInt(min, max)} and
     * {@code nextLong(min, max)}, whose range includes both ends; the double
     * streams draw through {@code nextDouble(min, max)}, which is half-open.
     */
    @Test
    public void testIntegralStreamsIncludeBothEndsAndDoublesDoNot() {
        int[] counts = new int[4];
        int[] v = new Stc64(9L).ints(100000L, 0, 3).toArray();
        for (int i = 0; i < v.length; ++i) {
            assertTrue("out of range: " + v[i], v[i] >= 0 && v[i] <= 3);
            ++counts[v[i]];
        }
        for (int i = 0; i < counts.length; ++i) {
            assertTrue("value " + i + " never drawn", counts[i] > 0);
        }
        long lo = Long.MAX_VALUE;
        long hi = Long.MIN_VALUE;
        long[] w = new Stc64(9L).longs(100000L, -2L, 2L).toArray();
        for (int i = 0; i < w.length; ++i) {
            assertTrue("out of range: " + w[i], w[i] >= -2L && w[i] <= 2L);
            lo = Math.min(lo, w[i]);
            hi = Math.max(hi, w[i]);
        }
        assertEquals("lower end never drawn", -2L, lo);
        assertEquals("upper end never drawn", 2L, hi);
        double[] d = new Stc64(9L).doubles(200000L, -1.0, 1.0).toArray();
        for (int i = 0; i < d.length; ++i) {
            assertTrue("out of range: " + d[i], d[i] >= -1.0 && d[i] < 1.0);
        }
    }

    /**
     * The round trip through the error function and its inverse is not exact,
     * so a bound the sampler cannot resolve used to leave the interval
     * entirely: on {@code (8.3, 8.5)} every one of 20000 draws landed outside
     * and 9867 of them were {@code +Infinity}, because the error function has
     * saturated to one there and its inverse answers infinity. An equality
     * test caught neither, so the clamp compares with {@code <=} and
     * {@code >=} now.
     * <p>
     * What the clamp restores is the range, not the resolution: in the
     * saturated window every value comes back as the same clamped bound. That
     * is a limit of the inverse error function, and a caller asking for a
     * window past 8.3 sigma is asking for something a double cannot resolve.
     */
    @Test
    public void testTruncatedNormalStaysInsideItsBounds() {
        double[][] windows = { { -3.0, 3.0 }, { -8.5, 8.5 }, { 5.0, 5.2 }, { 8.2, 8.4 }, { 8.3, 8.5 },
                { -8.5, -8.3 } };
        for (int k = 0; k < windows.length; ++k) {
            double lower = windows[k][0];
            double upper = windows[k][1];
            double[] v = new Stc64(42L).truncatedStandardNormal(20000L, lower, upper).toArray();
            for (int i = 0; i < v.length; ++i) {
                assertTrue("(" + lower + ", " + upper + ") produced " + v[i],
                        v[i] > lower && v[i] < upper);
            }
        }
    }

    /**
     * The Student-t sampler rejects {@code q == 0} and {@code q == 1} now, the
     * same acceptance region {@code AbstractRng64.nextGaussian()} uses for the
     * same polar method. Neither endpoint is reachable by sampling -- both
     * need two consecutive {@code nextDouble()} values to be exactly one half,
     * a 2^-106 event -- so this asserts the contract rather than the repair:
     * finite values, and the distribution they are supposed to have.
     * <p>
     * The KS bound was calibrated over 16 seeds at 50000 draws, where the worst
     * statistic was 0.00635 against a 5 % critical value of 0.00607.
     */
    @Test
    public void testStudentTIsFiniteAndCorrectlyDistributed() {
        double[] dfs = { 1.0, 3.0, 8.0, 30.0 };
        for (int k = 0; k < dfs.length; ++k) {
            double df = dfs[k];
            double[] v = new Stc64(500L).studentT(50000L, df).toArray();
            for (int i = 0; i < v.length; ++i) {
                assertTrue("df " + df + " produced " + v[i], !Double.isNaN(v[i]) && !Double.isInfinite(v[i]));
            }
            Arrays.sort(v);
            double ks = 0.0;
            for (int i = 0; i < v.length; ++i) {
                double f = ProbabilityFuncs.studentT(df, v[i]);
                ks = Math.max(ks, Math.max(Math.abs((i + 1.0) / v.length - f), Math.abs(f - i / (double) v.length)));
            }
            assertTrue("df " + df + " KS statistic " + ks, ks < 0.01);
        }
    }

    /**
     * {@code beta} and {@code fisherF} draw from two generators and take the
     * second from {@code detach}, which is the path {@code trySplit} shares. A
     * shape small enough to underflow used to make the ratio {@code NaN}; the
     * combination is done in log space now.
     */
    @Test
    public void testBetaAndFisherFAreFiniteForSmallShapes() {
        double[] shapes = { 0.001, 0.005, 0.01, 0.5, 2.0 };
        for (int k = 0; k < shapes.length; ++k) {
            double s = shapes[k];
            double[] v = new Stc64(9L).beta(100000L, s, s).toArray();
            for (int i = 0; i < v.length; ++i) {
                assertTrue("beta(" + s + ", " + s + ") produced " + v[i], !Double.isNaN(v[i]));
                assertTrue("beta out of [0, 1]: " + v[i], v[i] >= 0.0 && v[i] <= 1.0);
            }
        }
        double[] f = new Stc64(9L).fisherF(100000L, 3, 5).toArray();
        for (int i = 0; i < f.length; ++i) {
            assertTrue("fisherF produced " + f[i], !Double.isNaN(f[i]) && f[i] >= 0.0);
        }
    }

    /**
     * {@code intsSampledWithoutReplacement} has to return distinct values
     * inside the range, and reject what it cannot deliver.
     */
    @Test
    public void testSamplingWithoutReplacementIsDistinct() {
        PseudoRandom prng = new Stc64(42L);
        int[] v = prng.intsSampledWithoutReplacement(10, 59, 50);
        assertEquals("count", 50, v.length);
        int[] sorted = v.clone();
        Arrays.sort(sorted);
        for (int i = 0; i < sorted.length; ++i) {
            assertTrue("out of range: " + sorted[i], sorted[i] >= 10 && sorted[i] <= 59);
            if (i > 0) {
                assertTrue("duplicate " + sorted[i], sorted[i] != sorted[i - 1]);
            }
        }
        try {
            prng.intsSampledWithoutReplacement(0, 9, 11);
            org.junit.Assert.fail("asked for more distinct values than the range holds");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
    }
}
