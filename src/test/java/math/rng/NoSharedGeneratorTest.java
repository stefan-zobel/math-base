package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Test;

/**
 * Eleven classes of this package used to hand out a static singleton through
 * {@code getDefault()}, and none of the generators is thread-safe. Sharing one
 * cost more than statistical quality: eight threads drawing 200000
 * {@code nextLong()} each off the same instance produced 79.6 % duplicate
 * values on {@code MersenneTwister64} -- the instance
 * {@code DefaultRng.getGlobalPseudoRandom()} handed out -- along with 120
 * {@code ArrayIndexOutOfBoundsException} out of {@code mt[mti++]}, whose index
 * is re-read after the {@code mti >= NN} check. {@code nextGaussian()} answered
 * {@code NaN} 752 times, because the cache in {@link AbstractRng64} is a plain
 * field. Even {@code getDefault().split()} was no way out: 418 of 32000
 * children collided on their first draw, since {@code split()} mutates the
 * parent it is called on.
 * <p>
 * A generator is now something a caller owns. The seed comes from
 * {@link SplitMix64Seed}, which is synchronized, and no class of this package
 * keeps a generator in a static field any more.
 */
public class NoSharedGeneratorTest {

    /** the classes that used to carry a {@code defaultRng} singleton */
    private static final Class<?>[] GENERATORS = { Interleaved4Stc64.class, Lcg64Xor1024Mix.class,
            MersenneTwister64.class, Sfc64.class, SplitMix64.class, Stc64.class, XorShift1024Star.class,
            XorShift128Plus.class, XorShift64Star.class, XorShiftRot256PlusPlus.class,
            XorShiftRot256StarStar.class };

    private static final int THREADS = 8;
    private static final int PER_THREAD = 50000;

    /**
     * Eight threads, each building its own generator and drawing from it. The
     * draws have to be free of duplicates and free of exceptions.
     * <p>
     * The generators are deliberately built without a seed -- that is the case
     * under test, since the seed then comes from the shared
     * {@link SplitMix64Seed}. The run is therefore not reproducible, but the
     * assertion is not delicate: 400000 draws out of 2^64 collide with
     * probability about 4e-9 when the streams are independent, while the
     * shared singleton produced between 5.7 % and 82.3 % duplicates depending
     * on the generator.
     */
    @Test
    public void testConcurrentConstructionYieldsDistinctStreams() throws Exception {
        for (int i = 0; i < GENERATORS.length; ++i) {
            final Class<?> c = GENERATORS[i];
            final long[][] out = new long[THREADS][PER_THREAD];
            final AtomicReference<Throwable> failure = new AtomicReference<Throwable>();
            final CountDownLatch start = new CountDownLatch(1);
            Thread[] threads = new Thread[THREADS];
            for (int t = 0; t < THREADS; ++t) {
                final int me = t;
                threads[t] = new Thread(new Runnable() {
                    @Override
                    public void run() {
                        try {
                            PseudoRandom prng = (PseudoRandom) c.getDeclaredConstructor().newInstance();
                            start.await();
                            for (int k = 0; k < PER_THREAD; ++k) {
                                out[me][k] = prng.nextLong();
                            }
                        } catch (Throwable e) {
                            failure.compareAndSet(null, e);
                        }
                    }
                });
                threads[t].start();
            }
            start.countDown();
            for (int t = 0; t < THREADS; ++t) {
                threads[t].join();
            }
            assertEquals(c.getSimpleName() + " threw " + failure.get(), null, failure.get());
            assertEquals(c.getSimpleName() + " produced duplicate values", 0, duplicates(out));
        }
    }

    /**
     * The seed source is the one piece of shared mutable state left, and it is
     * synchronized. Its Weyl sequence never repeats a value inside its period,
     * so 800_000 concurrent draws must all be distinct.
     */
    @Test
    public void testConcurrentSeedingIsUnique() throws Exception {
        final long[][] out = new long[THREADS][100000];
        final AtomicReference<Throwable> failure = new AtomicReference<Throwable>();
        final CountDownLatch start = new CountDownLatch(1);
        Thread[] threads = new Thread[THREADS];
        for (int t = 0; t < THREADS; ++t) {
            final int me = t;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        start.await();
                        for (int k = 0; k < out[me].length; ++k) {
                            out[me][k] = SplitMix64Seed.seed();
                        }
                    } catch (Throwable e) {
                        failure.compareAndSet(null, e);
                    }
                }
            });
            threads[t].start();
        }
        start.countDown();
        for (int t = 0; t < THREADS; ++t) {
            threads[t].join();
        }
        assertEquals("seeding threw " + failure.get(), null, failure.get());
        assertEquals("SplitMix64Seed handed out the same seed twice", 0, duplicates(out));
    }

    /**
     * A sampler that owns its generator is not thread-safe either, and the
     * gaussian cache is where that shows first: eight threads on one instance
     * produced 752 {@code NaN} in 1.6 million draws. One generator per thread
     * must not.
     */
    @Test
    public void testGaussiansAreFiniteWithOneGeneratorPerThread() throws Exception {
        final int n = 100000;
        final double[][] out = new double[THREADS][n];
        final CountDownLatch start = new CountDownLatch(1);
        Thread[] threads = new Thread[THREADS];
        for (int t = 0; t < THREADS; ++t) {
            final int me = t;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    PseudoRandom prng = new Stc64();
                    try {
                        start.await();
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        return;
                    }
                    for (int k = 0; k < n; ++k) {
                        out[me][k] = prng.nextGaussian();
                    }
                }
            });
            threads[t].start();
        }
        start.countDown();
        for (int t = 0; t < THREADS; ++t) {
            threads[t].join();
        }
        double sum = 0.0;
        double sumSq = 0.0;
        for (int t = 0; t < THREADS; ++t) {
            for (int k = 0; k < n; ++k) {
                double v = out[t][k];
                if (Double.isNaN(v)) {
                    fail("nextGaussian() answered NaN at thread " + t + ", draw " + k);
                }
                sum += v;
                sumSq += v * v;
            }
        }
        int total = THREADS * n;
        double mean = sum / total;
        double sd = Math.sqrt(sumSq / total - mean * mean);
        assertTrue("mean " + mean, Math.abs(mean) < 0.02);
        assertTrue("standard deviation " + sd, Math.abs(sd - 1.0) < 0.02);
    }

    private static int duplicates(long[][] out) {
        int total = 0;
        for (int t = 0; t < out.length; ++t) {
            total += out[t].length;
        }
        long[] all = new long[total];
        int at = 0;
        for (int t = 0; t < out.length; ++t) {
            System.arraycopy(out[t], 0, all, at, out[t].length);
            at += out[t].length;
        }
        Arrays.sort(all);
        int dup = 0;
        for (int i = 1; i < all.length; ++i) {
            if (all[i] == all[i - 1]) {
                ++dup;
            }
        }
        return dup;
    }
}
