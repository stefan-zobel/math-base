package math.dl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Test;

import math.rng.PseudoRandom;
import math.rng.Stc64;

/**
 * {@link Softmax#sampleClass(double[])} drew from a static {@code Stc64} shared
 * by every caller. A sampling method has nowhere to keep state, so the
 * generator now comes either from the caller, through the overloads added for
 * it, or from one instance per thread.
 * <p>
 * The tolerance on the empirical frequencies was calibrated over 24 seeds at
 * 200000 draws; the worst deviation seen was 0.0027.
 */
public class SoftmaxTest {

    private static final double[] P = { 0.1, 0.2, 0.3, 0.4 };

    private static final float[] PF = { 0.1f, 0.2f, 0.3f, 0.4f };

    private static final int N = 200000;

    /** With an explicit generator the sequence of indices is reproducible. */
    @Test
    public void testAnExplicitGeneratorReproducesTheSameIndices() {
        PseudoRandom first = new Stc64(20260823L);
        PseudoRandom second = new Stc64(20260823L);
        for (int i = 0; i < 5000; ++i) {
            assertEquals("draw " + i, Softmax.sampleClass(P, first), Softmax.sampleClass(P, second));
        }
        PseudoRandom firstF = new Stc64(20260824L);
        PseudoRandom secondF = new Stc64(20260824L);
        for (int i = 0; i < 5000; ++i) {
            assertEquals("draw " + i, Softmax.sampleClassF(PF, firstF), Softmax.sampleClassF(PF, secondF));
        }
    }

    /**
     * The roulette wheel has to hit each class at its own probability. This is
     * the guard that routing the draw through a parameter did not change what
     * is sampled.
     */
    @Test
    public void testTheEmpiricalFrequenciesMatchTheDistribution() {
        PseudoRandom prng = new Stc64(20260825L);
        int[] count = new int[P.length];
        for (int i = 0; i < N; ++i) {
            ++count[Softmax.sampleClass(P, prng)];
        }
        for (int i = 0; i < P.length; ++i) {
            assertEquals("class " + i, P[i], count[i] / (double) N, 0.01);
        }
        PseudoRandom prngF = new Stc64(20260826L);
        int[] countF = new int[PF.length];
        for (int i = 0; i < N; ++i) {
            ++countF[Softmax.sampleClassF(PF, prngF)];
        }
        for (int i = 0; i < PF.length; ++i) {
            assertEquals("class " + i, PF[i], countF[i] / (double) N, 0.01);
        }
    }

    /**
     * The convenience form without a generator has to survive being called
     * from several threads at once, and still sample the right distribution
     * when the counts of all of them are put together.
     */
    @Test
    public void testTheConvenienceFormIsSafeAcrossThreads() throws Exception {
        final int threads = 8;
        final int perThread = 50000;
        final int[][] count = new int[threads][P.length];
        final AtomicReference<Throwable> failure = new AtomicReference<Throwable>();
        final CountDownLatch start = new CountDownLatch(1);
        Thread[] t = new Thread[threads];
        for (int i = 0; i < threads; ++i) {
            final int me = i;
            t[i] = new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        start.await();
                        for (int k = 0; k < perThread; ++k) {
                            ++count[me][Softmax.sampleClass(P)];
                        }
                    } catch (Throwable e) {
                        failure.compareAndSet(null, e);
                    }
                }
            });
            t[i].start();
        }
        start.countDown();
        for (int i = 0; i < threads; ++i) {
            t[i].join();
        }
        assertEquals("sampleClass threw " + failure.get(), null, failure.get());
        int total = threads * perThread;
        for (int c = 0; c < P.length; ++c) {
            int sum = 0;
            for (int i = 0; i < threads; ++i) {
                sum += count[i][c];
            }
            assertEquals("class " + c, P[c], sum / (double) total, 0.01);
        }
    }

    /** Every index the roulette wheel returns has to be a class. */
    @Test
    public void testTheIndexIsAlwaysInRange() {
        PseudoRandom prng = new Stc64(20260827L);
        for (int i = 0; i < 20000; ++i) {
            int c = Softmax.sampleClass(P, prng);
            assertTrue("index " + c, c >= 0 && c < P.length);
            int cf = Softmax.sampleClassF(PF, prng);
            assertTrue("index " + cf, cf >= 0 && cf < PF.length);
        }
    }
}
