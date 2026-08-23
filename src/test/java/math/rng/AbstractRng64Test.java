package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The generator contract of {@link AbstractRng64}, one test per defect that was
 * repaired in version 1.5.2. Every assertion below fails on the tree as it stood
 * before those repairs, and the numbers in the comments are what that tree produced.
 */
public class AbstractRng64Test {

    /**
     * {@code doubles()} and {@code doubles(n)} used to be built on
     * {@code -Double.MAX_VALUE, Double.MAX_VALUE} instead of the unit interval,
     * so they returned values around {@code 1.8e308}.
     */
    @Test
    public void testTheUnboundedDoubleStreamsAreTheUnitInterval() {
        double[] v = new Stc64(42L).doubles(20000L).toArray();
        for (int i = 0; i < v.length; ++i) {
            assertTrue("doubles() returned " + v[i], v[i] >= 0.0 && v[i] < 1.0);
        }
        double[] w = new Stc64(42L).doubles().limit(20000L).toArray();
        for (int i = 0; i < w.length; ++i) {
            assertTrue("doubles() returned " + w[i], w[i] >= 0.0 && w[i] < 1.0);
        }
    }

    /**
     * {@code longs()} and {@code longs(n)} were bounded by
     * {@code Integer.MIN_VALUE, Integer.MAX_VALUE}, so a stream of longs never
     * left the int range -- measured min and max were {@code -2147454105} and
     * {@code 2147480748} over 200000 draws.
     */
    @Test
    public void testTheUnboundedLongStreamsCoverTheLongRange() {
        long[] v = new Stc64(42L).longs(200000L).toArray();
        long lo = Long.MAX_VALUE;
        long hi = Long.MIN_VALUE;
        for (int i = 0; i < v.length; ++i) {
            lo = Math.min(lo, v[i]);
            hi = Math.max(hi, v[i]);
        }
        assertTrue("smallest value was " + lo + ", still inside the int range", lo < Integer.MIN_VALUE);
        assertTrue("largest value was " + hi + ", still inside the int range", hi > Integer.MAX_VALUE);
    }

    /**
     * {@code nextLong(min, max)} computed {@code max - min + 1} and rejected
     * the result when it came out non-positive, which is every span wider than
     * {@code Long.MAX_VALUE} -- the whole range included. It threw
     * {@code IllegalArgumentException}.
     */
    @Test
    public void testTheFullLongRangeIsADrawableRange() {
        PseudoRandom prng = new Stc64(42L);
        long lo = Long.MAX_VALUE;
        long hi = Long.MIN_VALUE;
        for (int i = 0; i < 200000; ++i) {
            long v = prng.nextLong(Long.MIN_VALUE, Long.MAX_VALUE);
            lo = Math.min(lo, v);
            hi = Math.max(hi, v);
        }
        assertTrue("no negative value in 200000 draws", lo < 0L);
        assertTrue("no positive value in 200000 draws", hi > 0L);
        // and a span just below the overflow, which took the ordinary path
        PseudoRandom p2 = new Stc64(43L);
        for (int i = 0; i < 20000; ++i) {
            long v = p2.nextLong(Long.MIN_VALUE + 1L, Long.MAX_VALUE);
            assertTrue("out of range: " + v, v >= Long.MIN_VALUE + 1L);
        }
    }

    /**
     * For a power-of-two bound the shortcut used to return {@code x & (n - 1)},
     * the low bits of the draw, which are the weakest bits of several
     * generators in this package. It takes the high ones now, so the result is
     * exactly {@code nextLong() >>> (64 - k)} -- an exact identity, not a
     * statistical claim. The old form matched it 5 times in 5000 at
     * {@code n = 2^10}.
     */
    @Test
    public void testThePowerOfTwoShortcutTakesTheHighBits() {
        for (int k = 1; k <= 20; ++k) {
            long n = 1L << k;
            PseudoRandom bounded = new Stc64(42L);
            PseudoRandom raw = new Stc64(42L);
            for (int i = 0; i < 2000; ++i) {
                assertEquals("nextLong(2^" + k + ")", raw.nextLong() >>> (64 - k), bounded.nextLong(n));
            }
        }
        // n == 1 is the degenerate power of two: the shift is by 64, which Java
        // takes mod 64 and therefore does not shift at all. The mask is what
        // saves it
        PseudoRandom prng = new Stc64(42L);
        for (int i = 0; i < 20000; ++i) {
            assertEquals("nextLong(1)", 0L, prng.nextLong(1L));
        }
    }

    /**
     * {@code min + (max - min) * nextDouble()} can round up to {@code max}: on
     * the interval {@code [1, nextUp(1))} it did so in 999704 of 2000000 draws,
     * half of them, because every product below one ulp rounds to the same
     * value. The interval is half-open, so the result is stepped back.
     */
    @Test
    public void testBoundedDoublesAndFloatsAreHalfOpen() {
        double max = Math.nextUp(1.0);
        PseudoRandom prng = new Stc64(11L);
        for (int i = 0; i < 500000; ++i) {
            double d = prng.nextDouble(1.0, max);
            assertTrue("nextDouble reached its upper bound", d < max);
            assertTrue("nextDouble below its lower bound", d >= 1.0);
        }
        float fmax = Math.nextUp(1.0f);
        PseudoRandom p2 = new Stc64(12L);
        for (int i = 0; i < 500000; ++i) {
            float f = p2.nextFloat(1.0f, fmax);
            assertTrue("nextFloat reached its upper bound", f < fmax);
            assertTrue("nextFloat below its lower bound", f >= 1.0f);
        }
        // an ordinary interval has to stay half-open as well
        PseudoRandom p3 = new Stc64(13L);
        for (int i = 0; i < 200000; ++i) {
            double d = p3.nextDouble(-1.0, 1.0);
            assertTrue("out of range: " + d, d >= -1.0 && d < 1.0);
        }
    }

    /**
     * An empty range is legal and returns its single value; the half-open rule
     * does not apply to it, since there is nothing below the bound.
     */
    @Test
    public void testAnEmptyRangeReturnsItsOnlyValue() {
        PseudoRandom prng = new Stc64(42L);
        for (int i = 0; i < 1000; ++i) {
            assertEquals("nextDouble(5, 5)", 5.0, prng.nextDouble(5.0, 5.0), 0.0);
            assertEquals("nextFloat(5, 5)", 5.0f, prng.nextFloat(5.0f, 5.0f), 0.0f);
            assertEquals("nextInt(3, 3)", 3, prng.nextInt(3, 3));
            assertEquals("nextLong(-7, -7)", -7L, prng.nextLong(-7L, -7L));
        }
    }

    /**
     * Java takes a shift count modulo 64, so {@code nextLong() >>> (64 - 0)}
     * did not shift and {@code next(0)} returned a full random int instead of
     * zero. Out-of-range bit counts went unnoticed too.
     */
    @Test
    public void testNextBitsIsValidatedAndZeroBitsGiveZero() {
        PseudoRandom prng = new Stc64(42L);
        for (int i = 0; i < 1000; ++i) {
            assertEquals("next(0)", 0, prng.next(0));
        }
        // next(bits) hands out the top `bits` bits of a draw. Up to 31 of them
        // fit below the sign, so the value is non-negative and bounded; 32 of
        // them fill the int and it is an ordinary signed one
        for (int bits = 1; bits <= 31; ++bits) {
            for (int i = 0; i < 200; ++i) {
                int v = prng.next(bits);
                assertTrue("next(" + bits + ") returned " + v, v >= 0 && v < (1L << bits));
            }
        }
        boolean negative = false;
        boolean positive = false;
        for (int i = 0; i < 200; ++i) {
            int v = prng.next(32);
            negative |= v < 0;
            positive |= v > 0;
        }
        assertTrue("next(32) never returned a negative value", negative);
        assertTrue("next(32) never returned a positive value", positive);
        assertThrows("next(-1)", prng, -1);
        assertThrows("next(33)", prng, 33);
        assertThrows("next(Integer.MIN_VALUE)", prng, Integer.MIN_VALUE);
    }

    private static void assertThrows(String what, PseudoRandom prng, int bits) {
        try {
            prng.next(bits);
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
    }

    /** {@code getSeed()} handed out the array the generator keeps. */
    @Test
    public void testGetSeedIsADefensiveCopy() {
        PseudoRandom prng = new Stc64(42L);
        long[] first = prng.getSeed();
        long original = first[0];
        first[0] = 0xDEADBEEFL;
        long[] second = prng.getSeed();
        assertEquals("getSeed() handed out its own array", original, second[0]);
        assertNotEquals("the two calls returned the same array", first, second);
    }

    /**
     * The range checks used {@code min > max}, which is false for a {@code NaN}
     * bound, so the degenerate range went through and the draw came back
     * {@code NaN}. They read {@code !(min <= max)} now. A non-positive standard
     * deviation was unchecked as well.
     */
    @Test
    public void testDegenerateArgumentsAreRejected() {
        final PseudoRandom prng = new Stc64(42L);
        expectIae("nextDouble(NaN, 1)", new Block() {
            @Override
            public void run() {
                prng.nextDouble(Double.NaN, 1.0);
            }
        });
        expectIae("nextDouble(0, NaN)", new Block() {
            @Override
            public void run() {
                prng.nextDouble(0.0, Double.NaN);
            }
        });
        expectIae("nextFloat(NaN, 1)", new Block() {
            @Override
            public void run() {
                prng.nextFloat(Float.NaN, 1.0f);
            }
        });
        expectIae("doubles(NaN, 1)", new Block() {
            @Override
            public void run() {
                prng.doubles(Double.NaN, 1.0);
            }
        });
        expectIae("nextGaussian(0, -1)", new Block() {
            @Override
            public void run() {
                prng.nextGaussian(0.0, -1.0);
            }
        });
        expectIae("nextGaussian(0, 0)", new Block() {
            @Override
            public void run() {
                prng.nextGaussian(0.0, 0.0);
            }
        });
        expectIae("nextDouble(1, 0)", new Block() {
            @Override
            public void run() {
                prng.nextDouble(1.0, 0.0);
            }
        });
        expectIae("nextLong(0)", new Block() {
            @Override
            public void run() {
                prng.nextLong(0L);
            }
        });
        expectIae("longs(-1)", new Block() {
            @Override
            public void run() {
                prng.longs(-1L);
            }
        });
        expectIae("ints(5, 3)", new Block() {
            @Override
            public void run() {
                prng.ints(5, 3);
            }
        });
    }

    private interface Block {
        void run();
    }

    private static void expectIae(String what, Block block) {
        try {
            block.run();
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
    }

    /**
     * {@code newInstance()} derived the second generator from
     * {@code getSeed()[0]}, the seed the source started from rather than the
     * state it has reached, so every call on the same source returned the same
     * stream. {@code BetaSpliterator} and {@code FisherFSpliterator} fall back
     * to it whenever the source cannot split, which correlated repeated calls.
     */
    @Test
    public void testNewInstanceUsesTheCurrentState() {
        AbstractRng64 source = new MersenneTwister64(42L);
        long first = source.newInstance().nextLong();
        long second = source.newInstance().nextLong();
        long third = source.newInstance().nextLong();
        assertNotEquals("two derived generators started identically", first, second);
        assertNotEquals("two derived generators started identically", second, third);
    }

    /** A gaussian mean and standard deviation the polar method has to hit. */
    @Test
    public void testGaussiansHaveTheRequestedMoments() {
        PseudoRandom prng = new Stc64(20260823L);
        int n = 400000;
        double sum = 0.0;
        double sumSq = 0.0;
        for (int i = 0; i < n; ++i) {
            double v = prng.nextGaussian(3.0, 2.0);
            sum += v;
            sumSq += v * v;
        }
        double mean = sum / n;
        double sd = Math.sqrt(sumSq / n - mean * mean);
        assertEquals("mean", 3.0, mean, 0.02);
        assertEquals("standard deviation", 2.0, sd, 0.02);
    }
}
