package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.lang.reflect.Field;
import java.util.Arrays;

import org.junit.Test;

/**
 * How the generators take their seed, and the two seeds that used to break
 * them.
 */
public class GeneratorSeedingTest {

    /** the eleven concrete generators, by their {@code (long)} constructor */
    private static final Class<?>[] GENERATORS = { Interleaved4Stc64.class, Lcg64Xor1024Mix.class,
            MersenneTwister64.class, Sfc64.class, SplitMix64.class, Stc64.class, XorShift1024Star.class,
            XorShift128Plus.class, XorShift64Star.class, XorShiftRot256PlusPlus.class,
            XorShiftRot256StarStar.class };

    private static PseudoRandom seeded(Class<?> c, long seed) throws Exception {
        return (PseudoRandom) c.getDeclaredConstructor(long.class).newInstance(Long.valueOf(seed));
    }

    /**
     * The all-zero state is the fixed point of every xorshift: once there,
     * {@code nextLong()} returns zero forever and {@code nextGaussian()} never
     * terminates, because the polar method rejects every candidate. The
     * no-argument constructor looped until it had a non-zero state, the seeded
     * one did not -- and its seed goes through
     * {@code BitMix.rrxmrrxmsx(seed + GOLDEN)}, which is zero for exactly one
     * argument: {@code -GOLDEN}.
     */
    @Test
    public void testTheXorshiftFamilySurvivesTheDeadlySeed() throws Exception {
        long deadly = 7046029254386353131L; // -GOLDEN
        Class<?>[] marsaglia = { XorShift64Star.class, XorShift128Plus.class, XorShift1024Star.class,
                XorShiftRot256PlusPlus.class, XorShiftRot256StarStar.class, Lcg64Xor1024Mix.class };
        for (int i = 0; i < marsaglia.length; ++i) {
            PseudoRandom prng = seeded(marsaglia[i], deadly);
            long zeros = 0L;
            for (int k = 0; k < 1000; ++k) {
                if (prng.nextLong() == 0L) {
                    ++zeros;
                }
            }
            assertTrue(marsaglia[i].getSimpleName() + " is stuck at the zero state", zeros < 5L);
            // this is the call that used to hang forever
            double g = prng.nextGaussian();
            assertTrue(marsaglia[i].getSimpleName() + " gaussian " + g, !Double.isNaN(g));
        }
        // and the array form of the same seed
        PseudoRandom prng = new XorShift64Star(new long[] { deadly });
        assertTrue("the array constructor is stuck at zero", prng.nextLong() != 0L);
    }

    /**
     * {@code Lcg64Xor1024Mix} draws the addend of its linear congruential part
     * from the seeder, and an LCG modulo a power of two only reaches full
     * period when the addend is odd. It was taken as drawn, so about half of
     * all seeds got a degenerate sub-generator -- 1012 of the first 2000.
     * <p>
     * The addend is not observable through the public API, because the output
     * is mixed with a 1024-bit xorshift that keeps the period long either way.
     * The field is what the repair changed, so the field is what is asserted.
     */
    @Test
    public void testTheLcgAddendIsAlwaysOdd() throws Exception {
        Field a = Lcg64Xor1024Mix.class.getDeclaredField("a");
        a.setAccessible(true);
        for (long seed = 1L; seed <= 2000L; ++seed) {
            long addend = a.getLong(new Lcg64Xor1024Mix(seed));
            assertTrue("seed " + seed + " gave the even addend " + addend, (addend & 1L) != 0L);
        }
        for (int i = 0; i < 200; ++i) {
            long addend = a.getLong(new Lcg64Xor1024Mix());
            assertTrue("an unseeded instance gave the even addend " + addend, (addend & 1L) != 0L);
        }
    }

    /**
     * {@code MersenneTwister64} consumes the array as seed material itself, so
     * an empty one carries nothing to seed with and used to fail inside
     * {@code init_by_array64} with an {@code ArrayIndexOutOfBoundsException}.
     * It is the one generator of the package that rejects both {@code null} and
     * an empty array; the other ten map them onto a fixed seed through
     * {@link SplitMix64Seed}, which is the documented convention.
     */
    @Test
    public void testMersenneTwisterValidatesItsSeedArray() {
        try {
            new MersenneTwister64((long[]) null);
            fail("a null seed array was accepted");
        } catch (NullPointerException expected) {
            // as specified
        }
        try {
            new MersenneTwister64(new long[0]);
            fail("an empty seed array was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
        assertTrue("a one-element array is enough", new MersenneTwister64(new long[] { 1L }).nextLong() != 0L);
    }

    /** The same seed has to give the same stream, for every generator. */
    @Test
    public void testSeededGeneratorsAreReproducible() throws Exception {
        for (int i = 0; i < GENERATORS.length; ++i) {
            long[] first = new long[2000];
            long[] second = new long[2000];
            seeded(GENERATORS[i], 20260823L).nextLongs(first);
            seeded(GENERATORS[i], 20260823L).nextLongs(second);
            assertTrue(GENERATORS[i].getSimpleName() + " is not reproducible from its seed",
                    Arrays.equals(first, second));
            long[] other = new long[2000];
            seeded(GENERATORS[i], 20260824L).nextLongs(other);
            assertTrue(GENERATORS[i].getSimpleName() + " gave the same stream for a different seed",
                    !Arrays.equals(first, other));
        }
    }

    /**
     * A seeded generator has to produce a stream, not a constant: no generator
     * may collapse onto a short cycle for an ordinary seed.
     */
    @Test
    public void testNoSeedCollapsesTheStream() throws Exception {
        for (int i = 0; i < GENERATORS.length; ++i) {
            for (long seed = 0L; seed < 40L; ++seed) {
                long[] v = new long[2000];
                seeded(GENERATORS[i], seed).nextLongs(v);
                long[] sorted = v.clone();
                Arrays.sort(sorted);
                int distinct = 1;
                for (int k = 1; k < sorted.length; ++k) {
                    if (sorted[k] != sorted[k - 1]) {
                        ++distinct;
                    }
                }
                assertEquals(GENERATORS[i].getSimpleName() + " at seed " + seed + " repeated itself", 2000,
                        distinct);
            }
        }
    }

    /**
     * {@code getSeed()} has to report something for a seeded generator, and the
     * report has to survive being handed out.
     */
    @Test
    public void testEveryGeneratorReportsASeed() throws Exception {
        for (int i = 0; i < GENERATORS.length; ++i) {
            PseudoRandom prng = seeded(GENERATORS[i], 20260823L);
            long[] seed = prng.getSeed();
            assertTrue(GENERATORS[i].getSimpleName() + " reports no seed", seed != null && seed.length > 0);
            assertTrue(GENERATORS[i].getSimpleName() + " reports an all-zero seed", notAllZero(seed));
            assertEquals(GENERATORS[i].getSimpleName() + " has no algorithm name",
                    GENERATORS[i].getSimpleName(), prng.getAlgorithm());
        }
    }

    private static boolean notAllZero(long[] seed) {
        for (int i = 0; i < seed.length; ++i) {
            if (seed[i] != 0L) {
                return true;
            }
        }
        return false;
    }
}
