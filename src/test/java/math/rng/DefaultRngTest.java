package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

/**
 * {@link DefaultRng} is what a caller with no reason to pick a particular
 * algorithm gets, and what the rest of the library funnels through, but it had
 * no test of its own. These are the properties a caller is entitled to whatever
 * generator sits behind the factory: the same seed gives the same stream, a
 * seed of nothing is refused, generators handed out as independent are
 * independent, and the default can split.
 */
public class DefaultRngTest {

    private static final long SEED = 20260825L;

    private static final int COUNT = 8;

    private static final int PER_GENERATOR = 50000;

    /** What a caller who stored a seed relies on. */
    @Test
    public void theDefaultIsReproducibleFromItsSeed() {
        long[] first = new long[2000];
        long[] second = new long[2000];
        DefaultRng.newPseudoRandom(SEED).nextLongs(first);
        DefaultRng.newPseudoRandom(SEED).nextLongs(second);
        assertTrue("the default is not reproducible from its seed", Arrays.equals(first, second));

        long[] other = new long[2000];
        DefaultRng.newPseudoRandom(SEED + 1L).nextLongs(other);
        assertFalse("a different seed gave the same stream", Arrays.equals(first, other));

        long[] a = new long[2000];
        long[] b = new long[2000];
        DefaultRng.newPseudoRandom(new long[] { SEED, 7L }).nextLongs(a);
        DefaultRng.newPseudoRandom(new long[] { SEED, 7L }).nextLongs(b);
        assertTrue("the array form is not reproducible", Arrays.equals(a, b));

        // the array is hashed down to a single long, and that hash has to
        // depend on every element it is given
        DefaultRng.newPseudoRandom(new long[] { SEED, 8L }).nextLongs(b);
        assertFalse("the last element of the seed array is ignored", Arrays.equals(a, b));
    }

    /**
     * The generator behind the factory accepts {@code null} and an empty array,
     * mapping both onto one fixed seed. The factory does not, and that is a
     * guard of its own rather than something it inherits: without it, seeding
     * from nothing stops being an error and silently becomes that one
     * well-known stream.
     */
    @Test
    public void theSeedArrayIsValidated() {
        try {
            DefaultRng.newPseudoRandom((long[]) null);
            fail("a null seed array was accepted");
        } catch (NullPointerException expected) {
            // as specified
        }
        try {
            DefaultRng.newPseudoRandom(new long[0]);
            fail("an empty seed array was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
        assertTrue("a one-element array is enough", DefaultRng.newPseudoRandom(new long[] { 1L }).nextLong() != 0L);
    }

    /**
     * Two generators built without a seed have to be two generators. This one
     * is deliberately not reproducible -- the unseeded constructor is the case
     * under test -- but the assertion is not delicate: 40000 draws out of
     * 2<sup>64</sup> collide with probability about 4e-11 when the streams are
     * independent.
     */
    @Test
    public void unseededGeneratorsDoNotShareAStream() {
        long[] a = new long[20000];
        long[] b = new long[20000];
        DefaultRng.newPseudoRandom().nextLongs(a);
        DefaultRng.newPseudoRandom().nextLongs(b);
        assertEquals("two unseeded generators produced a common value", 0, common(a, b));
    }

    /**
     * {@code newIndepPseudoRandoms} promises sequences that do not overlap in
     * any practically reachable way. 400000 draws off eight generators, all
     * distinct; the probability of a collision between independent streams is
     * about 4e-9.
     */
    @Test
    public void independentGeneratorsDoNotShareAStream() {
        PseudoRandom[] prng = DefaultRng.newIndepPseudoRandoms(COUNT);
        assertEquals("wrong number of generators", COUNT, prng.length);

        Set<Long> seen = new HashSet<Long>(2 * COUNT * PER_GENERATOR);
        for (int i = 0; i < prng.length; ++i) {
            long[] v = new long[PER_GENERATOR];
            prng[i].nextLongs(v);
            for (int k = 0; k < v.length; ++k) {
                if (!seen.add(Long.valueOf(v[k]))) {
                    fail("generator " + i + " repeated a value that had already been drawn");
                }
            }
        }

        assertEquals("a single generator is a special case", 1, DefaultRng.newIndepPseudoRandoms(1).length);
        expectIae("count 0", 0);
        expectIae("a negative count", -1);
    }

    /**
     * A derived generator is derived from the state its source has reached, not
     * from the seed the source started at, so two calls on the same source give
     * two different generators -- and neither of them replays what the source
     * itself goes on to produce.
     */
    @Test
    public void theDerivedGeneratorIsIndependent() {
        PseudoRandom source = DefaultRng.newPseudoRandom(SEED);
        PseudoRandom first = DefaultRng.newIndepPseudoRandom(source);
        PseudoRandom second = DefaultRng.newIndepPseudoRandom(source);
        assertNotSame("the source itself was handed back", source, first);
        assertNotSame("the source itself was handed back", source, second);

        long[] a = new long[2000];
        long[] b = new long[2000];
        long[] rest = new long[2000];
        first.nextLongs(a);
        second.nextLongs(b);
        source.nextLongs(rest);
        assertFalse("two derived generators started identically", Arrays.equals(a, b));
        assertEquals("two derived generators overlap", 0, common(a, b));
        assertEquals("a derived generator replays its source", 0, common(a, rest));
        assertEquals("a derived generator replays its source", 0, common(b, rest));

        // the javadoc promises exactly one value, which is what makes the
        // method usable in the middle of a reproducible run: what the source
        // produces next has to be what it would have produced anyway
        PseudoRandom probe = DefaultRng.newPseudoRandom(SEED);
        PseudoRandom check = DefaultRng.newPseudoRandom(SEED);
        DefaultRng.newIndepPseudoRandom(probe);
        check.nextLong();
        assertEquals("the source was not advanced by exactly one value", check.nextLong(), probe.nextLong());
    }

    /**
     * The default has to be splittable. {@code BetaSpliterator} and
     * {@code FisherFSpliterator} draw from two generators and derive the second
     * one through {@code split()}; a default that cannot split sends every beta
     * and every F stream through the reflective fallback in
     * {@link AbstractRng64#newInstance()} instead.
     */
    @Test
    public void theDefaultCanSplit() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        assertTrue("the default generator " + prng.getAlgorithm() + " cannot split",
                prng instanceof SplittablePseudoRandom);

        PseudoRandom child = ((SplittablePseudoRandom) prng).split();
        long[] a = new long[2000];
        long[] b = new long[2000];
        child.nextLongs(a);
        prng.nextLongs(b);
        assertEquals("a split child replays its parent", 0, common(a, b));

        double[] v = DefaultRng.newPseudoRandom(SEED).beta(2000L, 2.0, 3.0).toArray();
        assertEquals("the beta stream is short", 2000, v.length);
        for (int i = 0; i < v.length; ++i) {
            assertTrue("beta variate out of range: " + v[i], v[i] > 0.0 && v[i] < 1.0);
        }
    }

    /** The number of values the two arrays have in common. */
    private static int common(long[] a, long[] b) {
        Set<Long> seen = new HashSet<Long>(2 * a.length);
        for (int i = 0; i < a.length; ++i) {
            seen.add(Long.valueOf(a[i]));
        }
        int count = 0;
        for (int i = 0; i < b.length; ++i) {
            if (seen.contains(Long.valueOf(b[i]))) {
                ++count;
            }
        }
        return count;
    }

    private static void expectIae(String what, int count) {
        try {
            DefaultRng.newIndepPseudoRandoms(count);
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
    }
}
