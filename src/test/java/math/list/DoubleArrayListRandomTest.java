package math.list;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import math.rng.DefaultRng;
import math.rng.PseudoRandom;

/**
 * Tests the three sources of randomness in {@link DoubleArrayList}:
 * {@code randomUniform}, {@code randomNormal} and {@code shuffle}. They used to
 * draw from {@link java.util.concurrent.ThreadLocalRandom}, which no caller can
 * seed, and the uniform draw formed {@code max - min} first, so an interval as
 * wide as the {@code double} range produced a list of infinities.
 */
public class DoubleArrayListRandomTest {

    private static final long SEED = 20260824L;

    // ----- reproducibility ----------------------------------------------

    @Test
    public void theSameSeedGivesTheSameList() {
        double[] first = DoubleArrayList.randomUniform(-3.0, 7.0, 64, DefaultRng.newPseudoRandom(SEED)).toArray();
        double[] second = DoubleArrayList.randomUniform(-3.0, 7.0, 64, DefaultRng.newPseudoRandom(SEED)).toArray();
        assertTrue(Arrays.equals(first, second));

        first = DoubleArrayList.randomNormal(2.0, 0.5, 64, DefaultRng.newPseudoRandom(SEED)).toArray();
        second = DoubleArrayList.randomNormal(2.0, 0.5, 64, DefaultRng.newPseudoRandom(SEED)).toArray();
        assertTrue(Arrays.equals(first, second));
    }

    @Test
    public void theSameSeedGivesTheSamePermutation() {
        double[] first = ascending(50).shuffle(DefaultRng.newPseudoRandom(SEED)).toArray();
        double[] second = ascending(50).shuffle(DefaultRng.newPseudoRandom(SEED)).toArray();
        assertTrue(Arrays.equals(first, second));

        // and through a sub-list, which permutes a window of the root
        first = ascending(50).subList(10, 30).shuffle(DefaultRng.newPseudoRandom(SEED)).toArray();
        second = ascending(50).subList(10, 30).shuffle(DefaultRng.newPseudoRandom(SEED)).toArray();
        assertTrue(Arrays.equals(first, second));
    }

    @Test
    public void differentSeedsGiveDifferentLists() {
        double[] first = DoubleArrayList.randomUniform(0.0, 1.0, 32, DefaultRng.newPseudoRandom(SEED)).toArray();
        double[] second = DoubleArrayList.randomUniform(0.0, 1.0, 32, DefaultRng.newPseudoRandom(SEED + 1L)).toArray();
        assertFalse(Arrays.equals(first, second));
    }

    @Test
    public void theUnseededFactoriesStillDraw() {
        // no seed, no reproducibility, but 32 draws that agree twice running
        // would be a broken generator rather than a coincidence
        assertFalse(Arrays.equals(DoubleArrayList.randomUniform(0.0, 1.0, 32).toArray(),
                DoubleArrayList.randomUniform(0.0, 1.0, 32).toArray()));
        assertFalse(Arrays.equals(DoubleArrayList.randomNormal(0.0, 1.0, 32).toArray(),
                DoubleArrayList.randomNormal(0.0, 1.0, 32).toArray()));
        assertFalse(Arrays.equals(ascending(32).shuffle().toArray(), ascending(32).shuffle().toArray()));
    }

    // ----- the interval -------------------------------------------------

    @Test
    public void theValuesStayInsideTheInterval() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        double[][] intervals = { { 0.0, 1.0 }, { -3.0, 7.0 }, { -9.5, -2.25 }, { 1.0e300, 1.7e308 },
                { -1.0e-300, 1.0e-300 } };
        for (int k = 0; k < intervals.length; ++k) {
            double min = intervals[k][0];
            double max = intervals[k][1];
            double[] values = DoubleArrayList.randomUniform(min, max, 20000, prng).toArray();
            for (int i = 0; i < values.length; ++i) {
                String message = min + " .. " + max + " produced " + values[i];
                assertTrue(message, values[i] >= min);
                assertTrue(message, values[i] < max);
            }
        }
    }

    @Test
    public void anEmptyIntervalIsTheConstant() {
        double[] values = DoubleArrayList.randomUniform(2.5, 2.5, 16, DefaultRng.newPseudoRandom(SEED)).toArray();
        double[] expected = new double[16];
        Arrays.fill(expected, 2.5);
        assertTrue(Arrays.equals(expected, values));
    }

    @Test
    public void aWideIntervalDoesNotOverflow() {
        // min + (max - min) * u reads min + infinity * u here, which is what
        // the old formula did: every value came back infinite
        double[] values = DoubleArrayList
                .randomUniform(-Double.MAX_VALUE, Double.MAX_VALUE, 4096, DefaultRng.newPseudoRandom(SEED)).toArray();
        boolean negative = false;
        boolean positive = false;
        for (int i = 0; i < values.length; ++i) {
            assertFalse("value " + i + " is " + values[i], Double.isInfinite(values[i]));
            assertFalse("value " + i + " is NaN", Double.isNaN(values[i]));
            negative |= values[i] < 0.0;
            positive |= values[i] > 0.0;
        }
        // and it really is the whole range, not a collapsed half of it
        assertTrue(negative);
        assertTrue(positive);
    }

    @Test
    public void theMomentsAreRight() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        DoubleList uniform = DoubleArrayList.randomUniform(-2.0, 5.0, 200000, prng);
        // 1.5 and 7 / sqrt(12); the worst deviation measured over the seeds
        // tried is 4.3e-4 of the mean and 2.4e-3 of the standard deviation
        assertEquals(1.5, uniform.avg(), 0.02);
        assertEquals(7.0 / Math.sqrt(12.0), uniform.stddev(), 0.02);

        DoubleList normal = DoubleArrayList.randomNormal(3.0, 2.0, 200000, prng);
        assertEquals(3.0, normal.avg(), 0.02);
        assertEquals(2.0, normal.stddev(), 0.02);
    }

    // ----- what is refused ----------------------------------------------

    @Test
    public void aReversedIntervalIsRefused() {
        // it used to draw from (max, min] instead, which is not what the two
        // parameters are called
        expectIllegalArgument(1.0, 0.0);
        expectIllegalArgument(0.0, -Double.MIN_VALUE);
        expectIllegalArgument(Double.NaN, 1.0);
        expectIllegalArgument(0.0, Double.NaN);
        expectIllegalArgument(Double.NaN, Double.NaN);
    }

    @Test
    public void aNegativeSizeIsRefused() {
        try {
            DoubleArrayList.randomUniform(0.0, 1.0, -1, DefaultRng.newPseudoRandom(SEED));
            fail("no exception for a negative size");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("-1"));
        }
        try {
            DoubleArrayList.randomNormal(0.0, 1.0, -1, DefaultRng.newPseudoRandom(SEED));
            fail("no exception for a negative size");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("-1"));
        }
        // zero is not negative and stays a legal, empty list
        assertTrue(DoubleArrayList.randomUniform(0.0, 1.0, 0, DefaultRng.newPseudoRandom(SEED)).isEmpty());
        assertTrue(DoubleArrayList.randomNormal(0.0, 1.0, 0, DefaultRng.newPseudoRandom(SEED)).isEmpty());
    }

    @Test
    public void aStandardDeviationThatIsNotPositiveIsRefused() {
        double[] bad = { 0.0, -0.0, -1.0, Double.NaN, Double.NEGATIVE_INFINITY };
        for (int i = 0; i < bad.length; ++i) {
            try {
                DoubleArrayList.randomNormal(0.0, bad[i], 4, DefaultRng.newPseudoRandom(SEED));
                fail("no exception for a standard deviation of " + bad[i]);
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage(), expected.getMessage().contains("positive"));
            }
        }
    }

    @Test
    public void aNullGeneratorIsRefused() {
        try {
            DoubleArrayList.randomUniform(0.0, 1.0, 4, null);
            fail("no exception for a null generator");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            DoubleArrayList.randomNormal(0.0, 1.0, 4, null);
            fail("no exception for a null generator");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            ascending(4).shuffle(null);
            fail("no exception for a null generator");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            ascending(8).subList(2, 6).shuffle(null);
            fail("no exception for a null generator");
        } catch (NullPointerException expected) {
            // expected
        }
    }

    // ----- shuffle ------------------------------------------------------

    @Test
    public void theShuffleIsAPermutation() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        for (int size = 0; size < 40; ++size) {
            DoubleArrayList list = ascending(size);
            double[] before = list.toArray();
            list.shuffle(prng);
            double[] after = list.toArray();
            assertEquals(before.length, after.length);
            Arrays.sort(before);
            Arrays.sort(after);
            assertArrayEquals(before, after, 0.0);
        }
    }

    @Test
    public void aSubListShuffleLeavesTheRestAlone() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        for (int trial = 0; trial < 200; ++trial) {
            DoubleArrayList list = ascending(20);
            list.subList(5, 15).shuffle(prng);
            double[] values = list.toArray();
            for (int i = 0; i < 5; ++i) {
                assertEquals(i, values[i], 0.0);
            }
            for (int i = 15; i < 20; ++i) {
                assertEquals(i, values[i], 0.0);
            }
            double[] window = Arrays.copyOfRange(values, 5, 15);
            Arrays.sort(window);
            for (int i = 0; i < window.length; ++i) {
                assertEquals(5 + i, window[i], 0.0);
            }
        }
    }

    @Test
    public void everyPermutationOfFourElementsOccurs() {
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        Map<String, Integer> counts = new HashMap<String, Integer>();
        final int trials = 24000;
        for (int t = 0; t < trials; ++t) {
            DoubleArrayList list = ascending(4);
            list.shuffle(prng);
            String key = "" + (int) list.get(0) + (int) list.get(1) + (int) list.get(2) + (int) list.get(3);
            Integer seen = counts.get(key);
            counts.put(key, seen == null ? Integer.valueOf(1) : Integer.valueOf(seen.intValue() + 1));
        }
        assertEquals(24, counts.size());
        // the worst relative deviation from 1000 measured over the four seeds
        // tried is 0.083
        for (Map.Entry<String, Integer> e : counts.entrySet()) {
            assertEquals(e.getKey(), trials / 24.0, e.getValue().intValue(), 0.15 * trials / 24.0);
        }
    }

    @Test
    public void theArrayShuffleAgreesWithTheGetSetRoute() {
        // the route the default implementation of
        // DoubleList.shuffle(PseudoRandom) takes. It is spelled out here
        // because that default cannot be called from a subclass of
        // DoubleArrayList, which overrides it, and both draw the same values
        // in the same order, so they have to agree exactly
        for (int size = 0; size < 40; ++size) {
            DoubleArrayList fast = ascending(size);
            DoubleArrayList slow = ascending(size);
            fast.shuffle(DefaultRng.newPseudoRandom(SEED + size));
            getSetShuffle(slow, DefaultRng.newPseudoRandom(SEED + size));
            assertTrue("size " + size, Arrays.equals(fast.toArray(), slow.toArray()));
        }
        for (int size = 1; size < 20; ++size) {
            DoubleArrayList fast = ascending(2 * size);
            DoubleArrayList slow = ascending(2 * size);
            fast.subList(size / 2, size).shuffle(DefaultRng.newPseudoRandom(SEED + size));
            getSetShuffle(slow.subList(size / 2, size), DefaultRng.newPseudoRandom(SEED + size));
            assertTrue("sub-list " + size, Arrays.equals(fast.toArray(), slow.toArray()));
        }
    }

    @Test
    public void twoListsFromOneGeneratorDoNotRepeat() {
        // the generator carries its state across the calls, so a caller who
        // draws twice gets two different lists out of one seed
        PseudoRandom prng = DefaultRng.newPseudoRandom(SEED);
        Set<Double> seen = new HashSet<Double>();
        for (int k = 0; k < 8; ++k) {
            double[] values = DoubleArrayList.randomUniform(0.0, 1.0, 100, prng).toArray();
            for (int i = 0; i < values.length; ++i) {
                assertTrue("value " + values[i] + " came up twice", seen.add(Double.valueOf(values[i])));
            }
        }
    }

    // ----- helpers ------------------------------------------------------

    private static DoubleArrayList ascending(int size) {
        DoubleArrayList list = new DoubleArrayList(size);
        for (int i = 0; i < size; ++i) {
            list.add(i);
        }
        return list;
    }

    private static void getSetShuffle(DoubleList list, PseudoRandom prng) {
        for (int i = list.size(); i > 1; --i) {
            int j = prng.nextInt(i);
            double tmp = list.get(i - 1);
            list.set(i - 1, list.get(j));
            list.set(j, tmp);
        }
    }

    private static void expectIllegalArgument(double min, double max) {
        try {
            DoubleArrayList.randomUniform(min, max, 4, DefaultRng.newPseudoRandom(SEED));
            fail("no exception for the interval " + min + " .. " + max);
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
