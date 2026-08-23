package math.list;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.List;
import java.util.function.DoublePredicate;

import org.junit.Test;

/**
 * Tests {@link DoubleList#removeIf(DoublePredicate)}, the in-place counterpart
 * of {@link DoubleList#filter(DoublePredicate)}. It is the only way to delete
 * elements that {@code remove(double)} cannot reach, since that method compares
 * values and {@code NaN != NaN}.
 */
public class DoubleArrayListRemoveIfTest {

    /** The values on which the two comparison rules disagree, plus ordinary ones. */
    private static final double[] AWKWARD = { 0.0, -0.0, Double.NaN, 1.0, -1.0, 2.0, 3.0, Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY, Double.MIN_VALUE };

    private static final DoublePredicate IS_NAN = new DoublePredicate() {
        @Override
        public boolean test(double value) {
            return Double.isNaN(value);
        }
    };

    private static final DoublePredicate NEGATIVELY_SIGNED = new DoublePredicate() {
        @Override
        public boolean test(double value) {
            return Double.doubleToLongBits(value) < 0L;
        }
    };

    private static final DoublePredicate EVERYTHING = new DoublePredicate() {
        @Override
        public boolean test(double value) {
            return true;
        }
    };

    private static final DoublePredicate NOTHING = new DoublePredicate() {
        @Override
        public boolean test(double value) {
            return false;
        }
    };

    private static final DoublePredicate ABOVE_ONE = new DoublePredicate() {
        @Override
        public boolean test(double value) {
            return value > 1.0;
        }
    };

    private static final DoublePredicate[] PREDICATES = { IS_NAN, NEGATIVELY_SIGNED, EVERYTHING, NOTHING, ABOVE_ONE };

    private long lcg = 20260824L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return (lcg >>> 11) * 0x1.0p-53;
    }

    @Test
    public void everyNaNIsRemoved() {
        for (int trial = 0; trial < 500; ++trial) {
            double[] values = awkwardValues(1 + (int) (next() * 20));
            DoubleArrayList list = DoubleArrayList.of(values);
            int expected = 0;
            for (int i = 0; i < values.length; ++i) {
                if (!Double.isNaN(values[i])) {
                    ++expected;
                }
            }
            boolean changed = list.removeIf(IS_NAN);
            assertEquals(Arrays.toString(values), expected != values.length, changed);
            assertEquals(Arrays.toString(values), expected, list.size());
            for (int i = 0; i < list.size(); ++i) {
                assertFalse(Arrays.toString(values), Double.isNaN(list.get(i)));
            }
        }
    }

    @Test
    public void removeIfReachesWhatRemoveCannot() {
        DoubleArrayList list = DoubleArrayList.of(1.0, Double.NaN, 2.0, Double.NaN);
        // remove(double) compares values, so it cannot touch a NaN at all
        assertFalse(list.remove(Double.NaN));
        assertEquals(4, list.size());
        assertTrue(list.removeIf(IS_NAN));
        assertArrayEquals(new double[] { 1.0, 2.0 }, list.toArray(), 0.0);

        // and remove(double) cannot tell the two zeroes apart, a predicate can
        DoubleArrayList zeroes = DoubleArrayList.of(0.0, -0.0, 0.0, -0.0);
        assertTrue(zeroes.removeIf(NEGATIVELY_SIGNED));
        assertEquals(2, zeroes.size());
        assertEquals(0L, Double.doubleToLongBits(zeroes.get(0)));
        assertEquals(0L, Double.doubleToLongBits(zeroes.get(1)));
    }

    @Test
    public void removeIfIsTheOppositeOfFilter() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] values = awkwardValues((int) (next() * 20));
            DoublePredicate predicate = PREDICATES[(int) (next() * PREDICATES.length)];
            DoubleArrayList list = DoubleArrayList.of(values);
            DoubleList kept = list.filter(predicate.negate());
            boolean changed = list.removeIf(predicate);
            String message = Arrays.toString(values);
            assertEquals(message, kept.size() != values.length, changed);
            assertTrue(message, kept.equals(list));
        }
    }

    @Test
    public void removeIfAgreesWithTheListOfBoxedDoubles() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] values = awkwardValues((int) (next() * 20));
            final DoublePredicate predicate = PREDICATES[(int) (next() * PREDICATES.length)];
            List<Double> boxed = boxed(values);
            boolean boxedChanged = boxed.removeIf(d -> predicate.test(d.doubleValue()));
            DoubleArrayList list = DoubleArrayList.of(values);
            boolean changed = list.removeIf(predicate);
            String message = Arrays.toString(values);
            assertEquals(message, boxedChanged, changed);
            // Arrays.equals compares the bits, so this does not let -0.0 pass for 0.0
            assertTrue(message, Arrays.equals(unboxed(boxed), list.toArray()));
        }
    }

    @Test
    public void theFastPathAgreesWithTheIteratorRoute() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] values = awkwardValues((int) (next() * 20));
            DoublePredicate predicate = PREDICATES[(int) (next() * PREDICATES.length)];
            DoubleArrayList naive = DoubleArrayList.of(values.clone());
            DoubleArrayList list = DoubleArrayList.of(values);
            String message = Arrays.toString(values);
            assertEquals(message, naiveRemoveIf(naive, predicate), list.removeIf(predicate));
            assertTrue(message, naive.equals(list));
        }
    }

    @Test
    public void nothingRemovedLeavesTheArrayAndTheModCountAlone() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0);
        double[] before = list.getArrayUnsafe().clone();
        // an iterator taken before the call must survive it
        DListIterator it = list.listIterator();
        assertFalse(list.removeIf(NOTHING));
        assertTrue(Arrays.equals(before, list.getArrayUnsafe()));
        int seen = 0;
        while (it.hasNext()) {
            it.next();
            ++seen;
        }
        assertEquals(4, seen);
    }

    @Test
    public void theEdgesAreCovered() {
        assertFalse(DoubleArrayList.of().removeIf(EVERYTHING));
        assertFalse(DoubleArrayList.of().removeIf(NOTHING));

        DoubleArrayList single = DoubleArrayList.of(7.0);
        assertTrue(single.removeIf(EVERYTHING));
        assertEquals(0, single.size());

        DoubleArrayList all = DoubleArrayList.of(1.0, 2.0, 3.0);
        assertTrue(all.removeIf(EVERYTHING));
        assertEquals(0, all.size());

        DoubleArrayList none = DoubleArrayList.of(1.0, 2.0, 3.0);
        assertFalse(none.removeIf(NOTHING));
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0 }, none.toArray(), 0.0);

        DoubleArrayList first = DoubleArrayList.of(1.0, 2.0, 3.0);
        assertTrue(first.removeIf(v -> v == 1.0));
        assertArrayEquals(new double[] { 2.0, 3.0 }, first.toArray(), 0.0);

        DoubleArrayList last = DoubleArrayList.of(1.0, 2.0, 3.0);
        assertTrue(last.removeIf(v -> v == 3.0));
        assertArrayEquals(new double[] { 1.0, 2.0 }, last.toArray(), 0.0);
    }

    @Test
    public void aSubListRemovesFromTheRoot() {
        DoubleArrayList root = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
        DoubleList sub = root.subList(2, 7);
        assertTrue(sub.removeIf(v -> v % 2.0 == 0.0));
        assertArrayEquals(new double[] { 3.0, 5.0, 7.0 }, sub.toArray(), 0.0);
        // the elements outside the range keep their places
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0, 5.0, 7.0, 8.0, 9.0, 10.0 }, root.toArray(), 0.0);
        assertEquals(8, root.size());
        assertEquals(3, sub.size());
        // and the sub-list is still usable afterwards
        assertFalse(sub.removeIf(IS_NAN));
        assertEquals(3, sub.size());
    }

    @Test
    public void aNestedSubListKeepsAllSizesInStep() {
        DoubleArrayList root = new DoubleArrayList(new double[] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 });
        DoubleList outer = root.subList(1, 9);
        DoubleList inner = outer.subList(1, 6);
        assertArrayEquals(new double[] { 2.0, 3.0, 4.0, 5.0, 6.0 }, inner.toArray(), 0.0);
        assertTrue(inner.removeIf(v -> v == 3.0 || v == 5.0));
        assertEquals(3, inner.size());
        assertEquals(6, outer.size());
        assertEquals(8, root.size());
        assertArrayEquals(new double[] { 0.0, 1.0, 2.0, 4.0, 6.0, 7.0, 8.0, 9.0 }, root.toArray(), 0.0);
        assertArrayEquals(new double[] { 1.0, 2.0, 4.0, 6.0, 7.0, 8.0 }, outer.toArray(), 0.0);
    }

    @Test
    public void aThrowingPredicateLeavesTheSurvivorsAndTheUnexamined() {
        DoubleArrayList list = new DoubleArrayList(30, 0.0);
        for (int i = 0; i < 30; ++i) {
            list.set(i, i);
        }
        try {
            list.removeIf(v -> {
                if (v == 10.0) {
                    throw new IllegalStateException("boom");
                }
                return v % 2.0 == 0.0;
            });
            fail("the predicate should have thrown");
        } catch (IllegalStateException expected) {
            assertEquals("boom", expected.getMessage());
        }
        // the elements the predicate never saw count as survivors, so the odd
        // ones below ten are gone and everything from ten upwards is untouched
        double[] expected = new double[25];
        expected[0] = 1.0;
        expected[1] = 3.0;
        expected[2] = 5.0;
        expected[3] = 7.0;
        expected[4] = 9.0;
        for (int i = 10; i < 30; ++i) {
            expected[i - 5] = i;
        }
        assertArrayEquals(expected, list.toArray(), 0.0);
        assertEquals(25, list.size());
    }

    @Test
    public void aPredicateThatChangesTheListIsDetected() {
        final DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        try {
            list.removeIf(v -> {
                if (v == 3.0) {
                    list.add(99.0);
                }
                return v == 1.0;
            });
            fail("a structural modification should have been detected");
        } catch (ConcurrentModificationException expected) {
            // the list is left in whatever state the predicate put it in,
            // which is why the javadoc calls this undefined
        }
    }

    @Test
    public void aNullPredicateIsRejected() {
        try {
            DoubleArrayList.of(1.0, 2.0).removeIf(null);
            fail("null should have been rejected");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            DoubleArrayList.of().removeIf(null);
            fail("null should have been rejected on an empty list too");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            DoubleArrayList.of(1.0, 2.0).subList(0, 1).removeIf(null);
            fail("null should have been rejected through a sub-list");
        } catch (NullPointerException expected) {
            // expected
        }
    }

    // ----- helpers ------------------------------------------------------

    /**
     * The route the default implementation of
     * {@link DoubleList#removeIf(DoublePredicate)} takes. It is spelled out
     * here because that default cannot be called from a subclass of
     * {@link DoubleArrayList}, which overrides it.
     */
    private static boolean naiveRemoveIf(DoubleList list, DoublePredicate filter) {
        boolean removed = false;
        DListIterator it = list.listIterator();
        while (it.hasNext()) {
            if (filter.test(it.next())) {
                it.remove();
                removed = true;
            }
        }
        return removed;
    }

    private double[] awkwardValues(int length) {
        double[] values = new double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = AWKWARD[(int) (next() * AWKWARD.length)];
        }
        return values;
    }

    private static List<Double> boxed(double[] values) {
        List<Double> list = new ArrayList<Double>(values.length);
        for (int i = 0; i < values.length; ++i) {
            list.add(Double.valueOf(values[i]));
        }
        return list;
    }

    private static double[] unboxed(List<Double> values) {
        double[] array = new double[values.size()];
        for (int i = 0; i < array.length; ++i) {
            array[i] = values.get(i).doubleValue();
        }
        return array;
    }
}
