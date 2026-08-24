package math.list;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

/**
 * Tests that {@code equals} and {@code hashCode} of {@link DoubleArrayList}
 * agree. They used to disagree on signed zero -- {@code {0.0}} and
 * {@code {-0.0}} were equal but hashed differently, which is the direction that
 * breaks a {@code HashMap} -- and a list holding a {@code NaN} was not even
 * equal to its own copy.
 */
public class DoubleArrayListEqualityTest {

    /** The values that make the two comparisons disagree, plus ordinary ones. */
    private static final double[] AWKWARD = { 0.0, -0.0, Double.NaN, 1.0, -1.0, 2.0, Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY, Double.MIN_VALUE, -Double.MIN_VALUE };

    private long lcg = 20260824L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return (lcg >>> 11) * 0x1.0p-53;
    }

    @Test
    public void signedZeroesAreNotEqualAndDoNotClaimToBe() {
        DoubleArrayList plus = DoubleArrayList.of(0.0);
        DoubleArrayList minus = DoubleArrayList.of(-0.0);
        assertFalse(plus.equals(minus));
        assertFalse(minus.equals(plus));
        assertNotEquals(plus.hashCode(), minus.hashCode());
    }

    @Test
    public void aNaNEqualsItself() {
        DoubleArrayList one = DoubleArrayList.of(1.0, Double.NaN, 3.0);
        DoubleArrayList other = DoubleArrayList.of(1.0, Double.NaN, 3.0);
        assertTrue(one.equals(other));
        assertEquals(one.hashCode(), other.hashCode());
    }

    @Test
    public void aListEqualsItsOwnCopy() {
        for (int trial = 0; trial < 500; ++trial) {
            double[] values = awkwardValues(1 + (int) (next() * 12));
            DoubleArrayList list = DoubleArrayList.of(values);
            assertTrue(Arrays.toString(values), list.equals(new DoubleArrayList(list)));
            assertTrue(Arrays.toString(values), list.equals(DoubleArrayList.of(values.clone())));
            assertTrue(Arrays.toString(values), list.equals(list));
        }
    }

    @Test
    public void equalListsHashAlike() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] first = awkwardValues(1 + (int) (next() * 6));
            double[] second = nearlyTheSame(first);
            DoubleArrayList a = DoubleArrayList.of(first);
            DoubleArrayList b = DoubleArrayList.of(second);
            if (a.equals(b)) {
                assertEquals(a + " vs " + b, a.hashCode(), b.hashCode());
            }
            // symmetry, which the two comparison rules could also break
            assertEquals(a + " vs " + b, a.equals(b), b.equals(a));
        }
    }

    @Test
    public void equalsAgreesWithArraysEquals() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] first = awkwardValues(1 + (int) (next() * 8));
            double[] second = nearlyTheSame(first);
            // Arrays.equals compares double values by their bits as well, and
            // is the independent statement of the same rule
            assertEquals(Arrays.toString(first) + " vs " + Arrays.toString(second), Arrays.equals(first, second),
                    DoubleArrayList.of(first).equals(DoubleArrayList.of(second)));
        }
    }

    @Test
    public void equalsAgreesWithAListOfBoxedDoubles() {
        for (int trial = 0; trial < 2000; ++trial) {
            double[] first = awkwardValues(1 + (int) (next() * 8));
            double[] second = nearlyTheSame(first);
            String message = Arrays.toString(first) + " vs " + Arrays.toString(second);
            assertEquals(message, boxed(first).equals(boxed(second)),
                    DoubleArrayList.of(first).equals(DoubleArrayList.of(second)));
        }
    }

    @Test
    public void aHashSetFindsAnEqualList() {
        Set<DoubleList> set = new HashSet<DoubleList>();
        for (int trial = 0; trial < 300; ++trial) {
            double[] values = awkwardValues(1 + (int) (next() * 8));
            set.add(DoubleArrayList.of(values));
            assertTrue(Arrays.toString(values), set.contains(DoubleArrayList.of(values.clone())));
        }
    }

    @Test
    public void aHashMapKeepsSignedZeroesApart() {
        Map<DoubleList, String> map = new HashMap<DoubleList, String>();
        map.put(DoubleArrayList.of(0.0), "zero");
        map.put(DoubleArrayList.of(-0.0), "minus zero");
        // two keys that are not equal, so two entries, and each one findable
        assertEquals(2, map.size());
        assertEquals("zero", map.get(DoubleArrayList.of(0.0)));
        assertEquals("minus zero", map.get(DoubleArrayList.of(-0.0)));
        map.put(DoubleArrayList.of(Double.NaN), "not a number");
        assertEquals("not a number", map.get(DoubleArrayList.of(Double.NaN)));
    }

    @Test
    public void aSubListEqualsARootListHoldingTheSameValues() {
        for (int trial = 0; trial < 300; ++trial) {
            int n = 2 + (int) (next() * 10);
            double[] values = awkwardValues(n);
            DoubleArrayList root = DoubleArrayList.of(values);
            for (int from = 0; from < n; ++from) {
                for (int to = from + 1; to <= n; ++to) {
                    DoubleList sub = root.subList(from, to);
                    DoubleArrayList same = DoubleArrayList.of(Arrays.copyOfRange(values, from, to));
                    assertTrue(from + ".." + to, sub.equals(same));
                    assertTrue(from + ".." + to, same.equals(sub));
                    assertEquals(from + ".." + to, same.hashCode(), sub.hashCode());
                }
            }
        }
    }

    @Test
    public void listsOfDifferentSizeAreNotEqual() {
        assertFalse(DoubleArrayList.of(1.0, 2.0).equals(DoubleArrayList.of(1.0)));
        assertFalse(DoubleArrayList.of(1.0).equals(DoubleArrayList.of(1.0, 2.0)));
        assertFalse(DoubleArrayList.of().equals(DoubleArrayList.of(0.0)));
        assertTrue(DoubleArrayList.of().equals(DoubleArrayList.of()));
        assertFalse(DoubleArrayList.of(1.0).equals("not a list"));
        assertFalse(DoubleArrayList.of(1.0).equals(null));
    }

    @Test
    public void theSearchMethodsStillCompareValues() {
        // deliberately the other rule, and documented as o == e: a search over
        // a list of numbers asks for a number, not for a bit pattern
        DoubleArrayList list = DoubleArrayList.of(1.0, Double.NaN, 0.0, 3.0);
        assertTrue(list.contains(-0.0));
        assertEquals(2, list.indexOf(-0.0));
        assertEquals(2, list.lastIndexOf(-0.0));
        assertFalse(list.contains(Double.NaN));
        assertEquals(-1, list.indexOf(Double.NaN));
        assertFalse(list.remove(Double.NaN));
        assertEquals(4, list.size());
    }

    // ----- helpers ------------------------------------------------------

    /** A short list drawn from the values that make the two rules differ. */
    private double[] awkwardValues(int length) {
        double[] values = new double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = AWKWARD[(int) (next() * AWKWARD.length)];
        }
        return values;
    }

    /**
     * A copy that is often identical and otherwise differs in a place where
     * the two comparison rules disagree. Drawing both lists independently
     * almost never produces such a pair, and the cross-checks then pass
     * against any implementation at all.
     */
    private double[] nearlyTheSame(double[] values) {
        double[] copy = values.clone();
        int mutations = (int) (next() * 3.0);
        for (int k = 0; k < mutations; ++k) {
            int i = (int) (next() * copy.length);
            if (copy[i] == 0.0) {
                copy[i] = -copy[i];
            } else {
                copy[i] = AWKWARD[(int) (next() * AWKWARD.length)];
            }
        }
        return copy;
    }

    private static List<Double> boxed(double[] values) {
        List<Double> list = new ArrayList<Double>(values.length);
        for (int i = 0; i < values.length; ++i) {
            list.add(Double.valueOf(values[i]));
        }
        return list;
    }
}
