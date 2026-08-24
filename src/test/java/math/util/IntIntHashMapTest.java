package math.util;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

/**
 * {@link IntIntHashMap} against {@link HashMap}, of which it is a primitive
 * translation: the same sequence of operations has to answer the same at every
 * step and leave both containers in the same state.
 */
public class IntIntHashMapTest {

    /** The value answered for an absent key, by contract not a value of the map. */
    private static final int NO_KEY = Integer.MIN_VALUE;

    /** Keys are drawn from a small pool so that hits and collisions are frequent. */
    private static final int[] KEYS = pool();

    private static final int STEPS = 20000;

    @Test
    public void randomOperationsAgreeWithTheJdkMap() {
        IntIntHashMap mine = new IntIntHashMap(NO_KEY);
        Map<Integer, Integer> ref = new HashMap<Integer, Integer>();
        long lcg = 192837465L;
        for (int step = 0; step < STEPS; ++step) {
            lcg = next(lcg);
            int op = (int) ((lcg >>> 32) % 8L);
            lcg = next(lcg);
            int k = KEYS[(int) ((lcg >>> 32) % KEYS.length)];
            lcg = next(lcg);
            int v = 1 + (int) ((lcg >>> 32) % 100L);
            Integer bk = Integer.valueOf(k);
            Integer bv = Integer.valueOf(v);
            String at = " (step " + step + ", key " + k + ")";
            switch (op) {
            case 0:
                assertEquals("putInt" + at, unbox(ref.put(bk, bv)), mine.putInt(k, v));
                break;
            case 1:
                assertEquals("put" + at, ref.put(bk, bv), mine.put(bk, bv));
                break;
            case 2:
                assertEquals("getInt" + at, unbox(ref.get(bk)), mine.getInt(k));
                break;
            case 3:
                assertEquals("get" + at, ref.get(bk), mine.get(bk));
                break;
            case 4:
                assertEquals("removeInt" + at, unbox(ref.remove(bk)), mine.removeInt(k));
                break;
            case 5:
                assertEquals("remove" + at, ref.remove(bk), mine.remove(bk));
                break;
            case 6:
                assertEquals("containsKeyInt" + at, ref.containsKey(bk), mine.containsKeyInt(k));
                break;
            default:
                assertEquals("containsKey" + at, ref.containsKey(bk), mine.containsKey(bk));
                break;
            }
            assertEquals("size" + at, ref.size(), mine.size());
            assertEquals("isEmpty" + at, ref.isEmpty(), mine.isEmpty());
        }
        assertSameContents("after " + STEPS + " steps", ref, mine);
    }

    @Test
    public void theViewsAgreeWithTheJdkMap() {
        IntIntHashMap mine = new IntIntHashMap(NO_KEY);
        Map<Integer, Integer> ref = new HashMap<Integer, Integer>();
        long lcg = 55555L;
        for (int i = 0; i < 400; ++i) {
            lcg = next(lcg);
            int k = KEYS[(int) ((lcg >>> 32) % KEYS.length)];
            lcg = next(lcg);
            int v = 1 + (int) ((lcg >>> 32) % 100L);
            mine.putInt(k, v);
            ref.put(Integer.valueOf(k), Integer.valueOf(v));
        }
        assertTrue("equals the JDK map", mine.equals(ref));
        assertTrue("the JDK map equals it", ref.equals(mine));
        assertEquals("hashCode", ref.hashCode(), mine.hashCode());
        assertArrayEquals("keySet", sorted(ref.keySet()), sorted(mine.keySet()));
        assertArrayEquals("values", sorted(ref.values()), sorted(mine.values()));
        assertEquals("keySet size", ref.size(), mine.keySet().size());
        assertEquals("values size", ref.size(), mine.values().size());
        assertEquals("entrySet size", ref.size(), mine.entrySet().size());
        for (Map.Entry<Integer, Integer> e : mine.entrySet()) {
            assertEquals("entry " + e.getKey(), ref.get(e.getKey()), e.getValue());
        }
        for (Integer k : mine.keySet()) {
            assertTrue("keySet contains " + k, mine.keySet().contains(k));
        }
        assertTrue("containsValue", mine.containsValueInt(ref.values().iterator().next().intValue()));
        assertFalse("containsValue of an absent value", mine.containsValueInt(-1));

        IntIntHashMap fromMap = new IntIntHashMap(ref, NO_KEY);
        assertSameContents("the map constructor", ref, fromMap);

        mine.clear();
        assertEquals("clear", 0, mine.size());
        assertTrue("clear", mine.isEmpty());
    }

    @Test
    public void cloneIsAnIndependentCopy() {
        // clone() used to dereference a result that the compiler cannot prove
        // non-null, because the unreachable catch left it at null
        IntIntHashMap mine = new IntIntHashMap(NO_KEY);
        Map<Integer, Integer> ref = new HashMap<Integer, Integer>();
        long lcg = 777L;
        for (int i = 0; i < 200; ++i) {
            lcg = next(lcg);
            int k = KEYS[(int) ((lcg >>> 32) % KEYS.length)];
            mine.putInt(k, i + 1);
            ref.put(Integer.valueOf(k), Integer.valueOf(i + 1));
        }
        IntIntHashMap copy = (IntIntHashMap) mine.clone();
        assertSameContents("the clone", ref, copy);

        copy.putInt(4711, 42);
        copy.putInt(KEYS[0], 999);
        assertFalse("the clone does not write into the original", mine.containsKeyInt(4711));
        assertEquals("the clone does not write into the original", ref.get(Integer.valueOf(KEYS[0])),
                mine.get(Integer.valueOf(KEYS[0])));
        assertSameContents("the original", ref, mine);

        IntIntHashMap empty = new IntIntHashMap(NO_KEY);
        IntIntHashMap emptyCopy = (IntIntHashMap) empty.clone();
        assertEquals("an empty map clones to an empty map", 0, emptyCopy.size());
        emptyCopy.putInt(1, 1);
        assertEquals("and the two stay independent", 0, empty.size());
    }

    @Test
    public void theEntriesSupportSetValue() {
        // entries obtained from entrySet() are required to support setValue,
        // and this map supports put, so refusing it was the one part of the
        // Map.Entry contract not met
        IntIntHashMap mine = new IntIntHashMap(NO_KEY);
        Map<Integer, Integer> ref = new HashMap<Integer, Integer>();
        for (int i = 1; i <= 20; ++i) {
            mine.putInt(i, i);
            ref.put(Integer.valueOf(i), Integer.valueOf(i));
        }
        for (Map.Entry<Integer, Integer> e : mine.entrySet()) {
            Integer old = e.setValue(Integer.valueOf(e.getValue().intValue() + 100));
            assertEquals("the old value comes back", ref.get(e.getKey()), old);
        }
        for (Map.Entry<Integer, Integer> e : ref.entrySet()) {
            e.setValue(Integer.valueOf(e.getValue().intValue() + 100));
        }
        assertSameContents("after setValue through the entry set", ref, mine);
        assertEquals("the size is untouched", 20, mine.size());
        assertEquals("the entry stays in its bucket", 120, mine.getInt(20));
        assertTrue("the new value is found", mine.containsValueInt(120));
        assertFalse("the old value is gone", mine.containsValueInt(20));

        // Map.Entry.hashCode() is key ^ value, so it follows the new value,
        // and the hash of the map follows with it. That is prescribed, not
        // incidental, and the JDK map does the same
        Map.Entry<Integer, Integer> twenty = null;
        for (Map.Entry<Integer, Integer> e : mine.entrySet()) {
            if (e.getKey().intValue() == 20) {
                twenty = e;
            }
        }
        assertEquals("the entry hash follows the value", 20 ^ 120, twenty.hashCode());
        assertEquals("the map hash follows too", ref.hashCode(), mine.hashCode());

        // what must not follow the value is the placement: the stored hash an
        // entry is bucketed by comes from its key alone. Carry the changed
        // entries through several rehashes and look every one of them up
        // again -- cover for transfer, which nothing else here reaches with
        // values that were changed after insertion
        for (int i = 21; i <= 400; ++i) {
            mine.putInt(i, i);
            ref.put(Integer.valueOf(i), Integer.valueOf(i));
        }
        assertSameContents("after several rehashes", ref, mine);
        assertEquals("the map hash still agrees", ref.hashCode(), mine.hashCode());

        // the map permits no null values, and Map.Entry provides for exactly
        // this refusal
        Map.Entry<Integer, Integer> first = mine.entrySet().iterator().next();
        int untouched = first.getValue().intValue();
        try {
            first.setValue(null);
            fail("a null value has to be refused");
        } catch (NullPointerException expected) {
            assertEquals("the entry is unchanged", untouched, first.getValue().intValue());
            assertEquals("the map is unchanged", ref.size(), mine.size());
        }
    }

    // ----- helpers -----------------------------------------------------------

    private static long next(long lcg) {
        return lcg * 6364136223846793005L + 1442695040888963407L;
    }

    private static int[] pool() {
        int range = 24;
        int[] p = new int[2 * range + 3];
        for (int i = 0; i <= 2 * range; ++i) {
            p[i] = i - range;
        }
        p[2 * range + 1] = Integer.MIN_VALUE;
        p[2 * range + 2] = Integer.MAX_VALUE;
        return p;
    }

    /** What the map answers where the JDK map answers {@code null}. */
    private static int unbox(Integer v) {
        return v == null ? NO_KEY : v.intValue();
    }

    private static void assertSameContents(String what, Map<Integer, Integer> ref, IntIntHashMap mine) {
        assertEquals(what + " size", ref.size(), mine.size());
        assertArrayEquals(what + " keys", sorted(ref.keySet()), sorted(mine.keySet()));
        for (Map.Entry<Integer, Integer> e : ref.entrySet()) {
            assertEquals(what + " value at " + e.getKey(), e.getValue().intValue(),
                    mine.getInt(e.getKey().intValue()));
        }
    }

    private static int[] sorted(Iterable<Integer> values) {
        List<Integer> l = new ArrayList<Integer>();
        for (Integer v : values) {
            l.add(v);
        }
        int[] a = new int[l.size()];
        for (int i = 0; i < a.length; ++i) {
            a[i] = l.get(i).intValue();
        }
        Arrays.sort(a);
        return a;
    }
}
