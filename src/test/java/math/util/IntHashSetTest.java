package math.util;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.junit.Test;

/**
 * {@link IntHashSet} against {@link HashSet}, of which it is a primitive
 * translation: the same sequence of operations has to answer the same at every
 * step and leave both containers in the same state.
 */
public class IntHashSetTest {

    /** Sentinels the set has to cope with, both extremes included. */
    private static final int[] SENTINELS = { Integer.MIN_VALUE, Integer.MAX_VALUE, 1000, -1000 };

    /** Values are drawn from a small pool so that hits and collisions are frequent. */
    private static final int[] VALUES = pool();

    private static final int STEPS = 20000;

    // ----- the differential test ---------------------------------------------

    @Test
    public void randomOperationsAgreeWithTheJdkSet() {
        for (int s = 0; s < SENTINELS.length; ++s) {
            int sentinel = SENTINELS[s];
            IntHashSet mine = new IntHashSet(sentinel);
            Set<Integer> ref = new HashSet<Integer>();
            long lcg = 987654321L + s;
            for (int step = 0; step < STEPS; ++step) {
                lcg = next(lcg);
                int op = (int) ((lcg >>> 32) % 6L);
                lcg = next(lcg);
                int v = VALUES[(int) ((lcg >>> 32) % VALUES.length)];
                if (v == sentinel) {
                    // the sentinel is by contract not a value of the set
                    continue;
                }
                Integer boxed = Integer.valueOf(v);
                String at = " (sentinel " + sentinel + ", step " + step + ", value " + v + ")";
                switch (op) {
                case 0:
                    assertEquals("addInt" + at, ref.add(boxed), mine.addInt(v));
                    break;
                case 1:
                    assertEquals("add" + at, ref.add(boxed), mine.add(boxed));
                    break;
                case 2:
                    assertEquals("removeInt" + at, ref.remove(boxed), mine.removeInt(v));
                    break;
                case 3:
                    assertEquals("remove" + at, ref.remove(boxed), mine.remove(boxed));
                    break;
                case 4:
                    assertEquals("containsInt" + at, ref.contains(boxed), mine.containsInt(v));
                    break;
                default:
                    assertEquals("contains" + at, ref.contains(boxed), mine.contains(boxed));
                    break;
                }
                assertEquals("size" + at, ref.size(), mine.size());
                assertEquals("isEmpty" + at, ref.isEmpty(), mine.isEmpty());
            }
            assertSameContents("sentinel " + sentinel, ref, mine);
        }
    }

    @Test
    public void theBulkOperationsAgreeWithTheJdkSet() {
        for (int s = 0; s < SENTINELS.length; ++s) {
            int sentinel = SENTINELS[s];
            IntHashSet mine = new IntHashSet(sentinel);
            Set<Integer> ref = new HashSet<Integer>();
            long lcg = 24680L + s;
            for (int round = 0; round < 200; ++round) {
                lcg = next(lcg);
                int op = (int) ((lcg >>> 32) % 4L);
                lcg = next(lcg);
                // both branches of removeAll depend on how the two sizes
                // compare, so the argument has to be short as well as long
                int n = (int) ((lcg >>> 32) % 24L);
                List<Integer> arg = new ArrayList<Integer>(n);
                for (int i = 0; i < n; ++i) {
                    lcg = next(lcg);
                    int v = VALUES[(int) ((lcg >>> 32) % VALUES.length)];
                    if (v != sentinel) {
                        arg.add(Integer.valueOf(v));
                    }
                }
                String at = " (sentinel " + sentinel + ", round " + round + ", |arg| " + arg.size() + ")";
                switch (op) {
                case 0:
                    assertEquals("addAll" + at, ref.addAll(arg), mine.addAll(arg));
                    break;
                case 1:
                    assertEquals("removeAll" + at, ref.removeAll(arg), mine.removeAll(arg));
                    break;
                case 2:
                    assertEquals("retainAll" + at, ref.retainAll(arg), mine.retainAll(arg));
                    break;
                default:
                    assertEquals("containsAll" + at, ref.containsAll(arg), mine.containsAll(arg));
                    break;
                }
                assertSameContents("after op " + op + at, ref, mine);
            }
        }
    }

    @Test
    public void theViewsAndCopiesAgreeWithTheJdkSet() {
        IntHashSet mine = new IntHashSet(Integer.MIN_VALUE);
        Set<Integer> ref = new HashSet<Integer>();
        long lcg = 13579L;
        for (int i = 0; i < 400; ++i) {
            lcg = next(lcg);
            int v = VALUES[(int) ((lcg >>> 32) % VALUES.length)];
            if (v != Integer.MIN_VALUE) {
                mine.addInt(v);
                ref.add(Integer.valueOf(v));
            }
        }
        assertTrue("equals the JDK set", mine.equals(ref));
        assertTrue("the JDK set equals it", ref.equals(mine));
        assertEquals("hashCode", ref.hashCode(), mine.hashCode());

        assertArrayEquals("toArray()", sorted(Arrays.asList(ref.toArray())), sorted(Arrays.asList(mine.toArray())));
        assertArrayEquals("toArray(T[0])", sorted(Arrays.asList(ref.toArray(new Integer[0]))),
                sorted(Arrays.asList(mine.toArray(new Integer[0]))));
        Integer[] oversized = mine.toArray(new Integer[mine.size() + 3]);
        assertNull("toArray fills the tail with null", oversized[mine.size()]);

        IntHashSet copy = (IntHashSet) mine.clone();
        assertSameContents("clone", ref, copy);
        // 4711 is outside the pool the two were filled from
        copy.addInt(4711);
        assertFalse("the clone is independent", mine.containsInt(4711));

        IntHashSet fromCollection = new IntHashSet(ref, Integer.MIN_VALUE);
        assertSameContents("the collection constructor", ref, fromCollection);

        // removal through the iterator, which retainAll and removeAll use
        Iterator<Integer> it = mine.iterator();
        Iterator<Integer> refIt = ref.iterator();
        while (it.hasNext()) {
            Integer v = it.next();
            if ((v.intValue() & 1) == 0) {
                it.remove();
            }
        }
        while (refIt.hasNext()) {
            if ((refIt.next().intValue() & 1) == 0) {
                refIt.remove();
            }
        }
        assertSameContents("after iterator removal", ref, mine);

        mine.clear();
        assertEquals("clear", 0, mine.size());
        assertTrue("clear", mine.isEmpty());
    }

    // ----- the defects the differential test would have caught ---------------

	@Test
    public void removingAnAbsentElementAnswersFalse() {
        IntHashSet set = new IntHashSet(Integer.MIN_VALUE);
        set.addInt(1);
        set.addInt(2);
        assertTrue("a present element", set.remove(Integer.valueOf(1)));
        // this used to unbox the null the map answers for an absent key
        assertFalse("an absent element", set.remove(Integer.valueOf(99)));
        assertFalse("an element of the wrong type", set.remove((Object) "99"));
        assertEquals(1, set.size());
    }

    @Test
    public void removeAllSurvivesAnAbsentElement() {
        // the branch taken when the set is larger than the argument routes
        // every element through remove(Object), so a single absent one used to
        // abort the loop with the set already half modified
        IntHashSet set = new IntHashSet(Integer.MIN_VALUE);
        for (int i = 0; i < 10; ++i) {
            set.addInt(i);
        }
        assertTrue(set.removeAll(Arrays.asList(Integer.valueOf(3), Integer.valueOf(77))));
        assertEquals(9, set.size());
        assertFalse(set.containsInt(3));
        assertFalse("nothing was removed twice", set.removeAll(Arrays.asList(Integer.valueOf(3))));
        assertEquals(9, set.size());
    }

    @Test
    public void theSentinelMayBeTheLargestInt() {
        // Integer.MAX_VALUE used to be the marker the set stores for every
        // key, so a caller choosing it as the sentinel got nonsense back from
        // the int taking methods. AbstractRng64 picks exactly this sentinel
        for (int i = 0; i < SENTINELS.length; ++i) {
            int sentinel = SENTINELS[i];
            IntHashSet set = new IntHashSet(sentinel);
            assertTrue("first insert, sentinel " + sentinel, set.addInt(7));
            assertFalse("second insert, sentinel " + sentinel, set.addInt(7));
            assertEquals("size, sentinel " + sentinel, 1, set.size());
            assertTrue("contains, sentinel " + sentinel, set.containsInt(7));
            assertTrue("removal, sentinel " + sentinel, set.removeInt(7));
            assertFalse("second removal, sentinel " + sentinel, set.removeInt(7));
            assertEquals("size after removal, sentinel " + sentinel, 0, set.size());
        }
    }

    // ----- helpers -----------------------------------------------------------

    private static long next(long lcg) {
        return lcg * 6364136223846793005L + 1442695040888963407L;
    }

    private static int[] pool() {
        int range = 24;
        int[] p = new int[2 * range + 4];
        for (int i = 0; i <= 2 * range; ++i) {
            p[i] = i - range;
        }
        p[2 * range + 1] = Integer.MIN_VALUE;
        p[2 * range + 2] = Integer.MAX_VALUE;
        p[2 * range + 3] = Integer.MIN_VALUE + 1;
        return p;
    }

    private static void assertSameContents(String what, Set<Integer> ref, IntHashSet mine) {
        assertEquals(what + " size", ref.size(), mine.size());
        assertArrayEquals(what + " contents", sortedInts(ref), sortedInts(mine));
    }

    private static int[] sortedInts(Iterable<Integer> values) {
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

    private static int[] sorted(List<?> values) {
        int[] a = new int[values.size()];
        for (int i = 0; i < a.length; ++i) {
            a[i] = ((Integer) values.get(i)).intValue();
        }
        Arrays.sort(a);
        return a;
    }
}
