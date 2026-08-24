package math.list;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;

/**
 * Tests that {@link DoubleArrayList#getArrayUnsafe()} leaves the list alone.
 * Handing out the backing array used to replace the empty sentinel, which cost
 * the list its default-capacity marker: the first {@code add} then grew the
 * array to one element instead of ten, and the whole growth sequence became
 * {@code 1 2 3 4 6 9 ...} instead of {@code 10 15 22 ...}.
 */
public class DoubleArrayListCapacityTest {

    /** What a default-sized list is inflated to by its first element. */
    private static final int DEFAULT_CAPACITY = 10;

    @Test
    public void aFreshListGrowsToTheDefaultCapacity() {
        DoubleArrayList list = new DoubleArrayList();
        list.add(1.0);
        assertEquals(DEFAULT_CAPACITY, list.getArrayUnsafe().length);
    }

    @Test
    public void getArrayUnsafeKeepsTheDefaultCapacityMarker() {
        DoubleArrayList list = new DoubleArrayList();
        assertNotNull(list.getArrayUnsafe());
        assertEquals(0, list.getArrayUnsafe().length);
        list.add(1.0);
        assertEquals(DEFAULT_CAPACITY, list.getArrayUnsafe().length);
    }

    @Test
    public void theArgumentOfAVectorOperationKeepsItsMarker() {
        // these reach for the argument's array, not the receiver's
        DoubleArrayList argument = new DoubleArrayList();
        new DoubleArrayList().dot(argument);
        argument.add(1.0);
        assertEquals("dot", DEFAULT_CAPACITY, argument.getArrayUnsafe().length);

        DoubleArrayList summand = new DoubleArrayList();
        new DoubleArrayList().plusn(summand);
        summand.add(1.0);
        assertEquals("plusn", DEFAULT_CAPACITY, summand.getArrayUnsafe().length);

        DoubleArrayList first = new DoubleArrayList();
        DoubleArrayList second = new DoubleArrayList();
        DoubleList.approxEqual(first, second, 1.0e-9, 1.0e-12);
        first.add(1.0);
        second.add(1.0);
        assertEquals("approxEqual, first", DEFAULT_CAPACITY, first.getArrayUnsafe().length);
        assertEquals("approxEqual, second", DEFAULT_CAPACITY, second.getArrayUnsafe().length);
    }

    @Test
    public void aSubListDoesNotCostTheRootListItsMarker() {
        DoubleArrayList root = new DoubleArrayList();
        root.subList(0, 0).getArrayUnsafe();
        root.add(1.0);
        assertEquals(DEFAULT_CAPACITY, root.getArrayUnsafe().length);
    }

    @Test
    public void anExplicitlyEmptyListIsStillNotADefaultSizedOne() {
        // the two empty sentinels mean different things, and the difference
        // has to survive: a capacity asked for is a capacity meant
        DoubleArrayList list = new DoubleArrayList(0);
        list.getArrayUnsafe();
        list.add(1.0);
        assertEquals(1, list.getArrayUnsafe().length);
    }

    @Test
    public void theArrayHandedOutIsTheOneTheListUses() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0);
        double[] backing = list.getArrayUnsafe();
        backing[1] = 7.0;
        assertEquals(7.0, list.get(1), 0.0);
    }
}
