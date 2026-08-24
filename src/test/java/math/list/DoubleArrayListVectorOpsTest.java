package math.list;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Pins the three contracts the element-wise vector operations of
 * {@link DoubleArrayList} keep, which are deliberately not the same one:
 * {@code plusn}, {@code minusn} and {@code muln} run over the shorter of the
 * two lists and leave the rest of the longer one alone, {@code dot} refuses a
 * dimension mismatch, and {@code cross} refuses anything but three dimensions.
 * The truncation is documented behavior, so it is not to be tightened into an
 * exception by way of a later clean-up.
 */
public class DoubleArrayListVectorOpsTest {

    @Test
    public void plusnTruncatesToTheShorterList() {
        // this list is the longer one: the tail beyond the argument stays
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        assertSame(list, list.plusn(DoubleArrayList.of(10.0, 20.0)));
        assertArrayEquals(new double[] { 11.0, 22.0, 3.0, 4.0, 5.0 }, list.toArray(), 0.0);

        // and the argument is the longer one: what it carries beyond this
        // list is dropped
        list = DoubleArrayList.of(1.0, 2.0);
        assertSame(list, list.plusn(DoubleArrayList.of(10.0, 20.0, 30.0, 40.0)));
        assertArrayEquals(new double[] { 11.0, 22.0 }, list.toArray(), 0.0);
    }

    @Test
    public void minusnTruncatesToTheShorterList() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        assertSame(list, list.minusn(DoubleArrayList.of(10.0, 20.0)));
        assertArrayEquals(new double[] { -9.0, -18.0, 3.0, 4.0, 5.0 }, list.toArray(), 0.0);

        list = DoubleArrayList.of(1.0, 2.0);
        assertSame(list, list.minusn(DoubleArrayList.of(10.0, 20.0, 30.0, 40.0)));
        assertArrayEquals(new double[] { -9.0, -18.0 }, list.toArray(), 0.0);
    }

    @Test
    public void mulnTruncatesToTheShorterList() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        assertSame(list, list.muln(DoubleArrayList.of(2.0, 3.0)));
        assertArrayEquals(new double[] { 2.0, 6.0, 3.0, 4.0, 5.0 }, list.toArray(), 0.0);

        list = DoubleArrayList.of(1.0, 2.0);
        assertSame(list, list.muln(DoubleArrayList.of(2.0, 3.0, 4.0, 5.0)));
        assertArrayEquals(new double[] { 2.0, 6.0 }, list.toArray(), 0.0);
    }

    @Test
    public void theArrayOverloadsTruncateTheSameWay() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        list.plusn(new double[] { 10.0, 20.0 });
        assertArrayEquals(new double[] { 11.0, 22.0, 3.0, 4.0, 5.0 }, list.toArray(), 0.0);

        list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        list.minusn(new double[] { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0 });
        assertArrayEquals(new double[] { -9.0, -18.0, -27.0, -36.0, -45.0 }, list.toArray(), 0.0);

        list = DoubleArrayList.of(1.0, 2.0, 3.0);
        list.muln(new double[] { 2.0 });
        assertArrayEquals(new double[] { 2.0, 2.0, 3.0 }, list.toArray(), 0.0);
    }

    @Test
    public void anEmptyArgumentChangesNothing() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0);
        list.plusn(DoubleArrayList.of());
        list.minusn(new double[0]);
        list.muln(DoubleArrayList.of());
        assertArrayEquals(new double[] { 1.0, 2.0, 3.0 }, list.toArray(), 0.0);

        // the other way round as well: an empty receiver has nothing to write
        DoubleList empty = DoubleArrayList.of();
        empty.plusn(DoubleArrayList.of(1.0, 2.0));
        assertEquals(0, empty.size());
    }

    @Test
    public void aSubListTruncatesAndLeavesTheRootAlone() {
        DoubleArrayList root = ascending(10);
        DoubleList window = root.subList(3, 7);
        assertSame(window, window.plusn(DoubleArrayList.of(100.0, 200.0)));
        assertArrayEquals(new double[] { 0.0, 1.0, 2.0, 103.0, 204.0, 5.0, 6.0, 7.0, 8.0, 9.0 }, root.toArray(), 0.0);

        // an argument longer than the window stops at the window's end, it
        // does not run on into the root behind it
        root = ascending(10);
        window = root.subList(3, 7);
        window.minusn(new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
        assertArrayEquals(new double[] { 0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 7.0, 8.0, 9.0 }, root.toArray(), 0.0);

        root = ascending(10);
        window = root.subList(3, 7);
        window.muln(DoubleArrayList.of(2.0, 2.0, 2.0, 2.0, 2.0, 2.0));
        assertArrayEquals(new double[] { 0.0, 1.0, 2.0, 6.0, 8.0, 10.0, 12.0, 7.0, 8.0, 9.0 }, root.toArray(), 0.0);
    }

    @Test
    public void aSubListArgumentIsReadFromItsOwnOffset() {
        DoubleList argument = ascending(10).subList(5, 8);
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0);
        list.plusn(argument);
        assertArrayEquals(new double[] { 6.0, 8.0, 10.0, 4.0 }, list.toArray(), 0.0);
    }

    @Test
    public void dotRefusesADimensionMismatch() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0);
        assertEquals(32.0, list.dot(DoubleArrayList.of(4.0, 5.0, 6.0)), 0.0);
        assertEquals(32.0, list.dot(new double[] { 4.0, 5.0, 6.0 }), 0.0);
        // a bilinear form over two vectors is undefined unless the dimensions
        // agree, which is why this one throws where plusn truncates
        expectIllegalArgument(list, DoubleArrayList.of(4.0, 5.0));
        expectIllegalArgument(list, DoubleArrayList.of(4.0, 5.0, 6.0, 7.0));
        try {
            list.dot(new double[] { 4.0, 5.0 });
            fail("no exception for a dot product of length 3 with length 2");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        // and through a sub-list, which shares the rule
        DoubleList window = ascending(10).subList(2, 5);
        assertEquals(2.0 * 1.0 + 3.0 * 2.0 + 4.0 * 3.0, window.dot(DoubleArrayList.of(1.0, 2.0, 3.0)), 0.0);
        expectIllegalArgument(window, DoubleArrayList.of(1.0, 2.0));
    }

    @Test
    public void crossRefusesAnythingButThree() {
        DoubleArrayList x = DoubleArrayList.of(1.0, 0.0, 0.0);
        DoubleArrayList y = DoubleArrayList.of(0.0, 1.0, 0.0);
        assertArrayEquals(new double[] { 0.0, 0.0, 1.0 }, x.cross(y).toArray(), 0.0);
        assertArrayEquals(new double[] { 0.0, 0.0, 1.0 }, x.cross(new double[] { 0.0, 1.0, 0.0 }).toArray(), 0.0);
        // which side is wrong decides the exception, and the interface
        // documents both: a receiver that is not three-dimensional does not
        // have the operation, a wrong argument is a wrong argument
        try {
            DoubleArrayList.of(1.0, 2.0).cross(y);
            fail("no exception for a cross product of a vector of length 2");
        } catch (UnsupportedOperationException expected) {
            // expected
        }
        try {
            DoubleArrayList.of(1.0, 2.0).cross(new double[] { 0.0, 1.0, 0.0 });
            fail("no exception for a cross product of a vector of length 2");
        } catch (UnsupportedOperationException expected) {
            // expected
        }
        try {
            x.cross(DoubleArrayList.of(1.0, 2.0, 3.0, 4.0));
            fail("no exception for a cross product with a vector of length 4");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            x.cross(new double[] { 1.0, 2.0 });
            fail("no exception for a cross product with an array of length 2");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    @Test
    public void aNullArgumentIsRefused() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0);
        expectNullPointer(list, (DoubleList) null);
        try {
            list.plusn((double[]) null);
            fail("no exception for a null array");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            list.minusn((double[]) null);
            fail("no exception for a null array");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            list.muln((double[]) null);
            fail("no exception for a null array");
        } catch (NullPointerException expected) {
            // expected
        }
        expectNullPointer(ascending(6).subList(1, 4), (DoubleList) null);
    }

    // ----- helpers ------------------------------------------------------

    private static DoubleArrayList ascending(int size) {
        DoubleArrayList list = new DoubleArrayList(size);
        for (int i = 0; i < size; ++i) {
            list.add(i);
        }
        return list;
    }

    private static void expectIllegalArgument(DoubleList list, DoubleList argument) {
        try {
            list.dot(argument);
            fail("no exception for a dot product of length " + list.size() + " with length " + argument.size());
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    private static void expectNullPointer(DoubleList list, DoubleList argument) {
        try {
            list.plusn(argument);
            fail("no exception for a null list in plusn");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            list.minusn(argument);
            fail("no exception for a null list in minusn");
        } catch (NullPointerException expected) {
            // expected
        }
        try {
            list.muln(argument);
            fail("no exception for a null list in muln");
        } catch (NullPointerException expected) {
            // expected
        }
    }
}
