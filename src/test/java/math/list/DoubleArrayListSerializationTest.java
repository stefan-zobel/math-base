package math.list;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.ConcurrentModificationException;

import org.junit.Test;

/**
 * Tests the {@link java.io.Externalizable} half of {@link DoubleArrayList}: that
 * the stream carries the elements and nothing else, and that a stream which does
 * not agree with itself is refused instead of producing a list whose every
 * accessor is out of bounds.
 */
public class DoubleArrayListSerializationTest {

    private static final double[] AWKWARD = { 0.0, -0.0, Double.NaN, 1.0, -1.0, 2.0, 3.0, Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY, Double.MIN_VALUE };

    private long lcg = 20260824L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return (lcg >>> 11) * 0x1.0p-53;
    }

    @Test
    public void aRoundTripPreservesTheContents() throws Exception {
        for (int trial = 0; trial < 500; ++trial) {
            double[] values = awkwardValues((int) (next() * 30));
            DoubleArrayList list = DoubleArrayList.of(values);
            DoubleArrayList back = roundTrip(list);
            String message = Arrays.toString(values);
            assertEquals(message, list.size(), back.size());
            // bit for bit, so a NaN and a signed zero have to come back as
            // themselves
            assertTrue(message, list.equals(back));
            assertTrue(message, Arrays.equals(values, back.toArray()));
        }
    }

    @Test
    public void anEmptyListRoundTrips() throws Exception {
        DoubleArrayList back = roundTrip(new DoubleArrayList());
        assertEquals(0, back.size());
        assertTrue(back.isEmpty());
        assertEquals(0, back.toArray().length);
    }

    @Test
    public void aClearedListDoesNotShipWhatItHeld() throws Exception {
        DoubleArrayList list = new DoubleArrayList();
        for (int i = 0; i < 1000; ++i) {
            list.add(1.0e6 + i);
        }
        list.clear();
        // clear() only resets the size, so the values are still in the array
        assertTrue(list.getArrayUnsafe().length >= 1000);
        byte[] payload = payload(list);
        assertFalse(contains(payload, 1000042.0));
        assertFalse(contains(payload, 1000999.0));
        assertTrue("payload was " + payload.length + " bytes", payload.length < 200);
    }

    @Test
    public void aRemovedElementDoesNotShipEither() throws Exception {
        DoubleArrayList list = new DoubleArrayList();
        for (int i = 0; i < 10; ++i) {
            list.add(i);
        }
        assertTrue(list.remove(9.0));
        assertTrue(list.remove(8.0));
        byte[] payload = payload(list);
        assertTrue(contains(payload, 7.0));
        assertFalse(contains(payload, 8.0));
        assertFalse(contains(payload, 9.0));
    }

    @Test
    public void theStreamGrowsWithTheSizeNotTheCapacity() throws Exception {
        DoubleArrayList roomy = new DoubleArrayList(100000);
        roomy.add(1.0);
        roomy.add(2.0);
        roomy.add(3.0);
        assertEquals(100000, roomy.getArrayUnsafe().length);
        assertTrue(payload(roomy).length < 200);

        DoubleArrayList full = new DoubleArrayList();
        for (int i = 0; i < 1000; ++i) {
            full.add(i);
        }
        assertTrue(payload(full).length > 8000);
    }

    @Test
    public void anEmptyListComesBackAbleToGrow() throws Exception {
        DoubleArrayList back = roundTrip(new DoubleArrayList());
        back.add(1.0);
        // the same jump a freshly constructed list makes, not one element at
        // a time
        assertEquals(10, back.getArrayUnsafe().length);
    }

    @Test
    public void aSizeBeyondTheArrayIsRejected() throws Exception {
        try {
            readPayload(rawPayload(100, new double[] { 1.0, 2.0, 3.0 }));
            fail("a size of 100 for an array of 3 should have been refused");
        } catch (InvalidObjectException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("100"));
        }
    }

    @Test
    public void aNegativeSizeIsRejected() throws Exception {
        try {
            readPayload(rawPayload(-1, new double[] { 1.0, 2.0, 3.0 }));
            fail("a negative size should have been refused");
        } catch (InvalidObjectException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("-1"));
        }
    }

    @Test
    public void aMissingArrayIsRejected() throws Exception {
        try {
            readPayload(rawPayload(0, null));
            fail("a missing array should have been refused");
        } catch (InvalidObjectException expected) {
            // expected
        }
    }

    @Test
    public void aStreamFromTheOlderWriterStillReads() throws Exception {
        // what the writer produced before it stopped emitting the spare
        // capacity: a size smaller than the array that carries it
        DoubleArrayList back = readPayload(rawPayload(2, new double[] { 1.0, 2.0, 3.0, 0.0, 0.0 }));
        assertEquals(2, back.size());
        assertTrue(Arrays.equals(new double[] { 1.0, 2.0 }, back.toArray()));
        back.add(4.0);
        assertTrue(Arrays.equals(new double[] { 1.0, 2.0, 4.0 }, back.toArray()));
    }

    @Test
    public void readExternalResetsTheModificationCount() throws Exception {
        DoubleArrayList list = new DoubleArrayList();
        list.add(1.0);
        list.add(2.0);
        list.add(3.0);
        DListIterator it = list.listIterator();
        readInto(list, rawPayload(2, new double[] { 7.0, 8.0 }));
        assertTrue(Arrays.equals(new double[] { 7.0, 8.0 }, list.toArray()));
        try {
            it.next();
            fail("the iterator should not have survived a wholly replaced list");
        } catch (ConcurrentModificationException expected) {
            // expected
        }
    }

    // ----- helpers ------------------------------------------------------

    private double[] awkwardValues(int length) {
        double[] values = new double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = AWKWARD[(int) (next() * AWKWARD.length)];
        }
        return values;
    }

    /** A full object round trip through the serialization machinery. */
    private static DoubleArrayList roundTrip(DoubleArrayList list) throws Exception {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream(bytes);
        out.writeObject(list);
        out.close();
        ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()));
        DoubleArrayList back = (DoubleArrayList) in.readObject();
        in.close();
        return back;
    }

    /** Just what writeExternal puts on the wire, without the object header. */
    private static byte[] payload(DoubleArrayList list) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream(bytes);
        list.writeExternal(out);
        out.close();
        return bytes.toByteArray();
    }

    /** The same shape, but with a size and an array that need not agree. */
    private static byte[] rawPayload(int size, double[] array) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream(bytes);
        out.writeInt(size);
        out.writeObject(array);
        out.close();
        return bytes.toByteArray();
    }

    private static DoubleArrayList readPayload(byte[] data) throws Exception {
        DoubleArrayList list = new DoubleArrayList();
        readInto(list, data);
        return list;
    }

    private static void readInto(DoubleArrayList list, byte[] data) throws Exception {
        ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(data));
        list.readExternal(in);
        in.close();
    }

    private static boolean contains(byte[] data, double value) {
        long bits = Double.doubleToLongBits(value);
        for (int i = 0; i + 8 <= data.length; ++i) {
            long v = 0L;
            for (int k = 0; k < 8; ++k) {
                v = (v << 8) | (data[i + k] & 0xffL);
            }
            if (v == bits) {
                return true;
            }
        }
        return false;
    }
}
