package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * Tests for {@link Ilaenv}, the block size oracle.
 * <p>
 * It decomposes the routine name as {@code name[1..2]} and {@code name[3..5]},
 * so a name shorter than six characters used to leave the method with a
 * {@code StringIndexOutOfBoundsException}. Every current call site passes a
 * six-character literal, which is why the trap stayed latent.
 */
public class IlaenvTest {

    @Test
    public void testAShortNameIsTreatedAsUnrecognized() {
        // 1 is what the routine already answers for any name whose first
        // letter is not one of S, D, C or Z, so a name it cannot even parse
        // reaches the same place instead of throwing
        for (String name : new String[] { "", "D", "DG", "DGE", "DGEL", "DGELS" }) {
            assertEquals("name \"" + name + "\"", 1, Ilaenv.ilaenv(1, name, " ", 100, 100, -1, -1));
        }
    }

    @Test
    public void testANullNameIsTreatedAsUnrecognized() {
        assertEquals(1, Ilaenv.ilaenv(1, null, " ", 100, 100, -1, -1));
    }

    @Test
    public void testAnUnrecognizedPrefixIsStillUnrecognized() {
        // the pre-existing path, unchanged by the length guard
        assertEquals(1, Ilaenv.ilaenv(1, "XGEQRF", " ", 100, 100, -1, -1));
    }

    @Test
    public void testTheBlockSizesForTheNamesThisLibraryActuallyUses() {
        // pinned because the whole blocked/unblocked split in Dgeqrf, Dgelqf,
        // Dormqr and Dormlq hangs off these numbers
        assertEquals(32, Ilaenv.ilaenv(1, "DGEQRF", " ", 200, 64, -1, -1));
        assertEquals(32, Ilaenv.ilaenv(1, "DGELQF", " ", 64, 200, -1, -1));
        assertEquals(32, Ilaenv.ilaenv(1, "DORMQR", "LN", 200, 5, 64, -1));
        assertEquals(32, Ilaenv.ilaenv(1, "DORMLQ", "LT", 200, 5, 64, -1));
        assertEquals(64, Ilaenv.ilaenv(1, "DGETRF", " ", 200, 200, -1, -1));
        // the minimum block size and the crossover point
        assertEquals(2, Ilaenv.ilaenv(2, "DGEQRF", " ", 200, 64, -1, -1));
        assertEquals(128, Ilaenv.ilaenv(3, "DGEQRF", " ", 200, 64, -1, -1));
    }

    @Test
    public void testTheOptsArgumentIsIgnoredForTheBlockSizeQueries() {
        // Dormqr and Dormlq used to build this argument as side + trans with
        // trans an enum, yielding "LTRANS" rather than the two characters
        // LAPACK specifies. That was harmless only because ilaenv never reads
        // opts for ispec 1, 2 or 3, and this is what pins that
        for (String opts : new String[] { "LN", "LT", "RN", "RT", "LTRANS", " ", "" }) {
            assertEquals("opts \"" + opts + "\"", 32, Ilaenv.ilaenv(1, "DORMQR", opts, 200, 5, 64, -1));
        }
    }

    @Test
    public void testAnUnknownSpecificationIsRejectedWithMinusOne() {
        assertEquals(-1, Ilaenv.ilaenv(0, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(-1, Ilaenv.ilaenv(17, "DGEQRF", " ", 100, 100, -1, -1));
    }

    @Test
    public void testTheConstantSpecificationsAreUnchanged() {
        assertEquals(6, Ilaenv.ilaenv(4, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(2, Ilaenv.ilaenv(5, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(1, Ilaenv.ilaenv(7, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(50, Ilaenv.ilaenv(8, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(25, Ilaenv.ilaenv(9, "DGEQRF", " ", 100, 100, -1, -1));
        assertTrue(Ilaenv.ilaenv(6, "DGEQRF", " ", 100, 200, -1, -1) > 0);
    }
}
