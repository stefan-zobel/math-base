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
 * <p>
 * The name table below is the whole specification of this routine: LAPACK's
 * {@code SUBNAM} is a one-character type prefix, a two-character matrix group
 * and a three-character routine, and {@code Ilaenv} answers by matching those
 * three pieces. There is no invariant to assert instead of the table, so these
 * are golden values by necessity rather than by preference. Only five of the
 * roughly twenty routine names it knows are ones this library ever asks about;
 * the rest were untested until this file covered them.
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
        // whether IEEE NaN and infinity arithmetic can be trusted, which on a
        // JVM they always can
        assertEquals(1, Ilaenv.ilaenv(10, "DGEQRF", " ", 100, 100, -1, -1));
        assertEquals(1, Ilaenv.ilaenv(11, "DGEQRF", " ", 100, 100, -1, -1));
    }

    // =========================================================================
    // the name table, one assertion per row
    // =========================================================================

    /** ispec 1, the block size. Anything the table does not name answers 1. */
    @Test
    public void testTheBlockSizeTableForEveryNameItKnows() {
        expect(1, "DGETRF", 64);
        expect(1, "DGETRI", 64);
        expect(1, "DGEQRF", 32);
        expect(1, "DGERQF", 32);
        expect(1, "DGELQF", 32);
        expect(1, "DGEQLF", 32);
        expect(1, "DGEHRD", 32);
        expect(1, "DGEBRD", 32);
        expect(1, "DPOTRF", 64);
        expect(1, "DSYTRF", 64);
        expect(1, "DSYTRD", 32);
        expect(1, "DSYGST", 64);
        expect(1, "ZHETRF", 64);
        expect(1, "ZHETRD", 32);
        expect(1, "ZHEGST", 64);
        expect(1, "DTRTRI", 64);
        expect(1, "DLAUUM", 64);
        expect(1, "DSTEBZ", 1);

        // the generate and multiply forms of the orthogonal and unitary
        // factors, over every two-letter factorization they accept
        for (String pair : new String[] { "QR", "RQ", "LQ", "QL", "HR", "TR", "BR" }) {
            expect(1, "DORG" + pair, 32);
            expect(1, "DORM" + pair, 32);
            expect(1, "ZUNG" + pair, 32);
            expect(1, "ZUNM" + pair, 32);
        }
    }

    /** The two banded rows are the only ones that read the dimension arguments. */
    @Test
    public void testTheBandedRowsDependOnTheBandwidth() {
        // DGBTRF reads n4, DPBTRF reads n2, and both switch at 64
        assertEquals(1, Ilaenv.ilaenv(1, "DGBTRF", " ", 200, 200, 200, 64));
        assertEquals(32, Ilaenv.ilaenv(1, "DGBTRF", " ", 200, 200, 200, 65));
        assertEquals(1, Ilaenv.ilaenv(1, "DPBTRF", " ", 200, 64, 200, 200));
        assertEquals(32, Ilaenv.ilaenv(1, "DPBTRF", " ", 200, 65, 200, 200));
    }

    /**
     * ispec 2, the minimum block size. Every row in the table computes 2, which
     * is also the default, so exactly one name answers anything else. It reads
     * like a transcription mistake and is not.
     */
    @Test
    public void testTheMinimumBlockSizeTableHasExactlyOneNonDefaultRow() {
        expect(2, "DSYTRF", 8);
        for (String name : new String[] { "DGETRF", "DGEQRF", "DGERQF", "DGELQF", "DGEQLF", "DGEHRD", "DGEBRD",
                "DGETRI", "DPOTRF", "DSYTRD", "ZHETRD", "DORGQR", "DORMQR", "ZUNGQR", "ZUNMQR", "DTRTRI",
                "DLAUUM", "DSTEBZ", "DGBTRF", "DPBTRF", "DXXXXX" }) {
            expect(2, name, 2);
        }
    }

    /**
     * ispec 3, the crossover point, default 0.
     * <p>
     * The table is not the one ispec 1 and 2 use: there is <b>no multiply form</b>
     * here, so {@code DORGQR} crosses over at 128 while {@code DORMQR} answers
     * 0. That asymmetry is in reference LAPACK and is why {@code Dormqr} never
     * consults ispec 3.
     */
    @Test
    public void testTheCrossoverTableAndItsMissingMultiplyForm() {
        expect(3, "DGEQRF", 128);
        expect(3, "DGERQF", 128);
        expect(3, "DGELQF", 128);
        expect(3, "DGEQLF", 128);
        expect(3, "DGEHRD", 128);
        expect(3, "DGEBRD", 128);
        expect(3, "DSYTRD", 32);
        expect(3, "ZHETRD", 32);
        for (String pair : new String[] { "QR", "RQ", "LQ", "QL", "HR", "TR", "BR" }) {
            expect(3, "DORG" + pair, 128);
            expect(3, "ZUNG" + pair, 128);
            expect(3, "DORM" + pair, 0);
            expect(3, "ZUNM" + pair, 0);
        }
        // and the rows that have no crossover at all
        for (String name : new String[] { "DGETRF", "DGETRI", "DPOTRF", "DSYTRF", "DTRTRI", "DLAUUM", "DSTEBZ" }) {
            expect(3, name, 0);
        }
    }

    /**
     * The type prefix gates whole rows. {@code S} and {@code D} are the real
     * types, {@code C} and {@code Z} the complex ones, and a row that belongs
     * to one pair answers the default for the other.
     */
    @Test
    public void testThePrefixGate() {
        for (int ispec : new int[] { 1, 2, 3 }) {
            for (String rest : new String[] { "GEQRF", "GETRF", "SYTRF", "SYTRD", "HETRD", "ORGQR", "UNGQR",
                    "ORMQR", "POTRF", "TRTRI" }) {
                assertEquals("ispec " + ispec + ", S and D disagree on " + rest,
                        Ilaenv.ilaenv(ispec, "S" + rest, " ", 200, 200, 200, 200),
                        Ilaenv.ilaenv(ispec, "D" + rest, " ", 200, 200, 200, 200));
                assertEquals("ispec " + ispec + ", C and Z disagree on " + rest,
                        Ilaenv.ilaenv(ispec, "C" + rest, " ", 200, 200, 200, 200),
                        Ilaenv.ilaenv(ispec, "Z" + rest, " ", 200, 200, 200, 200));
            }
        }

        // the real-only rows, asked with a complex prefix
        expect(1, "ZSYTRD", 1);
        expect(1, "ZSYGST", 1);
        expect(1, "ZORGQR", 1);
        expect(1, "ZSTEBZ", 1);
        // the complex-only rows, asked with a real prefix
        expect(1, "DHETRF", 1);
        expect(1, "DHETRD", 1);
        expect(1, "DHEGST", 1);
        expect(1, "DUNGQR", 1);
        // and the one row that is not gated at all: SY + TRF answers 64 for
        // every prefix, unlike SY + TRD beside it
        expect(1, "ZSYTRF", 64);
    }

    @Test
    public void testAnUnknownGroupOrRoutineFallsThroughToTheDefault() {
        expect(1, "DXXQRF", 1);
        expect(1, "DGEXXX", 1);
        expect(1, "DORGXX", 1);
        expect(1, "DXXXXX", 1);
        expect(2, "DGEXXX", 2);
        expect(3, "DGEXXX", 0);
    }

    /**
     * The routine uppercases the name before it matches, which is the one
     * branch of that code left after the EBCDIC and Prime encodings were
     * removed and which nothing in this library reaches.
     */
    @Test
    public void testALowercaseNameGivesTheUppercaseAnswer() {
        expect(1, "dgeqrf", 32);
        expect(1, "dGeQrF", 32);
        expect(1, "dgetrf", 64);
        expect(3, "dgeqrf", 128);
        expect(1, "dormqr", 32);
    }

    private static void expect(int ispec, String name, int expected) {
        assertEquals("ilaenv(" + ispec + ", \"" + name + "\")", expected,
                Ilaenv.ilaenv(ispec, name, " ", 200, 200, 200, 200));
    }
}
