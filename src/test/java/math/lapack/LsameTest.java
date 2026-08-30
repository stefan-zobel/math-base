package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * {@code Lsame} decides every {@code side}, {@code uplo}, {@code trans} and
 * {@code diag} branch in this package, so its semantics are load bearing well
 * out of proportion to its two lines.
 * <p>
 * It was rewritten to compare the two leading characters directly instead of
 * taking a substring of each, which removed an allocation from a method called
 * at the entry of every routine and inside the loops of {@code Dlarft} and
 * {@code Dlarfb}. These assertions pin the semantics that rewrite had to
 * preserve.
 */
public class LsameTest {

    /** exactly what the old form, s1.substring(0, 1).equalsIgnoreCase(...), computed */
    private static boolean reference(String s, String s1) {
        return s1.substring(0, 1).equalsIgnoreCase(s.substring(0, 1));
    }

    private static void agrees(String s, String s1) {
        assertEquals("lsame(\"" + s + "\", \"" + s1 + "\")", Boolean.valueOf(reference(s, s1)),
                Boolean.valueOf(Lsame.lsame(s, s1)));
    }

    /**
     * Every argument pair this package actually forms, taken from the call
     * sites, against the form that was replaced.
     */
    @Test
    public void testItAgreesWithTheSubstringFormOnEveryPairThePackageUses() {
        String[] arguments = { "L", "R", "N", "T", "C", "U", "F", "B", "M", "I", "O", "E", "G", "H", "Q", "Z", "S",
                "P", "1", "Left", "Right", "No transpose", "Transpose", "Non-unit", "Unit", "Upper", "Lower",
                "Forward", "Backward", "Columnwise", "Rowwise", "Full", " " };
        for (int i = 0; i < arguments.length; ++i) {
            for (int j = 0; j < arguments.length; ++j) {
                agrees(arguments[i], arguments[j]);
            }
        }
    }

    @Test
    public void testTheComparisonIgnoresCase() {
        assertTrue(Lsame.lsame("Left", "L"));
        assertTrue(Lsame.lsame("left", "L"));
        assertTrue(Lsame.lsame("LEFT", "l"));
        assertTrue(Lsame.lsame("l", "L"));
        assertTrue(Lsame.lsame("No transpose", "n"));
        assertTrue(Lsame.lsame("Non-unit", "N"));
    }

    @Test
    public void testDifferentLettersDoNotMatch() {
        assertFalse(Lsame.lsame("Left", "R"));
        assertFalse(Lsame.lsame("Upper", "L"));
        assertFalse(Lsame.lsame("Transpose", "N"));
        // "Non-unit" and "No transpose" both begin with N, which is why the
        // routines that need to tell them apart never ask Lsame to
        assertTrue(Lsame.lsame("Non-unit", "No transpose"));
    }

    @Test
    public void testOnlyTheFirstCharacterIsRead() {
        assertTrue(Lsame.lsame("Lower", "Left"));
        assertTrue(Lsame.lsame("Transpose", "T is all that counts"));
        assertTrue(Lsame.lsame("Upper", "upside down, but it starts with a u"));
    }

    @Test
    public void testNonLettersCompareByIdentity() {
        assertTrue(Lsame.lsame("1", "1"));
        assertFalse(Lsame.lsame("1", "2"));
        assertTrue(Lsame.lsame(" ", " "));
        assertFalse(Lsame.lsame(" ", "N"));
    }
}
