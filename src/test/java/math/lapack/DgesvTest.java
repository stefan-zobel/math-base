package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The boundary of {@code Dgesv}, which is one of the three public entry points
 * of this package and had no test of its own.
 * <p>
 * The arithmetic underneath is the same {@code Dgetrf} and {@code Dgetrs} that
 * {@link LuFactorization} runs, and {@code LuFactorizationTest} already asserts
 * the two agree bit for bit. What is asserted here is everything around that:
 * the argument checks, the offsets, and what a matrix that cannot be solved
 * reports.
 */
public class DgesvTest {

    // =========================================================================
    // the answer
    // =========================================================================

    @Test
    public void testTheSolutionSatisfiesTheSystem() {
        for (int n : new int[] { 1, 2, 5, 12, 40 }) {
            double[] a = wellConditioned(n, 90210L + n);
            double[] b = new double[n];
            for (int i = 0; i < n; ++i) {
                b[i] = 1.0 + 0.25 * i;
            }
            double[] x = b.clone();
            assertTrue("order " + n, Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, x, 0, n));

            double scale = 0.0;
            double residual = 0.0;
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int j = 0; j < n; ++j) {
                    sum += a[j * n + i] * x[j];
                }
                residual = Math.max(residual, Math.abs(sum - b[i]));
                scale = Math.max(scale, Math.abs(b[i]));
            }
            assertTrue("order " + n + " left a residual of " + residual, residual <= 1.0e-13 * (scale + 1.0));
        }
    }

    @Test
    public void testSeveralRightHandSidesAtOnceAreTheColumnsOneAtATime() {
        int n = 9;
        int nrhs = 4;
        double[] a = wellConditioned(n, 13579L);
        double[] b = new double[n * nrhs];
        for (int i = 0; i < b.length; ++i) {
            b[i] = 1.0 + (i % 7);
        }

        double[] together = b.clone();
        assertTrue(Dgesv.dgesv(n, nrhs, a.clone(), 0, n, new int[n], 0, together, 0, n));

        for (int j = 0; j < nrhs; ++j) {
            double[] column = new double[n];
            System.arraycopy(b, j * n, column, 0, n);
            assertTrue(Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, column, 0, n));
            for (int i = 0; i < n; ++i) {
                assertEquals("column " + j + ", entry " + i, column[i], together[j * n + i], 0.0);
            }
        }
    }

    @Test
    public void testAnExactlySingularMatrixIsReportedRatherThanThrown() {
        // the second column is twice the first
        double[] a = { 1.0, 2.0, 2.0, 4.0 };
        assertFalse(Dgesv.dgesv(2, 1, a, 0, 2, new int[2], 0, new double[] { 1.0, 1.0 }, 0, 2));
    }

    @Test
    public void testAProblemOfOrderZeroIsNotAnError() {
        assertTrue(Dgesv.dgesv(0, 1, new double[0], 0, 1, new int[0], 0, new double[0], 0, 1));
        // no right hand side still factorizes, so the matrix has to be one
        // that can be factorized -- an all zero matrix is exactly singular
        assertTrue(Dgesv.dgesv(3, 0, wellConditioned(3, 17L), 0, 3, new int[3], 0, new double[0], 0, 3));
    }

    // =========================================================================
    // the offsets, which the pivot length check used to ignore
    // =========================================================================

    /**
     * The check compared {@code ipiv.length} against {@code n} and never looked
     * at the offset, so a pivot array with an offset that put most of its slots
     * past the end passed the check and failed inside {@code Dgetf2} with an
     * {@code ArrayIndexOutOfBoundsException}.
     */
    @Test
    public void testThePivotOffsetIsHonored() {
        int n = 5;
        double[] a = wellConditioned(n, 24680L);
        double[] b = new double[n];

        refuses("ipiv", () -> Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 4, b.clone(), 0, n));
        refuses("ipiv", () -> Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n - 1], 0, b.clone(), 0, n));
        refuses("_ipiv_offset", () -> Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], -1, b.clone(), 0, n));

        // an offset that does leave room is accepted, and gives the same answer
        double[] plain = b.clone();
        for (int i = 0; i < n; ++i) {
            plain[i] = 1.0 + i;
        }
        double[] offsetSolution = plain.clone();
        assertTrue(Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, plain, 0, n));
        assertTrue(Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n + 3], 3, offsetSolution, 0, n));
        for (int i = 0; i < n; ++i) {
            assertEquals("entry " + i, plain[i], offsetSolution[i], 0.0);
        }
    }

    @Test
    public void testOffsetsAndALargerLeadingDimensionGiveTheSameAnswer() {
        int n = 6;
        int lda = n + 4;
        int offset = 7;
        double[] a = wellConditioned(n, 111213L);
        double[] b = new double[n];
        for (int i = 0; i < n; ++i) {
            b[i] = 2.0 - 0.5 * i;
        }

        double[] plain = b.clone();
        assertTrue(Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, plain, 0, n));

        double[] padded = new double[offset + lda * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                padded[offset + j * lda + i] = a[j * n + i];
            }
        }
        double[] paddedB = new double[offset + n];
        System.arraycopy(b, 0, paddedB, offset, n);
        assertTrue(Dgesv.dgesv(n, 1, padded, offset, lda, new int[n + 2], 2, paddedB, offset, n));
        for (int i = 0; i < n; ++i) {
            assertEquals("entry " + i, plain[i], paddedB[offset + i], 0.0);
        }
    }

    // =========================================================================
    // the array lengths
    // =========================================================================

    @Test
    public void testTheArrayLengthsAreChecked() {
        int n = 5;
        double[] a = wellConditioned(n, 1414L);
        double[] b = new double[n];

        refuses("a", () -> Dgesv.dgesv(n, 1, new double[n * n - 1], 0, n, new int[n], 0, b.clone(), 0, n));
        refuses("b", () -> Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, new double[n - 1], 0, n));
        refuses("_a_offset", () -> Dgesv.dgesv(n, 1, a.clone(), -1, n, new int[n], 0, b.clone(), 0, n));
        refuses("_b_offset", () -> Dgesv.dgesv(n, 1, a.clone(), 0, n, new int[n], 0, b.clone(), -1, n));
    }

    /**
     * The bound is {@code lda * (n - 1) + n}, not {@code lda * n}: the last
     * column needs {@code n} entries and not a whole leading dimension of them,
     * so an array packed exactly to the last entry it uses must be accepted.
     */
    @Test
    public void testTheTightBoundDoesNotRejectAPackedArgument() {
        int n = 4;
        int lda = 9;
        int needed = lda * (n - 1) + n;
        double[] a = new double[needed];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                a[j * lda + i] = (i == j) ? 4.0 : 0.5;
            }
        }
        assertTrue("an array of exactly lda * (n - 1) + n must be accepted",
                Dgesv.dgesv(n, 1, a, 0, lda, new int[n], 0, new double[n], 0, n));
        refuses("a", () -> Dgesv.dgesv(n, 1, new double[needed - 1], 0, lda, new int[n], 0, new double[n], 0, n));
    }

    // =========================================================================
    // non-finite input
    // =========================================================================

    /**
     * {@code Idamax} compares with {@code >} against a running maximum, so a
     * {@code NaN} is never selected as a pivot and never displaced as one;
     * {@code Dgetf2}'s singularity test is {@code != 0.0}, which a {@code NaN}
     * passes. The factorization used to report success and hand back a
     * solution of {@code NaN}, and for an infinity a solution of ordinary
     * looking numbers whose residual cannot even be formed.
     */
    @Test
    public void testANonFiniteEntryIsReportedRatherThanPropagated() {
        int n = 4;
        double[] clean = wellConditioned(n, 2718L);
        double[] b = new double[n];
        for (int i = 0; i < n; ++i) {
            b[i] = 1.0;
        }
        assertTrue("the clean problem must succeed",
                Dgesv.dgesv(n, 1, clean.clone(), 0, n, new int[n], 0, b.clone(), 0, n));

        double[] bad = { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
        int[] positions = { 0, 1, n * n - 1 };
        for (int t = 0; t < bad.length; ++t) {
            for (int p = 0; p < positions.length; ++p) {
                double[] a = clean.clone();
                a[positions[p]] = bad[t];
                double[] x = b.clone();
                assertFalse("entry " + positions[p] + " set to " + bad[t] + " must be reported",
                        Dgesv.dgesv(n, 1, a, 0, n, new int[n], 0, x, 0, n));
                for (int i = 0; i < n; ++i) {
                    assertEquals("a rejected call must not have touched the right hand side", b[i], x[i], 0.0);
                }
            }
        }
    }

    @Test
    public void testTheSameHoldsForLuFactorization() {
        int n = 4;
        double[] clean = wellConditioned(n, 3141L);
        assertTrue(LuFactorization.factor(n, clean.clone(), new int[n]));

        double[] bad = { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
        for (int t = 0; t < bad.length; ++t) {
            double[] a = clean.clone();
            a[n + 2] = bad[t];
            assertFalse("a matrix holding " + bad[t] + " must be reported", LuFactorization.factor(n, a, new int[n]));
        }
    }

    // =========================================================================
    // helpers
    // =========================================================================

    /** a diagonally dominant matrix of order n, column-major, from a seeded LCG */
    private static double[] wellConditioned(int n, long seed) {
        double[] a = new double[n * n];
        long state = seed;
        for (int i = 0; i < a.length; ++i) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            a[i] = 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
        }
        for (int i = 0; i < n; ++i) {
            a[i * n + i] += n + 1;
        }
        return a;
    }

    private static void refuses(String expected, Runnable call) {
        try {
            call.run();
            fail("expected a rejection naming " + expected);
        } catch (IllegalArgumentException e) {
            assertTrue("the message should name " + expected + ", but is " + e.getMessage(),
                    e.getMessage().contains(expected));
        }
    }
}
