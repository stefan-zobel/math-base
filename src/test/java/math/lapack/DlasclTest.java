package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import org.junit.Test;

import math.gemm.Trans;

/**
 * Tests for {@link Dlascl}, whose scaling loop used not to terminate.
 * <p>
 * The loop shrinks {@code cfromc} by the smallest normal number, or
 * {@code ctoc} by the largest, until the two are close enough for a single
 * multiplication to carry the scaling. That terminates only while the
 * magnitudes strictly decrease, which fails for an infinite argument: the
 * product stays infinite and the same branch is taken forever. Every test here
 * carries a timeout, so a regression shows up as a failure rather than as a
 * stuck build.
 */
public class DlasclTest {

    private static final long TIMEOUT = 10000L;

    /** Xerbla writes to stderr as well as throwing; keep the output quiet. */
    private static void assertRejects(double cfrom, double cto, int expectedArgument) {
        PrintStream saved = System.err;
        System.setErr(new PrintStream(new ByteArrayOutputStream()));
        try {
            double[] a = { 1.0, 2.0, 3.0, 4.0 };
            Dlascl.dlascl("G", 0, 0, cfrom, cto, 2, 2, a, 0, 2, new intW(0));
            fail("expected a rejection for cfrom = " + cfrom + ", cto = " + cto);
        } catch (IllegalArgumentException expected) {
            assertTrue("the message should name argument " + expectedArgument + ", was: " + expected.getMessage(),
                    expected.getMessage().contains("'" + expectedArgument + "'"));
        } finally {
            System.setErr(saved);
        }
    }

    @Test(timeout = TIMEOUT)
    public void testANonFiniteScaleSourceIsRejected() {
        // argument 4 is CFROM
        assertRejects(Double.POSITIVE_INFINITY, 1.0, 4);
        assertRejects(Double.NEGATIVE_INFINITY, 1.0, 4);
        assertRejects(Double.NaN, 1.0, 4);
        // the pre-existing zero test still applies
        assertRejects(0.0, 1.0, 4);
    }

    @Test(timeout = TIMEOUT)
    public void testANonFiniteScaleTargetIsRejected() {
        // argument 5 is CTO
        assertRejects(1.0, Double.POSITIVE_INFINITY, 5);
        assertRejects(1.0, Double.NEGATIVE_INFINITY, 5);
        assertRejects(1.0, Double.NaN, 5);
    }

    @Test(timeout = TIMEOUT)
    public void testTheScalingItselfIsUnchanged() {
        // the guard must not disturb what the routine is for. A plain scaling
        // multiplies every entry by cto/cfrom
        double[] a = { 1.0, -2.0, 3.5, 4.0, -0.25, 8.0 };
        double[] expected = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            expected[i] = a[i] * (7.0 / 2.0);
        }
        Dlascl.dlascl("G", 0, 0, 2.0, 7.0, 3, 2, a, 0, 3, new intW(0));
        for (int i = 0; i < a.length; ++i) {
            assertEquals(expected[i], a[i], 0.0);
        }
    }

    @Test(timeout = TIMEOUT)
    public void testTheMultiPassScalingStillTerminatesAndIsAccurate() {
        // cfrom and cto far enough apart that cto/cfrom would itself overflow,
        // which is the whole reason the loop exists. It has to take several
        // passes and still land on the right answer
        double[] a = { 1.0, 2.0, 4.0, 8.0 };
        Dlascl.dlascl("G", 0, 0, 1.0e-300, 1.0e300, 2, 2, a, 0, 2, new intW(0));
        // 1 * 1e600 overflows, so the entries saturate at infinity; the point
        // is that the call returns at all
        for (int i = 0; i < a.length; ++i) {
            assertTrue("entry " + i + " = " + a[i], a[i] > 0.0);
        }

        // and a round trip that stays in range has to come back where it began
        double[] b = { 3.0, -5.0, 7.0, 11.0 };
        double[] original = b.clone();
        Dlascl.dlascl("G", 0, 0, 1.0, 1.0e-290, 2, 2, b, 0, 2, new intW(0));
        Dlascl.dlascl("G", 0, 0, 1.0e-290, 1.0, 2, 2, b, 0, 2, new intW(0));
        for (int i = 0; i < b.length; ++i) {
            assertEquals(original[i], b[i], 1.0e-13 * Math.abs(original[i]));
        }
    }

    // =========================================================================
    // through the public entry points, which is how the hang was reachable
    // =========================================================================

    private static double[] design(int m, int n) {
        double[] a = new double[m * n];
        long state = 314159265L;
        for (int i = 0; i < a.length; ++i) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            a[i] = 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
        }
        return a;
    }

    private static void runDgels(double[] a, double[] b, int m, int n) {
        PrintStream saved = System.err;
        System.setErr(new PrintStream(new ByteArrayOutputStream()));
        try {
            int ldb = Math.max(1, Math.max(m, n));
            double[] work = new double[m * 64 + 4096];
            Dgels.dgels(Trans.NO_TRANS, m, n, 1, a, 0, Math.max(1, m), b, 0, ldb, work, 0, work.length);
        } finally {
            System.setErr(saved);
        }
    }

    @Test(timeout = TIMEOUT)
    public void testDgelsTerminatesOnAnInfiniteDesignEntry() {
        for (double bad : new double[] { Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY }) {
            int m = 12;
            int n = 4;
            double[] a = design(m, n);
            a[7 + 2 * m] = bad;
            double[] b = new double[m];
            for (int i = 0; i < m; ++i) {
                b[i] = 1.0;
            }
            try {
                runDgels(a, b, m, n);
                // returning at all is the assertion; a wrong answer is a far
                // smaller problem than a call that never comes back
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage().contains("DLASCL"));
            }
        }
    }

    @Test(timeout = TIMEOUT)
    public void testDgelsTerminatesOnAnInfiniteRightHandSide() {
        int m = 12;
        int n = 4;
        double[] a = design(m, n);
        double[] b = new double[Math.max(m, n)];
        for (int i = 0; i < m; ++i) {
            b[i] = 1.0;
        }
        b[5] = Double.POSITIVE_INFINITY;
        try {
            runDgels(a, b, m, n);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("DLASCL"));
        }
    }

    @Test(timeout = TIMEOUT)
    public void testDgelsIsUnchangedOnExtremeButFiniteScales() {
        // the scaling path is the reason Dlascl is called at all, so this is
        // the regression that matters: the solution has to stay invariant when
        // A and b are scaled together
        int m = 40;
        int n = 6;
        double[] a0 = design(m, n);
        double[] b0 = new double[m];
        for (int i = 0; i < m; ++i) {
            b0[i] = 0.5 + 0.01 * i;
        }
        double[] reference = solveWith(a0, b0, m, n, 1.0);
        for (double s : new double[] { 1.0e-300, 1.0e-150, 1.0e150, 1.0e300 }) {
            double[] x = solveWith(a0, b0, m, n, s);
            for (int j = 0; j < n; ++j) {
                assertEquals("scale " + s + ", coefficient " + j, reference[j], x[j],
                        1.0e-12 * Math.abs(reference[j]));
            }
        }
    }

    private static double[] solveWith(double[] a0, double[] b0, int m, int n, double scale) {
        int ldb = Math.max(1, Math.max(m, n));
        double[] a = new double[a0.length];
        for (int i = 0; i < a0.length; ++i) {
            a[i] = scale * a0[i];
        }
        double[] b = new double[ldb];
        for (int i = 0; i < m; ++i) {
            b[i] = scale * b0[i];
        }
        double[] work = new double[m * 64 + 4096];
        Dgels.dgels(Trans.NO_TRANS, m, n, 1, a, 0, m, b, 0, ldb, work, 0, work.length);
        double[] x = new double[n];
        System.arraycopy(b, 0, x, 0, n);
        return x;
    }
}
