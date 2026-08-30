package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.gemm.Trans;

/**
 * The QR and LQ chain that {@code Dgels} runs on, blocked and unblocked.
 * <p>
 * Until this class existed no test anywhere in the repository drove a matrix
 * large enough to take a blocked path: {@code Dgetrf} blocks above order 64 and
 * {@code Dgeqrf} above {@code min(m, n) = 128}, and the largest matrix any test
 * handed to this package was 50 by 7. {@code Dlarft} and {@code Dlarfb}, 450
 * lines between them, had therefore never been executed at all.
 * <p>
 * There is no reference LAPACK here to compare against, so every assertion is
 * on a structural invariant -- orthonormality, reconstruction, the defining
 * equation of a least squares solution -- or on the blocked path agreeing with
 * the unblocked one on the same input.
 */
public class DgelsTest {

    /** the tolerances the invariants actually meet, measured over these shapes */
    private static final double ORTHONORMAL = 1.0e-14;
    private static final double RECONSTRUCT = 1.0e-13;
    private static final double DEFINING = 1.0e-12;
    /** blocked and unblocked sum in a different order, so they agree closely but not exactly */
    private static final double PATHS_AGREE = 1.0e-11;

    // =========================================================================
    // the coverage assumption itself
    // =========================================================================

    /**
     * The sizes every other test here picks are only meaningful while these
     * thresholds hold. If the tuning in {@code Ilaenv} moves, this fails first
     * and says why, rather than the blocked path quietly going untested again.
     */
    @Test
    public void testTheSizesBelowStillReachTheBlockedPaths() {
        assertEquals("DGEQRF block size", 32, Ilaenv.ilaenv(1, "DGEQRF", " ", 200, 200, -1, -1));
        assertEquals("DGEQRF crossover", 128, Ilaenv.ilaenv(3, "DGEQRF", " ", 200, 200, -1, -1));
        assertEquals("DGETRF block size", 64, Ilaenv.ilaenv(1, "DGETRF", " ", 200, 200, -1, -1));
        assertTrue("160 must be past the QR crossover", 160 > 128);
        assertTrue("100 must be past the LU block size", 100 > 64);
    }

    // =========================================================================
    // QR and LQ: orthonormality and reconstruction
    // =========================================================================

    @Test
    public void testTheQrFactorsReconstructTheMatrix() {
        qr(20, 20, false);
        qr(60, 25, false);
        qr(160, 160, true);
        qr(200, 90, false);
    }

    @Test
    public void testTheLqFactorsReconstructTheMatrix() {
        lq(20, 20, false);
        lq(25, 60, false);
        lq(160, 160, true);
    }

    /** Factor with Dgeqrf, rebuild Q by applying the reflectors to the identity. */
    private static void qr(int m, int n, boolean expectBlocked) {
        String what = m + " by " + n;
        assertEquals(what + " should be " + (expectBlocked ? "blocked" : "unblocked"), Boolean.valueOf(expectBlocked),
                Boolean.valueOf(isQrBlocked(m, n)));

        int k = Math.min(m, n);
        double[] original = matrix(m, n, 20260830L + m);
        double[] a = original.clone();
        double[] tau = new double[k];
        factorQr(m, n, a, tau);

        double[] q = identity(m);
        applyQ(m, k, a, tau, q);

        assertTrue(what + " Q is not orthonormal", orthonormality(q, m) <= ORTHONORMAL);
        double[] qr = multiply(false, false, m, n, m, q, m, upper(m, n, a), m);
        assertTrue(what + " QR does not reconstruct A",
                maxDifference(qr, original) <= RECONSTRUCT * maxMagnitude(original));
    }

    /** Factor with Dgelqf, rebuild Q by applying the reflectors to the identity. */
    private static void lq(int m, int n, boolean expectBlocked) {
        String what = m + " by " + n;
        int k = Math.min(m, n);
        int nb = Ilaenv.ilaenv(1, "DGELQF", " ", m, n, -1, -1);
        int nx = Math.max(0, Ilaenv.ilaenv(3, "DGELQF", " ", m, n, -1, -1));
        assertEquals(what + " should be " + (expectBlocked ? "blocked" : "unblocked"), Boolean.valueOf(expectBlocked),
                Boolean.valueOf(nb > 1 && nb < k && nx < k));

        double[] original = matrix(m, n, 20260901L + n);
        double[] a = original.clone();
        double[] tau = new double[k];
        intW info = new intW(0);
        double[] query = new double[1];
        Dgelqf.dgelqf(m, n, a, 0, m, tau, 0, query, 0, -1, info);
        int lwork = Math.max((int) query[0], Math.max(1, m));
        Dgelqf.dgelqf(m, n, a, 0, m, tau, 0, new double[lwork], 0, lwork, info);

        double[] q = identity(n);
        double[] query2 = new double[1];
        Dormlq.dormlq("Right", "No transpose", n, n, k, a, 0, m, tau, 0, q, 0, n, query2, 0, -1, new intW(0));
        int lw2 = Math.max((int) query2[0], Math.max(1, n));
        Dormlq.dormlq("Right", "No transpose", n, n, k, a, 0, m, tau, 0, q, 0, n, new double[lw2], 0, lw2,
                new intW(0));

        assertTrue(what + " Q is not orthonormal", orthonormality(q, n) <= ORTHONORMAL);
        double[] lq = multiply(false, false, m, n, n, lower(m, n, a), m, q, n);
        assertTrue(what + " LQ does not reconstruct A",
                maxDifference(lq, original) <= RECONSTRUCT * maxMagnitude(original));
    }

    // =========================================================================
    // the blocked path against the unblocked one, on the same input
    // =========================================================================

    @Test
    public void testTheBlockedQrAgreesWithTheUnblockedOne() {
        int m = 160;
        int n = 160;
        assertTrue("160 by 160 must take the blocked path", isQrBlocked(m, n));

        double[] original = matrix(m, n, 4242L);
        double[] blocked = original.clone();
        double[] tauBlocked = new double[n];
        factorQr(m, n, blocked, tauBlocked);

        double[] unblocked = original.clone();
        double[] tauUnblocked = new double[n];
        Dgeqr2.dgeqr2(m, n, unblocked, 0, m, tauUnblocked, 0, new double[n], 0, new intW(0));

        assertTrue("the two QR paths disagree on the factor",
                maxDifference(blocked, unblocked) <= PATHS_AGREE);
        assertTrue("the two QR paths disagree on tau", maxDifference(tauBlocked, tauUnblocked) <= PATHS_AGREE);
    }

    @Test
    public void testTheBlockedLuAgreesWithTheUnblockedOne() {
        int n = 100;
        assertTrue("100 must take the blocked LU path", Ilaenv.ilaenv(1, "DGETRF", " ", n, n, -1, -1) < n);

        double[] original = matrix(n, n, 777L);
        double[] blocked = original.clone();
        int[] pivotsBlocked = new int[n];
        Dgetrf.dgetrf(n, n, blocked, 0, n, pivotsBlocked, 0, new intW(0));

        double[] unblocked = original.clone();
        int[] pivotsUnblocked = new int[n];
        Dgetf2.dgetf2(n, n, unblocked, 0, n, pivotsUnblocked, 0, new intW(0));

        for (int i = 0; i < n; ++i) {
            assertEquals("pivot " + i, pivotsUnblocked[i], pivotsBlocked[i]);
        }
        assertTrue("the two LU paths disagree on the factor",
                maxDifference(blocked, unblocked) <= PATHS_AGREE);
    }

    // =========================================================================
    // Dgels itself: each of its four cases against its own defining equation
    // =========================================================================

    @Test
    public void testTheOverdeterminedCasesSolveTheLeastSquaresProblem() {
        gels(Trans.NO_TRANS, 40, 15, 3);
        gels(Trans.TRANS, 15, 40, 3);
        gels(Trans.NO_TRANS, 300, 150, 2);
    }

    @Test
    public void testTheUnderdeterminedCasesSolveTheSystemExactly() {
        gels(Trans.NO_TRANS, 15, 40, 3);
        gels(Trans.TRANS, 40, 15, 3);
    }

    /**
     * For an overdetermined problem the residual must be orthogonal to the
     * column space, which is the defining property of the least squares
     * solution; for an underdetermined one the system must hold exactly.
     */
    private static void gels(Trans trans, int m, int n, int nrhs) {
        boolean transposed = (trans != Trans.NO_TRANS);
        int rows = transposed ? n : m;
        int columns = transposed ? m : n;
        int ldb = Math.max(1, Math.max(m, n));
        String what = (transposed ? "TRANS " : "NO_TRANS ") + m + " by " + n;

        double[] a = matrix(m, n, 31337L + m + n);
        double[] b = new double[ldb * nrhs];
        long state = 909090L + nrhs;
        for (int j = 0; j < nrhs; ++j) {
            for (int i = 0; i < rows; ++i) {
                state = state * 6364136223846793005L + 1442695040888963407L;
                b[j * ldb + i] = 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
            }
        }

        double[] factor = a.clone();
        double[] solution = b.clone();
        int lwork = queryWorkspace(trans, m, n, nrhs, ldb);
        assertTrue(what + " did not succeed", Dgels.dgels(trans, m, n, nrhs, factor, 0, Math.max(1, m), solution, 0,
                ldb, new double[lwork], 0, lwork));

        double[] x = new double[columns * nrhs];
        for (int j = 0; j < nrhs; ++j) {
            for (int i = 0; i < columns; ++i) {
                x[j * columns + i] = solution[j * ldb + i];
            }
        }
        double[] product = multiply(transposed, false, rows, nrhs, columns, a, Math.max(1, m), x, columns);
        double[] residual = new double[rows * nrhs];
        for (int j = 0; j < nrhs; ++j) {
            for (int i = 0; i < rows; ++i) {
                residual[j * rows + i] = product[j * rows + i] - b[j * ldb + i];
            }
        }

        if (rows >= columns) {
            double[] normal = multiply(!transposed, false, columns, nrhs, rows, a, Math.max(1, m), residual, rows);
            double scaled = maxMagnitude(normal) / (maxMagnitude(a) * maxMagnitude(b));
            assertTrue(what + " residual is not orthogonal to the column space: " + scaled, scaled <= DEFINING);
        } else {
            double scaled = maxMagnitude(residual) / maxMagnitude(b);
            assertTrue(what + " does not solve the system exactly: " + scaled, scaled <= DEFINING);
        }
    }

    /**
     * The scaling branches fire only for a matrix whose largest entry leaves
     * {@code [1.0e-292, 9.98e291]}, which nothing else here reaches. Scaling A
     * by a constant must scale the solution by its reciprocal, and scaling b
     * must scale the solution by the same factor.
     */
    @Test
    public void testTheScalingPathIsUndoneAgain() {
        scaled(1.0e-300, 1.0, "A driven down to 1e-300");
        scaled(1.0e300, 1.0, "A driven up to 1e300");
        scaled(1.0, 1.0e-300, "b driven down to 1e-300");
        scaled(1.0, 1.0e300, "b driven up to 1e300");
    }

    private static void scaled(double scaleA, double scaleB, String what) {
        int m = 30;
        int n = 12;
        int ldb = Math.max(1, Math.max(m, n));
        double[] a = matrix(m, n, 5150L);
        double[] b = new double[ldb];
        for (int i = 0; i < m; ++i) {
            b[i] = 1.0 + i;
        }

        double[] plain = solveScaled(a, b, m, n, ldb, 1.0, 1.0);
        double[] moved = solveScaled(a, b, m, n, ldb, scaleA, scaleB);
        double expected = scaleB / scaleA;
        for (int i = 0; i < n; ++i) {
            double want = plain[i] * expected;
            assertEquals(what + ", entry " + i, want, moved[i], 1.0e-12 * Math.abs(want));
        }
    }

    private static double[] solveScaled(double[] a0, double[] b0, int m, int n, int ldb, double scaleA,
            double scaleB) {
        double[] a = new double[a0.length];
        for (int i = 0; i < a0.length; ++i) {
            a[i] = a0[i] * scaleA;
        }
        double[] b = new double[b0.length];
        for (int i = 0; i < b0.length; ++i) {
            b[i] = b0[i] * scaleB;
        }
        int lwork = queryWorkspace(Trans.NO_TRANS, m, n, 1, ldb);
        assertTrue(Dgels.dgels(Trans.NO_TRANS, m, n, 1, a, 0, m, b, 0, ldb, new double[lwork], 0, lwork));
        double[] x = new double[n];
        System.arraycopy(b, 0, x, 0, n);
        return x;
    }

    // =========================================================================
    // the boundary
    // =========================================================================

    /**
     * The workspace query reads nothing but {@code work}, so empty arrays for
     * {@code a} and {@code b} must stay legal. {@code LinearEquationsSolver}
     * passes exactly that, and it is why the length checks cannot run when
     * {@code lwork} is {@code -1}.
     */
    @Test
    public void testAWorkspaceQueryReadsNeitherMatrix() {
        double[] work = new double[1];
        assertTrue(Dgels.dgels(Trans.NO_TRANS, 40, 15, 3, new double[0], 0, 40, new double[0], 0, 40, work, 0, -1));
        assertTrue("the query must report a usable workspace size", work[0] >= 15.0);
    }

    @Test
    public void testTheArrayLengthsAreChecked() {
        int m = 5;
        int n = 3;
        int ldb = Math.max(1, Math.max(m, n));
        int lwork = queryWorkspace(Trans.NO_TRANS, m, n, 1, ldb);
        double[] a = matrix(m, n, 6L);
        double[] b = new double[ldb];

        refuses("a", () -> Dgels.dgels(Trans.NO_TRANS, m, n, 1, new double[m * n - 1], 0, m, b.clone(), 0, ldb,
                new double[lwork], 0, lwork));
        refuses("b", () -> Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), 0, m, new double[ldb - 1], 0, ldb,
                new double[lwork], 0, lwork));
        refuses("work", () -> Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), 0, m, b.clone(), 0, ldb,
                new double[lwork - 1], 0, lwork));
        refuses("_a_offset", () -> Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), -1, m, b.clone(), 0, ldb,
                new double[lwork], 0, lwork));

        // and the exactly sufficient arrays are accepted: the bound is
        // lda * (n - 1) + m, not lda * n, so a packed matrix is not rejected
        assertTrue("an exactly packed argument must be accepted", Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), 0,
                m, b.clone(), 0, ldb, new double[lwork], 0, lwork));
    }

    @Test
    public void testANonFiniteEntryIsReportedRatherThanPropagated() {
        int m = 8;
        int n = 4;
        int ldb = Math.max(1, Math.max(m, n));
        int lwork = queryWorkspace(Trans.NO_TRANS, m, n, 1, ldb);
        double[] a = matrix(m, n, 8080L);
        double[] b = new double[ldb];
        for (int i = 0; i < m; ++i) {
            b[i] = 1.0;
        }
        assertTrue("the clean problem must succeed", Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), 0, m, b.clone(),
                0, ldb, new double[lwork], 0, lwork));

        double[] bad = { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
        for (int t = 0; t < bad.length; ++t) {
            double[] spoiledA = a.clone();
            spoiledA[m + 1] = bad[t];
            assertFalse("a design matrix holding " + bad[t] + " must be reported, not propagated",
                    Dgels.dgels(Trans.NO_TRANS, m, n, 1, spoiledA, 0, m, b.clone(), 0, ldb, new double[lwork], 0,
                            lwork));

            double[] spoiledB = b.clone();
            spoiledB[2] = bad[t];
            assertFalse("a right hand side holding " + bad[t] + " must be reported, not propagated",
                    Dgels.dgels(Trans.NO_TRANS, m, n, 1, a.clone(), 0, m, spoiledB, 0, ldb, new double[lwork], 0,
                            lwork));
        }
    }

    // =========================================================================
    // helpers
    // =========================================================================

    private static boolean isQrBlocked(int m, int n) {
        int k = Math.min(m, n);
        int nb = Ilaenv.ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
        int nx = Math.max(0, Ilaenv.ilaenv(3, "DGEQRF", " ", m, n, -1, -1));
        return nb > 1 && nb < k && nx < k;
    }

    private static void factorQr(int m, int n, double[] a, double[] tau) {
        intW info = new intW(0);
        double[] query = new double[1];
        Dgeqrf.dgeqrf(m, n, a, 0, m, tau, 0, query, 0, -1, info);
        int lwork = Math.max((int) query[0], Math.max(1, n));
        Dgeqrf.dgeqrf(m, n, a, 0, m, tau, 0, new double[lwork], 0, lwork, info);
        assertEquals("Dgeqrf reported " + info.val, 0, info.val);
    }

    /** overwrites c, which is m by m, with Q * c */
    private static void applyQ(int m, int k, double[] a, double[] tau, double[] c) {
        double[] query = new double[1];
        Dormqr.dormqr("Left", Trans.NO_TRANS, m, m, k, a, 0, m, tau, 0, c, 0, m, query, 0, -1, new intW(0));
        int lwork = Math.max((int) query[0], Math.max(1, m));
        Dormqr.dormqr("Left", Trans.NO_TRANS, m, m, k, a, 0, m, tau, 0, c, 0, m, new double[lwork], 0, lwork,
                new intW(0));
    }

    private static int queryWorkspace(Trans trans, int m, int n, int nrhs, int ldb) {
        double[] query = new double[1];
        Dgels.dgels(trans, m, n, nrhs, new double[0], 0, Math.max(1, m), new double[0], 0, ldb, query, 0, -1);
        return Math.max((int) query[0], 1);
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

    /** an m by n matrix with entries in [-1, 1], column-major, from a seeded LCG */
    private static double[] matrix(int m, int n, long seed) {
        double[] a = new double[m * n];
        long state = seed;
        for (int i = 0; i < a.length; ++i) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            a[i] = 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
        }
        return a;
    }

    private static double[] identity(int n) {
        double[] a = new double[n * n];
        for (int i = 0; i < n; ++i) {
            a[i * n + i] = 1.0;
        }
        return a;
    }

    /** max |Q'Q - I| over an n by n matrix */
    private static double orthonormality(double[] q, int n) {
        double[] product = multiply(true, false, n, n, n, q, n, q, n);
        double worst = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double expected = (i == j) ? 1.0 : 0.0;
                worst = Math.max(worst, Math.abs(product[j * n + i] - expected));
            }
        }
        return worst;
    }

    /** the upper triangle of an m by n factor, everything below the diagonal zeroed */
    private static double[] upper(int m, int n, double[] a) {
        double[] r = new double[m * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i <= Math.min(j, m - 1); ++i) {
                r[j * m + i] = a[j * m + i];
            }
        }
        return r;
    }

    /** the lower triangle of an m by n factor, everything above the diagonal zeroed */
    private static double[] lower(int m, int n, double[] a) {
        double[] r = new double[m * n];
        for (int j = 0; j < n; ++j) {
            for (int i = j; i < m; ++i) {
                r[j * m + i] = a[j * m + i];
            }
        }
        return r;
    }

    /** a plain triple loop, column-major, deliberately independent of math.gemm */
    private static double[] multiply(boolean transA, boolean transB, int m, int n, int k, double[] a, int lda,
            double[] b, int ldb) {
        double[] c = new double[m * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int p = 0; p < k; ++p) {
                    double left = transA ? a[i * lda + p] : a[p * lda + i];
                    double right = transB ? b[p * ldb + j] : b[j * ldb + p];
                    sum += left * right;
                }
                c[j * m + i] = sum;
            }
        }
        return c;
    }

    private static double maxMagnitude(double[] a) {
        double worst = 0.0;
        for (int i = 0; i < a.length; ++i) {
            worst = Math.max(worst, Math.abs(a[i]));
        }
        return worst;
    }

    private static double maxDifference(double[] x, double[] y) {
        double worst = 0.0;
        for (int i = 0; i < x.length; ++i) {
            worst = Math.max(worst, Math.abs(x[i] - y[i]));
        }
        return worst;
    }
}
