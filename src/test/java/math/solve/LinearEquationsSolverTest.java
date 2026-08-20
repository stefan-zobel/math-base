package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.DMatrix;

/**
 * Tests for {@link LinearEquationsSolver}.
 * <p>
 * The assertions are on structural invariants rather than on tabulated values:
 * the residual of the defining equation for a square system, the normal
 * equation residual for a least squares problem, and the rank guarantee the
 * class documents.
 */
public class LinearEquationsSolverTest {

    private static final double TIGHT = 1.0e-12;

    // =========================================================================
    // helpers
    // =========================================================================

    /** Deterministic in-test LCG, per the test conventions. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            this.state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }

        /** Uniform on (-1, 1). */
        double symmetric() {
            return 2.0 * next() - 1.0;
        }
    }

    private static DMatrix random(int rows, int cols, Lcg rng) {
        DMatrix m = new DMatrix(rows, cols);
        for (int j = 0; j < cols; ++j) {
            for (int i = 0; i < rows; ++i) {
                m.setUnsafe(i, j, rng.symmetric());
            }
        }
        return m;
    }

    /**
     * The Hilbert matrix, the textbook near singular family. Its condition
     * number grows roughly like {@code e^(3.5*n)}, and it needs no seed.
     */
    private static DMatrix hilbert(int n) {
        DMatrix m = new DMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                m.setUnsafe(i, j, 1.0 / (i + j + 1.0));
            }
        }
        return m;
    }

    /** The Vandermonde matrix on an equispaced grid over [0, 1]. */
    private static DMatrix vandermonde(int rows, int cols) {
        DMatrix m = new DMatrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            double x = i / (double) (rows - 1);
            double p = 1.0;
            for (int j = 0; j < cols; ++j) {
                m.setUnsafe(i, j, p);
                p *= x;
            }
        }
        return m;
    }

    private static DMatrix solve(DMatrix A, DMatrix B) {
        return LinearEquationsSolver.solve(A, B, new DMatrix(A.numColumns(), B.numColumns()));
    }

    /** {@code ||A*X - B||_inf}, the residual of the system itself. */
    private static double systemResidual(DMatrix A, DMatrix X, DMatrix B) {
        double worst = 0.0;
        for (int j = 0; j < B.numColumns(); ++j) {
            for (int i = 0; i < A.numRows(); ++i) {
                double s = 0.0;
                for (int k = 0; k < A.numColumns(); ++k) {
                    s += A.getUnsafe(i, k) * X.getUnsafe(k, j);
                }
                worst = Math.max(worst, Math.abs(s - B.getUnsafe(i, j)));
            }
        }
        return worst;
    }

    /**
     * {@code ||A'(A*X - B)||_inf}. This vanishes exactly at the least squares
     * solution, so it is the defining property rather than a golden value.
     */
    private static double normalEquationResidual(DMatrix A, DMatrix X, DMatrix B) {
        int m = A.numRows();
        int n = A.numColumns();
        double worst = 0.0;
        for (int c = 0; c < B.numColumns(); ++c) {
            double[] r = new double[m];
            for (int i = 0; i < m; ++i) {
                double s = 0.0;
                for (int k = 0; k < n; ++k) {
                    s += A.getUnsafe(i, k) * X.getUnsafe(k, c);
                }
                r[i] = s - B.getUnsafe(i, c);
            }
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int i = 0; i < m; ++i) {
                    s += A.getUnsafe(i, j) * r[i];
                }
                worst = Math.max(worst, Math.abs(s));
            }
        }
        return worst;
    }

    /** The largest absolute entry, used to scale the tolerances. */
    private static double scaleOf(DMatrix m) {
        double s = 0.0;
        for (int j = 0; j < m.numColumns(); ++j) {
            for (int i = 0; i < m.numRows(); ++i) {
                s = Math.max(s, Math.abs(m.getUnsafe(i, j)));
            }
        }
        return s;
    }

    // =========================================================================
    // the square path
    // =========================================================================

    @Test
    public void testASquareSystemIsSolvedToMachinePrecision() {
        Lcg rng = new Lcg(11223344L);
        for (int n : new int[] { 1, 2, 5, 17, 60 }) {
            DMatrix A = random(n, n, rng);
            DMatrix B = random(n, 1, rng);
            DMatrix X = solve(A, B);
            assertEquals(n, X.numRows());
            assertEquals(1, X.numColumns());
            assertTrue("n = " + n + ", residual " + systemResidual(A, X, B),
                    systemResidual(A, X, B) <= TIGHT * n * scaleOf(A) * scaleOf(X));
        }
    }

    @Test
    public void testASquareSystemWithSeveralRightHandSides() {
        Lcg rng = new Lcg(55667788L);
        int n = 24;
        DMatrix A = random(n, n, rng);
        DMatrix B = random(n, 5, rng);
        DMatrix X = solve(A, B);
        assertEquals(5, X.numColumns());
        assertTrue(systemResidual(A, X, B) <= TIGHT * n * scaleOf(A) * scaleOf(X));
    }

    @Test
    public void testTheIdentityReproducesTheRightHandSide() {
        int n = 9;
        DMatrix B = random(n, 3, new Lcg(99001122L));
        DMatrix X = solve(DMatrix.identity(n), B);
        // I * X == B has to hold bit for bit, no arithmetic is involved
        assertEquals(0.0, systemResidual(DMatrix.identity(n), X, B), 0.0);
    }

    @Test
    public void testAnExactlySingularSquareMatrixIsReported() {
        // A duplicated column sometimes, but only sometimes, produces an
        // exact zero pivot in Dgetf2, which is the single case the LU path
        // was ever able to reject. Over 400 seeds of this construction it
        // fires 100 times; this seed is one of them, and the case that it
        // misses is testAnExactlyDuplicatedColumnInASquareMatrixIsReported.
        Lcg rng = new Lcg(2027808446L);
        int n = 6;
        DMatrix A = random(n, n, rng);
        for (int i = 0; i < n; ++i) {
            A.setUnsafe(i, n - 1, A.getUnsafe(i, 0));
        }
        assertThrowsRuntime(A, random(n, 1, rng));
    }

    // =========================================================================
    // the least squares path, m > n
    // =========================================================================

    @Test
    public void testAnOverdeterminedSystemMeetsTheNormalEquations() {
        Lcg rng = new Lcg(24682468L);
        // the sizes straddle the nb = 32 blocking threshold of the LAPACK QR,
        // so both the blocked and the unblocked path stay covered
        int[][] sizes = { { 5, 3 }, { 40, 10 }, { 200, 33 }, { 200, 64 }, { 400, 120 } };
        for (int[] size : sizes) {
            int m = size[0];
            int n = size[1];
            DMatrix A = random(m, n, rng);
            DMatrix B = random(m, 1, rng);
            DMatrix X = solve(A, B);
            assertEquals(n, X.numRows());
            double residual = normalEquationResidual(A, X, B);
            assertTrue(m + "x" + n + ", normal equation residual " + residual,
                    residual <= TIGHT * m * scaleOf(A) * scaleOf(A) * scaleOf(X));
        }
    }

    @Test
    public void testAConsistentOverdeterminedSystemIsSolvedExactly() {
        // b lies in the column space, so the least squares solution is the
        // exact one and the residual of the system itself has to vanish
        Lcg rng = new Lcg(36903690L);
        int m = 50;
        int n = 7;
        DMatrix A = random(m, n, rng);
        DMatrix xTrue = random(n, 1, rng);
        DMatrix B = new DMatrix(m, 1);
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int k = 0; k < n; ++k) {
                s += A.getUnsafe(i, k) * xTrue.getUnsafe(k, 0);
            }
            B.setUnsafe(i, 0, s);
        }
        DMatrix X = solve(A, B);
        for (int k = 0; k < n; ++k) {
            assertEquals(xTrue.getUnsafe(k, 0), X.getUnsafe(k, 0), 1.0e-11);
        }
    }

    @Test
    public void testAnOverdeterminedSystemWithSeveralRightHandSides() {
        Lcg rng = new Lcg(48154815L);
        int m = 80;
        int n = 12;
        DMatrix A = random(m, n, rng);
        DMatrix B = random(m, 4, rng);
        DMatrix X = solve(A, B);
        assertEquals(4, X.numColumns());
        assertTrue(normalEquationResidual(A, X, B) <= TIGHT * m * scaleOf(A) * scaleOf(A) * scaleOf(X));
    }

    // =========================================================================
    // the minimum norm path, m < n
    // =========================================================================

    @Test
    public void testAnUnderdeterminedSystemSatisfiesTheEquationsExactly() {
        Lcg rng = new Lcg(51625162L);
        int[][] sizes = { { 3, 5 }, { 10, 40 }, { 33, 200 } };
        for (int[] size : sizes) {
            int m = size[0];
            int n = size[1];
            DMatrix A = random(m, n, rng);
            DMatrix B = random(m, 1, rng);
            DMatrix X = solve(A, B);
            assertEquals(n, X.numRows());
            // an underdetermined consistent system is solved exactly, the
            // remaining freedom is spent on minimizing ||x||
            assertTrue(m + "x" + n + ", residual " + systemResidual(A, X, B),
                    systemResidual(A, X, B) <= TIGHT * n * scaleOf(A) * scaleOf(X));
        }
    }

    // =========================================================================
    // dimension checking, unchanged behaviour
    // =========================================================================

    @Test(expected = IndexOutOfBoundsException.class)
    public void testMismatchedRowCountIsRejected() {
        LinearEquationsSolver.solve(new DMatrix(4, 4), new DMatrix(3, 1), new DMatrix(4, 1));
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testMismatchedSolutionRowCountIsRejected() {
        LinearEquationsSolver.solve(new DMatrix(4, 3), new DMatrix(4, 1), new DMatrix(4, 1));
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testMismatchedSolutionColumnCountIsRejected() {
        LinearEquationsSolver.solve(new DMatrix(4, 4), new DMatrix(4, 2), new DMatrix(4, 1));
    }

    // =========================================================================
    // the one internal caller
    // =========================================================================

    @Test
    public void testMatrixInverseOfAWellConditionedMatrixIsUnaffected() {
        Lcg rng = new Lcg(60066006L);
        int n = 12;
        DMatrix A = random(n, n, rng);
        DMatrix inv = A.inverse();
        DMatrix product = A.mul(inv);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                assertEquals((i == j) ? 1.0 : 0.0, product.getUnsafe(i, j), 1.0e-10);
            }
        }
    }

    // =========================================================================
    // shared assertion
    // =========================================================================

    private static void assertThrowsRuntime(DMatrix A, DMatrix B) {
        try {
            solve(A, B);
            fail("expected a RuntimeException for a rank deficient " + A.numRows() + "x" + A.numColumns() + " matrix");
        } catch (RuntimeException expected) {
            assertTrue("the message should name the problem, was: " + expected.getMessage(),
                    expected.getMessage() != null && expected.getMessage().length() > 0);
        }
    }

    private static DMatrix ones(int rows) {
        DMatrix b = new DMatrix(rows, 1);
        for (int i = 0; i < rows; ++i) {
            b.setUnsafe(i, 0, 1.0);
        }
        return b;
    }

    // =========================================================================
    // the rank guard, the point of this class
    // =========================================================================

    @Test
    public void testAnExactlyDuplicatedColumnIsReportedOnTheLeastSquaresPath() {
        // the case that motivated the guard: without it Dgels returns success
        // and a solution of norm 5.5e14, because the Householder QR leaves
        // R[3][3] = 2.2e-16 rather than the exact zero Dtrtrs looks for
        Lcg rng = new Lcg(70017001L);
        int m = 20;
        int n = 4;
        DMatrix A = random(m, n, rng);
        for (int i = 0; i < m; ++i) {
            A.setUnsafe(i, n - 1, A.getUnsafe(i, 0));
        }
        assertThrowsRuntime(A, random(m, 1, rng));
    }

    @Test
    public void testALinearlyDependentColumnIsReported() {
        // subtler than an outright duplicate: one column is a combination of
        // two others, so no two columns are equal
        Lcg rng = new Lcg(80028002L);
        int m = 30;
        int n = 5;
        DMatrix A = random(m, n, rng);
        for (int i = 0; i < m; ++i) {
            A.setUnsafe(i, n - 1, 0.5 * A.getUnsafe(i, 0) - 2.0 * A.getUnsafe(i, 1));
        }
        assertThrowsRuntime(A, random(m, 1, rng));
    }

    @Test
    public void testRankDeficiencyIsReportedForTheUnderdeterminedCaseToo() {
        // m < n takes the LQ branch of Dgels, so the diagonal being tested is
        // that of L rather than of R
        Lcg rng = new Lcg(90039003L);
        int m = 5;
        int n = 30;
        DMatrix A = random(m, n, rng);
        for (int j = 0; j < n; ++j) {
            A.setUnsafe(m - 1, j, A.getUnsafe(0, j));
        }
        assertThrowsRuntime(A, random(m, 1, rng));
    }

    @Test
    public void testAnExactlyDuplicatedColumnInASquareMatrixIsReported() {
        // the counterpart on the LU path, and the seed is chosen so that
        // Dgesv returns success: the exact zero pivot that Dgetf2 looks for
        // does not occur, so only the new guard can catch this. Over 400
        // seeds of this construction the old code was silently wrong 300
        // times
        Lcg rng = new Lcg(13571357L);
        int n = 6;
        DMatrix A = random(n, n, rng);
        for (int i = 0; i < n; ++i) {
            A.setUnsafe(i, n - 1, A.getUnsafe(i, 0));
        }
        assertThrowsRuntime(A, random(n, 1, rng));
    }

    @Test
    public void testANumericallySingularSquareMatrixIsReported() {
        // Hilbert of order 16 has a condition number far beyond what double
        // precision carries, yet produces no exact zero pivot: Dgesv returns
        // success. min|d| / max|d| is 7.9e-17 against a threshold of 1.8e-15
        assertThrowsRuntime(hilbert(16), ones(16));
    }

    @Test
    public void testEveryDeficientPathThrowsTheSameWay() {
        // the invariant this whole guard exists for: square or not, exactly
        // singular or only numerically so, the caller sees one behaviour
        Lcg rng = new Lcg(11111111L);
        DMatrix square = random(7, 7, rng);
        for (int i = 0; i < 7; ++i) {
            square.setUnsafe(i, 6, square.getUnsafe(i, 1));
        }
        DMatrix tall = random(25, 6, rng);
        for (int i = 0; i < 25; ++i) {
            tall.setUnsafe(i, 5, tall.getUnsafe(i, 2));
        }
        DMatrix wide = random(6, 25, rng);
        for (int j = 0; j < 25; ++j) {
            wide.setUnsafe(5, j, wide.getUnsafe(2, j));
        }
        assertThrowsRuntime(square, ones(7));
        assertThrowsRuntime(tall, ones(25));
        assertThrowsRuntime(wide, ones(6));
        assertThrowsRuntime(hilbert(16), ones(16));
    }

    // =========================================================================
    // the false positive guard, the most important tests in this file
    // =========================================================================

    @Test
    public void testAnIllConditionedButFullRankSquareSystemStillSolves() {
        // Hilbert of order 10 has a condition number near 1e13 and is still
        // worth solving; its ratio is 2.6e-12 against a threshold of 1.1e-15,
        // a margin of more than three orders of magnitude
        DMatrix h = hilbert(10);
        DMatrix b = ones(10);
        DMatrix x = solve(h, b);
        assertEquals(10, x.numRows());
        // the answer is poor, as it must be at this conditioning, but the
        // system was solved rather than rejected
        assertTrue(systemResidual(h, x, b) <= 1.0e-9);
    }

    @Test
    public void testAnIllConditionedButFullRankLeastSquaresProblemStillSolves() {
        // a degree 12 Vandermonde fit on [0, 1] is badly conditioned by any
        // standard, and its ratio of 8.3e-8 sits seven orders of magnitude
        // above the threshold
        DMatrix v = vandermonde(60, 13);
        DMatrix b = ones(60);
        DMatrix x = solve(v, b);
        assertEquals(13, x.numRows());
        assertTrue(normalEquationResidual(v, x, b) <= 1.0e-10);
    }

    @Test
    public void testTheGuardLeavesTheWellConditionedCasesAlone() {
        // every shape, well conditioned, must pass without complaint
        Lcg rng = new Lcg(22222222L);
        for (int[] size : new int[][] { { 8, 8 }, { 50, 8 }, { 8, 50 }, { 200, 64 } }) {
            DMatrix A = random(size[0], size[1], rng);
            solve(A, ones(size[0]));
        }
    }

    // =========================================================================
    // scale invariance and the degenerate inputs
    // =========================================================================

    @Test
    public void testDetectionIsScaleInvariant() {
        // Dgels rescales A through Dlascl when its norm leaves the safe range,
        // so the diagonal being tested is the scaled one. A ratio of two
        // diagonal entries is invariant under that, and this pins it
        Lcg rng = new Lcg(33333333L);
        int m = 20;
        int n = 4;
        DMatrix base = random(m, n, rng);
        for (int i = 0; i < m; ++i) {
            base.setUnsafe(i, n - 1, base.getUnsafe(i, 0));
        }
        DMatrix good = random(m, n, rng);
        for (double s : new double[] { 1.0e-200, 1.0e-100, 1.0e100, 1.0e200 }) {
            assertThrowsRuntime(scaled(base, s), ones(m));
            // and the well conditioned matrix keeps being accepted
            solve(scaled(good, s), ones(m));
        }
    }

    private static DMatrix scaled(DMatrix m, double s) {
        DMatrix out = new DMatrix(m.numRows(), m.numColumns());
        for (int j = 0; j < m.numColumns(); ++j) {
            for (int i = 0; i < m.numRows(); ++i) {
                out.setUnsafe(i, j, s * m.getUnsafe(i, j));
            }
        }
        return out;
    }

    @Test
    public void testTheZeroMatrixIsReported() {
        // Dgels short circuits an all zero A before factorizing and returns a
        // zero solution, leaving min == max == 0 and hence a NaN ratio. The
        // zero matrix does not have full rank, so it belongs in the throwing
        // branch
        assertThrowsRuntime(new DMatrix(12, 4), ones(12));
        assertThrowsRuntime(new DMatrix(4, 4), ones(4));
    }

    @Test
    public void testANonFiniteEntryIsReported() {
        Lcg rng = new Lcg(44444444L);
        for (double bad : new double[] { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY }) {
            DMatrix tall = random(10, 3, rng);
            tall.setUnsafe(4, 1, bad);
            assertThrowsRuntime(tall, ones(10));

            DMatrix square = random(6, 6, rng);
            square.setUnsafe(2, 5, bad);
            assertThrowsRuntime(square, ones(6));
        }
    }

    @Test
    public void testANonFiniteRightHandSideIsReported() {
        Lcg rng = new Lcg(55555555L);
        DMatrix A = random(10, 3, rng);
        DMatrix b = ones(10);
        b.setUnsafe(3, 0, Double.POSITIVE_INFINITY);
        assertThrowsRuntime(A, b);
    }

    /**
     * An infinity used to make the call hang rather than fail. Dlascl, which
     * Dgels reaches when the norm of its argument leaves the safe range, scales
     * by multiplying cfrom by the smallest normal number until the product is
     * small enough; for an infinite cfrom that product stays infinite and the
     * loop never terminates. The timeout is what makes a regression show up as
     * a failure instead of a stuck build.
     */
    @Test(timeout = 10000L)
    public void testAnInfiniteEntryTerminatesInsteadOfHanging() {
        Lcg rng = new Lcg(66666666L);
        DMatrix A = random(12, 4, rng);
        A.setUnsafe(7, 2, Double.POSITIVE_INFINITY);
        assertThrowsRuntime(A, ones(12));
    }
}
