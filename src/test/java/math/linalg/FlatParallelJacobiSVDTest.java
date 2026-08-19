package math.linalg;

import static org.junit.Assert.*;

import org.junit.Test;

import math.linalg.FlatParallelJacobiSVD.Result;

public class FlatParallelJacobiSVDTest {

    private static final double EPSILON = 1e-10;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /** Packs a row-wise written matrix into the column-major flat layout. */
    private static double[] pack(double[][] rows) {
        int m = rows.length;
        int n = rows[0].length;
        double[] a = new double[m * n];
        for (int i = 0; i < m; i++) {
            assertEquals("not rectangular", n, rows[i].length);
            for (int j = 0; j < n; j++) {
                a[j * m + i] = rows[i][j];
            }
        }
        return a;
    }

    /** Deterministic m x n matrix (column-major) with entries in [-0.5, 0.5). */
    private static double[] random(int m, int n, long seed) {
        double[] a = new double[m * n];
        long lcg = seed;
        for (int idx = 0; idx < a.length; idx++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            a[idx] = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
        }
        return a;
    }

    /** {@code A^T A} as an n x n column-major matrix. */
    private static double[] gramian(double[] a, int m, int n) {
        double[] g = new double[n * n];
        for (int p = 0; p < n; p++) {
            for (int q = 0; q < n; q++) {
                double dot = 0.0;
                for (int i = 0; i < m; i++) {
                    dot += a[p * m + i] * a[q * m + i];
                }
                g[q * n + p] = dot;
            }
        }
        return g;
    }

    private static void assertGoodDecomposition(double[] a, int m, int n, Result r, double tol) {
        assertTrue("did not converge", r.converged);
        assertEquals(m, r.m);
        assertEquals(n, r.n);
        assertEquals(n, r.sigma.length);
        assertEquals(m * n, r.U.length);
        assertEquals(n * n, r.V.length);
        assertEquals("A = U S V^T", 0.0, FlatParallelJacobiSVD.reconstructionError(a, r), tol);
        assertEquals("U orthonormal", 0.0, FlatParallelJacobiSVD.orthonormalityError(r.U, m, n), tol);
        assertEquals("V orthonormal", 0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, n, n), tol);
        for (int j = 0; j < n; j++) {
            assertTrue("singular values must be non-negative", r.sigma[j] >= 0.0);
        }
        for (int j = 1; j < n; j++) {
            assertTrue("singular values must be descending", r.sigma[j - 1] >= r.sigma[j]);
        }
    }

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testIdentity() {
        int n = 6;
        double[] a = new double[n * n];
        for (int i = 0; i < n; i++) {
            a[i * n + i] = 1.0;
        }
        Result r = new FlatParallelJacobiSVD().decompose(a, n, n);
        for (int j = 0; j < n; j++) {
            assertEquals(1.0, r.sigma[j], EPSILON);
        }
        assertGoodDecomposition(a, n, n, r, EPSILON);
    }

    @Test
    public void testDiagonalMatrixIsSortedDescendingAndNonNegative() {
        // negative diagonal entries must show up as positive singular values
        double[] a = pack(new double[][] {
            { 2.0,  0.0, 0.0 },
            { 0.0, -5.0, 0.0 },
            { 0.0,  0.0, 3.0 },
        });
        Result r = new FlatParallelJacobiSVD().decompose(a, 3, 3);
        assertArrayEquals(new double[] { 5.0, 3.0, 2.0 }, r.sigma, EPSILON);
        assertGoodDecomposition(a, 3, 3, r, EPSILON);
    }

    @Test
    public void testZeroMatrix() {
        int m = 4, n = 3;
        double[] a = new double[m * n];
        Result r = new FlatParallelJacobiSVD().decompose(a, m, n);
        assertTrue(r.converged);
        for (int j = 0; j < n; j++) {
            assertEquals(0.0, r.sigma[j], 0.0);
        }
        // V is untouched, so it is still the identity; U stays zero (no direction is
        // defined), which is why orthonormality of U cannot be asserted here
        assertEquals(0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, n, n), EPSILON);
    }

    @Test
    public void testSingleColumn() {
        // n = 1: sigma is just the column norm, U the normalized column, V = [1]
        double[] a = { 3.0, 4.0 };
        Result r = new FlatParallelJacobiSVD().decompose(a, 2, 1);
        assertEquals(5.0, r.sigma[0], EPSILON);
        assertArrayEquals(new double[] { 0.6, 0.8 }, r.U, EPSILON);
        assertArrayEquals(new double[] { 1.0 }, r.V, EPSILON);
        assertGoodDecomposition(a, 2, 1, r, EPSILON);
    }

    @Test
    public void testKnownSingularValuesOfTwoByTwo() {
        // [[1, 1], [0, 1]] has the singular values (sqrt(5) +- 1) / 2
        double[] a = pack(new double[][] {
            { 1.0, 1.0 },
            { 0.0, 1.0 },
        });
        Result r = new FlatParallelJacobiSVD().decompose(a, 2, 2);
        double sqrt5 = Math.sqrt(5.0);
        assertEquals((sqrt5 + 1.0) / 2.0, r.sigma[0], 1e-12);
        assertEquals((sqrt5 - 1.0) / 2.0, r.sigma[1], 1e-12);
        assertGoodDecomposition(a, 2, 2, r, 1e-12);
    }

    @Test
    public void testTallAndSquareRandomMatrices() {
        for (int n = 1; n <= 12; n++) {
            for (int m = n; m <= n + 7; m++) {
                double[] a = random(m, n, 1000L + 31L * m + n);
                Result r = new FlatParallelJacobiSVD().decompose(a, m, n);
                assertGoodDecomposition(a, m, n, r, 1e-12);
            }
        }
    }

    @Test
    public void testRankDeficientMatrix() {
        // column 2 = 2 * column 1 -> rank 1, so the second singular value vanishes
        double[] a = pack(new double[][] {
            { 1.0, 2.0 },
            { 2.0, 4.0 },
            { 3.0, 6.0 },
        });
        Result r = new FlatParallelJacobiSVD().decompose(a, 3, 2);
        assertTrue(r.converged);
        assertEquals(Math.sqrt(14.0 * 5.0), r.sigma[0], 1e-12); // ||A||_F for a rank-1 matrix
        assertEquals(0.0, r.sigma[1], 1e-14 * r.sigma[0]);
        assertEquals(0.0, FlatParallelJacobiSVD.reconstructionError(a, r), 1e-13);
        // V is a product of rotations and stays orthogonal even for a singular A
        assertEquals(0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, 2, 2), EPSILON);
    }

    @Test
    public void testSingularValuesMatchEigenvaluesOfTheGramian() {
        // sigma_k(A) == sqrt(lambda_k(A^T A)) - cross-checked against an independent
        // algorithm (two-sided Jacobi on the Gramian)
        int m = 30, n = 14;
        double[] a = random(m, n, 777L);
        Result svd = new FlatParallelJacobiSVD().decompose(a, m, n);
        SymmetricJacobiEigen.Result eig =
                new SymmetricJacobiEigen().decompose(gramian(a, m, n), n);
        for (int k = 0; k < n; k++) {
            double expected = Math.sqrt(Math.max(0.0, eig.lambda[k]));
            assertEquals(expected, svd.sigma[k], 1e-9 * Math.max(1.0, svd.sigma[0]));
        }
    }

    @Test
    public void testScalingInvariance() {
        // sigma(c*A) == |c| * sigma(A); the internal power-of-two scaling must not
        // change anything beyond that
        int m = 11, n = 7;
        double[] a = random(m, n, 31415L);
        final double c = 1.0e-11;
        double[] scaled = new double[a.length];
        for (int idx = 0; idx < a.length; idx++) {
            scaled[idx] = c * a[idx];
        }
        FlatParallelJacobiSVD svd = new FlatParallelJacobiSVD();
        Result r = svd.decompose(a, m, n);
        Result rs = svd.decompose(scaled, m, n);
        for (int k = 0; k < n; k++) {
            assertEquals(c * r.sigma[k], rs.sigma[k], 1e-9 * c);
        }
    }

    @Test
    public void testInputIsNotModified() {
        int m = 9, n = 5;
        double[] a = random(m, n, 4711L);
        double[] copy = a.clone();
        new FlatParallelJacobiSVD().decompose(a, m, n);
        assertArrayEquals("decompose() must not touch its input", copy, a, 0.0);
    }

    @Test
    public void testDecomposeInPlaceAliasesTheInputAndMatchesDecompose() {
        int m = 9, n = 5;
        double[] a = random(m, n, 20260814L);
        FlatParallelJacobiSVD svd = new FlatParallelJacobiSVD();
        Result r1 = svd.decompose(a, m, n);
        double[] work = a.clone();
        Result r2 = svd.decomposeInPlace(work, m, n);
        assertSame("Result.U must alias the working array", work, r2.U);
        assertArrayEquals(r1.sigma, r2.sigma, 0.0);
        assertArrayEquals(r1.U, r2.U, 0.0);
        assertArrayEquals(r1.V, r2.V, 0.0);
    }

    @Test
    public void testParallelPathMatchesSequentialPath() {
        int m = 60, n = 40;
        double[] a = random(m, n, 13579L);
        // same eps and maxSweeps, only the threshold differs
        Result seq = new FlatParallelJacobiSVD(1e-15, 60, Long.MAX_VALUE).decompose(a, m, n);
        Result par = new FlatParallelJacobiSVD(1e-15, 60, 1L).decompose(a, m, n);
        assertGoodDecomposition(a, m, n, par, 1e-12);
        // rotations of one round act on disjoint columns, so the arithmetic is
        // identical regardless of the order in which they are applied
        assertArrayEquals(seq.sigma, par.sigma, 0.0);
        assertArrayEquals(seq.U, par.U, 0.0);
        assertArrayEquals(seq.V, par.V, 0.0);
    }

    @Test
    public void testLargeMatrixOnTheParallelPath() {
        int m = 300, n = 220; // m*n = 66000 -> above the default parallel threshold
        double[] a = random(m, n, 24680L);
        Result r = new FlatParallelJacobiSVD().decompose(a, m, n);
        assertGoodDecomposition(a, m, n, r, 1e-11);
    }

    @Test
    public void testReconstructRoundTrip() {
        int m = 7, n = 4;
        double[] a = random(m, n, 999L);
        Result r = new FlatParallelJacobiSVD().decompose(a, m, n);
        double[] rec = FlatParallelJacobiSVD.reconstruct(r);
        assertEquals(a.length, rec.length);
        assertArrayEquals(a, rec, 1e-13);
    }

    @Test
    public void testOrthonormalityErrorHelper() {
        // a 3 x 2 matrix with orthonormal columns
        double[] q = pack(new double[][] {
            { 1.0, 0.0 },
            { 0.0, 1.0 },
            { 0.0, 0.0 },
        });
        assertEquals(0.0, FlatParallelJacobiSVD.orthonormalityError(q, 3, 2), 0.0);

        // columns of equal direction -> off-diagonal entry 1 instead of 0
        double[] parallelCols = pack(new double[][] {
            { 1.0, 1.0 },
            { 0.0, 0.0 },
            { 0.0, 0.0 },
        });
        assertEquals(1.0, FlatParallelJacobiSVD.orthonormalityError(parallelCols, 3, 2), 0.0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWideMatrix() {
        new FlatParallelJacobiSVD().decompose(new double[6], 2, 3);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWrongLength() {
        new FlatParallelJacobiSVD().decompose(new double[5], 3, 2);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testInPlaceRejectsWideMatrix() {
        new FlatParallelJacobiSVD().decomposeInPlace(new double[6], 2, 3);
    }

    @Test
    public void testColumnPairOfExactlyEqualNormIsRotated() {
        // Regression: the rotation angle came from Math.signum(tau), which is
        // 0.0 at tau == 0 rather than +1. Two columns of exactly equal norm
        // give tau == 0, so the rotation collapsed to the identity, the pair
        // was never orthogonalized, the sweeps ran to the limit and the second
        // singular value came back equal to the first instead of zero.
        int m = 40;
        double[] a = new double[m * 2];
        double[] col = random(m, 1, 4711L);
        for (int i = 0; i < m; i++) {
            a[i] = col[i];
            a[m + i] = col[i];
        }

        Result r = new FlatParallelJacobiSVD().decompose(a, m, 2);

        assertTrue("a duplicated column pair has to converge", r.converged);
        assertEquals("the second singular value of a rank one matrix", 0.0, r.sigma[1], EPSILON);
        assertTrue("the first singular value must not vanish", r.sigma[0] > 0.0);
        assertEquals("A = U S V^T", 0.0, FlatParallelJacobiSVD.reconstructionError(a, r), EPSILON);
        // as in testRankDeficientMatrix, only V stays orthogonal for a singular A
        assertEquals("V orthonormal", 0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, 2, 2), EPSILON);
    }
}
