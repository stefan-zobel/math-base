package math.linalg;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

import math.linalg.SymmetricJacobiEigen.Result;

public class SymmetricJacobiEigenTest {

    private static final double EPSILON = 1e-10;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /** Packs a row-wise written matrix into the column-major flat layout. */
    private static double[] pack(double[][] rows) {
        int n = rows.length;
        double[] a = new double[n * n];
        for (int i = 0; i < n; i++) {
            assertEquals("not a square matrix", n, rows[i].length);
            for (int j = 0; j < n; j++) {
                a[j * n + i] = rows[i][j];
            }
        }
        return a;
    }

    /** Deterministic symmetric matrix with entries in [-0.5, 0.5). */
    private static double[] randomSymmetric(int n, long seed) {
        double[] a = new double[n * n];
        long lcg = seed;
        for (int j = 0; j < n; j++) {
            for (int i = j; i < n; i++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double x = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
                a[j * n + i] = x;
                a[i * n + j] = x;
            }
        }
        return a;
    }

    /** Deterministic symmetric positive definite matrix {@code B^T B + n*I}. */
    private static double[] randomSpd(int n, long seed) {
        double[] b = new double[n * n];
        long lcg = seed;
        for (int idx = 0; idx < b.length; idx++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            b[idx] = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
        }
        double[] a = new double[n * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += b[i * n + k] * b[j * n + k]; // (B^T B)(i,j)
                }
                a[j * n + i] = sum;
            }
        }
        for (int i = 0; i < n; i++) {
            a[i * n + i] += n;
        }
        return a;
    }

    private static double trace(double[] a, int n) {
        double t = 0.0;
        for (int i = 0; i < n; i++) {
            t += a[i * n + i];
        }
        return t;
    }

    private static void assertGoodDecomposition(double[] a, int n, Result r, double tol) {
        assertTrue("did not converge", r.converged);
        assertEquals(n, r.n);
        assertEquals(n, r.lambda.length);
        assertEquals(n * n, r.V.length);
        assertEquals("A = V L V^T", 0.0, SymmetricJacobiEigen.reconstructionError(a, r), tol);
        assertEquals("A v = lambda v", 0.0, SymmetricJacobiEigen.residualError(a, r), tol);
        assertEquals("V orthonormal", 0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, n, n), tol);
        assertEquals("trace == sum of eigenvalues", trace(a, n), sum(r.lambda), tol);
        for (int j = 1; j < n; j++) {
            assertTrue("eigenvalues must be descending", r.lambda[j - 1] >= r.lambda[j]);
        }
    }

    private static double sum(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s;
    }

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testOneByOne() {
        double[] a = { 42.0 };
        Result r = new SymmetricJacobiEigen().decompose(a, 1);
        assertTrue(r.converged);
        assertArrayEquals(new double[] { 42.0 }, r.lambda, EPSILON);
        assertArrayEquals(new double[] { 1.0 }, r.V, EPSILON);
    }

    @Test
    public void testIdentity() {
        int n = 5;
        double[] a = new double[n * n];
        for (int i = 0; i < n; i++) {
            a[i * n + i] = 1.0;
        }
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        for (int j = 0; j < n; j++) {
            assertEquals(1.0, r.lambda[j], EPSILON);
        }
        assertGoodDecomposition(a, n, r, EPSILON);
    }

    @Test
    public void testZeroMatrix() {
        int n = 4;
        double[] a = new double[n * n];
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        assertTrue(r.converged);
        for (int j = 0; j < n; j++) {
            assertEquals(0.0, r.lambda[j], EPSILON);
        }
        // no rotation should have happened at all -> V is still the identity
        assertEquals(0.0, FlatParallelJacobiSVD.orthonormalityError(r.V, n, n), EPSILON);
    }

    @Test
    public void testDiagonalMatrixIsSortedDescending() {
        double[] a = pack(new double[][] {
            { 3.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 2.0 },
        });
        Result r = new SymmetricJacobiEigen().decompose(a, 3);
        assertArrayEquals(new double[] { 3.0, 2.0, 1.0 }, r.lambda, EPSILON);
        assertGoodDecomposition(a, 3, r, EPSILON);
    }

    @Test
    public void testTwoByTwoAntiDiagonal() {
        // [[0, 1], [1, 0]] has the eigenvalues +1 and -1. This is the case where
        // both diagonal entries vanish, i.e. theta == 0 and t must become 1.
        double[] a = pack(new double[][] {
            { 0.0, 1.0 },
            { 1.0, 0.0 },
        });
        Result r = new SymmetricJacobiEigen().decompose(a, 2);
        assertArrayEquals(new double[] { 1.0, -1.0 }, r.lambda, EPSILON);
        assertGoodDecomposition(a, 2, r, EPSILON);
    }

    @Test
    public void testIndefiniteMatrix() {
        // symmetric, with eigenvalues of both signs
        double[] a = pack(new double[][] {
            {  4.0,  1.0, -2.0,  2.0 },
            {  1.0,  2.0,  0.0,  1.0 },
            { -2.0,  0.0,  3.0, -2.0 },
            {  2.0,  1.0, -2.0, -1.0 },
        });
        Result r = new SymmetricJacobiEigen().decompose(a, 4);
        assertGoodDecomposition(a, 4, r, EPSILON);
        assertTrue("expected a positive eigenvalue", r.lambda[0] > 0.0);
        assertTrue("expected a negative eigenvalue", r.lambda[3] < 0.0);
    }

    @Test
    public void testKnownSpectrum() {
        // A = Q * diag(7, -3, 0.5) * Q^T for an explicitly orthogonal Q built from
        // a rotation in the (0,1) plane composed with one in the (1,2) plane.
        int n = 3;
        double c1 = Math.cos(0.7), s1 = Math.sin(0.7);
        double c2 = Math.cos(-1.3), s2 = Math.sin(-1.3);
        double[][] q1 = { { c1, -s1, 0 }, { s1, c1, 0 }, { 0, 0, 1 } };
        double[][] q2 = { { 1, 0, 0 }, { 0, c2, -s2 }, { 0, s2, c2 } };
        double[][] q = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    s += q1[i][k] * q2[k][j];
                }
                q[i][j] = s;
            }
        }
        double[] expected = { 7.0, -3.0, 0.5 };
        double[][] rows = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    s += q[i][k] * expected[k] * q[j][k];
                }
                rows[i][j] = s;
            }
        }
        double[] a = pack(rows);
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        Arrays.sort(expected); // ascending
        double[] descending = { expected[2], expected[1], expected[0] };
        assertArrayEquals(descending, r.lambda, 1e-12);
        assertGoodDecomposition(a, n, r, 1e-12);
    }

    @Test
    public void testRepeatedEigenvalues() {
        // diag(5, 2, 2, 2) rotated by a Householder reflection stays symmetric and
        // keeps its spectrum, but the eigenvectors of the triple eigenvalue are no
        // longer unique - only the eigenvalue equation itself may be checked.
        int n = 4;
        double[] w = { 0.5, -0.5, 0.5, 0.5 }; // unit vector
        double[] d = { 5.0, 2.0, 2.0, 2.0 };
        double[] a = new double[n * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    double qik = (i == k ? 1.0 : 0.0) - 2.0 * w[i] * w[k];
                    double qjk = (j == k ? 1.0 : 0.0) - 2.0 * w[j] * w[k];
                    s += qik * d[k] * qjk;
                }
                a[j * n + i] = s;
            }
        }
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        assertArrayEquals(new double[] { 5.0, 2.0, 2.0, 2.0 }, r.lambda, 1e-12);
        assertGoodDecomposition(a, n, r, 1e-12);
    }

    @Test
    public void testRandomSymmetricOddAndEvenDimensions() {
        for (int n = 2; n <= 33; n++) {
            double[] a = randomSymmetric(n, 987654321L + n);
            Result r = new SymmetricJacobiEigen().decompose(a, n);
            assertGoodDecomposition(a, n, r, 1e-12);
        }
    }

    @Test
    public void testInputIsNotModified() {
        int n = 12;
        double[] a = randomSymmetric(n, 4711L);
        double[] copy = a.clone();
        new SymmetricJacobiEigen().decompose(a, n);
        assertArrayEquals("decompose() must not touch its input", copy, a, 0.0);
    }

    @Test
    public void testDecomposeInPlaceMatchesDecompose() {
        int n = 10;
        double[] a = randomSymmetric(n, 20260814L);
        SymmetricJacobiEigen eig = new SymmetricJacobiEigen();
        Result r1 = eig.decompose(a, n);
        Result r2 = eig.decomposeInPlace(a.clone(), n);
        assertArrayEquals(r1.lambda, r2.lambda, 0.0);
        assertArrayEquals(r1.V, r2.V, 0.0);
    }

    @Test
    public void testParallelPathMatchesSequentialPath() {
        int n = 40;
        double[] a = randomSymmetric(n, 13579L);
        // same eps and maxSweeps, only the threshold differs
        Result seq = new SymmetricJacobiEigen(1e-15, 60, Long.MAX_VALUE).decompose(a, n);
        Result par = new SymmetricJacobiEigen(1e-15, 60, 1L).decompose(a, n);
        assertGoodDecomposition(a, n, par, 1e-12);
        // the arithmetic per rotation is identical, so the results have to agree bit
        // for bit - only the order in which independent rotations are applied differs
        assertArrayEquals(seq.lambda, par.lambda, 0.0);
        assertArrayEquals(seq.V, par.V, 0.0);
    }

    @Test
    public void testLargeMatrixOnTheParallelPath() {
        int n = 260; // n*n = 67600 -> above the default parallel threshold
        double[] a = randomSymmetric(n, 24680L);
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        assertGoodDecomposition(a, n, r, 1e-10);
    }

    @Test
    public void testSpdEigenvaluesArePositiveAndMatchTheSvd() {
        int n = 24;
        double[] a = randomSpd(n, 112233L);
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        assertGoodDecomposition(a, n, r, 1e-9);
        for (int j = 0; j < n; j++) {
            assertTrue("SPD matrix must have positive eigenvalues", r.lambda[j] > 0.0);
        }
        // for an SPD matrix the singular values are exactly the eigenvalues
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a, n, n);
        for (int j = 0; j < n; j++) {
            assertEquals(svd.sigma[j], r.lambda[j], 1e-9 * Math.max(1.0, Math.abs(svd.sigma[j])));
        }
    }

    @Test
    public void testScalingInvariance() {
        // eig(c*A) == c*eig(A); the internal power-of-two scaling must not change
        // anything beyond that
        int n = 9;
        double[] a = randomSymmetric(n, 31415L);
        double[] scaled = new double[a.length];
        final double c = 1.0e12;
        for (int idx = 0; idx < a.length; idx++) {
            scaled[idx] = c * a[idx];
        }
        SymmetricJacobiEigen eig = new SymmetricJacobiEigen();
        Result r = eig.decompose(a, n);
        Result rs = eig.decompose(scaled, n);
        for (int j = 0; j < n; j++) {
            assertEquals(c * r.lambda[j], rs.lambda[j], 1e-9 * c);
        }
    }

    @Test
    public void testEigenvectorAccessor() {
        int n = 6;
        double[] a = randomSymmetric(n, 271828L);
        Result r = new SymmetricJacobiEigen().decompose(a, n);
        for (int j = 0; j < n; j++) {
            double[] vj = r.eigenvector(j);
            assertEquals(n, vj.length);
            for (int i = 0; i < n; i++) {
                assertEquals(r.V[j * n + i], vj[i], 0.0);
            }
            // A * v_j == lambda_j * v_j
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    s += a[k * n + i] * vj[k];
                }
                assertEquals(r.lambda[j] * vj[i], s, 1e-12);
            }
        }
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testEigenvectorIndexOutOfBounds() {
        double[] a = randomSymmetric(3, 1L);
        new SymmetricJacobiEigen().decompose(a, 3).eigenvector(3);
    }

    @Test
    public void testSymmetryErrorHelper() {
        double[] symmetric = pack(new double[][] {
            { 1.0, 2.0 },
            { 2.0, 3.0 },
        });
        assertEquals(0.0, SymmetricJacobiEigen.symmetryError(symmetric, 2), 0.0);

        double[] asymmetric = pack(new double[][] {
            { 1.0, 2.0 },
            { 2.5, 3.0 },
        });
        assertEquals(0.5, SymmetricJacobiEigen.symmetryError(asymmetric, 2), 0.0);
    }

    @Test
    public void testSlightlyAsymmetricInputIsSymmetrized() {
        // the working copy is symmetrized, so the result belongs to (A + A^T)/2
        double[] a = pack(new double[][] {
            { 1.0, 2.0,  0.0 },
            { 2.4, 3.0, -1.0 },
            { 0.0, -1.0, 5.0 },
        });
        double[] mean = pack(new double[][] {
            { 1.0, 2.2,  0.0 },
            { 2.2, 3.0, -1.0 },
            { 0.0, -1.0, 5.0 },
        });
        Result r = new SymmetricJacobiEigen().decompose(a, 3);
        assertGoodDecomposition(mean, 3, r, 1e-12);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWrongLength() {
        new SymmetricJacobiEigen().decompose(new double[5], 3);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsNonPositiveDimension() {
        new SymmetricJacobiEigen().decompose(new double[0], 0);
    }
}
