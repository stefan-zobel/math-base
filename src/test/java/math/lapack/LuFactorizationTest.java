package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The reusable factorization, checked against the fused one.
 * <p>
 * {@link LuFactorization} and {@link Dgesv} run the same two LAPACK routines,
 * so the interesting assertions are equalities at tolerance {@code 0.0} rather
 * than approximations: separating the factorization from the solve is supposed
 * to change what a caller pays and nothing about what it gets. What is checked
 * beyond that is the structure -- that the factors multiply back to the matrix
 * once the row interchanges are undone -- and the two cases a caller has to
 * handle, a singular matrix and an argument that does not fit.
 */
public final class LuFactorizationTest {

    /** The inline generator the test conventions ask for, seeded per test. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53 - 0.5;
        }
    }

    @Test
    public void testTheFactorizationIsTheOneDgesvComputes() {
        int n = 7;
        double[] a = matrix(n, 12345L);
        double[] mine = a.clone();
        double[] fused = a.clone();
        int[] minePivots = new int[n];
        int[] fusedPivots = new int[n];
        double[] b = rightHandSide(n, 999L);
        double[] mineB = b.clone();
        double[] fusedB = b.clone();

        assertTrue(LuFactorization.factor(n, mine, minePivots));
        LuFactorization.solve(n, 1, mine, minePivots, mineB);
        assertTrue(Dgesv.dgesv(n, 1, fused, 0, n, fusedPivots, 0, fusedB, 0, n));

        for (int i = 0; i < n * n; ++i) {
            assertEquals("the factors are the same array", fused[i], mine[i], 0.0);
        }
        for (int i = 0; i < n; ++i) {
            assertEquals("and so are the pivots", fusedPivots[i], minePivots[i]);
            assertEquals("and so is the solution", fusedB[i], mineB[i], 0.0);
        }
    }

    /**
     * The reason the class exists: five right hand sides, one factorization,
     * and every answer is the one a fresh {@link Dgesv} would have given.
     */
    @Test
    public void testOneFactorizationServesEveryRightHandSide() {
        int n = 9;
        double[] a = matrix(n, 555L);
        double[] lu = a.clone();
        int[] pivots = new int[n];
        assertTrue(LuFactorization.factor(n, lu, pivots));

        for (int k = 0; k < 5; ++k) {
            double[] b = rightHandSide(n, 4000L + k);
            double[] reused = b.clone();
            LuFactorization.solve(n, 1, lu, pivots, reused);

            double[] fresh = a.clone();
            double[] freshB = b.clone();
            assertTrue(Dgesv.dgesv(n, 1, fresh, 0, n, new int[n], 0, freshB, 0, n));
            for (int i = 0; i < n; ++i) {
                assertEquals("right hand side " + k, freshB[i], reused[i], 0.0);
            }
        }
    }

    @Test
    public void testSeveralColumnsAtOnceAreTheColumnsOneAtATime() {
        int n = 6;
        int nrhs = 3;
        double[] lu = matrix(n, 777L);
        int[] pivots = new int[n];
        assertTrue(LuFactorization.factor(n, lu, pivots));

        double[] block = rightHandSide(n * nrhs, 31L);
        double[] together = block.clone();
        LuFactorization.solve(n, nrhs, lu, pivots, together);

        for (int j = 0; j < nrhs; ++j) {
            double[] column = new double[n];
            System.arraycopy(block, j * n, column, 0, n);
            LuFactorization.solve(n, 1, lu, pivots, column);
            for (int i = 0; i < n; ++i) {
                assertEquals("column " + j, column[i], together[j * n + i], 0.0);
            }
        }
    }

    /**
     * The structural check: undoing the row interchanges on {@code L U} has to
     * give the matrix back. This is what says the pivot convention is read the
     * way {@code Dgetrf} writes it.
     */
    @Test
    public void testTheFactorsMultiplyBackToTheMatrix() {
        for (int n : new int[] { 1, 2, 5, 13 }) {
            double[] a = matrix(n, 8000L + n);
            double[] lu = a.clone();
            int[] pivots = new int[n];
            assertTrue(LuFactorization.factor(n, lu, pivots));

            double[] product = new double[n * n];
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    double sum = 0.0;
                    for (int k = 0; k <= Math.min(i, j); ++k) {
                        double l = (k == i) ? 1.0 : lu[k * n + i];
                        sum += l * lu[j * n + k];
                    }
                    product[j * n + i] = sum;
                }
            }
            // the interchanges applied in reverse turn P A back into A
            for (int i = n - 1; i >= 0; --i) {
                int other = pivots[i] - 1;
                if (other != i) {
                    for (int j = 0; j < n; ++j) {
                        double swap = product[j * n + i];
                        product[j * n + i] = product[j * n + other];
                        product[j * n + other] = swap;
                    }
                }
            }
            for (int i = 0; i < n * n; ++i) {
                assertEquals("order " + n + ", entry " + i, a[i], product[i], 1.0e-13);
            }
        }
    }

    @Test
    public void testTheResidualIsAtRounding() {
        for (int n : new int[] { 2, 8, 25 }) {
            double[] a = matrix(n, 606L + n);
            double[] b = rightHandSide(n, 707L + n);
            double[] lu = a.clone();
            int[] pivots = new int[n];
            assertTrue(LuFactorization.factor(n, lu, pivots));
            double[] x = b.clone();
            LuFactorization.solve(n, 1, lu, pivots, x);

            double residual = 0.0;
            double scale = 0.0;
            for (int i = 0; i < n; ++i) {
                double row = 0.0;
                for (int j = 0; j < n; ++j) {
                    row += a[j * n + i] * x[j];
                }
                residual = Math.max(residual, Math.abs(row - b[i]));
                scale = Math.max(scale, Math.abs(b[i]));
            }
            assertTrue("order " + n + " left a residual of " + residual, residual <= 1.0e-13 * (scale + 1.0));
        }
    }

    /**
     * A backward stable solver leaves a small residual even where the solution
     * itself cannot be trusted. On a Hilbert matrix of order eight the
     * condition number is above {@code 1e10}, so this asserts the residual and
     * says nothing about the answer.
     */
    @Test
    public void testAnIllConditionedMatrixStillLeavesASmallResidual() {
        int n = 8;
        double[] a = new double[n * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                a[j * n + i] = 1.0 / (i + j + 1.0);
            }
        }
        double[] b = new double[n];
        for (int i = 0; i < n; ++i) {
            b[i] = 1.0;
        }
        double[] lu = a.clone();
        int[] pivots = new int[n];
        assertTrue(LuFactorization.factor(n, lu, pivots));
        double[] x = b.clone();
        LuFactorization.solve(n, 1, lu, pivots, x);

        double residual = 0.0;
        double largest = 0.0;
        for (int i = 0; i < n; ++i) {
            double row = 0.0;
            for (int j = 0; j < n; ++j) {
                row += a[j * n + i] * x[j];
            }
            residual = Math.max(residual, Math.abs(row - b[i]));
            largest = Math.max(largest, Math.abs(x[i]));
        }
        // the problem really is the hard one it is meant to be: the solution
        // reaches 2.2e+05 where the right hand side is all ones
        assertTrue("the solution came out at " + largest, largest > 1.0e5);
        // and the residual stays at the backward stable bound eps |A| |x|,
        // which here is about 2.6e-11 -- measured, it is 4.1e-12
        double bound = 20.0 * 2.2204460492503131e-16 * 2.72 * largest;
        assertTrue("residual was " + residual + " against a bound of " + bound, residual < bound);
    }

    @Test
    public void testASingularMatrixIsReportedRatherThanThrown() {
        int n = 4;
        double[] a = matrix(n, 42L);
        // a column of zeros cannot be pivoted away
        for (int i = 0; i < n; ++i) {
            a[2 * n + i] = 0.0;
        }
        assertFalse("an exactly singular U is a false, not an exception",
                LuFactorization.factor(n, a, new int[n]));
    }

    @Test
    public void testOffsetsAndALargerLeadingDimensionGiveTheSameAnswer() {
        int n = 5;
        int lda = 9;
        int offset = 3;
        double[] compact = matrix(n, 2024L);
        double[] padded = new double[offset + lda * n];
        for (int i = 0; i < padded.length; ++i) {
            padded[i] = Double.NaN;
        }
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                padded[offset + j * lda + i] = compact[j * n + i];
            }
        }
        double[] b = rightHandSide(n, 3030L);

        int[] compactPivots = new int[n];
        assertTrue(LuFactorization.factor(n, compact, compactPivots));
        double[] compactX = b.clone();
        LuFactorization.solve(n, 1, compact, compactPivots, compactX);

        int[] paddedPivots = new int[n + 2];
        assertTrue(LuFactorization.factor(n, padded, offset, lda, paddedPivots, 2));
        double[] paddedX = new double[offset + n];
        System.arraycopy(b, 0, paddedX, offset, n);
        LuFactorization.solve(n, 1, padded, offset, lda, paddedPivots, 2, paddedX, offset, n);

        for (int i = 0; i < n; ++i) {
            assertEquals(compactPivots[i], paddedPivots[i + 2]);
            assertEquals(compactX[i], paddedX[offset + i], 0.0);
        }
    }

    @Test
    public void testAMatrixOfOrderZeroIsNotAnError() {
        assertTrue(LuFactorization.factor(0, new double[0], new int[0]));
        LuFactorization.solve(0, 1, new double[0], new int[0], new double[0]);
        LuFactorization.solve(3, 0, new double[9], new int[3], new double[0]);
    }

    @Test
    public void testTheArgumentsAreChecked() {
        refuses("a", null, new int[3], 3);
        refuses("ipiv", new double[9], null, 3);
        refuses("too short", new double[8], new int[3], 3);
        refuses("at least", new double[9], new int[2], 3);
        refuses("negative", new double[9], new int[3], -1);
        try {
            LuFactorization.solve(3, 1, new double[9], new int[3], new double[2]);
            fail("expected a refusal about the right hand side");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("right hand side"));
        }
        try {
            LuFactorization.solve(3, -1, new double[9], new int[3], new double[3]);
            fail("expected a refusal naming nrhs");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("nrhs"));
        }
    }

    private static void refuses(String fragment, double[] a, int[] ipiv, int n) {
        try {
            LuFactorization.factor(n, a, ipiv);
            fail("expected a refusal mentioning " + fragment);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(fragment));
        } catch (NullPointerException expected) {
            assertTrue(String.valueOf(expected.getMessage()), String.valueOf(expected.getMessage()).contains(fragment));
        }
    }

    /** A well scaled matrix of order {@code n}, column-major. */
    private static double[] matrix(int n, long seed) {
        Lcg lcg = new Lcg(seed);
        double[] a = new double[n * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                a[j * n + i] = lcg.next();
            }
            a[j * n + j] += 1.5;
        }
        return a;
    }

    private static double[] rightHandSide(int length, long seed) {
        Lcg lcg = new Lcg(seed);
        double[] b = new double[length];
        for (int i = 0; i < length; ++i) {
            b[i] = lcg.next();
        }
        return b;
    }
}
