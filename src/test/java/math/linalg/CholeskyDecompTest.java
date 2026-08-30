package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * {@link CholeskyDecomp}, and in particular the three members that make its
 * factor usable: until they existed the class produced an {@code L} that
 * nothing in the tree consumed.
 * <p>
 * Every check is against an oracle that shares no line of code with what it
 * checks -- {@link DMatrix#inverse()}, which is an LU factorization, and
 * {@link SymmetricJacobiEigen}, which is a rotation sweep -- or against the
 * residual of the defining equation.
 */
public final class CholeskyDecompTest {

    private static final long[] SEEDS = { 1L, 42L, -7L, 987654321L };
    private static final int[] SIZES = { 1, 2, 3, 5, 12, 30 };

    private static DMatrix column(double[] v) {
        DMatrix c = new DMatrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            c.set(i, 0, v[i]);
        }
        return c;
    }

    /** A right hand side that is not a special case, from the in-test LCG. */
    private static double[] rhs(int n, long seed) {
        double[] b = new double[n];
        long lcg = seed;
        for (int i = 0; i < n; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            b[i] = ((lcg >>> 11) * 0x1.0p-53) * 4.0 - 2.0;
        }
        return b;
    }

    // ------------------------------------------------------------------
    // the solve
    // ------------------------------------------------------------------

    @Test
    public void theSolveIsTheInverseTimesTheVector() {
        // against DMatrix.inverse(), which goes through LinearEquationsSolver
        // and shares nothing with a triangular substitution
        double worst = 0.0;
        for (int n : SIZES) {
            for (long seed : SEEDS) {
                DMatrix a = PosDefiniteMatrixGenerator.generate(n, seed);
                DMatrix l = CholeskyDecomp.cholesky(a);
                double[] b = rhs(n, seed);
                double[] x = new double[n];
                CholeskyDecomp.solve(l, b, x);

                DMatrix want = a.inverse().mul(column(b));
                for (int i = 0; i < n; i++) {
                    double scale = Math.max(1.0, Math.abs(want.get(i, 0)));
                    worst = Math.max(worst, Math.abs(x[i] - want.get(i, 0)) / scale);
                }
            }
        }
        assertTrue("the solve and the inverse disagree by " + worst, worst < 1.0e-8);
    }

    @Test
    public void theSolveReproducesTheRightHandSide() {
        // the residual of the defining equation, which needs no oracle at all
        double worst = 0.0;
        for (int n : SIZES) {
            for (long seed : SEEDS) {
                DMatrix a = PosDefiniteMatrixGenerator.generate(n, seed);
                DMatrix l = CholeskyDecomp.cholesky(a);
                double[] b = rhs(n, seed);
                double[] x = new double[n];
                CholeskyDecomp.solve(l, b, x);

                DMatrix back = a.mul(column(x));
                for (int i = 0; i < n; i++) {
                    worst = Math.max(worst, Math.abs(back.get(i, 0) - b[i]) / Math.max(1.0, Math.abs(b[i])));
                }
            }
        }
        assertTrue("A x does not come back as b, worst residual " + worst, worst < 1.0e-10);
    }

    @Test
    public void theForwardSubstitutionInvertsTheFactor() {
        double worst = 0.0;
        for (int n : SIZES) {
            for (long seed : SEEDS) {
                DMatrix l = CholeskyDecomp.cholesky(PosDefiniteMatrixGenerator.generate(n, seed));
                double[] b = rhs(n, seed);
                double[] y = new double[n];
                CholeskyDecomp.solveLower(l, b, y);

                DMatrix back = l.mul(column(y));
                for (int i = 0; i < n; i++) {
                    worst = Math.max(worst, Math.abs(back.get(i, 0) - b[i]) / Math.max(1.0, Math.abs(b[i])));
                }
            }
        }
        assertTrue("L y does not come back as b, worst residual " + worst, worst < 1.0e-10);
    }

    @Test
    public void theSolveMayWriteIntoItsOwnRightHandSide() {
        for (int n : SIZES) {
            DMatrix l = CholeskyDecomp.cholesky(PosDefiniteMatrixGenerator.generate(n, 5L));
            double[] b = rhs(n, 5L);
            double[] separate = new double[n];
            CholeskyDecomp.solve(l, b, separate);
            double[] inPlace = b.clone();
            CholeskyDecomp.solve(l, inPlace, inPlace);
            for (int i = 0; i < n; i++) {
                assertEquals("n=" + n + ": aliasing changed the answer at " + i, separate[i], inPlace[i], 0.0);
            }

            double[] lower = new double[n];
            CholeskyDecomp.solveLower(l, b, lower);
            double[] lowerInPlace = b.clone();
            CholeskyDecomp.solveLower(l, lowerInPlace, lowerInPlace);
            for (int i = 0; i < n; i++) {
                assertEquals("n=" + n + ": aliasing changed the forward answer at " + i, lower[i],
                        lowerInPlace[i], 0.0);
            }
        }
    }

    @Test
    public void onlyTheLowerTriangleIsRead() {
        // the documented contract, and what makes it exact rather than checked
        int n = 8;
        DMatrix l = CholeskyDecomp.cholesky(PosDefiniteMatrixGenerator.generate(n, 11L));
        double[] b = rhs(n, 11L);
        double[] clean = new double[n];
        CholeskyDecomp.solve(l, b, clean);

        DMatrix dirty = l.copy();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                dirty.set(i, j, 1.0e6 * (i + j + 1));
            }
        }
        double[] fromDirty = new double[n];
        CholeskyDecomp.solve(dirty, b, fromDirty);
        for (int i = 0; i < n; i++) {
            assertEquals("the upper triangle was read at " + i, clean[i], fromDirty[i], 0.0);
        }
        assertEquals("the upper triangle was read by logDeterminant", CholeskyDecomp.logDeterminant(l),
                CholeskyDecomp.logDeterminant(dirty), 0.0);
    }

    // ------------------------------------------------------------------
    // the log determinant
    // ------------------------------------------------------------------

    @Test
    public void theLogDeterminantIsTheSumOfTheEigenvalues() {
        // SymmetricJacobiEigen is a rotation sweep and shares nothing with a
        // product of diagonal entries. The matrix is symmetric, so row- and
        // column-major coincide and the flat array can be handed over as it is
        double worst = 0.0;
        for (int n : SIZES) {
            for (long seed : SEEDS) {
                DMatrix a = PosDefiniteMatrixGenerator.generate(n, seed);
                double got = CholeskyDecomp.logDeterminant(CholeskyDecomp.cholesky(a));

                SymmetricJacobiEigen.Result r = new SymmetricJacobiEigen()
                        .decompose(a.copy().getArrayUnsafe(), n);
                assertTrue("the eigensolver did not converge at n=" + n, r.converged);
                double want = 0.0;
                for (int i = 0; i < n; i++) {
                    want += Math.log(r.lambda[i]);
                }
                worst = Math.max(worst, Math.abs(got - want) / Math.max(1.0, Math.abs(want)));
            }
        }
        assertTrue("the two routes to the log determinant disagree by " + worst, worst < 1.0e-10);
    }

    @Test
    public void theLogDeterminantAnswersWhereTheDeterminantHasUnderflowed() {
        // fifty dimensions scaled down: the determinant itself is below the
        // smallest double and the logarithm is an ordinary number
        int n = 50;
        DMatrix a = PosDefiniteMatrixGenerator.generate(n, 3L).scaleInplace(1.0e-14);
        DMatrix l = CholeskyDecomp.cholesky(a);
        double logDet = CholeskyDecomp.logDeterminant(l);
        assertTrue("the log determinant is not finite : " + logDet, !Double.isInfinite(logDet)
                && !Double.isNaN(logDet));
        assertEquals("the determinant should have underflowed", 0.0, Math.exp(logDet), 0.0);

        // scaling by s multiplies the determinant by s^n, so the two answers
        // must differ by exactly n log s
        double unscaled = CholeskyDecomp.logDeterminant(
                CholeskyDecomp.cholesky(PosDefiniteMatrixGenerator.generate(n, 3L)));
        assertEquals("the scaling law of a determinant", n * Math.log(1.0e-14), logDet - unscaled,
                1.0e-9 * Math.abs(n * Math.log(1.0e-14)));
    }

    @Test
    public void theLogDeterminantIsOfTheProductAndNotOfTheFactor() {
        // the one mistake that is a silent factor of two everywhere downstream
        DMatrix a = DMatrix.diag(3, 4.0);
        DMatrix l = CholeskyDecomp.cholesky(a);
        // det(A) = 64, det(L) = 8
        assertEquals("logDeterminant must be of L L^T", Math.log(64.0), CholeskyDecomp.logDeterminant(l),
                1.0e-14);
    }

    // ------------------------------------------------------------------
    // the factorization itself, which had no test at all
    // ------------------------------------------------------------------

    @Test
    public void theFactorReproducesTheMatrix() {
        double worst = 0.0;
        for (int n : SIZES) {
            for (long seed : SEEDS) {
                DMatrix a = PosDefiniteMatrixGenerator.generate(n, seed);
                DMatrix l = CholeskyDecomp.cholesky(a);
                DMatrix back = l.mulBTrans(l);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        worst = Math.max(worst, Math.abs(back.get(i, j) - a.get(i, j))
                                / Math.max(1.0, Math.abs(a.get(i, j))));
                    }
                }
                for (int i = 0; i < n; i++) {
                    for (int j = i + 1; j < n; j++) {
                        assertEquals("the factor is not lower triangular at (" + i + ", " + j + ")", 0.0,
                                l.get(i, j), 0.0);
                    }
                }
            }
        }
        assertTrue("L L^T does not reproduce A, worst residual " + worst, worst < 1.0e-12);
    }

    @Test
    public void theSymmetryCheckIsAbsoluteAndSaysSoHere() {
        // recorded rather than fixed: cholesky compares entries with an
        // absolute tolerance of 1e-10, so a matrix whose entries are large is
        // accepted with an asymmetry that would be enormous in relative terms,
        // and one that is symmetric only to relative round-off at that scale is
        // refused. Callers at large scale have to symmetrize first, which is
        // what MultivariateNormal does
        DMatrix big = new DMatrix(2, 2);
        big.set(0, 0, 1.0e6);
        big.set(1, 1, 1.0e6);
        big.set(0, 1, 1.0);
        big.set(1, 0, 1.0 + 1.0e-9);
        try {
            CholeskyDecomp.cholesky(big);
            fail("the absolute tolerance has changed; this test records the behaviour, so update it");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("symmetric"));
        }
    }

    // ------------------------------------------------------------------
    // refusals
    // ------------------------------------------------------------------

    @Test
    public void theyRejectWhatTheyCannotSolve() {
        final DMatrix l = CholeskyDecomp.cholesky(PosDefiniteMatrixGenerator.generate(3, 1L));
        final double[] three = new double[3];

        refuses("a null factor", null, three, three);
        refuses("a null right hand side", l, null, three);
        refuses("a null output", l, three, null);
        refuses("a non-square factor", new DMatrix(3, 2), three, three);
        refuses("a short right hand side", l, new double[2], three);
        refuses("a short output", l, three, new double[2]);

        DMatrix singular = new DMatrix(3, 3);
        singular.set(0, 0, 1.0);
        singular.set(1, 1, 0.0);
        singular.set(2, 2, 1.0);
        refuses("a singular factor", singular, three, three);

        try {
            CholeskyDecomp.logDeterminant(null);
            fail("logDeterminant accepted null");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage() != null);
        }
        try {
            CholeskyDecomp.logDeterminant(new DMatrix(3, 2));
            fail("logDeterminant accepted a non-square factor");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage() != null);
        }
        try {
            CholeskyDecomp.logDeterminant(singular);
            fail("logDeterminant accepted a zero on the diagonal");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage() != null);
        }
    }

    private static void refuses(String what, DMatrix l, double[] b, double[] out) {
        try {
            CholeskyDecomp.solve(l, b, out);
            fail("solve accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": solve threw without a message", expected.getMessage() != null);
        }
        try {
            CholeskyDecomp.solveLower(l, b, out);
            fail("solveLower accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": solveLower threw without a message", expected.getMessage() != null);
        }
    }
}
