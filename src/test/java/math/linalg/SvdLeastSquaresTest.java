package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The checks here are against closed forms and against the equations the
 * results are supposed to satisfy, not against {@link OLS} or {@link Ridge}:
 * both of those call this class, so agreeing with them would only say that a
 * method agrees with itself.
 */
public class SvdLeastSquaresTest {

    /** A design of {@code n} rows with an intercept and {@code p} columns from an inline LCG. */
    private static double[] design(int n, int p, long seed) {
        double[] a = new double[n * (p + 1)];
        long lcg = seed;
        for (int i = 0; i < n; ++i) {
            a[i] = 1.0;
        }
        for (int j = 1; j <= p; ++j) {
            for (int i = 0; i < n; ++i) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                a[j * n + i] = ((lcg >>> 11) * 0x1.0p-53) * 4.0 - 2.0;
            }
        }
        return a;
    }

    private static double[] regressand(double[] a, int n, int cols, double[] truth, long seed) {
        double[] y = new double[n];
        long lcg = seed;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += a[j * n + i] * truth[j];
            }
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            y[i] = sum + (((lcg >>> 11) * 0x1.0p-53) - 0.5) * 0.1;
        }
        return y;
    }

    /** {@code X' X beta} against {@code X' y}, the equations least squares solves. */
    private static double normalEquationResidual(double[] a, int n, int cols, double[] y, double[] beta,
            double lambda) {
        double worst = 0.0;
        for (int j = 0; j < cols; ++j) {
            double lhs = lambda * beta[j];
            for (int k = 0; k < cols; ++k) {
                double xtx = 0.0;
                for (int i = 0; i < n; ++i) {
                    xtx += a[j * n + i] * a[k * n + i];
                }
                lhs += xtx * beta[k];
            }
            double rhs = 0.0;
            for (int i = 0; i < n; ++i) {
                rhs += a[j * n + i] * y[i];
            }
            worst = Math.max(worst, Math.abs(lhs - rhs) / Math.max(1.0, Math.abs(rhs)));
        }
        return worst;
    }

    /** The ordinary fit satisfies the normal equations it is defined by. */
    @Test
    public void solveSatisfiesTheNormalEquations() {
        int n = 40;
        int cols = 4;
        double[] a = design(n, 3, 11L);
        double[] y = regressand(a, n, cols, new double[] { 2.0, -1.0, 0.5, 3.0 }, 13L);

        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), n, cols);
        double[] beta = SvdLeastSquares.solve(svd, y, 0.0);

        assertEquals(cols, beta.length);
        assertTrue("residual " + normalEquationResidual(a, n, cols, y, beta, 0.0),
                normalEquationResidual(a, n, cols, y, beta, 0.0) < 1.0e-12);
    }

    /** The ridge fit satisfies its own, penalized normal equations. */
    @Test
    public void solveSatisfiesThePenalizedNormalEquations() {
        int n = 40;
        int cols = 4;
        double[] a = design(n, 3, 17L);
        double[] y = regressand(a, n, cols, new double[] { 1.0, 2.0, -2.0, 0.25 }, 19L);
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), n, cols);

        double[] lambdas = { 0.01, 1.0, 100.0 };
        double previous = Double.MAX_VALUE;
        for (int k = 0; k < lambdas.length; ++k) {
            double[] beta = SvdLeastSquares.solve(svd, y, lambdas[k]);
            assertTrue("lambda " + lambdas[k] + ", residual "
                    + normalEquationResidual(a, n, cols, y, beta, lambdas[k]),
                    normalEquationResidual(a, n, cols, y, beta, lambdas[k]) < 1.0e-12);
            double norm = 0.0;
            for (int j = 0; j < beta.length; ++j) {
                norm += beta[j] * beta[j];
            }
            assertTrue("a heavier penalty must not lengthen the coefficients", norm < previous);
            previous = norm;
        }
    }

    /** A line through three points, where the fit can be written down by hand. */
    @Test
    public void solveReproducesAClosedForm() {
        int n = 3;
        double[] a = { 1.0, 1.0, 1.0, 0.0, 1.0, 2.0 };
        double[] y = { 1.0, 3.0, 5.0 };

        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), n, 2);
        double[] beta = SvdLeastSquares.solve(svd, y, 0.0);

        assertEquals("intercept", 1.0, beta[0], 1.0e-12);
        assertEquals("slope", 2.0, beta[1], 1.0e-12);
    }

    /** The diagonal of the covariance is the diagonal the other method returns. */
    @Test
    public void varianceDiagonalIsTheDiagonalOfTheMatrix() {
        int n = 30;
        int cols = 3;
        double[] a = design(n, 2, 23L);
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), n, cols);

        double[] lambdas = { 0.0, 0.5, 5.0 };
        for (int k = 0; k < lambdas.length; ++k) {
            double[] diagonal = SvdLeastSquares.varianceDiagonal(svd, lambdas[k]);
            double[] matrix = SvdLeastSquares.varianceMatrix(svd, lambdas[k]);
            assertEquals(cols, diagonal.length);
            assertEquals(cols * cols, matrix.length);
            for (int j = 0; j < cols; ++j) {
                assertEquals("lambda " + lambdas[k] + ", entry " + j, matrix[j * cols + j], diagonal[j], 1.0e-15);
                assertTrue("a variance cannot be negative", diagonal[j] >= 0.0);
            }
            for (int c = 0; c < cols; ++c) {
                for (int r = 0; r < cols; ++r) {
                    assertEquals("the matrix must be symmetric", matrix[c * cols + r], matrix[r * cols + c], 0.0);
                }
            }
        }
    }

    /** At no penalty the effective degrees of freedom are the columns; they fall from there. */
    @Test
    public void effectiveDfFallsFromTheRankTowardsZero() {
        int n = 25;
        int cols = 4;
        double[] a = design(n, 3, 29L);
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), n, cols);

        assertEquals(cols, SvdLeastSquares.effectiveDf(svd, 0.0), 1.0e-12);
        double previous = cols;
        double[] lambdas = { 0.1, 1.0, 10.0, 1000.0, 1.0e8 };
        for (int k = 0; k < lambdas.length; ++k) {
            double df = SvdLeastSquares.effectiveDf(svd, lambdas[k]);
            assertTrue("df must fall with lambda: " + df + " after " + previous, df < previous);
            assertTrue("df cannot go negative", df >= 0.0);
            previous = df;
        }
        assertTrue("a huge penalty leaves nothing", previous < 0.01);
    }

    /** A repeated column is a rank deficiency, and the class says where. */
    @Test
    public void rankDeficiencyIsFound() {
        int n = 20;
        double[] full = design(n, 2, 31L);
        FlatParallelJacobiSVD.Result healthy = new FlatParallelJacobiSVD().decompose(full.clone(), n, 3);
        assertEquals("this design has full rank", -1, SvdLeastSquares.rankDeficientAt(healthy));

        double[] repeated = new double[n * 4];
        System.arraycopy(full, 0, repeated, 0, n * 3);
        System.arraycopy(full, n, repeated, n * 3, n); // column 1 again
        FlatParallelJacobiSVD.Result deficient = new FlatParallelJacobiSVD().decompose(repeated, n, 4);
        int at = SvdLeastSquares.rankDeficientAt(deficient);
        assertEquals("the fourth singular value must be the negligible one", 3, at);
    }

    /** Where the ordinary solve is unbounded, the truncated one stays finite. */
    @Test
    public void truncationSurvivesARankDeficientDesign() {
        int n = 20;
        double[] full = design(n, 2, 37L);
        double[] repeated = new double[n * 4];
        System.arraycopy(full, 0, repeated, 0, n * 3);
        System.arraycopy(full, n, repeated, n * 3, n);
        double[] y = regressand(repeated, n, 4, new double[] { 1.0, 1.0, 1.0, 0.0 }, 41L);

        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(repeated.clone(), n, 4);
        double tol = svd.sigma[0] * Math.max(n, 4) * 2.3e-16;
        double[] beta = SvdLeastSquares.solveTruncated(svd, y, tol);

        for (int j = 0; j < beta.length; ++j) {
            assertTrue("beta[" + j + "] = " + beta[j], !Double.isNaN(beta[j]) && Math.abs(beta[j]) < 1.0e3);
        }
        // the fit still reproduces the data, which is what the minimum norm solution promises
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double fitted = 0.0;
            for (int j = 0; j < 4; ++j) {
                fitted += repeated[j * n + i] * beta[j];
            }
            worst = Math.max(worst, Math.abs(fitted - y[i]));
        }
        assertTrue("worst residual " + worst, worst < 0.1);
    }

    /** A public method has to refuse what it cannot answer. */
    @Test
    public void argumentsAreValidated() {
        double[] a = design(10, 1, 43L);
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a.clone(), 10, 2);
        double[] y = new double[10];

        expect("null decomposition", null, y, 0.0);
        expect("null regressand", svd, null, 0.0);
        expect("wrong length", svd, new double[9], 0.0);
        expect("negative lambda", svd, y, -1.0);
        expect("infinite lambda", svd, y, Double.POSITIVE_INFINITY);
        expect("lambda not a number", svd, y, Double.NaN);

        try {
            SvdLeastSquares.solveTruncated(svd, y, -1.0);
            fail("expected IllegalArgumentException for a negative tolerance");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("tol"));
        }
        try {
            SvdLeastSquares.effectiveDf(null, 0.0);
            fail("expected IllegalArgumentException for a null decomposition");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("svd"));
        }
        try {
            SvdLeastSquares.rankDeficientAt(null);
            fail("expected IllegalArgumentException for a null decomposition");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("svd"));
        }
        try {
            SvdLeastSquares.varianceMatrix(svd, -0.5);
            fail("expected IllegalArgumentException for a negative penalty");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("lambda"));
        }
    }

    private static void expect(String what, FlatParallelJacobiSVD.Result svd, double[] y, double lambda) {
        try {
            SvdLeastSquares.solve(svd, y, lambda);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("the message must name the argument: " + expected.getMessage(),
                    expected.getMessage() != null && expected.getMessage().length() > 0);
        }
    }
}
