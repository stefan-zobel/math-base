package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Tests for the regularization path and its cross-validation. The design used
 * throughout carries signal in three of twenty columns, which is the situation
 * the method exists for and the one that makes a wrong answer visible.
 */
public class LassoPathTest {

    /** Deterministic uniform noise in {@code [-0.5, 0.5]}. */
    private static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) * 0x1.0p-53) - 0.5;
        }
    }

    private static final int N = 300;
    private static final int P = 20;
    /** The only columns that carry signal. */
    private static final int[] SIGNAL = { 1, 7, 14 };

    private static DMatrix sparseDesign(Lcg rng) {
        DMatrix X = new DMatrix(N, P);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < P; j++) {
                X.set(i, j, rng.next());
            }
        }
        return X;
    }

    private static DMatrix sparseResponse(DMatrix X, Lcg rng) {
        DMatrix y = new DMatrix(N, 1);
        for (int i = 0; i < N; i++) {
            double v = 2.0 + 4.0 * X.get(i, 1) - 3.0 * X.get(i, 7) + 1.5 * X.get(i, 14);
            y.set(i, 0, v + 0.3 * rng.next());
        }
        return y;
    }

    private static double[] columnScales(DMatrix X) {
        int n = X.numRows();
        int p = X.numColumns();
        double[] scale = new double[p];
        for (int j = 0; j < p; j++) {
            double m = 0.0;
            for (int i = 0; i < n; i++) {
                m += X.get(i, j);
            }
            m /= n;
            double ss = 0.0;
            for (int i = 0; i < n; i++) {
                ss += (X.get(i, j) - m) * (X.get(i, j) - m);
            }
            scale[j] = Math.sqrt(ss / n);
        }
        return scale;
    }

    @Test
    public void testPathDescendsAndStartsAtTheEmptyModel() {
        Lcg rng = new Lcg(101L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        Lasso.Path path = Lasso.path(X, y, 1.0);

        assertEquals("the default grid", 100, path.lambdas.length);
        assertTrue(path.converged);
        assertEquals("alpha is echoed", 1.0, path.alpha, 0.0);
        assertEquals("it starts at lambdaMax", Lasso.lambdaMax(X, y, 1.0), path.lambdas[0], 0.0);
        assertEquals("and nothing is selected there", 0, path.nonZeros[0]);
        for (int j = 0; j < P; j++) {
            assertEquals("beta[0][" + j + "]", 0.0, path.beta[0][j], 0.0);
        }
        for (int l = 1; l < path.lambdas.length; l++) {
            assertTrue("the grid must descend at " + l, path.lambdas[l] < path.lambdas[l - 1]);
        }
        assertTrue("the last point must select something", path.nonZeros[path.lambdas.length - 1] > 0);
        assertEquals("every point has an intercept", path.lambdas.length, path.intercepts.length);
    }

    @Test
    public void testL1NormIsMonotoneAlongThePath() {
        // Non-increasing in lambda is a theorem; the number of non-zeros is
        // not, since one variable can leave the active set as another enters,
        // so that is only asserted at the two ends.
        Lcg rng = new Lcg(107L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);
        double[] scale = columnScales(X);

        Lasso.Path path = Lasso.path(X, y, 1.0);

        double previous = Double.MAX_VALUE;
        for (int l = path.lambdas.length - 1; l >= 0; l--) {
            double l1 = 0.0;
            for (int j = 0; j < P; j++) {
                l1 += Math.abs(path.beta[l][j] * scale[j]);
            }
            assertTrue("||beta||_1 must not grow towards larger lambda at " + l, l1 <= previous + 1.0e-9);
            previous = l1;
        }
        assertTrue("more is selected at the weak end than at the strong end",
                path.nonZeros[path.lambdas.length - 1] > path.nonZeros[0]);
    }

    @Test
    public void testWarmStartedPathAgreesWithIndependentFits() {
        // The warm start is only a starting point, so it must not change where
        // the sweeps end up.
        Lcg rng = new Lcg(109L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        Lasso.Path path = Lasso.path(X, y, 0.8);

        for (int l : new int[] { 20, 55, 90 }) {
            Lasso.Result cold = Lasso.estimate(X, y, path.lambdas[l], 0.8);
            assertEquals("non-zeros at grid point " + l, cold.nonZeros, path.nonZeros[l]);
            assertEquals("intercept at grid point " + l, cold.intercept, path.intercepts[l], 1.0e-8);
            for (int j = 0; j < P; j++) {
                assertEquals("beta[" + l + "][" + j + "]", cold.beta[j], path.beta[l][j], 1.0e-8);
            }
        }
    }

    @Test
    public void testCustomGrid() {
        Lcg rng = new Lcg(113L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        Lasso.Path path = Lasso.path(X, y, 1.0, 25, 1.0e-3);

        assertEquals(25, path.lambdas.length);
        assertEquals("the grid ends where it was asked to", 1.0e-3 * path.lambdas[0],
                path.lambdas[24], 1.0e-12 * path.lambdas[0]);
        // log spacing: the ratio between neighbours is constant
        double ratio = path.lambdas[1] / path.lambdas[0];
        for (int l = 2; l < 25; l++) {
            assertEquals("spacing at " + l, ratio, path.lambdas[l] / path.lambdas[l - 1], 1.0e-9);
        }

        Lasso.Path single = Lasso.path(X, y, 1.0, 1, 0.5);
        assertEquals(1, single.lambdas.length);
        assertEquals(0, single.nonZeros[0]);
    }

    @Test
    public void testCrossValidationRecoversTheTruePredictors() {
        Lcg rng = new Lcg(127L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        Lasso.CvResult cv = Lasso.cv(X, y, 1.0, 5);

        assertEquals(5, cv.folds);
        assertEquals("one error per grid point", cv.path.lambdas.length, cv.meanError.length);
        assertEquals(cv.path.lambdas.length, cv.standardError.length);
        for (int l = 0; l < cv.meanError.length; l++) {
            assertTrue("meanError[" + l + "] must be a finite non-negative number",
                    cv.meanError[l] >= 0.0 && !Double.isNaN(cv.meanError[l]) && !Double.isInfinite(cv.meanError[l]));
            assertTrue("standardError[" + l + "]", cv.standardError[l] >= 0.0 && !Double.isNaN(cv.standardError[l]));
        }

        for (int j = 0; j < SIGNAL.length; j++) {
            assertTrue("the true predictor " + SIGNAL[j] + " must survive selection",
                    cv.best.beta[SIGNAL[j]] != 0.0);
        }
        assertTrue("and most of the noise must not, got " + cv.best.nonZeros + " of " + P, cv.best.nonZeros < P / 2);

        assertTrue("lambda1se is the more conservative choice", cv.lambda1se >= cv.lambdaMin);
        Lasso.Result direct = Lasso.estimate(X, y, cv.lambdaMin, 1.0);
        assertEquals("best is the fit at lambdaMin", direct.intercept, cv.best.intercept, 0.0);
        for (int j = 0; j < P; j++) {
            assertEquals("best.beta[" + j + "]", direct.beta[j], cv.best.beta[j], 0.0);
        }
    }

    @Test
    public void testCrossValidationIsReproducibleAndSeedDependent() {
        Lcg rng = new Lcg(131L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        Lasso.CvResult a = Lasso.cv(X, y, 1.0, 5, 42L);
        Lasso.CvResult b = Lasso.cv(X, y, 1.0, 5, 42L);
        for (int l = 0; l < a.meanError.length; l++) {
            assertEquals("the same seed must give the same partition", a.meanError[l], b.meanError[l], 0.0);
        }
        assertEquals(a.lambdaMin, b.lambdaMin, 0.0);

        Lasso.CvResult other = Lasso.cv(X, y, 1.0, 5, 4711L);
        boolean differs = false;
        for (int l = 0; l < a.meanError.length; l++) {
            if (a.meanError[l] != other.meanError[l]) {
                differs = true;
                break;
            }
        }
        assertTrue("a different seed must give a different partition", differs);
        // but both must find the signal
        for (int j = 0; j < SIGNAL.length; j++) {
            assertTrue("predictor " + SIGNAL[j] + " under the second seed", other.best.beta[SIGNAL[j]] != 0.0);
        }

        // the grid is shared with the path fitted on all rows
        assertEquals(a.path.lambdas.length, a.meanError.length);
        assertEquals(a.path.lambdas[0], Lasso.lambdaMax(X, y, 1.0), 0.0);
    }

    @Test
    public void testMoreColumnsThanRows() {
        // Neither OLS nor Ridge accepts this shape at all; here it is the
        // ordinary case, and the L1 penalty is what makes the problem
        // identified again.
        Lcg rng = new Lcg(137L);
        int n = 40;
        int p = 60;
        DMatrix X = new DMatrix(n, p);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                X.set(i, j, rng.next());
            }
        }
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            y.set(i, 0, 1.0 + 2.0 * X.get(i, 3) - 1.0 * X.get(i, 40) + 0.1 * rng.next());
        }

        Lasso.Path path = Lasso.path(X, y, 1.0);
        assertTrue(path.converged);
        assertEquals("the grid stops short of the unidentified end", 1.0e-2 * path.lambdas[0],
                path.lambdas[path.lambdas.length - 1], 1.0e-12 * path.lambdas[0]);

        Lasso.CvResult cv = Lasso.cv(X, y, 1.0, 5);
        assertTrue("column 3 must be found", cv.best.beta[3] != 0.0);
        assertTrue("column 40 must be found", cv.best.beta[40] != 0.0);
        // the lasso can never select more than there are rows, and here it
        // stays far below even that
        assertTrue("out of 60 columns only a third may survive, got " + cv.best.nonZeros, cv.best.nonZeros <= 20);
    }

    @Test
    public void testArgumentValidation() {
        Lcg rng = new Lcg(139L);
        DMatrix X = sparseDesign(rng);
        DMatrix y = sparseResponse(X, rng);

        // a path has no largest useful penalty as alpha goes to zero
        expectIae("path with alpha = 0", X, y, 0.0, 100, 1.0e-4, 5, true, false);
        expectIae("cv with alpha = 0", X, y, 0.0, 100, 1.0e-4, 5, false, true);
        try {
            Lasso.lambdaMax(X, y, 0.0);
            fail("expected IllegalArgumentException for lambdaMax with alpha = 0");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message must point at Ridge", expected.getMessage().contains("Ridge"));
        }

        expectIae("nLambda below 1", X, y, 1.0, 0, 1.0e-4, 5, true, false);
        expectIae("ratio of 0", X, y, 1.0, 100, 0.0, 5, true, false);
        expectIae("ratio of 1", X, y, 1.0, 100, 1.0, 5, true, false);
        expectIae("negative ratio", X, y, 1.0, 100, -0.5, 5, true, false);
        expectIae("a single fold", X, y, 1.0, 100, 1.0e-4, 1, false, true);
        expectIae("more folds than rows", X, y, 1.0, 100, 1.0e-4, N + 1, false, true);

        // leave one out is legal, if slow, so only the boundary is checked here
        assertEquals(2, Lasso.cv(X, y, 1.0, 2).folds);
    }

    private static void expectIae(String what, DMatrix X, DMatrix y, double alpha, int nLambda, double ratio,
            int folds, boolean callPath, boolean callCv) {
        try {
            if (callPath) {
                Lasso.path(X, y, alpha, nLambda, ratio);
            }
            if (callCv) {
                Lasso.cv(X, y, alpha, folds);
            }
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
