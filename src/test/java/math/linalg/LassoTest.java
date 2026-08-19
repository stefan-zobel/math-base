package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Tests for the lasso and the elastic net at a single penalty. Two of these
 * verify optimality outright rather than by comparison: the KKT conditions at
 * the returned point, and the closed form the solution takes when the design is
 * orthogonal.
 */
public class LassoTest {

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

    private static DMatrix predictors(int n, int p, Lcg rng) {
        DMatrix X = new DMatrix(n, p);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                X.set(i, j, rng.next());
            }
        }
        return X;
    }

    private static DMatrix response(DMatrix X, double intercept, double[] beta, double noise, Lcg rng) {
        DMatrix y = new DMatrix(X.numRows(), 1);
        for (int i = 0; i < X.numRows(); i++) {
            double v = intercept;
            for (int j = 0; j < X.numColumns(); j++) {
                v += beta[j] * X.get(i, j);
            }
            y.set(i, 0, v + noise * rng.next());
        }
        return y;
    }

    private static double mean(DMatrix y) {
        double m = 0.0;
        for (int i = 0; i < y.numRows(); i++) {
            m += y.get(i, 0);
        }
        return m / y.numRows();
    }

    /**
     * Largest violation of the stationarity conditions of the penalized
     * problem, recomputed from the reported coefficients without touching the
     * solver's own bookkeeping.
     */
    private static double kktViolation(DMatrix X, DMatrix y, Lasso.Result r, double lambda, double alpha) {
        int n = X.numRows();
        int p = X.numColumns();
        double[] xBar = new double[p];
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
            xBar[j] = m;
            scale[j] = Math.sqrt(ss / n);
        }
        double yBar = mean(y);
        double[] residual = new double[n];
        for (int i = 0; i < n; i++) {
            double v = y.get(i, 0) - yBar;
            for (int j = 0; j < p; j++) {
                v -= r.beta[j] * (X.get(i, j) - xBar[j]);
            }
            residual[i] = v;
        }

        double l1 = lambda * alpha;
        double l2 = lambda * (1.0 - alpha);
        double worst = 0.0;
        for (int j = 0; j < p; j++) {
            double g = 0.0;
            for (int i = 0; i < n; i++) {
                g += ((X.get(i, j) - xBar[j]) / scale[j]) * residual[i];
            }
            g /= n;
            // the solver works on the standardized coefficient
            double b = r.beta[j] * scale[j];
            double violation;
            if (b != 0.0) {
                violation = Math.abs(g - (((b > 0.0) ? l1 : -l1) + l2 * b));
            } else {
                violation = Math.max(0.0, Math.abs(g) - l1);
            }
            if (violation > worst) {
                worst = violation;
            }
        }
        return worst;
    }

    @Test
    public void testKarushKuhnTuckerConditionsHoldAtTheSolution() {
        // Optimality checked directly: the subgradient of the objective must
        // contain zero. Where a coefficient is non-zero the smooth gradient has
        // to balance the penalty exactly; where it is zero the gradient has to
        // be small enough for the penalty to hold it there.
        Lcg rng = new Lcg(17L);
        DMatrix X = predictors(200, 12, rng);
        double[] beta = new double[12];
        beta[0] = 3.0;
        beta[5] = -2.0;
        beta[9] = 0.5;
        DMatrix y = response(X, 1.0, beta, 0.1, rng);

        for (double alpha : new double[] { 1.0, 0.75, 0.25 }) {
            for (double lambda : new double[] { 0.5, 0.1, 0.01, 0.001 }) {
                Lasso.Result r = Lasso.estimate(X, y, lambda, alpha);
                assertTrue("must converge at lambda = " + lambda + ", alpha = " + alpha, r.converged);
                assertTrue("KKT violation at lambda = " + lambda + ", alpha = " + alpha + " : "
                        + kktViolation(X, y, r, lambda, alpha), kktViolation(X, y, r, lambda, alpha) < 1.0e-8);
            }
        }
    }

    @Test
    public void testOrthogonalDesignReproducesTheClosedForm() {
        // With orthogonal columns the coordinates decouple and the solution is
        // known in closed form, so this compares against an identity rather
        // than against a stored number.
        int n = 64;
        int p = 6;
        DMatrix X = new DMatrix(n, p);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                // Walsh pattern: pairwise orthogonal, zero mean, unit RMS
                X.set(i, j, (((i >> j) & 1) == 0) ? 1.0 : -1.0);
            }
        }
        Lcg rng = new Lcg(23L);
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            y.set(i, 0, 3.0 + 2.0 * X.get(i, 0) - 1.2 * X.get(i, 2) + 0.4 * X.get(i, 4) + 0.3 * rng.next());
        }

        double yBar = mean(y);
        double[] unpenalized = new double[p];
        for (int j = 0; j < p; j++) {
            double d = 0.0;
            for (int i = 0; i < n; i++) {
                d += X.get(i, j) * (y.get(i, 0) - yBar);
            }
            unpenalized[j] = d / n;
        }

        for (double alpha : new double[] { 1.0, 0.6 }) {
            for (double lambda : new double[] { 0.05, 0.3, 1.0 }) {
                Lasso.Result r = Lasso.estimate(X, y, lambda, alpha);
                for (int j = 0; j < p; j++) {
                    double z = unpenalized[j];
                    double g = lambda * alpha;
                    double soft = (z > g) ? (z - g) : ((z < -g) ? (z + g) : 0.0);
                    double expected = soft / (1.0 + lambda * (1.0 - alpha));
                    assertEquals("beta[" + j + "] at lambda = " + lambda + ", alpha = " + alpha, expected, r.beta[j],
                            1.0e-12);
                }
            }
        }
    }

    @Test
    public void testZeroPenaltyReproducesOrdinaryLeastSquares() {
        Lcg rng = new Lcg(13L);
        double[] beta = { 1.5, -2.0, 0.25, 3.0 };
        DMatrix X = predictors(200, 4, rng);
        DMatrix y = response(X, 4.0, beta, 0.05, rng);

        Lasso.Result r = Lasso.estimate(X, y, 0.0, 1.0);

        DMatrix Xi = new DMatrix(200, 5);
        for (int i = 0; i < 200; i++) {
            Xi.set(i, 0, 1.0);
            for (int j = 0; j < 4; j++) {
                Xi.set(i, j + 1, X.get(i, j));
            }
        }
        LSSummary ols = OLS.estimate(0.05, Xi, y);

        assertTrue(r.converged);
        assertEquals("intercept", ols.getBeta().get(0), r.intercept, 1.0e-9);
        for (int j = 0; j < 4; j++) {
            assertEquals("beta[" + j + "]", ols.getBeta().get(j + 1), r.beta[j], 1.0e-9);
        }
    }

    @Test
    public void testPureL2PenaltyAgreesWithRidge() {
        // An independent algorithm on the same problem: Ridge goes through the
        // singular values, this goes through sweeps. The conventions differ by
        // the 1/(2n) in front of the residual sum, so lambdaRidge = n * lambda
        // at alpha = 0.
        Lcg rng = new Lcg(31L);
        int n = 150;
        double[] beta = { 1.0, -0.5, 0.0, 0.0, 0.0 };
        DMatrix X = predictors(n, 5, rng);
        DMatrix y = response(X, 2.0, beta, 0.1, rng);

        double lambda = 0.37;
        Lasso.Result r = Lasso.estimate(X, y, lambda, 0.0);
        Ridge.Result ridge = Ridge.estimate(X, y, n * lambda);

        assertTrue(r.converged);
        assertEquals("intercept", ridge.intercept, r.intercept, 1.0e-8);
        for (int j = 0; j < 5; j++) {
            assertEquals("beta[" + j + "]", ridge.beta[j], r.beta[j], 1.0e-8);
        }
    }

    @Test
    public void testEverythingIsExactlyZeroAtLambdaMax() {
        // Sparsity is the point of the method, so it is asserted exactly rather
        // than against a cutoff.
        Lcg rng = new Lcg(37L);
        double[] beta = { 2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0 };
        DMatrix X = predictors(120, 8, rng);
        DMatrix y = response(X, 7.0, beta, 0.2, rng);
        double yBar = mean(y);

        for (double alpha : new double[] { 1.0, 0.5, 0.1 }) {
            double lambdaMax = Lasso.lambdaMax(X, y, alpha);
            Lasso.Result at = Lasso.estimate(X, y, lambdaMax, alpha);
            assertEquals("non-zeros at lambdaMax, alpha = " + alpha, 0, at.nonZeros);
            for (int j = 0; j < 8; j++) {
                assertEquals("beta[" + j + "] at lambdaMax", 0.0, at.beta[j], 0.0);
            }
            assertEquals("the intercept is the mean of y", yBar, at.intercept, 0.0);

            Lasso.Result below = Lasso.estimate(X, y, 0.99 * lambdaMax, alpha);
            assertTrue("just below lambdaMax something must enter, alpha = " + alpha, below.nonZeros >= 1);
        }
    }

    @Test
    public void testFitIsInvariantUnderRescalingAColumn() {
        // What the standardization buys. If centering, scaling or the back
        // transformation is wrong, this breaks loudly.
        Lcg rng = new Lcg(29L);
        double[] beta = { 1.0, -1.0, 2.0 };
        DMatrix X = predictors(120, 3, rng);
        DMatrix y = response(X, 2.0, beta, 0.1, rng);

        Lasso.Result plain = Lasso.estimate(X, y, 0.05, 0.8);

        double factor = 1000.0;
        DMatrix scaled = new DMatrix(120, 3);
        for (int i = 0; i < 120; i++) {
            for (int j = 0; j < 3; j++) {
                scaled.set(i, j, (j == 1 ? factor : 1.0) * X.get(i, j));
            }
        }
        Lasso.Result rescaled = Lasso.estimate(scaled, y, 0.05, 0.8);

        for (int i = 0; i < 120; i++) {
            assertEquals("fitted[" + i + "]", plain.fitted[i], rescaled.fitted[i], 1.0e-8);
        }
        assertEquals("the rescaled coefficient", plain.beta[1] / factor, rescaled.beta[1], 1.0e-11);
        assertEquals("untouched coefficient", plain.beta[0], rescaled.beta[0], 1.0e-8);
        assertEquals("intercept", plain.intercept, rescaled.intercept, 1.0e-8);
        assertEquals("the same variables are selected", plain.nonZeros, rescaled.nonZeros);
    }

    @Test
    public void testShiftingTheResponseMovesOnlyTheIntercept() {
        Lcg rng = new Lcg(41L);
        double[] beta = { 0.5, 1.25 };
        DMatrix X = predictors(80, 2, rng);
        DMatrix y = response(X, 1.0, beta, 0.05, rng);

        Lasso.Result base = Lasso.estimate(X, y, 0.02, 1.0);

        double shift = 17.5;
        DMatrix shifted = new DMatrix(80, 1);
        for (int i = 0; i < 80; i++) {
            shifted.set(i, 0, y.get(i, 0) + shift);
        }
        Lasso.Result moved = Lasso.estimate(X, shifted, 0.02, 1.0);

        assertEquals("intercept", base.intercept + shift, moved.intercept, 1.0e-10);
        for (int j = 0; j < 2; j++) {
            assertEquals("beta[" + j + "]", base.beta[j], moved.beta[j], 1.0e-11 * Math.abs(base.beta[j]));
        }
    }

    @Test
    public void testElasticNetKeepsACorrelatedPairThatTheLassoSplits() {
        // The reason the mixing parameter exists. Given two nearly identical
        // columns the L1 penalty is indifferent between loading one or both, so
        // it picks one; the L2 part is not indifferent and spreads the weight.
        Lcg rng = new Lcg(59L);
        int n = 200;
        DMatrix X = new DMatrix(n, 4);
        for (int i = 0; i < n; i++) {
            double base = rng.next();
            X.set(i, 0, base);
            X.set(i, 1, base + 1.0e-3 * rng.next());
            X.set(i, 2, rng.next());
            X.set(i, 3, rng.next());
        }
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            y.set(i, 0, 1.0 + 2.0 * X.get(i, 0) + 2.0 * X.get(i, 1) + 0.05 * rng.next());
        }

        Lasso.Result pure = Lasso.estimate(X, y, 0.05, 1.0);
        Lasso.Result net = Lasso.estimate(X, y, 0.05, 0.5);

        double pureGap = Math.abs(pure.beta[0] - pure.beta[1]);
        double pureMax = Math.max(Math.abs(pure.beta[0]), Math.abs(pure.beta[1]));
        assertTrue("the lasso must load the pair lopsidedly, got " + pure.beta[0] + " and " + pure.beta[1],
                pureGap > 0.5 * pureMax);

        assertTrue("the elastic net must keep both, got " + net.beta[0] + " and " + net.beta[1],
                net.beta[0] != 0.0 && net.beta[1] != 0.0);
        assertEquals("and split the weight evenly", net.beta[0], net.beta[1], 0.05 * Math.abs(net.beta[0]));

        // the disagreement is about how to distribute the effect, not about
        // how large it is: both explain the data, and the pair sums to the same
        double pureSum = pure.beta[0] + pure.beta[1];
        double netSum = net.beta[0] + net.beta[1];
        assertEquals("the combined effect of the pair", pureSum, netSum, 0.1 * Math.abs(pureSum));
        assertTrue(Math.abs(pure.rSquared - net.rSquared) < 0.01);
    }

    @Test
    public void testFittedValuesAndResidualsAreConsistent() {
        Lcg rng = new Lcg(89L);
        double[] beta = { 1.0, -0.5, 0.25, 0.0, 0.0 };
        DMatrix X = predictors(70, 5, rng);
        DMatrix y = response(X, 3.0, beta, 0.2, rng);

        Lasso.Result r = Lasso.estimate(X, y, 0.03, 0.9);

        double rss = 0.0;
        double sst = 0.0;
        double yBar = mean(y);
        int nonZeros = 0;
        for (int j = 0; j < 5; j++) {
            if (r.beta[j] != 0.0) {
                ++nonZeros;
            }
        }
        assertEquals("nonZeros counts what is in beta", nonZeros, r.nonZeros);

        for (int i = 0; i < 70; i++) {
            double v = r.intercept;
            for (int j = 0; j < 5; j++) {
                v += r.beta[j] * X.get(i, j);
            }
            assertEquals("fitted[" + i + "]", v, r.fitted[i], 1.0e-12);
            assertEquals("residual[" + i + "]", y.get(i, 0) - v, r.residuals[i], 1.0e-12);
            rss += r.residuals[i] * r.residuals[i];
            sst += (y.get(i, 0) - yBar) * (y.get(i, 0) - yBar);
        }
        assertEquals("R^2", 1.0 - rss / sst, r.rSquared, 1.0e-12);
        assertEquals("lambda is echoed", 0.03, r.lambda, 0.0);
        assertEquals("alpha is echoed", 0.9, r.alpha, 0.0);
        assertTrue("sweeps are counted", r.sweeps > 0);
        assertTrue(r.converged);
    }

    @Test
    public void testL1NormNeverGrowsWithThePenalty() {
        // A theorem about the standardized coefficients, so the reported ones
        // are scaled back up before the sum is taken.
        Lcg rng = new Lcg(97L);
        int n = 150;
        int p = 10;
        double[] beta = new double[p];
        beta[2] = 3.0;
        beta[6] = -2.0;
        DMatrix X = predictors(n, p, rng);
        DMatrix y = response(X, 0.0, beta, 0.1, rng);

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

        double previous = Double.MAX_VALUE;
        for (double lambda : new double[] { 0.001, 0.01, 0.05, 0.2, 0.5 }) {
            Lasso.Result r = Lasso.estimate(X, y, lambda, 1.0);
            double l1 = 0.0;
            for (int j = 0; j < p; j++) {
                l1 += Math.abs(r.beta[j] * scale[j]);
            }
            assertTrue("||beta||_1 must not grow at lambda = " + lambda, l1 <= previous);
            previous = l1;
        }
    }

    @Test
    public void testArgumentValidation() {
        Lcg rng = new Lcg(103L);
        DMatrix X = predictors(40, 3, rng);
        DMatrix y = response(X, 1.0, new double[] { 1.0, 1.0, 1.0 }, 0.1, rng);

        expectIae(X, new DMatrix(39, 1), 1.0, 1.0, "row count mismatch");
        expectIae(X, new DMatrix(40, 2), 1.0, 1.0, "y with two columns");
        expectIae(X, y, -1.0, 1.0, "negative lambda");
        expectIae(X, y, Double.NaN, 1.0, "lambda not a number");
        expectIae(X, y, Double.POSITIVE_INFINITY, 1.0, "infinite lambda");
        expectIae(X, y, 1.0, -0.1, "alpha below 0");
        expectIae(X, y, 1.0, 1.1, "alpha above 1");
        expectIae(X, y, 1.0, Double.NaN, "alpha not a number");

        DMatrix constant = new DMatrix(40, 2);
        for (int i = 0; i < 40; i++) {
            constant.set(i, 0, rng.next());
            constant.set(i, 1, 7.0);
        }
        expectIae(constant, y, 1.0, 1.0, "constant column");

        // alpha = 0 is a legal fit, unlike on a path where lambdaMax is unbounded
        assertTrue(Lasso.estimate(X, y, 0.5, 0.0).converged);
        // an exactly duplicated column is accepted here, which is what OLS and
        // Ridge refuse outright: the penalty makes the problem identified again
        DMatrix duplicated = new DMatrix(40, 2);
        for (int i = 0; i < 40; i++) {
            double v = rng.next();
            duplicated.set(i, 0, v);
            duplicated.set(i, 1, v);
        }
        Lasso.Result r = Lasso.estimate(duplicated, y, 0.1, 1.0);
        assertTrue(r.converged);
        assertFalse("the fit must be a number", Double.isNaN(r.rSquared));
    }

    private static void expectIae(DMatrix X, DMatrix y, double lambda, double alpha, String what) {
        try {
            Lasso.estimate(X, y, lambda, alpha);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
