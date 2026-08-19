package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Tests for ridge regression. Ridge has an unusual number of exactly checkable
 * properties, so nothing here compares against a golden value.
 */
public class RidgeTest {

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

    private static double norm(double[] a) {
        double s = 0.0;
        for (double v : a) {
            s += v * v;
        }
        return Math.sqrt(s);
    }

    @Test
    public void testZeroPenaltyReproducesOrdinaryLeastSquares() {
        Lcg rng = new Lcg(13L);
        double[] beta = { 1.5, -2.0, 0.25 };
        DMatrix X = predictors(150, 3, rng);
        DMatrix y = response(X, 4.0, beta, 0.05, rng);

        Ridge.Result r = Ridge.estimate(X, y, 0.0);

        // same problem for OLS, with an explicit intercept column in front
        DMatrix Xi = new DMatrix(150, 4);
        for (int i = 0; i < 150; i++) {
            Xi.set(i, 0, 1.0);
            for (int j = 0; j < 3; j++) {
                Xi.set(i, j + 1, X.get(i, j));
            }
        }
        LSSummary ols = OLS.estimate(0.05, Xi, y);

        assertEquals("intercept", ols.getBeta().get(0), r.intercept, 1.0e-10);
        for (int j = 0; j < 3; j++) {
            assertEquals("beta[" + j + "]", ols.getBeta().get(j + 1), r.beta[j], 1.0e-10);
        }
        assertEquals("effective df equals the column count at lambda = 0", 3.0, r.effectiveDf, 1.0e-12);
    }

    @Test
    public void testFitIsInvariantUnderRescalingAColumn() {
        // This is what the standardization inside estimate() buys. If centering,
        // scaling or the back transformation is wrong, this breaks loudly.
        Lcg rng = new Lcg(29L);
        double[] beta = { 1.0, -1.0, 2.0 };
        DMatrix X = predictors(120, 3, rng);
        DMatrix y = response(X, 2.0, beta, 0.1, rng);

        Ridge.Result plain = Ridge.estimate(X, y, 0.7);

        double factor = 1000.0;
        DMatrix scaled = new DMatrix(120, 3);
        for (int i = 0; i < 120; i++) {
            for (int j = 0; j < 3; j++) {
                scaled.set(i, j, (j == 1 ? factor : 1.0) * X.get(i, j));
            }
        }
        Ridge.Result rescaled = Ridge.estimate(scaled, y, 0.7);

        for (int i = 0; i < 120; i++) {
            assertEquals("fitted[" + i + "]", plain.fitted[i], rescaled.fitted[i], 1.0e-10);
        }
        assertEquals("the rescaled coefficient", plain.beta[1] / factor, rescaled.beta[1], 1.0e-12);
        assertEquals("untouched coefficient", plain.beta[0], rescaled.beta[0], 1.0e-10);
        assertEquals("intercept", plain.intercept, rescaled.intercept, 1.0e-10);
    }

    @Test
    public void testShiftingTheResponseMovesOnlyTheIntercept() {
        Lcg rng = new Lcg(41L);
        double[] beta = { 0.5, 1.25 };
        DMatrix X = predictors(80, 2, rng);
        DMatrix y = response(X, 1.0, beta, 0.05, rng);

        Ridge.Result base = Ridge.estimate(X, y, 0.3);

        double shift = 17.5;
        DMatrix shifted = new DMatrix(80, 1);
        for (int i = 0; i < 80; i++) {
            shifted.set(i, 0, y.get(i, 0) + shift);
        }
        Ridge.Result moved = Ridge.estimate(X, shifted, 0.3);

        assertEquals("intercept", base.intercept + shift, moved.intercept, 1.0e-10);
        for (int j = 0; j < 2; j++) {
            assertEquals("beta[" + j + "]", base.beta[j], moved.beta[j], 1.0e-12 * Math.abs(base.beta[j]));
        }
    }

    @Test
    public void testShrinkageIsMonotoneAndEndsAtTheMean() {
        Lcg rng = new Lcg(53L);
        double[] beta = { 2.0, -3.0, 1.0, 0.5 };
        DMatrix X = predictors(100, 4, rng);
        DMatrix y = response(X, 5.0, beta, 0.1, rng);

        double previousNorm = Double.MAX_VALUE;
        double previousDf = Double.MAX_VALUE;
        for (double lambda : new double[] { 0.0, 0.1, 1.0, 10.0, 100.0, 1.0e4 }) {
            Ridge.Result r = Ridge.estimate(X, y, lambda);

            assertTrue("||beta|| must fall at lambda = " + lambda, norm(r.beta) < previousNorm);
            assertTrue("effective df must fall at lambda = " + lambda, r.effectiveDf < previousDf);
            previousNorm = norm(r.beta);
            previousDf = r.effectiveDf;
        }

        Ridge.Result huge = Ridge.estimate(X, y, 1.0e14);
        assertEquals("beta vanishes under an overwhelming penalty", 0.0, norm(huge.beta), 1.0e-6);
        assertEquals("effective df vanishes with it", 0.0, huge.effectiveDf, 1.0e-8);
        double yBar = 0.0;
        for (int i = 0; i < 100; i++) {
            yBar += y.get(i, 0);
        }
        yBar /= 100.0;
        assertEquals("the intercept falls back to the mean of y", yBar, huge.intercept, 1.0e-6);
    }

    @Test
    public void testRidgeBeatsLeastSquaresOnACollinearDesign() {
        // The reason the method exists: under collinearity the unbiased
        // estimator has a variance so large that a biased one is closer to the
        // truth.
        Lcg rng = new Lcg(71L);
        int n = 60;
        double[] beta = { 1.0, 1.0, 1.0 };
        DMatrix X = new DMatrix(n, 3);
        for (int i = 0; i < n; i++) {
            double base = rng.next();
            X.set(i, 0, base);
            X.set(i, 1, base + 1.0e-3 * rng.next());
            X.set(i, 2, base + 1.0e-3 * rng.next());
        }
        DMatrix y = response(X, 0.0, beta, 0.05, rng);

        Ridge.Result ols = Ridge.estimate(X, y, 0.0);
        Ridge.Result ridge = Ridge.estimate(X, y, 1.0);

        double olsError = 0.0;
        double ridgeError = 0.0;
        for (int j = 0; j < 3; j++) {
            olsError += (ols.beta[j] - beta[j]) * (ols.beta[j] - beta[j]);
            ridgeError += (ridge.beta[j] - beta[j]) * (ridge.beta[j] - beta[j]);
        }
        assertTrue("ridge error " + Math.sqrt(ridgeError) + " is not below the least squares error "
                + Math.sqrt(olsError), ridgeError < olsError);
    }

    @Test
    public void testFittedValuesAndResidualsAreConsistent() {
        Lcg rng = new Lcg(89L);
        double[] beta = { 1.0, -0.5, 0.25 };
        DMatrix X = predictors(70, 3, rng);
        DMatrix y = response(X, 3.0, beta, 0.2, rng);

        Ridge.Result r = Ridge.estimate(X, y, 2.5);

        double rss = 0.0;
        double sst = 0.0;
        double yBar = 0.0;
        for (int i = 0; i < 70; i++) {
            yBar += y.get(i, 0);
        }
        yBar /= 70.0;
        for (int i = 0; i < 70; i++) {
            double v = r.intercept;
            for (int j = 0; j < 3; j++) {
                v += r.beta[j] * X.get(i, j);
            }
            assertEquals("fitted[" + i + "]", v, r.fitted[i], 1.0e-12);
            assertEquals("residual[" + i + "]", y.get(i, 0) - v, r.residuals[i], 1.0e-12);
            rss += r.residuals[i] * r.residuals[i];
            sst += (y.get(i, 0) - yBar) * (y.get(i, 0) - yBar);
        }
        assertEquals("R^2", 1.0 - rss / sst, r.rSquared, 1.0e-12);
        assertEquals("sigmaHat^2", rss / (70.0 - r.effectiveDf), r.sigmaHatSquared, 1.0e-12);
        for (int j = 0; j < 3; j++) {
            assertTrue("se[" + j + "]", r.standardErrors[j] > 0.0 && !Double.isNaN(r.standardErrors[j]));
        }
        assertTrue(r.converged);
    }

    @Test
    public void testArgumentValidation() {
        Lcg rng = new Lcg(103L);
        DMatrix X = predictors(40, 3, rng);
        DMatrix y = response(X, 1.0, new double[] { 1.0, 1.0, 1.0 }, 0.1, rng);

        expectIae(X, new DMatrix(39, 1), 1.0, "row count mismatch");
        expectIae(X, y, -1.0, "negative lambda");
        expectIae(X, y, Double.NaN, "lambda not a number");
        expectIae(predictors(3, 5, rng), new DMatrix(3, 1), 1.0, "more columns than rows");

        DMatrix constant = new DMatrix(40, 2);
        for (int i = 0; i < 40; i++) {
            constant.set(i, 0, rng.next());
            constant.set(i, 1, 7.0);
        }
        expectIae(constant, y, 1.0, "constant column");

        DMatrix duplicated = new DMatrix(40, 2);
        for (int i = 0; i < 40; i++) {
            double v = rng.next();
            duplicated.set(i, 0, v);
            duplicated.set(i, 1, v);
        }
        expectIae(duplicated, y, 0.0, "rank deficient at lambda = 0");
        // the same design is fine once the penalty is positive
        assertTrue(Ridge.estimate(duplicated, y, 1.0).converged);
    }

    private static void expectIae(DMatrix X, DMatrix y, double lambda, String what) {
        try {
            Ridge.estimate(X, y, lambda);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
