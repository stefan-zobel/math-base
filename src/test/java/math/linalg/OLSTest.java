package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.StudentT;

/**
 * Tests for ordinary least squares. The class had none before it was moved off
 * the normal equations, so these cover both the estimator and the properties
 * that the previous implementation could not hold.
 */
public class OLSTest {

    private static final double ALPHA = 0.05;

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

    /** Design with an intercept column and {@code p - 1} random predictors. */
    private static DMatrix design(int n, int p, Lcg rng) {
        DMatrix X = new DMatrix(n, p);
        for (int i = 0; i < n; i++) {
            X.set(i, 0, 1.0);
            for (int j = 1; j < p; j++) {
                X.set(i, j, rng.next());
            }
        }
        return X;
    }

    private static DMatrix response(DMatrix X, double[] beta, double noise, Lcg rng) {
        DMatrix y = new DMatrix(X.numRows(), 1);
        for (int i = 0; i < X.numRows(); i++) {
            double v = 0.0;
            for (int j = 0; j < X.numColumns(); j++) {
                v += beta[j] * X.get(i, j);
            }
            y.set(i, 0, v + noise * rng.next());
        }
        return y;
    }

    @Test
    public void testExactRecoveryWithoutNoise() {
        Lcg rng = new Lcg(11L);
        double[] beta = { 2.0, -1.5, 0.75, 3.25 };
        DMatrix X = design(200, 4, rng);
        DMatrix y = response(X, beta, 0.0, rng);

        LSSummary s = OLS.estimate(ALPHA, X, y);

        for (int j = 0; j < beta.length; j++) {
            assertEquals("beta[" + j + "]", beta[j], s.getBeta().get(j), 1.0e-12);
        }
        assertEquals("R^2 on a perfect fit", 1.0, s.getRSquared(), 1.0e-12);
        for (int i = 0; i < X.numRows(); i++) {
            assertEquals("residual " + i, 0.0, s.getResiduals().get(i), 1.0e-12);
        }
    }

    @Test
    public void testAgreesWithTheNormalEquationsWhereThoseStillWork() {
        // Guard for the move off (X'X)^-1 X'y: on a well conditioned design the
        // two routes have to give the same answer.
        Lcg rng = new Lcg(23L);
        double[] beta = { 1.0, 2.0, -3.0 };
        DMatrix X = design(150, 3, rng);
        DMatrix y = response(X, beta, 0.01, rng);

        LSSummary s = OLS.estimate(ALPHA, X, y);

        DMatrix Xt = X.transpose();
        DMatrix normalEquations = Xt.mul(X).inverse().mul(Xt).mul(y);

        for (int j = 0; j < beta.length; j++) {
            assertEquals("beta[" + j + "]", normalEquations.get(j, 0), s.getBeta().get(j), 1.0e-10);
        }
    }

    @Test
    public void testIllConditionedDesignThatTheNormalEquationsCannotHandle() {
        // Vandermonde of degree 7 on [0, 1]. Squaring its condition number puts
        // X'X beyond what an LU inversion can do: the previous implementation
        // refused the problem with "Matrix A may be (close to) singular".
        int n = 60;
        int p = 8;
        DMatrix X = new DMatrix(n, p);
        double[] beta = new double[p];
        for (int j = 0; j < p; j++) {
            beta[j] = 1.0 + 0.5 * j;
        }
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            double x = i / (n - 1.0);
            double power = 1.0;
            double v = 0.0;
            for (int j = 0; j < p; j++) {
                X.set(i, j, power);
                v += beta[j] * power;
                power *= x;
            }
            y.set(i, 0, v);
        }

        try {
            X.transpose().mul(X).inverse();
            fail("the normal equations were expected to be beyond an LU inversion here");
        } catch (RuntimeException expected) {
            // that is the point of this test
        }

        LSSummary s = OLS.estimate(ALPHA, X, y);

        for (int j = 0; j < p; j++) {
            assertEquals("beta[" + j + "]", beta[j], s.getBeta().get(j), 1.0e-9);
        }
    }

    @Test
    public void testStandardErrorsAreNeverNegative() {
        // V diag(1/d^2) V' has a diagonal of sums of squares. The normal
        // equations had no such guarantee, and the old code answered a negative
        // variance with Double.MIN_NORMAL, which turned a broken computation
        // into a coefficient of overwhelming significance.
        Lcg rng = new Lcg(37L);
        double[] beta = { 1.0, 1.0, 1.0, 1.0 };
        for (int e = 0; e <= 8; e++) {
            double scale = Math.pow(10.0, e);
            DMatrix X = design(120, 4, rng);
            for (int i = 0; i < X.numRows(); i++) {
                X.set(i, 3, scale * X.get(i, 3));
            }
            DMatrix y = response(X, beta, 0.05, rng);

            LSSummary s = OLS.estimate(ALPHA, X, y);

            for (int j = 0; j < 4; j++) {
                double se = s.getCoefficientStandardErrors().get(j);
                assertTrue("se[" + j + "] at scale 1e" + e + " is " + se, se >= 0.0 && !Double.isNaN(se));
            }
        }
    }

    @Test
    public void testRankDeficientDesignIsRejected() {
        Lcg rng = new Lcg(51L);
        DMatrix X = new DMatrix(30, 3);
        DMatrix y = new DMatrix(30, 1);
        for (int i = 0; i < 30; i++) {
            double x = rng.next();
            X.set(i, 0, 1.0);
            X.set(i, 1, x);
            X.set(i, 2, x);
            y.set(i, 0, 1.0 + x);
        }

        try {
            OLS.estimate(ALPHA, X, y);
            fail("a duplicated column has to be rejected");
        } catch (RuntimeException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("rank deficient"));
        }
    }

    /**
     * The overload that takes the tolerance is the same code path as the one
     * that does not: at the tolerance the default uses, every number in the
     * summary comes back with the same bits.
     */
    @Test
    public void testTheRankToleranceOverloadChangesNothingAtTheDefault() {
        Lcg rng = new Lcg(97L);
        DMatrix X = design(40, 4, rng);
        DMatrix y = response(X, new double[] { 2.0, -1.0, 0.5, 3.0 }, 0.1, rng);

        LSSummary a = OLS.estimate(ALPHA, X, y);
        LSSummary b = OLS.estimate(ALPHA, X, y, OLS.defaultRankTolerance(X));
        for (int j = 0; j < X.numColumns(); j++) {
            assertEquals("beta " + j, a.getBeta().get(j), b.getBeta().get(j), 0.0);
            assertEquals("standard error " + j, a.getCoefficientStandardErrors().get(j),
                    b.getCoefficientStandardErrors().get(j), 0.0);
            assertEquals("t value " + j, a.getTValues().get(j), b.getTValues().get(j), 0.0);
            assertEquals("p value " + j, a.getPValues().get(j), b.getPValues().get(j), 0.0);
        }
        assertEquals(a.getRSquared(), b.getRSquared(), 0.0);
        assertEquals(a.getSigmaHatSquared(), b.getSigmaHatSquared(), 0.0);
        assertEquals(a.getConditionNumber(), b.getConditionNumber(), 0.0);
    }

    /**
     * The knob turns both ways. A tolerance far above the design's smallest
     * relative singular value refuses a design the default accepts, which is
     * what says the number is being used rather than ignored.
     */
    @Test
    public void testALooseRankToleranceRefusesAnOtherwiseUsableDesign() {
        Lcg rng = new Lcg(23L);
        DMatrix X = design(40, 4, rng);
        DMatrix y = response(X, new double[] { 1.0, 1.0, 1.0, 1.0 }, 0.1, rng);

        LSSummary fit = OLS.estimate(ALPHA, X, y);
        double smallestRelative = 1.0 / fit.getConditionNumber();
        assertTrue("this design is supposed to be well conditioned", smallestRelative > 1.0e-6);
        assertTrue("and its smallest relative singular value below the tolerance used here",
                smallestRelative < 0.5);

        try {
            OLS.estimate(ALPHA, X, y, 0.5);
            fail("a tolerance above the smallest relative singular value has to refuse");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("rank deficient"));
        }
    }

    /**
     * And a tolerance of zero, the loosest there is, still refuses a design
     * that is exactly singular -- there is no answer to reach there, and the
     * message says so rather than pointing at a tolerance that would help.
     */
    @Test
    public void testZeroRankToleranceStillRefusesADuplicatedColumn() {
        Lcg rng = new Lcg(51L);
        DMatrix X = new DMatrix(30, 3);
        DMatrix y = new DMatrix(30, 1);
        for (int i = 0; i < 30; i++) {
            double x = rng.next();
            X.set(i, 0, 1.0);
            X.set(i, 1, x);
            X.set(i, 2, x);
            y.set(i, 0, 1.0 + x);
        }

        try {
            OLS.estimate(ALPHA, X, y, 0.0);
            fail("an exactly singular design has to be rejected at any tolerance");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("exactly singular"));
        }
    }

    /** The condition number the summary reports is the one the design has. */
    @Test
    public void testTheSummaryReportsTheConditionNumberOfTheDesign() {
        Lcg rng = new Lcg(1009L);
        DMatrix X = design(50, 5, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 0.05, rng);

        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(X.copy().getArrayUnsafe(),
                X.numRows(), X.numColumns());
        double expected = svd.sigma[0] / svd.sigma[X.numColumns() - 1];
        assertEquals(expected, OLS.estimate(ALPHA, X, y).getConditionNumber(), 0.0);

        // scaling the whole design leaves the conditioning where it was
        DMatrix scaled = X.copy().scaleInplace(1.0e6);
        assertEquals(expected, OLS.estimate(ALPHA, scaled, y).getConditionNumber(), 1.0e-9 * expected);
    }

    @Test
    public void testPValueDoesNotCancelToZero() {
        // Written as 2 * (1 - cdf(|t|)) the p value underflows to exactly 0.0
        // from about |t| = 8 upwards, because cdf(|t|) rounds to 1.0. The
        // regularized incomplete beta gives the same quantity without the
        // cancellation. A p value of exactly zero is always wrong.
        Lcg rng = new Lcg(67L);
        double[] beta = { 0.5, 4.0, -2.0 };
        DMatrix X = design(120, 3, rng);
        DMatrix y = response(X, beta, 0.02, rng);

        LSSummary s = OLS.estimate(ALPHA, X, y);

        int df = s.getDegreesOfFreedom();
        StudentT tDist = new StudentT(df);
        boolean sawCancellation = false;
        for (int j = 1; j < beta.length; j++) {
            double t = Math.abs(s.getTValues().get(j));
            double p = s.getPValues().get(j);
            assertTrue("p[" + j + "] = " + p + " at |t| = " + t, p > 0.0);
            assertTrue("p[" + j + "] = " + p, p < 1.0);
            if (2.0 * (1.0 - tDist.cdf(t)) == 0.0) {
                sawCancellation = true;
            }
        }
        assertTrue("the test data must reach the range where the old expression cancels", sawCancellation);
    }

    @Test
    public void testPValueFallsAsTheStatisticGrows() {
        Lcg rng = new Lcg(83L);
        DMatrix X = design(100, 3, rng);
        double previous = Double.MAX_VALUE;
        for (double effect : new double[] { 0.01, 0.05, 0.2, 1.0 }) {
            DMatrix y = response(X, new double[] { 0.0, effect, 0.0 }, 0.1, new Lcg(97L));

            LSSummary s = OLS.estimate(ALPHA, X, y);
            double p = s.getPValues().get(1);

            assertTrue("p must fall with the effect size, saw " + p + " after " + previous, p < previous);
            previous = p;
        }
    }

    @Test
    public void testArgumentValidation() {
        Lcg rng = new Lcg(101L);
        DMatrix X = design(20, 3, rng);
        DMatrix y = response(X, new double[] { 1.0, 1.0, 1.0 }, 0.1, rng);

        expectIae(X, new DMatrix(19, 1), ALPHA, "row count mismatch");
        expectIae(design(2, 3, rng), new DMatrix(2, 1), ALPHA, "degrees of freedom below one");
        expectIae(X, y, 0.0, "alpha at zero");
        expectIae(X, y, 1.0, "alpha at one");
    }

    private static void expectIae(DMatrix X, DMatrix y, double alpha, String what) {
        try {
            OLS.estimate(alpha, X, y);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
