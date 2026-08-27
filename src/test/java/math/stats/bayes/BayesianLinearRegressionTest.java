package math.stats.bayes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.List;

import org.junit.Test;

import math.distribution.MultivariateNormal;
import math.fun.DMultiFunction;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.Ridge;
import math.linalg.ScaledDesign;
import math.list.DoubleList;

/**
 * {@link BayesianLinearRegression}, checked against three oracles that share no
 * line of code with it: {@link OLS} at a flat prior, {@link Ridge} at a ridge
 * prior, and {@link BivariatePosterior} -- adaptive quadrature over the plane --
 * for the evidence when there is one coefficient.
 * <p>
 * The degrees of freedom and the constants of the evidence are the things to
 * get wrong by one here, and none of them is asserted from a derivation: the
 * oracles pin them. Measured, the posterior mean is <b>bit for bit</b> the
 * ridge estimate, the credible intervals agree with the least squares
 * confidence intervals to {@code 9.8e-12}, and the evidence agrees with
 * quadrature to {@code 2.8e-14}.
 */
public final class BayesianLinearRegressionTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private long lcg = 20260827L;

    private double uniform() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) + 0.5) * 0x1.0p-53;
    }

    private double gauss() {
        double u1 = uniform();
        double u2 = uniform();
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    }

    /**
     * A design whose columns sit at deliberately different scales -- the case
     * the standardization exists for -- and a response built from it.
     */
    private DMatrix[] data(int n, int p, double noise) {
        lcg = 20260827L;
        DMatrix x = new DMatrix(n, p);
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            double v = 3.0;
            for (int j = 0; j < p; j++) {
                double xij = (j + 1) * 2.0 + Math.pow(10.0, j - 1) * gauss();
                x.set(i, j, xij);
                v += (1.0 + j) * xij;
            }
            y.set(i, 0, v + noise * gauss());
        }
        return new DMatrix[] { x, y };
    }

    /** {@code [1 | X]}, which is the shape {@link OLS} wants. */
    private static DMatrix withIntercept(DMatrix x) {
        int n = x.numRows();
        int p = x.numColumns();
        DMatrix z = new DMatrix(n, p + 1);
        for (int i = 0; i < n; i++) {
            z.set(i, 0, 1.0);
            for (int j = 0; j < p; j++) {
                z.set(i, j + 1, x.get(i, j));
            }
        }
        return z;
    }

    // ------------------------------------------------------------------
    // the oracles
    // ------------------------------------------------------------------

    @Test
    public void aFlatPriorIsOrdinaryLeastSquares() {
        // The identity that pins the degrees of freedom. A flat prior on the
        // coefficients is NOT the limit of lambda -> 0 of a proper one -- a
        // proper prior contributes its own power of sigma^2, which cancels the
        // one integrating beta out produces, and a flat prior does not. So
        // lambda == 0 is the reference analysis, and it is the reference
        // analysis that reproduces least squares
        for (int[] shape : new int[][] { { 40, 2 }, { 60, 3 }, { 200, 5 } }) {
            int n = shape[0];
            int p = shape[1];
            String what = "n=" + n + " p=" + p;
            DMatrix[] d = data(n, p, 2.0);
            BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], 0.0);
            LSSummary ols = OLS.estimate(0.05, withIntercept(d[0]), d[1]);
            DoubleList olsBeta = ols.getBeta();

            double[] mean = new double[p];
            b.mean(mean);
            assertEquals(what + ": the intercept", olsBeta.get(0), b.intercept(),
                    1.0e-9 * Math.max(1.0, Math.abs(olsBeta.get(0))));
            for (int j = 0; j < p; j++) {
                assertEquals(what + ": coefficient " + j, olsBeta.get(j + 1), mean[j],
                        1.0e-9 * Math.max(1.0, Math.abs(olsBeta.get(j + 1))));
            }

            // exactly equal, and integers: n - p - 1 either way
            assertEquals(what + ": the degrees of freedom", ols.getDegreesOfFreedom(),
                    (int) Math.round(b.marginalDegreesOfFreedom()));
            assertEquals(what + ": the degrees of freedom are not an integer", 0.0,
                    b.marginalDegreesOfFreedom() - ols.getDegreesOfFreedom(), 0.0);

            // the posterior mean of an inverse gamma is rate/(shape-1) where a
            // least squares fit reports rate/shape, so the two differ by
            // exactly df/(df-2) and by nothing else
            double df = b.marginalDegreesOfFreedom();
            assertEquals(what + ": the variance is not the least squares one times df/(df-2)",
                    ols.getSigmaHatSquared() * df / (df - 2.0), b.variance().mean(),
                    1.0e-10 * ols.getSigmaHatSquared());

            // and the credible interval is the confidence interval, because
            // under the reference prior they are the same integral
            List<DoubleList> ci = ols.getConfidenceIntervals();
            for (int j = 0; j < p; j++) {
                BayesianLinearRegression.Interval interval = b.credibleInterval(j, 0.95);
                double lower = ci.get(j + 1).get(0);
                double upper = ci.get(j + 1).get(1);
                assertEquals(what + ": the lower end at " + j, lower, interval.lower,
                        1.0e-9 * Math.max(1.0, Math.abs(lower)));
                assertEquals(what + ": the upper end at " + j, upper, interval.upper,
                        1.0e-9 * Math.max(1.0, Math.abs(upper)));
                assertEquals(what + ": the level", 0.95, interval.level, 0.0);
                assertEquals(what + ": the width", interval.upper - interval.lower, interval.width(), 0.0);
            }
        }
    }

    @Test
    public void aRidgePriorIsRidgeRegression() {
        // the same solve, so this is bit for bit and not to a tolerance
        for (int[] shape : new int[][] { { 60, 3 }, { 120, 6 } }) {
            DMatrix[] d = data(shape[0], shape[1], 2.0);
            for (double lambda : new double[] { 1.0e-4, 0.1, 1.0, 10.0, 1000.0 }) {
                String what = "n=" + shape[0] + " p=" + shape[1] + " lambda=" + lambda;
                BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], lambda);
                Ridge.Result r = Ridge.estimate(d[0], d[1], lambda);
                double[] mean = new double[shape[1]];
                b.mean(mean);
                for (int j = 0; j < shape[1]; j++) {
                    assertEquals(what + ": coefficient " + j, r.beta[j], mean[j], 0.0);
                }
                assertEquals(what + ": the intercept", r.intercept, b.intercept(), 0.0);
                assertEquals(what + ": the effective degrees of freedom", r.effectiveDf,
                        b.effectiveDegreesOfFreedom(), 0.0);
                assertEquals(what + ": convergence", r.converged, b.converged());
                assertEquals(what + ": the penalty", lambda, b.penalty(), 0.0);
                assertEquals(what + ": the coefficient count", shape[1], b.coefficients());
            }
        }
    }

    @Test
    public void theEvidenceIsWhatQuadratureSays() {
        // The sharpest check in the file, and the one that pins every constant
        // of the closed form -- the (n-1) exponent, the -0.5 log n from the
        // flat intercept, the lambda^(p/2) from the prior. Adaptive quadrature
        // over (beta, sigma^2) shares nothing with a ratio of gamma functions
        for (int n : new int[] { 20, 50 }) {
            for (final double lambda : new double[] { 0.5, 5.0, 50.0 }) {
                DMatrix[] d = oneColumn(n);
                final ScaledDesign std = ScaledDesign.of(d[0].getArrayUnsafe(), d[1].getArrayUnsafe(), n, 1, null);
                final int rows = n;
                DMultiFunction logPosterior = new DMultiFunction() {
                    @Override
                    public double apply(double[] q) {
                        double beta = q[0];
                        double sigmaSquared = q[1];
                        if (!(sigmaSquared > 0.0)) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        double rss = 0.0;
                        for (int i = 0; i < rows; i++) {
                            double r = std.y[i] - beta * std.x[i];
                            rss += r * r;
                        }
                        return -0.5 * (rows - 1) * Math.log(2.0 * Math.PI) - 0.5 * Math.log(rows)
                                - 0.5 * (rows - 1) * Math.log(sigmaSquared) - rss / (2.0 * sigmaSquared)
                                + 0.5 * Math.log(lambda / (2.0 * Math.PI * sigmaSquared))
                                - lambda * beta * beta / (2.0 * sigmaSquared) - Math.log(sigmaSquared);
                    }
                };
                BivariatePosterior exact = BivariatePosterior.of(logPosterior,
                        new double[] { Double.NEGATIVE_INFINITY, 0.0 }, new double[] { INF, INF },
                        new double[] { 2.0, 4.0 }, 1.0e-10);
                double closed = BayesianLinearRegression.logEvidence(d[0], d[1], lambda);
                assertEquals("n=" + n + " lambda=" + lambda + ": the evidence", exact.logEvidence(), closed,
                        1.0e-11 * Math.abs(closed));
            }
        }
    }

    @Test
    public void theLaplaceApproximationIsExactHere() {
        // the posterior of beta at a known sigma^2 is Gaussian, so an
        // approximation that fits a Gaussian at the mode must return it -- a
        // third route, and the one that would catch a sign
        DMatrix[] d = data(80, 3, 2.0);
        final double lambda = 4.0;
        final double sigmaSquared = 3.0;
        final BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], lambda);
        final ScaledDesign std = ScaledDesign.of(d[0].getArrayUnsafe(), d[1].getArrayUnsafe(), 80, 3, null);

        DMultiFunction logPosterior = new DMultiFunction() {
            @Override
            public double apply(double[] q) {
                double rss = 0.0;
                for (int i = 0; i < 80; i++) {
                    double fit = 0.0;
                    for (int j = 0; j < 3; j++) {
                        fit += q[j] * std.x[j * 80 + i];
                    }
                    double r = std.y[i] - fit;
                    rss += r * r;
                }
                double penalty = 0.0;
                for (int j = 0; j < 3; j++) {
                    penalty += q[j] * q[j];
                }
                return -(rss + lambda * penalty) / (2.0 * sigmaSquared);
            }
        };
        double[] lo = { Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        double[] hi = { INF, INF, INF };
        LaplaceApproximation lap = LaplaceApproximation.of(logPosterior, lo, hi, new double[] { 0.0, 0.0, 0.0 });

        // the Laplace mode is the posterior mean, in the standardized scale
        double[] mode = new double[3];
        lap.mode(mode);
        double[] mean = new double[3];
        b.mean(mean);
        for (int j = 0; j < 3; j++) {
            assertEquals("the Laplace mode at " + j, mean[j] * std.scale[j], mode[j], 1.0e-6);
        }
        // and its covariance is sigma^2 (Z'Z + lambda I)^-1, which is what
        // covariance() holds once the scaling is undone
        DMatrix laplace = lap.asNormal().covariance();
        DMatrix mine = b.covariance();
        // scaled by the largest entry of the matrix and not by each entry:
        // an off-diagonal that is small against the diagonal is a difference
        // of nearly equal numbers once the curvature is inverted, and asserting
        // it at its own relative precision asks for what a difference quotient
        // cannot give. Measured worst against the matrix scale: 3.4e-6
        double scaleOf = 0.0;
        for (int i = 0; i < 3; i++) {
            scaleOf = Math.max(scaleOf, Math.abs(sigmaSquared * mine.get(i, i) * std.scale[i] * std.scale[i]));
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double want = sigmaSquared * mine.get(i, j) * std.scale[i] * std.scale[j];
                assertEquals("the Laplace covariance at (" + i + ", " + j + ")", want, laplace.get(i, j),
                        1.0e-4 * scaleOf);
            }
        }
    }

    // ------------------------------------------------------------------
    // the posterior itself
    // ------------------------------------------------------------------

    @Test
    public void theCovarianceIsThePosteriorOneAndNotTheEstimatorVariance() {
        // (Z'Z + lambda I)^-1 and not V diag(d^2/(d^2+lambda)^2) V', which is
        // what SvdLeastSquares.varianceMatrix returns and which would be a
        // silent substitution: the two coincide only at lambda == 0
        DMatrix[] d = data(50, 3, 2.0);
        DMatrix atZero = BayesianLinearRegression.of(d[0], d[1], 0.0).covariance();
        DMatrix atOne = BayesianLinearRegression.of(d[0], d[1], 1.0).covariance();
        boolean differ = false;
        for (int i = 0; i < 3 && !differ; i++) {
            for (int j = 0; j < 3 && !differ; j++) {
                differ = Math.abs(atZero.get(i, j) - atOne.get(i, j)) > 1.0e-9 * Math.abs(atZero.get(i, j));
            }
        }
        assertTrue("a penalty should change the posterior covariance", differ);

        for (double lambda : new double[] { 0.0, 1.0, 100.0 }) {
            DMatrix cov = BayesianLinearRegression.of(d[0], d[1], lambda).covariance();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    assertEquals("lambda=" + lambda + ": not symmetric at (" + i + ", " + j + ")", cov.get(i, j),
                            cov.get(j, i), 0.0);
                }
            }
            // positive definite, so Cholesky may not refuse it
            CholeskyDecomp.cholesky(cov);
        }
    }

    @Test
    public void theConditionalPosteriorIsAMultivariateNormal() {
        DMatrix[] d = data(50, 3, 2.0);
        BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], 2.0);
        double[] mean = new double[3];
        b.mean(mean);
        for (double sigmaSquared : new double[] { 0.5, 1.0, 7.0 }) {
            MultivariateNormal normal = b.given(sigmaSquared);
            double[] centre = new double[3];
            normal.mean(centre);
            for (int j = 0; j < 3; j++) {
                assertEquals("the conditional is not centered on the posterior mean", mean[j], centre[j], 0.0);
            }
            DMatrix cov = normal.covariance();
            DMatrix unit = b.covariance();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    assertEquals("the conditional covariance does not scale with sigma^2",
                            sigmaSquared * unit.get(i, j), cov.get(i, j), 1.0e-12 * Math.abs(cov.get(i, j)));
                }
            }
        }
    }

    @Test
    public void thePredictiveIsCenteredOnTheFit() {
        DMatrix[] d = data(60, 3, 2.0);
        double lambda = 3.0;
        BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], lambda);
        Ridge.Result r = Ridge.estimate(d[0], d[1], lambda);
        for (int i = 0; i < 60; i += 7) {
            double[] row = { d[0].get(i, 0), d[0].get(i, 1), d[0].get(i, 2) };
            assertEquals("the predictive mean is not the fitted value at row " + i, r.fitted[i],
                    b.predictiveMean(row), 1.0e-10 * Math.max(1.0, Math.abs(r.fitted[i])));
            // and the interval is centered on it and wider than the coefficient
            // interval, because it carries the error variance as well
            BayesianLinearRegression.Interval pi = b.predictionInterval(row, 0.95);
            assertEquals("the prediction interval is not centered on the fit",
                    b.predictiveMean(row), 0.5 * (pi.lower + pi.upper), 1.0e-9 * Math.max(1.0, Math.abs(r.fitted[i])));
            assertTrue("the predictive scale is not positive", b.predictiveScale(row) > 0.0);
            assertTrue("the prediction interval is not wider than the error alone",
                    pi.width() > 2.0 * Math.sqrt(b.variance().mean()));
        }
    }

    // ------------------------------------------------------------------
    // the evidence and the penalty
    // ------------------------------------------------------------------

    @Test
    public void theChosenPenaltyMaximizesTheEvidence() {
        DMatrix[] d = data(80, 4, 3.0);
        double best = BayesianLinearRegression.selectPenalty(d[0], d[1]);
        assertTrue("the chosen penalty is not positive : " + best, best > 0.0);
        double atBest = BayesianLinearRegression.logEvidence(d[0], d[1], best);
        assertTrue("the evidence at the chosen penalty is not finite", !Double.isNaN(atBest)
                && !Double.isInfinite(atBest));
        for (double factor : new double[] { 0.01, 0.1, 0.5, 0.9, 1.1, 2.0, 10.0, 100.0 }) {
            double other = BayesianLinearRegression.logEvidence(d[0], d[1], best * factor);
            assertTrue("the evidence at " + factor + " times the chosen penalty is higher : " + other + " > "
                    + atBest, other <= atBest + 1.0e-9);
        }
        // and the instance reports the same number as the static entry point
        assertEquals("the instance and the static evidence disagree", atBest,
                BayesianLinearRegression.of(d[0], d[1], best).logEvidence(), 1.0e-12 * Math.abs(atBest));
    }

    @Test
    public void theEvidenceIsRefusedWhereItWouldBeMeaningless() {
        // a flat prior on the coefficients leaves a constant that depends on
        // how many of them there are, so the one thing an evidence would be
        // used for -- comparing designs of different width -- is what it cannot
        // do. NaN, not a number that looks usable
        DMatrix[] d = data(50, 3, 2.0);
        assertTrue("the evidence at lambda == 0 must be NaN",
                Double.isNaN(BayesianLinearRegression.of(d[0], d[1], 0.0).logEvidence()));
        assertTrue("the static evidence at lambda == 0 must be NaN",
                Double.isNaN(BayesianLinearRegression.logEvidence(d[0], d[1], 0.0)));
    }

    @Test
    public void theTwoPriorsAreDifferentModelsAndNotAConstantApart() {
        // worth pinning because the opposite is the natural guess: the proper
        // prior on sigma^2 changes the shape and the rate of the posterior, and
        // the rate depends on lambda, so the two evidences are not a fixed
        // distance apart. Measured at 0.1, 1 and 100 the gap runs -7.09,
        // -11.56, -19.13
        DMatrix[] d = data(80, 4, 3.0);
        double[] gaps = new double[3];
        double[] penalties = { 0.1, 1.0, 100.0 };
        for (int k = 0; k < penalties.length; k++) {
            gaps[k] = BayesianLinearRegression.logEvidence(d[0], d[1], penalties[k], 2.0, 5.0)
                    - BayesianLinearRegression.logEvidence(d[0], d[1], penalties[k]);
        }
        assertTrue("the two priors differ by a constant, so one of them is wrong",
                Math.abs(gaps[0] - gaps[1]) > 1.0 && Math.abs(gaps[1] - gaps[2]) > 1.0);
        for (int k = 0; k < gaps.length; k++) {
            assertTrue("the proper prior evidence is not finite at " + penalties[k], !Double.isNaN(gaps[k]));
        }
    }

    @Test
    public void aTighterPriorUsesFewerEffectiveParameters() {
        // The honest statement of monotone shrinkage. It is NOT true of the
        // coefficient norm in the original scale: the columns sit at different
        // scales, so unscaling reweights the shrunk directions and the norm
        // measured 5.70, 7.74, 139.7, 4.32 over increasing lambda. What falls
        // monotonically is the effective number of parameters
        DMatrix[] d = data(80, 4, 3.0);
        double previous = Double.POSITIVE_INFINITY;
        for (double lambda : new double[] { 1.0e-3, 1.0, 100.0, 10000.0 }) {
            double df = BayesianLinearRegression.of(d[0], d[1], lambda).effectiveDegreesOfFreedom();
            assertTrue("the effective degrees of freedom rose at lambda=" + lambda, df < previous);
            assertTrue("the effective degrees of freedom left (0, p]", df > 0.0 && df <= 4.0 + 1.0e-9);
            previous = df;
        }
    }

    // ------------------------------------------------------------------
    // refusals
    // ------------------------------------------------------------------

    @Test
    public void theRegressionRejectsWhatItCannotFit() {
        DMatrix[] d = data(30, 2, 1.0);
        refuses("a null design", null, d[1], 1.0);
        refuses("a null response", d[0], null, 1.0);
        refuses("mismatched rows", new DMatrix(29, 2), d[1], 1.0);
        refuses("a response with two columns", d[0], new DMatrix(30, 2), 1.0);
        refuses("more columns than rows", new DMatrix(2, 3), new DMatrix(2, 1), 1.0);
        refuses("a negative penalty", d[0], d[1], -1.0);
        refuses("an infinite penalty", d[0], d[1], INF);
        refuses("a NaN penalty", d[0], d[1], Double.NaN);
        refuses("a constant column", constantColumn(), new DMatrix(10, 1), 1.0);

        for (double bad : new double[] { 0.0, -1.0, Double.NaN, INF }) {
            try {
                BayesianLinearRegression.of(d[0], d[1], 1.0, bad, 1.0);
                fail("a0 of " + bad + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
            try {
                BayesianLinearRegression.of(d[0], d[1], 1.0, 1.0, bad);
                fail("s0 of " + bad + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }

        BayesianLinearRegression b = BayesianLinearRegression.of(d[0], d[1], 1.0);
        for (int j : new int[] { -1, 2, Integer.MAX_VALUE }) {
            try {
                b.marginal(j);
                fail("coefficient " + j + " was accepted");
            } catch (IndexOutOfBoundsException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
        for (double level : new double[] { 0.0, 1.0, -0.5, 1.5, Double.NaN }) {
            try {
                b.credibleInterval(0, level);
                fail("a level of " + level + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
        for (double[] bad : new double[][] { null, new double[1], new double[3] }) {
            try {
                b.predictiveMean(bad);
                fail("a point of the wrong length was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
        try {
            b.given(0.0);
            fail("a zero variance was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without a message", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------

    private DMatrix[] oneColumn(int n) {
        lcg = 20260827L;
        DMatrix x = new DMatrix(n, 1);
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            double xi = 5.0 + 2.0 * gauss();
            x.set(i, 0, xi);
            y.set(i, 0, 1.0 + 2.0 * xi + 1.5 * gauss());
        }
        return new DMatrix[] { x, y };
    }

    private static DMatrix constantColumn() {
        DMatrix x = new DMatrix(10, 2);
        for (int i = 0; i < 10; i++) {
            x.set(i, 0, 3.0);
            x.set(i, 1, i);
        }
        return x;
    }

    private static void refuses(String what, DMatrix X, DMatrix y, double lambda) {
        try {
            BayesianLinearRegression.of(X, y, lambda);
            fail("the factory accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
