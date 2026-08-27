package math.stats.bayes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.cern.FastGamma;
import math.fun.DBiFunction;
import math.fun.DMultiFunction;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.distribution.InverseGamma;
import math.distribution.StudentT;
import math.solve.Quadrature;

/**
 * {@link BivariatePosterior}, which normalizes a two-parameter log posterior by
 * quadrature after whitening the parameters at the mode.
 * <p>
 * Checked the way {@link ScalarPosteriorTest} checks its class: hand it a law
 * the library already knows and demand the law back. The correlated normal pair
 * is the workhorse, because every part of the answer -- evidence, mode,
 * curvature, both moments and both marginals -- is closed form for it and the
 * correlation is exactly the thing the class exists to survive. The
 * Normal-InverseGamma at the end is a posterior it has to <em>derive</em>.
 * <p>
 * Tolerances were measured before they were written. Over the eight pairs below
 * the worst figures seen were {@code 1.3e-13} on the log evidence,
 * {@code 9.1e-15} on the agreement of the two axis orders and {@code 2.7e-12}
 * relative on a marginal density; what is asserted is about two orders of
 * magnitude above each.
 */
public final class BivariatePosteriorTest {

    private static final double INF = Double.POSITIVE_INFINITY;
    private static final double[] PLANE_LO = { Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
    private static final double[] PLANE_HI = { INF, INF };

    /** correlation, sd of x, sd of y, mean of x, mean of y */
    private static final double[][] PAIRS = { { 0.0, 1.0, 1.0, 0.0, 0.0 }, { 0.5, 1.0, 1.0, 0.0, 0.0 },
            { 0.9, 1.0, 1.0, 0.0, 0.0 }, { 0.99, 1.0, 1.0, 0.0, 0.0 }, { 0.999, 1.0, 1.0, 0.0, 0.0 },
            { 0.9999, 1.0, 1.0, 0.0, 0.0 }, { 0.999, 100.0, 0.001, 7.0, -3.0 },
            { 0.9, 1.0, 1.0, 1.0e4, -1.0e4 } };

    /** The log density of a correlated normal pair, already normalized. */
    private static DMultiFunction pair(final double rho, final double sx, final double sy, final double mx,
            final double my) {
        return new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double dx = (p[0] - mx) / sx;
                double dy = (p[1] - my) / sy;
                double q = (dx * dx - 2.0 * rho * dx * dy + dy * dy) / (1.0 - rho * rho);
                return -0.5 * q - Math.log(2.0 * Math.PI * sx * sy * Math.sqrt(1.0 - rho * rho));
            }
        };
    }

    private static String name(double[] c) {
        return "rho=" + c[0] + " sd=(" + c[1] + ", " + c[2] + ") at (" + c[3] + ", " + c[4] + ")";
    }

    /** Built once: a construction is seven integrations and about 130 ms. */
    private static final BivariatePosterior[] POSTERIORS = new BivariatePosterior[PAIRS.length];

    static {
        for (int i = 0; i < PAIRS.length; i++) {
            double[] c = PAIRS[i];
            POSTERIORS[i] = BivariatePosterior.of(pair(c[0], c[1], c[2], c[3], c[4]), PLANE_LO, PLANE_HI);
        }
    }

    private static BivariatePosterior over(int i) {
        return POSTERIORS[i];
    }

    // ------------------------------------------------------------------
    // the evidence and the density
    // ------------------------------------------------------------------

    @Test
    public void anAlreadyNormalizedDensityHasNoEvidenceToFind() {
        // a density that already integrates to one leaves nothing for the
        // normalizer to do, at every correlation and at every pair of scales.
        // Worst figure measured: 1.3e-13, at rho = 0.9999
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            BivariatePosterior p = over(idx);
            assertEquals(name(c) + " is already normalized", 0.0, p.logEvidence(), 1.0e-11);
        }
    }

    @Test
    public void theCorrelatedCaseTheAxisAlignedRouteLoses() {
        // The regression for the whole class. This is the posterior the
        // caller's own axes cannot integrate: the substitution samples its mass
        // and fails to resolve it, so it returns a number with no exception and
        // no warning. If the whitening is ever removed as an optimization, the
        // second half of this test is what fails
        final double rho = 0.999;
        final double sx = 100.0;
        final double sy = 0.001;
        final double mx = 7.0;
        final double my = -3.0;
        BivariatePosterior p = BivariatePosterior.of(pair(rho, sx, sy, mx, my), PLANE_LO, PLANE_HI);
        assertEquals("the whitened route lost the evidence", 0.0, p.logEvidence(), 1.0e-11);

        final DMultiFunction law = pair(rho, sx, sy, mx, my);
        double onTheCallersAxes = Quadrature.integrate(new DBiFunction() {
            @Override
            public double apply(double x, double y) {
                return Math.exp(law.apply(new double[] { x, y }));
            }
        }, Double.NEGATIVE_INFINITY, INF, Double.NEGATIVE_INFINITY, INF, 1.0e-10, mx, my);
        assertTrue("the axis-aligned route did not fail, so this test proves nothing : " + onTheCallersAxes,
                Math.abs(onTheCallersAxes - 1.0) > 0.5);
    }

    @Test
    public void theEvidenceRecoversAConstantThatWasAddedToIt() {
        // what makes the evidence usable for model comparison, including an
        // offset large enough that the density it produces is not representable
        for (final double offset : new double[] { 12.345, -12.345, -5000.0 }) {
            final DMultiFunction base = pair(0.9, 1.0, 1.0, 0.0, 0.0);
            BivariatePosterior p = BivariatePosterior.of(new DMultiFunction() {
                @Override
                public double apply(double[] x) {
                    return base.apply(x) + offset;
                }
            }, PLANE_LO, PLANE_HI);
            assertEquals("an offset of " + offset + " was not recovered", offset, p.logEvidence(),
                    1.0e-10 * Math.max(1.0, Math.abs(offset)));
        }
    }

    @Test
    public void outsideTheBoxThereIsNoDensity() {
        BivariatePosterior p = BivariatePosterior.of(pair(0.5, 1.0, 1.0, 0.0, 0.0), new double[] { -2.0, -2.0 },
                new double[] { 2.0, 2.0 });
        double[][] outside = { { -3.0, 0.0 }, { 3.0, 0.0 }, { 0.0, -3.0 }, { 0.0, 3.0 },
                { Double.NaN, 0.0 }, { 0.0, Double.NaN } };
        for (double[] x : outside) {
            assertEquals("density at (" + x[0] + ", " + x[1] + ")", 0.0, p.pdf(x), 0.0);
            assertEquals("log density at (" + x[0] + ", " + x[1] + ")", Double.NEGATIVE_INFINITY, p.logPdf(x), 0.0);
        }
    }

    @Test
    public void theModeIsWhereTheMassIs() {
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            BivariatePosterior p = over(idx);
            double[] mode = new double[2];
            p.mode(mode);
            double scale = Math.max(1.0, Math.max(Math.abs(c[3]), Math.abs(c[4])));
            assertEquals(name(c) + ": the mode is not at the center", c[3], mode[0], 1.0e-5 * scale);
            assertEquals(name(c) + ": the mode is not at the center", c[4], mode[1], 1.0e-5 * scale);
        }
    }

    @Test
    public void theCurvatureIsTheCurvature() {
        // the difference-quotient Hessian against the exact one. It only has to
        // be good enough to precondition with, and it is far better than that
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            double rho = c[0];
            double sx = c[1];
            double sy = c[2];
            double det = 1.0 - rho * rho;
            DMatrix k = over(idx).curvature();
            assertEquals(name(c) + ": curvature(0,0)", 1.0 / (sx * sx * det), k.get(0, 0),
                    1.0e-6 * Math.abs(1.0 / (sx * sx * det)));
            assertEquals(name(c) + ": curvature(1,1)", 1.0 / (sy * sy * det), k.get(1, 1),
                    1.0e-6 * Math.abs(1.0 / (sy * sy * det)));
            assertEquals(name(c) + ": curvature(0,1)", -rho / (sx * sy * det), k.get(0, 1),
                    1.0e-6 * Math.max(1.0, Math.abs(rho / (sx * sy * det))));
            assertEquals(name(c) + ": the curvature is not symmetric", k.get(0, 1), k.get(1, 0), 0.0);
        }
    }

    @Test
    public void theLaplaceValueIsExactForANormalAndNotOtherwise() {
        // the whitening had to factor the curvature anyway, so the Laplace
        // evidence is free -- and the gap between it and the quadrature answer
        // is a diagnostic only if it really is zero where the posterior is
        // Gaussian. Measured there: about 1e-8, which is the accuracy of the
        // difference-quotient curvature and not a real difference
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            BivariatePosterior p = over(idx);
            double gap = Math.abs(p.logEvidence() - p.logEvidenceLaplace());
            assertTrue(name(c) + ": a Gaussian should have no Laplace gap, this one has " + gap, gap < 1.0e-3);
        }

        // and a posterior that is genuinely not Gaussian must show one.
        // Measured: 9.2e-03 at a hundred observations, shrinking to 4.7e-04 at
        // two thousand as the likelihood takes over from the skew -- which is
        // itself the right behaviour and the reason the small sample is used
        // here
        BivariatePosterior skew = scaleAndLocation(100);
        double gap = Math.abs(skew.logEvidence() - skew.logEvidenceLaplace());
        assertTrue("a skewed posterior should show a Laplace gap, this one shows " + gap, gap > 5.0e-3);
    }

    @Test
    public void theTwoAxisOrdersAgree() {
        // exchanging the axes gives a different Cholesky factor, a different
        // substitution and a different subdivision of the same integral, so
        // this owes nothing to the tolerance either was asked for. Worst figure
        // measured: 9.1e-15
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            double gap = over(idx).errorEstimate();
            assertTrue(name(c) + ": the two axis orders disagree by " + gap, gap < 1.0e-11);
        }
    }

    @Test
    public void aBoundedAxisIsStillIntegrated() {
        // three shapes of support, and the route differs for each: neither axis
        // bounded is the native rule over the plane, one bounded is the native
        // rule over a box in u because the second whitened coordinate depends
        // on the second parameter alone, and both bounded needs the nested form
        // whose inner limits run between two lines
        BivariatePosterior half = BivariatePosterior.of(pair(0.9, 1.0, 1.0, 0.0, 1.0),
                new double[] { Double.NEGATIVE_INFINITY, 0.0 }, PLANE_HI);
        // P(Y > 0) for a standard normal centered at one
        assertEquals("a half-bounded axis", Math.log(0.8413447460685429), half.logEvidence(), 1.0e-11);
        assertTrue("the two axis orders disagree : " + half.errorEstimate(), half.errorEstimate() < 1.0e-11);

        BivariatePosterior box = BivariatePosterior.of(pair(0.9, 1.0, 1.0, 0.0, 0.0), new double[] { -3.0, -3.0 },
                new double[] { 3.0, 3.0 });
        assertEquals("a doubly bounded axis", Math.log(0.995821216644), box.logEvidence(), 1.0e-11);
        assertTrue("the two axis orders disagree : " + box.errorEstimate(), box.errorEstimate() < 1.0e-11);

        // the bounded axis first, which the permutation has to move
        BivariatePosterior swapped = BivariatePosterior.of(pair(0.9, 1.0, 1.0, 1.0, 0.0), new double[] { 0.0,
                Double.NEGATIVE_INFINITY }, PLANE_HI);
        assertEquals("the bounded axis given first", Math.log(0.8413447460685429), swapped.logEvidence(), 1.0e-11);
    }

    @Test
    public void theWholeThingWouldBeZeroWithoutTheShift() {
        // the test the design exists for, in two parameters: the mean and the
        // logarithm of the standard deviation of two thousand observations
        int n = 2000;
        double[] data = normalSample(n, 7.0, 1.5);
        DMultiFunction logPosterior = meanAndLogScale(data);
        assertEquals("the unshifted density has not underflowed, so this test proves nothing", 0.0,
                Math.exp(logPosterior.apply(new double[] { 7.0, Math.log(1.5) })), 0.0);

        BivariatePosterior p = BivariatePosterior.of(logPosterior, PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 },
                1.0e-10);
        assertTrue("the log evidence is not finite : " + p.logEvidence(), isFinite(p.logEvidence()));

        double xbar = average(data);
        double ss = 0.0;
        for (int i = 0; i < data.length; i++) {
            ss += (data[i] - xbar) * (data[i] - xbar);
        }
        double[] mode = new double[2];
        p.mode(mode);
        assertEquals("the mode of the mean is not the sample average", xbar, mode[0], 1.0e-6);
        assertEquals("the mode of the log scale is wrong", Math.log(Math.sqrt(ss / n)), mode[1], 1.0e-6);

        double[] mean = new double[2];
        p.mean(mean);
        assertEquals("the posterior mean of the mean is not the sample average", xbar, mean[0], 1.0e-8);
        assertEquals("the posterior spread of the mean is wrong", Math.sqrt(ss / n / n),
                Math.sqrt(p.covariance().get(0, 0)), 1.0e-3 * Math.sqrt(ss / n / n));
        assertTrue("the two axis orders disagree : " + p.errorEstimate(), p.errorEstimate() < 1.0e-9);
    }

    @Test
    public void theBivariatePosteriorRejectsWhatItCannotIntegrate() {
        final DMultiFunction ok = pair(0.5, 1.0, 1.0, 0.0, 0.0);
        refuses("a null integrand", null, PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 }, 1.0e-10);
        refuses("a null lower box", ok, null, PLANE_HI, new double[] { 0.0, 0.0 }, 1.0e-10);
        refuses("a one-entry box", ok, new double[] { 0.0 }, PLANE_HI, new double[] { 0.0, 0.0 }, 1.0e-10);
        refuses("a reversed axis", ok, new double[] { 1.0, -1.0 }, new double[] { 0.0, 1.0 },
                new double[] { 0.5, 0.0 }, 1.0e-10);
        refuses("a NaN bound", ok, new double[] { Double.NaN, -1.0 }, new double[] { 1.0, 1.0 },
                new double[] { 0.0, 0.0 }, 1.0e-10);
        refuses("a null start", ok, PLANE_LO, PLANE_HI, null, 1.0e-10);
        refuses("a start outside the box", ok, new double[] { -1.0, -1.0 }, new double[] { 1.0, 1.0 },
                new double[] { 5.0, 0.0 }, 1.0e-10);
        refuses("an infinite start", ok, PLANE_LO, PLANE_HI, new double[] { INF, 0.0 }, 1.0e-10);
        refuses("a zero tolerance", ok, PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 }, 0.0);

        // and the two that are about the posterior rather than the arguments
        refuses("an improper posterior", new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                return 0.0;
            }
        }, PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 }, 1.0e-10);
        refuses("a saddle rather than a peak", new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                return -(x[0] * x[0] - x[1] * x[1]);
            }
        }, PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 }, 1.0e-10);
    }

    // ------------------------------------------------------------------
    // the moments
    // ------------------------------------------------------------------

    @Test
    public void theMomentsAreTheMomentsItWasGiven() {
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            double rho = c[0];
            double sx = c[1];
            double sy = c[2];
            BivariatePosterior p = over(idx);
            double[] mean = new double[2];
            p.mean(mean);
            double scale = Math.max(1.0, Math.max(Math.abs(c[3]), Math.abs(c[4])));
            assertEquals(name(c) + ": the mean of x", c[3], mean[0], 1.0e-8 * scale);
            assertEquals(name(c) + ": the mean of y", c[4], mean[1], 1.0e-8 * scale);

            DMatrix cov = p.covariance();
            // a second moment over a pair whose scales are 1e5 apart is the
            // hardest thing here and reaches 3.3e-7 relative; everything with
            // matched scales reaches 1e-9 or better, and asserting one loose
            // tolerance over both would stop testing the easy cases at all
            double tol = sx == sy ? 1.0e-8 : 1.0e-5;
            assertEquals(name(c) + ": var(x)", sx * sx, cov.get(0, 0), tol * sx * sx);
            assertEquals(name(c) + ": var(y)", sy * sy, cov.get(1, 1), tol * sy * sy);
            // scaled by sx*sy and not by the covariance itself, because at
            // rho == 0 the covariance is zero and what is measured there is
            // 1.8e-17, which no relative tolerance can express
            assertEquals(name(c) + ": cov(x,y)", rho * sx * sy, cov.get(0, 1), tol * sx * sy);
            assertEquals(name(c) + ": the covariance is not symmetric", cov.get(0, 1), cov.get(1, 0), 0.0);
        }
    }

    @Test
    public void theCovarianceIsAPositiveDefiniteMatrix() {
        for (int idx = 0; idx < PAIRS.length; idx++) {
            DMatrix cov = over(idx).covariance();
            // exactly symmetric by construction, so Cholesky may not refuse it
            CholeskyDecomp.cholesky(cov);
        }
    }

    @Test
    public void theCovarianceIsTheInverseCurvatureForANormalAndNotOtherwise() {
        // that identity is exactly what the Laplace approximation asserts, so
        // it has to hold where the posterior is Gaussian and fail where it is
        // not -- otherwise the comparison says nothing.
        //
        // How well it holds depends on the conditioning, and that is arithmetic
        // rather than a defect: the product amplifies the 1e-8 accuracy of the
        // difference-quotient curvature by the condition number. Measured over
        // the total deviation from the identity: 7.9e-09 at rho = 0.5, 2.4e-07
        // at 0.9, 5.7e-06 at 0.999, 3.1e-04 at 0.9999 and 1.8e-03 for the pair
        // whose scales are 1e5 apart. So the identity is asserted over the
        // moderately correlated pairs, where the separation from a genuinely
        // skewed posterior is four orders of magnitude
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            if (c[0] > 0.99 || c[1] != c[2]) {
                continue;
            }
            BivariatePosterior p = over(idx);
            DMatrix product = p.curvature().mul(p.covariance());
            assertEquals(name(c) + ": curvature times covariance is not the identity", 1.0, product.get(0, 0),
                    1.0e-5);
            assertEquals(name(c) + ": curvature times covariance is not the identity", 1.0, product.get(1, 1),
                    1.0e-5);
            assertEquals(name(c) + ": curvature times covariance is not the identity", 0.0, product.get(0, 1),
                    1.0e-5);
        }

        // and a hundred observations of an unknown mean and log scale deviate
        // by 5.1e-02, four orders of magnitude clear of the Gaussians above
        BivariatePosterior skew = scaleAndLocation(100);
        DMatrix product = skew.curvature().mul(skew.covariance());
        double offDiagonal = Math.abs(product.get(0, 1)) + Math.abs(product.get(1, 0));
        double diagonal = Math.abs(product.get(0, 0) - 1.0) + Math.abs(product.get(1, 1) - 1.0);
        assertTrue("a skewed posterior should not have the curvature invert its covariance : "
                + (diagonal + offDiagonal), diagonal + offDiagonal > 1.0e-3);
    }

    // ------------------------------------------------------------------
    // the marginal density
    // ------------------------------------------------------------------

    @Test
    public void theMarginalIsTheMarginalItWasGiven() {
        // the marginal of a correlated normal pair is a normal that does not
        // depend on the correlation at all, which makes it a sharp check
        for (int idx = 0; idx < PAIRS.length; idx++) {
            double[] c = PAIRS[idx];
            BivariatePosterior p = over(idx);
            for (int axis = 0; axis < 2; axis++) {
                double sd = c[1 + axis];
                double center = c[3 + axis];
                for (int k = -3; k <= 3; k++) {
                    double x = center + k * sd;
                    double want = Math.exp(-0.5 * k * k) / (Math.sqrt(2.0 * Math.PI) * sd);
                    assertEquals(name(c) + ": the marginal of axis " + axis + " at " + k + " sd", want,
                            p.marginalDensity(axis, x), 1.0e-8 * want);
                }
            }
        }
    }

    @Test
    public void theMarginalSurvivesAConditionalWidthTheCenterAloneCannotFind() {
        // the regression for the scaling inside marginalDensity. At these
        // scales the conditional width of the free axis is 4.5e-5, and the
        // substitution has derivative one at its center -- so naming the center
        // is not enough and the inner integral came back as exactly zero
        BivariatePosterior p = BivariatePosterior.of(pair(0.999, 100.0, 0.001, 7.0, -3.0), PLANE_LO, PLANE_HI);
        double want = 1.0 / (Math.sqrt(2.0 * Math.PI) * 100.0);
        assertEquals("the marginal of the wide axis", want, p.marginalDensity(0, 7.0), 1.0e-6 * want);
        double wantNarrow = 1.0 / (Math.sqrt(2.0 * Math.PI) * 0.001);
        assertEquals("the marginal of the narrow axis", wantNarrow, p.marginalDensity(1, -3.0),
                1.0e-6 * wantNarrow);
    }

    @Test
    public void theMarginalIsZeroOutsideItsRange() {
        BivariatePosterior p = BivariatePosterior.of(pair(0.5, 1.0, 1.0, 0.0, 0.0), new double[] { -2.0, -2.0 },
                new double[] { 2.0, 2.0 });
        assertEquals(0.0, p.marginalDensity(0, -3.0), 0.0);
        assertEquals(0.0, p.marginalDensity(0, 3.0), 0.0);
        assertEquals(0.0, p.marginalDensity(1, Double.NaN), 0.0);
    }

    @Test
    public void theMarginalRejectsAnAxisThatIsNotThere() {
        BivariatePosterior p = over(2);
        for (int i : new int[] { -1, 2, Integer.MAX_VALUE }) {
            try {
                p.marginalDensity(i, 0.0);
                fail("the marginal accepted axis " + i);
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------
    // the conjugate family, where the answer is arithmetic
    // ------------------------------------------------------------------

    @Test
    public void theNormalInverseGammaPosteriorIsTheOneTheArithmeticGives() {
        // a genuine two-parameter conjugate posterior: the mean and the
        // variance of a normal sample, both unknown. Every part of the answer
        // is closed form, and the marginals are laws this library can write
        // down -- a Student t for the mean, an inverse gamma for the variance
        final double mu0 = 5.0;
        final double kappa0 = 2.0;
        final double a0 = 3.0;
        final double b0 = 4.0;
        final double[] data = normalSample(60, 7.0, 1.5);
        final int n = data.length;
        double xbar = average(data);
        double ss = 0.0;
        for (int i = 0; i < n; i++) {
            ss += (data[i] - xbar) * (data[i] - xbar);
        }
        double kappaN = kappa0 + n;
        double muN = (kappa0 * mu0 + n * xbar) / kappaN;
        double aN = a0 + n / 2.0;
        double bN = b0 + 0.5 * ss + kappa0 * n * (xbar - mu0) * (xbar - mu0) / (2.0 * kappaN);

        BivariatePosterior p = BivariatePosterior.of(normalInverseGamma(data, mu0, kappa0, a0, b0), new double[] {
                Double.NEGATIVE_INFINITY, 0.0 }, PLANE_HI, new double[] { xbar, ss / n }, 1.0e-10);

        double[] mean = new double[2];
        p.mean(mean);
        assertEquals("the posterior mean of mu", muN, mean[0], 1.0e-8);
        assertEquals("the posterior mean of sigma squared", bN / (aN - 1.0), mean[1], 1.0e-8);

        DMatrix cov = p.covariance();
        // Var(mu) = b_n / ((a_n - 1) kappa_n), Var(sigma2) = b_n^2 / ((a_n-1)^2 (a_n-2))
        assertEquals("the posterior variance of mu", bN / ((aN - 1.0) * kappaN), cov.get(0, 0), 1.0e-9);
        assertEquals("the posterior variance of sigma squared",
                bN * bN / ((aN - 1.0) * (aN - 1.0) * (aN - 2.0)), cov.get(1, 1), 1.0e-7);

        // the evidence in closed form
        double expected = -0.5 * n * Math.log(2.0 * Math.PI) + 0.5 * Math.log(kappa0 / kappaN)
                + a0 * Math.log(b0) - aN * Math.log(bN) + FastGamma.logGamma(aN) - FastGamma.logGamma(a0);
        assertEquals("the marginal likelihood of the data", expected, p.logEvidence(),
                1.0e-9 * Math.abs(expected));

        // the marginal of mu is a scaled and shifted Student t on 2 a_n
        StudentT t = new StudentT(2.0 * aN);
        double tScale = Math.sqrt(bN / (aN * kappaN));
        for (double z : new double[] { -2.0, -1.0, 0.0, 1.0, 2.0 }) {
            double at = muN + z * tScale;
            assertEquals("the marginal of mu at " + z, Math.exp(t.logPdf(z)) / tScale, p.marginalDensity(0, at),
                    1.0e-6 * Math.exp(t.logPdf(z)) / tScale);
        }

        // and the marginal of sigma squared is an inverse gamma
        InverseGamma ig = new InverseGamma(aN, bN);
        for (double q : new double[] { 0.2, 0.5, 0.8 }) {
            double at = ig.inverseCdf(q);
            assertEquals("the marginal of sigma squared at " + q, ig.pdf(at), p.marginalDensity(1, at),
                    1.0e-6 * ig.pdf(at));
        }
    }

    // ------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------

    /** The Normal-InverseGamma posterior of a normal sample, up to a constant. */
    private static DMultiFunction normalInverseGamma(final double[] data, final double mu0, final double kappa0,
            final double a0, final double b0) {
        return new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double mu = p[0];
                double sigma2 = p[1];
                if (!(sigma2 > 0.0)) {
                    return Double.NEGATIVE_INFINITY;
                }
                double value = -(a0 + 1.0) * Math.log(sigma2) - b0 / sigma2;
                value += -0.5 * Math.log(sigma2) - kappa0 * (mu - mu0) * (mu - mu0) / (2.0 * sigma2);
                double sum = 0.0;
                for (int i = 0; i < data.length; i++) {
                    double d = data[i] - mu;
                    sum += d * d;
                }
                value += -0.5 * data.length * Math.log(2.0 * Math.PI * sigma2) - sum / (2.0 * sigma2);
                // the two prior normalizing constants, so the evidence comes
                // out as the marginal likelihood and not as it times a factor
                value += a0 * Math.log(b0) - FastGamma.logGamma(a0);
                value += 0.5 * Math.log(kappa0 / (2.0 * Math.PI));
                return value;
            }
        };
    }

    /** The mean and the log standard deviation of a normal sample, flat priors. */
    private static DMultiFunction meanAndLogScale(final double[] data) {
        return new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double mu = p[0];
                double logSigma = p[1];
                double sigma = Math.exp(logSigma);
                double sum = 0.0;
                for (int i = 0; i < data.length; i++) {
                    double d = data[i] - mu;
                    sum += d * d;
                }
                return -data.length * logSigma - sum / (2.0 * sigma * sigma);
            }
        };
    }

    /** A posterior that is genuinely not Gaussian: the log scale is skew. */
    private static BivariatePosterior scaleAndLocation(int n) {
        double[] data = normalSample(n, 7.0, 1.5);
        return BivariatePosterior.of(meanAndLogScale(data), PLANE_LO, PLANE_HI, new double[] { 0.0, 0.0 },
                1.0e-10);
    }

    private static double average(double[] data) {
        double sum = 0.0;
        for (int i = 0; i < data.length; i++) {
            sum += data[i];
        }
        return sum / data.length;
    }

    /** Box-Muller over the LCG the tests in this repository use. */
    private static double[] normalSample(int n, double mean, double sd) {
        long lcg = 20260827L;
        double[] data = new double[n];
        for (int i = 0; i < n; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double u1 = ((lcg >>> 11) + 0.5) * 0x1.0p-53;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double u2 = ((lcg >>> 11) + 0.5) * 0x1.0p-53;
            data[i] = mean + sd * Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        }
        return data;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static void refuses(String what, DMultiFunction f, double[] lo, double[] hi, double[] start,
            double epsTol) {
        try {
            BivariatePosterior.of(f, lo, hi, start, epsTol);
            fail("the factory accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
