package math.stats.bayes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.Beta;
import math.distribution.Gamma;
import math.distribution.LogNormal;
import math.distribution.MultivariateNormal;
import math.distribution.Normal;
import math.distribution.StudentT;
import math.fun.DFunction;
import math.fun.DMultiFunction;
import math.linalg.DMatrix;

/**
 * {@link LaplaceApproximation}: the Gaussian that matches a posterior at its
 * mode.
 * <p>
 * Two kinds of check. Against an equicorrelated normal in {@code d} dimensions
 * the approximation is <em>exact</em> and every part of the answer has a closed
 * form, so the measurement is of the machinery. Against
 * {@link ScalarPosterior} and {@link BivariatePosterior} the approximation is
 * <em>wrong</em> by an amount only quadrature can supply, and that is the point:
 * one and two parameters are where the exact answer is affordable, so they are
 * what licenses the class in three and above.
 */
public final class LaplaceApproximationTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static double[] fill(int d, double v) {
        double[] a = new double[d];
        java.util.Arrays.fill(a, v);
        return a;
    }

    /**
     * An equicorrelated normal in {@code d} dimensions, already normalized, so
     * its log evidence is zero and its covariance is one on the diagonal and
     * {@code r} off it. Written through the closed form of the inverse of an
     * equicorrelation matrix.
     */
    private static DMultiFunction equicorrelated(final int d, final double r) {
        final double det = Math.pow(1.0 - r, d - 1) * (1.0 + (d - 1) * r);
        final double a = (1.0 + (d - 2) * r) / ((1.0 - r) * (1.0 + (d - 1) * r));
        final double b = -r / ((1.0 - r) * (1.0 + (d - 1) * r));
        final double[] center = new double[d];
        for (int i = 0; i < d; i++) {
            center[i] = (i + 1.0) / d;
        }
        return new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double q = 0.0;
                double s = 0.0;
                for (int i = 0; i < d; i++) {
                    double z = p[i] - center[i];
                    q += z * z;
                    s += z;
                }
                return -0.5 * ((a - b) * q + b * s * s) - 0.5 * (d * Math.log(2.0 * Math.PI) + Math.log(det));
            }
        };
    }

    // ------------------------------------------------------------------
    // where the approximation is exact
    // ------------------------------------------------------------------

    @Test
    public void aGaussianIsItsOwnLaplaceApproximation() {
        // measured worst figures: the log evidence -7.4e-08 at d=3 growing to
        // -1.8e-06 at d=20, the mode 3.1e-08 throughout, the covariance 2.8e-06
        for (double r : new double[] { 0.0, 0.5, 0.9 }) {
            for (int d : new int[] { 1, 2, 3, 5, 8, 12, 20 }) {
                String what = "d=" + d + " r=" + r;
                LaplaceApproximation l = LaplaceApproximation.of(equicorrelated(d, r), fill(d, -INF),
                        fill(d, INF));
                assertEquals(what + ": dimension", d, l.dimension());
                assertEquals(what + ": the evidence of a normalized law", 0.0, l.logEvidence(), 1.0e-4);

                double[] mode = new double[d];
                l.mode(mode);
                for (int i = 0; i < d; i++) {
                    assertEquals(what + ": the mode at " + i, (i + 1.0) / d, mode[i], 1.0e-6);
                }

                // the approximation of a Gaussian is that Gaussian
                MultivariateNormal n = l.asNormal();
                DMatrix cov = n.covariance();
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        assertEquals(what + ": the covariance at (" + i + ", " + j + ")", i == j ? 1.0 : r,
                                cov.get(i, j), 1.0e-4);
                    }
                }
                assertEquals(what + ": nothing should fall outside an unbounded box", 0.0,
                        l.massOutsideSupport(), 0.0);
                assertTrue(what + ": the mode search did not converge", l.modeConverged());
            }
        }
    }

    @Test
    public void theDensityOfTheApproximationIsThePosteriorItself() {
        // for a Gaussian posterior the two agree everywhere, not only at the
        // mode, which is a stronger statement than matching the evidence
        int d = 4;
        double r = 0.5;
        DMultiFunction law = equicorrelated(d, r);
        LaplaceApproximation l = LaplaceApproximation.of(law, fill(d, -INF), fill(d, INF));
        MultivariateNormal n = l.asNormal();
        double[] x = new double[d];
        for (int t = -2; t <= 2; t++) {
            for (int i = 0; i < d; i++) {
                x[i] = (i + 1.0) / d + t * 0.6;
            }
            assertEquals("the approximation is not the law at t=" + t, law.apply(x.clone()), n.logPdf(x), 1.0e-4);
        }
    }

    @Test
    public void theCurvatureIsTheCurvature() {
        // the difference-quotient Hessian against the closed-form precision
        for (double r : new double[] { 0.0, 0.5, 0.9 }) {
            for (int d : new int[] { 2, 5, 12 }) {
                double a = (1.0 + (d - 2) * r) / ((1.0 - r) * (1.0 + (d - 1) * r));
                double b = -r / ((1.0 - r) * (1.0 + (d - 1) * r));
                DMatrix k = LaplaceApproximation.of(equicorrelated(d, r), fill(d, -INF), fill(d, INF))
                        .curvature();
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        double want = i == j ? a : b;
                        assertEquals("d=" + d + " r=" + r + " at (" + i + ", " + j + ")", want, k.get(i, j),
                                1.0e-5 * Math.max(1.0, Math.abs(want)));
                        assertEquals("the curvature is not symmetric", k.get(i, j), k.get(j, i), 0.0);
                    }
                }
            }
        }
    }

    @Test
    public void aPosteriorOverManyObservationsIsStillFound() {
        // three parameters over five thousand observations, where the
        // unnormalized density is zero at every point a rule would evaluate
        final double[] data = normalSample(5000, 7.0, 1.5);
        DMultiFunction posterior = new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                // mean, log scale, and a third parameter that only shifts the
                // whole thing, so the answer stays known
                double mu = p[0];
                double sigma = Math.exp(p[1]);
                double nuisance = p[2];
                double sum = 0.0;
                for (int i = 0; i < data.length; i++) {
                    double e = data[i] - mu;
                    sum += e * e;
                }
                return -data.length * p[1] - sum / (2.0 * sigma * sigma) - 0.5 * nuisance * nuisance;
            }
        };
        assertEquals("the unshifted density has not underflowed, so this test proves nothing", 0.0,
                Math.exp(posterior.apply(new double[] { 7.0, Math.log(1.5), 0.0 })), 0.0);

        LaplaceApproximation l = LaplaceApproximation.of(posterior, fill(3, -INF), fill(3, INF),
                new double[] { 0.0, 0.0, 0.0 });
        assertTrue("the log evidence is not finite : " + l.logEvidence(), isFinite(l.logEvidence()));

        double xbar = average(data);
        double ss = 0.0;
        for (int i = 0; i < data.length; i++) {
            ss += (data[i] - xbar) * (data[i] - xbar);
        }
        double[] mode = new double[3];
        l.mode(mode);
        assertEquals("the mode of the mean", xbar, mode[0], 1.0e-5);
        assertEquals("the mode of the log scale", Math.log(Math.sqrt(ss / data.length)), mode[1], 1.0e-5);
        assertEquals("the mode of the nuisance parameter", 0.0, mode[2], 1.0e-5);
        // the posterior spread of the mean is sigma / sqrt(n)
        assertEquals("the spread of the mean", Math.sqrt(ss / data.length / data.length),
                Math.sqrt(l.asNormal().covariance().get(0, 0)), 1.0e-3 * Math.sqrt(ss / data.length / data.length));
    }

    @Test
    public void theLaplaceApproximationRejectsWhatItCannotApproximate() {
        final DMultiFunction ok = equicorrelated(2, 0.0);
        refuses("a null integrand", null, fill(2, -INF), fill(2, INF), fill(2, 0.0));
        refuses("a null box", ok, null, fill(2, INF), fill(2, 0.0));
        refuses("boxes of different length", ok, fill(2, -INF), fill(3, INF), fill(2, 0.0));
        refuses("an empty box", ok, new double[0], new double[0], new double[0]);
        refuses("a reversed axis", ok, new double[] { 1.0, -1.0 }, new double[] { 0.0, 1.0 },
                new double[] { 0.5, 0.0 });
        refuses("a NaN bound", ok, new double[] { Double.NaN, -1.0 }, fill(2, 1.0), fill(2, 0.0));
        refuses("a null start", ok, fill(2, -INF), fill(2, INF), null);
        refuses("a start outside the box", ok, fill(2, -1.0), fill(2, 1.0), new double[] { 5.0, 0.0 });
        refuses("an infinite start", ok, fill(2, -INF), fill(2, INF), new double[] { INF, 0.0 });

        refuses("a saddle rather than a peak", new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                return -(x[0] * x[0] - x[1] * x[1]);
            }
        }, fill(2, -INF), fill(2, INF), fill(2, 0.0));
        refuses("an improper posterior", new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                return 0.0;
            }
        }, fill(2, -INF), fill(2, INF), fill(2, 0.0));

        LaplaceApproximation l = LaplaceApproximation.of(ok, fill(2, -INF), fill(2, INF));
        for (double[] bad : new double[][] { null, new double[1], new double[3] }) {
            try {
                l.mode(bad);
                fail("mode accepted an output of the wrong length");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------
    // the mass a Gaussian puts where the posterior has none
    // ------------------------------------------------------------------

    @Test
    public void theMassOutsideIsTheTailItReallyIs() {
        // exact in one dimension, which is what makes it checkable at all:
        // measured 0.420740 against 0.420740 for a mode a fifth of a standard
        // deviation from the boundary
        for (final double mu : new double[] { 0.2, 1.0, 3.0, 8.0 }) {
            LaplaceApproximation l = LaplaceApproximation.of(new DMultiFunction() {
                @Override
                public double apply(double[] p) {
                    double z = p[0] - mu;
                    return -0.5 * z * z;
                }
            }, new double[] { 0.0 }, new double[] { INF }, new double[] { mu });
            double[] mode = new double[1];
            l.mode(mode);
            double sd = Math.sqrt(l.asNormal().covariance().get(0, 0));
            double exact = new Normal(mode[0], sd).cdf(0.0);
            assertEquals("mode near " + mu + ": the reported leak", exact, l.massOutsideSupport(), 1.0e-12);
        }
    }

    @Test
    public void aPosteriorWellInsideItsBoxLeaksNothing() {
        LaplaceApproximation l = LaplaceApproximation.of(equicorrelated(3, 0.5), fill(3, -50.0), fill(3, 50.0));
        assertTrue("a mode fifty widths from the wall leaks " + l.massOutsideSupport(),
                l.massOutsideSupport() < 1.0e-12);
    }

    @Test
    public void theLeakIsAUnionBoundAndSoNeverBelowTheTruth() {
        // two axes bounded at once: the axis tails are summed, so the reported
        // number counts a corner twice and can only be too large
        LaplaceApproximation l = LaplaceApproximation.of(equicorrelated(2, 0.0), new double[] { 0.0, 0.0 },
                fill(2, INF), new double[] { 0.5, 1.0 });
        double[] mode = new double[2];
        l.mode(mode);
        DMatrix cov = l.asNormal().covariance();
        double first = new Normal(mode[0], Math.sqrt(cov.get(0, 0))).cdf(0.0);
        double second = new Normal(mode[1], Math.sqrt(cov.get(1, 1))).cdf(0.0);
        assertEquals("the leak is not the sum of the axis tails", first + second, l.massOutsideSupport(),
                1.0e-12);
        // and the true mass outside the quadrant is at most that sum
        assertTrue("a union bound cannot be below either of its terms",
                l.massOutsideSupport() >= Math.max(first, second));
        assertTrue("the leak left [0, 1]", l.massOutsideSupport() >= 0.0 && l.massOutsideSupport() <= 1.0);
    }

    // ------------------------------------------------------------------
    // where the approximation is wrong, and by how much
    // ------------------------------------------------------------------

    @Test
    public void theApproximationAgreesWithQuadratureInOneDimension() {
        // The honest statement of when it holds, measured rather than asserted:
        //
        //   Normal(3, 2)      8.6e-09   exact, as it must be
        //   Gamma(3, 2)       4.1e-02
        //   Beta(2, 5)        9.7e-02
        //   StudentT(5)       1.4e-01
        //   StudentT(2)       3.2e-01
        //   LogNormal(0, 1)   5.0e-01   exactly sigma^2 / 2
        //
        // Each is pinned within a factor of two, so the test fails both if the
        // class gets worse and if these figures were ever wrong
        assertLaplaceGap("Normal(3, 2)", new Normal(3.0, 2.0)::logPdf, -INF, INF, 3.0, 0.0, 1.0e-6);
        assertLaplaceGap("Gamma(3, 2)", new Gamma(3.0, 2.0)::logPdf, 0.0, INF, 4.0, 4.13e-2, 4.0e-3);
        assertLaplaceGap("Beta(2, 5)", new Beta(2.0, 5.0)::logPdf, 0.0, 1.0, 0.2, 9.71e-2, 1.0e-2);
        assertLaplaceGap("StudentT(5)", new StudentT(5.0)::logPdf, -INF, INF, 0.0, 1.408e-1, 1.0e-2);
        assertLaplaceGap("StudentT(2)", new StudentT(2.0)::logPdf, -INF, INF, 0.0, 3.235e-1, 3.0e-2);
        // a log-normal is a normal on the logarithm, and the Laplace error in
        // the log evidence is then exactly sigma^2 / 2 -- a closed form, so it
        // is asserted tightly
        assertLaplaceGap("LogNormal(0, 1)", new LogNormal(0.0, 1.0)::logPdf, 0.0, INF, Math.exp(-1.0), 0.5,
                1.0e-6);
    }

    @Test
    public void theErrorShrinksWithTheSampleSize() {
        // and this is what licenses the class where quadrature cannot follow:
        // the approximation is not merely wrong by some amount, it is wrong by
        // an amount that falls like 1/n. Measured: 9.2e-03, 2.3e-03, 4.7e-04
        double[] gaps = new double[3];
        int[] sizes = { 100, 400, 2000 };
        for (int k = 0; k < sizes.length; k++) {
            DMultiFunction posterior = meanAndLogScale(normalSample(sizes[k], 7.0, 1.5));
            BivariatePosterior exact = BivariatePosterior.of(posterior, fill(2, -INF), fill(2, INF),
                    new double[] { 0.0, 0.0 }, 1.0e-10);
            LaplaceApproximation lap = LaplaceApproximation.of(posterior, fill(2, -INF), fill(2, INF),
                    new double[] { 0.0, 0.0 });
            gaps[k] = Math.abs(exact.logEvidence() - lap.logEvidence());
        }
        assertEquals("the gap at a hundred observations", 9.22e-3, gaps[0], 1.0e-3);
        assertEquals("the gap at four hundred", 2.30e-3, gaps[1], 3.0e-4);
        assertEquals("the gap at two thousand", 4.69e-4, gaps[2], 1.0e-4);
        assertTrue("the error does not fall with the sample size", gaps[0] > gaps[1] && gaps[1] > gaps[2]);
        // four times the observations, about four times less error
        assertEquals("the error does not fall like 1/n", 4.0, gaps[0] / gaps[1], 0.6);
        assertEquals("the error does not fall like 1/n", 5.0, gaps[1] / gaps[2], 1.0);
    }

    // ------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------

    private static void assertLaplaceGap(String what, final DFunction law, double lo, double hi, double start,
            double expected, double tolerance) {
        ScalarPosterior exact = ScalarPosterior.of(law, lo, hi, start, 1.0e-10);
        LaplaceApproximation lap = LaplaceApproximation.of(new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                return law.apply(p[0]);
            }
        }, new double[] { lo }, new double[] { hi }, new double[] { start });
        double gap = Math.abs(exact.logEvidence() - lap.logEvidence());
        assertEquals(what + ": the Laplace error against quadrature", expected, gap, tolerance);
    }

    private static DMultiFunction meanAndLogScale(final double[] data) {
        return new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double mu = p[0];
                double sigma = Math.exp(p[1]);
                double sum = 0.0;
                for (int i = 0; i < data.length; i++) {
                    double e = data[i] - mu;
                    sum += e * e;
                }
                return -data.length * p[1] - sum / (2.0 * sigma * sigma);
            }
        };
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

    private static void refuses(String what, DMultiFunction f, double[] lo, double[] hi, double[] start) {
        try {
            LaplaceApproximation.of(f, lo, hi, start);
            fail("the factory accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
