package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.linalg.PosDefiniteMatrixGenerator;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;
import math.solve.Quadrature;

/**
 * {@link MultivariateNormal}, the density half of a pair whose sampler was
 * built first.
 * <p>
 * Checked against oracles that share no line of code: {@link Normal} in one
 * dimension and coordinate by coordinate, {@link Quadrature} over the plane for
 * the normalizing constant, and {@link ChiSquare} for the law of the squared
 * Mahalanobis distance. Tolerances were measured first -- the worst figures
 * seen were {@code 8.9e-16} relative on the log density against {@code Normal},
 * {@code 2.2e-16} against a product of marginals and {@code 3.8e-12} on the
 * integral over the plane.
 */
public final class MultivariateNormalTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static DMatrix scalar(double variance) {
        DMatrix s = new DMatrix(1, 1);
        s.set(0, 0, variance);
        return s;
    }

    private static DMatrix pair(double sx, double sy, double rho) {
        DMatrix s = new DMatrix(2, 2);
        s.set(0, 0, sx * sx);
        s.set(1, 1, sy * sy);
        s.set(0, 1, rho * sx * sy);
        s.set(1, 0, rho * sx * sy);
        return s;
    }

    // ------------------------------------------------------------------
    // the density
    // ------------------------------------------------------------------

    @Test
    public void atOneDimensionItIsANormal() {
        // not bit for bit: the Cholesky factor of a 1x1 is sqrt(sd*sd), which
        // is not sd in the last bit. It matched exactly in 140 of 272 points
        // and the worst relative gap over the rest was 8.9e-16
        for (double mu : new double[] { 0.0, 3.0, -7.5, 1.0e4 }) {
            for (double sd : new double[] { 1.0, 0.25, 12.0, 1.0e-3 }) {
                MultivariateNormal m = new MultivariateNormal(new double[] { mu }, scalar(sd * sd));
                Normal n = new Normal(mu, sd);
                assertEquals("dimension", 1, m.dimension());
                for (double z = -4.0; z <= 4.0; z += 0.5) {
                    double x = mu + z * sd;
                    assertEquals("mu=" + mu + " sd=" + sd + ": the log density at " + z + " sd", n.logPdf(x),
                            m.logPdf(new double[] { x }), 1.0e-13 * Math.max(1.0, Math.abs(n.logPdf(x))));
                    assertEquals("mu=" + mu + " sd=" + sd + ": the density at " + z + " sd", n.pdf(x),
                            m.pdf(new double[] { x }), 1.0e-13 * Math.max(1.0e-300, n.pdf(x)));
                }
                assertEquals("the log determinant of a 1x1", Math.log(sd * sd), m.logDeterminant(),
                        1.0e-13 * Math.max(1.0, Math.abs(Math.log(sd * sd))));
            }
        }
    }

    @Test
    public void aDiagonalCovarianceIsAProductOfMarginals() {
        // independence written out: with no correlation the log density is the
        // sum of the one dimensional ones, which Normal supplies
        for (int d : new int[] { 2, 5, 20 }) {
            double[] mu = new double[d];
            DMatrix s = new DMatrix(d, d);
            for (int i = 0; i < d; i++) {
                mu[i] = i - d / 2.0;
                s.set(i, i, (i + 1) * 0.5);
            }
            MultivariateNormal m = new MultivariateNormal(mu, s);
            double[] x = new double[d];
            double want = 0.0;
            for (int i = 0; i < d; i++) {
                x[i] = mu[i] + 0.7 * (i % 3 - 1);
                want += new Normal(mu[i], Math.sqrt(s.get(i, i))).logPdf(x[i]);
            }
            assertEquals("d=" + d + ": the log density is not the sum of the marginals", want, m.logPdf(x),
                    1.0e-13 * Math.abs(want));
        }
    }

    @Test
    public void theDensityIntegratesToOne() {
        // an independent route to the normalizing constant, and the only check
        // here that does not go through a closed form
        for (double rho : new double[] { 0.0, 0.5, 0.9, -0.7 }) {
            final MultivariateNormal m = new MultivariateNormal(new double[] { 1.0, -2.0 }, pair(1.0, 2.0, rho));
            double value = Quadrature.integrate((x, y) -> m.pdf(new double[] { x, y }), -INF, INF, -INF, INF,
                    1.0e-10, 1.0, -2.0);
            assertEquals("rho=" + rho + ": the density does not integrate to one", 1.0, value, 1.0e-9);
        }
    }

    @Test
    public void theMarginalIsTheMarginalItWasGiven() {
        for (double rho : new double[] { 0.0, 0.5, 0.9 }) {
            MultivariateNormal m = new MultivariateNormal(new double[] { 1.0, -2.0 }, pair(1.0, 2.0, rho));
            Normal first = m.marginal(0);
            Normal second = m.marginal(1);
            assertEquals("the first marginal mean", 1.0, first.mean(), 0.0);
            assertEquals("the first marginal variance", 1.0, first.variance(), 1.0e-14);
            assertEquals("the second marginal mean", -2.0, second.mean(), 0.0);
            assertEquals("the second marginal variance", 4.0, second.variance(), 1.0e-14);
            // the marginal of a multivariate normal does not depend on the
            // correlation, which is what makes this a sharp statement
            for (double x = -6.0; x <= 8.0; x += 1.0) {
                assertEquals("the first marginal density at " + x, new Normal(1.0, 1.0).pdf(x), first.pdf(x),
                        1.0e-14);
            }
        }
    }

    @Test
    public void theLogDensityAnswersWhereTheDensityHasUnderflowed() {
        // fifty coordinates, two and a half units out in each: the density is
        // exactly zero and the logarithm is an ordinary number
        int d = 50;
        MultivariateNormal m = new MultivariateNormal(new double[d],
                PosDefiniteMatrixGenerator.generate(d, 3L));
        double[] x = new double[d];
        for (int i = 0; i < d; i++) {
            x[i] = 2.5;
        }
        assertEquals("the density should have underflowed", 0.0, m.pdf(x), 0.0);
        double logPdf = m.logPdf(x);
        assertTrue("the log density is not finite : " + logPdf, !Double.isInfinite(logPdf)
                && !Double.isNaN(logPdf));
        assertTrue("the log density is not below the underflow threshold : " + logPdf, logPdf < -745.0);

        // and where both are defined the two routes agree
        double[] near = new double[d];
        assertEquals("the mass and its logarithm disagree at the mean", Math.log(m.pdf(near)), m.logPdf(near),
                1.0e-12);
    }

    @Test
    public void aCovarianceFarFromTheOriginIsStillANormal() {
        // The regression for the symmetrization. CholeskyDecomp compares the
        // two halves with an ABSOLUTE tolerance of 1e-10, so a covariance whose
        // entries are of order 1e6 -- symmetric to within its own last bit --
        // is refused outright. This class compares relatively and averages
        // within that, so it accepts what its own factorizer will not
        int d = 6;
        DMatrix big = PosDefiniteMatrixGenerator.generate(d, 7L).scaleInplace(1.0e6);
        big.set(1, 0, Math.nextAfter(Math.nextAfter(big.get(0, 1), INF), INF));

        boolean refused = false;
        try {
            CholeskyDecomp.cholesky(big);
        } catch (IllegalArgumentException expected) {
            refused = true;
        }
        assertTrue("CholeskyDecomp accepted it, so this test proves nothing", refused);

        MultivariateNormal m = new MultivariateNormal(new double[d], big);
        assertTrue("the log determinant is not finite : " + m.logDeterminant(),
                !Double.isInfinite(m.logDeterminant()) && !Double.isNaN(m.logDeterminant()));
        assertTrue("the density at the mean is not positive", m.pdf(new double[d]) > 0.0);
        // and what it hands back is exactly symmetric whatever it was given
        DMatrix used = m.covariance();
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                assertEquals("the covariance it uses is not symmetric at (" + i + ", " + j + ")", used.get(i, j),
                        used.get(j, i), 0.0);
            }
        }
    }

    @Test
    public void theInputIsCopiedRatherThanKept() {
        double[] mu = { 1.0, -2.0 };
        DMatrix s = pair(1.0, 2.0, 0.5);
        MultivariateNormal m = new MultivariateNormal(mu, s);
        double before = m.logPdf(new double[] { 0.0, 0.0 });
        mu[0] = 99.0;
        s.set(0, 0, 99.0);
        assertEquals("the density followed the caller's arrays", before, m.logPdf(new double[] { 0.0, 0.0 }), 0.0);
        m.covariance().set(0, 0, 99.0);
        assertEquals("the covariance it hands out is not a copy", before, m.logPdf(new double[] { 0.0, 0.0 }), 0.0);
    }

    @Test
    public void theMultivariateNormalRejectsWhatItCannotBe() {
        refuses("a null mean", null, pair(1.0, 1.0, 0.0));
        refuses("a null covariance", new double[] { 0.0, 0.0 }, null);
        // DMatrix refuses a zero-size matrix of its own accord, so the empty
        // mean is offered with a covariance that is otherwise valid
        refuses("an empty mean", new double[0], scalar(1.0));
        refuses("a non-finite mean", new double[] { 0.0, Double.NaN }, pair(1.0, 1.0, 0.0));
        refuses("a non-square covariance", new double[] { 0.0, 0.0 }, new DMatrix(2, 3));
        refuses("a covariance of the wrong order", new double[] { 0.0, 0.0, 0.0 }, pair(1.0, 1.0, 0.0));
        refuses("a genuinely asymmetric covariance", new double[] { 0.0, 0.0 }, asymmetric());
        refuses("a singular covariance", new double[] { 0.0, 0.0 }, pair(1.0, 1.0, 1.0));
        refuses("an indefinite covariance", new double[] { 0.0, 0.0 }, pair(1.0, 1.0, 2.0));
        refuses("a zero covariance", new double[] { 0.0, 0.0 }, new DMatrix(2, 2));

        MultivariateNormal m = new MultivariateNormal(new double[] { 0.0, 0.0 }, pair(1.0, 1.0, 0.0));
        for (double[] bad : new double[][] { null, new double[1], new double[3] }) {
            try {
                m.logPdf(bad);
                fail("the density accepted a point of the wrong length");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
        try {
            m.marginal(2);
            fail("the marginal accepted coordinate 2");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without a message", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------------------
    // the squared Mahalanobis distance and the draw
    // ------------------------------------------------------------------

    @Test
    public void theSquaredDistanceIsTheQuadraticForm() {
        // against the inverse, computed the way the class deliberately does not
        for (double rho : new double[] { 0.0, 0.5, 0.9 }) {
            DMatrix s = pair(1.0, 2.0, rho);
            MultivariateNormal m = new MultivariateNormal(new double[] { 1.0, -2.0 }, s);
            DMatrix inv = s.inverse();
            for (double x = -2.0; x <= 4.0; x += 1.5) {
                for (double y = -5.0; y <= 1.0; y += 1.5) {
                    double dx = x - 1.0;
                    double dy = y + 2.0;
                    double want = dx * dx * inv.get(0, 0) + 2.0 * dx * dy * inv.get(0, 1)
                            + dy * dy * inv.get(1, 1);
                    assertEquals("rho=" + rho + " at (" + x + ", " + y + ")", want,
                            m.squaredMahalanobisDistance(new double[] { x, y }), 1.0e-11 * Math.max(1.0, want));
                }
            }
            assertEquals("the distance at the mean", 0.0, m.squaredMahalanobisDistance(new double[] { 1.0, -2.0 }),
                    1.0e-25);
        }
    }

    @Test
    public void theDrawIsTheFactorTimesAStandardNormal() {
        // the construction written out in the test, from the same seed
        int d = 4;
        DMatrix s = PosDefiniteMatrixGenerator.generate(d, 23L);
        double[] mu = { 1.0, -1.0, 0.5, 3.0 };
        MultivariateNormal m = new MultivariateNormal(mu, s);
        DMatrix l = CholeskyDecomp.cholesky(m.covariance());

        PseudoRandom one = DefaultRng.newPseudoRandom(99L);
        PseudoRandom two = DefaultRng.newPseudoRandom(99L);
        double[] out = new double[d];
        for (int t = 0; t < 200; t++) {
            m.sample(one, out);
            double[] z = new double[d];
            for (int i = 0; i < d; i++) {
                z[i] = two.nextGaussian();
            }
            for (int i = 0; i < d; i++) {
                double want = mu[i];
                for (int k = 0; k <= i; k++) {
                    want += l.get(i, k) * z[k];
                }
                assertEquals("draw " + t + " coordinate " + i, want, out[i], 0.0);
            }
        }
    }

    @Test
    public void theDrawHasTheMomentsAndTheShellItShould() {
        // 200000 draws: the worst mean gap measured was 0.015 and the worst
        // covariance gap 0.10, and the fraction beyond the 95 per cent shell
        // came out 0.0501 and 0.0503 -- which is the chi-squared law of the
        // squared distance, confirmed to three digits
        for (int d : new int[] { 2, 5 }) {
            DMatrix s = PosDefiniteMatrixGenerator.generate(d, 13L);
            double[] mu = new double[d];
            for (int i = 0; i < d; i++) {
                mu[i] = i + 1.0;
            }
            MultivariateNormal m = new MultivariateNormal(mu, s);
            ChiSquare chi = new ChiSquare(d);
            PseudoRandom prng = DefaultRng.newPseudoRandom(20260827L);

            int n = 200000;
            double[] out = new double[d];
            double[] sum = new double[d];
            double[][] cross = new double[d][d];
            int beyond = 0;
            for (int t = 0; t < n; t++) {
                m.sample(prng, out);
                if (chi.cdf(m.squaredMahalanobisDistance(out)) > 0.95) {
                    beyond++;
                }
                for (int i = 0; i < d; i++) {
                    sum[i] += out[i];
                    for (int j = 0; j < d; j++) {
                        cross[i][j] += out[i] * out[j];
                    }
                }
            }
            for (int i = 0; i < d; i++) {
                assertEquals("d=" + d + ": the drawn mean of coordinate " + i, mu[i], sum[i] / n, 0.05);
                for (int j = 0; j < d; j++) {
                    double c = cross[i][j] / n - (sum[i] / n) * (sum[j] / n);
                    assertEquals("d=" + d + ": the drawn covariance at (" + i + ", " + j + ")", s.get(i, j), c,
                            0.3);
                }
            }
            assertEquals("d=" + d + ": the 95 per cent shell does not hold 5 per cent", 0.05, (double) beyond / n,
                    0.005);
        }
    }

    @Test
    public void theDrawRejectsWhatItCannotDrawInto() {
        MultivariateNormal m = new MultivariateNormal(new double[] { 0.0, 0.0 }, pair(1.0, 1.0, 0.0));
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        try {
            m.sample(null, new double[2]);
            fail("the draw accepted a null generator");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without a message", expected.getMessage() != null);
        }
        for (double[] bad : new double[][] { null, new double[1], new double[3] }) {
            try {
                m.sample(prng, bad);
                fail("the draw accepted an output of the wrong length");
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------

    private static DMatrix asymmetric() {
        DMatrix s = new DMatrix(2, 2);
        s.set(0, 0, 2.0);
        s.set(1, 1, 3.0);
        s.set(0, 1, 0.5);
        s.set(1, 0, -0.5);
        return s;
    }

    private static void refuses(String what, double[] mean, DMatrix covariance) {
        try {
            new MultivariateNormal(mean, covariance);
            fail("the constructor accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
