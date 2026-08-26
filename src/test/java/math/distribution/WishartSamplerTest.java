package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;
import math.stats.Alternative;
import math.stats.HypothesisTests;

/**
 * The Wishart sampler, checked against the two moments of the law it claims and
 * then, where the scale matrix is diagonal, against the chi-squared its
 * diagonal must be.
 */
public final class WishartSamplerTest {

    /** A scale matrix with the same value down the diagonal and off it. */
    private static DMatrix symmetric(int p, double diagonal, double offDiagonal) {
        DMatrix m = new DMatrix(p, p);
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                m.set(i, j, i == j ? diagonal : offDiagonal);
            }
        }
        return m;
    }

    @Test
    public void testTheMeanIsTheScaleTimesTheDegreesOfFreedom() {
        // E[W] = nu Sigma is exact, and it is the first thing a Bartlett
        // decomposition with the degrees of freedom counted wrongly breaks
        int draws = 100000;
        for (double[] shape : new double[][] { { 2, 5 }, { 3, 10 }, { 4, 4 } }) {
            int p = (int) shape[0];
            double nu = shape[1];
            DMatrix scale = symmetric(p, 2.0, 0.5);
            WishartSampler sampler = WishartSampler.of(scale, nu);
            assertEquals("dimension", p, sampler.dimension());
            assertEquals("degrees of freedom", nu, sampler.degreesOfFreedom(), 0.0);

            PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
            double[][] total = new double[p][p];
            DMatrix out = new DMatrix(p, p);
            for (int d = 0; d < draws; d++) {
                sampler.sample(prng, out);
                for (int i = 0; i < p; i++) {
                    for (int j = 0; j < p; j++) {
                        total[i][j] += out.getUnsafe(i, j);
                    }
                }
            }
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    // measured over 200000 draws: 0.0295 at worst
                    assertEquals("p=" + p + " nu=" + nu + " at (" + i + "," + j + ")",
                            nu * scale.getUnsafe(i, j), total[i][j] / draws, 0.1);
                }
            }
        }
    }

    @Test
    public void testTheVarianceIsTheOneTheLawPrescribes() {
        // Var(W_ij) = nu (Sigma_ij^2 + Sigma_ii Sigma_jj), which the mean
        // cannot see and which catches a diagonal whose degrees of freedom are
        // off by one
        int p = 3;
        double nu = 8.0;
        int draws = 200000;
        DMatrix scale = symmetric(p, 2.0, 0.5);
        WishartSampler sampler = WishartSampler.of(scale, nu);
        PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
        double[][] sum = new double[p][p];
        double[][] sumOfSquares = new double[p][p];
        DMatrix out = new DMatrix(p, p);
        for (int d = 0; d < draws; d++) {
            sampler.sample(prng, out);
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    double v = out.getUnsafe(i, j);
                    sum[i][j] += v;
                    sumOfSquares[i][j] += v * v;
                }
            }
        }
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                double mean = sum[i][j] / draws;
                double variance = sumOfSquares[i][j] / draws - mean * mean;
                double wanted = nu * (scale.getUnsafe(i, j) * scale.getUnsafe(i, j)
                        + scale.getUnsafe(i, i) * scale.getUnsafe(j, j));
                // measured over 500000 draws: 0.4 percent at worst
                assertEquals("at (" + i + "," + j + ")", wanted, variance, 0.03 * wanted);
            }
        }
    }

    @Test
    public void testADiagonalScaleMakesTheDiagonalChiSquared() {
        // with no correlation to carry, W_ii / Sigma_ii is exactly a
        // chi-squared with nu degrees of freedom -- a continuous check, not a
        // binned one
        DMatrix scale = new DMatrix(3, 3);
        scale.set(0, 0, 1.0);
        scale.set(1, 1, 4.0);
        scale.set(2, 2, 0.25);
        for (double nu : new double[] { 3.0, 10.0 }) {
            WishartSampler sampler = WishartSampler.of(scale, nu);
            ChiSquare law = new ChiSquare(nu);
            for (int component = 0; component < 3; component++) {
                int reps = 80;
                int draws = 4000;
                double[] pValues = new double[reps];
                for (int r = 0; r < reps; r++) {
                    PseudoRandom prng = DefaultRng.newPseudoRandom(r * 7919L + 1);
                    double[] drawn = new double[draws];
                    DMatrix out = new DMatrix(3, 3);
                    for (int d = 0; d < draws; d++) {
                        sampler.sample(prng, out);
                        drawn[d] = out.getUnsafe(component, component) / scale.getUnsafe(component, component);
                    }
                    pValues[r] = HypothesisTests.kolmogorovSmirnov(drawn, law, Alternative.TWO_SIDED).pValue;
                }
                double uniformity = HypothesisTests
                        .kolmogorovSmirnov(pValues, new Uniform(0.0, 1.0), Alternative.TWO_SIDED).pValue;
                // measured: 0.1645 .. 0.4951
                assertTrue("nu=" + nu + " component " + component + ": uniformity p = " + uniformity,
                        uniformity > 0.001);
            }
        }
    }

    @Test
    public void testEveryDrawIsSymmetricAndPositiveDefinite() {
        // the second is an independent statement about the result rather than
        // about the algorithm: a Cholesky of the draw either succeeds or does
        // not
        for (int p : new int[] { 2, 3, 5 }) {
            DMatrix scale = symmetric(p, 3.0, 0.7);
            WishartSampler sampler = WishartSampler.of(scale, p + 2.0);
            PseudoRandom prng = DefaultRng.newPseudoRandom(7L);
            for (int d = 0; d < 1000; d++) {
                DMatrix drawn = sampler.sample(prng);
                for (int i = 0; i < p; i++) {
                    for (int j = 0; j < i; j++) {
                        // measured over 6000 draws: the two halves agree bit
                        // for bit, since both come out of the same product
                        assertEquals("p=" + p + " draw " + d + " at (" + i + "," + j + ")",
                                drawn.getUnsafe(i, j), drawn.getUnsafe(j, i), 0.0);
                    }
                }
                CholeskyDecomp.cholesky(drawn);
            }
        }
    }

    @Test
    public void testOneDimensionIsAScaledChiSquared() {
        DMatrix scale = new DMatrix(1, 1);
        scale.set(0, 0, 3.0);
        WishartSampler sampler = WishartSampler.of(scale, 6.0);
        PseudoRandom prng = DefaultRng.newPseudoRandom(9L);
        double[] drawn = new double[200000];
        DMatrix out = new DMatrix(1, 1);
        for (int d = 0; d < drawn.length; d++) {
            sampler.sample(prng, out);
            drawn[d] = out.getUnsafe(0, 0) / 3.0;
        }
        double pValue = HypothesisTests
                .kolmogorovSmirnov(drawn, new ChiSquare(6.0), Alternative.TWO_SIDED).pValue;
        // measured: 0.7109
        assertTrue("a one by one Wishart is not a scaled chi-squared, p = " + pValue, pValue > 0.001);
    }

    @Test
    public void testTheScaleMatrixIsUsedAsGiven() {
        // the multivariate normal sampler in this package perturbs its
        // covariance by 1e-7 times the identity before factorizing it. This one
        // does not, and a scale matrix that only just fails to be positive
        // definite is refused rather than nudged
        DMatrix singular = new DMatrix(2, 2);
        singular.set(0, 0, 1.0);
        singular.set(0, 1, 1.0);
        singular.set(1, 0, 1.0);
        singular.set(1, 1, 1.0);
        try {
            WishartSampler.of(singular, 5.0);
            fail("of accepted a singular scale matrix");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }

        // and a scale matrix that is positive definite is reproduced in the
        // mean without a shift: a perturbation of 1e-7 on a scale of 1e-5 would
        // be one percent
        DMatrix tiny = new DMatrix(2, 2);
        tiny.set(0, 0, 1.0e-5);
        tiny.set(1, 1, 1.0e-5);
        WishartSampler sampler = WishartSampler.of(tiny, 4.0);
        PseudoRandom prng = DefaultRng.newPseudoRandom(13L);
        DMatrix out = new DMatrix(2, 2);
        double total = 0.0;
        int draws = 100000;
        for (int d = 0; d < draws; d++) {
            sampler.sample(prng, out);
            total += out.getUnsafe(0, 0);
        }
        assertEquals("the mean of a tiny scale", 4.0 * 1.0e-5, total / draws, 0.02 * 4.0 * 1.0e-5);
    }

    @Test
    public void testTheSamplerIsReproducibleFromTheSeed() {
        WishartSampler sampler = WishartSampler.of(symmetric(3, 2.0, 0.4), 7.0);
        PseudoRandom one = DefaultRng.newPseudoRandom(23L);
        PseudoRandom two = DefaultRng.newPseudoRandom(23L);
        DMatrix first = new DMatrix(3, 3);
        DMatrix again = new DMatrix(3, 3);
        for (int d = 0; d < 200; d++) {
            sampler.sample(one, first);
            sampler.sample(two, again);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    assertEquals("draw " + d + " at (" + i + "," + j + ")", first.getUnsafe(i, j),
                            again.getUnsafe(i, j), 0.0);
                }
            }
        }
    }

    @Test
    public void testTheSamplerRejectsWhatItCannotDraw() {
        DMatrix good = symmetric(3, 2.0, 0.3);
        rejectsBuild("a null scale matrix", null, 5.0);
        rejectsBuild("a non-square scale matrix", new DMatrix(2, 3), 5.0);
        rejectsBuild("degrees of freedom at the dimension less one", good, 2.0);
        rejectsBuild("degrees of freedom below the dimension less one", good, 1.0);
        rejectsBuild("degrees of freedom that are not a number", good, Double.NaN);

        DMatrix asymmetric = new DMatrix(2, 2);
        asymmetric.set(0, 0, 2.0);
        asymmetric.set(0, 1, 1.0);
        asymmetric.set(1, 0, 0.5);
        asymmetric.set(1, 1, 2.0);
        rejectsBuild("an asymmetric scale matrix", asymmetric, 5.0);

        DMatrix indefinite = new DMatrix(2, 2);
        indefinite.set(0, 0, 1.0);
        indefinite.set(0, 1, 2.0);
        indefinite.set(1, 0, 2.0);
        indefinite.set(1, 1, 1.0);
        rejectsBuild("an indefinite scale matrix", indefinite, 5.0);

        WishartSampler sampler = WishartSampler.of(good, 5.0);
        try {
            sampler.sample(null, new DMatrix(3, 3));
            fail("sample accepted a null generator");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            sampler.sample(DefaultRng.newPseudoRandom(1L), (DMatrix) null);
            fail("sample accepted a null output matrix");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            sampler.sample(DefaultRng.newPseudoRandom(1L), new DMatrix(2, 2));
            fail("sample accepted an output matrix of the wrong shape");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsBuild(String what, DMatrix scale, double nu) {
        try {
            WishartSampler.of(scale, nu);
            fail("of accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
