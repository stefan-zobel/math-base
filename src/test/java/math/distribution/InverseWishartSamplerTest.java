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
 * The inverse-Wishart sampler. The one thing that can go wrong in a wrapper
 * this thin is inverting the scale matrix in the wrong place, or not at all,
 * and the two tests that matter here are the two that would catch it.
 */
public final class InverseWishartSamplerTest {

    private static DMatrix symmetric(int p, double diagonal, double offDiagonal) {
        DMatrix m = new DMatrix(p, p);
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                m.set(i, j, i == j ? diagonal : offDiagonal);
            }
        }
        return m;
    }

    private static void symmetrize(DMatrix m) {
        for (int i = 0; i < m.numRows(); i++) {
            for (int j = 0; j < i; j++) {
                double average = 0.5 * (m.getUnsafe(i, j) + m.getUnsafe(j, i));
                m.setUnsafe(i, j, average);
                m.setUnsafe(j, i, average);
            }
        }
    }

    @Test
    public void testTheMeanIsTheScaleOverTheDegreesOfFreedomLessTheDimension() {
        // E[W] = Psi / (nu - p - 1). A sampler that forgot to invert the scale
        // matrix would land here on Psi's inverse instead, which is not a
        // subtle difference
        int draws = 150000;
        for (double[] shape : new double[][] { { 2, 6 }, { 3, 9 }, { 2, 12 } }) {
            int p = (int) shape[0];
            double nu = shape[1];
            DMatrix psi = symmetric(p, 2.0, 0.5);
            InverseWishartSampler sampler = InverseWishartSampler.of(psi, nu);
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
            double denominator = nu - p - 1.0;
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    double wanted = psi.getUnsafe(i, j) / denominator;
                    // measured over 300000 draws: 0.58 percent at worst
                    assertEquals("p=" + p + " nu=" + nu + " at (" + i + "," + j + ")", wanted,
                            total[i][j] / draws, 0.03 * Math.abs(wanted));
                }
            }
        }
    }

    @Test
    public void testInvertingADrawGivesTheWishartWithTheInvertedScale() {
        // the definition, taken the other way round. From the same seeded
        // generator the two have to be the same numbers, not merely the same
        // distribution -- which is what makes the wrapper honest rather than
        // plausible
        int p = 3;
        double nu = 8.0;
        DMatrix psi = symmetric(p, 2.0, 0.5);
        InverseWishartSampler inverse = InverseWishartSampler.of(psi, nu);
        DMatrix invertedScale = psi.inverse();
        symmetrize(invertedScale);
        WishartSampler direct = WishartSampler.of(invertedScale, nu);

        PseudoRandom one = DefaultRng.newPseudoRandom(17L);
        PseudoRandom two = DefaultRng.newPseudoRandom(17L);
        double worst = 0.0;
        for (int d = 0; d < 1000; d++) {
            DMatrix back = inverse.sample(one).inverse();
            DMatrix straight = direct.sample(two);
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    double scale = Math.max(Math.abs(straight.getUnsafe(i, j)), 1.0e-9);
                    worst = Math.max(worst,
                            Math.abs(back.getUnsafe(i, j) - straight.getUnsafe(i, j)) / scale);
                }
            }
        }
        // measured over 2000 draws: 1.4e-12, which is two matrix inversions
        // worth of round-off and nothing else
        assertTrue("the round trip differs by " + worst, worst < 1.0e-9);
    }

    @Test
    public void testOneDimensionIsAnInverseGamma() {
        DMatrix psi = new DMatrix(1, 1);
        psi.set(0, 0, 3.0);
        for (double nu : new double[] { 5.0, 9.0 }) {
            InverseWishartSampler sampler = InverseWishartSampler.of(psi, nu);
            PseudoRandom prng = DefaultRng.newPseudoRandom(9L);
            int draws = 200000;
            DMatrix out = new DMatrix(1, 1);
            double total = 0.0;
            double[] reciprocal = new double[draws];
            for (int d = 0; d < draws; d++) {
                sampler.sample(prng, out);
                total += out.getUnsafe(0, 0);
                // the reciprocal of an inverse-Wishart draw scaled by psi is
                // the chi-squared the Wishart behind it was built from
                reciprocal[d] = 3.0 / out.getUnsafe(0, 0);
            }
            // measured: 0.99982 against 1.0, and 0.42867 against 0.42857
            assertEquals("nu=" + nu + ": the mean", 3.0 / (nu - 2.0), total / draws,
                    0.02 * 3.0 / (nu - 2.0));
            double pValue = HypothesisTests
                    .kolmogorovSmirnov(reciprocal, new ChiSquare(nu), Alternative.TWO_SIDED).pValue;
            // measured: 0.0870 and 0.8532
            assertTrue("nu=" + nu + ": the reciprocal is not chi-squared, p = " + pValue, pValue > 0.001);
        }
    }

    @Test
    public void testEveryDrawIsSymmetricAndPositiveDefinite() {
        for (int p : new int[] { 2, 3, 5 }) {
            InverseWishartSampler sampler = InverseWishartSampler.of(symmetric(p, 3.0, 0.7), p + 3.0);
            PseudoRandom prng = DefaultRng.newPseudoRandom(7L);
            for (int d = 0; d < 600; d++) {
                DMatrix drawn = sampler.sample(prng);
                for (int i = 0; i < p; i++) {
                    for (int j = 0; j < i; j++) {
                        // the inversion behind a draw is a solve, which does
                        // not come back exactly symmetric; the sampler averages
                        // the two halves, so this is exact after all
                        assertEquals("p=" + p + " draw " + d + " at (" + i + "," + j + ")",
                                drawn.getUnsafe(i, j), drawn.getUnsafe(j, i), 0.0);
                    }
                }
                CholeskyDecomp.cholesky(drawn);
            }
        }
    }

    @Test
    public void testTooFewDegreesOfFreedomStillDraw() {
        // the mean exists only for nu > p + 1, but the distribution is defined
        // for nu > p - 1 and there is no reason to refuse it
        InverseWishartSampler sampler = InverseWishartSampler.of(symmetric(3, 1.0, 0.0), 3.0);
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        for (int d = 0; d < 200; d++) {
            DMatrix drawn = sampler.sample(prng);
            double trace = 0.0;
            for (int i = 0; i < 3; i++) {
                trace += drawn.getUnsafe(i, i);
            }
            assertTrue("a draw at nu = p has trace " + trace, trace > 0.0 && !Double.isNaN(trace));
        }
    }

    @Test
    public void testTheSamplerIsReproducibleFromTheSeed() {
        InverseWishartSampler sampler = InverseWishartSampler.of(symmetric(3, 2.0, 0.4), 8.0);
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

        DMatrix singular = new DMatrix(2, 2);
        singular.set(0, 0, 1.0);
        singular.set(0, 1, 1.0);
        singular.set(1, 0, 1.0);
        singular.set(1, 1, 1.0);
        rejectsBuild("a singular scale matrix", singular, 5.0);

        InverseWishartSampler sampler = InverseWishartSampler.of(good, 6.0);
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
            InverseWishartSampler.of(scale, nu);
            fail("of accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
