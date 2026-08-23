package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.linalg.DMatrix;

/**
 * {@link MultivariateNormalSampler} used to draw from a static
 * {@code Lcg64Xor1024Mix} shared by every instance, which made it both unsafe
 * to use from several threads -- the gaussian cache in {@code AbstractRng64}
 * answers {@code NaN} when two threads race on it -- and impossible to
 * reproduce, since there was no way to say which stream to draw from. It owns
 * its generator now and takes a seed.
 * <p>
 * The tolerances below were calibrated over 24 seeds at 200000 samples: the
 * worst deviation seen was 0.0123 on a mean and 0.0134 on a covariance entry.
 */
public class MultivariateNormalSamplerTest {

    private static final double[] MU = { 1.0, -2.0, 0.5 };

    private static final double[][] COV = { { 2.0, 0.3, -0.4 }, { 0.3, 1.0, 0.2 }, { -0.4, 0.2, 1.5 } };

    private static DMatrix mean() {
        DMatrix m = new DMatrix(3, 1);
        for (int i = 0; i < 3; ++i) {
            m.set(i, 0, MU[i]);
        }
        return m;
    }

    private static DMatrix covariance() {
        DMatrix c = new DMatrix(3, 3);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                c.set(i, j, COV[i][j]);
            }
        }
        return c;
    }

    /** The same seed has to reproduce the same samples to the last bit. */
    @Test
    public void testTheSameSeedReproducesTheSameSamples() {
        DMatrix a = new MultivariateNormalSampler(mean(), covariance(), 42L).sample(500);
        DMatrix b = new MultivariateNormalSampler(mean(), covariance(), 42L).sample(500);
        for (int i = 0; i < a.numRows(); ++i) {
            for (int j = 0; j < a.numColumns(); ++j) {
                assertEquals("sample(" + i + ", " + j + ")", a.get(i, j), b.get(i, j), 0.0);
            }
        }
    }

    /** And a different seed must not. */
    @Test
    public void testADifferentSeedGivesDifferentSamples() {
        DMatrix a = new MultivariateNormalSampler(mean(), covariance(), 42L).sample(500);
        DMatrix b = new MultivariateNormalSampler(mean(), covariance(), 43L).sample(500);
        int equal = 0;
        for (int j = 0; j < a.numColumns(); ++j) {
            for (int i = 0; i < a.numRows(); ++i) {
                if (a.get(i, j) == b.get(i, j)) {
                    ++equal;
                }
            }
        }
        assertEquals("the two seeds produced identical samples", 0, equal);
    }

    /**
     * The distribution itself: empirical mean and covariance have to converge
     * on what was asked for. This is the guard that giving the sampler its own
     * generator did not change what it samples.
     */
    @Test
    public void testTheEmpiricalMomentsMatchTheParameters() {
        int n = 200000;
        DMatrix x = new MultivariateNormalSampler(mean(), covariance(), 20260823L).sample(n);
        double[] m = new double[3];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < 3; ++i) {
                m[i] += x.get(i, j);
            }
        }
        for (int i = 0; i < 3; ++i) {
            m[i] /= n;
            assertEquals("mean[" + i + "]", MU[i], m[i], 0.04);
        }
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double s = 0.0;
                for (int j = 0; j < n; ++j) {
                    s += (x.get(a, j) - m[a]) * (x.get(b, j) - m[b]);
                }
                assertEquals("cov[" + a + "][" + b + "]", COV[a][b], s / (n - 1), 0.05);
            }
        }
    }

    /** The sampled matrix has the shape the contract promises. */
    @Test
    public void testTheSampleHasTheRequestedShape() {
        DMatrix x = new MultivariateNormalSampler(mean(), covariance(), 7L).sample(17);
        assertEquals("rows", 3, x.numRows());
        assertEquals("columns", 17, x.numColumns());
        for (int i = 0; i < x.numRows(); ++i) {
            for (int j = 0; j < x.numColumns(); ++j) {
                assertTrue("sample(" + i + ", " + j + ") is not finite",
                        !Double.isNaN(x.get(i, j)) && !Double.isInfinite(x.get(i, j)));
            }
        }
    }

    /** The argument checks are on the seeded constructor as well. */
    @Test
    public void testTheArgumentChecksStillHold() {
        try {
            new MultivariateNormalSampler(new DMatrix(2, 1), covariance(), 1L);
            org.junit.Assert.fail("inconsistent dimensions were accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
        try {
            new MultivariateNormalSampler(new DMatrix(3, 2), covariance(), 1L);
            org.junit.Assert.fail("a mean that is not a column vector was accepted");
        } catch (IllegalArgumentException expected) {
            // as specified
        }
    }
}
