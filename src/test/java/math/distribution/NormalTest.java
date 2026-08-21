package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.RootFinder;

/**
 * The density is checked against the closed form and against its own
 * distribution function, which is what nothing did before: `Normal` had no
 * caller anywhere in the test tree, and the normalizing constant was
 * {@code 1 / (variance * sqrt(2 pi))} where its own javadoc said
 * {@code 1 / (stdDev * sqrt(2 pi))}. The density was therefore wrong by a
 * factor of {@code 1 / sigma} everywhere except at {@code sigma == 1}, where
 * the two agree and a casual check would have seen nothing.
 */
public class NormalTest {

    private static final double[] SIGMAS = { 0.079011, 0.5, 1.0, 2.0, 10.0 };
    private static final double[] MEANS = { 0.0, 3.0, 299.8524 };

    private static DFunction density(final Normal n) {
        return new DFunction() {
            @Override
            public double apply(double x) {
                return n.pdf(x);
            }
        };
    }

    /** A density integrates to one. That is what makes it a density. */
    @Test
    public void densityIntegratesToOne() {
        for (int m = 0; m < MEANS.length; ++m) {
            for (int s = 0; s < SIGMAS.length; ++s) {
                Normal normal = new Normal(MEANS[m], SIGMAS[s]);
                double mass = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, density(normal),
                        MEANS[m] - 12.0 * SIGMAS[s], MEANS[m] + 12.0 * SIGMAS[s], 1.0e-13, 40);
                assertEquals("mean " + MEANS[m] + ", sigma " + SIGMAS[s], 1.0, mass, 1.0e-9);
            }
        }
    }

    /** The peak of the density is one over sigma root two pi, not one over the variance. */
    @Test
    public void densityMatchesTheClosedForm() {
        for (int s = 0; s < SIGMAS.length; ++s) {
            double sigma = SIGMAS[s];
            Normal normal = new Normal(3.0, sigma);
            for (int k = -30; k <= 30; ++k) {
                double x = 3.0 + 0.1 * k * sigma;
                double z = (x - 3.0) / sigma;
                double expected = Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2.0 * Math.PI));
                assertEquals("sigma " + sigma + " at z = " + (0.1 * k), expected, normal.pdf(x),
                        1.0e-13 * Math.max(1.0, expected));
            }
        }
    }

    /** The distribution function is the integral of the density it belongs to. */
    @Test
    public void cdfIsTheIntegralOfThePdf() {
        for (int s = 0; s < SIGMAS.length; ++s) {
            double sigma = SIGMAS[s];
            Normal normal = new Normal(3.0, sigma);
            double from = 3.0 - 12.0 * sigma;
            double below = normal.cdf(from);
            for (int k = -30; k <= 30; k += 3) {
                double x = 3.0 + 0.1 * k * sigma;
                double integrated = below + AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, density(normal),
                        from, x, 1.0e-13, 40);
                assertEquals("sigma " + sigma + " at z = " + (0.1 * k), normal.cdf(x), integrated, 1.0e-9);
            }
        }
    }

    /** The closed-form quantile against an independent root of the same equation. */
    @Test
    public void inverseCdfAgreesWithAnIndependentRoot() {
        final Normal normal = new Normal(3.0, 2.0);
        double[] ps = { 0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999 };
        for (int k = 0; k < ps.length; ++k) {
            final double p = ps[k];
            double byRoot = RootFinder.brentDekker(3.0 - 40.0, 3.0 + 40.0, new DFunction() {
                @Override
                public double apply(double x) {
                    return normal.cdf(x) - p;
                }
            }, 1.0e-14);
            assertEquals("p = " + p, byRoot, normal.inverseCdf(p), 1.0e-9);
            assertEquals("round trip at p = " + p, p, normal.cdf(normal.inverseCdf(p)), 1.0e-12);
        }
    }

    /** Mean and variance are the first two moments of the density. */
    @Test
    public void momentsAreTheIntegralsTheyClaimToBe() {
        final double mu = 3.0;
        final double sigma = 2.0;
        final Normal normal = new Normal(mu, sigma);
        double first = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                return x * normal.pdf(x);
            }
        }, mu - 14.0 * sigma, mu + 14.0 * sigma, 1.0e-13, 40);
        double second = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                return (x - mu) * (x - mu) * normal.pdf(x);
            }
        }, mu - 14.0 * sigma, mu + 14.0 * sigma, 1.0e-13, 40);

        assertEquals(normal.mean(), first, 1.0e-8);
        assertEquals(normal.variance(), second, 1.0e-8);
        assertEquals(sigma * sigma, normal.variance(), 0.0);
    }

    /** The no-argument constructor is the standard normal. */
    @Test
    public void theDefaultIsStandard() {
        Normal standard = new Normal();
        assertEquals(0.0, standard.mean(), 0.0);
        assertEquals(1.0, standard.variance(), 0.0);
        assertEquals(0.3989422804014327, standard.pdf(0.0), 1.0e-15);
        assertEquals(0.5, standard.cdf(0.0), 1.0e-15);
        assertEquals(0.8413447460685429, standard.cdf(1.0), 1.0e-12);
    }

    /** A standard deviation of zero or less is not a distribution. */
    @Test
    public void aNonPositiveScaleIsRejected() {
        for (int k = 0; k < 2; ++k) {
            try {
                new Normal(0.0, k == 0 ? 0.0 : -1.0);
                org.junit.Assert.fail("expected IllegalArgumentException");
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage().contains("Standard deviation"));
            }
        }
    }
}
