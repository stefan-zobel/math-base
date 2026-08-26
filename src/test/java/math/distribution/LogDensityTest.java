package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * The log density against the logarithm of the density, and then against the
 * places where there is no density left to take the logarithm of -- which is
 * the whole reason the method exists.
 */
public class LogDensityTest {

    /** Every continuous implementation, at settings that separate the branches. */
    private static ContinuousDistribution[] continuous() {
        return new ContinuousDistribution[] { new Normal(0.0, 1.0), new Normal(-3.0, 7.0), new Normal(0.0, 1.0e-9),
                new Normal(0.0, 1.0e9), new Exponential(1.0), new Exponential(1.0e6), new Cauchy(0.0, 1.0),
                new Cauchy(2.0, 0.5), new Uniform(0.0, 1.0), new Uniform(-3.0, 8.0), new LogNormal(0.0, 1.0),
                new LogNormal(2.0, 0.25), new Weibull(1.0, 2.0), new Weibull(3.0, 0.5), new Weibull(1.0, 1.0),
                new Gamma(2.0, 1.0), new Gamma(0.5, 3.0), new Gamma(1.0, 2.0), new ChiSquare(3.0),
                new ChiSquare(1.0), new Beta(2.0, 3.0), new Beta(0.5, 0.5), new Beta(1.0, 1.0),
                new Beta(400.0, 400.0), new StudentT(1.0), new StudentT(3.0), new StudentT(1000.0),
                new FisherF(1.0, 1.0), new FisherF(3.0, 4.0), new FisherF(10.0, 20.0), new InverseGamma(2.5, 3.0),
                new InverseGamma(0.5, 0.25), new InverseGamma(1.0, 1.0) };
    }

    /** Every discrete implementation, likewise. */
    private static DiscreteDistribution[] discrete() {
        return new DiscreteDistribution[] { new Binomial(20, 0.3), new Binomial(2000, 0.5), new Binomial(5, 1.0),
                new Binomial(5, 0.0), new Poisson(1.0), new Poisson(1000.0), new Geometric(0.5),
                new Geometric(1.0e-6), new Geometric(1.0), new NegativeBinomial(5, 0.5),
                new NegativeBinomial(1, 0.9), new Hypergeometric(50, 20, 10),
                new Hypergeometric(4000, 2000, 2000), new DiscreteUniform(0, 9), new DiscreteUniform(-5, 5) };
    }

    /**
     * The five continuous and four discrete classes whose density was already
     * being formed as the exponential of a logarithm before this slice, and
     * whose {@code pdf} now evaluates that logarithm and exponentiates it. They
     * have to return exactly the bits they returned before, which is what the
     * bit-for-bit assertion below is for; the rest add a method and cannot
     * change anything.
     */
    private static boolean delegates(Object d) {
        return d instanceof Beta || d instanceof InverseGamma || d instanceof Binomial
                || d instanceof NegativeBinomial || d instanceof Poisson || d instanceof Hypergeometric;
    }

    /** A grid inside the support, taken from the quantiles. */
    private static double[] grid(ContinuousDistribution d) {
        double[] p = { 1.0e-6, 1.0e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1.0 - 1.0e-6 };
        double[] x = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            x[i] = d.inverseCdf(p[i]);
        }
        return x;
    }

    // ------------------------------------------- the same function, in logs --

    @Test
    public void theLogDensityIsTheLogarithmOfTheDensity() {
        // wherever the density is comfortably inside the double range, taking
        // its logarithm and computing the logarithm directly have to agree.
        // Denormal densities are left out on purpose: there the logarithm of
        // the density has already lost most of its bits and the closed form is
        // the more accurate of the two, so an agreement test would be asserting
        // that the better answer matches the worse one
        double worst = 0.0;
        String at = "";
        for (ContinuousDistribution d : continuous()) {
            for (double x : grid(d)) {
                double pdf = d.pdf(x);
                if (!(pdf > 1.0e-300) || Double.isInfinite(pdf)) {
                    continue;
                }
                double wanted = Math.log(pdf);
                double got = d.logPdf(x);
                double relative = Math.abs(got - wanted) / Math.max(1.0, Math.abs(wanted));
                if (relative > worst) {
                    worst = relative;
                    at = d.getClass().getSimpleName() + " at x = " + x;
                }
            }
        }
        // measured over 33 settings: 4.4e-16, at FisherF(3, 4)
        assertTrue("the two differ by " + worst + " relative, " + at, worst < 1.0e-13);
    }

    @Test
    public void theLogMassIsTheLogarithmOfTheMass() {
        double worst = 0.0;
        String at = "";
        for (DiscreteDistribution d : discrete()) {
            long lo = d.supportLowerBound();
            long hi = d.supportUpperBound();
            for (int t = 0; t <= 40; t++) {
                long k = lo + Math.round((hi - lo) * (t / 40.0));
                if (k < Integer.MIN_VALUE || k > Integer.MAX_VALUE) {
                    continue;
                }
                double pmf = d.pmf((int) k);
                if (!(pmf > 1.0e-300)) {
                    continue;
                }
                double relative = Math.abs(d.logPmf((int) k) - Math.log(pmf)) / Math.max(1.0, Math.abs(Math.log(pmf)));
                if (relative > worst) {
                    worst = relative;
                    at = d.getClass().getSimpleName() + " at k = " + k;
                }
            }
        }
        // measured over 15 settings: 1.9e-16, at DiscreteUniform(0, 9)
        assertTrue("the two differ by " + worst + " relative, " + at, worst < 1.0e-13);
    }

    @Test
    public void theDensityThatDelegatesIsUnchangedBitForBit() {
        // the refactor these classes went through has to be invisible: pdf now
        // exponentiates logPdf, and every branch that used to stand in it --
        // outside the support, the two poles, the normalizing factor at a shape
        // of exactly one -- has to come back out of the exponential unchanged
        int checked = 0;
        for (ContinuousDistribution d : continuous()) {
            if (!delegates(d)) {
                continue;
            }
            for (double x : new double[] { Double.NEGATIVE_INFINITY, -1.0, 0.0, Double.MIN_VALUE, 1.0e-30, 0.25,
                    0.5, 0.75, 1.0, 2.0, 1.0e30, Double.MAX_VALUE, Double.POSITIVE_INFINITY }) {
                assertEquals(d.getClass().getSimpleName() + " at x = " + x, d.pdf(x), Math.exp(d.logPdf(x)), 0.0);
                checked++;
            }
        }
        for (DiscreteDistribution d : discrete()) {
            if (!delegates(d)) {
                continue;
            }
            int lo = d.supportLowerBound();
            for (int k = lo - 2; k < lo + 80; k++) {
                assertEquals(d.getClass().getSimpleName() + " at k = " + k, d.pmf(k), Math.exp(d.logPmf(k)), 0.0);
                checked++;
            }
        }
        assertTrue("nothing was checked", checked > 400);
    }

    // ------------------------------------------------ where the density ends --

    @Test
    public void theLogDensityAnswersWhereTheDensityHasUnderflowed() {
        // the point of the method, one row per family. Every argument below was
        // measured: the density is exactly zero there and the logarithm is an
        // ordinary negative number
        far(new Normal(0.0, 1.0), 39.0, -761.0, -762.0);
        far(new Normal(-3.0, 7.0), 322.0, -1080.0, -1081.5);
        far(new Exponential(1.0), 827.0, -826.5, -827.5);
        far(new Cauchy(0.0, 1.0), 8.7e153, -710.0, -710.5);
        far(new LogNormal(0.0, 1.0), 2.6e16, -752.5, -753.5);
        far(new Weibull(1.0, 2.0), 33.0, -1084.0, -1085.5);
        far(new Gamma(2.0, 1.0), 912.0, -905.0, -905.5);
        far(new ChiSquare(3.0), 1670.0, -832.0, -832.5);
        far(new Beta(400.0, 400.0), 1.0e-4, -3118.0, -3119.5);
        far(new StudentT(3.0), 1.35e81, -746.0, -746.5);
        far(new InverseGamma(2.5, 3.0), 6.0e92, -745.0, -745.5);

        // and not every family has such a point. A StudentT(1) and a FisherF
        // both form x * x (or d1 * x) before the logarithm, so the density and
        // its logarithm run out of range at the same argument; a Uniform has a
        // constant density that cannot underflow at all. Those are covered by
        // the agreement test above rather than here
        assertEquals("a uniform density does not underflow", -Math.log(1.0),
                new Uniform(0.0, 1.0).logPdf(0.5), 0.0);
    }

    /** The density is zero at {@code x} and the log density lies in the band. */
    private static void far(ContinuousDistribution d, double x, double upper, double lower) {
        String at = d.getClass().getSimpleName() + " at x = " + x;
        assertEquals(at + ": the density has not underflowed yet", 0.0, d.pdf(x), 0.0);
        double log = d.logPdf(x);
        assertTrue(at + ": the log density is " + log + ", not in (" + lower + ", " + upper + ")",
                log < upper && log > lower);
        assertEquals(at + ": the logarithm of the density is not what it should be", Double.NEGATIVE_INFINITY,
                Math.log(d.pdf(x)), 0.0);
    }

    @Test
    public void theLogMassAnswersWhereTheMassHasUnderflowed() {
        // two thousand fair coin flips landing on no heads at all is the
        // clearest of these: the mass is 2^-2000, which is not a double, and
        // the logarithm is 2000 log(0.5) exactly
        Binomial coins = new Binomial(2000, 0.5);
        assertEquals("the mass has not underflowed", 0.0, coins.pmf(0), 0.0);
        assertEquals("the log mass", 2000.0 * Math.log(0.5), coins.logPmf(0), 1.0e-9);
        assertEquals("2000 log(0.5)", -1386.2943611198906, coins.logPmf(0), 1.0e-9);

        gone(new Poisson(1.0), 178, -748.5, -749.5);
        gone(new Geometric(0.5), 1074, -745.0, -745.5);
        gone(new NegativeBinomial(5, 0.5), 1106, -745.0, -745.5);
        gone(new Hypergeometric(4000, 2000, 2000), 0, -2768.0, -2769.0);
    }

    /** The mass is zero at {@code k} and the log mass lies in the band. */
    private static void gone(DiscreteDistribution d, int k, double upper, double lower) {
        String at = d.getClass().getSimpleName() + " at k = " + k;
        assertEquals(at + ": the mass has not underflowed yet", 0.0, d.pmf(k), 0.0);
        double log = d.logPmf(k);
        assertTrue(at + ": the log mass is " + log + ", not in (" + lower + ", " + upper + ")",
                log < upper && log > lower);
    }

    @Test
    public void aLikelihoodSurvivesInLogsWhereTheProductCannot() {
        // the reason the whole slice exists, in one comparison. A thousand
        // standard normal draws, their densities multiplied and their
        // logarithms added
        Normal standard = new Normal(0.0, 1.0);
        long lcg = 12345L;
        double product = 1.0;
        double logSum = 0.0;
        for (int i = 0; i < 1000; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double draw = standard.inverseCdf(((lcg >>> 11) + 0.5) * 0x1.0p-53);
            product *= standard.pdf(draw);
            logSum += standard.logPdf(draw);
        }
        assertEquals("the product survived, so this sample proves nothing", 0.0, product, 0.0);
        assertEquals("and its logarithm is no answer either", Double.NEGATIVE_INFINITY, Math.log(product), 0.0);
        // measured: -1434.1670
        assertTrue("the log-likelihood is " + logSum, logSum < -1400.0 && logSum > -1470.0);
    }

    // ------------------------------------------------------ the totality rule --

    @Test
    public void theLogDensityIsMinusInfinityOutsideTheSupportAndNaNAtNaN() {
        for (ContinuousDistribution d : continuous()) {
            String name = d.getClass().getSimpleName();
            double below = d.supportLowerBound();
            if (below > Double.NEGATIVE_INFINITY) {
                double outside = below - 1.0;
                assertEquals(name + " below its support", Double.NEGATIVE_INFINITY, d.logPdf(outside), 0.0);
                assertEquals(name + " below its support, in the density", 0.0, d.pdf(outside), 0.0);
            }
            double above = d.supportUpperBound();
            if (above < Double.POSITIVE_INFINITY) {
                double outside = above + 1.0;
                assertEquals(name + " above its support", Double.NEGATIVE_INFINITY, d.logPdf(outside), 0.0);
            }
            assertTrue(name + " answered a question that is not one", Double.isNaN(d.logPdf(Double.NaN)));
        }
        for (DiscreteDistribution d : discrete()) {
            String name = d.getClass().getSimpleName();
            long below = d.supportLowerBound() - 1L;
            if (below >= Integer.MIN_VALUE) {
                assertEquals(name + " below its support", Double.NEGATIVE_INFINITY, d.logPmf((int) below), 0.0);
            }
            long above = d.supportUpperBound() + 1L;
            if (above <= Integer.MAX_VALUE) {
                assertEquals(name + " above its support", Double.NEGATIVE_INFINITY, d.logPmf((int) above), 0.0);
            }
        }
    }

    @Test
    public void theDefaultIsTheLogarithmOfWhateverTheImplementationGives() {
        // an implementation outside this package inherits the fallback, which
        // has to be correct rather than merely present. This is also what makes
        // the two new methods a source- and binary-compatible addition
        ContinuousDistribution triangular = new ContinuousDistribution() {
            @Override
            public double pdf(double x) {
                return (x >= 0.0 && x <= 2.0) ? (1.0 - Math.abs(x - 1.0)) : 0.0;
            }

            @Override
            public double cdf(double x) {
                throw new UnsupportedOperationException();
            }

            @Override
            public double inverseCdf(double p) {
                throw new UnsupportedOperationException();
            }

            @Override
            public double mean() {
                return 1.0;
            }

            @Override
            public double variance() {
                return 1.0 / 6.0;
            }
        };
        assertEquals("at the peak", 0.0, triangular.logPdf(1.0), 0.0);
        assertEquals("halfway up", Math.log(0.5), triangular.logPdf(0.5), 0.0);
        assertEquals("outside", Double.NEGATIVE_INFINITY, triangular.logPdf(3.0), 0.0);

        DiscreteDistribution coin = new DiscreteDistribution() {
            @Override
            public double pmf(int k) {
                return (k == 0 || k == 1) ? 0.5 : 0.0;
            }

            @Override
            public double cdf(int k) {
                throw new UnsupportedOperationException();
            }

            @Override
            public double mean() {
                return 0.5;
            }

            @Override
            public double variance() {
                return 0.25;
            }

            @Override
            public int supportLowerBound() {
                return 0;
            }

            @Override
            public int supportUpperBound() {
                return 1;
            }
        };
        assertEquals("heads", Math.log(0.5), coin.logPmf(0), 0.0);
        assertEquals("neither", Double.NEGATIVE_INFINITY, coin.logPmf(7), 0.0);
    }
}
