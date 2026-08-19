/*
 * Copyright 2026 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.distribution;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

import org.junit.Test;

import math.cern.ProbabilityFuncs;

public class DiscreteDistributionTest {

    private static final double EPSILON = 1.0e-11;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /**
     * Wraps a distribution so that the generic quantile search of
     * {@link DiscreteDistribution} is used instead of a closed form override.
     */
    private static DiscreteDistribution generic(final DiscreteDistribution d) {
        return new DiscreteDistribution() {
            @Override
            public double pmf(int k) {
                return d.pmf(k);
            }

            @Override
            public double cdf(int k) {
                return d.cdf(k);
            }

            @Override
            public double mean() {
                return d.mean();
            }

            @Override
            public double variance() {
                return d.variance();
            }

            @Override
            public int supportLowerBound() {
                return d.supportLowerBound();
            }

            @Override
            public int supportUpperBound() {
                return d.supportUpperBound();
            }
        };
    }

    /** Last index worth summing to: the support end, capped for unbounded ones. */
    private static int upperSummationLimit(DiscreteDistribution d) {
        return Math.max(d.quantile(1.0 - 1.0e-15), d.supportLowerBound());
    }

    /**
     * The invariants every discrete distribution has to satisfy: a non-negative
     * mass function whose running sum is the CDF, a total mass of one, and
     * first and second moments that agree with the closed forms.
     */
    private static void assertDistributionInvariants(String name, DiscreteDistribution d) {
        int lo = d.supportLowerBound();
        int hi = upperSummationLimit(d);
        double sum = 0.0;
        double mean = 0.0;
        for (int k = lo; k <= hi; k++) {
            double f = d.pmf(k);
            assertTrue(name + ": negative pmf at " + k, f >= 0.0);
            sum += f;
            mean += (double) k * f;
            assertEquals(name + ": running pmf sum must equal the cdf at " + k, sum, d.cdf(k), EPSILON);
        }
        assertEquals(name + ": total mass", 1.0, sum, EPSILON);
        assertEquals(name + ": mean", d.mean(), mean, 1.0e-8 * (1.0 + Math.abs(mean)));

        double variance = 0.0;
        for (int k = lo; k <= hi; k++) {
            variance += (k - mean) * (k - mean) * d.pmf(k);
        }
        assertEquals(name + ": variance", d.variance(), variance, 1.0e-7 * (1.0 + Math.abs(variance)));

        assertEquals(name + ": pmf below the support", 0.0, d.pmf(lo - 1), 0.0);
        assertEquals(name + ": cdf below the support", 0.0, d.cdf(lo - 1), 0.0);
        assertEquals(name + ": cdf at the support end", 1.0, d.cdf(d.supportUpperBound()), 0.0);
    }

    /**
     * The quantile is the smallest k with cdf(k) &gt;= p, and a closed form
     * override has to agree with the generic search of the interface.
     */
    private static void assertQuantileInvariants(String name, DiscreteDistribution d) {
        int lo = d.supportLowerBound();
        DiscreteDistribution searched = generic(d);
        for (int i = 1; i < 500; i++) {
            double p = i / 500.0;
            int q = d.quantile(p);
            assertTrue(name + ": quantile below the support at p=" + p, q >= lo);
            assertTrue(name + ": quantile above the support at p=" + p, q <= d.supportUpperBound());
            assertTrue(name + ": cdf(q) >= p at p=" + p, d.cdf(q) >= p);
            if (q > lo) {
                assertTrue(name + ": cdf(q-1) < p at p=" + p, d.cdf(q - 1) < p);
            }
            assertEquals(name + ": closed form vs generic search at p=" + p, searched.quantile(p), q);
        }
        assertEquals(name + ": quantile(0)", lo, d.quantile(0.0));
        assertEquals(name + ": quantile(1)", d.supportUpperBound(), d.quantile(1.0));

        // round trip on the values that carry enough mass to be separable
        int hi = Math.min(upperSummationLimit(d), lo + 500);
        for (int k = lo; k <= hi; k++) {
            double c = d.cdf(k);
            if (c > 1.0e-12 && c < 1.0 - 1.0e-12 && d.pmf(k) > 1.0e-12) {
                assertEquals(name + ": quantile(cdf(" + k + "))", k, d.quantile(c));
            }
        }
    }

    // -----------------------------------------------------------------
    // per distribution invariants
    // -----------------------------------------------------------------

    @Test
    public void testDiscreteUniformInvariants() {
        int[][] params = { { 0, 9 }, { -5, 5 }, { 7, 7 }, { -13, -3 }, { 100, 400 } };
        for (int[] ab : params) {
            DiscreteUniform d = new DiscreteUniform(ab[0], ab[1]);
            String name = "U(" + ab[0] + "," + ab[1] + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    @Test
    public void testGeometricInvariants() {
        double[] probabilities = { 0.5, 0.01, 0.25, 0.999, 1.0 };
        for (double p : probabilities) {
            Geometric d = new Geometric(p);
            String name = "Geo(" + p + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    @Test
    public void testBinomialInvariants() {
        double[][] params = { { 20, 0.3 }, { 1, 0.5 }, { 500, 0.97 }, { 64, 0.5 }, { 200, 0.004 } };
        for (double[] np : params) {
            Binomial d = new Binomial((int) np[0], np[1]);
            String name = "Bin(" + (int) np[0] + "," + np[1] + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    @Test
    public void testPoissonInvariants() {
        double[] lambdas = { 4.0, 0.01, 1.0, 250.0 };
        for (double lambda : lambdas) {
            Poisson d = new Poisson(lambda);
            String name = "Poi(" + lambda + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    @Test
    public void testNegativeBinomialInvariants() {
        double[][] params = { { 1, 0.3 }, { 7, 0.4 }, { 3, 0.8 }, { 25, 0.5 } };
        for (double[] rp : params) {
            NegativeBinomial d = new NegativeBinomial((int) rp[0], rp[1]);
            String name = "NB(" + (int) rp[0] + "," + rp[1] + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    @Test
    public void testHypergeometricInvariants() {
        int[][] params = { { 50, 5, 10 }, { 50, 25, 25 }, { 20, 7, 12 }, { 10, 8, 5 }, { 200, 3, 199 } };
        for (int[] p : params) {
            Hypergeometric d = new Hypergeometric(p[0], p[1], p[2]);
            String name = "HG(" + p[0] + "," + p[1] + "," + p[2] + ")";
            assertDistributionInvariants(name, d);
            assertQuantileInvariants(name, d);
        }
    }

    // -----------------------------------------------------------------
    // cross checks against independent implementations
    // -----------------------------------------------------------------

    @Test
    public void testCdfAgreesWithTheCernComplements() {
        Binomial binomial = new Binomial(30, 0.42);
        for (int k = 0; k < 30; k++) {
            assertEquals("binomial complement at " + k, 1.0,
                    binomial.cdf(k) + ProbabilityFuncs.binomialComplemented(k, 30, 0.42), 1.0e-13);
        }
        Poisson poisson = new Poisson(6.5);
        for (int k = 0; k < 40; k++) {
            assertEquals("poisson complement at " + k, 1.0,
                    poisson.cdf(k) + ProbabilityFuncs.poissonComplemented(k, 6.5), 1.0e-13);
        }
        NegativeBinomial negativeBinomial = new NegativeBinomial(5, 0.35);
        for (int k = 0; k < 60; k++) {
            assertEquals("negative binomial complement at " + k, 1.0,
                    negativeBinomial.cdf(k) + ProbabilityFuncs.negativeBinomialComplemented(k, 5, 0.35), 1.0e-13);
        }
    }

    @Test
    public void testHypergeometricPmfAgreesWithExactBinomialCoefficients() {
        int[][] params = { { 50, 5, 10 }, { 20, 7, 12 }, { 30, 15, 15 }, { 10, 8, 5 } };
        for (int[] p : params) {
            Hypergeometric d = new Hypergeometric(p[0], p[1], p[2]);
            for (int k = d.supportLowerBound(); k <= d.supportUpperBound(); k++) {
                double exact = exactHypergeometricPmf(p[0], p[1], p[2], k);
                assertEquals("HG(" + p[0] + "," + p[1] + "," + p[2] + ") pmf at " + k, exact, d.pmf(k), 1.0e-12);
            }
        }
    }

    /** C(K,k) * C(N-K,n-k) / C(N,n) in exact integer arithmetic. */
    private static double exactHypergeometricPmf(int population, int successes, int draws, int k) {
        BigInteger numerator = binomialCoefficient(successes, k)
                .multiply(binomialCoefficient(population - successes, draws - k));
        if (numerator.signum() == 0) {
            return 0.0;
        }
        return new BigDecimal(numerator)
                .divide(new BigDecimal(binomialCoefficient(population, draws)), MathContext.DECIMAL64).doubleValue();
    }

    private static BigInteger binomialCoefficient(int a, int b) {
        if (b < 0 || b > a) {
            return BigInteger.ZERO;
        }
        BigInteger r = BigInteger.ONE;
        for (int i = 0; i < b; i++) {
            r = r.multiply(BigInteger.valueOf(a - i)).divide(BigInteger.valueOf(i + 1));
        }
        return r;
    }

    // -----------------------------------------------------------------
    // relations between the families
    // -----------------------------------------------------------------

    @Test
    public void testGeometricIsNegativeBinomialWithOneSuccess() {
        double[] probabilities = { 0.05, 0.5, 0.93, 1.0 };
        for (double p : probabilities) {
            Geometric geometric = new Geometric(p);
            NegativeBinomial negativeBinomial = new NegativeBinomial(1, p);
            for (int k = 0; k < 200; k++) {
                assertEquals("pmf at " + k + " for p=" + p, negativeBinomial.pmf(k), geometric.pmf(k), 1.0e-15);
                assertEquals("cdf at " + k + " for p=" + p, negativeBinomial.cdf(k), geometric.cdf(k), 1.0e-13);
            }
            assertEquals("mean for p=" + p, negativeBinomial.mean(), geometric.mean(), 0.0);
            assertEquals("variance for p=" + p, negativeBinomial.variance(), geometric.variance(), 0.0);
        }
    }

    @Test
    public void testBinomialApproachesPoisson() {
        Poisson limit = new Poisson(3.0);
        Binomial approximation = new Binomial(200000, 3.0 / 200000.0);
        for (int k = 0; k <= 15; k++) {
            assertEquals("pmf at " + k, limit.pmf(k), approximation.pmf(k), 1.0e-4);
        }
    }

    @Test
    public void testHypergeometricApproachesBinomial() {
        Binomial limit = new Binomial(20, 0.25);
        Hypergeometric approximation = new Hypergeometric(4000000, 1000000, 20);
        for (int k = 0; k <= 20; k++) {
            assertEquals("pmf at " + k, limit.pmf(k), approximation.pmf(k), 1.0e-5);
        }
    }

    // -----------------------------------------------------------------
    // numerics that a naive implementation would get wrong
    // -----------------------------------------------------------------

    @Test
    public void testHypergeometricWithAnUnderflowingSupportEnd() {
        // this population is peaked enough that the pmf at the ends of the
        // support underflows to zero, so a cdf summed from the far end would
        // return garbage in the centre. Cross checked against the normal limit.
        Hypergeometric d = new Hypergeometric(1000000000, 500000000, 1000000);
        double sd = Math.sqrt(d.variance());
        assertEquals("pmf at the support end must underflow", 0.0, d.pmf(0), 0.0);
        assertEquals("cdf at the mean", 0.5, d.cdf((int) Math.round(d.mean())), 1.0e-3);
        assertEquals("cdf two sigma below the mean", 0.02275, d.cdf((int) Math.round(d.mean() - 2.0 * sd)), 5.0e-4);
        assertEquals("cdf two sigma above the mean", 0.97725, d.cdf((int) Math.round(d.mean() + 2.0 * sd)), 5.0e-4);
        assertTrue("far lower tail", d.cdf((int) Math.round(d.mean() - 40.0 * sd)) < 1.0e-300);
        assertEquals("median", d.mean(), d.quantile(0.5), 3.0 * sd);
    }

    @Test
    public void testGeometricQuantileFarOutInTheTail() {
        Geometric d = new Geometric(1.0e-6);
        double[] probabilities = { 0.001, 0.5, 0.9, 0.999, 1.0 - 1.0e-9 };
        for (double p : probabilities) {
            assertEquals("closed form vs generic search at p=" + p, generic(d).quantile(p), d.quantile(p));
        }
    }

    @Test
    public void testExtremeSupportBoundsDoNotOverflow() {
        DiscreteUniform d = new DiscreteUniform(Integer.MIN_VALUE, Integer.MAX_VALUE);
        assertEquals(Integer.MIN_VALUE, d.quantile(0.0));
        assertEquals(Integer.MAX_VALUE, d.quantile(1.0));
        assertEquals(-0.5, d.mean(), 0.0);
        assertTrue("variance must not wrap", d.variance() > 0.0);
        assertEquals(1.0, new Poisson(2.0).cdf(Integer.MAX_VALUE), 0.0);
        assertEquals(1.0, new NegativeBinomial(2, 0.5).cdf(Integer.MAX_VALUE), 0.0);
        assertEquals(0.0, new NegativeBinomial(2, 0.5).pmf(Integer.MAX_VALUE), 0.0);
    }

    // -----------------------------------------------------------------
    // degenerate parameters
    // -----------------------------------------------------------------

    @Test
    public void testDegenerateParameters() {
        Binomial certainFailure = new Binomial(10, 0.0);
        assertEquals(1.0, certainFailure.pmf(0), 0.0);
        assertEquals(0.0, certainFailure.pmf(1), 0.0);
        assertEquals(1.0, certainFailure.cdf(0), 0.0);
        assertEquals(0, certainFailure.quantile(0.5));

        Binomial certainSuccess = new Binomial(10, 1.0);
        assertEquals(1.0, certainSuccess.pmf(10), 0.0);
        assertEquals(0.0, certainSuccess.cdf(9), 0.0);
        assertEquals(10, certainSuccess.quantile(0.5));

        Binomial noTrials = new Binomial(0, 0.4);
        assertEquals(1.0, noTrials.pmf(0), 0.0);
        assertEquals(0, noTrials.quantile(0.5));

        Poisson noEvents = new Poisson(0.0);
        assertEquals(1.0, noEvents.pmf(0), 0.0);
        assertEquals(0.0, noEvents.pmf(1), 0.0);
        assertEquals(0, noEvents.quantile(0.5));

        NegativeBinomial noFailures = new NegativeBinomial(3, 1.0);
        assertEquals(1.0, noFailures.pmf(0), 0.0);
        assertEquals(1.0, noFailures.cdf(0), 0.0);

        Hypergeometric allSuccesses = new Hypergeometric(30, 30, 10);
        assertEquals(10, allSuccesses.supportLowerBound());
        assertEquals(1.0, allSuccesses.pmf(10), 0.0);

        Hypergeometric noSuccesses = new Hypergeometric(30, 0, 10);
        assertEquals(1.0, noSuccesses.pmf(0), 0.0);

        Hypergeometric empty = new Hypergeometric(0, 0, 0);
        assertEquals(1.0, empty.pmf(0), 0.0);
        assertEquals(0.0, empty.mean(), 0.0);
        assertEquals(0.0, empty.variance(), 0.0);

        DiscreteUniform single = new DiscreteUniform(7, 7);
        assertEquals(1.0, single.pmf(7), 0.0);
        assertEquals(0.0, single.variance(), 0.0);
        assertEquals(7, single.quantile(0.5));
    }

    @Test
    public void testProbabilityOverAnInterval() {
        DiscreteUniform d = new DiscreteUniform(0, 9);
        assertEquals(0.3, d.probability(2, 5), 1.0e-15);
        assertEquals(0.0, d.probability(4, 4), 0.0);
        assertEquals(1.0, d.probability(-1, 9), 1.0e-15);
    }

    // -----------------------------------------------------------------
    // argument validation
    // -----------------------------------------------------------------

    @Test(expected = IllegalArgumentException.class)
    public void testProbabilityRejectsInvertedInterval() {
        new DiscreteUniform(0, 9).probability(5, 2);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testQuantileRejectsProbabilityAboveOne() {
        new Poisson(1.0).quantile(1.5);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testQuantileRejectsNaN() {
        new Geometric(0.5).quantile(Double.NaN);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testBinomialRejectsNegativeTrials() {
        new Binomial(-1, 0.5);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testPoissonRejectsNegativeRate() {
        new Poisson(-1.0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testGeometricRejectsZeroProbability() {
        new Geometric(0.0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testNegativeBinomialRejectsZeroSuccesses() {
        new NegativeBinomial(0, 0.5);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testHypergeometricRejectsTooManySuccesses() {
        new Hypergeometric(10, 11, 5);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testHypergeometricRejectsTooManyDraws() {
        new Hypergeometric(10, 5, 11);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testDiscreteUniformRejectsInvertedSupport() {
        new DiscreteUniform(5, 4);
    }

    // -----------------------------------------------------------------
    // randomized sweep
    // -----------------------------------------------------------------

    @Test
    public void testRandomizedParameterSweep() {
        long lcg = 20260819L;
        for (int trial = 0; trial < 40; trial++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            int population = 1 + (int) ((lcg >>> 33) % 120L);
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            int successes = (int) ((lcg >>> 33) % (population + 1L));
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            int draws = (int) ((lcg >>> 33) % (population + 1L));
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double p = 0.02 + 0.96 * ((lcg >>> 11) / (double) (1L << 53));

            DiscreteDistribution[] tested = { new Hypergeometric(population, successes, draws),
                    new Binomial(population, p), new Poisson(0.5 + population / 4.0),
                    new NegativeBinomial(1 + successes, p), new Geometric(p),
                    new DiscreteUniform(-successes, draws) };
            for (DiscreteDistribution d : tested) {
                String name = d.getClass().getSimpleName() + " trial " + trial;
                assertDistributionInvariants(name, d);
                int lo = d.supportLowerBound();
                DiscreteDistribution searched = generic(d);
                for (int i = 1; i < 40; i++) {
                    double probability = i / 40.0;
                    int q = d.quantile(probability);
                    assertTrue(name + ": cdf(q) >= p at p=" + probability, d.cdf(q) >= probability);
                    if (q > lo) {
                        assertTrue(name + ": cdf(q-1) < p at p=" + probability, d.cdf(q - 1) < probability);
                    }
                    assertEquals(name + ": closed form vs generic search at p=" + probability,
                            searched.quantile(probability), q);
                }
            }
        }
    }
}
