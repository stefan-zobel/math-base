package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.cern.Arithmetic;
import math.cern.FastGamma;
import math.linalg.DMatrix;
import math.rng.DefaultRng;
import math.rng.DirichletSampler;
import math.rng.PseudoRandom;
import math.solve.Quadrature;

/**
 * The Dirichlet against classes that know nothing about it: {@link Beta}, which
 * it reduces to at two components and which every marginal is, {@link Binomial}
 * and {@link Quadrature}, which reproduce the evidence by integrating a
 * likelihood, and {@link DirichletSampler}, which draws from it.
 * <p>
 * The fit of the sampler against this class lives in
 * {@code math.rng.SingleDrawTest}.
 */
public class DirichletTest {

    /**
     * The density integrated over the simplex, which has to be one.
     * <p>
     * With {@code x = u} and {@code y = (1 - u) v} the triangle becomes the
     * unit square and the Jacobian is {@code 1 - u}, so an ordinary rectangular
     * rule can take it.
     */
    private static double integrateOverTheSimplex(final Dirichlet d) {
        return Quadrature.integrate((u, v) -> {
            double x = u;
            double y = (1.0 - u) * v;
            double z = 1.0 - x - y;
            if (!(x > 0.0) || !(y > 0.0) || !(z > 0.0)) {
                return 0.0;
            }
            // the three coordinates are formed by subtraction and need not sum
            // to one closely enough for logPdf, which is entitled to insist
            double sum = x + y + z;
            return (1.0 - u) * d.pdf(new double[] { x / sum, y / sum, z / sum });
        }, 0.0, 1.0, 0.0, 1.0, 1.0e-12);
    }

    @Test
    public void theDensityIntegratesToOneOverTheSimplex() {
        // every concentration is either exactly one, where the factor is
        // constant, or at least two, where the exponent is at least one and the
        // integrand is smooth to the edge. In between it is not the class that
        // fails but the rule: at alpha = 1.5 the factor is a square root, whose
        // derivative is infinite at the edge, and the integral reaches only
        // 6.1e-07. Below one there is an integrable pole and it reaches 1.8e-01.
        // Those shapes are covered by the beta comparison below instead
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 2.0, 3.0, 4.0 }, { 20.0, 20.0, 20.0 }, { 60.0, 5.0, 1.0 },
                { 1.0, 1.0, 8.0 }, { 3.0, 3.0, 3.0 }, { 200.0, 200.0, 200.0 } };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            double integral = integrateOverTheSimplex(new Dirichlet(settings[s]));
            double gap = Math.abs(integral - 1.0);
            if (gap > worst) {
                worst = gap;
                at = Arrays.toString(settings[s]) + " integrates to " + integral;
            }
        }
        // measured: 1.1e-13, at the most concentrated of the seven
        assertTrue("the density integrates to " + worst + " off one, " + at, worst < 1.0e-12);
    }

    private static final double[][] PAIRS = { { 0.5, 0.5 }, { 1.0, 1.0 }, { 2.0, 5.0 }, { 0.1, 3.0 }, { 30.0, 30.0 },
            { 1.0, 400.0 }, { 1.0e-3, 1.0e-3 }, { 500.0, 500.0 } };

    @Test
    public void atTwoComponentsItIsABeta() {
        // the strongest check of the normalizing factor that does not need an
        // integral, and it reaches the shapes below one that the integral
        // cannot
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < PAIRS.length; s++) {
            double a = PAIRS[s][0];
            double b = PAIRS[s][1];
            Dirichlet d = new Dirichlet(new double[] { a, b });
            Beta beta = new Beta(a, b);
            for (double x = 1.0e-6; x < 1.0; x *= 1.7) {
                double viaDirichlet = d.logPdf(new double[] { x, 1.0 - x });
                double viaBeta = beta.logPdf(x);
                // a log density crosses zero, so the scale has to be bounded
                // below or the comparison divides by nothing
                double relative = Math.abs(viaDirichlet - viaBeta) / Math.max(1.0, Math.abs(viaBeta));
                if (relative > worst) {
                    worst = relative;
                    at = "alpha=" + a + "," + b + " x=" + x + " : " + viaDirichlet + " vs " + viaBeta;
                }
            }
        }
        // measured over 216 points: 8.3e-15, and 164 of them agree bit for bit.
        // The rest cannot: Beta forms the second factor as log1p(-x) where this
        // class takes the logarithm of the second component
        assertTrue("the two differ by " + worst + ", " + at, worst < 1.0e-13);
    }

    @Test
    public void mergingComponentsGivesADirichlet() {
        // the marginal of component zero is the two component Dirichlet with
        // the other concentrations added up, which is checked here through the
        // class itself rather than against a formula restated in the test
        double worst = 0.0;
        String at = "";
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 2.0, 3.0, 4.0 }, { 0.5, 0.5, 0.5 }, { 30.0, 1.0, 9.0 } };
        for (int s = 0; s < settings.length; s++) {
            double[] a = settings[s];
            Beta marginal = new Dirichlet(a).marginal(0);
            Dirichlet merged = new Dirichlet(new double[] { a[0], a[1] + a[2] });
            for (double x = 0.02; x < 1.0; x += 0.07) {
                double viaMerge = merged.logPdf(new double[] { x, 1.0 - x });
                double relative = Math.abs(viaMerge - marginal.logPdf(x)) / Math.max(1.0, Math.abs(marginal.logPdf(x)));
                if (relative > worst) {
                    worst = relative;
                    at = Arrays.toString(a) + " x=" + x;
                }
            }
        }
        // measured: 8.6e-16
        assertTrue("the two differ by " + worst + ", " + at, worst < 1.0e-13);
    }

    @Test
    public void theMarginalAgreesWithTheMomentsAndTheRowsSumToZero() {
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 2.0, 3.0, 4.0, 5.0 }, { 0.1, 0.2, 100.0 }, { 7.0, 7.0 },
                { 1.0e-6, 1.0e6, 3.0 }, { 1.0e8, 1.0, 1.0, 1.0, 1.0 }, { 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 } };
        double worstMean = 0.0;
        double worstVariance = 0.0;
        double worstRow = 0.0;
        for (int s = 0; s < settings.length; s++) {
            Dirichlet d = new Dirichlet(settings[s]);
            int m = d.components();
            double[] mean = new double[m];
            d.mean(mean);
            DMatrix c = d.covariance();
            for (int i = 0; i < m; i++) {
                Beta marginal = d.marginal(i);
                worstMean = Math.max(worstMean, Math.abs(marginal.mean() - mean[i]) / mean[i]);
                worstVariance = Math.max(worstVariance, Math.abs(marginal.variance() - c.get(i, i)) / c.get(i, i));
                double row = 0.0;
                for (int j = 0; j < m; j++) {
                    row += c.get(i, j);
                }
                worstRow = Math.max(worstRow, Math.abs(row));
            }
        }
        // the two means came out bit identical at all 30 components measured,
        // but that is not owed: the marginal's denominator is the sum of the
        // other concentrations plus this one, and alpha_0 is the sum taken in
        // its own order
        assertTrue("the marginal mean differs by " + worstMean, worstMean < 1.0e-15);
        // measured: 3.4e-16
        assertTrue("the marginal variance differs by " + worstVariance, worstVariance < 1.0e-14);
        // the components sum to one, so the covariance matrix is singular and
        // every row of it sums to zero. Measured: 1.7e-18
        assertTrue("a row of the covariance sums to " + worstRow, worstRow < 1.0e-15);
    }

    @Test
    public void theEvidenceSumsToOneOverEveryCountVector() {
        // the multinomial coefficient is part of the evidence, so what it
        // returns is the probability of a count vector and the probabilities of
        // all of them add up
        double[][] settings = { { 1.0, 1.0, 1.0 }, { 0.5, 2.0, 7.0 }, { 100.0, 1.0, 3.0 } };
        int[] totals = { 5, 12 };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < settings.length; s++) {
            Dirichlet d = new Dirichlet(settings[s]);
            for (int t = 0; t < totals.length; t++) {
                int n = totals[t];
                double sum = 0.0;
                for (int i = 0; i <= n; i++) {
                    for (int j = 0; j <= n - i; j++) {
                        sum += Math.exp(d.logMarginalLikelihood(new int[] { i, j, n - i - j }));
                    }
                }
                if (Math.abs(sum - 1.0) > worst) {
                    worst = Math.abs(sum - 1.0);
                    at = Arrays.toString(settings[s]) + " n=" + n + " sums to " + sum;
                }
            }
        }
        // and once with four components, where there are 84 of them
        Dirichlet four = new Dirichlet(new double[] { 1.5, 0.3, 4.0, 2.0 });
        int n = 6;
        double sum = 0.0;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                for (int k = 0; k <= n - i - j; k++) {
                    sum += Math.exp(four.logMarginalLikelihood(new int[] { i, j, k, n - i - j - k }));
                }
            }
        }
        if (Math.abs(sum - 1.0) > worst) {
            worst = Math.abs(sum - 1.0);
            at = "four components, n=6, sums to " + sum;
        }
        // measured over up to 91 count vectors: 6.2e-15
        assertTrue("the evidence sums to " + worst + " off one, " + at, worst < 1.0e-13);
    }

    @Test
    public void theEvidenceIsTheIntegralOfTheLikelihoodAgainstThePrior() {
        // at two components the evidence is a beta-binomial, so it must be what
        // integrating a Binomial mass against a Beta density gives -- two other
        // classes and a quadrature rule, none of which knows this one.
        // Concentrations below one are left out because the Beta then has a pole
        // at either end and it is the rule that fails, not the identity
        double worst = 0.0;
        String at = "";
        double[][] settings = { { 1.0, 1.0 }, { 2.0, 5.0 }, { 20.0, 3.0 }, { 1.0, 8.0 } };
        int[] totals = { 3, 10 };
        for (int s = 0; s < settings.length; s++) {
            Dirichlet d = new Dirichlet(settings[s]);
            final Beta prior = new Beta(settings[s][0], settings[s][1]);
            for (int t = 0; t < totals.length; t++) {
                final int n = totals[t];
                for (int a = 0; a <= n; a++) {
                    final int hits = a;
                    double direct = Math.exp(d.logMarginalLikelihood(new int[] { hits, n - hits }));
                    double integral = Quadrature.integrate(p -> new Binomial(n, p).pmf(hits) * prior.pdf(p), 0.0, 1.0,
                            1.0e-13);
                    double relative = Math.abs(direct - integral) / integral;
                    if (relative > worst) {
                        worst = relative;
                        at = Arrays.toString(settings[s]) + " n=" + n + " a=" + hits + " : " + direct + " vs "
                                + integral;
                    }
                }
            }
        }
        // measured over 64 count vectors: 6.9e-15
        assertTrue("the two differ by " + worst + " relative, " + at, worst < 1.0e-12);
    }

    @Test
    public void theEvidenceSurvivesAConcentrationTheSubtractionWouldNot() {
        // with both concentrations enormous and equal the proportions are
        // pinned at one half, so the evidence is a plain binomial probability
        Dirichlet d = new Dirichlet(new double[] { 5.0e14, 5.0e14 });
        int[] counts = { 30, 10 };
        double evidence = Math.exp(d.logMarginalLikelihood(counts));
        double binomial = new Binomial(40, 0.5).pmf(30);
        assertEquals("the evidence is not the binomial probability", binomial, evidence, 1.0e-11 * binomial);

        // the same quantity written as a difference of two log gammas, which is
        // the subtraction FastGamma.logGammaRatio exists to avoid
        double naive = Arithmetic.logFactorial(40) - (FastGamma.logGamma(1.0e15 + 40) - FastGamma.logGamma(1.0e15));
        for (int i = 0; i < 2; i++) {
            naive += FastGamma.logGamma(5.0e14 + counts[i]) - FastGamma.logGamma(5.0e14);
            naive -= Arithmetic.logFactorial(counts[i]);
        }
        // measured: 7.93e-05 where 7.71e-04 is right, an order of magnitude out.
        // At this size the difference of the two logarithms is below one ulp of
        // either of them and comes back quantized
        assertTrue("the subtraction should have failed here, it gave " + Math.exp(naive),
                Math.exp(naive) < 0.5 * binomial);
    }

    @Test
    public void thePosteriorComposes() {
        Dirichlet prior = new Dirichlet(new double[] { 0.5, 2.0, 7.0, 1.0 });
        int[] first = { 3, 0, 11, 4 };
        int[] second = { 1, 7, 0, 2 };
        int[] both = new int[4];
        for (int i = 0; i < 4; i++) {
            both[i] = first[i] + second[i];
        }
        Dirichlet stepwise = prior.posterior(first).posterior(second);
        Dirichlet once = prior.posterior(both);
        for (int i = 0; i < 4; i++) {
            // both routes add the same integers to the same double, so this is
            // an identity and not an approximation
            assertEquals("component " + i, stepwise.concentration(i), once.concentration(i), 0.0);
        }
        assertEquals("the total", stepwise.totalConcentration(), once.totalConcentration(), 0.0);

        // and no counts at all changes nothing
        Dirichlet unchanged = prior.posterior(new int[] { 0, 0, 0, 0 });
        for (int i = 0; i < 4; i++) {
            assertEquals("component " + i, prior.concentration(i), unchanged.concentration(i), 0.0);
        }
    }

    @Test
    public void theCornersOfTheSimplexAreTheShapesAlone() {
        // a concentration of one contributes nothing whatever the component is.
        // Taken as written the term would read 0.0 * -Infinity there, which is
        // NaN, and the uniform density on the simplex would have no value at
        // its own corner
        Dirichlet uniform = new Dirichlet(new double[] { 1.0, 1.0, 1.0 });
        assertEquals("the uniform density at a corner", Math.log(2.0), uniform.logPdf(new double[] { 0.0, 0.0, 1.0 }),
                0.0);
        assertEquals("and in the middle", Math.log(2.0),
                uniform.logPdf(new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 }), 1.0e-15);

        Dirichlet mixed = new Dirichlet(new double[] { 0.5, 2.0, 1.0 });
        assertEquals("a concentration below one is a pole", Double.POSITIVE_INFINITY,
                mixed.logPdf(new double[] { 0.0, 0.5, 0.5 }), 0.0);
        assertEquals("a concentration above one is a zero", Double.NEGATIVE_INFINITY,
                mixed.logPdf(new double[] { 0.5, 0.0, 0.5 }), 0.0);
        // where the two meet the answer depends on the direction the corner is
        // approached from, and there is no number that says so
        assertTrue("a pole meeting a zero should not have a value",
                Double.isNaN(mixed.logPdf(new double[] { 0.0, 0.0, 1.0 })));

        assertTrue("a NaN component", Double.isNaN(uniform.logPdf(new double[] { Double.NaN, 0.5, 0.5 })));
    }

    @Test
    public void theLogDensityAnswersWhereTheDensityHasUnderflowed() {
        // the reason the class answers in logarithms at all
        Dirichlet d = new Dirichlet(new double[] { 400.0, 400.0, 400.0 });
        double[] far = { 0.01, 0.01, 0.98 };
        assertEquals("the density is not zero out there", 0.0, d.pdf(far), 0.0);
        // measured: -2361.05
        double logDensity = d.logPdf(far);
        assertTrue("the log density is " + logDensity, logDensity < -2000.0 && logDensity > -3000.0);

        double[] centre = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
        assertTrue("and the peak is an ordinary number", d.pdf(centre) > 900.0 && d.pdf(centre) < 1100.0);
    }

    @Test
    public void theCovarianceSurvivesAConcentrationTheSquareWouldNot() {
        // alpha_0^2 leaves the double range from about 1.3e154 upwards, and
        // forming the quotients first is what reaches an answer that is
        // perfectly representable
        double[] alpha = { 3.0e300, 1.0e300, 6.0e300 };
        Dirichlet d = new Dirichlet(alpha);
        double a0 = d.totalConcentration();
        assertEquals("the square of the total is representable after all", Double.POSITIVE_INFINITY, a0 * a0, 0.0);

        DMatrix c = d.covariance();
        double m0 = alpha[0] / a0;
        assertEquals("the leading variance", m0 * (1.0 - m0) / (a0 + 1.0), c.get(0, 0), 0.0);
        assertTrue("the variance came back as " + c.get(0, 0), c.get(0, 0) > 0.0 && c.get(0, 0) < 1.0e-300);
        assertTrue("the covariance came back as " + c.get(0, 1), c.get(0, 1) < 0.0);
    }

    @Test
    public void theMarginalSurvivesACancellationTheDifferenceWouldNot() {
        // alpha_0 - alpha_i is exactly zero when one component carries all of
        // the concentration, and Beta refuses a shape of zero. Summing the
        // others instead is what keeps the marginal available
        Dirichlet d = new Dirichlet(new double[] { 1.0e300, 1.0, 2.0 });
        assertEquals("the difference should have cancelled", 0.0, d.totalConcentration() - 1.0e300, 0.0);
        // Beta(1e300, 0.0) would have been refused outright; built from the sum
        // of the others the second shape is three
        Beta marginal = d.marginal(0);
        assertEquals("the marginal mean", 1.0, marginal.mean(), 1.0e-15);
        // its variance is 3e-600 and underflows, which is a fact about the
        // double and not about the shapes it was built from
        assertEquals("the marginal variance", 0.0, marginal.variance(), 0.0);
        // and the small components still have marginals of their own, which is
        // where the concentration ratio survives
        assertTrue("the second component's mean is " + d.marginal(1).mean(),
                d.marginal(1).mean() > 0.0 && d.marginal(1).mean() < 1.0e-299);
        assertEquals("and twice that for the third", 2.0 * d.marginal(1).mean(), d.marginal(2).mean(), 1.0e-314);
    }

    @Test
    public void theSimplexToleranceAcceptsRoundingAndRefusesCounts() {
        Dirichlet d = new Dirichlet(new double[] { 2.0, 3.0 });
        // the sampler in this library leaves a residual of at most 3.6e-15 over
        // component counts up to a thousand, so anything of that size has to be
        // accepted
        double reference = d.logPdf(new double[] { 0.5, 0.5 });
        assertEquals("a residual of 1e-12", reference, d.logPdf(new double[] { 0.5, 0.5 + 1.0e-12 }), 1.0e-11);
        assertEquals("a residual of 1e-10", reference, d.logPdf(new double[] { 0.5, 0.5 + 1.0e-10 }), 1.0e-9);
        // and the mistake it is here to catch misses by a factor
        refusesArgument("a vector summing to 0.9", d, new double[] { 0.45, 0.45 });
        refusesArgument("a vector of counts", d, new double[] { 3.0, 5.0 });
        refusesArgument("a residual of 1e-8", d, new double[] { 0.5, 0.5 + 1.0e-8 });
        refusesArgument("a negative component", d, new double[] { -0.1, 1.1 });
        refusesArgument("the wrong length", d, new double[] { 0.3, 0.3, 0.4 });

        // the draws of the library's own sampler are all accepted
        DirichletSampler sampler = DirichletSampler.of(new double[] { 0.01, 5.0, 100.0, 0.5 });
        Dirichlet law = new Dirichlet(new double[] { 0.01, 5.0, 100.0, 0.5 });
        long lcg = 20260827L;
        for (int i = 0; i < 500; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double[] x = new double[4];
            sampler.sample(DefaultRng.newPseudoRandom(lcg), x);
            double value = law.logPdf(x);
            assertTrue("a draw of the sampler was not on the simplex : " + Arrays.toString(x), !Double.isNaN(value));
        }
    }

    @Test
    public void aSingleComponentIsOneWithProbabilityOne() {
        Dirichlet d = new Dirichlet(new double[] { 2.5 });
        // the density with respect to a measure on a point is one, and the
        // formula gives that on its own: the normalizing factor cancels and the
        // only component is one
        assertEquals("the log density at the only point there is", 0.0, d.logPdf(new double[] { 1.0 }), 0.0);
        double[] mean = new double[1];
        d.mean(mean);
        assertEquals("the mean", 1.0, mean[0], 0.0);
        assertEquals("the variance", 0.0, d.covariance().get(0, 0), 0.0);
        try {
            d.marginal(0);
            fail("a single component was given a beta marginal");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    @Test
    public void theDrawIsTheSamplersDraw() {
        double[] alpha = { 0.7, 3.0, 12.0 };
        Dirichlet d = new Dirichlet(alpha);
        DirichletSampler sampler = DirichletSampler.of(alpha);
        for (long seed = 1L; seed <= 20L; seed++) {
            double[] viaDistribution = new double[3];
            double[] viaSampler = new double[3];
            double[] viaGenerator = new double[3];
            d.sample(DefaultRng.newPseudoRandom(seed), viaDistribution);
            sampler.sample(DefaultRng.newPseudoRandom(seed), viaSampler);
            PseudoRandom prng = DefaultRng.newPseudoRandom(seed);
            prng.nextDirichlet(alpha, viaGenerator);
            for (int i = 0; i < 3; i++) {
                // all three routes run the same code over the same generator
                assertEquals("seed " + seed + " component " + i, viaSampler[i], viaDistribution[i], 0.0);
                assertEquals("seed " + seed + " component " + i, viaSampler[i], viaGenerator[i], 0.0);
            }
        }
    }

    @Test
    public void permutingTheConcentrationsPermutesEverything() {
        double[] alpha = { 0.5, 4.0, 2.0 };
        double[] swapped = { 2.0, 4.0, 0.5 };
        Dirichlet d = new Dirichlet(alpha);
        Dirichlet s = new Dirichlet(swapped);
        double[] x = { 0.2, 0.3, 0.5 };
        double[] reversed = { 0.5, 0.3, 0.2 };
        assertEquals("the log density", d.logPdf(x), s.logPdf(reversed), 1.0e-14);
        assertEquals("the evidence", d.logMarginalLikelihood(new int[] { 2, 5, 1 }),
                s.logMarginalLikelihood(new int[] { 1, 5, 2 }), 1.0e-14);
        assertEquals("the marginal mean", d.marginal(0).mean(), s.marginal(2).mean(), 0.0);
    }

    @Test
    public void theInputIsCopiedRatherThanKept() {
        double[] alpha = { 1.0, 2.0, 3.0 };
        Dirichlet d = new Dirichlet(alpha);
        alpha[0] = 99.0;
        assertEquals("the constructor kept the caller's array", 1.0, d.concentration(0), 0.0);
        assertEquals("the total moved with it", 6.0, d.totalConcentration(), 0.0);
    }

    @Test
    public void theDirichletRejectsWhatItCannotBe() {
        refuses("null", null);
        refuses("an empty array", new double[0]);
        refuses("a zero concentration", new double[] { 1.0, 0.0 });
        refuses("a negative concentration", new double[] { 1.0, -1.0 });
        refuses("an infinite concentration", new double[] { 1.0, Double.POSITIVE_INFINITY });
        refuses("a NaN concentration", new double[] { 1.0, Double.NaN });

        Dirichlet d = new Dirichlet(new double[] { 1.0, 2.0 });
        refusesCounts("null counts", d, null);
        refusesCounts("counts of the wrong length", d, new int[] { 1, 2, 3 });
        refusesCounts("a negative count", d, new int[] { 1, -1 });
        try {
            d.mean(new double[3]);
            fail("the mean accepted an array of the wrong length");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            d.concentration(2);
            fail("a component out of range was accepted");
        } catch (IndexOutOfBoundsException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refuses(String what, double[] alpha) {
        try {
            new Dirichlet(alpha);
            fail("Dirichlet accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refusesArgument(String what, Dirichlet d, double[] x) {
        try {
            d.logPdf(x);
            fail("logPdf accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void refusesCounts(String what, Dirichlet d, int[] counts) {
        try {
            d.posterior(counts);
            fail("posterior accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
        try {
            d.logMarginalLikelihood(counts);
            fail("logMarginalLikelihood accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
