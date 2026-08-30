package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.distribution.Beta;
import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;

/**
 * The Dirichlet sampler, checked against its marginals -- which are Beta, and
 * therefore checkable exactly and continuously rather than through bins.
 */
public final class DirichletSamplerTest {

    private static double concentrationSum(double[] alpha) {
        double sum = 0.0;
        for (int i = 0; i < alpha.length; i++) {
            sum += alpha[i];
        }
        return sum;
    }

    @Test
    public void testEveryMarginalIsTheBetaItShouldBe() {
        // the marginal of a Dirichlet is Beta(alpha_i, alpha_0 - alpha_i), so
        // this is a continuous exact cross-check, and the p-values of several
        // hundred of them have to be uniform
        Uniform uniform = new Uniform(0.0, 1.0);
        double[][] cases = { { 1.0, 1.0, 1.0 }, { 0.5, 0.5 }, { 2.0, 3.0, 5.0, 10.0 },
                { 0.1, 0.1, 0.1 }, { 100.0, 1.0 } };
        for (double[] alpha : cases) {
            DirichletSampler sampler = DirichletSampler.of(alpha);
            double total = concentrationSum(alpha);
            for (int component = 0; component < alpha.length; component++) {
                int reps = 150;
                int draws = 2000;
                double[] pValues = new double[reps];
                Beta marginal = new Beta(alpha[component], total - alpha[component]);
                for (int r = 0; r < reps; r++) {
                    PseudoRandom prng = DefaultRng.newPseudoRandom(r * 7919L + 1);
                    double[] out = new double[alpha.length];
                    double[] drawn = new double[draws];
                    for (int d = 0; d < draws; d++) {
                        sampler.sample(prng, out);
                        drawn[d] = out[component];
                    }
                    pValues[r] = HypothesisTests
                            .kolmogorovSmirnov(drawn, marginal, Alternative.TWO_SIDED).pValue;
                }
                double uniformity = HypothesisTests
                        .kolmogorovSmirnov(pValues, uniform, Alternative.TWO_SIDED).pValue;
                // measured over the five vectors: 0.1450 .. 0.6919
                assertTrue(Arrays.toString(alpha) + " component " + component + ": uniformity p = "
                        + uniformity, uniformity > 0.001);
            }
        }
    }

    @Test
    public void testTwoComponentsAreABeta() {
        DirichletSampler sampler = DirichletSampler.of(new double[] { 2.0, 5.0 });
        PseudoRandom prng = DefaultRng.newPseudoRandom(9L);
        double[] out = new double[2];
        double[] first = new double[200000];
        for (int d = 0; d < first.length; d++) {
            sampler.sample(prng, out);
            first[d] = out[0];
        }
        double pValue = HypothesisTests.kolmogorovSmirnov(first, new Beta(2.0, 5.0), Alternative.TWO_SIDED)
                .pValue;
        // measured: 0.9837
        assertTrue("the first component does not fit Beta(2, 5), p = " + pValue, pValue > 0.001);
    }

    @Test
    public void testTheComponentsAlwaysSumToOne() {
        // and never come back as NaN, which is what the logarithmic form is
        // for: the normalization subtracts the largest before exponentiating,
        // so one term of the sum is exactly one and the sum cannot underflow
        for (double concentration : new double[] { 100.0, 1.0, 0.01, 0.001, 1.0e-6 }) {
            DirichletSampler sampler = DirichletSampler
                    .of(new double[] { concentration, concentration, concentration });
            PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
            double[] out = new double[3];
            double worst = 0.0;
            for (int d = 0; d < 100000; d++) {
                sampler.sample(prng, out);
                double sum = 0.0;
                for (int i = 0; i < 3; i++) {
                    assertTrue("alpha=" + concentration + ": a component came back as " + out[i],
                            out[i] >= 0.0 && out[i] <= 1.0);
                    sum += out[i];
                }
                worst = Math.max(worst, Math.abs(sum - 1.0));
            }
            // measured at every concentration from 1e-6 up: 3.3e-16
            assertTrue("alpha=" + concentration + ": the sum was off by " + worst, worst < 1.0e-14);
        }
    }

    @Test
    public void testTheShapeAtWhichTheDirectRouteWouldHaveFailed() {
        // the direct route -- draw Gamma, sum, divide -- underflows: at a shape
        // of 0.001 the gamma sampler returns exactly zero for 47562 of 100000
        // uniforms, so all three components of a draw underflow together often
        // enough to matter. It is measured here rather than asserted, because
        // it is the whole reason the sampler works in logarithms
        PseudoRandom prng = DefaultRng.newPseudoRandom(3L);
        int allUnderflowed = 0;
        int trials = 20000;
        for (int d = 0; d < trials; d++) {
            double sum = 0.0;
            for (int i = 0; i < 3; i++) {
                sum += GammaSpliterator.sample(prng, 0.001, 1.0);
            }
            if (sum == 0.0) {
                allUnderflowed++;
            }
        }
        // measured: 2175 of 20000, which the direct route would have returned
        // as three NaN
        assertTrue("the direct route did not underflow at all, so the log form "
                + "would be pointless: " + allUnderflowed, allUnderflowed > 200);

        // and the sampler itself, at the same shape, returns proportions
        DirichletSampler sampler = DirichletSampler.of(new double[] { 0.001, 0.001, 0.001 });
        double[] out = new double[3];
        double[] mean = new double[3];
        long[] largest = new long[3];
        int draws = 200000;
        for (int d = 0; d < draws; d++) {
            sampler.sample(prng, out);
            int best = 0;
            for (int i = 0; i < 3; i++) {
                assertTrue("a component came back as " + out[i], out[i] >= 0.0 && out[i] <= 1.0);
                mean[i] += out[i];
                if (out[i] > out[best]) {
                    best = i;
                }
            }
            largest[best]++;
        }
        // at this concentration the marginal is so nearly a two-point law at
        // zero and one that a Kolmogorov-Smirnov test says nothing: a third of
        // the true mass sits below 1e-300, which a double cannot tell from
        // zero. What still characterizes it is the mean and the symmetry, and
        // both are exact statements
        for (int i = 0; i < 3; i++) {
            // measured: 0.33362, 0.33319, 0.33319
            assertEquals("the mean of component " + i, 1.0 / 3.0, mean[i] / draws, 0.01);
            // measured: 66714, 66647, 66639 of 200000
            assertEquals("how often component " + i + " is the largest", draws / 3.0, largest[i],
                    0.02 * draws);
        }
    }

    @Test
    public void testTheMeanIsTheNormalizedConcentration() {
        double[][] cases = { { 2.0, 3.0, 5.0 }, { 0.5, 0.5, 0.5, 0.5 }, { 50.0, 1.0 } };
        for (double[] alpha : cases) {
            DirichletSampler sampler = DirichletSampler.of(alpha);
            double total = concentrationSum(alpha);
            PseudoRandom prng = DefaultRng.newPseudoRandom(5L);
            double[] out = new double[alpha.length];
            double[] mean = new double[alpha.length];
            int draws = 200000;
            for (int d = 0; d < draws; d++) {
                sampler.sample(prng, out);
                for (int i = 0; i < alpha.length; i++) {
                    mean[i] += out[i];
                }
            }
            for (int i = 0; i < alpha.length; i++) {
                // measured: 0.19986 / 0.30007 / 0.50007 against 0.2 / 0.3 / 0.5
                assertEquals(Arrays.toString(alpha) + " component " + i, alpha[i] / total, mean[i] / draws,
                        0.005);
            }
        }
    }

    @Test
    public void testTheSymmetricCaseIsUniformOnTheSimplex() {
        // alpha all one is the uniform distribution over the simplex, whose
        // marginals are Beta(1, m - 1)
        int m = 4;
        double[] alpha = new double[m];
        Arrays.fill(alpha, 1.0);
        DirichletSampler sampler = DirichletSampler.of(alpha);
        PseudoRandom prng = DefaultRng.newPseudoRandom(13L);
        double[] out = new double[m];
        double[] first = new double[200000];
        for (int d = 0; d < first.length; d++) {
            sampler.sample(prng, out);
            first[d] = out[0];
        }
        double pValue = HypothesisTests
                .kolmogorovSmirnov(first, new Beta(1.0, m - 1.0), Alternative.TWO_SIDED).pValue;
        assertTrue("the marginal is not Beta(1, 3), p = " + pValue, pValue > 0.001);
    }

    @Test
    public void testOneComponentTakesEverything() {
        DirichletSampler sampler = DirichletSampler.of(new double[] { 2.5 });
        assertEquals("components", 1, sampler.components());
        PseudoRandom prng = DefaultRng.newPseudoRandom(17L);
        double[] out = new double[1];
        for (int d = 0; d < 1000; d++) {
            sampler.sample(prng, out);
            assertEquals("the only component", 1.0, out[0], 0.0);
        }
    }

    @Test
    public void testTheSamplerIsReproducibleFromTheSeed() {
        DirichletSampler sampler = DirichletSampler.of(new double[] { 0.7, 1.3, 4.0 });
        double[] first = new double[3];
        double[] again = new double[3];
        PseudoRandom one = DefaultRng.newPseudoRandom(23L);
        PseudoRandom two = DefaultRng.newPseudoRandom(23L);
        for (int d = 0; d < 500; d++) {
            sampler.sample(one, first);
            sampler.sample(two, again);
            for (int i = 0; i < 3; i++) {
                assertEquals("draw " + d + " at " + i, first[i], again[i], 0.0);
            }
        }
    }

    @Test
    public void testTheInputIsCopiedRatherThanKept() {
        double[] alpha = { 1.0, 1.0 };
        DirichletSampler sampler = DirichletSampler.of(alpha);
        alpha[0] = 1000.0;
        assertEquals("the sampler followed a change to the caller's array", 1.0,
                sampler.concentration(0), 0.0);
    }

    @Test
    public void testTheSamplerRejectsWhatItCannotDraw() {
        rejectsBuild("null concentrations", null);
        rejectsBuild("no concentrations", new double[0]);
        rejectsBuild("a concentration of zero", new double[] { 1.0, 0.0 });
        rejectsBuild("a negative concentration", new double[] { 1.0, -0.5 });
        rejectsBuild("a NaN concentration", new double[] { 1.0, Double.NaN });
        rejectsBuild("an infinite concentration", new double[] { 1.0, Double.POSITIVE_INFINITY });

        DirichletSampler sampler = DirichletSampler.of(new double[] { 1.0, 1.0, 1.0 });
        PseudoRandom prng = DefaultRng.newPseudoRandom(1L);
        rejectsDraw("a null generator", sampler, null, new double[3]);
        rejectsDraw("a null output array", sampler, prng, null);
        rejectsDraw("an output array of the wrong length", sampler, prng, new double[2]);
    }

    private static void rejectsBuild(String what, double[] alpha) {
        try {
            DirichletSampler.of(alpha);
            fail("of accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    private static void rejectsDraw(String what, DirichletSampler sampler, PseudoRandom prng,
            double[] proportions) {
        try {
            sampler.sample(prng, proportions);
            fail("sample accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
