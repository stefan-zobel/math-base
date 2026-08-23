package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.Test;

import math.cern.BetaFun;

/**
 * {@code beta(alpha, beta)} used to draw two {@code Gamma(k, 1)} variates and
 * divide. The identity {@code X / (X + Y)} is exact, but forming {@code X} and
 * {@code Y} first is not: below a shape of about {@code 0.01} a gamma variate
 * spans hundreds of decades and falls off the bottom of the {@code double}
 * range. One underflowed draw dragged the quotient to an exact zero where the
 * true value was representable -- the rate of exact zeros came out at twice the
 * rate the distribution has -- and two underflowed draws made it {@code 0/0},
 * which is where 22 per cent of {@code NaN} at {@code alpha = beta = 0.001}
 * came from. The sampler now takes the difference of the logarithms instead.
 * <p>
 * The oracle throughout is {@link BetaFun#incompleteBeta(double, double,
 * double)}, the regularized incomplete beta function, which is the exact
 * distribution function.
 * <p>
 * What this does <em>not</em> repair, and what the last test pins down so that
 * nobody tries: for small shapes a large share of the result is exactly
 * {@code 1.0}. A {@code double} resolves {@code 4.9e-324} above zero but only
 * {@code 1.1e-16} below one, and the distribution genuinely puts that much mass
 * within half an ulp of one. Those exact ones are the correctly rounded answer.
 */
public class BetaSpliteratorTest {

    /** half an ulp of one: everything above this rounds to 1.0 */
    private static final double HALF_ULP = Math.ulp(1.0) / 2.0;

    private static double[] draw(long seed, long count, double alpha, double beta) {
        return new XorShiftRot256StarStar(seed).beta(count, alpha, beta).toArray();
    }

    /**
     * The outright bug: both gamma draws underflowing gave {@code 0/0}. This
     * failed before the change with 227 NaN at {@code 0.005} and 90225 at
     * {@code 0.001} out of 400000.
     */
    @Test
    public void testSmallShapesProduceNoNaN() {
        for (double shape : new double[] { 0.05, 0.01, 0.005, 0.001 }) {
            double[] sample = draw(99L, 400000L, shape, shape);
            int nan = 0;
            for (int i = 0; i < sample.length; ++i) {
                if (Double.isNaN(sample[i])) {
                    ++nan;
                }
            }
            assertEquals("NaN at alpha = beta = " + shape, 0, nan);
        }
    }

    /**
     * A single underflowed gamma draw used to force the quotient to zero, which
     * roughly doubled the share of exact zeros: 0.000595 measured against a
     * true 0.000292 at a shape of {@code 0.01}, and 0.023598 against 0.012091
     * at {@code 0.005}. The tolerance is a quarter of the true share, which the
     * repaired sampler stays inside over every seed tried (worst deviation 20
     * per cent of the share at {@code 0.01}, where only about 117 of 400000
     * draws are expected to underflow at all).
     */
    @Test
    public void testTheZeroRateMatchesTheExactCdf() {
        for (double shape : new double[] { 0.01, 0.005, 0.001 }) {
            double[] sample = draw(99L, 400000L, shape, shape);
            int zeros = 0;
            for (int i = 0; i < sample.length; ++i) {
                if (sample[i] == 0.0) {
                    ++zeros;
                }
            }
            double truth = BetaFun.incompleteBeta(shape, shape, Double.MIN_VALUE);
            double share = zeros / (double) sample.length;
            assertEquals("share of exact zeros at alpha = beta = " + shape, truth, share, 0.25 * truth);
        }
    }

    /**
     * The part of the lower tail that a {@code double} can hold. Between
     * {@code Double.MIN_VALUE} and {@code 1e-300} the distribution puts about
     * 208 of a million draws at a shape of {@code 0.01}; the old sampler
     * delivered fewer than half of them, because the quotient needed a
     * numerator far below the smallest positive {@code double}. The tolerance
     * is 30 per cent, against a measured spread of 13 per cent over five seeds.
     */
    @Test
    public void testTheLowerTailIsRepresentable() {
        final double shape = 0.01;
        final int n = 1000000;
        double[] sample = draw(3L, n, shape, shape);
        int inside = 0;
        for (int i = 0; i < sample.length; ++i) {
            if (sample[i] > 0.0 && sample[i] < 1.0e-300) {
                ++inside;
            }
        }
        double expected = n * (BetaFun.incompleteBeta(shape, shape, 1.0e-300)
                - BetaFun.incompleteBeta(shape, shape, Double.MIN_VALUE));
        assertEquals("draws strictly inside (MIN_VALUE, 1e-300)", expected, inside, 0.3 * expected);
    }

    /**
     * The guard that the repair did not move the distribution: a
     * Kolmogorov-Smirnov statistic against the exact distribution function, for
     * shapes where the old arithmetic worked. This has to pass before and after
     * the change. The threshold is well above the one per cent critical value
     * of {@code 1.63 / sqrt(n) = 0.0073} because it is a regression guard and
     * not a hypothesis test; the measured statistics run from 0.0024 to 0.0043.
     */
    @Test
    public void testTheDistributionIsUnchangedForOrdinaryShapes() {
        double[][] shapes = { { 2.0, 3.0 }, { 0.5, 0.5 }, { 1.0, 1.0 }, { 0.25, 0.25 }, { 5.0, 0.5 },
                { 0.3, 7.0 } };
        for (int s = 0; s < shapes.length; ++s) {
            double alpha = shapes[s][0];
            double beta = shapes[s][1];
            double[] sample = draw(7L, 50000L, alpha, beta);
            Arrays.sort(sample);
            double d = 0.0;
            for (int i = 0; i < sample.length; ++i) {
                double cdf = BetaFun.incompleteBeta(alpha, beta, sample[i]);
                d = Math.max(d, Math.abs(cdf - (i + 1.0) / sample.length));
                d = Math.max(d, Math.abs(cdf - i / (double) sample.length));
            }
            assertTrue("Kolmogorov-Smirnov distance " + d + " for beta(" + alpha + ", " + beta + ")", d < 0.008);
        }
    }

    /**
     * The exact ones are correct and must stay. Their share matches the true
     * probability of landing within half an ulp of one, so a future reader who
     * takes them for the same defect as the zeros can see here that they are
     * not. Passes before and after the change alike. The tolerance is eight per
     * cent of the share against a measured spread of 3.6 per cent over four
     * seeds; at a shape of {@code 0.125} only about 2000 of the 400000 draws
     * round to one, so a couple of per cent is counting noise.
     */
    @Test
    public void testExactOnesAreCorrectRounding() {
        for (double shape : new double[] { 0.125, 0.05, 0.01 }) {
            double[] sample = draw(99L, 400000L, shape, shape);
            int ones = 0;
            for (int i = 0; i < sample.length; ++i) {
                if (sample[i] == 1.0) {
                    ++ones;
                }
            }
            double truth = 1.0 - BetaFun.incompleteBeta(shape, shape, 1.0 - HALF_ULP);
            double share = ones / (double) sample.length;
            assertEquals("share rounding to 1.0 at alpha = beta = " + shape, truth, share, 0.08 * truth);
        }
    }

    /**
     * {@code fisherF} used to reach its value through a beta variate and then
     * divide by {@code d1 - d1 * y}, which is zero as soon as {@code y} rounds
     * to one. It now takes the same difference of logarithms the beta sampler
     * does. The means are checked against {@code d2 / (d2 - 2)}, which exists
     * only above two denominator degrees of freedom.
     */
    @Test
    public void testFisherFIsFiniteAndCentered() {
        int[][] dfs = { { 1, 1 }, { 1, 5 }, { 5, 2 }, { 10, 20 }, { 100, 100 } };
        for (int i = 0; i < dfs.length; ++i) {
            int d1 = dfs[i][0];
            int d2 = dfs[i][1];
            double[] sample = new XorShiftRot256StarStar(99L).fisherF(200000L, d1, d2).toArray();
            double sum = 0.0;
            for (int j = 0; j < sample.length; ++j) {
                assertTrue("non-finite F(" + d1 + ", " + d2 + ") variate " + sample[j],
                        !Double.isNaN(sample[j]) && !Double.isInfinite(sample[j]));
                sum += sample[j];
            }
            if (d2 > 2) {
                double mean = d2 / (d2 - 2.0);
                assertEquals("mean of F(" + d1 + ", " + d2 + ")", mean, sum / sample.length, 0.05 * mean);
            }
        }
    }
}
