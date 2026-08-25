package math.stats.gof;

import java.util.Arrays;
import java.util.stream.IntStream;

import math.distribution.ContinuousDistribution;
import math.distribution.Exponential;
import math.distribution.LogNormal;
import math.distribution.Normal;
import math.rng.BitMix;
import math.rng.PseudoRandom;
import math.rng.XorShiftRot256StarStar;
import math.stats.mle.MLE;
import math.stats.mle.ParExponential;
import math.stats.mle.ParLogNormal;
import math.stats.mle.ParNormal;

/**
 * The null distribution of the Kolmogorov-Smirnov statistic when the
 * distribution it is measured against was fitted to the same sample, which
 * makes the ordinary null distribution far too generous.
 * <p>
 * The fit is the maximum likelihood one that {@link math.stats.mle.MLE}
 * supplies, which for the normal family means the standard deviation with
 * {@code n} in the denominator. <b>The published Lilliefors tables use the
 * one with {@code n - 1}</b>, a larger number, so their critical values are
 * smaller than the ones simulated here by exactly
 * {@code sqrt(n / (n - 1))}; rescaled by that factor they agree. What makes
 * a test correct is that the estimator which produced the statistic is the
 * estimator the null distribution was drawn with, and here that is the one
 * the caller would have used anyway.
 * <p>
 * https://en.wikipedia.org/wiki/Lilliefors_test
 *
 * @since 1.5.3
 */
public final class Lilliefors {

    /**
     * The families this null distribution is available for.
     * <p>
     * All three are location-scale families, {@link #LOGNORMAL} after taking
     * logarithms, and their maximum likelihood estimators carry that structure:
     * scaling the sample scales the fit with it. The standardized residuals,
     * and therefore the statistic, are left where they were, so the null
     * distribution depends on the sample size alone and can be simulated
     * without knowing the parameters. Gamma, Beta and Weibull are absent
     * because a shape parameter breaks exactly that argument.
     */
    public enum Family {
        /** The normal distribution, fitted in location and scale. */
        NORMAL,
        /** The exponential distribution, fitted in scale. */
        EXPONENTIAL,
        /** The log-normal distribution, the normal one on the logarithms. */
        LOGNORMAL
    }

    /** The number of replications the short form of the test uses. */
    public static final int DEFAULT_REPLICATIONS = 10_000;

    /**
     * The seed the short form of the test uses. Its value is arbitrary; what
     * matters is that it does not change between two runs of the same code, so
     * that a p-value near a decision threshold cannot move across it on its
     * own.
     */
    public static final long DEFAULT_SEED = 42L;

    /** The smallest sample the test says anything about. */
    public static final int MINIMUM_SAMPLE = 3;

    /** The odd golden-ratio increment SplitMix64 spaces its seeds by. */
    private static final long GOLDEN_GAMMA = 0x9E3779B97F4A7C15L;

    /**
     * Replications times sample size from which the replications are spread
     * over the available cores. The threads were measured to pay for
     * themselves well below this -- 2000 replications of five observations
     * already run eight times faster -- but the sequential cost there is a
     * millisecond either way, so the threshold sits where the saving starts
     * to be worth having rather than where it starts.
     */
    private static final long PARALLEL_THRESHOLD = 10_000L;

    /**
     * The Kolmogorov-Smirnov distance between {@code x} and the member of
     * {@code family} fitted to {@code x} by maximum likelihood.
     * <p>
     * The statistic lives here rather than with the tests because the
     * simulation below has to measure the same thing on its own samples, and
     * two implementations of one definition would be free to drift apart.
     *
     * @param family
     *            the family that is fitted and tested against
     * @param x
     *            the sample, at least {@link #MINIMUM_SAMPLE} finite
     *            observations, strictly positive for
     *            {@link Family#EXPONENTIAL} and {@link Family#LOGNORMAL}
     * @return {@code sup |F_n(t) - F(t)|} against the fitted member
     * @throws IllegalArgumentException
     *             if {@code family} or {@code x} is {@code null}, if {@code x}
     *             is shorter than {@link #MINIMUM_SAMPLE}, or if the fit
     *             rejects an observation
     */
    public static double statistic(Family family, double[] x) {
        requireFamily(family);
        if (x == null) {
            throw new IllegalArgumentException("x must not be null");
        }
        requireSampleSize(x.length);
        return distance(fit(family, x), x);
    }

    /**
     * The simulated upper tail {@code P[D_n >= d]} of the fitted statistic.
     * <p>
     * There is no closed form, so the null is drawn: {@code replications}
     * samples of size {@code n} from the standard member of the family, each
     * fitted and measured the way {@link #statistic(Family, double[])} measures
     * the real one. The tail is
     * {@code (1 + #{D_b >= d}) / (replications + 1)}, which is never zero and
     * keeps the test valid at any number of replications. Its own uncertainty
     * is {@code sqrt(p (1 - p) / replications)}.
     * <p>
     * The answer is reproducible from {@code seed} alone: every replication
     * derives its generator from the replication index rather than from shared
     * state, so it does not matter whether the work was spread over threads.
     *
     * @param family
     *            the family that was fitted
     * @param n
     *            the sample size, at least {@link #MINIMUM_SAMPLE}
     * @param d
     *            the value of the statistic
     * @param replications
     *            how many samples to draw, {@code 1} or more
     * @param seed
     *            the seed the replications are derived from
     * @return {@code P[D_n >= d]} as estimated from the replications
     * @throws IllegalArgumentException
     *             if {@code family} is {@code null}, if {@code n} is below
     *             {@link #MINIMUM_SAMPLE}, if {@code replications} is not
     *             strictly positive, or if {@code d} is {@code NaN}
     */
    public static double barF(Family family, int n, double d, int replications, long seed) {
        requireFamily(family);
        requireSampleSize(n);
        if (replications < 1) {
            throw new IllegalArgumentException("replications must be strictly positive : " + replications);
        }
        if (Double.isNaN(d)) {
            throw new IllegalArgumentException("d must not be NaN");
        }
        if (d <= 0.0) {
            return 1.0;
        }
        IntStream replicationIndices = IntStream.range(0, replications);
        if ((long) replications * n >= PARALLEL_THRESHOLD) {
            replicationIndices = replicationIndices.parallel();
        }
        long atLeast = replicationIndices.filter(b -> replicate(family, n, seed, b) >= d).count();
        return (1.0 + atLeast) / (replications + 1.0);
    }

    /** One replication: draw, fit, measure. */
    private static double replicate(Family family, int n, long masterSeed, int replication) {
        // the generator comes from the replication index and never from shared
        // state, which is what makes the run reproducible whatever the threads
        // decide to do
        PseudoRandom rng = new XorShiftRot256StarStar(
                BitMix.staffordMix13(masterSeed + (replication + 1L) * GOLDEN_GAMMA));
        double[] sample = new double[n];
        switch (family) {
        case EXPONENTIAL:
            for (int i = 0; i < n; i++) {
                // the inverse distribution function of the unit exponential,
                // written so that a draw close to zero keeps its digits
                sample[i] = -Math.log1p(-rng.nextDouble());
            }
            break;
        case LOGNORMAL:
            for (int i = 0; i < n; i++) {
                sample[i] = Math.exp(rng.nextGaussian());
            }
            break;
        default:
            for (int i = 0; i < n; i++) {
                sample[i] = rng.nextGaussian();
            }
            break;
        }
        return distance(fit(family, sample), sample);
    }

    private static ContinuousDistribution fit(Family family, double[] x) {
        switch (family) {
        case EXPONENTIAL:
            ParExponential exponential = MLE.getExponentialMLE(x);
            return new Exponential(exponential.lambda);
        case LOGNORMAL:
            ParLogNormal logNormal = MLE.getLogNormalMLE(x);
            return new LogNormal(logNormal.mu, logNormal.sigma);
        default:
            ParNormal normal = MLE.getNormalMLE(x);
            return new Normal(normal.mean, normal.stdDev);
        }
    }

    /** The largest gap between the empirical distribution function and {@code fitted}. */
    private static double distance(ContinuousDistribution fitted, double[] x) {
        int n = x.length;
        double[] uniform = new double[n];
        for (int i = 0; i < n; i++) {
            uniform[i] = fitted.cdf(x[i]);
        }
        Arrays.sort(uniform);
        double dPlus = 0.0;
        double dMinus = 0.0;
        for (int i = 0; i < n; i++) {
            dPlus = Math.max(dPlus, (i + 1.0) / n - uniform[i]);
            dMinus = Math.max(dMinus, uniform[i] - i / (double) n);
        }
        return Math.max(dPlus, dMinus);
    }

    private static void requireFamily(Family family) {
        if (family == null) {
            throw new IllegalArgumentException("family must not be null");
        }
    }

    private static void requireSampleSize(int n) {
        if (n < MINIMUM_SAMPLE) {
            throw new IllegalArgumentException("at least " + MINIMUM_SAMPLE
                    + " observations are needed for a fitted test, got " + n);
        }
    }

    private Lilliefors() {
        throw new AssertionError();
    }
}
