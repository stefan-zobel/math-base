/*
 * Copyright 2013, 2026 Stefan Zobel
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
package math.rng;

/**
 * A generator of uniform pseudorandom values.
 * <p>
 * Ranges given as {@code [min, max)} are half-open, those given as
 * {@code [min, max]} include both ends. A seeded generator is reproducible: the
 * same seed yields the same sequence, and the same sequential stream.
 * <p>
 * Implementations are <b>not</b> thread-safe and must not be shared between
 * threads. Give each thread its own instance, or use
 * {@link SplittablePseudoRandom#split()} where it is available. A parallel
 * stream is safe -- each split draws from its own generator -- but its values
 * depend on the parallelism of the common
 * {@link java.util.concurrent.ForkJoinPool} and are therefore not
 * reproducible from the seed.
 */
public interface PseudoRandom extends PseudoRandomStream {

    /**
     * Returns the next pseudorandom {@code long}, uniform over all
     * 2<sup>64</sup> values.
     *
     * @return a pseudorandom {@code long}
     */
    long nextLong();

    /**
     * Returns the next pseudorandom {@code double} in {@code [0, 1)}.
     *
     * @return a pseudorandom {@code double} in {@code [0, 1)}
     */
    double nextDouble();

    /**
     * Returns the next pseudorandom {@code double} in {@code [min, max)}, or
     * {@code min} itself if the range is empty.
     *
     * @param min the smallest value that can be returned
     * @param max the exclusive upper bound
     * @return a pseudorandom {@code double} in {@code [min, max)}
     * @throws IllegalArgumentException if {@code min > max} or either bound is
     *             {@code NaN}
     */
    double nextDouble(double min, double max);

    /**
     * Returns the next standard normal {@code double}, drawn by
     * <a href="https://en.wikipedia.org/wiki/Marsaglia_polar_method">
     * Marsaglia's polar method</a>.
     *
     * @return a pseudorandom {@code double} with mean {@code 0.0} and standard
     *         deviation {@code 1.0}
     */
    double nextGaussian();

    /**
     * Returns the next normally distributed {@code double} with the given mean
     * and standard deviation.
     *
     * @param mean the mean of the distribution
     * @param stdDeviation the standard deviation of the distribution, positive
     * @return a pseudorandom {@code double} from
     *         {@code N(mean, stdDeviation^2)}
     * @throws IllegalArgumentException if {@code stdDeviation} is not
     *             greater than zero
     */
    double nextGaussian(double mean, double stdDeviation);

    /**
     * Returns the next pseudorandom {@code float} in {@code [0, 1)}.
     *
     * @return a pseudorandom {@code float} in {@code [0, 1)}
     */
    float nextFloat();

    /**
     * Returns the next pseudorandom {@code float} in {@code [min, max)}, or
     * {@code min} itself if the range is empty.
     *
     * @param min the smallest value that can be returned
     * @param max the exclusive upper bound
     * @return a pseudorandom {@code float} in {@code [min, max)}
     * @throws IllegalArgumentException if {@code min > max} or either bound is
     *             {@code NaN}
     */
    float nextFloat(float min, float max);

    /**
     * Returns the next pseudorandom {@code int}, uniform over all
     * 2<sup>32</sup> values.
     *
     * @return a pseudorandom {@code int}
     */
    int nextInt();

    /**
     * Fills the given array with pseudorandom bytes.
     *
     * @param bytes the array to fill
     */
    void nextBytes(byte[] bytes);

    /**
     * Returns {@code true} and {@code false} with equal probability.
     *
     * @return a pseudorandom {@code boolean}
     */
    boolean nextBoolean();

    /**
     * Returns a pseudorandom {@code long} in {@code [0, n)}.
     *
     * @param n the exclusive upper bound, positive
     * @return a pseudorandom {@code long} in {@code [0, n)}
     * @throws IllegalArgumentException if {@code n <= 0}
     */
    long nextLong(long n);

    /**
     * Returns a pseudorandom {@code int} in {@code [0, n)}.
     *
     * @param n the exclusive upper bound, positive
     * @return a pseudorandom {@code int} in {@code [0, n)}
     * @throws IllegalArgumentException if {@code n <= 0}
     */
    int nextInt(int n);

    /**
     * Returns a pseudorandom {@code int} in {@code [min, max]}.
     *
     * @param min the smallest value that can be returned
     * @param max the largest value that can be returned
     * @return a pseudorandom {@code int} in {@code [min, max]}
     * @throws IllegalArgumentException if {@code min > max}
     */
    int nextInt(int min, int max);

    /**
     * Returns a pseudorandom {@code long} in {@code [min, max]}. The whole
     * {@code long} range is a legal argument.
     *
     * @param min the smallest value that can be returned
     * @param max the largest value that can be returned
     * @return a pseudorandom {@code long} in {@code [min, max]}
     * @throws IllegalArgumentException if {@code min > max}
     */
    long nextLong(long min, long max);

    /**
     * Returns the top {@code bits} bits of a draw: {@code 0} for
     * {@code bits == 0}, a non-negative value for {@code bits} up to
     * {@code 31}, and an ordinary signed {@code int} for {@code bits == 32}.
     *
     * @param bits the number of random bits to return, in {@code [0, 32]}
     * @return an {@code int} carrying {@code bits} random bits
     * @throws IllegalArgumentException if {@code bits} lies outside
     *             {@code [0, 32]}
     */
    int next(int bits);

    /**
     * Fills the given array with the values of {@link #nextLong()}.
     *
     * @param longs the array to fill
     */
    void nextLongs(long[] longs);

    /**
     * Fills the given array with the values of {@link #nextDouble()}.
     *
     * @param doubles the array to fill
     */
    void nextDoubles(double[] doubles);

    /**
     * Returns {@code count} distinct pseudorandom {@code int} values from
     * {@code [min, max]}, in the order they were drawn.
     *
     * @param min the smallest value that can be returned
     * @param max the largest value that can be returned
     * @param count the number of distinct values to draw
     * @return {@code count} distinct values from {@code [min, max]}
     * @throws IllegalArgumentException if {@code min > max}, if {@code count}
     *             is negative, if the range holds fewer than {@code count}
     *             values, or if {@code [min, max]} is the whole {@code int}
     *             range, which is rejected outright
     */
    int[] intsSampledWithoutReplacement(int min, int max, int count);

    /**
     * Returns the name of the generator algorithm.
     *
     * @return the algorithm name
     */
    String getAlgorithm();

    /**
     * Returns a copy of the seed this generator was created from.
     *
     * @return a copy of the initial seed, or {@code null} if none was recorded
     */
    long[] getSeed();

    /**
     * Returns one pseudorandom value from the Cauchy distribution with the
     * given location and scale.
     * <p>
     * The same value {@link #cauchy(long, double, double)} would put first in
     * its stream, without building one.
     *
     * @param location the location of the distribution, its median
     * @param scale the scale of the distribution, strictly positive
     * @return a {@code Cauchy(location, scale)} variate
     * @throws IllegalArgumentException if {@code scale} is not greater than
     *             zero
     * @since 1.5.3
     */
    default double nextCauchy(double location, double scale) {
        CauchySpliterator.checkScale(scale);
        return CauchySpliterator.sample(this, location, scale);
    }

    /**
     * Returns one pseudorandom value from the exponential distribution with
     * the given rate.
     *
     * @param lambda the rate of the distribution, strictly positive. Its mean
     *            is {@code 1 / lambda}
     * @return an {@code Exponential(lambda)} variate
     * @throws IllegalArgumentException if {@code lambda} is not greater
     *             than zero
     * @since 1.5.3
     */
    default double nextExponential(double lambda) {
        ExponentialSpliterator.checkRate(lambda);
        return ExponentialSpliterator.sample(this, lambda);
    }

    /**
     * Returns one pseudorandom value from the gamma distribution with the
     * given shape and scale.
     *
     * @param k the shape of the distribution, strictly positive
     * @param theta the scale of the distribution, strictly positive. The mean
     *            is {@code k * theta}
     * @return a {@code Gamma(k, theta)} variate
     * @throws IllegalArgumentException if {@code k} or {@code theta} is not
     *             greater than zero
     * @since 1.5.3
     */
    default double nextGamma(double k, double theta) {
        GammaSpliterator.checkParameters(k, theta);
        return GammaSpliterator.sample(this, k, theta);
    }

    /**
     * Returns one pseudorandom value from the beta distribution with the given
     * shapes.
     * <p>
     * A beta variate is built from two gamma variates. {@link #beta(long,
     * double, double)} draws the second from an independent generator split off
     * this one, and refuses to run where that is impossible; a single draw takes
     * both from this generator in sequence, which is equally independent and
     * works for every generator. The two therefore do <b>not</b> return the same
     * value from the same seed.
     *
     * @param alpha the first shape of the distribution, strictly positive
     * @param beta the second shape of the distribution, strictly positive
     * @return a {@code Beta(alpha, beta)} variate in {@code [0, 1]}
     * @throws IllegalArgumentException if {@code alpha} or {@code beta} is
     *             not greater than zero
     * @since 1.5.3
     */
    default double nextBeta(double alpha, double beta) {
        BetaSpliterator.checkParameters(alpha, beta);
        return BetaSpliterator.sample(this, this, alpha, beta);
    }

    /**
     * Returns one pseudorandom value from the chi-squared distribution with the
     * given degrees of freedom.
     *
     * @param k the degrees of freedom, strictly positive. It need not be a
     *            whole number
     * @return a chi-squared variate with {@code k} degrees of freedom
     * @throws IllegalArgumentException if {@code k} is not greater than zero
     * @since 1.5.3
     */
    default double nextChiSquare(double k) {
        ChiSquareSpliterator.checkDegreesOfFreedom(k);
        return ChiSquareSpliterator.sample(this, k);
    }

    /**
     * Returns one pseudorandom value from the F distribution with the given
     * degrees of freedom.
     * <p>
     * As with {@link #nextBeta(double, double)} both underlying variates come
     * from this generator in sequence, where {@link #fisherF(long, int, int)}
     * splits a second generator off, so the two do not agree from the same seed.
     *
     * @param numeratorDF the numerator degrees of freedom, at least one
     * @param denominatorDF the denominator degrees of freedom, at least one
     * @return an {@code F(numeratorDF, denominatorDF)} variate
     * @throws IllegalArgumentException if either argument is smaller than one
     * @since 1.5.3
     */
    default double nextFisherF(int numeratorDF, int denominatorDF) {
        FisherFSpliterator.checkParameters(numeratorDF, denominatorDF);
        return FisherFSpliterator.sample(this, this, numeratorDF, denominatorDF);
    }

    /**
     * Returns one pseudorandom value from the log-normal distribution whose
     * logarithm has the given mean and standard deviation.
     *
     * @param mu the mean of the logarithm of the distribution
     * @param sigma the standard deviation of the logarithm of the
     *            distribution, strictly positive
     * @return a {@code LogNormal(mu, sigma)} variate, strictly positive
     * @throws IllegalArgumentException if {@code sigma} is not greater than
     *             zero
     * @since 1.5.3
     */
    default double nextLogNormal(double mu, double sigma) {
        LogNormalSpliterator.checkSigma(sigma);
        return LogNormalSpliterator.sample(this, mu, sigma);
    }

    /**
     * Returns one pseudorandom value from Student's t distribution with the
     * given degrees of freedom.
     *
     * @param df the degrees of freedom, strictly positive. It need not be a
     *            whole number
     * @return a {@code t(df)} variate
     * @throws IllegalArgumentException if {@code df} is not greater than zero
     * @since 1.5.3
     */
    default double nextStudentT(double df) {
        StudentTSpliterator.checkDegreesOfFreedom(df);
        return StudentTSpliterator.sample(this, df);
    }

    /**
     * Returns one pseudorandom value from the Weibull distribution with the
     * given scale and shape.
     *
     * @param scale the scale of the distribution, strictly positive
     * @param shape the shape of the distribution, strictly positive
     * @return a {@code Weibull(scale, shape)} variate
     * @throws IllegalArgumentException if {@code scale} or {@code shape} is
     *             not greater than zero
     * @since 1.5.3
     */
    default double nextWeibull(double scale, double shape) {
        WeibullSpliterator.checkParameters(scale, shape);
        return WeibullSpliterator.sample(this, scale, shape);
    }

    /**
     * Returns one pseudorandom value from the standard normal distribution
     * truncated to {@code (min, max)}.
     * <p>
     * The result lies strictly between the two bounds. For an untruncated
     * normal use {@link #nextGaussian(double, double)}, which is what a
     * {@code nextNormal} would be.
     *
     * @param min the lower bound, exclusive
     * @param max the upper bound, exclusive, strictly greater than {@code min}
     * @return a standard normal variate conditioned on {@code (min, max)}
     * @throws IllegalArgumentException if {@code min} is not smaller than
     *             {@code max}
     * @since 1.5.3
     */
    default double nextTruncatedStandardNormal(double min, double max) {
        return TruncatedNormalSpliterator.sampleFor(this, min, max);
    }

    /**
     * Returns one pseudorandom value from the LeCun normal distribution with
     * the given standard deviation.
     *
     * @param sigma the standard deviation of the distribution, strictly
     *            positive
     * @return a LeCun normal variate
     * @throws IllegalArgumentException if {@code sigma} is not greater than
     *             zero
     * @since 1.5.3
     */
    default double nextLeCunNormal(double sigma) {
        return LeCunNormalSpliterator.sampleFor(this, sigma);
    }

    /**
     * Returns one pseudorandom value from the inverse gamma distribution with
     * the given shape and scale.
     *
     * @param alpha the shape of the distribution, strictly positive
     * @param beta the scale of the distribution, strictly positive
     * @return an {@code InverseGamma(alpha, beta)} variate
     * @throws IllegalArgumentException if {@code alpha} or {@code beta} is
     *             not greater than zero
     * @since 1.5.3
     */
    default double nextInverseGamma(double alpha, double beta) {
        return InverseGammaSpliterator.sampleFor(this, alpha, beta);
    }

    /**
     * Returns one pseudorandom count from the Poisson distribution with the
     * given mean.
     *
     * @param lambda the mean of the distribution, strictly positive
     * @return a {@code Poisson(lambda)} count, zero or more
     * @throws IllegalArgumentException if {@code lambda} is not strictly
     *             positive or is too large to draw from
     * @since 1.5.3
     */
    default int nextPoisson(double lambda) {
        PoissonSpliterator.checkMean(lambda);
        return PoissonSpliterator.sample(this, lambda);
    }

    /**
     * Returns one pseudorandom count from the binomial distribution with the
     * given number of trials and success probability.
     *
     * @param n the number of trials, zero or more
     * @param p the success probability of a single trial, in {@code [0, 1]}
     * @return a {@code Binomial(n, p)} count, from zero to {@code n}
     * @throws IllegalArgumentException if {@code n} is negative or {@code p}
     *             lies outside {@code [0, 1]}
     * @since 1.5.3
     */
    default int nextBinomial(int n, double p) {
        BinomialSpliterator.checkParameters(n, p);
        return BinomialSpliterator.sample(this, n, p);
    }

    /**
     * Returns one pseudorandom count from the geometric distribution with the
     * given success probability: the number of failures before the first
     * success.
     *
     * @param p the success probability, in {@code (0, 1]}
     * @return a {@code Geometric(p)} count, zero or more
     * @throws IllegalArgumentException if {@code p} is not in {@code (0, 1]}
     * @since 1.5.3
     */
    default int nextGeometric(double p) {
        GeometricSpliterator.checkProbability(p);
        return GeometricSpliterator.sampleFor(this, p);
    }

    /**
     * Returns one pseudorandom count from the negative binomial distribution:
     * the number of failures before the {@code r}th success.
     *
     * @param r the number of successes to wait for, at least one
     * @param p the success probability of a single trial, in {@code (0, 1]}
     * @return a {@code NegativeBinomial(r, p)} count, zero or more
     * @throws IllegalArgumentException if {@code r} is smaller than one or
     *             {@code p} is not in {@code (0, 1]}
     * @since 1.5.3
     */
    default int nextNegativeBinomial(int r, double p) {
        NegativeBinomialSpliterator.checkParameters(r, p);
        return NegativeBinomialSpliterator.sampleFor(this, r, p);
    }

    /**
     * Returns one pseudorandom count from the hypergeometric distribution: the
     * number of successes among {@code draws} taken without replacement from a
     * population holding {@code successes} of them.
     *
     * @param population the size of the population, zero or more
     * @param successes the number of successes in it, from zero to
     *            {@code population}
     * @param draws the number of items drawn, from zero to {@code population}
     * @return the number of successes drawn
     * @throws IllegalArgumentException if any argument is negative, or if
     *             {@code successes} or {@code draws} exceeds {@code population}
     * @since 1.5.3
     */
    default int nextHypergeometric(int population, int successes, int draws) {
        HypergeometricSpliterator.checkParameters(population, successes, draws);
        return HypergeometricSpliterator.sample(this, population, successes, draws);
    }
}
