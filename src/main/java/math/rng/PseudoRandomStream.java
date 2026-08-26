/*
 * Copyright 2021, 2026 Stefan Zobel
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

import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

/**
 * Each method that returns a stream produces a stream of values each of which
 * is chosen in the same manner as for a method that returns a single
 * pseudorandomly chosen value.
 */
public interface PseudoRandomStream {

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code int} values.
     *
     * @return a stream of pseudorandomly chosen {@code int} values
     */
    IntStream ints();

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code int} values.
     *
     * @param streamSize
     *            the number of values to generate
     *
     * @return a stream of pseudorandomly chosen {@code int} values
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero
     */
    IntStream ints(long streamSize);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code int} values, where each value is between the
     * specified min (inclusive) and the specified max (inclusive).
     *
     * @param streamSize
     *            the number of values to generate
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code int} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code min} is
     *             greater than {@code max}
     */
    IntStream ints(long streamSize, int min, int max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code int} values, where each value is between the specified min
     * (inclusive) and the specified max (inclusive).
     *
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code int} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
     *
     * @throws IllegalArgumentException
     *             if {@code min} is greater than {@code max}
     */
    IntStream ints(int min, int max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code long} values.
     *
     * @return a stream of pseudorandomly chosen {@code long} values
     */
    LongStream longs();

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code long} values.
     *
     * @param streamSize
     *            the number of values to generate
     *
     * @return a stream of pseudorandomly chosen {@code long} values
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero
     */
    LongStream longs(long streamSize);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code long} values, where each value is between
     * the specified min (inclusive) and the specified max (inclusive).
     *
     * @param streamSize
     *            the number of values to generate
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code long} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code min} is
     *             greater than {@code max}
     */
    LongStream longs(long streamSize, long min, long max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code long} values, where each value is between the specified min
     * (inclusive) and the specified max (inclusive).
     *
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code long} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
     *
     * @throws IllegalArgumentException
     *             if {@code min} is greater than {@code max}
     */
    LongStream longs(long min, long max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code double} values.
     *
     * @return a stream of pseudorandomly chosen {@code double} values
     */
    DoubleStream doubles();

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code double} values.
     *
     * @param streamSize
     *            the number of values to generate
     *
     * @return a stream of pseudorandomly chosen {@code double} values
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero
     */
    DoubleStream doubles(long streamSize);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen {@code double} values, where each value is between
     * the specified min (inclusive) and the specified max (inclusive).
     *
     * @param streamSize
     *            the number of values to generate
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code double} values, each
     *         between the specified min and the specified max
     *
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code min} is
     *             not finite, or {@code max} is not finite, or {@code min} is
     *             greater than {@code max}
     */
    DoubleStream doubles(long streamSize, double min, double max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * {@code double} values, where each value is between the specified min
     * (inclusive) and the specified max (inclusive).
     *
     * @param min
     *            the smallest value that can be produced
     * @param max
     *            the largest value that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code double} values, each
     *         between the specified min and the specified max
     *
     * @throws IllegalArgumentException
     *             if {@code min} is not finite, or {@code max} is not finite,
     *             or {@code min} is greater than {@code max}
     */
    DoubleStream doubles(double min, double max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen normal
     * variates with expectation {@code mu} and standard deviation
     * {@code sigma}.
     * 
     * @param mu
     *            the expectation of the normal variate
     * @param sigma
     *            the standard deviation of the normal variate
     * @return a stream of pseudorandomly chosen normal variates with the
     *         specified expecation and standard deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code sigma} is not greater than zero
     */
    DoubleStream normal(double mu, double sigma);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen normal variates with expectation {@code mu} and
     * standard deviation {@code sigma}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param mu
     *            the expectation of the normal variate
     * @param sigma
     *            the standard deviation of the normal variate
     * @return a stream of pseudorandomly chosen normal variates with the
     *         specified expecation and standard deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code sigma} is
     *             not greater than zero
     */
    DoubleStream normal(long streamSize, double mu, double sigma);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Cauchy-distributed variates with location parameter {@code location} and
     * scale parameter {@code scale}.
     * 
     * @param location
     *            the location parameter of the Cauchy distribution
     * @param scale
     *            the scale parameter of the Cauchy distribution
     * @return a stream of pseudorandomly chosen Cauchy-distributed variates
     *         with the specified location and scale
     * @throws IllegalArgumentException
     *             if {@code scale} is not greater than zero
     */
    DoubleStream cauchy(double location, double scale);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Cauchy-distributed variates with location parameter
     * {@code location} and scale parameter {@code scale}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param location
     *            the location parameter of the Cauchy distribution
     * @param scale
     *            the scale parameter of the Cauchy distribution
     * @return a stream of pseudorandomly chosen Cauchy-distributed variates
     *         with the specified location and scale
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code scale} is
     *             not greater than zero
     */
    DoubleStream cauchy(long streamSize, double location, double scale);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * exponentially distributed variates with rate parameter {@code lambda}.
     * 
     * @param lambda
     *            the rate parameter {@code lambda} of the Exponential
     *            distribution
     * @return a stream of pseudorandomly chosen exponentially distributed
     *         variates with parameter {@code lambda}
     * @throws IllegalArgumentException
     *             if {@code lambda} is not greater than zero
     */
    DoubleStream exponential(double lambda);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * exponentially distributed variates with rate parameter {@code lambda}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param lambda
     *            the rate parameter {@code lambda} of the Exponential
     *            distribution
     * @return a stream of pseudorandomly chosen exponentially distributed
     *         variates with parameter {@code lambda}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code lambda} is
     *             not greater than zero
     */
    DoubleStream exponential(long streamSize, double lambda);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Gamma-distributed variates with shape parameter {@code k} and scale
     * parameter {@code theta}.
     * 
     * @param k
     *            the shape parameter of the Gamma distribution
     * @param theta
     *            the scale parameter of the Gamma distribution
     * @return a stream of pseudorandomly chosen Gamma-distributed variates with
     *         with shape parameter {@code k} and scale parameter {@code theta}
     * @throws IllegalArgumentException
     *             if {@code k} is not greater than zero, or {@code theta} is
     *             not greater than zero
     */
    DoubleStream gamma(double k, double theta);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Gamma-distributed variates with shape parameter
     * {@code k} and scale parameter {@code theta}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param k
     *            the shape parameter of the Gamma distribution
     * @param theta
     *            the scale parameter of the Gamma distribution
     * @return a stream of pseudorandomly chosen Gamma-distributed variates with
     *         with shape parameter {@code k} and scale parameter {@code theta}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code k} is not
     *             greater than zero, or {@code theta} is not greater than zero
     */
    DoubleStream gamma(long streamSize, double k, double theta);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Beta-distributed variates with shape parameters {@code alpha} and
     * {@code beta}.
     * 
     * @param alpha
     *            the shape parameter {@code alpha} of the Beta distribution
     * @param beta
     *            the shape parameter {@code beta} of the Beta distribution
     * @return a stream of pseudorandomly chosen Beta-distributed variates with
     *         shape parameters {@code alpha} and {@code beta}
     * @throws IllegalArgumentException
     *             if {@code alpha} is not greater than zero, or {@code beta} is
     *             not greater than zero
     */
    DoubleStream beta(double alpha, double beta);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Beta-distributed variates with shape parameters
     * {@code alpha} and {@code beta}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param alpha
     *            the shape parameter {@code alpha} of the Beta distribution
     * @param beta
     *            the shape parameter {@code beta} of the Beta distribution
     * @return a stream of pseudorandomly chosen Beta-distributed variates with
     *         shape parameters {@code alpha} and {@code beta}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code alpha} is
     *             not greater than zero, or {@code beta} is not greater than
     *             zero
     */
    DoubleStream beta(long streamSize, double alpha, double beta);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Chi-squared distributed variates with {@code k} degrees of freedom.
     * 
     * @param k
     *            the number of degrees of freedom
     * @return a stream of pseudorandomly chosen Chi-squared distributed
     *         variates with {@code k} degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code k} is not greater than zero
     */
    DoubleStream chiSquare(double k);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Chi-squared distributed variates with {@code k}
     * degrees of freedom.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param k
     *            the number of degrees of freedom
     * @return a stream of pseudorandomly chosen Chi-squared distributed
     *         variates with {@code k} degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code k} is not
     *             greater than zero
     */
    DoubleStream chiSquare(long streamSize, double k);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen FisherF
     * distributed variates with {@code numeratorDF} and {@code denominatorDF}
     * degrees of freedom.
     * 
     * @param numeratorDF
     *            degrees of freedom of the numerator
     * @param denominatorDF
     *            degrees of freedom of the denominator
     * @return a stream of pseudorandomly chosen FisherF distributed variates
     *         with {@code numeratorDF} and {@code denominatorDF} degrees of
     *         freedom
     * @throws IllegalArgumentException
     *             if {@code numeratorDF} is less than one, or
     *             {@code denominatorDF} is less than one
     */
    DoubleStream fisherF(int numeratorDF, int denominatorDF);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen FisherF distributed variates with
     * {@code numeratorDF} and {@code denominatorDF} degrees of freedom.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param numeratorDF
     *            degrees of freedom of the numerator
     * @param denominatorDF
     *            degrees of freedom of the denominator
     * @return a stream of pseudorandomly chosen FisherF distributed variates
     *         with {@code numeratorDF} and {@code denominatorDF} degrees of
     *         freedom
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or
     *             {@code numeratorDF} is less than one, or
     *             {@code denominatorDF} is less than one
     */
    DoubleStream fisherF(long streamSize, int numeratorDF, int denominatorDF);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Log-normal variates whose natural logarithm has expectation {@code mu}
     * and standard deviation {@code sigma}.
     * 
     * @param mu
     *            the expectation of the natural logarithm of the Log-normal
     *            variate
     * @param sigma
     *            the standard deviation of the natural logarithm of the
     *            Log-normal variate
     * @return a stream of pseudorandomly chosen Log-normal variates whose
     *         natural logarithm has the specified expecation and standard
     *         deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code sigma} is not greater than zero
     */
    DoubleStream logNormal(double mu, double sigma);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Log-normal variates whose natural logarithm has
     * expectation {@code mu} and standard deviation {@code sigma}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param mu
     *            the expectation of the natural logarithm of the Log-normal
     *            variate
     * @param sigma
     *            the standard deviation of the natural logarithm of the
     *            Log-normal variate
     * @return a stream of pseudorandomly chosen Log-normal variates whose
     *         natural logarithm has the specified expecation and standard
     *         deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code sigma} is
     *             not greater than zero
     */
    DoubleStream logNormal(long streamSize, double mu, double sigma);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen StudentT
     * distributed variates with {@code df} degrees of freedom.
     * 
     * @param df
     *            degrees of freedom
     * @return a stream of pseudorandomly chosen StudentT distributed variates
     *         with {@code df} degrees of freedom
     * 
     * @throws IllegalArgumentException
     *             if {@code df} is not greater than zero
     */
    DoubleStream studentT(double df);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen StudentT distributed variates with {@code df}
     * degrees of freedom.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param df
     *            degrees of freedom
     * @return a stream of pseudorandomly chosen StudentT distributed variates
     *         with {@code df} degrees of freedom
     * 
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code df} is not
     *             greater than zero
     */
    DoubleStream studentT(long streamSize, double df);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen Weibull
     * distributed variates with parameters {@code scale} and {@code shape}.
     * 
     * @param scale
     *            the scale parameter of the Weibull distribution
     * @param shape
     *            the shape parameter of the Weibull distribution
     * @return a stream of pseudorandomly chosen Weibull distributed variates
     *         with parameters {@code scale} and {@code shape}
     * 
     * @throws IllegalArgumentException
     *             if {@code scale} and {@code shape} are not greater than zero
     */
    DoubleStream weibull(double scale, double shape);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Weibull distributed variates with parameters
     * {@code scale} and {@code shape}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param scale
     *            the scale parameter of the Weibull distribution
     * @param shape
     *            the shape parameter of the Weibull distribution
     * @return a stream of pseudorandomly chosen Weibull distributed variates
     *         with parameters {@code scale} and {@code shape}
     * 
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code scale} and
     *             {@code shape} are not greater than zero
     */
    DoubleStream weibull(long streamSize, double scale, double shape);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * truncated standard normal random variates on the domain
     * {@code (min, max)}.
     * 
     * @param min
     *            the lower bound for truncation
     * @param max
     *            the upper bound for truncation
     * @return a stream of pseudorandomly chosen truncated standard normal
     *         samples on the domain {@code (min, max)}
     * @throws IllegalArgumentException
     *             if {@code min} is not finite, or {@code max} is not finite,
     *             or {@code min} is greater than {@code max}
     */
    DoubleStream truncatedStandardNormal(double min, double max);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen truncated standard normal random variates on the
     * domain {@code (min, max)}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param min
     *            the lower bound for truncation
     * @param max
     *            the upper bound for truncation
     * @return a stream of pseudorandomly chosen truncated standard normal
     *         samples on the domain {@code (min, max)}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code min} is
     *             not finite, or {@code max} is not finite, or {@code min} is
     *             greater than {@code max}
     */
    DoubleStream truncatedStandardNormal(long streamSize, double min, double max);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen LeCun
     * normal variates with standard deviation {@code sigma} (usually
     * {@code Math.sqrt(1.0 / fan_in)} in a deep learning setting).
     * 
     * @param sigma
     *            the standard deviation of the LeCun normal variate
     * @return a stream of pseudorandomly chosen LeCun normal variates with the
     *         specified standard deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code sigma} is not greater than zero
     */
    DoubleStream leCunNormal(double sigma);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen LeCun normal variates with standard deviation
     * {@code sigma} (usually {@code Math.sqrt(1.0 / fan_in)} in a deep learning
     * setting).
     * 
     * @param streamSize
     *            the number of values to generate
     * @param sigma
     *            the standard deviation of the LeCun normal variate
     * @return a stream of pseudorandomly chosen LeCun normal variates with the
     *         specified standard deviation
     * 
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code sigma} is
     *             not greater than zero
     */
    DoubleStream leCunNormal(long streamSize, double sigma);

    /**
     * Returns an effectively unlimited stream of pseudorandomly chosen
     * Inverse-gamma distributed variates with shape parameter {@code alpha} and
     * scale parameter {@code beta}.
     * 
     * @param alpha
     *            the shape parameter of the Inverse-gamma distribution
     * @param beta
     *            the scale parameter of the Inverse-gamma distribution
     * @return a stream of pseudorandomly chosen Inverse-gamma distributed
     *         variates with with shape parameter {@code alpha} and scale
     *         parameter {@code beta}
     * @throws IllegalArgumentException
     *             if {@code alpha} is not greater than zero, or {@code beta} is
     *             not greater than zero
     */
    DoubleStream inverseGamma(double alpha, double beta);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * pseudorandomly chosen Inverse-gamma distributed variates with shape
     * parameter {@code alpha} and scale parameter {@code beta}.
     * 
     * @param streamSize
     *            the number of values to generate
     * @param alpha
     *            the shape parameter of the Inverse-gamma distribution
     * @param beta
     *            the scale parameter of the Inverse-gamma distribution
     * @return a stream of pseudorandomly chosen Inverse-gamma distributed
     *         variates with with shape parameter {@code alpha} and scale
     *         parameter {@code beta}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or {@code alpha} is
     *             not greater than zero, or {@code beta} is not greater than
     *             zero
     */
    DoubleStream inverseGamma(long streamSize, double alpha, double beta);

    /**
     * Returns an effectively unlimited stream of outcomes drawn from the
     * categorical distribution the given weights describe.
     * <p>
     * The weights need not sum to one; they are normalized, so counts serve as
     * well as probabilities. The table behind the stream is built once, not
     * once per value, and a caller drawing repeatedly from the same
     * distribution should build an {@link AliasTable} and use the overload that
     * takes one.
     *
     * @param weights
     *            the weight of each outcome, at least one of them, each finite
     *            and not negative, and not all zero
     * @return a stream of outcomes in {@code 0 .. weights.length - 1}
     * @throws IllegalArgumentException
     *             if {@code weights} is {@code null}, is empty, holds a value
     *             that is negative or not finite, or sums to zero
     */
    IntStream categorical(double[] weights);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * outcomes drawn from the categorical distribution the given weights
     * describe.
     *
     * @param streamSize
     *            the number of values to generate
     * @param weights
     *            the weight of each outcome, at least one of them, each finite
     *            and not negative, and not all zero
     * @return a stream of outcomes in {@code 0 .. weights.length - 1}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or if
     *             {@code weights} is {@code null}, is empty, holds a value that
     *             is negative or not finite, or sums to zero
     */
    IntStream categorical(long streamSize, double[] weights);

    /**
     * Returns an effectively unlimited stream of outcomes drawn from a table
     * that was built once and can be drawn from by any number of streams.
     *
     * @param table
     *            the prepared categorical distribution
     * @return a stream of outcomes in
     *         {@code 0 .. table.outcomes() - 1}
     * @throws IllegalArgumentException
     *             if {@code table} is {@code null}
     */
    IntStream categorical(AliasTable table);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * outcomes drawn from a prepared table.
     *
     * @param streamSize
     *            the number of values to generate
     * @param table
     *            the prepared categorical distribution
     * @return a stream of outcomes in
     *         {@code 0 .. table.outcomes() - 1}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or if {@code table}
     *             is {@code null}
     */
    IntStream categorical(long streamSize, AliasTable table);

    /**
     * Returns an effectively unlimited stream of Poisson distributed counts
     * with mean {@code lambda}.
     *
     * @param lambda
     *            the mean of the Poisson distribution, not negative and at most
     *            {@code 1e9}, above which a count need not fit an {@code int}
     * @return a stream of pseudorandomly chosen Poisson distributed counts with
     *         mean {@code lambda}
     * @throws IllegalArgumentException
     *             if {@code lambda} is negative, is not a number, or exceeds
     *             {@code 1e9}
     */
    IntStream poisson(double lambda);

    /**
     * Returns a stream producing the given {@code streamSize} number of Poisson
     * distributed counts with mean {@code lambda}.
     *
     * @param streamSize
     *            the number of values to generate
     * @param lambda
     *            the mean of the Poisson distribution, not negative and at most
     *            {@code 1e9}, above which a count need not fit an {@code int}
     * @return a stream of pseudorandomly chosen Poisson distributed counts with
     *         mean {@code lambda}
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or if {@code lambda}
     *             is negative, is not a number, or exceeds {@code 1e9}
     */
    IntStream poisson(long streamSize, double lambda);

    /**
     * Returns an effectively unlimited stream of binomial counts: the number of
     * successes in {@code n} independent trials that each succeed with
     * probability {@code p}.
     *
     * @param n
     *            the number of trials, not negative
     * @param p
     *            the probability of success in one trial, in {@code [0, 1]}
     * @return a stream of pseudorandomly chosen binomial counts
     * @throws IllegalArgumentException
     *             if {@code n} is negative, or {@code p} does not lie in
     *             {@code [0, 1]}
     */
    IntStream binomial(int n, double p);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * binomial counts.
     *
     * @param streamSize
     *            the number of values to generate
     * @param n
     *            the number of trials, not negative
     * @param p
     *            the probability of success in one trial, in {@code [0, 1]}
     * @return a stream of pseudorandomly chosen binomial counts
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, if {@code n} is
     *             negative, or if {@code p} does not lie in {@code [0, 1]}
     */
    IntStream binomial(long streamSize, int n, double p);

    /**
     * Returns an effectively unlimited stream of geometric counts: the number of
     * failures before the first success, each trial succeeding with probability
     * {@code p}.
     *
     * @param p
     *            the probability of success in one trial, in
     *            {@code [1e-7, 1]}. Below that bound a count need not fit an
     *            {@code int}
     * @return a stream of pseudorandomly chosen geometric counts
     * @throws IllegalArgumentException
     *             if {@code p} does not lie in {@code [1e-7, 1]} or is not a
     *             number
     */
    IntStream geometric(double p);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * geometric counts.
     *
     * @param streamSize
     *            the number of values to generate
     * @param p
     *            the probability of success in one trial, in
     *            {@code [1e-7, 1]}
     * @return a stream of pseudorandomly chosen geometric counts
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or if {@code p} does
     *             not lie in {@code [1e-7, 1]} or is not a number
     */
    IntStream geometric(long streamSize, double p);

    /**
     * Returns an effectively unlimited stream of negative binomial counts: the
     * number of failures before the {@code r}-th success, each trial succeeding
     * with probability {@code p}.
     *
     * @param r
     *            the number of successes to wait for, {@code 1} or more
     * @param p
     *            the probability of success in one trial, in {@code (0, 1]}
     * @return a stream of pseudorandomly chosen negative binomial counts
     * @throws IllegalArgumentException
     *             if {@code r} is less than one, if {@code p} does not lie in
     *             {@code (0, 1]} or is not a number, or if the mean
     *             {@code r (1-p) / p} exceeds {@code 1e9} and so need not fit an
     *             {@code int}
     */
    IntStream negativeBinomial(int r, double p);

    /**
     * Returns a stream producing the given {@code streamSize} number of negative
     * binomial counts.
     *
     * @param streamSize
     *            the number of values to generate
     * @param r
     *            the number of successes to wait for, {@code 1} or more
     * @param p
     *            the probability of success in one trial, in {@code (0, 1]}
     * @return a stream of pseudorandomly chosen negative binomial counts
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or for any of the
     *             reasons {@link #negativeBinomial(int, double)} throws
     */
    IntStream negativeBinomial(long streamSize, int r, double p);

    /**
     * Returns an effectively unlimited stream of hypergeometric counts: the
     * number of successes among {@code draws} taken without replacement from a
     * population of {@code population} that holds {@code successes} of them.
     * <p>
     * A draw costs time linear in the width of the support,
     * {@code min(draws, successes) - max(0, draws + successes - population)}:
     * measured, 51 nanoseconds where that width is 20 and 35 microseconds where
     * it is 20000.
     *
     * @param population
     *            the size of the population, not negative
     * @param successes
     *            how many of the population count as a success, in
     *            {@code [0, population]}
     * @param draws
     *            how many are taken, in {@code [0, population]}
     * @return a stream of pseudorandomly chosen hypergeometric counts
     * @throws IllegalArgumentException
     *             if {@code population} is negative, or {@code successes} or
     *             {@code draws} lies outside {@code [0, population]}
     */
    IntStream hypergeometric(int population, int successes, int draws);

    /**
     * Returns a stream producing the given {@code streamSize} number of
     * hypergeometric counts.
     *
     * @param streamSize
     *            the number of values to generate
     * @param population
     *            the size of the population, not negative
     * @param successes
     *            how many of the population count as a success, in
     *            {@code [0, population]}
     * @param draws
     *            how many are taken, in {@code [0, population]}
     * @return a stream of pseudorandomly chosen hypergeometric counts
     * @throws IllegalArgumentException
     *             if {@code streamSize} is less than zero, or for any of the
     *             reasons {@link #hypergeometric(int, int, int)} throws
     */
    IntStream hypergeometric(long streamSize, int population, int successes, int draws);
}
