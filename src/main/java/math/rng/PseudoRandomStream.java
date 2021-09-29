/*
 * Copyright 2021 Stefan Zobel
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
}
