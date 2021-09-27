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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code double} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
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
     *            the least value that can be produced
     * @param max
     *            the largest value (inclusive) that can be produced
     *
     * @return a stream of pseudorandomly chosen {@code double} values, each
     *         between the specified min (inclusive) and the specified max
     *         (inclusive)
     *
     * @throws IllegalArgumentException
     *             if {@code min} is not finite, or {@code max} is not finite,
     *             or {@code min} is greater than {@code max}
     */
    DoubleStream doubles(double min, double max);
}
