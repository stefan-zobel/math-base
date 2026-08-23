/*
 * Copyright 2013, 2024 Stefan Zobel
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
     * @throws IllegalArgumentException if {@code stdDeviation <= 0.0}
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
}
