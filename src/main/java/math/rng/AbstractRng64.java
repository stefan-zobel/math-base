/*
 * Copyright 2013, 2021 Stefan Zobel
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

import java.util.Arrays;
import java.util.Spliterator;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.StreamSupport;

/**
 * Abstract base class for 64-bit pseudo RNGs.
 * <p>
 * Derived classes need to supply an implementation of {@link #nextLong()}.
 * </p>
 * This base class implementation is efficient for {@link #nextDouble()},
 * {@link #nextGaussian()} and {@link #nextBytes(byte[])} but somehow wasteful
 * for the other methods because it dissipates valuable random bits piled up in
 * the call to {@link #nextLong()} whenever less than {@code 33} random bits are
 * needed for the result type.
 */
public abstract class AbstractRng64 implements PseudoRandom {

    protected static final double DOUBLE_NORM = 1.0 / (1L << 53);
    protected static final float FLOAT_NORM = 1.0F / (1 << 24);

    private static final String BAD_RANGE = "max must be >= min";

    // not used, but a potential target for a store
    protected static byte unused;

    private long[] initialSeed;

    /** cache for the next gaussian */
    protected double nextGaussian = Double.NaN;

    @Override
    public abstract long nextLong();

    @Override
    public double nextDouble() {
        return (nextLong() >>> 11) * DOUBLE_NORM;
    }

    @Override
    public double nextDouble(double min, double max) {
        return min + (max - min) * nextDouble();
    }

    @Override
    public final double nextGaussian() {
        final double rndVal;
        if (Double.isNaN(nextGaussian)) {
            // Marsaglia's polar method
            double u1, u2, q;
            do {
                u1 = 2.0 * nextDouble() - 1.0; // between -1 and 1
                u2 = 2.0 * nextDouble() - 1.0; // between -1 and 1
                q = u1 * u1 + u2 * u2;
            } while (q >= 1 || q == 0.0);
            final double p = Math.sqrt(-2.0 * Math.log(q) / q);
            rndVal = u1 * p;
            nextGaussian = u2 * p;
        } else {
            rndVal = nextGaussian;
            nextGaussian = Double.NaN;
        }
        return rndVal;
    }

    @Override
    public double nextGaussian(double mean, double stdDeviation) {
        return mean + stdDeviation * nextGaussian();
    }

    @Override
    public float nextFloat() {
        return (nextLong() >>> 40) * FLOAT_NORM;
    }

    @Override
    public float nextFloat(float min, float max) {
        return min + (max - min) * nextFloat();
    }

    @Override
    public int nextInt() {
        return (int) (nextLong() >>> 32);
    }

    @Override
    public void nextBytes(byte[] bytes) {
        // awful code (adapted from java.util.Random)
        for (int i = 0, len = bytes.length; i < len; /**/) {
            for (long rnd = nextLong(), n = Math.min(len - i, Long.SIZE / Byte.SIZE); n-- > 0; rnd >>= Byte.SIZE) {
                bytes[i++] = (byte) rnd;
            }
        }
    }

    @Override
    public void nextLongs(long[] longs) {
        for (int i = 0; i < longs.length; ++i) {
            longs[i] = nextLong();
        }
    }

    @Override
    public void nextDoubles(double[] doubles) {
        for (int i = 0; i < doubles.length; ++i) {
            doubles[i] = nextDouble();
        }
    }

    @Override
    public boolean nextBoolean() {
        return nextLong() < 0L;
    }

    @Override
    public long nextLong(long n) {
        if (n <= 0L) {
            throw new IllegalArgumentException("n must be positive");
        }
        final long nMinus1 = n - 1L;
        long x = nextLong();
        if ((n & nMinus1) == 0L) {
            // power of two shortcut
            return x & nMinus1;
        }
        // rejection-based algorithm to get uniform longs
        for (long y = x >>> 1; y + nMinus1 - (x = y % n) < 0L; y = nextLong() >>> 1) {
            ;
        }
        return x;
    }

    @Override
    public int nextInt(int n) {
        return (int) nextLong(n);
    }

    @Override
    public int nextInt(int min, int max) {
        return (int) nextLong(min, max);
    }

    @Override
    public long nextLong(long min, long max) {
        return min + nextLong((max - min) + 1L);
    }

    @Override
    public int next(int bits) {
        return (int) (nextLong() >>> (64 - bits));
    }

    @Override
    public String getAlgorithm() {
        return getClass().getSimpleName();
    }

    @Override
    public long[] getSeed() {
        return initialSeed;
    }

    public IntStream ints() {
        return intStream(
                new PseudoRandomIntSpliterator(this, 0L, Long.MAX_VALUE, Integer.MIN_VALUE, Integer.MAX_VALUE));
    }

    public IntStream ints(long streamSize) {
        checkStreamSize(streamSize);
        return intStream(new PseudoRandomIntSpliterator(this, 0L, streamSize, Integer.MIN_VALUE, Integer.MAX_VALUE));
    }

    public IntStream ints(int min, int max) {
        checkRange(min, max);
        return intStream(new PseudoRandomIntSpliterator(this, 0L, Long.MAX_VALUE, min, max));
    }

    public IntStream ints(long streamSize, int min, int max) {
        checkStreamSize(streamSize);
        checkRange(min, max);
        return intStream(new PseudoRandomIntSpliterator(this, 0L, streamSize, min, max));
    }

    protected void saveSeed(long[] seed) {
        initialSeed = Arrays.copyOf(seed, seed.length);
    }

    protected void saveSeed(long seed) {
        initialSeed = new long[] { seed };
    }

    private static IntStream intStream(Spliterator.OfInt spliterator) {
        return StreamSupport.intStream(spliterator, false);
    }

    private static LongStream longStream(Spliterator.OfLong spliterator) {
        return StreamSupport.longStream(spliterator, false);
    }

    private static DoubleStream doubleStream(Spliterator.OfDouble spliterator) {
        return StreamSupport.doubleStream(spliterator, false);
    }

    private static void checkStreamSize(long streamSize) {
        if (streamSize < 0L) {
            throw new IllegalArgumentException("stream size must be non-negative");
        }
    }

    private static void checkRange(double min, double max) {
        if (!(min <= max && (max - min) < Double.POSITIVE_INFINITY)) {
            throw new IllegalArgumentException(BAD_RANGE);
        }
    }

    private static void checkRange(int min, int max) {
        if (min > max) {
            throw new IllegalArgumentException(BAD_RANGE);
        }
    }

    private static void checkRange(long min, long max) {
        if (min > max) {
            throw new IllegalArgumentException(BAD_RANGE);
        }
    }
}
