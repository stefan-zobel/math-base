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

/**
 * The splitmix64() generator (here only used for seed generation) which is
 * based on <a href=http://gee.cs.oswego.edu/dl/papers/oopsla14.pdf>Guy Steele,
 * Doug Lea, Christine Flood (2014)</a> and has been included in Java 8 as
 * {@link java.util.SplittableRandom} (this implementation uses a mixing
 * function which is different from that used in {@code SplittableRandom}.
 */
public final class SplitMix64Seed {

    // the golden ratio scaled to 64 bits
    private static final long GOLDEN = 0x9e3779b97f4a7c15L;

    private static long state = Seed.seed();

    /**
     * Returns a statistically good long random seed.
     * 
     * @return a random long seed
     */
    public static synchronized long seed() {
        return BitMix.rrxmrrxmsx(state += GOLDEN);
    }

    /**
     * Returns a statistically good long[] random seed array.
     * 
     * @param numLongs
     *            length ({@code >= 0}) of the seed array
     * @return a random long[] seed array
     */
    public static synchronized long[] seed(int numLongs) {
        long[] seed = new long[numLongs];
        for (int i = 0; i < seed.length; ++i) {
            seed[i] = seed();
        }
        return seed;
    }

    /**
     * Computes a deterministic seed value from a given value.
     * 
     * @param seed
     *            the seed to start with
     * @return a deterministically computed seed value
     */
    public static long seed(long seed) {
        seed = (seed == 0L) ? -1L : seed;
        return BitMix.rrxmrrxmsx(seed + GOLDEN);
    }

    /**
     * Computes a deterministic seed value from a given sequence of seed values.
     * 
     * @param seed
     *            the sequence of seed values to start with
     * @return a deterministically computed seed value
     */
    public static long seed(long[] seed) {
        if (seed == null || seed.length == 0) {
            return seed(0L);
        }
        long s = seed(seed[0]);
        for (int i = 1; i < seed.length; ++i) {
            s += ((i & 1) != 0) ? flipMix64(seed[i] + GOLDEN) : seed(seed[i]);
        }
        return (seed.length > 1) ? seed(s) : s;
    }

    private static long flipMix64(long s) {
        s = (s ^ (s >>> 33)) * 0xff51afd7ed558ccdL;
        s = (s ^ (s >>> 33)) * 0xc4ceb9fe1a85ec53L;
        // force to be odd
        s = (s ^ (s >>> 33)) | 1L;
        // try to support enough transitions
        int n = Long.bitCount(s ^ (s >>> 1));
        return (n < 24) ? s ^ 0xaaaaaaaaaaaaaaaaL : s;
    }

    private SplitMix64Seed() {
    }
}
