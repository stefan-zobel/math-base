/*
 * Copyright 2013 Stefan Zobel
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

public final class Seed {

    // from java.util.Random
    private static long seedUniquifier = 0x1ed8b55fac9decL;

    private static long lastSeed = pseudoRandomSeed();

    // from java.util.Random
    private static long nextSeedUniquifier() {
        // Pierre L'Ecuyer: "Tables of Linear Congruential Generators
        // of Different Sizes and Good Lattice Structure"
        seedUniquifier *= 0x106689d45497fdb5L;
        return seedUniquifier;
    }

    /*
     * Returns a reasonably good (pseudo) random seed
     */
    private static long pseudoRandomSeed() {
        long seed = nextSeedUniquifier() ^ System.nanoTime();
        // apply Austin Appleby's fmix64() hash
        return fmix64(seed);
    }

    /*
     * Austin Appleby's fmix64() hash
     */
    private static long fmix64(long seed) {
        seed ^= seed >>> 33;
        seed *= 0xff51afd7ed558ccdL;
        seed ^= seed >>> 33;
        seed *= 0xc4ceb9fe1a85ec53L;
        seed ^= seed >>> 33;

        return seed;
    }

    /**
     * Returns a reasonably good long random seed.
     * 
     * @return a long random seed.
     */
    public static synchronized long seed() {
        long seed = pseudoRandomSeed();
        while (seed == lastSeed || seed == 0L) {
            seed = pseudoRandomSeed();
        }
        lastSeed = seed;
        return seed;
    }

    /**
     * Returns an array of length {@code 1} that contains an arbitrary but fixed
     * constant that can be used as a fixed seed.
     * 
     * @return an array that contains {@code 1} fixed seed
     */
    public static long[] get1Constants() {
        return Arrays.copyOf(SEED_CONSTS, 1);
    }

    /**
     * Returns an array of length {@code 2} that contains {@code 2} arbitrary
     * but fixed constants that can be used as fixed seed values.
     * 
     * @return an array that contains {@code 2} fixed seed values
     */
    public static long[] get2Constants() {
        return Arrays.copyOf(SEED_CONSTS, 2);
    }

    /**
     * Returns an array of length {@code 4} that contains {@code 4} arbitrary
     * but fixed constants that can be used as fixed seed values.
     * 
     * @return an array that contains {@code 4} fixed seed values
     */
    public static long[] get4Constants() {
        return Arrays.copyOf(SEED_CONSTS, 4);
    }

    /**
     * Returns an array of length {@code 8} that contains {@code 8} arbitrary
     * but fixed constants that can be used as fixed seed values.
     * 
     * @return an array that contains {@code 8} fixed seed values
     */
    public static long[] get8Constants() {
        return Arrays.copyOf(SEED_CONSTS, 8);
    }

    /**
     * Returns an array of length {@code 16} that contains {@code 16} arbitrary
     * but fixed constants that can be used as fixed seed values.
     * 
     * @return an array that contains {@code 16} fixed seed values
     */
    public static long[] get16Constants() {
        return Arrays.copyOf(SEED_CONSTS, 16);
    }

    // These are fractional parts of the third root of prime numbers where each
    // second number has been multiplied by -1 and then filtered for odd numbers
    // that have a bit count of 32
    private static final long[] SEED_CONSTS = new long[] { 0x4f92dedd83aec3L, 0xffe50940645b5395L, 0x1baceeb778af09L,
            0xffd610f20b54043fL, 0x6a8b2b6fb43e1fL, 0xff4edb408816c467L, 0x1ddf25ab6283f7L, 0xffe8b4dd40032167L,
            0x9be4fd33558f7L, 0xffc86ea752251419L, 0x1bbaaccf48bfe1L, 0xff688d3f400d4a8bL, 0x9a6eb695e11f3dL,
            0xff64b00aa59c9369L, 0x9df0e7c6b2d267L, 0xffea4f554090c1d1L };

    private Seed() {
    }
}
