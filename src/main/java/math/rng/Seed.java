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

    private static final long[] SEED_CONSTS = new long[] { 0x93bf7da7c7e7cL, 0xfff049c37a98faddL, 0x1939319101e278L,
            0xffdf90f07376f98fL, 0x4f92dedd83aec3L, 0xff832e4bda1cf731L, 0x144bc69ac76f4fL, 0xffe840ec239b3f07L,
            0x1dfaecebb7a3c6L, 0xffcdc5812b9c1c35L, 0x760766c165e061L, 0xff60c2e1247a67ffL, 0x11e260f1b32067L,
            0xffea5ec251a22640L, 0x1ade60fd4a549fL, 0xffe0463eb27636c4L, 0x2145646b51d0a3L, 0xffcdf8b346c2121aL,
            0x3fb6ce0d75f2b2L, 0xfff5aad25c9b98ceL, 0xcdd03798e7898L, 0xff5ae3ae06a86a18L, 0x1520c5ea46ad7aL,
            0xffe8a88888bba4b3L, 0x186d36c71196a6L, 0xffe571e578ab7baeL, 0x1b99773f1f2057L, 0xffe25977177f5d9dL,
            0x371d9a98f6f744L, 0xfff999fc9821f955L, 0xab5c160c2eb0cL, 0xfff473781028c89bL, 0xe05b0ecea14e2L,
            0xff5baad0dc7d4c29L, 0x120208596985f4L, 0xffebad3176a48a4fL, 0x1696076247c437L, 0xffe8abc6ffc2bc96L,
            0x1af6bf9ba4ac6bL, 0xffe45311488750f7L, 0x1d158f19620408L, 0xfffd05987bbf894eL, 0x21de97e1378030L,
            0xffdc6992c6d96453L, 0x29ef0df4abfbc1L, 0xffc97bd98febd85aL, 0x491f62192df0d5L, 0xffb0be9af944579cL };

    private Seed() {
    }
}
