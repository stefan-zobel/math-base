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

    private static final long[] SEED_CONSTS = new long[] { 0xfb63c85670523L, 0xffb06d21227c513dL, 0x7cd1b425e308cfL,
            0xffebb439653890b1L, 0x17bf13dc64c0f9L, 0xffcdc5812b9c1c35L, 0x760766c165e061L, 0xff60c2e1247a67ffL,
            0x1ade60fd4a549fL, 0xffe8a88888bba4b3L, 0x1b99773f1f2057L, 0xfff999fc9821f955L, 0xb8c87efd73765L,
            0xff5baad0dc7d4c29L, 0x1452ce895b75b1L, 0xffe969f89db83bc9L, 0x1af6bf9ba4ac6bL, 0xffe45311488750f7L,
            0x23966d39269badL, 0xffd610f20b54043fL, 0x491f62192df0d5L, 0xff80dcb71e4dd83bL, 0xe754446507ed7L,
            0xff5891e9cf4d3f97L, 0x1d7bd8283bf89L, 0xffa2dc029f151763L, 0x6a8b2b6fb43e1fL, 0xff7b0ee546175dd9L,
            0x963fd0afe8e43bL, 0xff546e4ec30df0a9L, 0x13abc5414dadd3L, 0xffeaadc5a9f14463L, 0x1a279681925db3L,
            0xffe50e5f20621bb5L, 0x1ed2e4d7050c75L, 0xfff5d25d0d9a3887L, 0xa85d5a975eb3bL, 0xfffe92ee3892e807L,
            0x1c2061e9eade7L, 0xffe8f0fb46d751cdL, 0x18016e82575aefL, 0xffe6bda6502d6f19L, 0x1b1ee8a24225e9L,
            0xffe4437f8aa9baabL, 0x2dd92bafa1cf3L, 0xffe220c306f2da2dL, 0x1f61423134717bL, 0xffdcc3fa2171603bL };

    private Seed() {
    }
}
