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
 * Sebastiano Vigna's splitmix64() generator.
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
        return mix64(state += GOLDEN);
    }

    private static long mix64(long seed) {
        seed = (seed ^ (seed >>> 30)) * 0xbf58476d1ce4e5b9L;
        seed = (seed ^ (seed >>> 27)) * 0x94d049bb133111ebL;
        return seed ^ (seed >>> 31);
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
        return mix64(seed + GOLDEN);
    }

    private SplitMix64Seed() {
    }
}
