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

/**
 * Factory methods for the default {@link PseudoRandom} algorithm which is
 * currently a {@link MersenneTwister64}.
 * <p>
 * The generators returned are not thread-safe; see {@link PseudoRandom} for
 * what that rules out. To hand one to each of several threads, take them from
 * {@link #newIndepPseudoRandoms(int)} rather than sharing a single one.
 */
public final class DefaultRng {

    /**
     * Returns a new generator seeded unpredictably, and therefore not
     * reproducible.
     *
     * @return a new {@link PseudoRandom}
     */
    public static PseudoRandom newPseudoRandom() {
        return new MersenneTwister64();
    }

    /**
     * Returns a new generator for the given seed. The same seed always yields
     * the same sequence.
     *
     * @param seed the seed to start from
     * @return a new {@link PseudoRandom}
     */
    public static PseudoRandom newPseudoRandom(long seed) {
        return new MersenneTwister64(seed);
    }

    /**
     * Returns a new generator for the given seed values, all of which are used
     * as seed material. The same array always yields the same sequence.
     *
     * @param seed the seed values to start from, neither {@code null} nor
     *            empty
     * @return a new {@link PseudoRandom}
     * @throws NullPointerException if {@code seed} is {@code null}
     * @throws IllegalArgumentException if {@code seed} is empty
     */
    public static PseudoRandom newPseudoRandom(long[] seed) {
        return new MersenneTwister64(seed);
    }

    /**
     * Returns {@code count} generators seeded so that their sequences do not
     * overlap in any practically reachable way. All of them are seeded
     * unpredictably and are therefore not reproducible.
     *
     * @param count the number of generators to create, positive
     * @return an array of {@code count} independently seeded generators
     * @throws IllegalArgumentException if {@code count <= 0}
     */
    public static PseudoRandom[] newIndepPseudoRandoms(int count) {
        if (count <= 0) {
            throw new IllegalArgumentException("count <= 0 : " + count);
        }
        if (count == 1) {
            return new PseudoRandom[] { newPseudoRandom() };
        }
        final int NN = 312;
        final long[] seed = SplitMix64Seed.seed(NN);
        final PseudoRandom[] multiplePrng = new PseudoRandom[count];
        for (int i = 0; i < multiplePrng.length; ++i) {
            final PseudoRandom prng = newPseudoRandom(seed);
            // now change the seed
            reseed(NN, seed, prng);
            multiplePrng[i] = prng;
        }
        return multiplePrng;
    }

    /**
     * Returns a new generator seeded from the given one. Note that this
     * <b>advances {@code prng}</b>: the seed is 312 {@code long} values drawn
     * from it, so the caller's own sequence moves on by that much.
     *
     * @param prng the generator to draw the seed from, advanced by 312 values
     * @return a new {@link PseudoRandom} independent of {@code prng}
     */
    public static PseudoRandom newIndepPseudoRandom(PseudoRandom prng) {
        final int NN = 312;
        final long[] seed = new long[NN];
        reseed(NN, seed, prng);
        return newPseudoRandom(seed);
    }

    private static void reseed(int len, long[] seed, PseudoRandom prng) {
        prng.nextLongs(seed);
        int j = 0;
        while (j < seed.length && seed[j] == 0L) {
            ++j;
        }
        final long nucleus = (j < seed.length) ? seed[j] : -1L;
        final long[] half_seed = new long[len / 2];
        new XorShift1024Star(nucleus ^ SplitMix64Seed.seed()).nextLongs(half_seed);
        System.arraycopy(half_seed, 0, seed, 0, half_seed.length);
    }

    private DefaultRng() {
        throw new AssertionError();
    }
}
