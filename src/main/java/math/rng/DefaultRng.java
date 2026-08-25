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

import java.util.Objects;

/**
 * Factory methods for the default {@link PseudoRandom} algorithm which is
 * currently a {@link XorShiftRot256StarStar}: period
 * 2<sup>256</sup>&nbsp;&minus;&nbsp;1, {@code 4}-dimensionally equidistributed,
 * and splittable.
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
        return new XorShiftRot256StarStar();
    }

    /**
     * Returns a new generator for the given seed. The same seed always yields
     * the same sequence.
     *
     * @param seed the seed to start from
     * @return a new {@link PseudoRandom}
     */
    public static PseudoRandom newPseudoRandom(long seed) {
        return new XorShiftRot256StarStar(seed);
    }

    /**
     * Returns a new generator for the given seed values. They are hashed down
     * to a single {@code long} through {@link SplitMix64Seed#seed(long[])}, so
     * a longer array does not carry more of the seed into the generator. The
     * same array always yields the same sequence.
     * <p>
     * That mapping accepts {@code null} and an empty array, sending both onto
     * one fixed seed; this factory rejects them instead, because having
     * nothing to seed with is far more likely a mistake than a request for
     * that one well-known stream.
     *
     * @param seed the seed values to start from, neither {@code null} nor
     *            empty
     * @return a new {@link PseudoRandom}
     * @throws NullPointerException if {@code seed} is {@code null}
     * @throws IllegalArgumentException if {@code seed} is empty
     */
    public static PseudoRandom newPseudoRandom(long[] seed) {
        if (Objects.requireNonNull(seed, "seed").length == 0) {
            throw new IllegalArgumentException("seed must not be empty");
        }
        return new XorShiftRot256StarStar(seed);
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
        final PseudoRandom[] multiplePrng = new PseudoRandom[count];
        for (int i = 0; i < multiplePrng.length; ++i) {
            // the no-argument constructor takes four independent values from
            // SplitMix64Seed, so every generator starts from 256 bits of its
            // own. All of them sit on the same cycle of length 2^256 - 1, but
            // the chance that one of them starts within reach of another is
            // not a number anybody has to plan around
            multiplePrng[i] = new XorShiftRot256StarStar();
        }
        return multiplePrng;
    }

    /**
     * Returns a new generator independent of the given one. It is derived from
     * the state {@code prng} has reached rather than from the seed it started
     * at, so repeated calls on the same generator do not repeat themselves --
     * at the price of <b>advancing {@code prng} by one value</b>.
     * <p>
     * A generator that can {@link SplittablePseudoRandom#split() split} is
     * split, and the result is then of the same algorithm as {@code prng};
     * one that cannot is re-created from its current state; and for anything
     * that is neither, a generator of the default algorithm is seeded from a
     * single draw.
     *
     * @param prng the generator to derive from, advanced by one value
     * @return a new {@link PseudoRandom} independent of {@code prng}
     * @throws NullPointerException if {@code prng} is {@code null}
     */
    public static PseudoRandom newIndepPseudoRandom(PseudoRandom prng) {
        Objects.requireNonNull(prng, "prng");
        PseudoRandom indep = PseudoRandomSpliterator.detach(prng);
        if (indep != null) {
            return indep;
        }
        return new XorShiftRot256StarStar(prng.nextLong());
    }

    private DefaultRng() {
        throw new AssertionError();
    }
}
