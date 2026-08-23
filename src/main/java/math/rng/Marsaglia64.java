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
 * Abstract base class for {@link XorShift64Star} and other 64-bit xorshift
 * generators.
 */
abstract class Marsaglia64 extends AbstractRng64 {

    protected long seed;

    Marsaglia64() {
        long seed = 0L;
        do {
            seed = SplitMix64Seed.seed();
        } while (seed == 0L);
        this.seed = seed;
        escape();
    }

    Marsaglia64(long seed) {
        this.seed = nonZero(SplitMix64Seed.seed(seed));
        escape();
    }

    Marsaglia64(long[] seed) {
        this.seed = nonZero(SplitMix64Seed.seed(seed));
        escape();
    }

    /*
     * The all-zero state is the fixed point of every xorshift: once there the
     * generator can never leave, nextLong() returns 0 forever and
     * nextGaussian() never terminates. It is reachable, because the mixer in
     * SplitMix64Seed maps 0 to 0 and therefore turns the seed -GOLDEN into a
     * zero state. XorShift1024Star, XorShift128Plus, Xoshiro256 and
     * Lcg64Xor1024Mix all seed themselves through a XorShift64Star, so they
     * would inherit it. Every other seed is left exactly as it was.
     */
    private static long nonZero(long state) {
        return (state == 0L) ? -1L : state;
    }

    /*
     * Attempt to protect against poor seeds
     */
    private void escape() {
        saveSeed(seed);
        long l = 0L;
        for (int i = 0; i < 10; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
