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
 * Abstract base class for {@link MarsagliaXOR64} and
 * {@link MarsagliaXOR64Star}.
 */
public abstract class Marsaglia64 extends AbstractRng64 {

    protected long seed;

    public Marsaglia64() {
        long seed = 0L;
        do {
            seed = SplitMix64Seed.seed();
        } while (seed == 0L);
        this.seed = seed;
        recover();
    }

    public Marsaglia64(long seed) {
        this.seed = (seed == 0L) ? -1L : seed;
        recover();
    }

    public Marsaglia64(long[] seed) {
        MersenneTwister64 seeder = new MersenneTwister64(seed);
        long seed_ = seeder.nextLong();
        this.seed = (seed_ == 0L) ? -1L : seed_;
        recover();
    }

    /*
     * Protect against poor seeds.
     */
    private void recover() {
        long l = 0L;
        for (int i = 0; i < 10; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            // this cannot happen
            throw new AssertionError("0L");
        }
    }
}