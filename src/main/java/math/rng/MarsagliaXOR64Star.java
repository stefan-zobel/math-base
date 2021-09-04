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
 * 64-bit Xorshift* pseudo random generator by George Marsaglia (2003).
 * <p>
 * A faster and statistically better pseudo RNG than {@link java.util.Random} or
 * {@link java.util.concurrent.ThreadLocalRandom}. This generator has a period
 * of 2<sup>64</sup>&nbsp;&minus;&nbsp;1.
 */
public class MarsagliaXOR64Star extends AbstractRng64 {

    private static final MarsagliaXOR64Star defaultRng = new MarsagliaXOR64Star();

    private long seed;

    public MarsagliaXOR64Star() {
        long seed = 0L;
        do {
            seed = Seed.seed();
        } while (seed == 0L);
        this.seed = seed;
        recover();
    }

    public MarsagliaXOR64Star(final long seed) {
        this.seed = (seed == 0L) ? -1L : seed;
        recover();
    }

    public MarsagliaXOR64Star(final long[] seed) {
        MersenneTwister64 seeder = new MersenneTwister64(seed);
        long seed_ = seeder.nextLong();
        this.seed = (seed_ == 0L) ? -1L : seed_;
        recover();
    }

    public final long nextLong() {
        long x = seed;
        x ^= (x >>> 12);
        x ^= (x << 25);
        x ^= (x >>> 27);
        seed = x;
        return x * 0x2545F4914F6CDD1DL;
    }

    public static MarsagliaXOR64Star getDefault() {
        return defaultRng;
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
