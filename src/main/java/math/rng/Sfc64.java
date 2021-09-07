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
 * Chris Doty-Humphrey's 256-bit "Small Fast Counting RNG" ({@code sfc64}).
 * <p>
 * This generator has a guaranteed period of at least 2<sup>64</sup> but the
 * average period can be larger.
 */
public class Sfc64 extends AbstractRng64 {

    private long a;
    private long b;
    private long c;
    private long counter;

    public Sfc64() {
        a = b = c = SplitMix64Seed.seed();
        escape();
    }

    public Sfc64(long seed) {
        a = b = c = SplitMix64Seed.seed(seed);
        escape();
    }

    public Sfc64(long[] seed) {
        a = b = c = SplitMix64Seed.seed(seed);
        escape();
    }

    @Override
    public final long nextLong() {
        long rnd = a + b + counter++;
        a = b ^ (b >>> 11);
        b = c + (c << 3);
        c = ((c << 24) | (c >>> 40)) + rnd;
        return rnd;
    }

    private void escape() {
        counter = 1L;
        saveSeed(new long[] { a, b, c, counter });
        long l = 0L;
        for (int i = 0; i < 12; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
