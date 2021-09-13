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
 * This generator has a guaranteed period of at least 2<sup>64</sup> and an
 * average period of 2<sup>255</sup>.
 */
public class Sfc64 extends AbstractRng64 implements SplittablePseudoRandom {

    private static final Sfc64 defaultRng = new Sfc64();

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

    protected Sfc64(long a, long b, long c, long counter) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.counter = counter;
        saveSeed(new long[] { a, b, c, counter });
    }

    @Override
    public final long nextLong() {
        long xa = a;
        long xb = b;
        long xc = c;
        long rnd = xa + xb + counter++;
        a = xb ^ (xb >>> 11);
        b = xc + (xc << 3);
        c = ((xc << 24) | (xc >>> 40)) + rnd;
        return rnd;
    }

    @Override
    public Sfc64 split() {
        long l = 0L;
        if ((l = nextLong()) == 0L) {
            unused = (byte) l;
        }
        long[] mix = Seed.get4Constants();
        SpookyMix.mix(new long[] { a, b, c, counter }, mix);
        return new Sfc64(mix[0], mix[1], mix[2], mix[3]);
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
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
