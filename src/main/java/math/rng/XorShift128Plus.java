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
 * 128-bit {@code xorshift128+} pseudo random generator suggested by
 * <a href=https://vigna.di.unimi.it/ftp/papers/xorshiftplus.pdf>Sebastiano
 * Vigna (2017)</a>. It is about 40% faster than {@link XorShift64Star}
 * despite its larger state size.
 * <p>
 * This generator has a period of 2<sup>128</sup>&nbsp;&minus;&nbsp;1.
 */
public class XorShift128Plus extends AbstractRng64 {

    private static final XorShift128Plus defaultRng = new XorShift128Plus();

    private long x0;
    private long x1;

    public XorShift128Plus() {
        x0 = SplitMix64Seed.seed();
        x1 = SplitMix64Seed.seed();
        escape();
    }

    public XorShift128Plus(long seed) {
        this(new XorShift64Star(seed));
    }

    public XorShift128Plus(long[] seed) {
        this(new XorShift64Star(seed));
    }

    private XorShift128Plus(XorShift64Star seeder) {
        long[] x = new long[2];
        seeder.nextLongs(x);
        x0 = x[0];
        x1 = x[1];
        escape();
    }

    @Override
    public final long nextLong() {
        long s0 = x1;
        long s1 = x0;
        long s = x0 + x1;
        s1 ^= (s1 << 23);
        x1 = s1 ^ s0 ^ (s1 >>> 18) ^ (s0 >>> 5);
        x0 = s0;
        return s;
    }

    public static PseudoRandom getDefault() {
        return defaultRng;
    }

    private void escape() {
        saveSeed(new long[] { x0, x1 });
        long l = 0L;
        for (int i = 0; i < 20; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            // this cannot happen
            throw new AssertionError("0L");
        }
    }
}
