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
 * The 256-bit generator {@code Stc64} is Tyge Løvset's improved variation of
 * {@link Sfc64}. See
 * "https://github.com/tylov/STC/blob/master/include/stc/crandom.h".
 * <p>
 * This generator has a guaranteed period of at least 2<sup>64</sup> and an
 * average period of 2<sup>255</sup>.
 */
public class Stc64 extends AbstractRng64 implements SplittablePseudoRandom {

    private static final Stc64 defaultRng = new Stc64();

    private long s0;
    private long s1;
    private long s2;
    private long s3;
    private long inc = 0L;
    private final long seq;

    public Stc64() {
        seq = init(SplitMix64Seed.seed());
        escape();
    }

    public Stc64(long seed) {
        seq = init(SplitMix64Seed.seed(seed));
        escape();
    }

    public Stc64(long[] seed) {
        seq = init(SplitMix64Seed.seed(seed));
        escape();
    }

    public Stc64(long seed, long seq) {
        init(SplitMix64Seed.seed(seed));
        this.seq = seq | 1L;
        escape();
    }

    private long init(long seed) {
        s0 = seed;
        s1 = seed + 0x26aa069ea2fb1a4dL;
        s2 = seed + 0x70c72c95cd592d04L;
        s3 = seed + 0x504f333d3aa0b359L;
        // seq must be odd
        return (((seed + 0x3504f333d3aa0b37L) << 1) | 1L);
    }

    protected Stc64(long s0, long s1, long s2, long s3, long seq) {
        this.s0 = s0;
        this.s1 = s1;
        this.s2 = s2;
        this.s3 = s3;
        this.seq = seq;
        saveSeed(new long[] { s0, s1, s2, s3, seq });
    }

    @Override
    public final long nextLong() {
        long xb = s1;
        long xc = s2;
        long rnd = (s0 ^ (s3 += seq)) + xb;

        s0 = xb ^ (xb >>> 11);
        s1 = xc + (xc << 3);
        s2 = ((xc << 24) | (xc >>> 40)) + rnd;

        return rnd;
    }

    @Override
    public Stc64 split() {
        long l = 0L;
        if ((l = nextLong()) == 0L) {
            unused = (byte) l;
        }
        long[] mix = Seed.get4Constants();
        SpookyMix.mix(new long[] { s0, s1, s2, s3 }, mix);
        return new Stc64(mix[0], mix[1], mix[2], mix[3], seq + (inc += 2L));
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
    }

    private void escape() {
        saveSeed(new long[] { s0, s1, s2, s3, seq });
        long l = 0L;
        for (int i = 0; i < 12; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
