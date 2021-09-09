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
 * 256-bit {@code xoshiro256++} pseudo random generator suggested by
 * <a href=https://arxiv.org/pdf/1805.01407.pdf>David Blackman and Sebastiano
 * Vigna (2019)</a>. This implementation is not significantly faster than
 * {@link XorShiftRot256StarStar} which has slightly better statistical
 * properties. It is about 40% faster than {@link XorShift64Star} despite having
 * a 4 times larger state space.
 * <p>
 * This generator has a period of 2<sup>256</sup>&nbsp;&minus;&nbsp;1.
 */
public class XorShiftRot256PlusPlus extends Xoshiro256 {

    private static final XorShiftRot256PlusPlus defaultRng = new XorShiftRot256PlusPlus();

    public XorShiftRot256PlusPlus() {
    }

    public XorShiftRot256PlusPlus(long seed) {
        super(seed);
    }

    public XorShiftRot256PlusPlus(long[] seed) {
        super(seed);
    }

    protected XorShiftRot256PlusPlus(long x0, long x1, long x2, long x3) {
        super(x0, x1, x2, x3);
    }

    @Override
    public final long nextLong() {
        long s0 = x0;
        long s1 = x1;
        long s3 = x3;
        long x = s0 + s3;
        long t = s1 << 17;
        long rnd = ((x << 23) | (x >>> 41)) + s0;

        long s2 = (x2 ^= s0);
        s3 = (x3 ^= s1);
        x1 ^= s2;
        x0 ^= s3;

        x2 ^= t;
        x3 = ((s3 << 45) | (s3 >>> 19));

        return rnd;
    }

    @Override
    public XorShiftRot256PlusPlus split() {
        long l = 0L;
        if ((l = nextLong()) == 0L) {
            unused = (byte) l;
        }
        long[] mix = Seed.get4Constants();
        SpookyMix.mix(new long[] { x0, x1, x2, x3 }, mix);
        return new XorShiftRot256PlusPlus(mix[0], mix[1], mix[2], mix[3]);
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
    }
}
