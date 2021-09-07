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
 * 256-bit {@code xoshiro256**} pseudo random generator suggested by
 * <a href=https://arxiv.org/pdf/1805.01407.pdf>D. Blackman and S. Vigna
 * (2019)</a>
 * <p>
 * This generator has a period of 2<sup>256</sup>&nbsp;&minus;&nbsp;1.
 */
public class XorShiftRot256StarStar extends AbstractRng64 {

    private long x0;
    private long x1;
    private long x2;
    private long x3;

    public XorShiftRot256StarStar() {
        x0 = SplitMix64Seed.seed();
        x1 = SplitMix64Seed.seed();
        x2 = SplitMix64Seed.seed();
        x3 = SplitMix64Seed.seed();
        escape();
    }

    public XorShiftRot256StarStar(long seed) {
        this(new XorShift64Star(seed));
    }

    public XorShiftRot256StarStar(long[] seed) {
        this(new XorShift64Star(seed));
    }

    private XorShiftRot256StarStar(XorShift64Star seeder) {
        x0 = seeder.nextLong();
        x1 = seeder.nextLong();
        x2 = seeder.nextLong();
        x3 = seeder.nextLong();
        escape();
    }

    @Override
    public final long nextLong() {
        long s1 = x1;
        long t = s1 << 17;
        long x = s1 + (s1 << 2);
        long rnd = ((x << 7) | (x >>> 57));
        rnd += rnd << 3;

        long s2 = (x2 ^= x0);
        long s3 = (x3 ^= s1);
        x1 ^= s2;
        x0 ^= s3;

        x2 ^= t;
        x3 = ((s3 << 45) | (s3 >>> 19));

        return rnd;
    }

    private void escape() {
        saveSeed(new long[] { x0, x1, x2, x3 });
        long l = 0L;
        for (int i = 0; i < 20; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
