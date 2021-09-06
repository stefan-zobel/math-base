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
 * 1024-bit Xorshift pseudo random generator suggested by
 * <a href=https://www.jstatsoft.org/article/view/v008i14>George Marsaglia
 * (2003)</a> as studied by
 * <a href=https://arxiv.org/pdf/1402.6246.pdf>Sebastiano Vigna (2016)</a>. It
 * has better statistical properties than {@link XorShift64Star} and other
 * 64-bit xorshift generators and is about 25% faster despite its larger state
 * size.
 * <p>
 * This generator has a period of 2<sup>1024</sup>&nbsp;&minus;&nbsp;1.
 */
public class XorShift1024Star extends AbstractRng64 {

    private static final XorShift1024Star defaultRng = new XorShift1024Star();

    private int pos = 0;
    private final long[] seed = new long[16];

    public XorShift1024Star() {
        XorShift64Star seeder = new XorShift64Star();
        seeder.nextLongs(this.seed);
        escape();
    }

    public XorShift1024Star(long seed) {
        XorShift64Star seeder = new XorShift64Star(seed);
        seeder.nextLongs(this.seed);
        escape();
    }

    public XorShift1024Star(long[] seed) {
        XorShift64Star seeder = new XorShift64Star(seed);
        seeder.nextLongs(this.seed);
        escape();
    }

    @Override
    public final long nextLong() {
        long[] x = seed;
        long s0 = x[pos];
        long s1 = x[pos = (pos + 1) & 15];
        s1 ^= (s1 << 31);
        long s = s1 ^ s0 ^ (s1 >>> 11) ^ (s0 >>> 30);
        x[pos] = s;
        return s * 0x106689D45497FDB5L;
    }

    public static PseudoRandom getDefault() {
        return defaultRng;
    }

    /*
     * Escape from "zeroland"
     */
    private void escape() {
        saveSeed(seed);
        long l = 0L;
        for (int i = 0; i < 200; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
