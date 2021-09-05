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
 * <a href=https://arxiv.org/pdf/1402.6246.pdf>Sebastiano Vigna (2016)</a>.
 * <p>
 * This generator has a period of 2<sup>1024</sup>&nbsp;&minus;&nbsp;1.
 */
public class MarsagliaXOR1024Star extends AbstractRng64 {

    private static final MarsagliaXOR1024Star defaultRng = new MarsagliaXOR1024Star();

    private int pos = 0;
    private final long[] seed = new long[16];

    public MarsagliaXOR1024Star() {
        MarsagliaXOR64Star seeder = new MarsagliaXOR64Star();
        seeder.nextLongs(this.seed);
        recover();
    }

    public MarsagliaXOR1024Star(long seed) {
        MarsagliaXOR64Star seeder = new MarsagliaXOR64Star(seed);
        seeder.nextLongs(this.seed);
        recover();
    }

    public MarsagliaXOR1024Star(long[] seed) {
        MersenneTwister64 seeder = new MersenneTwister64(seed);
        seeder.nextLongs(this.seed);
        recover();
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
     * Protect against poor seeds.
     */
    private void recover() {
        long l = 0L;
        for (int i = 0; i < 10; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            pos = 0;
        }
    }
}
