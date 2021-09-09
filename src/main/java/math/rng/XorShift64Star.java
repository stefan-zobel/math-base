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
 * 64-bit Xorshift* pseudo random generator by
 * <a href=https://www.jstatsoft.org/article/view/v008i14>George Marsaglia
 * (2003)</a>.
 * <p>
 * A faster and statistically better pseudo RNG than {@link java.util.Random}
 * This generator has a period of 2<sup>64</sup>&nbsp;&minus;&nbsp;1.
 */
public class XorShift64Star extends Marsaglia64 {

    private static final XorShift64Star defaultRng = new XorShift64Star();

    public XorShift64Star() {
    }

    public XorShift64Star(long seed) {
        super(seed);
    }

    public XorShift64Star(long[] seed) {
        super(seed);
    }

    public final long nextLong() {
        long x = seed;
        x ^= (x >>> 12);
        x ^= (x << 25);
        x ^= (x >>> 27);
        seed = x;
        return x * 0x2545F4914F6CDD1DL;
    }

    public static PseudoRandom getDefault() {
        return defaultRng;
    }
}
