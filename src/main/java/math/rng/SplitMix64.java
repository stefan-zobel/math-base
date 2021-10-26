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
 * 64-bit {@code SplitMix}-type pseudo random generator suggested by
 * <a href="http://gee.cs.oswego.edu/dl/papers/oopsla14.pdf">Guy Steele, Doug Lea,
 * Christine Flood (2014)</a>. This implementation is basically the same
 * algorithm as used in {@link java.util.SplittableRandom} but it employs
 * different mixing functions to (hopefully) protect better against the possible
 * problems with weak gamma values reported in
 * <a href="https://labs.oracle.com/pls/apex/f?p=LABS:0::APPLICATION_PROCESS%3DGETDOC_INLINE:::DOC_ID:1050">Steele
 * (2017)</a>. Performance should be about the same as for
 * {@code SplittableRandom}.
 * <p>
 * This generator has a period of 2<sup>64</sup>.
 * <p>
 * This generator is {@code 1}-dimensionally equidistributed.
 */
public class SplitMix64 extends AbstractRng64 implements SplittablePseudoRandom {

    /*
     * the golden ratio scaled to 64 bits
     */
    private static final long GOLDEN = 0x9e3779b97f4a7c15L;

    private static final SplitMix64 defaultRng = new SplitMix64();

    private long state;

    /*
     * Weyl generator step value
     */
    private final long gamma;

    public SplitMix64() {
        state = SplitMix64Seed.seed();
        gamma = mixGamma(state + GOLDEN);
    }

    public SplitMix64(long seed) {
        this(SplitMix64Seed.seed(seed), GOLDEN);
    }

    public SplitMix64(long[] seed) {
        this(SplitMix64Seed.seed(seed), GOLDEN);
    }

    protected SplitMix64(long seed, long gamma) {
        this.state = seed;
        this.gamma = gamma;
        saveSeed(new long[] { seed, gamma });
    }

    @Override
    public final long nextLong() {
        return BitMix.xnasam(state += gamma);
    }

    @Override
    public final int nextInt() {
        return BitMix.staffordMix04(state += gamma);
    }

    @Override
    public SplitMix64 split() {
        return new SplitMix64(nextLong(), mixGamma(state += gamma));
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
    }

    /*
     * compute the next gamma value from v
     */
    private static long mixGamma(long v) {
        // force v to be odd
        v = BitMix.rrxmrrxmsx(v) | 1L;
        // try to support enough 01 and 10 transitions
        int n = Long.bitCount(v ^ (v >>> 1));
        return (n < 24) ? v ^ 0xaaaaaaaaaaaaaaaaL : v;
    }
}
