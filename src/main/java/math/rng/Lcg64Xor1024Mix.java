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
 * The {@code L64X1024MixRandom} algorithm from JDK 17.
 */
public class Lcg64Xor1024Mix extends AbstractRng64 implements SplittablePseudoRandom {

    private static final Lcg64Xor1024Mix defaultRng = new Lcg64Xor1024Mix();

    /*
     * Comment taken from the Java 17 implementation: Multiplier used in the LCG
     * portion of the algorithm. Chosen based on research by Sebastiano Vigna
     * and Guy Steele (2019). The spectral scores for dimensions 2 through 8 for
     * the multiplier 0xd1342543de82ef95 are [0.958602, 0.937479, 0.870757,
     * 0.822326, 0.820405, 0.813065, 0.760215].
     */
    private static final long M = 0xd1342543de82ef95L;

    /*
     * Comment taken from the Java 17 implementation: The parameter that is used
     * as an additive constant for the LCG. Must be odd.
     */
    private final long a;

    /*
     * Comment taken from the Java 17 implementation: The per-instance state: s
     * for the LCG; the array seed for the XBG; pos is the rotating pointer into
     * the array seed. At least one of the 16 elements of the array seed must be
     * nonzero.
     */
    private long s;
    private int pos = 15;
    private final long[] seed;

    public Lcg64Xor1024Mix() {
        this(new XorShift64Star());
    }

    public Lcg64Xor1024Mix(long seed) {
        this(new XorShift64Star(seed));
    }

    public Lcg64Xor1024Mix(long[] seed) {
        this(new XorShift64Star(seed));
    }

    private Lcg64Xor1024Mix(XorShift64Star seeder) {
        seed = new long[16];
        a = seeder.nextLong();
        s = seeder.nextLong();
        seeder.nextLongs(seed);
        saveSeed();
    }

    protected Lcg64Xor1024Mix(long a, long s, long[] seed) {
        this.a = a;
        this.s = s;
        this.seed = seed;
        saveSeed();
    }

    private void saveSeed() {
        long[] state = new long[seed.length + 2];
        System.arraycopy(seed, 0, state, 2, seed.length);
        state[0] = a;
        state[1] = s;
        saveSeed(state);
    }

    @Override
    public Lcg64Xor1024Mix split() {
        long l = 0L;
        if ((l = nextLong()) == 0L) {
            unused = (byte) l;
        }
        long[] mix = Seed.get2Constants();
        long[] mixSeed = Seed.get16Constants();
        SpookyMix.mix(new long[] { a, s }, mix);
        SpookyMix.mix(seed, mixSeed);
        return new Lcg64Xor1024Mix(mix[0], mix[1], mixSeed);
    }

    // this is actually faster than the JDK 17 implementation
    @Override
    public final long nextLong() {
        // xoroshiro1024: part 1
        long[] x = seed;
        int p = pos;
        long s15 = x[p];
        long s0 = x[pos = (p + 1) & 15];
        // compute result
        long rnd = BitMix.leaMix64(s + s0);
        // update LCG sub-generator
        s = M * s + a;
        // xoroshiro1024: part 2
        s15 ^= s0;
        x[p] = ((s0 << 25) | (s0 >>> 39)) ^ s15 ^ (s15 << 27);
        x[pos] = ((s15 << 36) | (s15 >>> 28));
        return rnd;
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
    }
}
