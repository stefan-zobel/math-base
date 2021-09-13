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
 * A selection of 64-bit mixing functions for longs.
 */
public final class BitMix {

    /**
     * Pelle Evensen's better mixer (see
     * "http://mostlymangling.blogspot.com/2018/07/on-mixing-functions-in-fast-splittable.html").
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long rrmxmx(long v) {
        v ^= ((v >>> 49) | (v << 15)) ^ ((v >>> 24) | (v << 40));
        v *= 0x9fb21c651e98df25L;
        v ^= v >>> 28;
        v *= 0x9fb21c651e98df25L;
        return v ^ (v >>> 28);
    }

    /**
     * Pelle Evensen's even better mixer (a bit slower than
     * {@link #rrmxmx(long)}. See
     * "https://mostlymangling.blogspot.com/2019/01/better-stronger-mixer-and-test-procedure.html".
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long rrxmrrxmsx(long v) {
        v ^= ((v >>> 25) | (v << 39)) ^ ((v >>> 50) | (v << 14));
        v *= 0xa24baed4963ee407L;
        v ^= ((v >>> 24) | (v << 40)) ^ ((v >>> 49) | (v << 15));
        v *= 0x9fb21c651e98df25L;
        return v ^ (v >>> 28);
    }

    /**
     * Pelle Evensen's {@code xNASAM} ("Not Another Strange Acronym Mixer")
     * mixer (slower than the {@code rr...} mixers). See
     * "http://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html".
     * 
     * @param v
     *            long to mix
     * @param c
     *            a constant with which {@code v} becomes xor-ed before the
     *            mixing steps
     * @return the mixed long
     */
    public static long xnasam(long v, long c) {
        v ^= c;
        v ^= ((v >>> 25) | (v << 39)) ^ ((v >>> 47) | (v << 17));
        v *= 0x9e6c63d0676a9a99L;
        v ^= (v >>> 23) ^ (v >>> 51);
        v *= 0x9e6d62d06f6a9a9bL;
        v ^= (v >>> 23) ^ (v >>> 51);
        return v;
    }

    /**
     * A variant of {@link #xnasam(long, long)} where {@code 0xff64b00aa59c9369}
     * is used as a fixed value for the constant {@code c}.
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long xnasam(long v) {
        return xnasam(v, 0xff64b00aa59c9369L);
    }

    /**
     * David Stafford's variant 13 of his 64-bit mix functions. See
     * "http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html"
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long staffordMix13(long v) {
        v = (v ^ (v >>> 30)) * 0xbf58476d1ce4e5b9L;
        v = (v ^ (v >>> 27)) * 0x94d049bb133111ebL;
        return v ^ (v >>> 31);
    }

    /**
     * Austin Appleby's fmix64() hash function used in {@code MurmurHash3}. See
     * "https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp".
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long murmurhash3(long v) {
        v ^= v >>> 33;
        v *= 0xff51afd7ed558ccdL;
        v ^= v >>> 33;
        v *= 0xc4ceb9fe1a85ec53L;
        v ^= v >>> 33;
        return v;
    }

    private BitMix() {
        throw new AssertionError();
    }
}
