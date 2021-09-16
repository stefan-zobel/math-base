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
 * A selection of 64-bit mixing functions for longs and 32-bit mixing functions
 * for ints.
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
        return v ^ ((v >>> 23) ^ (v >>> 51));
    }

    /**
     * A variant of {@link #xnasam(long, long)} where
     * {@code 0x6a09e667f3bcc909L}
     * (<a href="https://dilbert.com/strip/2001-10-25">chosen at random</a>) is
     * used as a fixed value for the constant {@code c}.
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long xnasam(long v) {
        // first 64 bits of 1 + sqrt(2), forced to be odd
        return xnasam(v, 0x6a09e667f3bcc909L);
    }

    /**
     * Tommy Ettinger's {@code Pelican} mixer. See
     * "https://github.com/tommyettinger/sarong/blob/master/src/main/java/sarong/PelicanRNG.java".
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long pelican(long v) {
        return (v = ((v = (v ^ (v << 41 | v >>> 23) ^ (v << 17 | v >>> 47) ^ 0xd1b54a32d192ed03L) * 0xaef17502108ef2d9L)
                ^ v >>> 43 ^ v >>> 31 ^ v >>> 23) * 0xdb4f0b9175ae2165L) ^ v >>> 28;
    }

    /**
     * David Stafford's variant 13 of his 64-bit mixing functions. See
     * "http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html".
     * Note that if the argument {@code v} is 0, the result is 0.
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
     * Doug Lea's 64-bit mixing function. Note that if the argument {@code v} is
     * 0, the result is 0.
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long leaMix64(long v) {
        v = (v ^ (v >>> 32)) * 0xdaba0b6eb09322e3L;
        v = (v ^ (v >>> 32)) * 0xdaba0b6eb09322e3L;
        return v ^ (v >>> 32);
    }

    /**
     * Austin Appleby's fmix64() mix function used in {@code MurmurHash3}. See
     * "https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp".
     * Note that if the argument {@code v} is 0, the result is 0.
     * 
     * @param v
     *            long to mix
     * @return the mixed long
     */
    public static long applebyMix64(long v) {
        v ^= v >>> 33;
        v *= 0xff51afd7ed558ccdL;
        v ^= v >>> 33;
        v *= 0xc4ceb9fe1a85ec53L;
        return v ^ (v >>> 33);
    }

    /**
     * David Stafford's variant 4 of his 64-bit mixing functions adapted to
     * return the 32 high bits as an int. See
     * "http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html".
     * Note that if the argument {@code v} is 0, the result is 0.
     * 
     * @param v
     *            long to mix
     * @return a mixed int
     */
    public static int staffordMix04(long v) {
        v = (v ^ (v >>> 33)) * 0x62a9d9ed799705f5L;
        return (int) (((v ^ (v >>> 28)) * 0xcb24d0a5c88c35b3L) >>> 32);
    }

    /**
     * Doug Lea's 32-bit mixing function. Note that if the argument {@code v} is
     * 0, the result is 0.
     * 
     * @param v
     *            int to mix
     * @return the mixed int
     */
    public static int leaMix32(int v) {
        v = (v ^ (v >>> 16)) * 0xd36d884b;
        v = (v ^ (v >>> 16)) * 0xd36d884b;
        return v ^ (v >>> 16);
    }

    private BitMix() {
        throw new AssertionError();
    }
}
