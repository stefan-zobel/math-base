/*
 * Copyright 2015 Stefan Zobel
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

import java.util.Objects;

/**
 * Mixing of long[] arrays based on Bob Jenkin's
 * <a href="https://burtleburtle.net/bob/hash/spooky.html">SpookyHash V2</a>
 * algorithm.
 */
final class SpookyMix {

    /**
     * Constants which:
     * <ul>
     * <li>are not zero
     * <li>are odd
     * <li>are a not-very-regular mix of 1's and 0's
     * <li>don't need any other special mathematical properties
     * </ul>
     */
    private static final long SC_CONST = 0xdeadbeefdeadbeefL;
    private static final long SC_SEED = 0xbb67ae8584caa73bL;

    /**
     * Computes a mixing for the array {@code src}, using the seeds provided in
     * {@code inout} and stores the resulting mix in {@code inout}. The arrays
     * {@code src} and {@code inout} must have the same length which must be
     * {@code > 0}.
     * 
     * @param src
     *            contains the long values to mix
     * @param inout
     *            on entry, contains the seeds to use; on exit, contains the
     *            mixed longs
     */
    static void mix(long[] src, long[] inout) {
        int srcLen = Objects.requireNonNull(src, "src").length;
        int outLen = Objects.requireNonNull(inout, "inout").length;
        if (srcLen != outLen) {
            throw new IllegalArgumentException("src.length != inout.length : " + srcLen + " != " + outLen);
        }
        if (srcLen == 0) {
            throw new IllegalArgumentException("length must be > 0 : 0");
        }
        if (srcLen <= 4) {
            mix4(src, 0, srcLen, inout, 0);
            return;
        }
        if (srcLen <= 12) {
            mix12(src, 0, srcLen, inout, 0);
            return;
        }
        for (int start = 0; start < srcLen; /**/) {
            int length = Math.min(srcLen - start, 12);
            mix12(src, start, length, inout, start);
            start += length;
        }
    }

    /**
     * Computes a mixing for an array of longs up to length {@code 4}, using up
     * to two seeds provided in {@code inout} starting at position {@code off}
     * and stores the resulting mix in {@code inout} starting at position
     * {@code off}.
     * 
     * @param src
     *            contains the long values to mix
     * @param start
     *            start index of the first element in {@code src} to consider
     * @param length
     *            number of long values to consider
     * @param inout
     *            on entry, up to two seeds starting at position {@code off}; on
     *            exit, {@code length} mixed longs starting at position
     *            {@code off}
     * @param off
     *            start offset in the {@code inout} array
     */
    private static void mix4(long[] src, int start, int length, long[] inout, int off) {
        int remaining = length;
        int pos = start;

        long h0 = inout[off + 0];
        long h1;
        long h2 = SC_CONST;
        long h3 = SC_CONST;

        int outLen = inout.length;
        if (outLen >= off + 2) {
            h1 = inout[off + 1];
        } else {
            h1 = SC_SEED;
        }

        while (remaining >= 4) {
            h2 += src[pos];
            h3 += src[pos + 1];

            h2 = (h2 << 50) | (h2 >>> 14);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 52) | (h3 >>> 12);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 30) | (h0 >>> 34);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 41) | (h1 >>> 23);
            h1 += h2;
            h3 ^= h1;
            h2 = (h2 << 54) | (h2 >>> 10);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 48) | (h3 >>> 16);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 38) | (h0 >>> 26);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 37) | (h1 >>> 27);
            h1 += h2;
            h3 ^= h1;
            h2 = (h2 << 62) | (h2 >>> 2);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 34) | (h3 >>> 30);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 5) | (h0 >>> 59);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 36) | (h1 >>> 28);
            h1 += h2;
            h3 ^= h1;

            h0 += src[pos + 2];
            h1 += src[pos + 3];
            pos += 4;
            remaining -= 4;
        }

        if (remaining >= 2) {
            h2 += src[pos];
            h3 += src[pos + 1];
            pos += 2;
            remaining -= 2;

            h2 = (h2 << 50) | (h2 >>> 14);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 52) | (h3 >>> 12);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 30) | (h0 >>> 34);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 41) | (h1 >>> 23);
            h1 += h2;
            h3 ^= h1;
            h2 = (h2 << 54) | (h2 >>> 10);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 48) | (h3 >>> 16);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 38) | (h0 >>> 26);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 37) | (h1 >>> 27);
            h1 += h2;
            h3 ^= h1;
            h2 = (h2 << 62) | (h2 >>> 2);
            h2 += h3;
            h0 ^= h2;
            h3 = (h3 << 34) | (h3 >>> 30);
            h3 += h0;
            h1 ^= h3;
            h0 = (h0 << 5) | (h0 >>> 59);
            h0 += h1;
            h2 ^= h0;
            h1 = (h1 << 36) | (h1 >>> 28);
            h1 += h2;
            h3 ^= h1;
        }

        h3 += ((long) (length << 3)) << 56;

        if (remaining > 0) {
            h2 += src[pos];
        } else {
            h2 += SC_CONST;
            h3 += SC_CONST;
        }

        h3 ^= h2;
        h2 = (h2 << 15) | (h2 >>> 49);
        h3 += h2;
        h0 ^= h3;
        h3 = (h3 << 52) | (h3 >>> 12);
        h0 += h3;
        h1 ^= h0;
        h0 = (h0 << 26) | (h0 >>> 38);
        h1 += h0;
        h2 ^= h1;
        h1 = (h1 << 51) | (h1 >>> 13);
        h2 += h1;
        h3 ^= h2;
        h2 = (h2 << 28) | (h2 >>> 36);
        h3 += h2;
        h0 ^= h3;
        h3 = (h3 << 9) | (h3 >>> 55);
        h0 += h3;
        h1 ^= h0;
        h0 = (h0 << 47) | (h0 >>> 17);
        h1 += h0;
        h2 ^= h1;
        h1 = (h1 << 54) | (h1 >>> 10);
        h2 += h1;
        h3 ^= h2;
        h2 = (h2 << 32) | (h2 >>> 32);
        h3 += h2;
        h0 ^= h3;
        h3 = (h3 << 25) | (h3 >>> 39);
        h0 += h3;
        h1 ^= h0;
        h0 = (h0 << 63) | (h0 >>> 1);
        h1 += h0;

        if (outLen >= off + 1) {
            inout[off + 0] = h0;
        }
        if (outLen >= off + 2) {
            inout[off + 1] = h1;
        }
        if (outLen >= off + 3) {
            inout[off + 2] = h2;
        }
        if (outLen >= off + 4) {
            inout[off + 3] = h3;
        }
    }

    /**
     * Computes a mixing for an array of longs up to length {@code 12}, using up
     * to two seeds provided in {@code inout} starting at position {@code off}
     * and stores the resulting mix in {@code inout} starting at position
     * {@code off}.
     * 
     * @param src
     *            contains the long values to mix
     * @param start
     *            start index of the first element in {@code src} to consider
     * @param length
     *            number of long values to consider
     * @param inout
     *            on entry, up to two seeds starting at position {@code off}; on
     *            exit, {@code length} mixed longs starting at position
     *            {@code off}
     * @param off
     *            start offset in the {@code inout} array
     */
    private static void mix12(long[] src, int start, int length, long[] inout, int off) {
        if (length <= 4) {
            mix4(src, start, length, inout, off);
            return;
        }
        int remaining = length;
        int pos = start;

        long seed0 = inout[off + 0];
        long seed1;

        int outLen = inout.length;
        if (outLen >= off + 2) {
            seed1 = inout[off + 1];
        } else {
            seed1 = SC_SEED;
        }

        long h0 = seed0;
        long h1 = seed1;
        long h2 = SC_CONST;
        long h3 = seed0;
        long h4 = seed1;
        long h5 = SC_CONST;
        long h6 = seed0;
        long h7 = seed1;
        long h8 = SC_CONST;
        long h9 = seed0;
        long h10 = seed1;
        long h11 = SC_CONST;

        // handle whole blocks of up to 12 values
        while (remaining > 0) {

            h0 += src[pos];
            h2 ^= h10;
            h11 ^= h0;
            h0 = (h0 << 11) | (h0 >>> 53);
            h11 += h1;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h1 += src[pos];
            h3 ^= h11;
            h0 ^= h1;
            h1 = (h1 << 32) | (h1 >>> 32);
            h0 += h2;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h2 += src[pos];
            h4 ^= h0;
            h1 ^= h2;
            h2 = (h2 << 43) | (h2 >>> 21);
            h1 += h3;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h3 += src[pos];
            h5 ^= h1;
            h2 ^= h3;
            h3 = (h3 << 31) | (h3 >>> 33);
            h2 += h4;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h4 += src[pos];
            h6 ^= h2;
            h3 ^= h4;
            h4 = (h4 << 17) | (h4 >>> 47);
            h3 += h5;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h5 += src[pos];
            h7 ^= h3;
            h4 ^= h5;
            h5 = (h5 << 28) | (h5 >>> 36);
            h4 += h6;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h6 += src[pos];
            h8 ^= h4;
            h5 ^= h6;
            h6 = (h6 << 39) | (h6 >>> 25);
            h5 += h7;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h7 += src[pos];
            h9 ^= h5;
            h6 ^= h7;
            h7 = (h7 << 57) | (h7 >>> 7);
            h6 += h8;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h8 += src[pos];
            h10 ^= h6;
            h7 ^= h8;
            h8 = (h8 << 55) | (h8 >>> 9);
            h7 += h9;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h9 += src[pos];
            h11 ^= h7;
            h8 ^= h9;
            h9 = (h9 << 54) | (h9 >>> 10);
            h8 += h10;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h10 += src[pos];
            h0 ^= h8;
            h9 ^= h10;
            h10 = (h10 << 22) | (h10 >>> 42);
            h9 += h11;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }

            h11 += src[pos];
            h1 ^= h9;
            h10 ^= h11;
            h11 = (h11 << 46) | (h11 >>> 18);
            h10 += h0;
            ++pos;
            --remaining;
            if (remaining == 0) {
                break;
            }
        }

        h11 += ((long) (length << 3)) << 56;

        for (int i = 0; i < 3; i++) {
            h11 += h1;
            h2 ^= h11;
            h1 = (h1 << 44) | (h1 >>> 20);
            h0 += h2;
            h3 ^= h0;
            h2 = (h2 << 15) | (h2 >>> 49);
            h1 += h3;
            h4 ^= h1;
            h3 = (h3 << 34) | (h3 >>> 30);
            h2 += h4;
            h5 ^= h2;
            h4 = (h4 << 21) | (h4 >>> 43);
            h3 += h5;
            h6 ^= h3;
            h5 = (h5 << 38) | (h5 >>> 26);
            h4 += h6;
            h7 ^= h4;
            h6 = (h6 << 33) | (h6 >>> 31);
            h5 += h7;
            h8 ^= h5;
            h7 = (h7 << 10) | (h7 >>> 54);
            h6 += h8;
            h9 ^= h6;
            h8 = (h8 << 13) | (h8 >>> 51);
            h7 += h9;
            h10 ^= h7;
            h9 = (h9 << 38) | (h9 >>> 26);
            h8 += h10;
            h11 ^= h8;
            h10 = (h10 << 53) | (h10 >>> 11);
            h9 += h11;
            h0 ^= h9;
            h11 = (h11 << 42) | (h11 >>> 22);
            h10 += h0;
            h1 ^= h10;
            h0 = (h0 << 54) | (h0 >>> 10);
        }

        if (outLen >= off + 1) {
            inout[off + 0] = h0;
        }
        if (outLen >= off + 2) {
            inout[off + 1] = h1;
        }
        if (outLen >= off + 3) {
            inout[off + 2] = h2;
        }
        if (outLen >= off + 4) {
            inout[off + 3] = h3;
        }
        if (outLen >= off + 5) {
            inout[off + 4] = h4;
        }
        if (outLen >= off + 6) {
            inout[off + 5] = h5;
        }
        if (outLen >= off + 7) {
            inout[off + 6] = h6;
        }
        if (outLen >= off + 8) {
            inout[off + 7] = h7;
        }
        if (outLen >= off + 9) {
            inout[off + 8] = h8;
        }
        if (outLen >= off + 10) {
            inout[off + 9] = h9;
        }
        if (outLen >= off + 11) {
            inout[off + 10] = h10;
        }
        if (outLen >= off + 12) {
            inout[off + 11] = h11;
        }
    }

    private SpookyMix() {
        throw new AssertionError();
    }
}
