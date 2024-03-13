/*
 * Copyright 2024 Stefan Zobel
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
package math.dl;

/**
 * The softmax function in deep learning.
 */
public final class Softmax {

    /**
     * Determine the softmax distribution for the {@code in} array partition of
     * size {@code length} starting at {@code ioff} and store it in the
     * {@code out} array partition of size {@code length} starting at
     * {@code ooff}.
     * 
     * @param length
     *            the size of the array slices in {@code in} and {@code out}
     * @param ioff
     *            starting position in {@code in}
     * @param in
     *            the source array
     * @param ooff
     *            starting position in {@code out}
     * @param out
     *            the destination array
     * @return the destination array
     */
    public static double[] softmax(int length, int ioff, double[] in, int ooff, double[] out) {
        double max = max(length, ioff, in);
        double s = 0.0;
        for (int i = ioff; i < ioff + length; ++i) {
            double q = Math.exp(in[i] - max);
            s += q;
            out[ooff + i - ioff] = q;
        }
        s = 1.0 / s;
        for (int i = ooff; i < ooff + length; ++i) {
            out[i] *= s;
        }
        return out;
    }

    /**
     * Determine the softmax distribution for the {@code in} array partition of
     * size {@code length} starting at {@code ioff} and store it in the
     * {@code out} array partition of size {@code length} starting at
     * {@code ooff}.
     * 
     * @param length
     *            the size of the array slices in {@code in} and {@code out}
     * @param ioff
     *            starting position in {@code in}
     * @param in
     *            the source array
     * @param ooff
     *            starting position in {@code out}
     * @param out
     *            the destination array
     * @return the destination array
     */
    public static float[] softmaxF(int length, int ioff, float[] in, int ooff, float[] out) {
        float max = maxF(length, ioff, in);
        float s = 0.0f;
        for (int i = ioff; i < ioff + length; ++i) {
            float q = (float) Math.exp(in[i] - max);
            s += q;
            out[ooff + i - ioff] = q;
        }
        s = 1.0f / s;
        for (int i = ooff; i < ooff + length; ++i) {
            out[i] *= s;
        }
        return out;
    }

    static double max(int length, int aoff, double[] a) {
        double max = a[aoff];
        for (int i = aoff + 1; i < aoff + length; ++i) {
            max = Math.max(max, a[i]);
        }
        return max;
    }

    static float maxF(int length, int aoff, float[] a) {
        float max = a[aoff];
        for (int i = aoff + 1; i < aoff + length; ++i) {
            max = Math.max(max, a[i]);
        }
        return max;
    }

    private Softmax() {
        throw new AssertionError();
    }
}
