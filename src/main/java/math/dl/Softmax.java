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

import math.rng.Stc64;

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

    /**
     * Reweighs a softmax distribution by the positive temperature {@code t}. A
     * temperature of {@code 1.0} does nothing to the distribution. When
     * {@code t} gets above {@code 1.0} smaller probabilities get pushed up
     * while higher probabilities get pushed down (approaching a discrete
     * uniform distribution when {@code t} approaches Double.MAX_VALUE. When
     * {@code t} is below {@code 1.0} higher probabilities get pushed up and
     * lower probabilities get pushed down (approaching a one-hot encoded
     * distribution where the class with the highest probability gets all the
     * probability when {@code t} approaches {@code 0.0}). The latter extreme
     * case should be avoided since the computation becomes numerically unstable
     * when t approaches {@code 0.0} too closely - it is better to avoid values
     * for {@code t} below {@code 0.05}.
     * 
     * @param softmax
     *            a softmax (discrete probability) distribution
     * @param t
     *            temperature coefficient to use (t must be > 0.0)
     * @return a reweighted distribution (to a different temperature when
     *         {@code t != 1.0})
     * @throws IllegalArgumentException
     *             if t is zero or negative or if one of the probabilities in
     *             the softmax distribution is negative
     */
    public static double[] reweigh(double[] softmax, double t) {
        if (t == 1.0) {
            return softmax;
        }
        if (t <= 0.0) {
            throw new IllegalArgumentException("temperature must be strictly positive: " + t);
        }
        double[] b = new double[softmax.length];
        double sum = 0.0;
        double max = -Double.MAX_VALUE;
        int maxpos = 0;
        for (int i = 0; i < softmax.length; ++i) {
            double p = Math.log(clamp(softmax[i])) / t;
            if (p > max) {
                max = p;
                maxpos = i;
            }
            p = Math.exp(p);
            sum += p;
            b[i] = p;
        }
        if (sum == 0.0) {
            b[maxpos] = 1.0;
        } else {
            for (int i = 0; i < b.length; ++i) {
                b[i] /= sum;
            }
        }
        return b;
    }

    /**
     * Reweighs a softmax distribution by the positive temperature {@code t}. A
     * temperature of {@code 1.0} does nothing to the distribution. When
     * {@code t} gets above {@code 1.0} smaller probabilities get pushed up
     * while higher probabilities get pushed down (approaching a discrete
     * uniform distribution when {@code t} approaches Float.MAX_VALUE. When
     * {@code t} is below {@code 1.0} higher probabilities get pushed up and
     * lower probabilities get pushed down (approaching a one-hot encoded
     * distribution where the class with the highest probability gets all the
     * probability when {@code t} approaches {@code 0.0}). The latter extreme
     * case should be avoided since the computation becomes numerically unstable
     * when t approaches {@code 0.0} too closely - it is better to avoid values
     * for {@code t} below {@code 0.05}.
     * 
     * @param softmax
     *            a softmax (discrete probability) distribution
     * @param t
     *            temperature coefficient to use (t must be > 0.0)
     * @return a reweighted distribution (to a different temperature when
     *         {@code t != 1.0})
     * @throws IllegalArgumentException
     *             if t is zero or negative or if one of the probabilities in
     *             the softmax distribution is negative
     */
    public static float[] reweighF(float[] softmax, float t) {
        if (t == 1.0f) {
            return softmax;
        }
        if (t <= 0.0f) {
            throw new IllegalArgumentException("temperature must be strictly positive: " + t);
        }
        float[] b = new float[softmax.length];
        float sum = 0.0f;
        float max = -Float.MAX_VALUE;
        int maxpos = 0;
        for (int i = 0; i < softmax.length; ++i) {
            float p = (float) Math.log(clampF(softmax[i])) / t;
            if (p > max) {
                max = p;
                maxpos = i;
            }
            p = (float) Math.exp(p);
            sum += p;
            b[i] = p;
        }
        if (sum == 0.0) {
            b[maxpos] = 1.0f;
        } else {
            for (int i = 0; i < b.length; ++i) {
                b[i] /= sum;
            }
        }
        return b;
    }

    /**
     * Samples a class index from the softmax distribution randomly according to
     * the class probabilities in the {@code softmax} array.
     * 
     * @param softmax
     *            a softmax (discrete probability) distribution
     * @return class index {@code i} with probability {@code softmax[i]}
     */
    public static int sampleClass(double[] softmax) {
        int classIdx = 0;
        double p = Stc64.getDefault().nextDouble();
        // roulette wheel selection
        for (int i = 0; i < softmax.length; ++i) {
            p -= softmax[i];
            if (p <= 0.0) {
                classIdx = i;
                break;
            }
        }
        return classIdx;
    }

    /**
     * Samples a class index from the softmax distribution randomly according to
     * the class probabilities in the {@code softmax} array.
     * 
     * @param softmax
     *            a softmax (discrete probability) distribution
     * @return class index {@code i} with probability {@code softmax[i]}
     */
    public static int sampleClassF(float[] softmax) {
        int classIdx = 0;
        float p = Stc64.getDefault().nextFloat();
        // roulette wheel selection
        for (int i = 0; i < softmax.length; ++i) {
            p -= softmax[i];
            if (p <= 0.0f) {
                classIdx = i;
                break;
            }
        }
        return classIdx;
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

    static double clamp(double p) {
        if (p > 0.0) {
            return p;
        }
        // p should never be <= 0
        if (p == 0.0) {
            return Double.MIN_NORMAL;
        }
        throw new IllegalArgumentException("negative probability: " + p);
    }

    static float clampF(float p) {
        if (p > 0.0f) {
            return p;
        }
        // p should never be <= 0
        if (p == 0.0f) {
            return Float.MIN_NORMAL;
        }
        throw new IllegalArgumentException("negative probability: " + p);
    }

    private Softmax() {
        throw new AssertionError();
    }
}
