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

import math.MathConsts;

/**
 * Similarity calculation between two discrete probability distributions based
 * on the Hellinger distance.
 */
public final class Hellinger {

    private static final double ONE_OVER_SQRT_TWO = 1.0 / MathConsts.SQRT_TWO;

    /**
     * Computes the similarity between two discrete probability distributions
     * contained in the array slices {@code a} and {@code b} of size
     * {@code length} both starting at position {@code off}.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param a
     *            array containing the first probability distribution
     * @param off
     *            starting position in both {@code a} and {@code b}
     * @param b
     *            array containing the second probability distribution
     * @return a similarity measure in {@code [0..1]}
     */
    public static double similarity(int length, double[] a, int off, double[] b) {
        return similarity(length, off, a, off, b);
    }

    /**
     * Computes the similarity between two discrete probability distributions
     * contained in the array slices {@code a} and {@code b} of size
     * {@code length} both starting at position {@code off}.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param a
     *            array containing the first probability distribution
     * @param off
     *            starting position in both {@code a} and {@code b}
     * @param b
     *            array containing the second probability distribution
     * @return a similarity measure in {@code [0..1]}
     */
    public static float similarityF(int length, float[] a, int off, float[] b) {
        return similarityF(length, off, a, off, b);
    }

    /**
     * Computes the similarity between two discrete probability distributions
     * contained in the array slices {@code a} and {@code b} of size
     * {@code length} starting at positions {@code aoff} and {@code boff}
     * respectively.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param aoff
     *            starting position in array {@code a}
     * @param a
     *            array containing the first probability distribution
     * @param boff
     *            starting position in array {@code b}
     * @param b
     *            array containing the second probability distribution
     * @return a similarity measure in {@code [0..1]}
     */
    public static double similarity(int length, int aoff, double[] a, int boff, double[] b) {
        double dist = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            double x = Math.sqrt(a[i]) - Math.sqrt(b[boff + i - aoff]);
            x *= x;
            dist += x;
        }
        return 1.0 - ONE_OVER_SQRT_TWO * Math.sqrt(dist);
    }

    /**
     * Computes the similarity between two discrete probability distributions
     * contained in the array slices {@code a} and {@code b} of size
     * {@code length} starting at positions {@code aoff} and {@code boff}
     * respectively.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param aoff
     *            starting position in array {@code a}
     * @param a
     *            array containing the first probability distribution
     * @param boff
     *            starting position in array {@code b}
     * @param b
     *            array containing the second probability distribution
     * @return a similarity measure in {@code [0..1]}
     */
    public static float similarityF(int length, int aoff, float[] a, int boff, float[] b) {
        double dist = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            double x = Math.sqrt(a[i]) - Math.sqrt(b[boff + i - aoff]);
            x *= x;
            dist += x;
        }
        return (float) (1.0 - ONE_OVER_SQRT_TWO * Math.sqrt(dist));
    }

    private Hellinger() {
        throw new AssertionError();
    }
}
