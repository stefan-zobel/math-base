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
 * Determines the index of the largest element in a given part of a float[] or
 * double[] array.
 */
public final class ArrayMaxPos {

    /**
     * Finds the index of the largest element in the passed array slice.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param off
     *            starting position in {@code a}
     * @param a
     *            the array to inspect
     * @return index of the largest element in the given array slice
     */
    public static int maxPos(int length, int off, double[] a) {
        double max = a[off];
        int pos = off;
        for (int i = off + 1; i < off + length; ++i) {
            double b = a[i];
            if (b > max) {
                max = b;
                pos = i;
            }
        }
        return pos;
    }

    /**
     * Finds the index of the largest element in the passed array slice.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param off
     *            starting position in {@code a}
     * @param a
     *            the array to inspect
     * @return index of the largest element in the given array slice
     */
    public static int maxPosF(int length, int off, float[] a) {
        float max = a[off];
        int pos = off;
        for (int i = off + 1; i < off + length; ++i) {
            float b = a[i];
            if (b > max) {
                max = b;
                pos = i;
            }
        }
        return pos;
    }

    private ArrayMaxPos() {
        throw new AssertionError();
    }
}
