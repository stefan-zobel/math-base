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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import math.cern.Arithmetic;

/**
 * Calculation of the accuracy of a multi-class classifier.
 */
public class MultiClassAccuracy implements Externalizable {

    private static final long serialVersionUID = 1L;

    private long ok = 0L;
    private long total = 0L;

    public MultiClassAccuracy() {
    }

    /**
     * Get the average accuracy achieved so far.
     * 
     * @return the current accuracy measure
     */
    public double getAverage() {
        if (total == 0L) {
            return 0.0;
        }
        return Arithmetic.round(ok / (double) total, 4);
    }

    /**
     * Reset the accuracy measure to {@code 0.0}.
     */
    public void reset() {
        ok = 0L;
        total = 0L;
    }

    /**
     * Compare the two discrete probability distributions contained in the array
     * slices {@code a} and {@code b} of size {@code length} both starting at
     * position {@code off} for a match of the positions containing the highest
     * probability. If both positions match this is interpreted as agreement on
     * the class labels, otherwise the two distributions disagree.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param a
     *            array slice containing the first probability distribution
     * @param off
     *            starting position in both {@code a} and {@code b}
     * @param b
     *            array slice containing the second probability distribution
     */
    public void compare(int length, double[] a, int off, double[] b) {
        compare(length, off, a, off, b);
    }

    /**
     * Compare the two discrete probability distributions contained in the array
     * slices {@code a} and {@code b} of size {@code length} both starting at
     * position {@code off} for a match of the positions containing the highest
     * probability. If both positions match this is interpreted as agreement on
     * the class labels, otherwise the two distributions disagree.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param a
     *            array slice containing the first probability distribution
     * @param off
     *            starting position in both {@code a} and {@code b}
     * @param b
     *            array slice containing the second probability distribution
     */
    public void compareF(int length, float[] a, int off, float[] b) {
        compareF(length, off, a, off, b);
    }

    /**
     * Compare the two discrete probability distributions contained in the array
     * slices {@code a} and {@code b} of size {@code length} starting at
     * positions {@code aoff} and {@code boff} for a match of the positions
     * containing the highest probability. If both positions relative to the
     * offsets match this is interpreted as agreement on the class labels,
     * otherwise the two distributions disagree.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param aoff
     *            starting position in array {@code a}
     * @param a
     *            array slice containing the first probability distribution
     * @param boff
     *            starting position in array {@code b}
     * @param b
     *            array slice containing the second probability distribution
     */
    public void compare(int length, int aoff, double[] a, int boff, double[] b) {
        if (ArrayMaxPos.maxPos(length, aoff, a) - aoff == ArrayMaxPos.maxPos(length, boff, b) - boff) {
            ++ok;
        }
        ++total;
    }

    /**
     * Compare the two discrete probability distributions contained in the array
     * slices {@code a} and {@code b} of size {@code length} starting at
     * positions {@code aoff} and {@code boff} for a match of the positions
     * containing the highest probability. If both positions relative to the
     * offsets match this is interpreted as agreement on the class labels,
     * otherwise the two distributions disagree.
     * 
     * @param length
     *            the size of the array slice to consider
     * @param aoff
     *            starting position in array {@code a}
     * @param a
     *            array slice containing the first probability distribution
     * @param boff
     *            starting position in array {@code b}
     * @param b
     *            array slice containing the second probability distribution
     */
    public void compareF(int length, int aoff, float[] a, int boff, float[] b) {
        if (ArrayMaxPos.maxPosF(length, aoff, a) - aoff == ArrayMaxPos.maxPosF(length, boff, b) - boff) {
            ++ok;
        }
        ++total;
    }

    @Override
    public void writeExternal(ObjectOutput out) throws IOException {
        out.writeLong(ok);
        out.writeLong(total);
    }

    @Override
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        ok = in.readLong();
        total = in.readLong();
    }
}
