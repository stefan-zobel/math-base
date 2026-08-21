/*
 * Copyright 2018, 2026 Stefan Zobel
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
package math.fft;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import math.MathConsts;

/**
 * For computations with arrays of complex numbers.
 * <p>
 * Every operation here that produces an array flushes the round-off it
 * generated to exact zeros, so that a coefficient the mathematics makes zero
 * comes back as {@code 0.0}. The threshold for that is <em>relative</em> to the
 * largest magnitude in the same result -- an absolute one would mean something
 * only for data that happens to be of order one, and would silently delete the
 * whole result of a transform over data that is not. See
 * {@link #relativeFactor(int)}.
 * <p>
 * Note that indexes in {@link #set(int, double, double)} are 1-based!
 */
public final class ComplexArray {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    private final double[] re;
    private final double[] im;

    public ComplexArray(int size) {
        if (size < 0) {
            throw new IllegalArgumentException("size < 0 : " + size);
        }
        re = new double[size];
        im = new double[size];
    }

    public ComplexArray(double[] re) {
        this(re.clone(), new double[re.length], false);
    }

    public ComplexArray(double[] re, double[] im) {
        this(re, im, true);
    }

    public ComplexArray(double[] re, double[] im, boolean copy) {
        if (re.length != im.length) {
            throw new IllegalArgumentException(re.length + " != " + im.length);
        }
        if (copy) {
            this.re = re.clone();
            this.im = im.clone();
        } else {
            this.re = re;
            this.im = im;
        }
    }

    public void set(int index, double re, double im) {
        checkArg(index);
        this.re[index - 1] = re;
        this.im[index - 1] = im;
    }

    public ComplexArray naiveForwardDFT() {
        return naiveDFT(-1.0, re, im, 1.0);
    }

    public ComplexArray naiveInverseDFT() {
        return naiveDFT(1.0, re, im, (1.0 / re.length));
    }

    private static ComplexArray naiveDFT(double sign, double[] re, double[] im, double scale) {
        int N = re.length;
        double[] imag = new double[N];
        double[] real = new double[N];
        double[] cos = new double[N];
        double[] sin = new double[N];
        for (int i = 0; i < N; ++i) {
            double angle = (sign * MathConsts.TWO_PI * i) / N;
            cos[i] = Math.cos(angle);
            sin[i] = Math.sin(angle);
        }
        for (int i = 0; i < N; ++i) {
            double rZ = 0.0;
            double iZ = 0.0;
            for (long j = 0; j < N; ++j) {
                int idx = (int) ((i * j) % N);
                double cosine = cos[idx];
                double sine = sin[idx];
                double rX = re[(int) j];
                double iY = (im == null) ? 0.0 : im[(int) j];
                rZ += cosine * rX - sine * iY;
                iZ += sine * rX + cosine * iY;
            }
            real[i] = scale * rZ;
            imag[i] = scale * iZ;
        }
        zeroNegligible(real, imag, largestMagnitude(real, imag), N);
        return new ComplexArray(real, imag, false);
    }

    public static ComplexArray naiveForwarDFT(double[] data) {
        return naiveDFT(-1.0, data, null, 1.0);
    }

    public static ComplexArray naiveInverseDFT(ComplexArray freqs) {
        return naiveDFT(1.0, freqs.re, freqs.im, (1.0 / freqs.re.length));
    }

    public double[] absSquared() {
        return absSquaredScaled(false);
    }

    // for power density spectrum
    public double[] absSquaredScaled() {
        return absSquaredScaled(true);
    }

    private double[] absSquaredScaled(boolean withScaling) {
        double[] real = re;
        double[] imag = im;
        int N = real.length;
        double[] res = new double[N];
        double invScale = withScaling ? (1.0 / N) : 1.0;

        int i = 0;
        int upper = SPECIES.loopBound(N);
        DoubleVector vInvScale = DoubleVector.broadcast(SPECIES, invScale);

        for (; i < upper; i += SPECIES.length()) {
            DoubleVector vr = DoubleVector.fromArray(SPECIES, real, i);
            DoubleVector vi = DoubleVector.fromArray(SPECIES, imag, i);

            DoubleVector square = vr.fma(vr, vi.mul(vi)).mul(vInvScale);
            square.intoArray(res, i);
        }
        for (; i < N; i += SPECIES.length()) {
            VectorMask<Double> m = SPECIES.indexInRange(i, N);
            DoubleVector vr = DoubleVector.fromArray(SPECIES, real, i, m);
            DoubleVector vi = DoubleVector.fromArray(SPECIES, imag, i, m);

            DoubleVector square = vr.fma(vr, vi.mul(vi)).mul(vInvScale);
            square.intoArray(res, i, m);
        }
        // these are squares, so the threshold is the square of the one that
        // applies to a magnitude: a coefficient that is nothing but round-off
        // is at most factor * max|X|, and its square at most factor^2 * maxSquare
        double max = 0.0;
        for (int k = 0; k < N; ++k) {
            max = Math.max(max, res[k]);
        }
        if (max > 0.0 && !Double.isInfinite(max)) {
            double factor = relativeFactor(N);
            double cutoff = max * factor * factor;
            blendBelow(res, cutoff);
        }
        return res;
    }

    public ComplexArray fftshift() {
        return shift(false);
    }

    public ComplexArray ifftshift() {
        return shift(true);
    }

    private ComplexArray shift(boolean inverse) {
        final int length = re.length;
        int mid = -1;
        double[] re_this = re;
        double[] im_this = im;
        double[] re_shift = new double[length];
        double[] im_shift = new double[length];
        if (length % 2 == 0) {
            mid = (length / 2);
            System.arraycopy(re_this, 0, re_shift, mid, mid);
            System.arraycopy(re_this, mid, re_shift, 0, mid);
            System.arraycopy(im_this, 0, im_shift, mid, mid);
            System.arraycopy(im_this, mid, im_shift, 0, mid);
        } else {
            mid = (length - 1) / 2;
            if (inverse) {
                System.arraycopy(re_this, 0, re_shift, mid + 1, mid);
                System.arraycopy(re_this, mid, re_shift, 0, mid + 1);
                System.arraycopy(im_this, 0, im_shift, mid + 1, mid);
                System.arraycopy(im_this, mid, im_shift, 0, mid + 1);
            } else {
                System.arraycopy(re_this, 0, re_shift, mid, mid + 1);
                System.arraycopy(re_this, mid + 1, re_shift, 0, mid);
                System.arraycopy(im_this, 0, im_shift, mid, mid + 1);
                System.arraycopy(im_this, mid + 1, im_shift, 0, mid);
            }
        }
        return new ComplexArray(re_shift, im_shift, false);
    }

    public static double[] dot(ComplexArray a, ComplexArray b) {
        if (a.length() != b.length()) {
            throw new IllegalArgumentException("Unequal dimensions: " + a.length() + " != " + b.length());
        }
        if (a.length() == 0) {
            throw new IllegalArgumentException("Arrays are empty: length = 0");
        }

        double[] a_re_ = a.re;
        double[] b_re_ = b.re;
        double[] a_im_ = a.im;
        double[] b_im_ = b.im;

        int i = 0;
        int n = a_re_.length;
        int upper = SPECIES.loopBound(n);

        DoubleVector accRe = DoubleVector.zero(SPECIES);
        DoubleVector accIm = DoubleVector.zero(SPECIES);
        DoubleVector vMax = DoubleVector.zero(SPECIES);

        for (; i < upper; i += SPECIES.length()) {
            DoubleVector aRe = DoubleVector.fromArray(SPECIES, a_re_, i);
            DoubleVector bRe = DoubleVector.fromArray(SPECIES, b_re_, i);
            DoubleVector aIm = DoubleVector.fromArray(SPECIES, a_im_, i);
            DoubleVector bIm = DoubleVector.fromArray(SPECIES, b_im_, i);

            DoubleVector rePart = aRe.fma(bRe, aIm.mul(bIm).neg());
            DoubleVector imPart = aRe.fma(bIm, aIm.mul(bRe));

            // the terms themselves are not cleaned: zeroing an addend changes
            // the sum, and nobody ever sees an intermediate
            accRe = accRe.add(rePart);
            accIm = accIm.add(imPart);
            vMax = vMax.max(rePart.abs()).max(imPart.abs());
        }
        for (; i < n; i += SPECIES.length()) {
            VectorMask<Double> m = SPECIES.indexInRange(i, n);
            DoubleVector aRe = DoubleVector.fromArray(SPECIES, a_re_, i, m);
            DoubleVector bRe = DoubleVector.fromArray(SPECIES, b_re_, i, m);
            DoubleVector aIm = DoubleVector.fromArray(SPECIES, a_im_, i, m);
            DoubleVector bIm = DoubleVector.fromArray(SPECIES, b_im_, i, m);

            DoubleVector rePart = aRe.fma(bRe, aIm.mul(bIm).neg());
            DoubleVector imPart = aRe.fma(bIm, aIm.mul(bRe));

            accRe = accRe.add(rePart, m);
            accIm = accIm.add(imPart, m);
            // inactive lanes hold whatever the masked load left, so they are
            // blended away before they can raise the maximum
            vMax = vMax.max(rePart.abs().blend(0.0, m.not())).max(imPart.abs().blend(0.0, m.not()));
        }

        double res_re = accRe.reduceLanes(VectorOperators.ADD);
        double res_im = accIm.reduceLanes(VectorOperators.ADD);

        // the error of a sum of n terms is bounded by n times the largest of
        // them; a maximum rather than the sum of magnitudes, because it does
        // not depend on the order the terms were visited in
        double[] result = { res_re, res_im };
        zeroNegligible(result, new double[2], n * vMax.reduceLanes(VectorOperators.MAX), n);
        return result;
    }

    public static ComplexArray elementwiseProduct(ComplexArray a, ComplexArray b) {
        if (a.length() != b.length()) {
            throw new IllegalArgumentException("Unequal dimensions: " + a.length() + " != " + b.length());
        }

        int n = a.length();
        double[] real = new double[n];
        double[] imag = new double[n];
        double[] a_re_ = a.re;
        double[] b_re_ = b.re;
        double[] a_im_ = a.im;
        double[] b_im_ = b.im;

        int i = 0;
        int upper = SPECIES.loopBound(n);

        for (; i < upper; i += SPECIES.length()) {
            DoubleVector aRe = DoubleVector.fromArray(SPECIES, a_re_, i);
            DoubleVector bRe = DoubleVector.fromArray(SPECIES, b_re_, i);
            DoubleVector aIm = DoubleVector.fromArray(SPECIES, a_im_, i);
            DoubleVector bIm = DoubleVector.fromArray(SPECIES, b_im_, i);

            DoubleVector rePart = aRe.fma(bRe, aIm.mul(bIm).neg());
            DoubleVector imPart = aRe.fma(bIm, aIm.mul(bRe));

            rePart.intoArray(real, i);
            imPart.intoArray(imag, i);
        }
        for (; i < n; i += SPECIES.length()) {
            VectorMask<Double> m = SPECIES.indexInRange(i, n);
            DoubleVector aRe = DoubleVector.fromArray(SPECIES, a_re_, i, m);
            DoubleVector bRe = DoubleVector.fromArray(SPECIES, b_re_, i, m);
            DoubleVector aIm = DoubleVector.fromArray(SPECIES, a_im_, i, m);
            DoubleVector bIm = DoubleVector.fromArray(SPECIES, b_im_, i, m);

            DoubleVector rePart = aRe.fma(bRe, aIm.mul(bIm).neg());
            DoubleVector imPart = aRe.fma(bIm, aIm.mul(bRe));

            rePart.intoArray(real, i, m);
            imPart.intoArray(imag, i, m);
        }
        zeroNegligible(real, imag, largestMagnitude(real, imag), n);

        return new ComplexArray(real, imag, false);
    }

    public double[] re() {
        return re;
    }

    public double[] im() {
        return im;
    }

    public int length() {
        return re.length;
    }

    @Override
    public String toString() {
        int max = length() - 1;
        if (max == -1) {
            return "[]";
        }
        StringBuilder b = new StringBuilder(40 * (max + 1));
        b.append('[');
        for (int i = 0; ; i++) {
            b.append(re[i]).append("  ").append(im[i]).append('i');
            if (i == max) {
                return b.append(']').toString();
            }
            b.append(",\n ");
        }
    }

    private void checkArg(int idx) {
        if (idx < 1 || idx > re.length) {
            throw new IllegalArgumentException("Invalid index " + idx + " for [1.." + re.length + "] array");
        }
    }

    /**
     * The relative accuracy an {@code n} point transform or an {@code n} term
     * sum can be held to: anything at or below this fraction of the largest
     * magnitude in the same result is round-off and not information.
     * <p>
     * The round-off floor of both transform paths grows in proportion to
     * {@code n}. Measured over 24000 transforms, radix-2 and Bluestein, lengths
     * from 256 to 16384 and one to five spectral lines, the floor stayed below
     * {@code 1.6 * n * eps} of the maximum, with a median near
     * {@code 0.3 * n * eps}; the factor of four is the margin over the worst
     * of those. Making it larger would cost dynamic range: a genuine
     * coefficient below this fraction of the largest one is indistinguishable
     * from round-off here and is lost.
     *
     * @param n
     *            the number of points the result was computed from
     * @return the relative threshold, as a fraction of the largest magnitude
     */
    static double relativeFactor(int n) {
        return 4.0 * n * MathConsts.MACH_EPS_DBL;
    }

    /**
     * The largest magnitude across both components. It is a maximum and not a
     * sum, so it is exact and independent of the order of traversal, which is
     * what lets this class and its Java 8 counterpart agree. A {@code NaN}
     * anywhere propagates, which switches the cleanup off rather than letting
     * it act on a meaningless reference.
     *
     * @param re
     *            the real parts
     * @param im
     *            the imaginary parts, of the same length
     * @return the largest of {@code |re[i]|} and {@code |im[i]|}
     */
    static double largestMagnitude(double[] re, double[] im) {
        final int n = re.length;
        int i = 0;
        int upper = SPECIES.loopBound(n);
        DoubleVector vMax = DoubleVector.zero(SPECIES);
        for (; i < upper; i += SPECIES.length()) {
            vMax = vMax.max(DoubleVector.fromArray(SPECIES, re, i).abs());
            vMax = vMax.max(DoubleVector.fromArray(SPECIES, im, i).abs());
        }
        double max = vMax.reduceLanes(VectorOperators.MAX);
        // the tail is scalar so that it matches the Java 8 form exactly; a
        // masked reduction would have to decide what the inactive lanes
        // contribute, and this way there is nothing to decide
        for (; i < n; ++i) {
            max = Math.max(max, Math.abs(re[i]));
            max = Math.max(max, Math.abs(im[i]));
        }
        return max;
    }

    /**
     * Zeroes every entry at or below {@code reference * relativeFactor(n)}. A
     * reference that is zero or not finite leaves both arrays untouched: there
     * is then nothing to be relative to, and an infinite cutoff would zero
     * everything.
     *
     * @param re
     *            the real parts, modified in place
     * @param im
     *            the imaginary parts, modified in place
     * @param reference
     *            the largest magnitude in the result
     * @param n
     *            the number of points the result was computed from
     */
    static void zeroNegligible(double[] re, double[] im, double reference, int n) {
        if (!(reference > 0.0) || Double.isInfinite(reference)) {
            return;
        }
        final double cutoff = reference * relativeFactor(n);
        blendBelow(re, cutoff);
        blendBelow(im, cutoff);
    }

    private static void blendBelow(double[] a, double cutoff) {
        final int n = a.length;
        int i = 0;
        int upper = SPECIES.loopBound(n);
        for (; i < upper; i += SPECIES.length()) {
            DoubleVector v = DoubleVector.fromArray(SPECIES, a, i);
            v.blend(0.0, v.abs().compare(VectorOperators.LE, cutoff)).intoArray(a, i);
        }
        for (; i < n; ++i) {
            if (Math.abs(a[i]) <= cutoff) {
                a[i] = 0.0;
            }
        }
    }
}
