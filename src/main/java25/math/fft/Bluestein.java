/* 
 * Copyright (c) 2017, 2018 Project Nayuki (MIT License) and Stefan Zobel (Apache 2.0)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */
/*
 * Any changes, bugfixes or additions made by the maintainers
 * of the https://github.com/stefan-zobel/FFT library are
 * licensed under the Apache License, Version 2.0, as explained
 * at http://www.apache.org/licenses/LICENSE-2.0
 */
package math.fft;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Bluestein chirp-z transform
 */
final class Bluestein {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    static ComplexArray forwardDFT(double[] data, double[] imag) {
        // find a power of 2 convolution length m such that m >= n * 2 + 1
        int n = data.length;
        if (n >= 0x20000000) {
            throw new IllegalArgumentException("array too large: " + n);
        }

        double[] cos = new double[n];
        double[] sin = new double[n];
        for (int i = 0; i < n; ++i) {
            int j = (int) ((long) i * i % (n * 2));
            double angle = Math.PI * j / n;
            cos[i] = Math.cos(angle);
            sin[i] = Math.sin(angle);
        }

        int m = Integer.highestOneBit(n) * 4;

        // temporary arrays
        double[] a_re = new double[m];
        double[] a_im = new double[m];
        double[] b_re = new double[m];
        double[] b_im = new double[m];

        b_re[0] = cos[0];
        b_im[0] = sin[0];

        int i = 0;
        int upper = SPECIES.loopBound(n);

        if (imag != null) {
            for (; i < upper; i += SPECIES.length()) {
                DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i);
                DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i);
                DoubleVector vRe = DoubleVector.fromArray(SPECIES, data, i);
                DoubleVector vIm = DoubleVector.fromArray(SPECIES, imag, i);

                vRe.fma(vCos, vIm.mul(vSin)).intoArray(a_re, i);
                vIm.fma(vCos, vRe.mul(vSin).neg()).intoArray(a_im, i);
            }
            for (; i < n; i += SPECIES.length()) {
                VectorMask<Double> vMask = SPECIES.indexInRange(i, n);
                DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i, vMask);
                DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i, vMask);
                DoubleVector vRe = DoubleVector.fromArray(SPECIES, data, i, vMask);
                DoubleVector vIm = DoubleVector.fromArray(SPECIES, imag, i, vMask);

                vRe.fma(vCos, vIm.mul(vSin)).intoArray(a_re, i, vMask);
                vIm.fma(vCos, vRe.mul(vSin).neg()).intoArray(a_im, i, vMask);
            }
        } else {
            for (; i < upper; i += SPECIES.length()) {
                DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i);
                DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i);
                DoubleVector vRe = DoubleVector.fromArray(SPECIES, data, i);

                vRe.mul(vCos).intoArray(a_re, i);
                vRe.mul(vSin).neg().intoArray(a_im, i);
            }
            for (; i < n; i += SPECIES.length()) {
                VectorMask<Double> vMask = SPECIES.indexInRange(i, n);
                DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i, vMask);
                DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i, vMask);
                DoubleVector vRe = DoubleVector.fromArray(SPECIES, data, i, vMask);

                vRe.mul(vCos).intoArray(a_re, i, vMask);
                vRe.mul(vSin).neg().intoArray(a_im, i, vMask);
            }
        }

        for (int k = 0; k < n; ++k) {
            if (k != 0) {
                double cos_i = cos[k];
                double sin_i = sin[k];
                b_re[k] = b_re[m - k] = cos_i;
                b_im[k] = b_im[m - k] = sin_i;
            }
        }

        // convolution
        ComplexArray conv = convolve(new ComplexArray(a_re, a_im, false), new ComplexArray(b_re, b_im, false));
        double[] c_re = conv.re();
        double[] c_im = conv.im();

        // result
        double[] re = new double[n];
        double[] im = new double[n];

        // postprocessing
        i = 0;
        for (; i < upper; i += SPECIES.length()) {
            DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i);
            DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i);
            DoubleVector vCre = DoubleVector.fromArray(SPECIES, c_re, i);
            DoubleVector vCim = DoubleVector.fromArray(SPECIES, c_im, i);

            DoubleVector vRe = vCre.fma(vCos, vCim.mul(vSin));
            DoubleVector vIm = vCim.fma(vCos, vCre.mul(vSin).neg());

            zeroSmall(vRe).intoArray(re, i);
            zeroSmall(vIm).intoArray(im, i);
        }
        for (; i < n; i += SPECIES.length()) {
            VectorMask<Double> vMask = SPECIES.indexInRange(i, n);
            DoubleVector vSin = DoubleVector.fromArray(SPECIES, sin, i, vMask);
            DoubleVector vCos = DoubleVector.fromArray(SPECIES, cos, i, vMask);
            DoubleVector vCre = DoubleVector.fromArray(SPECIES, c_re, i, vMask);
            DoubleVector vCim = DoubleVector.fromArray(SPECIES, c_im, i, vMask);

            DoubleVector vRe = vCre.fma(vCos, vCim.mul(vSin));
            DoubleVector vIm = vCim.fma(vCos, vCre.mul(vSin).neg());

            zeroSmall(vRe).intoArray(re, i, vMask);
            zeroSmall(vIm).intoArray(im, i, vMask);
        }

        return new ComplexArray(re, im, false);
    }

    static ComplexArray inverseDFT(ComplexArray freqs) {
        ComplexArray inv = forwardDFT(freqs.re(), freqs.im());
        double[] re = inv.re();
        double[] im = inv.im();
        final int n = re.length;
        int i = 0;
        int upper = SPECIES.loopBound(n);
        DoubleVector vInvN = DoubleVector.broadcast(SPECIES, 1.0 / n);

        for (; i < upper; i += SPECIES.length()) {
            DoubleVector vRe = DoubleVector.fromArray(SPECIES, re, i).mul(vInvN);
            DoubleVector vIm = DoubleVector.fromArray(SPECIES, im, i).mul(vInvN);
            zeroSmall(vRe).intoArray(re, i);
            zeroSmall(vIm).intoArray(im, i);
        }
        for (; i < n; i += SPECIES.length()) {
            VectorMask<Double> m = SPECIES.indexInRange(i, n);
            DoubleVector vRe = DoubleVector.fromArray(SPECIES, re, i, m).mul(vInvN);
            DoubleVector vIm = DoubleVector.fromArray(SPECIES, im, i, m).mul(vInvN);
            zeroSmall(vRe).intoArray(re, i, m);
            zeroSmall(vIm).intoArray(im, i, m);
        }

        for (int k = 1; k <= n / 2; ++k) {
            double re_tmp = re[n - k];
            double im_tmp = im[n - k];
            re[n - k] = re[k];
            re[k] = re_tmp;
            im[n - k] = im[k];
            im[k] = im_tmp;
        }
        return inv;
    }

    private static ComplexArray convolve(ComplexArray x, ComplexArray y) {

        x = Fourier.forwardDFT(x.re(), x.im());
        y = Fourier.forwardDFT(y.re(), y.im());

        double[] x_re = x.re();
        double[] x_im = x.im();
        double[] y_re = y.re();
        double[] y_im = y.im();

        int i = 0;
        int n = x_re.length;
        int upper = SPECIES.loopBound(n);

        for (; i < upper; i += SPECIES.length()) {
            DoubleVector vxRe = DoubleVector.fromArray(SPECIES, x_re, i);
            DoubleVector vxIm = DoubleVector.fromArray(SPECIES, x_im, i);
            DoubleVector vyRe = DoubleVector.fromArray(SPECIES, y_re, i);
            DoubleVector vyIm = DoubleVector.fromArray(SPECIES, y_im, i);

            DoubleVector rePart = vxRe.fma(vyRe, vxIm.mul(vyIm).neg());
            DoubleVector imPart = vxIm.fma(vyRe, vxRe.mul(vyIm));

            rePart.intoArray(x_re, i);
            imPart.intoArray(x_im, i);
        }
        for (; i < n; ++i) {
            double x_re_i = x_re[i];
            double y_re_i = y_re[i];
            double x_im_i = x_im[i];
            double y_im_i = y_im[i];
            x_re[i] = x_re_i * y_re_i - x_im_i * y_im_i;
            x_im[i] = x_im_i * y_re_i + x_re_i * y_im_i;
        }

        return Fourier.inverseDFT(new ComplexArray(x_re, x_im, false));
    }

    private static DoubleVector zeroSmall(DoubleVector v) {
        return v.blend(0.0, v.abs().compare(VectorOperators.LE, ComplexArray.TOL));
    }

    private Bluestein() {
        throw new AssertionError();
    }
}
