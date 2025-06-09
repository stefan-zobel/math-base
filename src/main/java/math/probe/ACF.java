/*
 * Copyright 2025 Stefan Zobel
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
package math.probe;

import math.cern.Arithmetic;
import math.fft.ComplexArray;
import math.fft.Fourier;

public final class ACF {

    public static double[] acf(double[] data) {
        return acf(data, -1);
    }

    public static double[] acf(double[] data, int lags) {
        if (data == null) {
            return null;
        }
        final int length = data.length;
        if (length == 0) {
            return new double[0];
        }
        if (length == 1 || lags == 0) {
            return new double[] { 1.0 };
        }
        double sum = sum(data);
        double avg = sum / length;
        data = centerAndPadd(data, avg);

        ComplexArray fft = Fourier.forwardDFT(data);
        double[] absSquared = fft.absSquared();
        ComplexArray ifft = Fourier.inverseDFT(new ComplexArray(absSquared));
        double[] realParts = ifft.re();

        double div = realParts[0] / length;

        int size = length;
        if (lags > 0 && lags < length - 1) {
            size = lags + 1;
        }

        double[] result = new double[size];
        System.arraycopy(realParts, 0, result, 0, size);

        return divAndRound(result, length * div);
    }

    private static double sum(double[] a) {
        double sum = 0.0;
        for (int i = 0; i < a.length; ++i) {
            sum += a[i];
        }
        return sum;
    }

    // zero-padding is necessary because FFT is actually circular
    // cross-correlation, meaning the signal wraps around at the ends.
    // So add enough zeros to get rid of the overlap, to simulate a
    // signal that is zero out to infinity.
    // https://dsp.stackexchange.com/questions/741/why-should-i-zero-pad-a-signal-before-taking-the-discrete-fourier-transform
    private static double[] centerAndPadd(double[] data, double avg) {
        // double the length; increase if necessary to ensure a full
        // cache line is used
        int newLen = data.length << 1;
        int mod8 = newLen % Double.BYTES;
        if (mod8 != 0) {
            mod8 = Double.BYTES - mod8;
            newLen += mod8;
        }
        double[] zeroPadded = new double[newLen];
        for (int i = 0; i < data.length; ++i) {
            zeroPadded[i] = data[i] - avg;
        }
        return zeroPadded;
    }

    private static double[] divAndRound(double[] data, double divisor) {
        for (int i = 0; i < data.length; ++i) {
            data[i] = Arithmetic.round(data[i] / divisor, 8);
        }
        return data;
    }

    private ACF() {
        throw new AssertionError();
    }
}
