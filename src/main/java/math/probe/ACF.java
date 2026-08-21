/*
 * Copyright 2025, 2026 Stefan Zobel
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

import math.fft.ComplexArray;
import math.fft.Fourier;

/**
 * The sample autocorrelation function of a time series.
 * <p>
 * The estimator is the usual one, {@code r_k = gamma_k / gamma_0} with
 * {@code gamma_k = sum_t (x_t - m)(x_{t+k} - m)} over the {@code n - k} pairs
 * that exist, divided by the same {@code n} at every lag -- the biased form
 * that R's {@code acf} and statsmodels return by default, and the one whose
 * sequence is guaranteed to be positive semi-definite. {@code r_0} is exactly
 * {@code 1.0}.
 * <p>
 * A series whose values are all equal has no variation to correlate, so
 * {@code gamma_0} is zero and every {@code r_k} comes back as {@code NaN}.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Autocorrelation">Autocorrelation</a>
 */
public final class ACF {

    // Up to this many lags the plain double sum beats the transform, which
    // costs the same whether one lag is wanted or all of them. Measured, the
    // crossing moves from about 88 lags at n = 500 to about 195 at n = 10000
    // and the cost is flat either side of it, so one constant in the middle of
    // that range is worth more than a formula. The two paths agree to 5e-15.
    private static final int DIRECT_SUMMATION_MAX_LAGS = 96;

    /**
     * The autocorrelations of {@code data} at every lag from {@code 0} to
     * {@code data.length - 1}.
     * <p>
     * Note that the last lags rest on very few pairs -- the one at
     * {@code n - 1} on a single product -- and are correspondingly noisy.
     * Passing an explicit lag count is usually what a caller wants.
     *
     * @param data
     *            the time series, not {@code null}
     * @return the autocorrelations, {@code data.length} of them, starting at
     *         lag 0; an empty array for an empty series
     * @throws IllegalArgumentException
     *             if {@code data} is {@code null}
     */
    public static double[] acf(double[] data) {
        return acf(data, -1);
    }

    /**
     * The autocorrelations of {@code data} up to the given lag.
     *
     * @param data
     *            the time series, not {@code null}
     * @param lags
     *            the highest lag wanted. A value that is not positive, or one
     *            that reaches beyond the series, asks for every lag the series
     *            has rather than being refused.
     * @return the autocorrelations from lag 0 up to and including {@code lags},
     *         so {@code lags + 1} of them; an empty array for an empty series,
     *         and {@code {1.0}} for a single observation or {@code lags == 0}
     * @throws IllegalArgumentException
     *             if {@code data} is {@code null}
     */
    public static double[] acf(double[] data, int lags) {
        if (data == null) {
            throw new IllegalArgumentException("data must not be null");
        }
        final int length = data.length;
        if (length == 0) {
            return new double[0];
        }
        if (length == 1 || lags == 0) {
            return new double[] { 1.0 };
        }
        int size = length;
        if (lags > 0 && lags < length - 1) {
            size = lags + 1;
        }
        // the two paths compute the same estimator; which one runs is not
        // visible in the result
        if (size - 1 <= DIRECT_SUMMATION_MAX_LAGS) {
            return bySummation(data, size);
        }
        return byTransform(data, size);
    }

    private static double mean(double[] data) {
        double sum = 0.0;
        for (int i = 0; i < data.length; ++i) {
            sum += data[i];
        }
        return sum / data.length;
    }

    /** Few lags: centre once, then one inner loop per lag. */
    private static double[] bySummation(double[] data, int size) {
        final int length = data.length;
        final double avg = mean(data);
        double[] centered = new double[length];
        double gamma0 = 0.0;
        for (int i = 0; i < length; ++i) {
            double d = data[i] - avg;
            centered[i] = d;
            gamma0 += d * d;
        }
        double[] result = new double[size];
        for (int k = 0; k < size; ++k) {
            double sum = 0.0;
            for (int t = 0, end = length - k; t < end; ++t) {
                sum += centered[t] * centered[t + k];
            }
            result[k] = sum / gamma0;
        }
        return result;
    }

    /** Many lags: the Wiener-Khinchin route, one forward and one inverse transform. */
    private static double[] byTransform(double[] data, int size) {
        final int length = data.length;
        final double avg = mean(data);
        double[] padded = new double[paddedLength(length)];
        for (int i = 0; i < length; ++i) {
            padded[i] = data[i] - avg;
        }
        ComplexArray fft = Fourier.forwardDFT(padded);
        ComplexArray ifft = Fourier.inverseDFT(new ComplexArray(fft.absSquared()));
        double[] realParts = ifft.re();
        double gamma0 = realParts[0];
        double[] result = new double[size];
        for (int i = 0; i < size; ++i) {
            result[i] = realParts[i] / gamma0;
        }
        return result;
    }

    // The transform is a circular correlation, so the series has to be padded
    // with at least its own length in zeros for the wrap-around not to reach
    // the lags that are kept. Beyond that the length has to be a power of two,
    // or Fourier falls back on Bluestein and the transform costs several times
    // as much.
    // https://dsp.stackexchange.com/questions/741/why-should-i-zero-pad-a-signal-before-taking-the-discrete-fourier-transform
    private static int paddedLength(int length) {
        int twice = length << 1;
        return Integer.highestOneBit((twice - 1) << 1);
    }

    private ACF() {
        throw new AssertionError();
    }
}
