/*
 * Copyright 2026 Stefan Zobel
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

/**
 * Provides static methods for calculating the Partial Autocorrelation Function
 * (PACF). The partial autocorrelation at lag k is the correlation that results
 * after removing the effect of any correlations due to the terms at shorter
 * lags.
 * 
 * This implementation primarily uses the Durbin-Levinson recursion to solve the
 * Yule-Walker equations efficiently.
 *
 * @since 1.4.1
 */
public final class PartialACF {

    /**
     * Computes the partial autocorrelations for a given time series up to the
     * specified number of lags. This method uses the Durbin-Levinson algorithm,
     * which recursively solves the Yule-Walker equations based on the
     * Autocorrelation Function (ACF).
     * 
     * @param data
     *            the original time series data
     * @param lags
     *            the maximum number of lags to calculate
     * @return an array of partial autocorrelations from lag 1 to {@code lags}
     * @throws IllegalArgumentException
     *             if lags is greater than or equal to data length
     */
    public static double[] partialAutocorrelation(double[] data, int lags) {
        if (data == null || data.length == 0) {
            throw new IllegalArgumentException("data[] must not be null or empty");
        }
        if (lags < 0 || lags >= data.length) {
            throw new IllegalArgumentException("lags must be >= 0 and < data.length");
        }

        double[] r = ACF.acf(data, lags);
        double[] pacf = new double[lags];

        double[] phiPrev = new double[lags + 1];
        double[] phiCurr = new double[lags + 1];
        double[] err = new double[lags + 1];
        err[0] = r[0];

        for (int k = 1; k <= lags; k++) {
            // Stability check: Avoid division by zero if error variance becomes too small
            if (Math.abs(err[k - 1]) < 1e-15) {
                // If the error is zero, the remaining PACF values are set to 0.0
                break; 
            }

            double num = r[k];
            for (int j = 1; j < k; j++) {
                num -= phiPrev[j] * r[k - j];
            }
            phiCurr[k] = num / err[k - 1];

            // that's the PACF value for lag k
            pacf[k - 1] = phiCurr[k];

            for (int j = 1; j < k; j++) {
                phiCurr[j] = phiPrev[j] - phiCurr[k] * phiPrev[k - j];
            }
            err[k] = err[k - 1] * (1.0 - phiCurr[k] * phiCurr[k]);
            System.arraycopy(phiCurr, 1, phiPrev, 1, k);
        }
        return pacf;
    }

    /**
     * Computes the 95% confidence intervals for the Partial Autocorrelation
     * Function (PACF). According to Bartlett's formula, for a white noise
     * process, the standard error is approximately 1/sqrt(n). Values outside
     * +/- 1.96 * standard error are considered statistically significant (alpha
     * = 0.05).
     * 
     * @param n
     *            the number of observations in the original time series
     * @return an array containing [lower_bound, upper_bound]
     */
    public static double[] getConfidenceInterval(int n) {
        if (n <= 0) {
            return new double[] { 0.0, 0.0 };
        }
        // 1.96 is the approximate value of the 97.5 percentile point of the
        // normal distribution
        double limit = 1.96 / Math.sqrt(n);
        return new double[] { -limit, limit };
    }

    /**
     * Checks if a given PACF value is statistically significant based on the
     * sample size.
     * 
     * @param value
     *            the PACF value at a specific lag
     * @param n
     *            the number of observations in the original time series
     * @return true if the value is outside the 95% confidence interval
     */
    public static boolean isSignificant(double value, int n) {
        double limit = 1.96 / Math.sqrt(n);
        return Math.abs(value) > limit;
    }

    private PartialACF() {
        throw new AssertionError();
    }
}
