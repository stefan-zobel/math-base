/*
 * Copyright 2023, 2026 Stefan Zobel
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
package math.linalg;

/**
 * Ordinary least squares regression.
 * <p>
 * Solved through the singular values rather than through the normal equations:
 * forming {@code X'X} squares the condition number, which used to let the
 * variances come out negative on a merely awkward design. The body is shared
 * with {@link Wls}, which is the same estimator on weighted rows.
 */
public final class OLS {

    /**
     * Fits {@code y} on {@code X}.
     *
     * @param alpha
     *            significance level for the confidence intervals, strictly
     *            between 0 and 1
     * @param X
     *            the design matrix, {@code n x p} with {@code n > p}
     * @param y
     *            the regressand, {@code n x 1}
     * @return the summary of the fit
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} is outside
     *             {@code (0, 1)}, or if there are not more rows than columns
     */
    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y) {
        return LeastSquaresFit.estimate(alpha, X, y, null);
    }

    private OLS() {
        // no instances
    }
}
