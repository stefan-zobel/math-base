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
package math.linalg;

/**
 * Weighted least squares: minimizes {@code sum w_i (y_i - x_i'beta)^2} for
 * observations of unequal reliability, such as measurements with known
 * variances (where {@code w_i = 1 / sigma_i^2}) or aggregated data weighted by
 * group size.
 * <p>
 * The estimator is the ordinary one on the rows of {@code X} and {@code y}
 * scaled by {@code sqrt(w_i)}, and it goes through the same singular value core
 * as {@link OLS} and {@link Ridge} rather than through the normal equations.
 * That matters here more than it does for the unweighted fit: weights
 * concentrated on a few observations drive the condition number of the scaled
 * design up towards {@code cond(X) * sqrt(max w / min w)}, and at a spread of
 * {@code 1e12} an inversion of {@code X'WX} was measured to refuse a design
 * whose own condition number is 4.
 * <p>
 * It reports the same {@link LSSummary} as {@link OLS}, and the classical
 * distributions still apply, so the t values, p values and confidence
 * intervals mean what they say. {@link LSSummary#getResiduals()} returns the
 * raw {@code y - X beta}, the same quantity it returns for an unweighted fit;
 * {@link LSSummary#getSigmaHatSquared()}, {@link LSSummary#getYBar()} and
 * {@link LSSummary#getRSquared()} are weighted.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Weighted_least_squares">Weighted least
 * squares</a>
 *
 * @since 1.5.2
 */
public final class Wls {

    /**
     * Fits {@code y} on {@code X} under the given weights.
     *
     * @param alpha
     *            significance level for the confidence intervals, strictly
     *            between 0 and 1
     * @param X
     *            the design matrix, {@code n x p} with {@code n > p}
     * @param y
     *            the regressand, {@code n x 1}
     * @param weights
     *            one weight per observation, length {@code n}, each finite and
     *            strictly positive
     * @return the summary of the fit
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} is outside
     *             {@code (0, 1)}, if there are not more rows than columns, if a
     *             weight is missing, not finite or not positive, or if the
     *             scaled design is rank deficient at the customary tolerance
     */
    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y, double[] weights) {
        return estimate(alpha, X, y, weights, OLS.defaultRankTolerance(X));
    }

    /**
     * Fits {@code y} on {@code X} under the given weights, deciding for itself
     * what counts as rank deficient.
     * <p>
     * The tolerance applies to the singular values of the <em>scaled</em>
     * design {@code sqrt(W) X}, which is the matrix the fit goes through and
     * the one whose conditioning the weights drive up. See
     * {@link OLS#estimate(double, DMatrix, DMatrix, double)} for what the
     * tolerance means and what taking it into one's own hands implies.
     *
     * @param alpha
     *            significance level for the confidence intervals, strictly
     *            between 0 and 1
     * @param X
     *            the design matrix, {@code n x p} with {@code n > p}
     * @param y
     *            the regressand, {@code n x 1}
     * @param weights
     *            one weight per observation, length {@code n}, each finite and
     *            strictly positive
     * @param rankTolerance
     *            a singular value of {@code sqrt(W) X} at or below this
     *            multiple of the largest one makes the design unusable; finite
     *            and in {@code [0, 1)}
     * @return the summary of the fit
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} is outside
     *             {@code (0, 1)}, if there are not more rows than columns, if a
     *             weight is missing, not finite or not positive, if
     *             {@code rankTolerance} is outside {@code [0, 1)}, or if the
     *             scaled design is rank deficient at {@code rankTolerance}
     * @since 1.5.2
     */
    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y, double[] weights, double rankTolerance) {
        if (X == null || y == null) {
            throw new IllegalArgumentException("null argument");
        }
        if (weights == null) {
            throw new IllegalArgumentException("weights must not be null; use OLS for the unweighted fit");
        }
        if (weights.length != X.numRows()) {
            throw new IllegalArgumentException(
                    "weights.length != X.numRows : " + weights.length + " != " + X.numRows());
        }
        for (int i = 0; i < weights.length; i++) {
            double w = weights[i];
            if (Double.isNaN(w) || Double.isInfinite(w)) {
                throw new IllegalArgumentException("weights[" + i + "] is not finite : " + w);
            }
            if (w == 0.0) {
                // dropping the observation changes n and with it the degrees
                // of freedom, so it is the caller's edit to make, not ours
                throw new IllegalArgumentException("weights[" + i
                        + "] is 0, which is a request to drop that observation; remove the row instead");
            }
            if (w < 0.0) {
                throw new IllegalArgumentException("weights[" + i + "] is negative : " + w);
            }
        }
        return LeastSquaresFit.estimate(alpha, X, y, weights, rankTolerance);
    }

    private Wls() {
        throw new AssertionError();
    }
}
