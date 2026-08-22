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

import math.MathConsts;

/**
 * Ordinary least squares regression.
 * <p>
 * Solved through the singular values rather than through the normal equations:
 * forming {@code X'X} squares the condition number, which used to let the
 * variances come out negative on a merely awkward design. The body is shared
 * with {@link Wls}, which is the same estimator on weighted rows.
 * <p>
 * A design whose smallest singular value is negligible against its largest one
 * is refused rather than fitted. The criterion is the customary
 * {@code max(n, p) * eps}, and it does not separate a design that is rank
 * deficient from one that is merely ill conditioned -- the second still has an
 * answer, and {@link #estimate(double, DMatrix, DMatrix, double)} is how a
 * caller who knows that reaches it.
 */
public final class OLS {

    /**
     * Fits {@code y} on {@code X}, refusing a design that is rank deficient at
     * the customary tolerance.
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
     *             {@code (0, 1)}, if there are not more rows than columns, or
     *             if {@code X} is rank deficient at the customary tolerance
     */
    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y) {
        return estimate(alpha, X, y, defaultRankTolerance(X));
    }

    /**
     * Fits {@code y} on {@code X}, deciding for itself what counts as rank
     * deficient.
     * <p>
     * {@code rankTolerance} is relative to the largest singular value of
     * {@code X} and therefore dimensionless. The three-argument overload uses
     * {@code max(n, p) * eps}, the level at which rounding alone could have
     * produced the singular value; {@code 0.0} accepts every design that is
     * not exactly singular. Between the two lie the designs that are ill
     * conditioned rather than rank deficient, which the singular value route
     * still fits -- the NIST Filip benchmark is one, and it is reached to seven
     * digits at a tolerance the default rejects it at.
     * <p>
     * A fit obtained this way is a fit like any other and comes back as a full
     * {@link LSSummary}; what makes it different is that the caller has taken
     * responsibility for the conditioning, and
     * {@link LSSummary#getConditionNumber()} is the number to judge it by.
     *
     * @param alpha
     *            significance level for the confidence intervals, strictly
     *            between 0 and 1
     * @param X
     *            the design matrix, {@code n x p} with {@code n > p}
     * @param y
     *            the regressand, {@code n x 1}
     * @param rankTolerance
     *            a singular value at or below this multiple of the largest one
     *            makes {@code X} unusable; finite and in {@code [0, 1)}
     * @return the summary of the fit
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} is outside
     *             {@code (0, 1)}, if there are not more rows than columns, if
     *             {@code rankTolerance} is outside {@code [0, 1)}, or if
     *             {@code X} is rank deficient at {@code rankTolerance}
     * @since 1.5.2
     */
    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y, double rankTolerance) {
        return LeastSquaresFit.estimate(alpha, X, y, null, rankTolerance);
    }

    /**
     * The rank tolerance {@link #estimate(double, DMatrix, DMatrix)} applies to
     * a design of this shape, {@code max(n, p) * eps}.
     *
     * @param X
     *            the design matrix
     * @return the default rank tolerance, relative to the largest singular
     *         value
     * @throws IllegalArgumentException
     *             if {@code X} is null
     * @since 1.5.2
     */
    public static double defaultRankTolerance(DMatrix X) {
        if (X == null) {
            throw new IllegalArgumentException("X is null");
        }
        return Math.max(X.numRows(), X.numColumns()) * MathConsts.MACH_EPS_DBL;
    }

    private OLS() {
        // no instances
    }
}
