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

import java.util.ArrayList;

import math.cern.ProbabilityFuncs;
import math.distribution.StudentT;
import math.list.DoubleArrayList;
import math.list.DoubleList;

/**
 * Shared body of {@link OLS} and {@link Wls}. Weighted least squares is
 * ordinary least squares on the rows of {@code X} and {@code y} scaled by
 * {@code sqrt(w_i)}, so one routine covers both: {@code weights == null} is the
 * unweighted case and skips every scaling rather than multiplying by one, which
 * is what makes the two paths agree bit for bit.
 * <p>
 * The linear algebra is one level further down in {@link SvdLeastSquares},
 * which works on the singular values rather than on the normal equations.
 */
final class LeastSquaresFit {

    /**
     * Fits {@code y} on {@code X}, optionally weighted, and fills an
     * {@link LSSummary}.
     *
     * @param alpha
     *            significance level for the confidence intervals, strictly
     *            between 0 and 1
     * @param X
     *            the design matrix, {@code n x p} with {@code n > p}
     * @param y
     *            the regressand, {@code n x 1}
     * @param weights
     *            the weights, length {@code n} and all strictly positive, or
     *            {@code null} for the unweighted fit. The caller has already
     *            validated them.
     * @return the summary of the fit
     */
    static LSSummary estimate(double alpha, DMatrix X, DMatrix y, double[] weights) {
        if (X.numRows() != y.numRows()) {
            throw new IllegalArgumentException("X.numRows != y.numRows : " + X.numRows() + " != " + y.numRows());
        }
        if (X.numRows() - X.numColumns() < 1) {
            throw new IllegalArgumentException("degrees of freedom < 1 : " + (X.numRows() - X.numColumns()));
        }
        if (alpha <= 0.0) {
            throw new IllegalArgumentException("alpha <= 0 : " + alpha);
        }
        if (alpha >= 1.0) {
            throw new IllegalArgumentException("alpha >= 1 : " + alpha);
        }
        LSSummary smmry = new LSSummary(alpha, X, y);
        int n = X.numRows();
        int p = X.numColumns();

        // the square roots of the weights, or null, in which case every
        // scaling below hands back the matrix it was given
        double[] sqrtW = null;
        if (weights != null) {
            sqrtW = new double[n];
            for (int i = 0; i < n; ++i) {
                sqrtW[i] = Math.sqrt(weights[i]);
            }
            smmry.setWeights(weights);
        }
        DMatrix Xw = scaleRows(X, sqrtW);
        DMatrix yw = scaleRows(y, sqrtW);

        // Solved through the singular values rather than through the normal
        // equations: forming X'WX squares the condition number, which used to
        // let the variances come out negative on a merely awkward design.
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(Xw.getArrayUnsafe(), n, p);
        if (!svd.converged) {
            throw new RuntimeException("the singular value decomposition of X did not converge");
        }
        int deficient = SvdLeastSquares.rankDeficientAt(svd);
        if (deficient >= 0) {
            throw new RuntimeException(
                    "X is rank deficient: singular value " + deficient + " of " + p + " is negligible");
        }
        DMatrix beta = new DMatrix(p, 1, SvdLeastSquares.solve(svd, yw.getArrayUnsafe(), 0.0));
        smmry.setBeta(beta);
        // the fitted values live in the scale of the data, not of the weights
        DMatrix yHat = X.mul(beta);
        smmry.setYHat(yHat);

        // ybar = sum w_i y_i / sum w_i, written as onesW' yw over onesW' onesW;
        // with no weights the denominator is exactly n and this is the mean
        DMatrix onesRow = new DMatrix(1, n);
        for (int i = 0; i < n; ++i) {
            onesRow.setUnsafe(0, i, sqrtW == null ? 1.0 : sqrtW[i]);
        }
        double sumW = n;
        if (sqrtW != null) {
            DMatrix onesCol = new DMatrix(n, 1);
            for (int i = 0; i < n; ++i) {
                onesCol.setUnsafe(i, 0, sqrtW[i]);
            }
            sumW = onesRow.mul(onesCol).get(0, 0);
        }
        double ybar = onesRow.mul(yw).scaleInplace(1.0 / sumW).get(0, 0);
        smmry.setYBar(ybar);

        DMatrix ones = new DMatrix(n, 1);
        for (int i = 0; i < n; ++i) {
            ones.setUnsafe(i, 0, 1.0);
        }
        DMatrix yBarMat = ones.scaleInplace(ybar);
        DMatrix a = scaleRows(yHat.minus(yBarMat), sqrtW);
        DMatrix b = scaleRows(y.minus(yBarMat), sqrtW);
        double SQE = a.transpose().mul(a).get(0, 0);
        double SQT = b.transpose().mul(b).get(0, 0);
        double R_squared = SQE / SQT;
        smmry.setRSquared(R_squared > 1.0 ? 1.0 : R_squared);

        // the residuals reported are the raw ones, y - X beta, so that they
        // mean the same thing as in an unweighted fit; the weighted ones are
        // sqrt(w_i) times these and the weights come back with the summary
        DMatrix epsHat = y.minus(yHat);
        smmry.setResiduals(epsHat);
        int df = n - p;
        smmry.setDegreesOfFreedom(df);
        DMatrix epsW = scaleRows(epsHat, sqrtW);
        double sigmaHatSquared = epsW.transpose().mul(epsW).scaleInplace(1.0 / (df)).get(0, 0);
        smmry.setSigmaHatSquared(sigmaHatSquared);

        // Xw' Xw is X'WX, so the variance matrix needs no weighting of its own
        DMatrix varCov = new DMatrix(p, p, SvdLeastSquares.varianceMatrix(svd, 0.0)).scaleInplace(sigmaHatSquared);
        smmry.setVarianceCovarianceMatrix(varCov);
        DoubleList standardErrors = new DoubleArrayList(varCov.numRows());
        for (int i = 0; i < varCov.numRows(); ++i) {
            // V diag(1/d_i^2) V' has a diagonal of sums of squares, so no
            // guard against a negative variance is needed here any more
            standardErrors.add(Math.sqrt(varCov.get(i, i)));
        }
        smmry.setCoefficientStandardErrors(standardErrors);
        DoubleList tValues = new DoubleArrayList(varCov.numRows());
        DoubleList pValues = new DoubleArrayList(varCov.numRows());
        ArrayList<DoubleList> confIntervals = new ArrayList<>();
        StudentT tDist = new StudentT(df);
        double tval = tDist.inverseCdf(1.0 - (alpha / 2.0));
        for (int i = 0; i < varCov.numRows(); ++i) {
            double coeff = beta.get(i, 0);
            double se = standardErrors.get(i);
            double t = coeff / se;
            // P(|T| > |t|) as a regularized incomplete beta. Written as
            // 2 * (1 - cdf(|t|)) this cancels to exactly 0.0 from about
            // |t| = 8 upwards, and a p value of exactly zero is always wrong.
            double pv = ProbabilityFuncs.beta(df / 2.0, 0.5, df / (df + t * t));
            double min = coeff - tval * se;
            double max = coeff + tval * se;
            tValues.add(t);
            pValues.add(pv);
            confIntervals.add(DoubleList.of(min, max));
        }
        smmry.setTValues(tValues);
        smmry.setPValues(pValues);
        smmry.setConfidenceIntervals(confIntervals);
        return smmry;
    }

    /**
     * Row scaling by {@code sqrtW} into a fresh matrix, or the matrix itself
     * when there are no weights. The caller's {@code X} and {@code y} are never
     * modified.
     */
    private static DMatrix scaleRows(DMatrix M, double[] sqrtW) {
        if (sqrtW == null) {
            return M;
        }
        int rows = M.numRows();
        int cols = M.numColumns();
        DMatrix out = new DMatrix(rows, cols);
        double[] src = M.getArrayUnsafe();
        double[] dst = out.getArrayUnsafe();
        for (int j = 0; j < cols; j++) {
            int off = j * rows;
            for (int i = 0; i < rows; i++) {
                dst[off + i] = sqrtW[i] * src[off + i];
            }
        }
        return out;
    }

    private LeastSquaresFit() {
        throw new AssertionError();
    }
}
