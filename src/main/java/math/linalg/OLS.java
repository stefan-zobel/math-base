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
 * Poor man's naive ordinary least squares regression.
 */
public final class OLS {

    public static LSSummary estimate(double alpha, DMatrix X, DMatrix y) {
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
        // Solved through the singular values rather than through the normal
        // equations: forming X'X squares the condition number, which used to
        // let the variances come out negative on a merely awkward design.
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(X.getArrayUnsafe(), n, p);
        if (!svd.converged) {
            throw new RuntimeException("the singular value decomposition of X did not converge");
        }
        int deficient = SvdLeastSquares.rankDeficientAt(svd);
        if (deficient >= 0) {
            throw new RuntimeException(
                    "X is rank deficient: singular value " + deficient + " of " + p + " is negligible");
        }
        DMatrix beta = new DMatrix(p, 1, SvdLeastSquares.solve(svd, y.getArrayUnsafe(), 0.0));
        smmry.setBeta(beta);
        DMatrix yHat = X.mul(beta);
        smmry.setYHat(yHat);
        DMatrix ones = new DMatrix(1, y.numRows());
        for (int i = 0; i < y.numRows(); ++i) {
            ones.setUnsafe(0, i, 1.0);
        }
        double ybar = ones.mul(y).scaleInplace(1.0 / y.numRows()).get(0, 0);
        smmry.setYBar(ybar);
        ones = new DMatrix(y.numRows(), 1);
        for (int i = 0; i < y.numRows(); ++i) {
            ones.setUnsafe(i, 0, 1.0);
        }
        DMatrix yBarMat = ones.scaleInplace(ybar);
        DMatrix a = yHat.minus(yBarMat);
        DMatrix b = y.minus(yBarMat);
        double SQE = a.transpose().mul(a).get(0, 0);
        double SQT = b.transpose().mul(b).get(0, 0);
        double R_squared = SQE / SQT;
        smmry.setRSquared(R_squared > 1.0 ? 1.0 : R_squared);
        DMatrix epsHat = y.minus(yHat);
        smmry.setResiduals(epsHat);
        int df = epsHat.numRows() - X.numColumns();
        smmry.setDegreesOfFreedom(df);
        double sigmaHatSquared = epsHat.transpose().mul(epsHat).scaleInplace(1.0 / (df)).get(0, 0);
        smmry.setSigmaHatSquared(sigmaHatSquared);
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

    private OLS() {
        // no instances
    }
}
