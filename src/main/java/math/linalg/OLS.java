/*
 * Copyright 2023 Stefan Zobel
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
        DMatrix Xtrans = X.transpose();
        // Note: this may be numerically unstable!
        DMatrix XtransTimesXInverse = Xtrans.mul(X).inverse();
        DMatrix beta = XtransTimesXInverse.mul(Xtrans).mul(y);
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
        DMatrix varCov = XtransTimesXInverse.scaleInplace(sigmaHatSquared);
        smmry.setVarianceCovarianceMatrix(varCov);
        DoubleList standardErrors = new DoubleArrayList(varCov.numRows());
        for (int i = 0; i < varCov.numRows(); ++i) {
            double vari = varCov.get(i, i);
            if (vari < 0.0) {
                vari = Double.MIN_NORMAL;
                varCov.set(i, i, vari);
            }
            standardErrors.add(Math.sqrt(vari));
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
            double p = 2.0 * (1.0 - tDist.cdf(Math.abs(t)));
            double min = coeff - tval * se;
            double max = coeff + tval * se;
            tValues.add(t);
            pValues.add(p);
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
