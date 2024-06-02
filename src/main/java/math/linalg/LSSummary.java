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

import java.util.List;

import math.list.DoubleArrayList;
import math.list.DoubleList;

/**
 * Least squares regression summary.
 */
public class LSSummary {

    // significance level used
    private double alpha;

    // parameters (beta)
    private DoubleList coefficients;

    // designMatrix X
    private DMatrix designMatrix;

    // regressand y
    private DMatrix regressand;

    // average value of y
    private double yBar;

    // predicted values of y
    private DoubleList yHat;

    // y - yHat = epsilon
    private DoubleList residuals;

    // coefficient of determination
    private double rSquared;

    // population variance estimator
    private double sigmaHatSquared;

    // var-cov matrix of the coefficients
    private DMatrix varCovMatrix;

    // standard error estimates of the coefficient estimators 
    private DoubleList coefficientStandardErrors;

    // t values of the coefficients
    private DoubleList tValues;

    // p values of the coefficients
    private DoubleList pValues;

    // confidence intervals of the coefficients
    private List<DoubleList> confidenceIntervals;

    // df for t-distribution
    private int degreesOfFreedom;

    public LSSummary(double alpha, DMatrix designMatrix, DMatrix regressand) {
        this.alpha = alpha;
        this.designMatrix = designMatrix;
        this.regressand = regressand;
    }

    public void clearTemporaries() {
        designMatrix = null;
        regressand = null;
        yHat = null;
        residuals = null;
        tValues = null;
    }

    public DoubleList getBeta() {
        return coefficients;
    }

    void setBeta(DMatrix beta) {
        coefficients = new DoubleArrayList(beta.numRows());
        for (int i = 0; i < beta.numRows(); ++i) {
            coefficients.add(beta.get(i, 0));
        }
    }

    public double getYBar() {
        return yBar;
    }

    void setYBar(double yBar) {
        this.yBar = yBar;
    }

    public DoubleList getYHat() {
        return yHat;
    }

    void setYHat(DMatrix yEst) {
        yHat = new DoubleArrayList(yEst.numRows());
        for (int i = 0; i < yEst.numRows(); ++i) {
            yHat.add(yEst.get(i, 0));
        }
    }

    public DoubleList getResiduals() {
        return residuals;
    }

    void setResiduals(DMatrix epsilonHat) {
        residuals = new DoubleArrayList(epsilonHat.numRows());
        for (int i = 0; i < epsilonHat.numRows(); ++i) {
            residuals.add(epsilonHat.get(i, 0));
        }
    }

    public double getRSquared() {
        return rSquared;
    }

    void setRSquared(double rSquared) {
        this.rSquared = rSquared;
    }

    public double getSigmaHatSquared() {
        return sigmaHatSquared;
    }

    void setSigmaHatSquared(double sigmaHatSquared) {
        this.sigmaHatSquared = sigmaHatSquared;
    }

    public DMatrix getVarianceCovarianceMatrix() {
        return varCovMatrix;
    }

    void setVarianceCovarianceMatrix(DMatrix varianceCovarianceMatrix) {
        this.varCovMatrix = varianceCovarianceMatrix;
    }

    public DoubleList getCoefficientStandardErrors() {
        return coefficientStandardErrors;
    }

    void setCoefficientStandardErrors(DoubleList coefficientStandardErrors) {
        this.coefficientStandardErrors = coefficientStandardErrors;
    }

    public DoubleList getTValues() {
        return tValues;
    }

    void setTValues(DoubleList tValues) {
        this.tValues = tValues;
    }

    public DoubleList getPValues() {
        return pValues;
    }

    void setPValues(DoubleList pValues) {
        this.pValues = pValues;
    }

    public List<DoubleList> getConfidenceIntervals() {
        return confidenceIntervals;
    }

    void setConfidenceIntervals(List<DoubleList> confidenceIntervals) {
        this.confidenceIntervals = confidenceIntervals;
    }

    public int getDegreesOfFreedom() {
        return degreesOfFreedom;
    }

    void setDegreesOfFreedom(int degreesOfFreedom) {
        this.degreesOfFreedom = degreesOfFreedom;
    }

    public double getAlpha() {
        return alpha;
    }

    public DMatrix getXMatrix() {
        return designMatrix;
    }

    public DMatrix getYVector() {
        return regressand;
    }

    public int getCoefficientsCount() {
        return coefficients.size();
    }

    @Override
    public String toString() {
        return "Summary [alpha=" + alpha + ", numCoefficients=" + getCoefficientsCount() + ",\n coefficients="
                + coefficients + ",\n yBar=" + yBar + ", rSquared=" + rSquared + ", sigmaHatSquared=" + sigmaHatSquared
                + ",\n varCovMatrix=" + varCovMatrix + ", coefficientStandardErrors=" + coefficientStandardErrors
                + ",\n pValues=" + pValues + ",\n confidenceIntervals=" + confidenceIntervals + ",\n degreesOfFreedom="
                + degreesOfFreedom + "]";
    }
}
