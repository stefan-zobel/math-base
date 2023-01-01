/*
 * Copyright 2013 Stefan Zobel
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
package math.distribution;

import math.cern.FastGamma;
import math.cern.ProbabilityFuncs;

/**
 * StudentT distribution (a.k.a "T distribution").
 */
public class StudentT implements ContinuousDistribution {

    private final double df;
    private final double pdfConst;

    public StudentT(double df) {
        if (df <= 0.0) {
            throw new IllegalArgumentException("df <= 0.0 : " + df);
        }
        double tmp = FastGamma.logGamma((df + 1.0) / 2.0) - FastGamma.logGamma(df / 2.0);
        this.pdfConst = Math.exp(tmp) / Math.sqrt(Math.PI * df);
        this.df = df;
    }

    @Override
    public double pdf(double x) {
        return pdfConst * Math.pow((1.0 + x * x / df), -(df + 1.0) * 0.5);
    }

    @Override
    public double cdf(double x) {
        return ProbabilityFuncs.studentT(df, x);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (probability >= 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return findRoot(probability, mean(), -Double.MAX_VALUE, Double.MAX_VALUE);
    }

    @Override
    public double mean() {
        if (df <= 1.0) {
            return Double.NaN;
        }
        return 0.0;
    }

    @Override
    public double variance() {
        if (df > 2.0) {
            return df / ((double) df - 2.0);
        }
        if (df == 2.0) {
            return Double.POSITIVE_INFINITY;
        }
        return Double.NaN;
    }

    public double getDegreesOfFreedom() {
        return df;
    }
}
