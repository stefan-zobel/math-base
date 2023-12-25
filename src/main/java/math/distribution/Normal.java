/*
 * Copyright 2013, 2023 Stefan Zobel
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

import math.cern.ProbabilityFuncs;
import static math.MathConsts.SQRT_TWO_PI;

/**
 * Normal (a.k.a Gaussian) distribution.
 */
public class Normal implements ContinuousDistribution {

    /** Mean of this distribution */
    private final double mean;

    /** Standard deviation of this distribution */
    private final double stdDev;

    /** Variance of this distribution */
    private final double variance;

    /** 1.0 / (stdDev * sqrt(2 * PI)) */
    private final double factor;

    public Normal() {
        this(0.0, 1.0);
    }

    public Normal(double mean, double stdDev) {
        if (stdDev <= 0.0) {
            throw new IllegalArgumentException("Standard deviation must be positive (" + stdDev + ")");
        }
        this.mean = mean;
        this.stdDev = stdDev;
        this.variance = stdDev * stdDev;
        this.factor = (1.0 / (this.variance * SQRT_TWO_PI));
    }

    @Override
    public double pdf(double x) {
        double xMinusMu = (x - mean);
        return factor * Math.exp(-(xMinusMu * xMinusMu) / (2.0 * variance));
    }

    @Override
    public double cdf(double x) {
        return ProbabilityFuncs.normal(mean, variance, x);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (probability >= 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return mean + stdDev * ProbabilityFuncs.normalInverse(probability);
    }

    @Override
    public double mean() {
        return mean;
    }

    @Override
    public double variance() {
        return variance;
    }
}
