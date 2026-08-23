/*
 * Copyright 2013, 2026 Stefan Zobel
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
 * 
 * <pre>
 *                 1                       2
 *    pdf(x) = -----------   exp( - (x-mean) / 2v ) 
 *             sqrt(2pi*v)
 * 
 *                            x
 *                            -
 *                 1         | |                 2
 *    cdf(x) = -----------   |    exp( - (t-mean) / 2v ) dt
 *             sqrt(2pi*v) | |
 *                          -
 *                         -inf.
 * </pre>
 * 
 * where <tt>v = variance = standardDeviation^2</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Normal_distribution">Wikipedia Normal
 * distribution</a>.
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
        this.factor = (1.0 / (this.stdDev * SQRT_TWO_PI));
    }

    /**
     * {@inheritDoc}
     * <p>
     * The exponent is formed from {@code (x - mean) / stdDev} rather than from
     * the squared distance over the variance. {@code variance} is
     * {@code stdDev * stdDev} and leaves the {@code double} range on both
     * sides -- it underflows to zero below a standard deviation of about
     * {@code 1.5e-162} and overflows above {@code 1.3e154} -- and the quotient
     * then read {@code 0/0} or {@code infinity/infinity}. A
     * {@code Normal(0, 1e-170)} had density {@code NaN} at its own mean, where
     * the value is {@code 3.99e169}.
     */
    @Override
    public double pdf(double x) {
        double z = (x - mean) / stdDev;
        return factor * Math.exp(-0.5 * (z * z));
    }

    @Override
    public double cdf(double x) {
        return ProbabilityFuncs.normal(mean, variance, x);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
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
