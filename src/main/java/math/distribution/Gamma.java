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
 * The &Gamma;(x; k, &theta;) distribution for x &gt;= 0 with PDF:
 * <p>
 * <tt>f(x; k, &theta;) = (x^(k-1) * e^(-x/&theta;)) / (&theta;^k * &Gamma;(k)) </tt>
 * where <tt>&Gamma;()</tt> is the Gamma function.
 * <p>
 * Valid parameter ranges: <tt>k &gt; 0</tt>, <tt>&theta;  &gt; 0</tt>,
 * <tt>x &gt;= 0</tt>.
 * <p>
 * Note: For a Gamma distribution to have the mean <tt>E(X)</tt> and variance
 * <tt>Var(X)</tt>, set the parameters as follows:
 * 
 * <pre>
 * k = E(X) * E(X) / Var(X)
 * &theta; = Var(X) / E(X)
 * </pre>
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Gamma_distribution">Wikipedia Gamma
 * distribution</a>.
 */
public class Gamma implements ContinuousDistribution {

    private final double shape_k;
    private final double scale_theta;
    private final double rate_beta;

    public Gamma(double shape /* k */) {
        this(shape, 1.0 /* scale */);
    }

    public Gamma(double shape /* k */, double scale /* theta */) {
        if (shape <= 0.0) {
            throw new IllegalArgumentException("shape <= 0.0");
        }
        if (scale <= 0.0) {
            throw new IllegalArgumentException("scale <= 0.0");
        }
        this.shape_k = shape;
        this.scale_theta = scale;
        this.rate_beta = (1.0 / this.scale_theta);
    }

    /**
     * Returns the probability distribution function.
     * 
     * @param x
     *            Where to compute the density function.
     * 
     * @return The value of the gamma density at x.
     */
    @Override
    public double pdf(final double x) {
        if (x < 0.0) {
            throw new IllegalArgumentException("x < 0.0");
        }
        if (x == 0.0) {
            if (shape_k == 1.0) {
                return rate_beta;
            } else if (shape_k < 1.0) {
                return Double.POSITIVE_INFINITY;
            } else {
                return 0.0;
            }
        }
        if (shape_k == 1.0) {
            return rate_beta * Math.exp(-rate_beta * x);
        }

        return rate_beta
                * Math.exp((shape_k - 1.0) * Math.log(rate_beta * x) - (rate_beta * x) - FastGamma.logGamma(shape_k));
    }

    @Override
    public double cdf(double x) {
        return ProbabilityFuncs.gamma(shape_k, rate_beta, x);
    }

    @Override
    public double mean() {
        return shape_k * scale_theta; // k * theta
    }

    @Override
    public double variance() {
        return shape_k * scale_theta * scale_theta; // k * (theta^2)
    }

    /**
     * Inverse of the Gamma cumulative distribution function.
     * 
     * @return the value X for which P(x&lt;=X).
     */
    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return 0.0;
        }
        if (probability >= 1.0) {
            return Double.MAX_VALUE;
        }
        return findRoot(probability, mean(), 0.0, Double.MAX_VALUE);
    }

    /**
     * Returns the shape parameter of this distribution.
     * 
     * @return the shape parameter.
     */
    public double getShape() {
        return shape_k;
    }

    /**
     * Returns the scale parameter of this distribution.
     * 
     * @return the scale parameter.
     */
    public double getScale() {
        return scale_theta;
    }
}
