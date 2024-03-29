/*
 * Copyright 2018 Stefan Zobel
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

/**
 * The Weibull distribution with scale parameter &lambda; &gt; {@code 0} and
 * shape parameter {@code k} &gt; {@code 0}.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Weibull_distribution">Wikipedia
 * Weibull distribution</a>.
 */
public class Weibull implements ContinuousDistribution {

    private final double scale_lambda;
    private final double shape_k;
    // helper constants
    private final double inverse_scale;
    private final double inverse_shape;
    private final double shape_dividedby_scale;
    // cached mean
    private double cached_mean = Double.NaN;

    public Weibull(double scale /* lambda */, double shape /* k */) {
        if (scale <= 0.0) {
            throw new IllegalArgumentException("scale <= 0.0");
        }
        if (shape <= 0.0) {
            throw new IllegalArgumentException("shape <= 0.0");
        }
        this.scale_lambda = scale;
        this.shape_k = shape;
        // pre-computed constants
        this.inverse_scale = 1.0 / scale;
        this.inverse_shape = 1.0 / shape;
        this.shape_dividedby_scale = inverse_scale * shape_k;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double pdf(double x) {
        if (x < 0.0) {
            return 0.0;
        }
        double xscale = x / scale_lambda;
        double xscalepow = Math.pow(xscale, shape_k - 1);
        double xscalepowshape = xscalepow * xscale;
        return shape_dividedby_scale * xscalepow * Math.exp(-xscalepowshape);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double cdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        double y = Math.pow(inverse_scale * x, shape_k);
        return -Math.expm1(-y);
    }

    /**
     * Inverse of the Weibull cumulative distribution function.
     * 
     * @param probability
     *            a given probability
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
        return scale_lambda * Math.pow(-Math.log1p(-probability), inverse_shape);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double mean() {
        if (Double.isNaN(cached_mean)) {
            cached_mean = scale_lambda * Math.exp(FastGamma.logGamma(1.0 + inverse_shape));
        }
        return cached_mean;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double variance() {
        double mean = mean();
        return ((scale_lambda * scale_lambda) * Math.exp(FastGamma.logGamma(1.0 + (2.0 * inverse_shape))))
                - (mean * mean);
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
        return scale_lambda;
    }
}
