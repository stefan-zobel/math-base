/*
 * Copyright 2017, 2026 Stefan Zobel
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

/**
 * The Cauchy distribution.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Cauchy_distribution">Wikipedia Cauchy
 * distribution</a>.
 */
public class Cauchy implements ContinuousDistribution {

    private final double loc;
    private final double scale;

    public Cauchy(double location, double scale) {
        if (scale <= 0.0) {
            throw new IllegalArgumentException("scale <= 0.0");
        }
        this.loc = location;
        this.scale = scale;
    }

    public Cauchy() {
        this(0.0, 1.0);
    }

    @Override
    public double pdf(double x) {
        double y = (x - loc) / scale;
        return 1.0 / (Math.PI * scale * (1.0 + y * y));
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code Math.PI * scale} is not formed, since it overflows for a scale
     * above about {@code 5.7e307} where the density is still an ordinary
     * number.
     */
    @Override
    public double logPdf(double x) {
        double y = (x - loc) / scale;
        return -Math.log(Math.PI) - Math.log(scale) - Math.log1p(y * y);
    }

    @Override
    public double cdf(double x) {
        double y = (x - loc) / scale;
        if (y < -0.5) {
            return Math.atan(-1.0 / y) / Math.PI;
        }
        return 0.5 + (Math.atan(y) / Math.PI);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        return loc + scale * Math.tan(Math.PI * (probability - 0.5));
    }

    @Override
    public double mean() {
        return Double.NaN;
    }

    @Override
    public double variance() {
        return Double.NaN;
    }

    /**
     * Returns the location parameter of this distribution.
     * 
     * @return the location parameter
     */
    public double getLocation() {
        return loc;
    }

    /**
     * Returns the scale parameter of this distribution.
     * 
     * @return the scale parameter
     */
    public double getScale() {
        return scale;
    }
}
