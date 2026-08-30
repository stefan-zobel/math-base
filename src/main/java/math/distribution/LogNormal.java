/*
 * Copyright 2015, 2026 Stefan Zobel
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

/**
 * The log-normal distribution with parameters &mu; and &sigma; &gt; {@code 0}.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Log-normal_distribution">Wikipedia
 * Log-normal distribution</a>.
 */
public class LogNormal implements ContinuousDistribution {

    private final double mu;
    private final double sigma;

    public LogNormal(double mu, double sigma) {
        if (sigma <= 0.0) {
            throw new IllegalArgumentException("sigma <= 0.0 : " + sigma);
        }
        this.mu = mu;
        this.sigma = sigma;
    }

    @Override
    public double pdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        double d = Math.log(x) - mu;
        return Math.exp((-d * d) / (2.0 * (sigma * sigma))) / (Math.sqrt(2.0 * Math.PI) * sigma * x);
    }

    /**
     * {@inheritDoc}
     * <p>
     * The density underflows to zero from about {@code x = 2.1e16} for the
     * standard log-normal, which is well inside the range of the variable
     * itself; this goes on answering there.
     */
    @Override
    public double logPdf(double x) {
        if (x <= 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        double logX = Math.log(x);
        double z = (logX - mu) / sigma;
        return -0.5 * (z * z) - Math.log(sigma) - logX - Math.log(Math.sqrt(2.0 * Math.PI));
    }

    @Override
    public double cdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        return ProbabilityFuncs.normal((Math.log(x) - mu) / sigma);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        return Math.exp(mu + sigma * ProbabilityFuncs.normalInverse(probability));
    }

    @Override
    public double supportLowerBound() {
        return 0.0;
    }
    @Override
    public double mean() {
        return (Math.exp(mu + (sigma * sigma) / 2.0));
    }

    @Override
    public double variance() {
        double sigsig = sigma * sigma;
        return ((Math.exp(2.0 * mu + sigsig) * (Math.exp(sigsig) - 1.0)));
    }
}
