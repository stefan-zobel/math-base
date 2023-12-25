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

import math.cern.ProbabilityFuncs;

import static math.cern.FastGamma.logGamma;

/**
 * The Beta distribution.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Beta_distribution">Wikipedia Beta
 * distribution</a>.
 */
public class Beta implements ContinuousDistribution {

    private final double alpha;
    private final double beta;
    private final double pdfNormFactor;

    public Beta(double alpha, double beta) {
        if (alpha <= 0.0) {
            throw new IllegalArgumentException("alpha <= 0.0");
        }
        if (beta <= 0.0) {
            throw new IllegalArgumentException("beta <= 0.0");
        }
        this.alpha = alpha;
        this.beta = beta;
        this.pdfNormFactor = Math.exp(logGamma(alpha + beta) - (logGamma(alpha) + logGamma(beta)));
    }

    @Override
    public double pdf(double x) {
        if (x < 0.0 || x > 1.0) {
            return 0.0;
        }
        return pdfNormFactor * Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1);
    }

    @Override
    public double cdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        } else if (x >= 1.0) {
            return 1.0;
        }
        return ProbabilityFuncs.beta(alpha, beta, x);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return 0.0;
        }
        if (probability >= 1.0) {
            return 1.0;
        }
        return findRoot(probability, mean(), 0.0, 1.0);
    }

    @Override
    public double mean() {
        return alpha / (alpha + beta);
    }

    @Override
    public double variance() {
        double alphaPlusBeta = alpha + beta;
        return (alpha * beta) / (alphaPlusBeta * alphaPlusBeta * (alphaPlusBeta + 1.0));
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }
}
