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
    private final double logPdfNormFactor;

    public Beta(double alpha, double beta) {
        if (alpha <= 0.0) {
            throw new IllegalArgumentException("alpha <= 0.0");
        }
        if (beta <= 0.0) {
            throw new IllegalArgumentException("beta <= 0.0");
        }
        this.alpha = alpha;
        this.beta = beta;
        this.logPdfNormFactor = logGamma(alpha + beta) - (logGamma(alpha) + logGamma(beta));
    }

    /**
     * {@inheritDoc}
     * <p>
     * The normalizing factor is kept as its logarithm and the three factors are
     * combined in one exponential. Held as a plain number it leaves the
     * {@code double} range from {@code alpha = beta = 511} upwards, and the
     * product with the two powers then read {@code infinity} times something
     * tiny: the density of a {@code Beta(1000, 1000)} at one half is
     * {@code 35.7}, and it used to come back as {@code Infinity}, of a
     * {@code Beta(1e8, 1e8)} as {@code NaN}. Those are ordinary shapes -- an F
     * test with a few thousand degrees of freedom reaches them through
     * {@link FisherF#pdf(double)}.
     */
    @Override
    public double pdf(double x) {
        if (x < 0.0 || x > 1.0) {
            return 0.0;
        }
        if (x == 0.0) {
            // the shape alone decides the left end: a pole below one, the
            // normalizing factor at one, and zero above it. Computed,
            // (alpha - 1) * log(x) reads zero times minus infinity here
            if (alpha < 1.0) {
                return Double.POSITIVE_INFINITY;
            }
            return (alpha == 1.0) ? Math.exp(logPdfNormFactor) : 0.0;
        }
        if (x == 1.0) {
            if (beta < 1.0) {
                return Double.POSITIVE_INFINITY;
            }
            return (beta == 1.0) ? Math.exp(logPdfNormFactor) : 0.0;
        }
        return Math.exp(logPdfNormFactor + (alpha - 1.0) * Math.log(x) + (beta - 1.0) * Math.log1p(-x));
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
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        return findRoot(probability, mean(), 0.0, 1.0);
    }

    @Override
    public double supportLowerBound() {
        return 0.0;
    }

    @Override
    public double supportUpperBound() {
        return 1.0;
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
