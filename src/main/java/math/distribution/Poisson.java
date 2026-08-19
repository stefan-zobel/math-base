/*
 * Copyright 2026 Stefan Zobel
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

import math.cern.Arithmetic;
import math.cern.ProbabilityFuncs;

/**
 * The Poisson distribution of the number of events <tt>k &gt;= 0</tt> in a
 * fixed interval with PMF:
 * <p>
 * <tt>P(X = k) = (&lambda;<sup>k</sup> / k!) * e<sup>-&lambda;</sup></tt>
 * <p>
 * Valid parameter ranges: <tt>&lambda; &gt;= 0</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Poisson_distribution">Wikipedia
 * Poisson distribution</a>.
 *
 * @since 1.5.1
 */
public class Poisson implements DiscreteDistribution {

    private final double lambda;

    /**
     * Creates a Poisson distribution with the given rate.
     *
     * @param lambda
     *            the expected number of events, must be non-negative
     */
    public Poisson(double lambda) {
        if (Double.isNaN(lambda) || lambda < 0.0) {
            throw new IllegalArgumentException("lambda < 0.0 : " + lambda);
        }
        this.lambda = lambda;
    }

    @Override
    public double pmf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (lambda == 0.0) {
            return k == 0 ? 1.0 : 0.0;
        }
        return Math.exp(k * Math.log(lambda) - lambda - Arithmetic.logFactorial(k));
    }

    @Override
    public double cdf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (k == Integer.MAX_VALUE) {
            // ProbabilityFuncs.poisson() forms k + 1 internally
            return 1.0;
        }
        return ProbabilityFuncs.poisson(k, lambda);
    }

    @Override
    public double mean() {
        return lambda;
    }

    @Override
    public double variance() {
        return lambda;
    }

    @Override
    public int supportLowerBound() {
        return 0;
    }

    @Override
    public int supportUpperBound() {
        // a rate of zero degenerates to all mass on zero
        return lambda == 0.0 ? 0 : Integer.MAX_VALUE;
    }

    /**
     * @return the expected number of events
     */
    public double getLambda() {
        return lambda;
    }
}
