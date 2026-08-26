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
 * The binomial distribution of the number of successes <tt>k</tt> in <tt>n</tt>
 * independent Bernoulli trials with PMF:
 * <p>
 * <tt>P(X = k) = C(n, k) * p<sup>k</sup> * (1 - p)<sup>n-k</sup></tt>
 * <p>
 * Valid parameter ranges: <tt>n &gt;= 0</tt>; <tt>0 &lt;= p &lt;= 1</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Binomial_distribution">Wikipedia
 * Binomial distribution</a>.
 *
 * @since 1.5.1
 */
public class Binomial implements DiscreteDistribution {

    private final int n;
    private final double p;

    /**
     * Creates a binomial distribution for {@code n} trials with success
     * probability {@code p}.
     *
     * @param n
     *            the number of trials
     * @param p
     *            the probability of success of a single trial
     */
    public Binomial(int n, double p) {
        if (n < 0) {
            throw new IllegalArgumentException("n < 0 : " + n);
        }
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + p);
        }
        this.n = n;
        this.p = p;
    }

    @Override
    public double pmf(int k) {
        // the mass was already formed as the exponential of the logarithm, so
        // this returns the bits it always did
        return Math.exp(logPmf(k));
    }

    /**
     * {@inheritDoc}
     * <p>
     * The form {@link #pmf(int)} is the exponential of. It answers where the
     * mass itself cannot: {@code Binomial(2000, 0.5).pmf(0)} underflows to
     * zero and the logarithm is {@code -1386.29}.
     */
    @Override
    public double logPmf(int k) {
        if (k < 0 || k > n) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 0.0) {
            return k == 0 ? 0.0 : Double.NEGATIVE_INFINITY;
        }
        if (p == 1.0) {
            return k == n ? 0.0 : Double.NEGATIVE_INFINITY;
        }
        return Arithmetic.logFactorial(n) - Arithmetic.logFactorial(k) - Arithmetic.logFactorial(n - k)
                + k * Math.log(p) + (n - k) * Math.log1p(-p);
    }

    @Override
    public double cdf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (k >= n) {
            return 1.0;
        }
        // the incomplete beta integral that ProbabilityFuncs.binomial() uses is
        // undefined at these two endpoints
        if (p == 0.0) {
            return 1.0;
        }
        if (p == 1.0) {
            return 0.0;
        }
        return ProbabilityFuncs.binomial(k, n, p);
    }

    @Override
    public double mean() {
        return n * p;
    }

    @Override
    public double variance() {
        return n * p * (1.0 - p);
    }

    @Override
    public int supportLowerBound() {
        // p == 1 degenerates to all mass on n
        return p == 1.0 ? n : 0;
    }

    @Override
    public int supportUpperBound() {
        // p == 0 degenerates to all mass on zero
        return p == 0.0 ? 0 : n;
    }

    /**
     * @return the number of trials
     */
    public int getNumberOfTrials() {
        return n;
    }

    /**
     * @return the probability of success of a single trial
     */
    public double getProbability() {
        return p;
    }
}
