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
 * The negative binomial distribution of the number of failures
 * <tt>k &gt;= 0</tt> before the <tt>r</tt>-th success in a sequence of
 * Bernoulli trials with PMF:
 * <p>
 * <tt>P(X = k) = C(k + r - 1, k) * p<sup>r</sup> * (1 - p)<sup>k</sup></tt>
 * <p>
 * Valid parameter ranges: <tt>r &gt;= 1</tt>; <tt>0 &lt; p &lt;= 1</tt>. The
 * case <tt>r = 1</tt> is the {@link Geometric} distribution.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Negative_binomial_distribution">
 * Wikipedia Negative binomial distribution</a>.
 *
 * @since 1.5.1
 */
public class NegativeBinomial implements DiscreteDistribution {

    private final int r;
    private final double p;

    /**
     * Creates a negative binomial distribution counting the failures before the
     * {@code r}-th success.
     *
     * @param r
     *            the number of successes, must be positive
     * @param p
     *            the probability of success of a single trial
     */
    public NegativeBinomial(int r, double p) {
        if (r < 1) {
            throw new IllegalArgumentException("r < 1 : " + r);
        }
        if (!(p > 0.0) || p > 1.0) {
            throw new IllegalArgumentException("p must be in (0, 1] : " + p);
        }
        this.r = r;
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
     * The form {@link #pmf(int)} is the exponential of.
     */
    @Override
    public double logPmf(int k) {
        if (k < 0 || k > Integer.MAX_VALUE - r) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 1.0) {
            return k == 0 ? 0.0 : Double.NEGATIVE_INFINITY;
        }
        return Arithmetic.logFactorial(k + r - 1) - Arithmetic.logFactorial(k)
                - Arithmetic.logFactorial(r - 1) + r * Math.log(p) + k * Math.log1p(-p);
    }

    @Override
    public double cdf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (p == 1.0) {
            return 1.0;
        }
        if (k == Integer.MAX_VALUE) {
            // ProbabilityFuncs.negativeBinomial() forms k + 1 internally
            return 1.0;
        }
        return ProbabilityFuncs.negativeBinomial(k, r, p);
    }

    @Override
    public double mean() {
        return r * (1.0 - p) / p;
    }

    @Override
    public double variance() {
        return r * (1.0 - p) / (p * p);
    }

    @Override
    public int supportLowerBound() {
        return 0;
    }

    @Override
    public int supportUpperBound() {
        // p == 1 degenerates to all mass on zero
        return p == 1.0 ? 0 : Integer.MAX_VALUE;
    }

    /**
     * @return the number of successes
     */
    public int getNumberOfSuccesses() {
        return r;
    }

    /**
     * @return the probability of success of a single trial
     */
    public double getProbability() {
        return p;
    }
}
