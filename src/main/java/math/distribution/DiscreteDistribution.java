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

import math.cern.ProbabilityFuncs;

/**
 * A discrete distribution whose support is a contiguous range of integers.
 * <p>
 * This is the discrete counterpart of {@link ContinuousDistribution}. It is a
 * separate hierarchy because the two differ in more than a name: the mass
 * function is a probability, not a density, and the CDF is a step function, so
 * the Newton iteration of {@link ContinuousDistribution#findRoot} has nothing
 * to work on.
 *
 * @since 1.5.1
 */
public interface DiscreteDistribution {

    /**
     * Returns the probability mass function (PMF) of this distribution
     * evaluated at {@code k}, i.e. {@code P(X == k)}. Outside of the support
     * the result is {@code 0.0}.
     *
     * @param k
     *            the point at which the PMF is evaluated
     * @return the probability that a random variable with this distribution
     *         takes the value {@code k}
     */
    double pmf(int k);

    /**
     * For a random variable {@code X} whose values are distributed according to
     * this distribution, this method returns {@code P(X <= k)}, i.e. the
     * (cumulative) distribution function (CDF) of this distribution.
     *
     * @param k
     *            the point at which the CDF is evaluated
     * @return the probability that a random variable with this distribution
     *         takes a value less than or equal to {@code k}
     */
    double cdf(int k);

    /**
     * Use this method to get the mean of this distribution.
     *
     * @return the mean or {@code Double.NaN} if it is not defined
     */
    double mean();

    /**
     * Use this method to get the variance of this distribution.
     *
     * @return the variance or {@code Double.NaN} if it is not defined
     */
    double variance();

    /**
     * Returns the smallest value the support of this distribution contains.
     *
     * @return the lower end of the support (inclusive)
     */
    int supportLowerBound();

    /**
     * Returns the largest value the support of this distribution contains, or
     * {@code Integer.MAX_VALUE} if the support is unbounded above.
     *
     * @return the upper end of the support (inclusive)
     */
    int supportUpperBound();

    /**
     * For a random variable {@code X} whose values are distributed according to
     * this distribution, this method returns {@code P(k0 < X <= k1)}.
     *
     * @param k0
     *            lower bound (excluded)
     * @param k1
     *            upper bound (included)
     * @return the probability that a random variable with this distribution
     *         takes a value between {@code k0} and {@code k1}, excluding the
     *         lower and including the upper endpoint
     * @throws IllegalArgumentException
     *             if {@code k0 > k1}
     */
    default double probability(int k0, int k1) {
        if (k0 > k1) {
            throw new IllegalArgumentException(
                    "Lower endpoint (" + k0 + ") must be less than or equal to upper endpoint (" + k1 + ")");
        }
        return cdf(k1) - cdf(k0);
    }

    /**
     * Returns the smallest {@code k} in the support for which
     * {@code cdf(k) >= p}, i.e. the {@code p}-quantile of this distribution.
     * {@code p == 0.0} returns {@link #supportLowerBound()}, {@code p == 1.0}
     * returns {@link #supportUpperBound()}.
     * <p>
     * The default implementation starts from a first order Cornish-Fisher
     * estimate, brackets the answer by doubling steps outwards and then
     * bisects. Subclasses with a closed form should override it. Note that a
     * CDF evaluated in floating point is monotone only up to rounding, so what
     * is guaranteed is the bracket condition on the <em>computed</em> CDF.
     *
     * @param p
     *            the probability for which the quantile is evaluated, in
     *            {@code [0, 1]}
     * @return the smallest {@code k} with {@code cdf(k) >= p}
     * @throws IllegalArgumentException
     *             if {@code p} is {@code NaN} or lies outside of {@code [0, 1]}
     */
    default int quantile(double p) {
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + p);
        }
        long lo = supportLowerBound();
        long hi = supportUpperBound();
        if (p == 0.0) {
            return (int) lo;
        }
        if (p == 1.0) {
            return (int) hi;
        }
        // first order Cornish-Fisher starting point, clamped into the support
        long x;
        double guess = Math.floor(mean() + ProbabilityFuncs.normalInverse(p) * Math.sqrt(variance()));
        if (Double.isNaN(guess) || guess <= lo) {
            x = lo;
        } else if (guess >= hi) {
            x = hi;
        } else {
            x = (long) guess;
        }
        // bracket: cdf(a) < p <= cdf(b), with a == lo - 1 standing for "no such
        // a exists". All arithmetic is done in long so that the doubling steps
        // cannot wrap around; cdf is only called once the index is known to lie
        // within the support.
        long a;
        long b;
        if (cdf((int) x) >= p) {
            b = x;
            long step = 1L;
            a = x - 1L;
            while (a >= lo && cdf((int) a) >= p) {
                b = a;
                step += step;
                a = b - step;
            }
            if (a < lo) {
                a = lo - 1L;
            }
        } else {
            a = x;
            long step = 1L;
            b = x + 1L;
            while (b <= hi && cdf((int) b) < p) {
                a = b;
                step += step;
                b = a + step;
            }
            if (b > hi) {
                b = hi;
            }
        }
        while (b - a > 1L) {
            long m = a + (b - a) / 2L;
            if (cdf((int) m) >= p) {
                b = m;
            } else {
                a = m;
            }
        }
        return (int) b;
    }
}
