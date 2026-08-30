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

/**
 * The geometric distribution of the number of failures <tt>k &gt;= 0</tt>
 * before the first success in a sequence of Bernoulli trials with PMF:
 * <p>
 * <tt>P(X = k) = (1 - p)<sup>k</sup> * p</tt>
 * <p>
 * This is the {@code r = 1} case of {@link NegativeBinomial}. Valid parameter
 * ranges: <tt>0 &lt; p &lt;= 1</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Geometric_distribution">Wikipedia
 * Geometric distribution</a>.
 *
 * @since 1.5.1
 */
public class Geometric implements DiscreteDistribution {

    private final double p;
    private final double log1mp;

    /**
     * Creates a geometric distribution with success probability {@code p}.
     *
     * @param p
     *            the probability of success of a single trial
     */
    public Geometric(double p) {
        if (!(p > 0.0) || p > 1.0) {
            throw new IllegalArgumentException("p must be in (0, 1] : " + p);
        }
        this.p = p;
        this.log1mp = Math.log1p(-p);
    }

    @Override
    public double pmf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (p == 1.0) {
            return k == 0 ? 1.0 : 0.0;
        }
        return p * Math.exp(k * log1mp);
    }

    /**
     * {@inheritDoc}
     * <p>
     * The mass of a {@code Geometric(0.5)} underflows to zero from
     * {@code k = 1074} upwards, where this is still answering.
     */
    @Override
    public double logPmf(int k) {
        if (k < 0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 1.0) {
            return k == 0 ? 0.0 : Double.NEGATIVE_INFINITY;
        }
        return Math.log(p) + k * log1mp;
    }

    @Override
    public double cdf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (p == 1.0) {
            return 1.0;
        }
        return -Math.expm1((k + 1.0) * log1mp);
    }

    @Override
    public int quantile(double prob) {
        if (Double.isNaN(prob) || prob < 0.0 || prob > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + prob);
        }
        if (prob == 0.0 || p == 1.0) {
            return 0;
        }
        if (prob == 1.0) {
            return Integer.MAX_VALUE;
        }
        double x = Math.ceil(Math.log1p(-prob) / log1mp) - 1.0;
        long k;
        if (Double.isNaN(x) || x <= 0.0) {
            k = 0L;
        } else if (x >= Integer.MAX_VALUE) {
            k = Integer.MAX_VALUE;
        } else {
            k = (long) x;
        }
        // guard against the rounding of the two logarithms
        if (k > 0L && cdf((int) (k - 1L)) >= prob) {
            --k;
        } else if (k < Integer.MAX_VALUE && cdf((int) k) < prob) {
            ++k;
        }
        return (int) k;
    }

    @Override
    public double mean() {
        return (1.0 - p) / p;
    }

    @Override
    public double variance() {
        return (1.0 - p) / (p * p);
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
     * @return the probability of success of a single trial
     */
    public double getProbability() {
        return p;
    }
}
