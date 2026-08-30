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
 * The discrete uniform distribution on the integers <tt>a, a+1, ..., b</tt>
 * with PMF:
 * <p>
 * <tt>P(X = k) = 1 / (b - a + 1)</tt>
 * <p>
 * Valid parameter ranges: <tt>a &lt;= b</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Discrete_uniform_distribution">
 * Wikipedia Discrete uniform distribution</a>.
 *
 * @since 1.5.1
 */
public class DiscreteUniform implements DiscreteDistribution {

    private final int a;
    private final int b;
    private final long count;

    /**
     * Creates a discrete uniform distribution on <tt>{a, ..., b}</tt>.
     *
     * @param a
     *            the lower end of the support (inclusive)
     * @param b
     *            the upper end of the support (inclusive)
     */
    public DiscreteUniform(int a, int b) {
        if (b < a) {
            throw new IllegalArgumentException("b < a : " + b + " < " + a);
        }
        this.a = a;
        this.b = b;
        this.count = (long) b - (long) a + 1L;
    }

    @Override
    public double pmf(int k) {
        if (k < a || k > b) {
            return 0.0;
        }
        return 1.0 / count;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double logPmf(int k) {
        if (k < a || k > b) {
            return Double.NEGATIVE_INFINITY;
        }
        return -Math.log(count);
    }

    @Override
    public double cdf(int k) {
        if (k < a) {
            return 0.0;
        }
        if (k >= b) {
            return 1.0;
        }
        return ((double) ((long) k - (long) a + 1L)) / count;
    }

    @Override
    public int quantile(double p) {
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + p);
        }
        if (p == 0.0) {
            return a;
        }
        if (p == 1.0) {
            return b;
        }
        long k = (long) Math.ceil(p * count) - 1L;
        if (k < 0L) {
            k = 0L;
        }
        if (k > count - 1L) {
            k = count - 1L;
        }
        // guard against the rounding of p * count
        if (k > 0L && cdf((int) (a + k - 1L)) >= p) {
            --k;
        } else if (k < count - 1L && cdf((int) (a + k)) < p) {
            ++k;
        }
        return (int) (a + k);
    }

    @Override
    public double mean() {
        return ((double) a + (double) b) / 2.0;
    }

    @Override
    public double variance() {
        return ((double) count * (double) count - 1.0) / 12.0;
    }

    @Override
    public int supportLowerBound() {
        return a;
    }

    @Override
    public int supportUpperBound() {
        return b;
    }

    /**
     * @return the lower end of the support
     */
    public int getA() {
        return a;
    }

    /**
     * @return the upper end of the support
     */
    public int getB() {
        return b;
    }
}
