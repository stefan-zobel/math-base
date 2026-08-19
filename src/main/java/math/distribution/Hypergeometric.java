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

/**
 * The hypergeometric distribution of the number of successes <tt>k</tt> in
 * <tt>n</tt> draws without replacement from a population of <tt>N</tt> items
 * containing <tt>K</tt> successes, with PMF:
 * <p>
 * <tt>P(X = k) = (C(K, k) * C(N-K, n-k)) / C(N, n)</tt>
 * <p>
 * The support is <tt>max(0, n+K-N) &lt;= k &lt;= min(n, K)</tt>. Valid
 * parameter ranges: <tt>N &gt;= 0</tt>; <tt>0 &lt;= K &lt;= N</tt>;
 * <tt>0 &lt;= n &lt;= N</tt>.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Hypergeometric_distribution">Wikipedia
 * Hypergeometric distribution</a>.
 *
 * @since 1.5.1
 */
public class Hypergeometric implements DiscreteDistribution {

    /** Relative size below which a term no longer moves the partial sum. */
    private static final double CUTOFF = 1.0e-17;

    private final int population;
    private final int successes;
    private final int draws;
    private final int lo;
    private final int hi;

    /**
     * Creates a hypergeometric distribution for {@code draws} draws without
     * replacement.
     *
     * @param population
     *            the size of the population
     * @param successes
     *            the number of items in the population that count as a success
     * @param draws
     *            the number of items drawn
     */
    public Hypergeometric(int population, int successes, int draws) {
        if (population < 0) {
            throw new IllegalArgumentException("population < 0 : " + population);
        }
        if (successes < 0 || successes > population) {
            throw new IllegalArgumentException("successes not in [0, " + population + "] : " + successes);
        }
        if (draws < 0 || draws > population) {
            throw new IllegalArgumentException("draws not in [0, " + population + "] : " + draws);
        }
        this.population = population;
        this.successes = successes;
        this.draws = draws;
        this.lo = (int) Math.max(0L, (long) draws + successes - population);
        this.hi = Math.min(draws, successes);
    }

    @Override
    public double pmf(int k) {
        if (k < lo || k > hi) {
            return 0.0;
        }
        return Math.exp(logPmf(k));
    }

    @Override
    public double cdf(int k) {
        if (k < lo) {
            return 0.0;
        }
        if (k >= hi) {
            return 1.0;
        }
        // always sum the tail that does not contain the mode: its terms decrease
        // away from k, so the sum neither cancels nor loses its leading term to
        // underflow, which summing from the far end of the support would risk
        if (k <= mean()) {
            return Math.min(1.0, lowerTail(k));
        }
        return Math.max(0.0, 1.0 - upperTail(k + 1));
    }

    @Override
    public double mean() {
        if (population == 0) {
            return 0.0;
        }
        return draws * ((double) successes / (double) population);
    }

    @Override
    public double variance() {
        if (population <= 1) {
            return 0.0;
        }
        double N = population;
        double q = successes / N;
        return draws * q * (1.0 - q) * ((N - draws) / (N - 1.0));
    }

    @Override
    public int supportLowerBound() {
        return lo;
    }

    @Override
    public int supportUpperBound() {
        return hi;
    }

    /**
     * @return the size of the population
     */
    public int getPopulation() {
        return population;
    }

    /**
     * @return the number of successes in the population
     */
    public int getSuccesses() {
        return successes;
    }

    /**
     * @return the number of items drawn
     */
    public int getDraws() {
        return draws;
    }

    /** log of C(K,k) * C(N-K,n-k) / C(N,n), for k within the support. */
    private double logPmf(int k) {
        return logBinomial(successes, k) + logBinomial(population - successes, draws - k)
                - logBinomial(population, draws);
    }

    private static double logBinomial(int a, int b) {
        return Arithmetic.logFactorial(a) - Arithmetic.logFactorial(b) - Arithmetic.logFactorial(a - b);
    }

    /**
     * Sum of the PMF over {@code [lo, k]}, stepped downwards from {@code k} with
     * the recurrence {@code t(j-1) = t(j) * j (N-K-n+j) / ((K-j+1) (n-j+1))}.
     */
    private double lowerTail(int k) {
        long N = population;
        long K = successes;
        long n = draws;
        double t = Math.exp(logPmf(k));
        double s = t;
        for (int j = k; j > lo; j--) {
            t *= (j / (double) (K - j + 1L)) * ((N - K - n + j) / (double) (n - j + 1L));
            if (t == 0.0) {
                break;
            }
            s += t;
            if (t < s * CUTOFF) {
                break;
            }
        }
        return s;
    }

    /**
     * Sum of the PMF over {@code [k, hi]}, stepped upwards from {@code k} with
     * the recurrence {@code t(j+1) = t(j) * (K-j) (n-j) / ((j+1) (N-K-n+j+1))}.
     */
    private double upperTail(int k) {
        long N = population;
        long K = successes;
        long n = draws;
        double t = Math.exp(logPmf(k));
        double s = t;
        for (int j = k; j < hi; j++) {
            t *= ((K - j) / (double) (j + 1L)) * ((n - j) / (double) (N - K - n + j + 1L));
            if (t == 0.0) {
                break;
            }
            s += t;
            if (t < s * CUTOFF) {
                break;
            }
        }
        return s;
    }
}
