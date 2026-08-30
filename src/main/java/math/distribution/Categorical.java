package math.distribution;

import math.rng.AliasTable;
import math.rng.PseudoRandom;

/**
 * The categorical distribution over the categories {@code 0 .. m-1}: a single
 * draw landing in a single category, with mass {@code P(X = k) = p_k}.
 * <p>
 * It is the one draw case of {@link Multinomial} -- {@code Multinomial(1, p)}
 * is this law written as a vector of indicators -- and the likelihood
 * {@link Dirichlet} is the conjugate prior of when observations arrive one at a
 * time rather than already counted.
 * <p>
 * Unlike those two this is an ordinary {@link DiscreteDistribution}. An outcome
 * is a single integer, so there is a distribution function, a quantile, a mean
 * and a variance, and all four are here. What is unusual is only that its
 * masses are given rather than computed from a formula: it is the one law in
 * this package whose shape is an argument.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom)} cannot -- that restriction belongs to
 * the generator.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Categorical_distribution">Wikipedia
 * Categorical distribution</a>.
 *
 * @since 1.5.3
 */
public final class Categorical implements DiscreteDistribution {

    /** The weights as given, which need not sum to one. */
    private final double[] p;
    /** Their sum, the divisor that turns a weight into a probability. */
    private final double p0;
    /** {@code log(p_i / p0)}, which is {@code -Infinity} at a weight of zero. */
    private final double[] logP;
    /** The distribution function, {@code (p_0 + ... + p_i) / p0}. */
    private final double[] cum;
    private final double mean;
    private final double variance;
    /** The first and the last category of positive weight. */
    private final int lo;
    private final int hi;
    private final AliasTable table;

    /**
     * A categorical distribution over the given categories.
     * <p>
     * The weights need not sum to one; they are normalized, so counts serve as
     * well as probabilities. A weight of zero is allowed and its category is
     * then impossible rather than merely unlikely.
     *
     * @param weights
     *            the weight of each category, at least one of them, each finite
     *            and not negative, and not all zero. Not modified
     * @throws IllegalArgumentException
     *             if {@code weights} is {@code null}, is empty, holds a value
     *             that is negative or not finite, or does not sum to a finite
     *             positive number
     */
    public Categorical(double[] weights) {
        if (weights == null) {
            throw new IllegalArgumentException("weights must not be null");
        }
        int m = weights.length;
        if (m == 0) {
            throw new IllegalArgumentException("weights must not be empty");
        }
        double sum = 0.0;
        for (int i = 0; i < m; i++) {
            double w = weights[i];
            if (!(w >= 0.0) || w == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException("weights[" + i + "] is not a finite non-negative number : " + w);
            }
            sum += w;
        }
        if (!(sum > 0.0) || sum == Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException("the weights must sum to a finite positive number, not " + sum);
        }

        this.p = weights.clone();
        this.p0 = sum;
        this.logP = new double[m];
        double logP0 = Math.log(sum);
        for (int i = 0; i < m; i++) {
            double q = p[i] / sum;
            // the quotient is the accurate form wherever it is a normal number,
            // but weights spanning hundreds of orders of magnitude -- what
            // importance weights do just before a resampling step -- underflow
            // it to zero and lose an answer that is merely small. There the
            // difference of the logarithms is exact enough: it is only reached
            // below exp(-708), where neither logarithm exceeds 709 in magnitude
            // and there is no cancellation to fear
            logP[i] = q > Double.MIN_NORMAL ? Math.log(q) : Math.log(p[i]) - logP0;
        }

        this.cum = new double[m];
        double running = 0.0;
        for (int i = 0; i < m; i++) {
            running += p[i];
            cum[i] = running / sum;
        }
        // assigned rather than computed, so that the last category closes the
        // distribution function exactly and the search in quantile terminates
        cum[m - 1] = 1.0;

        int first = -1;
        int last = -1;
        double weighted = 0.0;
        for (int i = 0; i < m; i++) {
            if (p[i] > 0.0) {
                if (first < 0) {
                    first = i;
                }
                last = i;
            }
            weighted += i * p[i];
        }
        this.lo = first;
        this.hi = last;
        this.mean = weighted / sum;

        // two passes rather than E[X*X] - E[X]*E[X], which cancels for a law
        // that sits far from the origin. Here that is prudence rather than a
        // rescue: the categories are indices into this array, so they cannot
        // reach the magnitude where a square stops being exact, and both forms
        // in fact agree over every support a double[] of weights can span
        double spread = 0.0;
        for (int i = 0; i < m; i++) {
            double d = i - mean;
            spread += p[i] * d * d;
        }
        this.variance = spread / sum;

        this.table = AliasTable.of(p);
    }

    /**
     * The number of categories.
     *
     * @return the number of categories, one or more
     */
    public int categories() {
        return p.length;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The weight of category {@code k} over their sum, taken as the quotient
     * rather than as the exponential of {@link #logPmf(int)}: here the mass is
     * what was given and the logarithm is what is derived, which is the reverse
     * of every other law in this package.
     */
    @Override
    public double pmf(int k) {
        if (k < 0 || k >= p.length) {
            return 0.0;
        }
        return p[k] / p0;
    }

    /**
     * {@inheritDoc}
     * <p>
     * It answers where the mass cannot: a weight small enough that its share
     * underflows to zero still has a finite logarithm, and this returns it.
     */
    @Override
    public double logPmf(int k) {
        if (k < 0 || k >= p.length) {
            return Double.NEGATIVE_INFINITY;
        }
        return logP[k];
    }

    /**
     * {@inheritDoc}
     * <p>
     * The running sum of the mass, precomputed. It is exactly {@code 1.0} at
     * the last category whatever the rounding of the individual shares.
     */
    @Override
    public double cdf(int k) {
        if (k < 0) {
            return 0.0;
        }
        if (k >= p.length) {
            return 1.0;
        }
        return cum[k];
    }

    /**
     * {@inheritDoc}
     * <p>
     * A binary search over the precomputed distribution function rather than
     * the bracketing of the default, which starts from a normal approximation
     * that means nothing for an arbitrary finite law.
     * <p>
     * The answer always has positive weight. For {@code p > 0} the smallest
     * {@code k} with {@code cdf(k) >= p} sits at a jump of the distribution
     * function by construction, and {@code p == 0} is answered with
     * {@link #supportLowerBound()}. Its {@link #pmf(int)} may still be zero,
     * but only where the share has underflowed and {@link #logPmf(int)}
     * answers instead.
     */
    @Override
    public int quantile(double p) {
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + p);
        }
        if (p == 0.0) {
            return lo;
        }
        if (p == 1.0) {
            return hi;
        }
        int low = 0;
        int high = cum.length - 1;
        while (low < high) {
            int mid = (low + high) >>> 1;
            if (cum[mid] >= p) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return low;
    }

    @Override
    public double mean() {
        return mean;
    }

    @Override
    public double variance() {
        return variance;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The first category of positive weight, which need not be {@code 0}: a
     * weight of zero is not a small probability but no support at all.
     */
    @Override
    public int supportLowerBound() {
        return lo;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The last category of positive weight. Between it and
     * {@link #supportLowerBound()} the support may still have holes, one for
     * every interior weight of zero -- it is a set of categories, not a range.
     */
    @Override
    public int supportUpperBound() {
        return hi;
    }

    /**
     * The prepared table this distribution draws from, so that a caller who
     * wants many draws can stream them:
     * {@code prng.categorical(dist.aliasTable())} takes them from the table
     * that was built here rather than building a second one.
     *
     * @return the alias table over these categories
     */
    public AliasTable aliasTable() {
        return table;
    }

    /**
     * Draws one category.
     * <p>
     * The draw is the one {@link AliasTable} makes, which this distribution
     * holds so that a caller need not build one.
     *
     * @param prng
     *            the generator to draw from
     * @return a category in {@code 0 .. categories() - 1}
     * @throws IllegalArgumentException
     *             if {@code prng} is {@code null}
     */
    public int sample(PseudoRandom prng) {
        return table.sample(prng);
    }
}
