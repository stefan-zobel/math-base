package math.distribution;

import math.cern.Arithmetic;
import math.linalg.DMatrix;
import math.rng.MultinomialSampler;
import math.rng.PseudoRandom;

/**
 * The multinomial distribution: the joint law of the counts that fall into each
 * of {@code m} categories when {@code n} independent draws are placed according
 * to fixed probabilities, with mass
 * {@code P(k) = n! / (k_1! ... k_m!) * prod p_i^k_i} over the count vectors
 * that sum to {@code n}.
 * <p>
 * It is the likelihood {@link Dirichlet} is the conjugate prior of, and the two
 * are the halves of one pair: this mass is what
 * {@link Dirichlet#logMarginalLikelihood(int[])} integrates the proportions out
 * of and what {@link Dirichlet#posterior(int[])} updates against. At
 * {@code n == 1} it is the categorical distribution, a single draw landing in a
 * single category.
 * <p>
 * This class implements neither {@link DiscreteDistribution} nor
 * {@link ContinuousDistribution}, for the reason {@link Dirichlet} gives at
 * length: an outcome here is a vector of counts, so there is no distribution
 * function to write down and no quantile to ask for, and a mean and a variance
 * are a vector and a matrix. Of what those interfaces declare only the mass
 * survives the move to several dimensions; {@link #marginal(int)} brings the
 * rest back one category at a time.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, int[])} cannot -- that restriction
 * belongs to the generator.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Multinomial_distribution">Wikipedia
 * Multinomial distribution</a>.
 *
 * @since 1.5.3
 */
public final class Multinomial {

    private final int n;
    /** The weights as given, which need not sum to one. */
    private final double[] p;
    /** Their sum, the divisor that turns a weight into a probability. */
    private final double p0;
    /** {@code log(p_i / p0)}, which is {@code -Infinity} at a weight of zero. */
    private final double[] logP;
    private final MultinomialSampler sampler;

    /**
     * A multinomial distribution over the given categories for {@code n}
     * independent draws.
     * <p>
     * The weights need not sum to one; they are normalized, so counts serve as
     * well as probabilities. A weight of zero is allowed and its category is
     * then impossible rather than merely unlikely.
     *
     * @param n
     *            the number of draws to place, not negative
     * @param probabilities
     *            the weight of each category, at least one of them, each finite
     *            and not negative, and not all zero. Not modified
     * @throws IllegalArgumentException
     *             if {@code n} is negative, or if {@code probabilities} is
     *             {@code null}, is empty, holds a value that is negative or not
     *             finite, or does not sum to a finite positive number
     */
    public Multinomial(int n, double[] probabilities) {
        if (n < 0) {
            throw new IllegalArgumentException("n < 0 : " + n);
        }
        if (probabilities == null) {
            throw new IllegalArgumentException("probabilities must not be null");
        }
        int m = probabilities.length;
        if (m == 0) {
            throw new IllegalArgumentException("probabilities must not be empty");
        }
        double sum = 0.0;
        for (int i = 0; i < m; i++) {
            double w = probabilities[i];
            if (!(w >= 0.0) || w == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException(
                        "probabilities[" + i + "] is not a finite non-negative number : " + w);
            }
            sum += w;
        }
        if (!(sum > 0.0) || sum == Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException("the weights must sum to a finite positive number, not " + sum);
        }
        this.n = n;
        this.p = probabilities.clone();
        this.p0 = sum;
        this.logP = new double[m];
        double logP0 = Math.log(sum);
        for (int i = 0; i < m; i++) {
            double q = this.p[i] / sum;
            // the quotient is the accurate form wherever it is a normal number,
            // but weights spanning hundreds of orders of magnitude underflow it
            // to zero and lose an answer that is merely small. There the
            // difference of the logarithms is exact enough: it is only reached
            // below exp(-708), where neither logarithm exceeds 709 in magnitude
            // and there is no cancellation to fear. See Categorical, which
            // carries the same guard for the same reason
            logP[i] = q > Double.MIN_NORMAL ? Math.log(q) : Math.log(this.p[i]) - logP0;
        }
        this.sampler = MultinomialSampler.of(this.p);
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
     * The number of draws one outcome places.
     *
     * @return the number of trials
     */
    public int trials() {
        return n;
    }

    /**
     * The probability of one category, normalized whether the weights the
     * constructor was given were or not.
     *
     * @param i
     *            the category
     * @return {@code p_i}, in {@code [0, 1]}
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not a category
     */
    public double probability(int i) {
        checkCategory(i);
        return p[i] / p0;
    }

    /**
     * The natural logarithm of the probability of the outcome {@code counts}.
     * <p>
     * A count vector that does not sum to {@link #trials()} is not an outcome
     * of small probability but no outcome at all, and is refused rather than
     * answered -- the same line {@link Dirichlet#logPdf(double[])} draws beside
     * the simplex. Counts are integers, so the sum is compared exactly and
     * there is no tolerance to choose.
     * <p>
     * A category of weight zero is a different matter: it lies in the support
     * of nothing, so a positive count there gives
     * {@link Double#NEGATIVE_INFINITY}, and a count of zero contributes
     * nothing at all.
     *
     * @param counts
     *            how often each category was hit, one entry per category, none
     *            of them negative, summing to {@link #trials()}. Not modified
     * @return the natural logarithm of the probability of {@code counts}
     * @throws IllegalArgumentException
     *             if {@code counts} is {@code null}, is not as long as there
     *             are categories, holds a negative value, or does not sum to
     *             {@link #trials()}
     */
    public double logPmf(int[] counts) {
        int m = p.length;
        checkCounts(counts, m);
        long total = 0L;
        for (int i = 0; i < m; i++) {
            total += counts[i];
        }
        if (total != n) {
            throw new IllegalArgumentException("the counts sum to " + total + ", not to the " + n + " trials");
        }
        double logMass = logMultinomialCoefficient(counts, n);
        for (int i = 0; i < m; i++) {
            int k = counts[i];
            if (k != 0) {
                // the power is one at a count of zero, and taken as written the
                // term would read 0.0 * -Infinity for a category of weight
                // zero, which is NaN rather than the nothing it contributes
                logMass += k * logP[i];
            }
        }
        return logMass;
    }

    /**
     * The natural logarithm of the multinomial coefficient
     * {@code n! / (k_1! ... k_m!)}, the number of orderings of {@code n}
     * labelled draws that produce the given counts.
     * <p>
     * Package private and static because
     * {@link Dirichlet#logMarginalLikelihood(int[])} needs exactly this factor:
     * the Dirichlet-multinomial evidence is this coefficient times a ratio of
     * gamma functions, and the counting half of the two belongs here rather
     * than in both places.
     * <p>
     * The caller is responsible for the arguments.
     * {@link math.cern.Arithmetic#logFactorial(int)} has no guard of its own
     * and reads its lookup table out of bounds for a negative count.
     *
     * @param counts
     *            how often each category was hit, none of them negative
     * @param n
     *            their sum, not negative
     * @return the natural logarithm of the multinomial coefficient
     */
    static double logMultinomialCoefficient(int[] counts, int n) {
        double coefficient = Arithmetic.logFactorial(n);
        for (int i = 0; i < counts.length; i++) {
            coefficient -= Arithmetic.logFactorial(counts[i]);
        }
        return coefficient;
    }

    /**
     * The probability of the outcome {@code counts}.
     * <p>
     * The exponential of {@link #logPmf(int[])}, which is where the expression
     * lives and where the accuracy is: a few thousand draws over a handful of
     * categories already put every single outcome below the smallest
     * {@code double}.
     *
     * @param counts
     *            how often each category was hit, one entry per category, none
     *            of them negative, summing to {@link #trials()}. Not modified
     * @return the probability of {@code counts}
     * @throws IllegalArgumentException
     *             under the conditions {@link #logPmf(int[])} states
     */
    public double pmf(int[] counts) {
        return Math.exp(logPmf(counts));
    }

    /**
     * The law of the count in one category on its own, which is a
     * {@code Binomial(n, p_i)}: every draw either lands there or does not.
     * <p>
     * This is where everything a one dimensional distribution offers comes
     * back. The returned object has a distribution function and a quantile, so
     * {@code marginal(i).quantile(0.975)} is an upper bound on the count in
     * category {@code i}.
     * <p>
     * Unlike {@link Dirichlet#marginal(int)} this is defined for a single
     * category as well: a {@link Beta} has no shape of zero, but a
     * {@code Binomial(n, 1)} is an ordinary degenerate law that puts every
     * draw where the only category is.
     *
     * @param i
     *            the category
     * @return the marginal law of the count in category {@code i}
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not a category
     */
    public Binomial marginal(int i) {
        checkCategory(i);
        // the weights are not negative and their sum is accumulated from them,
        // so the quotient cannot exceed one and Binomial cannot refuse it
        return new Binomial(n, p[i] / p0);
    }

    /**
     * Writes the mean, {@code n p_i} in each category, into {@code out}.
     *
     * @param out
     *            where the result is written, one entry per category. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or is not as long as there are
     *             categories
     */
    public void mean(double[] out) {
        int m = p.length;
        if (out == null) {
            throw new IllegalArgumentException("out must not be null");
        }
        if (out.length != m) {
            throw new IllegalArgumentException("out must hold " + m + " categories, not " + out.length);
        }
        for (int i = 0; i < m; i++) {
            out[i] = n * (p[i] / p0);
        }
    }

    /**
     * The covariance matrix, {@code -n p_i p_j} off the diagonal and
     * {@code n p_i (1 - p_i)} on it.
     * <p>
     * <b>It is singular.</b> The counts sum to {@code n} whatever the outcome,
     * so every row of this matrix sums to zero and its rank is one less than
     * its order; a Cholesky factorization of it fails, and correctly so. The
     * negative off-diagonal is the whole of what makes a multinomial more than
     * a row of independent binomials: a draw counted in one category is a draw
     * the others did not get.
     *
     * @return a fresh {@code m} by {@code m} covariance matrix
     */
    public DMatrix covariance() {
        int m = p.length;
        DMatrix c = new DMatrix(m, m);
        for (int i = 0; i < m; i++) {
            double mi = p[i] / p0;
            // others(i) / p0 rather than 1 - mi, which cancels to zero when one
            // category carries almost all of the weight and would report a
            // variance of zero where there is one
            c.setUnsafe(i, i, n * mi * (others(i) / p0));
            for (int j = i + 1; j < m; j++) {
                double covariance = -n * mi * (p[j] / p0);
                c.setUnsafe(i, j, covariance);
                c.setUnsafe(j, i, covariance);
            }
        }
        return c;
    }

    /**
     * Draws the counts of one experiment into {@code counts}.
     * <p>
     * The draw is the one {@link MultinomialSampler} makes, which this
     * distribution holds so that a caller need not build one.
     *
     * @param prng
     *            the generator to draw from
     * @param counts
     *            where the result is written, one entry per category. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code prng} or {@code counts} is {@code null}, or if
     *             {@code counts} is not as long as there are categories
     */
    public void sample(PseudoRandom prng, int[] counts) {
        sampler.sample(prng, n, counts);
    }

    /**
     * The sum of every weight but the {@code i}-th.
     * <p>
     * Summed rather than taken as {@code p0 - p[i]}, which cancels to zero when
     * one category carries almost all of the weight.
     */
    private double others(int i) {
        double sum = 0.0;
        for (int j = 0; j < p.length; j++) {
            if (j != i) {
                sum += p[j];
            }
        }
        return sum;
    }

    private void checkCategory(int i) {
        if (i < 0 || i >= p.length) {
            throw new IndexOutOfBoundsException("no such category : " + i);
        }
    }

    private static void checkCounts(int[] counts, int m) {
        if (counts == null) {
            throw new IllegalArgumentException("counts must not be null");
        }
        if (counts.length != m) {
            throw new IllegalArgumentException("counts must hold " + m + " categories, not " + counts.length);
        }
        for (int i = 0; i < m; i++) {
            if (counts[i] < 0) {
                throw new IllegalArgumentException("counts[" + i + "] is negative : " + counts[i]);
            }
        }
    }
}
