package math.rng;

/**
 * A prepared multinomial distribution: the joint law of the counts that fall
 * into each of {@code m} categories when {@code n} independent draws are placed
 * according to fixed probabilities.
 * <p>
 * A draw is a chain of conditional binomials. The first category takes
 * {@code Binomial(n, q_0)}, the second what is left of {@code n} with
 * {@code q_1}, and so on, where
 * {@code q_i = p_i / (p_i + ... + p_(m-1))} is the probability that a draw
 * still unplaced after the first {@code i} categories lands in category
 * {@code i}. Those conditional probabilities depend only on the weights, so
 * they are computed once here rather than once per draw.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, int, int[])} cannot -- that
 * restriction belongs to the generator, not to this sampler.
 * <p>
 * https://en.wikipedia.org/wiki/Multinomial_distribution
 *
 * @since 1.5.3
 */
public final class MultinomialSampler {

    /**
     * The probability that a draw not yet placed lands in this category rather
     * than in any later one. The last is one by construction.
     */
    private final double[] conditional;

    private MultinomialSampler(double[] conditional) {
        this.conditional = conditional;
    }

    /**
     * Builds the sampler for the given weights.
     * <p>
     * The weights need not sum to one; they are normalized, so counts serve as
     * well as probabilities. A weight of zero is allowed and its category is
     * then never counted.
     *
     * @param probabilities
     *            the weight of each category, at least one of them, each finite
     *            and not negative, and not all zero. Not modified
     * @return a sampler over {@code 0 .. probabilities.length - 1}
     * @throws IllegalArgumentException
     *             if {@code probabilities} is {@code null}, is empty, holds a
     *             value that is negative or not finite, or sums to zero
     */
    public static MultinomialSampler of(double[] probabilities) {
        if (probabilities == null) {
            throw new IllegalArgumentException("probabilities must not be null");
        }
        int m = probabilities.length;
        if (m == 0) {
            throw new IllegalArgumentException("probabilities must not be empty");
        }
        for (int i = 0; i < m; i++) {
            double w = probabilities[i];
            if (!(w >= 0.0) || w == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException(
                        "probabilities[" + i + "] is not a finite non-negative number : " + w);
            }
        }

        // the suffix sums, taken from the back. Writing the denominator as
        // 1 - p_0 - ... - p_(i-1) instead would lose digits to repeated
        // subtraction, and would get the degenerate case wrong: where every
        // later weight is zero the suffix sum is exactly p_i, so the quotient
        // below is exactly one, every remaining count is placed here, and a
        // trailing category of weight zero cannot be handed a leftover
        double[] conditional = new double[m];
        double suffix = 0.0;
        for (int i = m - 1; i >= 0; i--) {
            suffix += probabilities[i];
            conditional[i] = suffix > 0.0 ? probabilities[i] / suffix : 0.0;
            if (conditional[i] > 1.0) {
                conditional[i] = 1.0;
            }
        }
        if (!(suffix > 0.0) || suffix == Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException(
                    "the weights must sum to a finite positive number, not " + suffix);
        }
        return new MultinomialSampler(conditional);
    }

    /**
     * Draws the counts of one multinomial experiment into {@code counts}.
     *
     * @param prng
     *            the generator to draw from
     * @param n
     *            the number of draws to place, not negative
     * @param counts
     *            where the result is written, one entry per category. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code prng} or {@code counts} is {@code null}, if
     *             {@code n} is negative, or if {@code counts} is not as long as
     *             there are categories
     */
    public void sample(PseudoRandom prng, int n, int[] counts) {
        if (prng == null) {
            throw new IllegalArgumentException("prng must not be null");
        }
        if (counts == null) {
            throw new IllegalArgumentException("counts must not be null");
        }
        if (n < 0) {
            throw new IllegalArgumentException("n < 0 : " + n);
        }
        int m = conditional.length;
        if (counts.length != m) {
            throw new IllegalArgumentException(
                    "counts must hold " + m + " categories, not " + counts.length);
        }

        int remaining = n;
        for (int i = 0; i < m - 1; i++) {
            if (remaining == 0) {
                counts[i] = 0;
                continue;
            }
            int k = BinomialSpliterator.sample(prng, remaining, conditional[i]);
            counts[i] = k;
            remaining -= k;
        }
        // the last category takes what is left, which is what makes the counts
        // sum to n exactly rather than nearly
        counts[m - 1] = remaining;
    }

    /**
     * The number of categories this sampler covers.
     *
     * @return the number of categories, one or more
     */
    public int categories() {
        return conditional.length;
    }

    /**
     * The probability that a draw still unplaced after the first {@code i}
     * categories lands in category {@code i}, exposed for the tests that check
     * the construction rather than its consequences.
     *
     * @param i
     *            the category
     * @return the conditional probability of category {@code i}
     */
    double conditionalProbability(int i) {
        return conditional[i];
    }
}
