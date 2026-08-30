package math.rng;

/**
 * A prepared Dirichlet distribution: the law of a vector of proportions that
 * sum to one, and the conjugate prior of the multinomial that
 * {@link MultinomialSampler} draws from.
 * <p>
 * A draw is one {@code Gamma(alpha_i, 1)} per component, normalized by their
 * sum. <b>It is taken in logarithms</b>, and that is not a refinement: for a
 * shape below one the gamma sampler returns {@code p^(1/alpha)}, which
 * underflows to exactly zero for an ever larger share of the uniforms as the
 * shape falls. Measured, that share is 45 in 100000 at {@code alpha = 0.01},
 * 2432 at {@code 0.005} and 47562 at {@code 0.001} -- and once every component
 * of a draw underflows, the sum is zero and the whole vector comes back
 * {@link Double#NaN}. At {@code alpha = 0.001} that happened to 11 percent of
 * three-component draws.
 * <p>
 * {@link GammaSpliterator} keeps a logarithmic form for exactly this reason and
 * {@code BetaSpliterator} already leans on it. Taken that way the draw is
 * sound at every shape tried down to {@code 1e-6}, where the logarithms reach
 * {@code -1.2e7} and nothing but a logarithm could hold them.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, double[])} cannot -- that restriction
 * belongs to the generator, not to this sampler.
 * <p>
 * https://en.wikipedia.org/wiki/Dirichlet_distribution
 *
 * @since 1.5.3
 */
public final class DirichletSampler {

    private final double[] alpha;

    private DirichletSampler(double[] alpha) {
        this.alpha = alpha;
    }

    /**
     * Builds the sampler for the given concentration parameters.
     * <p>
     * They are not normalized: their sum is the concentration and decides how
     * tightly the draws gather around their mean {@code alpha_i / alpha_0},
     * where {@code alpha_0} is that sum.
     *
     * @param alpha
     *            the concentration of each component, at least one of them,
     *            each finite and strictly positive. Not modified
     * @return a sampler over {@code alpha.length} components
     * @throws IllegalArgumentException
     *             if {@code alpha} is {@code null}, is empty, or holds a value
     *             that is not finite and strictly positive
     */
    public static DirichletSampler of(double[] alpha) {
        if (alpha == null) {
            throw new IllegalArgumentException("alpha must not be null");
        }
        if (alpha.length == 0) {
            throw new IllegalArgumentException("alpha must not be empty");
        }
        for (int i = 0; i < alpha.length; i++) {
            double a = alpha[i];
            if (!(a > 0.0) || a == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException(
                        "alpha[" + i + "] is not a finite positive number : " + a);
            }
        }
        return new DirichletSampler(alpha.clone());
    }

    /**
     * Draws one vector of proportions into {@code proportions}.
     * <p>
     * The result sums to one up to round-off and every entry is positive. It
     * cannot come back as {@link Double#NaN}: the normalization subtracts the
     * largest logarithm before exponentiating, so one of the terms of the sum
     * is exactly one and the sum can never underflow.
     *
     * @param prng
     *            the generator to draw from
     * @param proportions
     *            where the result is written, one entry per component. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code prng} or {@code proportions} is {@code null}, or if
     *             {@code proportions} is not as long as there are components
     */
    public void sample(PseudoRandom prng, double[] proportions) {
        if (prng == null) {
            throw new IllegalArgumentException("prng must not be null");
        }
        if (proportions == null) {
            throw new IllegalArgumentException("proportions must not be null");
        }
        int m = alpha.length;
        if (proportions.length != m) {
            throw new IllegalArgumentException(
                    "proportions must hold " + m + " components, not " + proportions.length);
        }

        // the logarithms go into the caller's array, which keeps this method
        // free of state and therefore free to be called from several threads
        double largest = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < m; i++) {
            double logGamma = GammaSpliterator.logSample(prng, alpha[i]);
            proportions[i] = logGamma;
            if (logGamma > largest) {
                largest = logGamma;
            }
        }
        double sum = 0.0;
        for (int i = 0; i < m; i++) {
            double weight = Math.exp(proportions[i] - largest);
            proportions[i] = weight;
            sum += weight;
        }
        for (int i = 0; i < m; i++) {
            proportions[i] /= sum;
        }
    }

    /**
     * The number of components this sampler covers.
     *
     * @return the number of components, one or more
     */
    public int components() {
        return alpha.length;
    }

    /**
     * The concentration of one component, exposed for the tests that need the
     * marginal law it implies.
     *
     * @param i
     *            the component
     * @return {@code alpha_i}
     */
    double concentration(int i) {
        return alpha[i];
    }
}
