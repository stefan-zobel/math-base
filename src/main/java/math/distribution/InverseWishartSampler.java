package math.distribution;

import math.linalg.DMatrix;
import math.rng.PseudoRandom;

/**
 * A sampler for the inverse-Wishart distribution, the conjugate prior for the
 * covariance matrix of a multivariate normal -- and therefore the form in which
 * a random covariance is usually wanted.
 * <p>
 * A draw is the inverse of a {@link WishartSampler} draw taken with the
 * <b>inverted</b> scale matrix. That inversion is the step this class exists
 * for: it is easy to leave out, it leaves a sampler that runs and produces
 * plausible matrices, and the result is a different distribution. Doing it in
 * one named place is cheaper than remembering it.
 * <p>
 * If {@code W} is drawn here with scale {@code Psi} and {@code nu} degrees of
 * freedom then {@code E[W] = Psi / (nu - p - 1)}, which exists only for
 * {@code nu > p + 1}. Smaller degrees of freedom are still drawn from -- the
 * distribution is perfectly well defined there and only its mean is not.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, DMatrix)} cannot.
 * <p>
 * https://en.wikipedia.org/wiki/Inverse-Wishart_distribution
 *
 * @since 1.5.3
 */
public final class InverseWishartSampler {

    private final WishartSampler inverseScale;

    private InverseWishartSampler(WishartSampler inverseScale) {
        this.inverseScale = inverseScale;
    }

    /**
     * Builds the sampler for the given scale matrix and degrees of freedom.
     *
     * @param scale
     *            the scale matrix, square, symmetric and positive definite.
     *            Not modified
     * @param nu
     *            the degrees of freedom, strictly greater than
     *            {@code scale.numRows() - 1}. It need not be a whole number
     * @return a sampler for the inverse-Wishart with that scale and those
     *         degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code scale} is {@code null}, is not square, is not
     *             symmetric, is not positive definite or cannot be inverted, or
     *             if {@code nu} is not greater than one less than the dimension
     */
    public static InverseWishartSampler of(DMatrix scale, double nu) {
        if (scale == null) {
            throw new IllegalArgumentException("scale must not be null");
        }
        if (!scale.isSquareMatrix()) {
            throw new IllegalArgumentException("the scale matrix is not square");
        }
        DMatrix inverted;
        try {
            inverted = scale.inverse();
        } catch (RuntimeException e) {
            if (e instanceof IllegalArgumentException) {
                throw e;
            }
            throw new IllegalArgumentException("the scale matrix cannot be inverted", e);
        }
        // the inverse of a symmetric matrix is symmetric, but it is computed by
        // a solve and comes back symmetric only to round-off. The Cholesky
        // inside the Wishart sampler tests for exact symmetry, so the two
        // halves are averaged first
        symmetrize(inverted);
        return new InverseWishartSampler(WishartSampler.of(inverted, nu));
    }

    private static void symmetrize(DMatrix m) {
        int p = m.numRows();
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < i; j++) {
                double average = 0.5 * (m.getUnsafe(i, j) + m.getUnsafe(j, i));
                m.setUnsafe(i, j, average);
                m.setUnsafe(j, i, average);
            }
        }
    }

    /**
     * Draws one matrix into {@code out}.
     *
     * @param prng
     *            the generator to draw from
     * @param out
     *            where the result is written, of the dimension of the scale
     *            matrix. Its previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code prng} or {@code out} is {@code null}, or if
     *             {@code out} is not of the right shape
     */
    public void sample(PseudoRandom prng, DMatrix out) {
        if (out == null) {
            throw new IllegalArgumentException("out must not be null");
        }
        int p = inverseScale.dimension();
        if (out.numRows() != p || out.numColumns() != p) {
            throw new IllegalArgumentException(
                    "out must be " + p + " by " + p + ", not " + out.numRows() + " by " + out.numColumns());
        }
        DMatrix drawn = sample(prng);
        for (int j = 0; j < p; j++) {
            for (int i = 0; i < p; i++) {
                out.setUnsafe(i, j, drawn.getUnsafe(i, j));
            }
        }
    }

    /**
     * Draws one matrix.
     *
     * @param prng
     *            the generator to draw from
     * @return a fresh matrix distributed as the inverse-Wishart with this
     *         sampler's scale and degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code prng} is {@code null}
     */
    public DMatrix sample(PseudoRandom prng) {
        DMatrix drawn = inverseScale.sample(prng).inverse();
        // as above: the solve behind the inversion is not exactly symmetric,
        // and a covariance matrix that is asymmetric in its last digit is a
        // nuisance to every consumer of it
        symmetrize(drawn);
        return drawn;
    }

    /**
     * The dimension of the matrices this sampler draws.
     *
     * @return the dimension, one or more
     */
    public int dimension() {
        return inverseScale.dimension();
    }

    /**
     * The degrees of freedom this sampler draws with.
     *
     * @return the degrees of freedom
     */
    public double degreesOfFreedom() {
        return inverseScale.degreesOfFreedom();
    }
}
