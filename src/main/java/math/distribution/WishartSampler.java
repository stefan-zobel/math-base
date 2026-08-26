package math.distribution;

import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.rng.PseudoRandom;

/**
 * A sampler for the Wishart distribution, the law of a random covariance
 * matrix: the scatter {@code X' X} of {@code nu} draws from a centered
 * multivariate normal with covariance {@code Sigma}.
 * <p>
 * A draw is Bartlett's decomposition rather than that sum. With
 * {@code Sigma = L L'} and {@code A} lower triangular carrying
 * {@code sqrt(chi-squared with nu - i degrees of freedom)} on its diagonal and
 * standard normals below it, {@code (L A)(L A)'} has exactly the Wishart law --
 * at a cost that does not grow with {@code nu}, where forming the scatter
 * directly would need {@code nu} normal vectors per draw.
 * <p>
 * <b>The scale matrix is used as given.</b> The multivariate normal sampler in
 * this package adds {@code 1e-7} times the identity to its covariance before
 * factorizing it; that is not done here, because it would shift every draw of a
 * distribution whose whole content is a covariance. A scale matrix that is not
 * positive definite is refused instead.
 * <p>
 * Instances are immutable and can be shared between threads -- unlike
 * {@link MultivariateNormalSampler}, which owns its generator. The generator
 * passed to {@link #sample(PseudoRandom, DMatrix)} is the one that cannot be
 * shared, which is a restriction belonging to the generator rather than to this
 * sampler.
 * <p>
 * https://en.wikipedia.org/wiki/Wishart_distribution
 *
 * @since 1.5.3
 */
public final class WishartSampler {

    /** The lower Cholesky factor of the scale matrix, taken once. */
    private final DMatrix choleskyL;

    /** The degrees of freedom of each diagonal entry, {@code nu - i}. */
    private final double[] diagonalDf;

    private final double nu;

    private WishartSampler(DMatrix choleskyL, double[] diagonalDf, double nu) {
        this.choleskyL = choleskyL;
        this.diagonalDf = diagonalDf;
        this.nu = nu;
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
     * @return a sampler for {@code Wishart(scale, nu)}
     * @throws IllegalArgumentException
     *             if {@code scale} is {@code null}, is not square, is not
     *             symmetric, is not positive definite, or if {@code nu} is not
     *             greater than one less than the dimension
     */
    public static WishartSampler of(DMatrix scale, double nu) {
        if (scale == null) {
            throw new IllegalArgumentException("scale must not be null");
        }
        if (!scale.isSquareMatrix()) {
            throw new IllegalArgumentException("the scale matrix is not square");
        }
        int p = scale.numRows();
        if (Double.isNaN(nu) || nu <= p - 1.0) {
            throw new IllegalArgumentException(
                    "nu must be greater than " + (p - 1) + " for a " + p + "-dimensional scale matrix : "
                            + nu);
        }
        DMatrix choleskyL;
        try {
            choleskyL = CholeskyDecomp.cholesky(scale);
        } catch (RuntimeException e) {
            // the decomposition reports an indefinite matrix by throwing a bare
            // RuntimeException, which is not what a caller error should look
            // like on the way out of here
            if (e instanceof IllegalArgumentException) {
                throw e;
            }
            throw new IllegalArgumentException("the scale matrix is not positive definite", e);
        }
        for (int i = 0; i < p; i++) {
            if (!(choleskyL.getUnsafe(i, i) > 0.0)) {
                throw new IllegalArgumentException(
                        "the scale matrix is not positive definite: its Cholesky factor has "
                                + choleskyL.getUnsafe(i, i) + " on the diagonal at " + i);
            }
        }
        double[] diagonalDf = new double[p];
        for (int i = 0; i < p; i++) {
            diagonalDf[i] = nu - i;
        }
        return new WishartSampler(choleskyL, diagonalDf, nu);
    }

    /**
     * Draws one matrix into {@code out}.
     * <p>
     * The working matrices of the decomposition are allocated per call, so this
     * saves the caller the result allocation rather than all of them.
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
        int p = choleskyL.numRows();
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
     * @return a fresh matrix distributed as {@code Wishart(scale, nu)}
     * @throws IllegalArgumentException
     *             if {@code prng} is {@code null}
     */
    public DMatrix sample(PseudoRandom prng) {
        if (prng == null) {
            throw new IllegalArgumentException("prng must not be null");
        }
        int p = choleskyL.numRows();
        DMatrix a = new DMatrix(p, p);
        for (int i = 0; i < p; i++) {
            // the chi-squared draws come from the stream rather than from
            // inverting ChiSquare. Measured, inversion costs 2.4 to 32
            // microseconds a draw against the stream's 87 to 285 nanoseconds,
            // and at one degree of freedom -- which is the last diagonal entry
            // whenever nu equals p -- its far lower tail is wrong by four
            // orders of magnitude
            a.setUnsafe(i, i, Math.sqrt(chiSquared(prng, diagonalDf[i])));
            for (int j = 0; j < i; j++) {
                a.setUnsafe(i, j, prng.nextGaussian());
            }
        }
        // L and A are both lower triangular, so this product could be formed in
        // a sixth of the operations. The general one is used because it is
        // tested and a triangular product written here would not be
        DMatrix m = choleskyL.mul(a);
        return m.mulBTrans(m);
    }

    private static double chiSquared(PseudoRandom prng, double df) {
        return prng.chiSquare(1L, df).toArray()[0];
    }

    /**
     * The dimension of the matrices this sampler draws.
     *
     * @return the dimension, one or more
     */
    public int dimension() {
        return choleskyL.numRows();
    }

    /**
     * The degrees of freedom this sampler draws with.
     *
     * @return the degrees of freedom
     */
    public double degreesOfFreedom() {
        return nu;
    }
}
