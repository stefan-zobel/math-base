package math.distribution;

import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.rng.PseudoRandom;

/**
 * The multivariate normal distribution, with density
 * {@code (2 pi)^(-d/2) det(S)^(-1/2) exp(-0.5 (x-mu)' S^-1 (x-mu))}.
 * <p>
 * The missing half of a pair: the previous "MultivariateNormalSampler"
 * could draw from this law and nothing in the library could evaluate it.
 * It is also what a Laplace approximation of a posterior <em>is</em>, so
 * {@code math.stats.bayes.LaplaceApproximation} hands one back.
 * <p>
 * Like {@link Dirichlet} and {@link Multinomial} it implements no interface:
 * an outcome here is a vector, so of what {@link ContinuousDistribution}
 * declares only the density survives the move to several dimensions -- there is
 * no distribution function in closed form, no quantile, and a mean and a
 * variance are a vector and a matrix. {@link #marginal(int)} brings the one
 * dimensional apparatus back one coordinate at a time.
 * <p>
 * <b>The covariance is used as given.</b>
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, double[])} cannot -- that restriction
 * belongs to the generator.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Multivariate_normal_distribution">
 * Wikipedia Multivariate normal distribution</a>.
 *
 * @since 1.5.3
 */
public final class MultivariateNormal {

    /**
     * How far the two halves of the covariance may differ before it is called
     * asymmetric, relative to the largest variance on its diagonal.
     */
    private static final double SYMMETRY_TOL = 1.0e-10;

    private final double[] mean;
    /** The covariance as used: the caller's, symmetrized within tolerance. */
    private final DMatrix covariance;
    /** Its Cholesky factor, lower triangular. */
    private final DMatrix chol;
    private final double logDeterminant;
    /** {@code -0.5 (d log(2 pi) + logDeterminant)}, the constant of the density. */
    private final double logNormalizer;

    /**
     * The multivariate normal with the given mean vector and covariance matrix.
     * <p>
     * The covariance has to be symmetric and positive definite. It is compared
     * against its own transpose <em>relatively</em>, to the largest entry on
     * its diagonal, and averaged within that tolerance -- a covariance
     * accumulated from data is symmetric only up to its own round-off, and at a
     * scale of {@code 1e6} that round-off is larger than the absolute
     * tolerance {@link CholeskyDecomp} applies.
     *
     * @param mean
     *            the mean of each coordinate, at least one, all finite. Not
     *            modified
     * @param covariance
     *            the covariance matrix, square, as wide as {@code mean} is
     *            long, symmetric and positive definite. Not modified
     * @throws IllegalArgumentException
     *             if either argument is {@code null}, if {@code mean} is empty
     *             or holds a value that is not finite, if {@code covariance} is
     *             not square or not of the same order, if it holds a value that
     *             is not finite, if it is not symmetric, or if it is not
     *             positive definite
     */
    public MultivariateNormal(double[] mean, DMatrix covariance) {
        if (mean == null) {
            throw new IllegalArgumentException("mean must not be null");
        }
        if (covariance == null) {
            throw new IllegalArgumentException("covariance must not be null");
        }
        int d = mean.length;
        if (d == 0) {
            throw new IllegalArgumentException("mean must not be empty");
        }
        for (int i = 0; i < d; i++) {
            if (!isFinite(mean[i])) {
                throw new IllegalArgumentException("mean[" + i + "] is not finite : " + mean[i]);
            }
        }
        if (!covariance.isSquareMatrix()) {
            throw new IllegalArgumentException("covariance is not square");
        }
        if (covariance.numRows() != d) {
            throw new IllegalArgumentException(
                    "covariance is of order " + covariance.numRows() + ", not " + d);
        }

        double scale = 0.0;
        for (int i = 0; i < d; i++) {
            double v = covariance.get(i, i);
            if (!isFinite(v)) {
                throw new IllegalArgumentException("covariance[" + i + "][" + i + "] is not finite : " + v);
            }
            scale = Math.max(scale, Math.abs(v));
        }
        if (!(scale > 0.0)) {
            throw new IllegalArgumentException("the covariance has no positive variance on its diagonal");
        }
        double tol = SYMMETRY_TOL * scale;

        DMatrix symmetric = new DMatrix(d, d);
        for (int i = 0; i < d; i++) {
            symmetric.set(i, i, covariance.get(i, i));
            for (int j = i + 1; j < d; j++) {
                double a = covariance.get(i, j);
                double b = covariance.get(j, i);
                if (!isFinite(a) || !isFinite(b)) {
                    throw new IllegalArgumentException(
                            "covariance[" + i + "][" + j + "] is not finite : " + a + ", " + b);
                }
                if (Math.abs(a - b) > tol) {
                    throw new IllegalArgumentException("the covariance is not symmetric at (" + i + ", " + j
                            + ") : " + a + " against " + b);
                }
                // averaged rather than taken from one side, so that the result
                // is exactly symmetric whichever half was the rounder one
                double m = 0.5 * (a + b);
                symmetric.set(i, j, m);
                symmetric.set(j, i, m);
            }
        }

        DMatrix factor;
        try {
            factor = CholeskyDecomp.cholesky(symmetric);
        } catch (RuntimeException e) {
            throw new IllegalArgumentException("the covariance is not positive definite : " + e.getMessage());
        }
        for (int i = 0; i < d; i++) {
            if (!(factor.get(i, i) > 0.0)) {
                throw new IllegalArgumentException(
                        "the covariance is singular: its factor has a zero at (" + i + ", " + i + ")");
            }
        }

        this.mean = mean.clone();
        this.covariance = symmetric;
        this.chol = factor;
        this.logDeterminant = CholeskyDecomp.logDeterminant(factor);
        this.logNormalizer = -0.5 * (d * Math.log(2.0 * Math.PI) + logDeterminant);
    }

    /**
     * The number of coordinates.
     *
     * @return the dimension, one or more
     */
    public int dimension() {
        return mean.length;
    }

    /**
     * The natural logarithm of the density at {@code x}.
     * <p>
     * Formed through the Cholesky factor and never through an inverse: the
     * quadratic form is {@code r'r} where {@code L r = x - mu}, which
     * {@link CholeskyDecomp#solveLower} supplies. Multiplying by
     * {@link DMatrix#inverse()} would be a second factorization and would
     * square the condition number of the covariance, for a quantity the factor
     * already holds.
     * <p>
     * It answers where the density cannot: fifty coordinates already put the
     * density of an ordinary point below the smallest {@code double}.
     * <p>
     * A multivariate normal has a density everywhere in {@code R^d}, so there
     * is no support to fall off; a {@code NaN} coordinate is outside the domain
     * rather than outside the support and comes back as {@code NaN}.
     *
     * @param x
     *            the point, as long as this law has coordinates. Not modified
     * @return the natural logarithm of the density there
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is not as long as this law
     *             has coordinates
     */
    public double logPdf(double[] x) {
        return logNormalizer - 0.5 * squaredMahalanobisDistance(x);
    }

    /**
     * The density at {@code x}.
     * <p>
     * The exponential of {@link #logPdf(double[])}, which is where the accuracy
     * is and where the answer survives.
     *
     * @param x
     *            the point, as long as this law has coordinates. Not modified
     * @return the density there
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is not as long as this law
     *             has coordinates
     */
    public double pdf(double[] x) {
        return Math.exp(logPdf(x));
    }

    /**
     * The <b>squared</b> Mahalanobis distance from the mean,
     * {@code (x-mu)' S^-1 (x-mu)}.
     * <p>
     * Under this law it is distributed as a {@link ChiSquare} on
     * {@link #dimension()} degrees of freedom, which is what makes it useful:
     * {@code new ChiSquare(d).cdf(mvn.squaredMahalanobisDistance(x))} is how
     * unusual {@code x} is, and the same quantity draws the confidence
     * ellipsoids.
     *
     * @param x
     *            the point, as long as this law has coordinates. Not modified
     * @return the squared Mahalanobis distance, not negative
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is not as long as this law
     *             has coordinates
     */
    public double squaredMahalanobisDistance(double[] x) {
        int d = mean.length;
        if (x == null) {
            throw new IllegalArgumentException("x must not be null");
        }
        if (x.length != d) {
            throw new IllegalArgumentException("x must hold " + d + " coordinates, not " + x.length);
        }
        // its own array, because this is evaluated from several threads at once
        // wherever the caller shares the instance, which is what immutability
        // is for
        double[] r = new double[d];
        for (int i = 0; i < d; i++) {
            r[i] = x[i] - mean[i];
        }
        CholeskyDecomp.solveLower(chol, r, r);
        double sum = 0.0;
        for (int i = 0; i < d; i++) {
            sum += r[i] * r[i];
        }
        return sum;
    }

    /**
     * Writes the mean into {@code out}.
     *
     * @param out
     *            where the result is written, one entry per coordinate. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or is not as long as this law
     *             has coordinates
     */
    public void mean(double[] out) {
        if (out == null || out.length != mean.length) {
            throw new IllegalArgumentException("out must hold " + mean.length + " coordinates");
        }
        System.arraycopy(mean, 0, out, 0, mean.length);
    }

    /**
     * The covariance matrix as this law uses it, which is the caller's
     * symmetrized and is exactly symmetric.
     *
     * @return a fresh copy of the covariance
     */
    public DMatrix covariance() {
        return covariance.copy();
    }

    /**
     * The natural logarithm of the determinant of the covariance.
     * <p>
     * Answers where the determinant cannot: fifty coordinates of a covariance
     * scaled by {@code 1e-14} have a determinant below the smallest
     * {@code double} and a logarithm of about {@code -1600}.
     *
     * @return {@code log det(S)}
     */
    public double logDeterminant() {
        return logDeterminant;
    }

    /**
     * The law of one coordinate on its own, which is a
     * {@code Normal(mu_i, sqrt(S_ii))}.
     * <p>
     * This is where everything a one dimensional distribution offers comes
     * back, as {@link Dirichlet#marginal(int)} and
     * {@link Multinomial#marginal(int)} do for their laws: the returned object
     * has a distribution function and a quantile.
     *
     * @param i
     *            the coordinate
     * @return the marginal law of coordinate {@code i}
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not a coordinate
     */
    public Normal marginal(int i) {
        if (i < 0 || i >= mean.length) {
            throw new IndexOutOfBoundsException("no such coordinate : " + i);
        }
        return new Normal(mean[i], Math.sqrt(covariance.get(i, i)));
    }

    /**
     * Draws one point into {@code out}.
     * <p>
     * The mean plus the Cholesky factor times a vector of independent standard
     * normals.
     *
     * @param prng
     *            the generator to draw from
     * @param out
     *            where the result is written, one entry per coordinate. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code prng} or {@code out} is {@code null}, or if
     *             {@code out} is not as long as this law has coordinates
     */
    public void sample(PseudoRandom prng, double[] out) {
        int d = mean.length;
        if (prng == null) {
            throw new IllegalArgumentException("prng must not be null");
        }
        if (out == null || out.length != d) {
            throw new IllegalArgumentException("out must hold " + d + " coordinates");
        }
        double[] z = new double[d];
        for (int i = 0; i < d; i++) {
            z[i] = prng.nextGaussian();
        }
        for (int i = 0; i < d; i++) {
            double sum = mean[i];
            for (int k = 0; k <= i; k++) {
                sum += chol.getUnsafe(i, k) * z[k];
            }
            out[i] = sum;
        }
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
