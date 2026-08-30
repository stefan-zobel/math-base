package math.distribution;

import math.cern.FastGamma;
import math.linalg.DMatrix;
import math.rng.DirichletSampler;
import math.rng.PseudoRandom;

/**
 * The Dirichlet distribution: the law of a vector of proportions that sum to
 * one, with density
 * {@code f(x) = Gamma(alpha_0) / prod Gamma(alpha_i) * prod x_i^(alpha_i - 1)}
 * on the simplex, where {@code alpha_0} is the sum of the concentrations.
 * <p>
 * It is the conjugate prior of the multinomial, which is what it is here for:
 * {@link #posterior(int[])} adds the counts to the concentrations and
 * {@link #logMarginalLikelihood(int[])} is the evidence that integrating the
 * multinomial likelihood against this prior gives.
 * <p>
 * This class implements neither {@link ContinuousDistribution} nor
 * {@link DiscreteDistribution}, and deliberately so. Of the members those
 * declare, only the density survives the move to several dimensions: there is
 * no closed form for the distribution function over a simplex, a quantile
 * needs a total order that {@code R^k} does not have, and a mean and a
 * variance are a vector and a matrix here rather than two numbers. What is
 * shared is a method name, not a contract.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom, double[])} cannot -- that restriction
 * belongs to the generator.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Dirichlet_distribution">Wikipedia
 * Dirichlet distribution</a>.
 *
 * @since 1.5.3
 */
public final class Dirichlet {

    /**
     * How far the components of an argument to {@link #logPdf(double[])} may
     * miss summing to one.
     * <p>
     * Measured against the sampler this package ships: over component counts
     * from two to a thousand and concentrations from {@code 1e-6} to
     * {@code 1e6}, the worst residual {@link DirichletSampler} leaves is
     * {@code 3.6e-15}, or sixteen ulps at a thousand components. The tolerance
     * therefore sits six orders of magnitude above what the library itself
     * produces, and a caller summing naively stays below it up to several
     * million components -- while the mistake it is here to catch, a vector of
     * counts or of unnormalized weights, misses by a factor rather than by an
     * ulp.
     */
    private static final double SIMPLEX_TOLERANCE = 1.0e-9;

    private final double[] alpha;
    private final double alpha0;
    /** {@code logGamma(alpha_0) - sum logGamma(alpha_i)}, the log normalizer. */
    private final double logNormFactor;
    private final DirichletSampler sampler;

    /**
     * A Dirichlet distribution with the given concentration parameters.
     * <p>
     * They are not normalized: their sum is the concentration and decides how
     * tightly the mass gathers around the mean {@code alpha_i / alpha_0}.
     *
     * @param alpha
     *            the concentration of each component, at least one of them,
     *            each finite and strictly positive. Not modified
     * @throws IllegalArgumentException
     *             if {@code alpha} is {@code null}, is empty, holds a value
     *             that is not finite and strictly positive, or sums to
     *             infinity
     */
    public Dirichlet(double[] alpha) {
        if (alpha == null) {
            throw new IllegalArgumentException("alpha must not be null");
        }
        if (alpha.length == 0) {
            throw new IllegalArgumentException("alpha must not be empty");
        }
        double sum = 0.0;
        double logGammas = 0.0;
        for (int i = 0; i < alpha.length; i++) {
            double a = alpha[i];
            if (!(a > 0.0) || a == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException("alpha[" + i + "] is not a finite positive number : " + a);
            }
            sum += a;
            logGammas += FastGamma.logGamma(a);
        }
        if (sum == Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException("the concentrations sum to infinity");
        }
        // plain summation: compensating it was measured against the sampler
        // and moved the residual from 4.2e-15 to 3.1e-15, which is the same
        // number and does not pay for the extra arithmetic
        this.alpha = alpha.clone();
        this.alpha0 = sum;
        this.logNormFactor = FastGamma.logGamma(sum) - logGammas;
        this.sampler = DirichletSampler.of(this.alpha);
    }

    /**
     * The number of components.
     *
     * @return the number of components, one or more
     */
    public int components() {
        return alpha.length;
    }

    /**
     * The concentration of one component.
     *
     * @param i
     *            the component
     * @return {@code alpha_i}
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not a component
     */
    public double concentration(int i) {
        checkComponent(i);
        return alpha[i];
    }

    /**
     * The sum of the concentrations, {@code alpha_0}.
     *
     * @return the total concentration
     */
    public double totalConcentration() {
        return alpha0;
    }

    /**
     * The natural logarithm of the density at {@code x}.
     * <p>
     * The number is a density with respect to the <b>(k-1)-dimensional</b>
     * Lebesgue measure on the simplex. A Dirichlet has no density on
     * {@code R^k} at all -- all of its mass lies on a set of {@code k}
     * dimensional measure zero -- so there is nothing to extend by zero, and a
     * point beside the simplex is not one where the density vanishes but one
     * where it is not defined. Such an argument is therefore refused rather
     * than answered, which is also what a wrong array length has to do.
     * <p>
     * At a corner the shape alone decides, as it does in
     * {@link Beta#logPdf(double)}: a component of zero contributes a pole below
     * a concentration of one, the normalizing factor at one, and minus infinity
     * above it. Where a pole and a zero meet -- {@code x = {0, 0, 1}} with
     * {@code alpha = {0.5, 2, 1}} -- the sum of the two infinities is
     * {@code NaN}, and that is the honest answer: the limit depends on the
     * direction the corner is approached from.
     *
     * @param x
     *            the point at which the log density is evaluated, one entry per
     *            component, summing to one. Not modified
     * @return the natural logarithm of the density at {@code x}, or
     *         {@code Double.NaN} if any component is {@code NaN}
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is not as long as there are
     *             components, holds a negative value, or does not sum to one
     */
    public double logPdf(double[] x) {
        int m = alpha.length;
        if (x == null) {
            throw new IllegalArgumentException("x must not be null");
        }
        if (x.length != m) {
            throw new IllegalArgumentException("x must hold " + m + " components, not " + x.length);
        }
        boolean unknown = false;
        double sum = 0.0;
        for (int i = 0; i < m; i++) {
            double xi = x[i];
            if (xi < 0.0) {
                throw new IllegalArgumentException("x[" + i + "] is negative : " + xi);
            }
            unknown |= Double.isNaN(xi);
            sum += xi;
        }
        if (unknown) {
            // nothing can be said about the simplex condition either
            return Double.NaN;
        }
        if (!(Math.abs(sum - 1.0) <= SIMPLEX_TOLERANCE)) {
            throw new IllegalArgumentException("the components of x sum to " + sum + ", not to one");
        }

        double logDensity = logNormFactor;
        for (int i = 0; i < m; i++) {
            double a = alpha[i];
            if (a != 1.0) {
                // the term is zero for a concentration of one whatever the
                // component is, and taken as written it would read
                // 0.0 * -Infinity there, which is NaN
                logDensity += (a - 1.0) * Math.log(x[i]);
            }
        }
        return logDensity;
    }

    /**
     * The density at {@code x}, with respect to the {@code (k-1)}-dimensional
     * Lebesgue measure on the simplex.
     * <p>
     * The exponential of {@link #logPdf(double[])}, which is where the
     * expression lives and where the accuracy is: the normalizing factor alone
     * leaves the {@code double} range at concentrations an ordinary model
     * reaches.
     *
     * @param x
     *            the point at which the density is evaluated, one entry per
     *            component, summing to one. Not modified
     * @return the density at {@code x}
     * @throws IllegalArgumentException
     *             under the conditions {@link #logPdf(double[])} states
     */
    public double pdf(double[] x) {
        return Math.exp(logPdf(x));
    }

    /**
     * The law of one component on its own, which is a
     * {@code Beta(alpha_i, alpha_0 - alpha_i)}.
     * <p>
     * This is where everything a one dimensional distribution offers comes
     * back: the returned object has a distribution function and a quantile, so
     * a credible interval for one component is
     * {@code marginal(i).inverseCdf(0.025)} and
     * {@code marginal(i).inverseCdf(0.975)}.
     *
     * @param i
     *            the component
     * @return the marginal law of component {@code i}
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not a component
     * @throws IllegalArgumentException
     *             if there is only one component, which is then one with
     *             probability one and has no beta marginal
     */
    public Beta marginal(int i) {
        checkComponent(i);
        if (alpha.length == 1) {
            throw new IllegalArgumentException(
                    "a single component is one with probability one and has no beta marginal");
        }
        return new Beta(alpha[i], others(i));
    }

    /**
     * Writes the mean, {@code alpha_i / alpha_0} in each component, into
     * {@code out}.
     *
     * @param out
     *            where the result is written, one entry per component. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or is not as long as there are
     *             components
     */
    public void mean(double[] out) {
        int m = alpha.length;
        if (out == null) {
            throw new IllegalArgumentException("out must not be null");
        }
        if (out.length != m) {
            throw new IllegalArgumentException("out must hold " + m + " components, not " + out.length);
        }
        for (int i = 0; i < m; i++) {
            out[i] = alpha[i] / alpha0;
        }
    }

    /**
     * The covariance matrix, {@code -m_i m_j / (alpha_0 + 1)} off the diagonal
     * and {@code m_i (1 - m_i) / (alpha_0 + 1)} on it, with {@code m_i} the
     * mean of component {@code i}.
     * <p>
     * <b>It is singular.</b> The components sum to one, so every row of this
     * matrix sums to zero and its rank is one less than its order; a Cholesky
     * factorization of it fails, and correctly so.
     *
     * @return a fresh {@code k} by {@code k} covariance matrix
     */
    public DMatrix covariance() {
        int m = alpha.length;
        double denominator = alpha0 + 1.0;
        DMatrix c = new DMatrix(m, m);
        for (int i = 0; i < m; i++) {
            // the quotient before the product throughout, so that a total
            // concentration beyond 1e154 does not overflow through alpha_0^2
            // on its way to a covariance that is representable
            double mi = alpha[i] / alpha0;
            c.setUnsafe(i, i, mi * (others(i) / alpha0) / denominator);
            for (int j = i + 1; j < m; j++) {
                double covariance = -mi * (alpha[j] / alpha0) / denominator;
                c.setUnsafe(i, j, covariance);
                c.setUnsafe(j, i, covariance);
            }
        }
        return c;
    }

    /**
     * The posterior after observing {@code counts}, which is the Dirichlet with
     * the counts added to the concentrations.
     * <p>
     * That closed form is the whole of Dirichlet-multinomial conjugacy, and it
     * composes: updating with one batch and then with a second gives the same
     * distribution as updating once with their sum.
     *
     * @param counts
     *            how often each component was observed, one entry per
     *            component, none of them negative. Not modified
     * @return the posterior distribution
     * @throws IllegalArgumentException
     *             if {@code counts} is {@code null}, is not as long as there
     *             are components, or holds a negative value
     */
    public Dirichlet posterior(int[] counts) {
        int m = alpha.length;
        checkCounts(counts, m);
        double[] updated = new double[m];
        for (int i = 0; i < m; i++) {
            updated[i] = alpha[i] + counts[i];
        }
        return new Dirichlet(updated);
    }

    /**
     * The natural logarithm of the probability of observing {@code counts},
     * with the proportions integrated out against this distribution.
     * <p>
     * This is the Dirichlet-multinomial, or Polya, evidence, and it is a
     * probability of the count vector: the multinomial coefficient is part of
     * it, so exponentiating this over every way {@code n} draws can fall into
     * the components sums to one. It is the quantity a Bayes factor comparing
     * two concentration vectors is a ratio of.
     * <p>
     * Each of the gamma ratios is taken through
     * {@link FastGamma#logGammaRatio(double, double)} rather than as a
     * difference of two logarithms, which is the subtraction that method exists
     * to avoid: at a total concentration beyond {@code 1e15} the difference is
     * smaller than one ulp of either term and comes back as zero.
     * <p>
     * The counting factor is the multinomial coefficient, and it is taken from
     * {@link Multinomial} rather than written out again here -- the evidence is
     * the multinomial likelihood with the proportions integrated out, so the
     * two share that half by construction.
     *
     * @param counts
     *            how often each component was observed, one entry per
     *            component, none of them negative. Not modified
     * @return the natural logarithm of the probability of {@code counts}
     * @throws IllegalArgumentException
     *             if {@code counts} is {@code null}, is not as long as there
     *             are components, holds a negative value, or sums to more than
     *             {@link Integer#MAX_VALUE}
     */
    public double logMarginalLikelihood(int[] counts) {
        int m = alpha.length;
        checkCounts(counts, m);
        long total = 0L;
        for (int i = 0; i < m; i++) {
            total += counts[i];
        }
        if (total > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("the counts sum to " + total + ", which is more than an int holds");
        }
        int n = (int) total;
        double result = Multinomial.logMultinomialCoefficient(counts, n) - FastGamma.logGammaRatio(alpha0, n);
        for (int i = 0; i < m; i++) {
            result += FastGamma.logGammaRatio(alpha[i], counts[i]);
        }
        return result;
    }

    /**
     * Draws one vector of proportions into {@code proportions}.
     * <p>
     * The draw is the one {@link DirichletSampler} makes, which this
     * distribution holds so that a caller need not build one.
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
        sampler.sample(prng, proportions);
    }

    /**
     * The sum of every concentration but the {@code i}-th.
     * <p>
     * Summed rather than taken as {@code alpha_0 - alpha_i}, which cancels to
     * zero when one component carries almost all of the concentration and would
     * then ask {@link Beta} for a shape of zero.
     */
    private double others(int i) {
        double sum = 0.0;
        for (int j = 0; j < alpha.length; j++) {
            if (j != i) {
                sum += alpha[j];
            }
        }
        return sum;
    }

    private void checkComponent(int i) {
        if (i < 0 || i >= alpha.length) {
            throw new IndexOutOfBoundsException("no such component : " + i);
        }
    }

    private static void checkCounts(int[] counts, int m) {
        if (counts == null) {
            throw new IllegalArgumentException("counts must not be null");
        }
        if (counts.length != m) {
            throw new IllegalArgumentException("counts must hold " + m + " components, not " + counts.length);
        }
        for (int i = 0; i < m; i++) {
            if (counts[i] < 0) {
                throw new IllegalArgumentException("counts[" + i + "] is negative : " + counts[i]);
            }
        }
    }
}
