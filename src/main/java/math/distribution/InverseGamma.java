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

import math.cern.FastGamma;
import math.cern.ProbabilityFuncs;

/**
 * The inverse gamma distribution: if {@code X ~ Gamma(alpha, rate = beta)} then
 * {@code 1/X ~ InverseGamma(alpha, beta)}, with density
 * {@code f(x) = beta^alpha / Gamma(alpha) * x^(-alpha-1) * exp(-beta/x)} on
 * {@code x > 0}.
 * <p>
 * It is the conjugate prior for the variance of a normal distribution, and the
 * marginal posterior of that variance when the mean is unknown as well, which
 * is what it is here for.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Inverse-gamma_distribution">Wikipedia
 * Inverse-gamma distribution</a>.
 *
 * @since 1.5.3
 */
public class InverseGamma implements ContinuousDistribution {

    private final double shape_alpha;
    private final double scale_beta;
    /** {@code alpha * log(beta) - logGamma(alpha)}, the log normalizer. */
    private final double logNormFactor;
    /** The law of {@code 1/X}, used to start the quantile search. */
    private final Gamma reciprocal;

    /**
     * An inverse gamma distribution with the given shape and scale.
     *
     * @param shape
     *            the shape {@code alpha}, strictly positive and finite
     * @param scale
     *            the scale {@code beta}, strictly positive and finite
     * @throws IllegalArgumentException
     *             if either argument is not strictly positive and finite
     */
    public InverseGamma(double shape /* alpha */, double scale /* beta */) {
        // written this way round so that NaN is refused with the rest
        if (!(shape > 0.0 && shape < Double.POSITIVE_INFINITY)) {
            throw new IllegalArgumentException("shape must be positive and finite : " + shape);
        }
        if (!(scale > 0.0 && scale < Double.POSITIVE_INFINITY)) {
            throw new IllegalArgumentException("scale must be positive and finite : " + scale);
        }
        this.shape_alpha = shape;
        this.scale_beta = scale;
        this.logNormFactor = shape * Math.log(scale) - FastGamma.logGamma(shape);
        // Gamma takes a scale, and the reciprocal of an InverseGamma(a, b) is a
        // Gamma with rate b, which is scale 1/b
        this.reciprocal = new Gamma(shape, 1.0 / scale);
    }

    /**
     * {@inheritDoc}
     * <p>
     * The exponential of {@link #logPdf(double)}, which is where the expression
     * lives: the density underflows to zero for a shape and scale that leave
     * the logarithm an ordinary number.
     */
    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double logPdf(double x) {
        if (!(x > 0.0)) {
            // NaN reaches here too, and Math.log below would return it anyway
            return Double.isNaN(x) ? Double.NaN : Double.NEGATIVE_INFINITY;
        }
        if (x == Double.POSITIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        return logNormFactor - (shape_alpha + 1.0) * Math.log(x) - scale_beta / x;
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code Q(alpha, beta/x)}, the upper regularized incomplete gamma, taken
     * directly rather than as {@code 1 - F} of the reciprocal gamma. The
     * subtraction is what a lower tail cannot afford: at
     * {@code InverseGamma(2.5, 3).cdf(0.05)} it returns exactly zero where this
     * returns {@code 3.14e-24}.
     * <p>
     * {@link ProbabilityFuncs#gammaComplemented(double, double, double)}
     * delegates to {@code incompleteGammaComplement(b, a * x)}, so the shape
     * goes in the second argument and the first two multiply into the
     * integration limit -- hence the order below.
     */
    @Override
    public double cdf(double x) {
        if (Double.isNaN(x)) {
            // the test below is written so that NaN fails it, and would then
            // answer zero to a question that is not one
            return Double.NaN;
        }
        if (!(x > 0.0)) {
            return 0.0;
        }
        if (x == Double.POSITIVE_INFINITY) {
            return 1.0;
        }
        return ProbabilityFuncs.gammaComplemented(scale_beta, shape_alpha, 1.0 / x);
    }

    /**
     * {@inheritDoc}
     * <p>
     * The search starts from the exact relation {@code X = 1 / Y} with
     * {@code Y} the reciprocal gamma, which is the quantile itself apart from
     * the cancellation in {@code 1 - probability}; the iteration then removes
     * that. This is the shape {@link FisherF} uses for the same reason, and it
     * is also what serves for {@code alpha <= 1}, where the mean does not exist
     * and cannot start anything.
     */
    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        double y = reciprocal.inverseCdf(1.0 - probability);
        double start = (y > 0.0) ? 1.0 / y : Double.MIN_NORMAL;
        if (!(start > 0.0) || start == Double.POSITIVE_INFINITY) {
            start = (shape_alpha > 1.0) ? mean() : scale_beta / (shape_alpha + 1.0);
        }
        return findRoot(probability, start, 0.0, Double.MAX_VALUE);
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code beta / (alpha - 1)}, which exists only for {@code alpha > 1}.
     */
    @Override
    public double mean() {
        if (shape_alpha <= 1.0) {
            return Double.NaN;
        }
        return scale_beta / (shape_alpha - 1.0);
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code beta^2 / ((alpha - 1)^2 (alpha - 2))}, which exists only for
     * {@code alpha > 2}.
     */
    @Override
    public double variance() {
        if (shape_alpha <= 2.0) {
            return Double.NaN;
        }
        double gap = shape_alpha - 1.0;
        // the ratio before the square, so that a scale beyond 1e154 does not
        // overflow on its way to a variance that is representable
        double ratio = scale_beta / gap;
        return (ratio * ratio) / (shape_alpha - 2.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double supportLowerBound() {
        return 0.0;
    }

    /**
     * The shape parameter {@code alpha}.
     *
     * @return the shape parameter
     */
    public double getShape() {
        return shape_alpha;
    }

    /**
     * The scale parameter {@code beta}.
     *
     * @return the scale parameter
     */
    public double getScale() {
        return scale_beta;
    }
}
