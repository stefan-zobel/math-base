/*
 * Copyright 2013, 2026 Stefan Zobel
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

/**
 * Snedecor's F distribution.
 * <p>
 * The degrees of freedom are not required to be whole numbers, as they are
 * not for {@link StudentT} and {@link ChiSquare} either: Welch's
 * approximation to a k-sample test of equal means puts a fractional value in
 * the denominator, and a distribution that refused it could not answer that
 * test.
 * <p>
 * <b>See</b> <a href="https://en.wikipedia.org/wiki/F-distribution">Wikipedia
 * F-distribution</a>.
 */
public class FisherF implements ContinuousDistribution {

    private final double d1; // numerator DF
    private final double d2; // denominator DF
    private final Beta beta;

    /**
     * An F distribution with the given degrees of freedom, neither of which has
     * to be a whole number.
     *
     * @param numeratorDF
     *            degrees of freedom of the numerator, strictly positive and
     *            finite
     * @param denominatorDF
     *            degrees of freedom of the denominator, strictly positive and
     *            finite
     * @throws IllegalArgumentException
     *             if either argument is not strictly positive and finite
     */
    public FisherF(double numeratorDF, double denominatorDF) {
        // written this way round so that NaN is refused with the rest
        if (!(numeratorDF > 0.0 && numeratorDF < Double.POSITIVE_INFINITY)) {
            throw new IllegalArgumentException("numeratorDF must be positive and finite : " + numeratorDF);
        }
        if (!(denominatorDF > 0.0 && denominatorDF < Double.POSITIVE_INFINITY)) {
            throw new IllegalArgumentException(
                    "denominatorDF must be positive and finite : " + denominatorDF);
        }
        this.d1 = numeratorDF;
        this.d2 = denominatorDF;
        this.beta = new Beta((d1 / 2.0), (d2 / 2.0));
    }

    @Override
    public double pdf(double x) {
        if (x < 0.0) {
            return 0.0;
        }
        if (x == 0.0) {
            if (d1 < 2.0) {
                return Double.POSITIVE_INFINITY;
            }
            if (d1 == 2.0) {
                return 1.0;
            }
            return 0.0;
        }
        // A quite clever variable substitution approach:
        final double w = d2 / (d2 + (d1 * x));
        // Fact: if X ~ F(d1, d1) then (1 - W) ~ Beta(d1/2, d2/2).
        //
        // Further note that (1): (1 - w) = (d1 / d2) * x * w
        // and (2): (1 / w) = 1 + (d1 / d2) * x
        //
        // First write out the density of the Beta((1-w); d1/2, d2/2)
        // and then substitute (1) into the resulting (1 - w) term.
        //
        // Then multiply the density by (w * w * (d1/d2)); finally
        // replace the remaining "w" term with the inverse of (2).
        //
        // Compare your result with the density of the F(x; d1, d2)
        // - both are identical! This proves that the following is
        // the correct transformation:
        final double scaled = (w * d1) * w;
        if (scaled == 0.0) {
            // The prefactor has underflowed. For d2 < 2 the Beta density has
            // a pole at one, so the product would read 0 * infinity; the
            // exponent of w wins it and the density is zero this far out.
            return 0.0;
        }
        return (scaled * beta.pdf(1.0 - w)) / d2;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The same substitution {@link #pdf(double)} uses, in logarithms. The
     * prefactor {@code d1 * w * w} underflows to zero far out in the tail --
     * which is the branch the density has to answer zero from -- while the sum
     * of its logarithms is an ordinary number there.
     */
    @Override
    public double logPdf(double x) {
        if (x < 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (x == 0.0) {
            if (d1 < 2.0) {
                return Double.POSITIVE_INFINITY;
            }
            return (d1 == 2.0) ? 0.0 : Double.NEGATIVE_INFINITY;
        }
        final double w = d2 / (d2 + (d1 * x));
        return Math.log(d1) + 2.0 * Math.log(w) + beta.logPdf(1.0 - w) - Math.log(d2);
    }

    @Override
    public double cdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        final double z = d1 * x;
        if (z > Double.MAX_VALUE) {
            // d1 * x overflowed, so z / (d2 + z) would read inf / inf. The
            // ratio is one out there, and so is the distribution function.
            return 1.0;
        }
        final double y = z / (d2 + z);
        return beta.cdf(y);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        return findRoot(probability, startingPoint(probability), 0.0, Double.MAX_VALUE);
    }

    /**
     * A finite point for {@link #inverseCdf(double)} to start its search from.
     * The mean serves whenever it exists, but for {@code d2 <= 2} it does not,
     * and seeding the search with {@code NaN} made every quantile {@code NaN} --
     * a quantile that exists refused because a moment that does not was asked
     * for it. The fallback inverts the Beta relation that {@link #cdf(double)}
     * is built on, which is exact apart from the cancellation in {@code 1 - y}
     * far out in the upper tail; the search then removes that.
     *
     * @param probability
     *            the probability whose quantile is being sought
     * @return a finite starting point inside the support
     */
    private double startingPoint(double probability) {
        double start = mean();
        if (!Double.isNaN(start)) {
            return start;
        }
        double y = beta.inverseCdf(probability);
        if (y >= 1.0) {
            return Double.MAX_VALUE;
        }
        start = (d2 * y) / (d1 * (1.0 - y));
        return (start > 0.0) ? start : Double.MIN_NORMAL;
    }

    @Override
    public double mean() {
        if (d2 <= 2.0) {
            return Double.NaN;
        }
        return d2 / (d2 - 2.0);
    }

    @Override
    public double variance() {
        if (d2 <= 4.0) {
            return Double.NaN;
        }
        final double z = d2 - 2.0;
        return 2.0 * d2 * d2 * (d1 + z) / (d1 * z * z * (d2 - 4.0));
    }

    @Override
    public double supportLowerBound() {
        return 0.0;
    }

    /**
     * Degrees of freedom of the numerator.
     *
     * @return degrees of freedom of the numerator
     */
    public double getNumeratorDegreesOfFreedom() {
        return d1;
    }

    /**
     * Degrees of freedom of the denominator.
     *
     * @return degrees of freedom of the denominator
     */
    public double getDenominatorDegreesOfFreedom() {
        return d2;
    }
}
