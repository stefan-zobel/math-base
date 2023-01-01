/*
 * Copyright 2013 Stefan Zobel
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
 * A continuous distribution.
 */
public interface ContinuousDistribution {

    /**
     * Returns the probability density function (PDF) of this distribution
     * evaluated at the specified point {@code x}. In general, the PDF is the
     * derivative of the {@link #cdf(double) CDF}. If the derivative does not
     * exist at {@code x}, then an appropriate replacement should be returned,
     * e.g. {@code Double.POSITIVE_INFINITY}, {@code Double.NaN}, or the limit
     * inferior or limit superior of the difference quotient.
     * 
     * @param x
     *            the point at which the PDF is evaluated
     * @return the value of the probability density function at point {@code x}
     */
    double pdf(double x);

    /**
     * For a random variable {@code X} whose values are distributed according to
     * this distribution, this method returns {@code P(X <= x)}. In other words,
     * this method represents the (cumulative) distribution function (CDF) for
     * this distribution.
     * 
     * @param x
     *            the point at which the CDF is evaluated
     * @return the probability that a random variable with this distribution
     *         takes a value less than or equal to {@code x}
     */
    double cdf(double x);

    /**
     * The inverse of {@link #cdf(double)}, i.e. a method that returns the value
     * X for which P(x&lt;=X) = {@code probability}.
     * 
     * @param probability
     *            the probability for which the inverse CDF is evaluated
     * @return the value X for which P(x&lt;=X) = {@code probability}.
     */
    double inverseCdf(double probability);

    /**
     * Use this method to get the the mean of this distribution.
     * 
     * @return the mean or {@code Double.NaN} if it is not defined
     */
    double mean();

    /**
     * Use this method to get the variance of this distribution.
     * 
     * @return the variance (possibly {@code Double.POSITIVE_INFINITY} as for
     *         certain cases in {@link StudentT}) or {@code Double.NaN} if it is
     *         not defined
     */
    double variance();

    /**
     * For a random variable {@code X} whose values are distributed according to
     * this distribution, this method returns {@code P(x0 < X <= x1)}.
     * 
     * @param x0
     *            Lower bound (excluded).
     * @param x1
     *            Upper bound (included).
     * @return the probability that a random variable with this distribution
     *         takes a value between {@code x0} and {@code x1}, excluding the
     *         lower and including the upper endpoint.
     * @throws IllegalArgumentException
     *             if {@code x0 > x1}.
     * 
     *             The default implementation uses the identity
     *             {@code P(x0 < X <= x1) = P(X <= x1) - P(X <= x0)}
     */
    default double probability(double x0, double x1) {
        if (x0 > x1) {
            throw new IllegalArgumentException(
                    "Lower endpoint (" + x0 + ") must be less than or equal to upper endpoint (" + x1 + ")");
        }
        return cdf(x1) - cdf(x0);
    }

    static final double FINDROOT_ACCURACY = 1.0e-15;
    static final int FINDROOT_MAX_ITERATIONS = 150;

    /**
     * This method approximates the value of {@code x} for which
     * {@code P(X <= x) = p} where {@code p} is a given probability.
     * <p>
     * It applies a combination of the Newton-Raphson algorithm and the
     * bisection method to the value {@code start} as a starting point.
     * <p>
     * Furthermore, to ensure convergence and stability, the caller must supply
     * an interval {@code [xMin, xMax]} in which the probability distribution
     * reaches the value {@code p}.
     * <p>
     * Caution: this method does not check its arguments! It will produce wrong
     * results if bad values for the parameters are supplied. To be used with
     * care!
     * 
     * @param p
     *            the given probability for which we want to find the
     *            corresponding value of {@code x} such that
     *            {@code P(X <= x) = p}
     * @param start
     *            an initial guess that must lie in the interval
     *            {@code [xMin, xMax]} as a starting point for the search for
     *            {@code x}
     * @param xMin
     *            lower bound for an interval that must contain the searched
     *            {@code x}
     * @param xMax
     *            upper bound for an interval that must contain the searched
     *            {@code x}
     * @return an approximation for the value of {@code x} for which
     *         {@code P(X <= x) = p}
     */
    default double findRoot(double p, double start, double xMin, double xMax) {
        double x = start;
        double xNew = start;
        double dx = 1.0;
        int i = 0;
        while (Math.abs(dx) > FINDROOT_ACCURACY && i++ < FINDROOT_MAX_ITERATIONS) {
            // apply Newton-Raphson step
            double error = cdf(x) - p;
            if (error < 0.0) {
                xMin = x;
            } else {
                xMax = x;
            }
            double density = pdf(x);
            if (density != 0.0) { // avoid division by zero
                dx = error / density;
                xNew = x - dx;
            }
            // If Newton-Raphson fails to converge (which, for example, may be
            // the case if the initial guess is too rough) we apply a bisection
            // step to determine a more narrow interval around the root
            if (xNew < xMin || xNew > xMax || density == 0.0) {
                xNew = (xMin + xMax) / 2.0;
                dx = xNew - x;
            }
            x = xNew;
        }
        return x;
    }
}
