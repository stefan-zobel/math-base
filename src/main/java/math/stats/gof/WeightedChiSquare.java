package math.stats.gof;

import math.cern.ProbabilityFuncs;
import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod;

/**
 * The distribution of a weighted sum of independent chi-squares with one degree
 * of freedom each, {@code Q = sum lambda_j z_j^2}, by Imhof's inversion of its
 * characteristic function.
 * <p>
 * Several goodness of fit statistics are exactly of this shape once their
 * eigenvalues are known, and the eigenvalues are the only thing that
 * distinguishes them: {@code 1 / (pi^2 j^2)} gives the limit of the
 * Cramer-von Mises statistic, {@code 1 / (j (j+1))} that of the
 * Anderson-Darling statistic, and the eigenvalues of a design matrix give the
 * Durbin-Watson distribution. A ratio {@code z' B z / z' z} is the same problem
 * at {@code x = 0}, since it is at most {@code t} exactly when
 * {@code sum (nu_i - t) z_i^2} is at most zero.
 * <p>
 * The weights may be of either sign, and a zero weight is dropped -- it
 * contributes nothing to the sum and nothing to the integral.
 * <p>
 * <b>What a call costs depends on the shift and on how many weights there
 * are, and the two interact.</b> Imhof's integrand carries a term
 * {@code -x v / 2} inside the sine, so for {@code x != 0} it oscillates at a
 * fixed frequency however far out one goes, while its envelope decays only as
 * {@code v^(-1-m/2)} in the number of weights {@code m}. Few weights therefore
 * mean a slowly decaying oscillation over a very long range, covered by panels
 * that double in width, so the outer panels hold thousands of periods each:
 * eight weights and a shift take about a fifth of a second, thirty weights
 * about twelve milliseconds, and {@code x = 0} -- where the oscillating term
 * is absent and the integrand settles -- about one. Both callers in this
 * package sit in the cheap corner; a caller who does not should know why.
 * <p>
 * https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)
 *
 * @since 1.5.3
 */
public final class WeightedChiSquare {

    /**
     * The accuracy the quadrature is asked for, per panel.
     * <p>
     * It is an <em>absolute</em> tolerance, and both tails read the integral as
     * {@code 1/2 - I/pi} and {@code 1/2 + I/pi}, so it also bounds how far down
     * a tail keeps its digits. Asking for less does not help: measured against
     * a closed form, tightening it from {@code 1e-11} to {@code 1e-16} changed
     * nothing at all, because what limits the answer is the subtraction and not
     * the rule.
     */
    public static final double QUADRATURE_TOLERANCE = 1.0e-13;

    /** The rule the panels are integrated with. */
    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;

    /**
     * How deep a single panel may be bisected before its value is taken as it
     * stands.
     * <p>
     * The rule halves its tolerance at every level, so by level {@code d} it is
     * asking for {@code tolerance / 2^d} on a panel -- a demand that cannot be
     * met in {@code double} arithmetic once the panel itself is small, and the
     * work spent trying is spent failing. Measured on the worst shape this
     * class meets, eight weights with a shift: depth {@code 14} costs 1513545
     * integrand evaluations, depth {@code 10} costs 247605 and lands within
     * {@code 2.5e-14} of it, and depth {@code 6} costs 25395 and lands within
     * {@code 2.7e-9}. Neither of the two callers in this package is affected at
     * all -- their panels converge long before the cap.
     */
    private static final int MAX_DEPTH = 10;

    /** The size of the tail of the integral that is dropped rather than integrated. */
    private static final double TRUNCATION_TOLERANCE = 1.0e-15;

    /** Bisections of the first panel taken before the error estimate is consulted. */
    private static final int FORCED_SUBDIVISIONS = 3;

    /**
     * {@code P[sum lambda_j z_j^2 <= x]}.
     * <p>
     * See {@link #barF} for how far down a small probability keeps its digits:
     * the two tails are the same integral read from either side of one half,
     * and the limit that puts on them is the same.
     *
     * @param lambda
     *            the weights, of either sign and at least one of them
     * @param x
     *            the value the sum is compared against
     * @return {@code P[Q <= x]}
     * @throws IllegalArgumentException
     *             if {@code lambda} is {@code null} or empty, holds a value
     *             that is not finite, or if {@code x} is not a number
     */
    public static double cdf(double[] lambda, double x) {
        return tail(lambda, x, false, QUADRATURE_TOLERANCE);
    }

    /**
     * {@code P[sum lambda_j z_j^2 <= x]}, computed to {@code tolerance}.
     * <p>
     * See {@link #barF(double[], double, double)} for what a looser tolerance
     * is worth.
     *
     * @param lambda
     *            the weights, of either sign and at least one of them
     * @param x
     *            the value the sum is compared against
     * @param tolerance
     *            the absolute accuracy asked of each panel, strictly positive
     * @return {@code P[Q <= x]}
     * @throws IllegalArgumentException
     *             if {@code lambda} is {@code null} or empty, holds a value that
     *             is not finite, if {@code x} is not a number, or if
     *             {@code tolerance} is not strictly positive
     */
    public static double cdf(double[] lambda, double x, double tolerance) {
        return tail(lambda, x, false, tolerance);
    }

    /**
     * {@code P[sum lambda_j z_j^2 >= x]}.
     * <p>
     * <b>Both tails are a number near one half plus or minus the same
     * integral</b>, so a small probability here is a difference of two numbers
     * near one half -- the one thing the rest of this package is written to
     * avoid, and here it is Imhof's formula rather than a choice. Measured
     * against an exact closed form, the relative error was below {@code 1e-11}
     * at {@code 1e-3}, below {@code 1e-9} at {@code 1e-4} and below
     * {@code 1e-5} at {@code 1e-6}; from about {@code 1e-8} downwards only the
     * order of magnitude survives, and further down the answer is clamped to
     * zero. Every conventional decision level is far above that.
     *
     * @param lambda
     *            the weights, of either sign and at least one of them
     * @param x
     *            the value the sum is compared against
     * @return {@code P[Q >= x]}
     * @throws IllegalArgumentException
     *             if {@code lambda} is {@code null} or empty, holds a value
     *             that is not finite, or if {@code x} is not a number
     */
    public static double barF(double[] lambda, double x) {
        return tail(lambda, x, true, QUADRATURE_TOLERANCE);
    }

    /**
     * {@code P[sum lambda_j z_j^2 >= x]}, computed to {@code tolerance}.
     * <p>
     * The default is tight enough that nothing is gained by tightening it
     * further, but plenty is gained by loosening it when the caller does not
     * need it. The integrand costs one arc tangent and one logarithm per
     * weight, so it is the number of evaluations that decides the running time:
     * measured on thirty weights spanning nine hundred to one, asking for
     * {@code 1e-9} instead of {@code 1e-13} cut 7380 evaluations to 2580 and the
     * time from 40 milliseconds to 10, and moved the answer by {@code 4.9e-11}.
     *
     * @param lambda
     *            the weights, of either sign and at least one of them
     * @param x
     *            the value the sum is compared against
     * @param tolerance
     *            the absolute accuracy asked of each panel, strictly positive
     * @return {@code P[Q >= x]}
     * @throws IllegalArgumentException
     *             if {@code lambda} is {@code null} or empty, holds a value that
     *             is not finite, if {@code x} is not a number, or if
     *             {@code tolerance} is not strictly positive
     */
    public static double barF(double[] lambda, double x, double tolerance) {
        return tail(lambda, x, true, tolerance);
    }

    /**
     * Either tail, from Imhof's inversion: {@code P[Q >= x]} is
     * {@code 1/2 + I/pi} with
     * <p>
     * {@code I = integral_0^inf sin(theta(v)) / (v rho(v)) dv},
     * {@code theta(v) = (1/2) sum atan(lambda_j v) - x v / 2} and
     * {@code rho(v) = prod (1 + lambda_j^2 v^2)^(1/4)}.
     * <p>
     * The integral runs over a half line and is cut where Imhof's own bound on
     * the remainder falls below {@link #TRUNCATION_TOLERANCE}. That point can
     * be very far out when a weight is close to zero, so the range is covered
     * by panels that double in width rather than by one interval.
     */
    private static double tail(double[] weights, double x, boolean upper, double tolerance) {
        if (weights == null || weights.length == 0) {
            throw new IllegalArgumentException("lambda must not be null or empty");
        }
        if (Double.isNaN(x)) {
            throw new IllegalArgumentException("x must not be NaN");
        }
        if (!(tolerance > 0.0)) {
            throw new IllegalArgumentException("tolerance must be strictly positive : " + tolerance);
        }
        final double[] lambda = nonZero(weights);
        if (lambda.length == 0) {
            // every weight is zero, so Q is identically zero
            if (upper) {
                return x <= 0.0 ? 1.0 : 0.0;
            }
            return x >= 0.0 ? 1.0 : 0.0;
        }
        if (lambda.length == 1) {
            // one term, so this is a chi-square with one degree of freedom
            // read through a scale that may be negative. Both tails of that
            // are error functions, which beats an inversion of a
            // characteristic function that has almost nothing left to invert
            return oneTerm(lambda[0], x, upper);
        }

        DFunction integrand = new DFunction() {
            @Override
            public double apply(double v) {
                if (v <= 0.0) {
                    // the limit at the origin, where sin(theta)/v is not yet
                    // an expression a machine can evaluate
                    double sum = 0.0;
                    for (int i = 0; i < lambda.length; i++) {
                        sum += lambda[i];
                    }
                    return 0.5 * (sum - x);
                }
                double theta = 0.0;
                for (int i = 0; i < lambda.length; i++) {
                    theta += Math.atan(lambda[i] * v);
                }
                return Math.sin(0.5 * theta - 0.5 * x * v) * Math.exp(-Math.log(v) - logRho(lambda, v));
            }
        };

        // one panel for the interval the integrand is O(1) on, then panels
        // that double in width. A single panel over the whole range would not
        // do: the rule would place its nodes so far apart that it never sees
        // the part that carries the integral, and its own error estimate would
        // agree that there was nothing there
        double integral = AdaptiveGaussKronrod.integrate1DAdaptive(RULE, integrand, 0.0, 1.0,
                tolerance, MAX_DEPTH, FORCED_SUBDIVISIONS);
        double upperLimit = truncationLimit(lambda);
        double from = 1.0;
        while (from < upperLimit) {
            double to = 2.0 * from;
            integral += AdaptiveGaussKronrod.integrate1DAdaptive(RULE, integrand, from, to,
                    tolerance, MAX_DEPTH, 0);
            from = to;
        }
        return clamp(upper ? 0.5 + integral / Math.PI : 0.5 - integral / Math.PI);
    }

    /**
     * {@code P[lambda z^2 <= x]} or its complement, from the error function.
     * <p>
     * {@code P[z^2 <= q]} is {@code erf(sqrt(q/2))} and {@code P[z^2 >= q]} is
     * {@code erfc(sqrt(q/2))}, so neither tail is formed as one minus the
     * other. A negative weight reflects the sum, which exchanges them.
     */
    private static double oneTerm(double lambda, double x, boolean upper) {
        double q = x / lambda;
        boolean atMost = upper == (lambda < 0.0);
        if (q < 0.0) {
            // the whole support lies on one side of x. At q = 0 the error
            // functions already say 0 and 1, so only a strictly negative
            // threshold needs answering here
            return atMost ? 0.0 : 1.0;
        }
        double root = Math.sqrt(0.5 * q);
        return clamp(atMost ? ProbabilityFuncs.errorFunction(root)
                : ProbabilityFuncs.errorFunctionComplemented(root));
    }

    /**
     * The logarithm of Imhof's {@code rho(v) = prod (1 + lambda_j^2 v^2)^(1/4)}.
     * <p>
     * Written as a sum because the product passes the largest {@code double}
     * for a few dozen weights, and each term avoids {@code log1p} of a square
     * that would overflow before its logarithm would.
     */
    private static double logRho(double[] lambda, double v) {
        double sum = 0.0;
        for (int i = 0; i < lambda.length; i++) {
            double lv = lambda[i] * v;
            double magnitude = Math.abs(lv);
            sum += magnitude > 1.0e8 ? 2.0 * Math.log(magnitude) : Math.log1p(lv * lv);
        }
        return 0.25 * sum;
    }

    /**
     * Where the rest of the integral is smaller than
     * {@link #TRUNCATION_TOLERANCE}.
     * <p>
     * Beyond {@code U} each factor of {@code rho} grows at least as
     * {@code (v/U)^(1/2)}, so the remainder is bounded by
     * {@code 2 / (pi m rho(U))} -- Imhof's own bound, which is why the half
     * line does not have to be substituted away.
     */
    private static double truncationLimit(double[] lambda) {
        double u = 1.0;
        for (int doublings = 0; doublings < 200; doublings++) {
            if (2.0 / (Math.PI * lambda.length * Math.exp(logRho(lambda, u))) < TRUNCATION_TOLERANCE) {
                return u;
            }
            u *= 2.0;
        }
        return u;
    }

    /** The weights with the exact zeros left out, which change neither tail. */
    private static double[] nonZero(double[] weights) {
        int kept = 0;
        for (int i = 0; i < weights.length; i++) {
            if (!isFinite(weights[i])) {
                throw new IllegalArgumentException("lambda[" + i + "] is not finite : " + weights[i]);
            }
            if (weights[i] != 0.0) {
                kept++;
            }
        }
        double[] lambda = new double[kept];
        int at = 0;
        for (int i = 0; i < weights.length; i++) {
            if (weights[i] != 0.0) {
                lambda[at++] = weights[i];
            }
        }
        return lambda;
    }

    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }

    private static double clamp(double p) {
        if (p <= 0.0) {
            return 0.0;
        }
        if (p >= 1.0) {
            return 1.0;
        }
        return p;
    }

    private WeightedChiSquare() {
        throw new AssertionError();
    }
}
