/*
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
/*
 * Copyright © 1999 CERN - European Organization for Nuclear Research.
 * Permission to use, copy, modify, distribute and sell this software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation. CERN makes no representations about the suitability of this
 * software for any purpose. It is provided "as is" without expressed or
 * implied warranty.
 */
package math.cern;


/**
 * Contains fast approximations for the &Gamma; (Gamma) family of functions.
 * <p>
 * <b>Implementation:</b> High performance implementation.
 * </p>
 * This is a port of <tt>gen_fun.cpp</tt> from the <A
 * HREF="http://www.cis.tu-graz.ac.at/stat/stadl/random.html">C-RAND /
 * WIN-RAND</A> library.
 * 
 * @author wolfgang.hoschek@cern.ch
 */
public final class FastGamma {

    private static final double c0 = 9.1893853320467274e-01;
    private static final double c1 = 8.3333333333333333e-02;
    private static final double c2 = -2.7777777777777777e-03;
    private static final double c3 = 7.9365079365079365e-04;
    private static final double c4 = -5.9523809523809524e-04;
    private static final double c5 = 8.4175084175084175e-04;
    private static final double c6 = -1.9175269175269175e-03;

    /** the six term Stirling series is accurate from here up */
    private static final double STIRLING_MIN = 11.0;

    /** below this the rational approximation in {@link GammaFun} takes over */
    private static final double RATIONAL_MAX = 13.0;

    /**
     * below this gamma is the reciprocal to well past the last bit, and the
     * reciprocal is what the rational approximation would have to form
     */
    private static final double RECIPROCAL_MAX = 1.0e-305;

    /** natural logarithm of pi */
    private static final double LN_PI = 1.14472988584940017414;

    /**
     * Returns a quick approximation of the gamma function <tt>gamma(x)</tt>.
     * <p>
     * Below zero the sign has to be put back on by hand, because
     * {@link #logGamma(double)} is the logarithm of the magnitude and gamma
     * alternates sign there: negative on <tt>(-1, 0)</tt>, positive on
     * <tt>(-2, -1)</tt>, and so on. At a pole there is no value to return --
     * the two one-sided limits are <tt>+infinity</tt> and <tt>-infinity</tt>
     * -- so that answers <tt>NaN</tt>.
     *
     * @param x the value
     * @return gamma(x)
     */
    public static double gamma(final double x) {
        if (x > 0.0 || Double.isNaN(x)) {
            return Math.exp(logGamma(x));
        }
        if (x == Math.floor(x)) {
            return Double.NaN;
        }
        final double magnitude = Math.exp(logGamma(x));
        return ((((long) Math.floor(-x)) & 1L) == 0L) ? -magnitude : magnitude;
    }

    /**
     * Returns <tt>log(|gamma(x)|)</tt>, the natural logarithm of the absolute
     * value of the gamma function.
     * <p>
     * Defined on the whole real line: <tt>+infinity</tt> at the poles, which
     * are zero and the negative integers, and the logarithm of the magnitude
     * elsewhere, since gamma alternates sign below zero. <tt>NaN</tt> answers
     * <tt>NaN</tt> and <tt>-infinity</tt>, where the magnitude has no limit --
     * it runs to zero between the poles and to infinity at each of them.
     * <p>
     * Above <tt>x = 13</tt> this is the Stirling series; below it the rational
     * approximation of {@link GammaFun#lnGamma(double)}, which is both more
     * accurate and faster there than shifting the argument up and subtracting
     * the logarithm of the shift. That subtraction used to cost about five
     * ulps of absolute error everywhere below 13, which near the two zeros of
     * <tt>log(gamma)</tt> at <tt>x = 1</tt> and <tt>x = 2</tt> is a relative
     * error of <tt>1e-10</tt>; <tt>logGamma(1.0)</tt> came out as
     * <tt>1.78e-15</tt> rather than zero.
     *
     * @param x the value
     * @return log(|gamma(x)|)
     */
    public static double logGamma(double x) {
        if (x > 0.0) {
            if (x < RATIONAL_MAX) {
                if (x < RECIPROCAL_MAX) {
                    // gamma(x) = 1/x - euler + O(x), and the correction is far
                    // below the last bit of log(x) this far down. The rational
                    // approximation would divide by x to shift the argument
                    // up, and that overflows below 1 / Double.MAX_VALUE
                    return -Math.log(x);
                }
                return GammaFun.lnGamma(x);
            }
            if (x == Double.POSITIVE_INFINITY) {
                return x;
            }
            // the shift loop the Stirling branch used to carry is gone: it ran
            // below x = 11, and nothing below 13 reaches here any more
            return (x - 0.5) * Math.log(x) - x + c0 + stirlingSeries(x);
        }
        if (Double.isNaN(x)) {
            return Double.NaN;
        }
        if (x == Double.NEGATIVE_INFINITY) {
            // |gamma| runs to zero between the poles and to infinity at each
            // of them, so there is nothing to converge to
            return Double.NaN;
        }
        if (x == Math.floor(x)) {
            return Double.POSITIVE_INFINITY;
        }
        // reflection: |gamma(x)| = pi / (|sin(pi x)| * |gamma(1 - x)|) and
        // gamma(1 - x) = q * gamma(q) with q = -x, so
        // log|gamma(x)| = log(pi) - log|sin(pi x)| - log(q) - logGamma(q).
        // The three logarithms are kept apart on purpose: q * sin(pi q)
        // underflows for a q below about 1e-162, while each factor on its own
        // is still a normal number
        final double q = -x;
        if (q < RECIPROCAL_MAX) {
            // the reciprocal again, only negative. Taken through the formula
            // below it would go through sin(pi q), and pi times a q this small
            // is a deep subnormal with a handful of bits left
            return -Math.log(q);
        }
        double f = q - Math.floor(q);
        if (f > 0.5) {
            // sin(pi f) is symmetric about one half; the smaller argument
            // is the one Math.sin can reduce without losing digits
            f = 1.0 - f;
        }
        return LN_PI - Math.log(Math.sin(Math.PI * f)) - Math.log(q) - logGamma(q);
    }

    /**
     * Returns <tt>log(gamma(x + delta) / gamma(x))</tt> without ever forming
     * the two logarithms it is the difference of.
     * <p>
     * Subtracting them is hopeless once <tt>x</tt> is large, because the
     * difference is tiny beside what it is taken from: at
     * <tt>x = 5&middot;10<sup>14</sup></tt> the answer is about 17 while one
     * ulp of <tt>logGamma(x)</tt> is already 2, so
     * <tt>logGamma(x + 0.5) - logGamma(x)</tt> comes out quantized to whole
     * units -- 20.0 where 16.92 is right -- and from
     * <tt>x = 5&middot;10<sup>15</sup></tt> it is 0. That is not a weakness of
     * {@link #logGamma(double)}, which is accurate to about one ulp
     * throughout; the information simply is not in the difference of two
     * doubles.
     * <p>
     * Written out from the Stirling series, the difference is
     * <tt>delta&middot;log(x) + (x + delta - 1/2)&middot;log1p(delta/x) -
     * delta</tt> plus the difference of the two series remainders, and every
     * term of that is computed at the size of the answer rather than at the
     * size of <tt>logGamma</tt>. Below <tt>x = 11</tt> the series does not
     * apply and the plain difference is taken, which is accurate there.
     *
     * @param x the value
     * @param delta the offset
     * @return log(|gamma(x + delta)| / |gamma(x)|)
     * @since 1.5.2
     */
    public static double logGammaRatio(double x, double delta) {
        final double y = x + delta;
        if (delta == 0.0) {
            return 0.0;
        }
        if (x < STIRLING_MIN || y < STIRLING_MIN) {
            // both logarithms are of order one here, so their difference keeps
            // every digit it has
            return logGamma(y) - logGamma(x);
        }
        return delta * Math.log(x) + (y - 0.5) * Math.log1p(delta / x) - delta + stirlingSeries(y)
                - stirlingSeries(x);
    }

    /**
     * The tail of the Stirling series, <tt>B(2n) / (2n(2n-1) x^(2n-1))</tt>
     * summed to six terms, which is what carries <tt>logGamma</tt> from
     * <tt>x = 11</tt> upwards.
     *
     * @param x the value, at least <tt>11</tt>
     * @return the series remainder
     */
    private static double stirlingSeries(double x) {
        final double r = 1.0 / (x * x);
        return (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * c6))))) / x;
    }

    private FastGamma() {
    }
}
