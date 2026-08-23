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

    /** below this the Stirling series is not used and the argument is shifted up */
    private static final double STIRLING_MIN = 11.0;

    /**
     * Returns a quick approximation of the gamma function <tt>gamma(x)</tt>.
     * 
     * @param x the value
     * @return gamma(x)
     */
    public static double gamma(final double x) {
        return Math.exp(logGamma(x));
    }

    /**
     * Returns a quick approximation of <tt>log(gamma(x))</tt>.
     * 
     * @param x the value
     * @return log(gamma(x))
     */
    public static double logGamma(double x) {
        if (x <= 0.0 /* || x > 1.3e19 */) {
            return -999;
        }

        double z;
        for (z = 1.0; x < STIRLING_MIN; x++) {
            z *= x;
        }

        final double g = (x - 0.5) * Math.log(x) - x + c0 + stirlingSeries(x);
        if (z == 1.0) {
            return g;
        }
        return g - Math.log(z);
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
     * @param x the value, strictly positive
     * @param delta the offset, with <tt>x + delta</tt> strictly positive
     * @return log(gamma(x + delta) / gamma(x)), or <tt>-999</tt> if either
     *         argument of the two gamma functions is not strictly positive
     * @since 1.5.2
     */
    public static double logGammaRatio(double x, double delta) {
        final double y = x + delta;
        if (x <= 0.0 || y <= 0.0) {
            return -999;
        }
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
