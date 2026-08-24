/*
 * Copyright (c) 1997 - 1998 by Visual Numerics, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this software is freely
 * granted by Visual Numerics, Inc., provided that the copyright notice
 * above and the following warranty disclaimer are preserved in human
 * readable form.
 *
 * Because this software is licensed free of charge, it is provided
 * "AS IS", with NO WARRANTY.  TO THE EXTENT PERMITTED BY LAW, VNI
 * DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
 * TO ITS PERFORMANCE, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * VNI WILL NOT BE LIABLE FOR ANY DAMAGES WHATSOEVER ARISING OUT OF THE USE
 * OF OR INABILITY TO USE THIS SOFTWARE, INCLUDING BUT NOT LIMITED TO DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, PUNITIVE, AND EXEMPLARY DAMAGES, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
 */
package math.imsl;

/**
 * Some ancient (but quick) trigonometric routines by Visual Numerics.
 */
public final class Trig {

    private static final double TANH_COEFF[] = { -0.25828756643634709, -0.11836106330053497, 0.0098694426480063976,
            -0.00083579866234458197, 7.0904321198942996e-5, -6.01642431812e-6, 5.1052419079999996e-7,
            -4.3320729076999997e-8, 3.6759990550000002e-9, -3.1192849599999999e-10, 2.6468827999999999e-11,
            -2.2460229999999999e-12, 1.9058700000000001e-13, -1.6171999999999999e-14, 1.372e-15,
            -1.1600000000000001e-16, 8.9999999999999999e-18 };

    // The point at which tanh reaches 1.0 in double precision:
    // 0.5 * ln(8 / eps), where 2 * e^(-2x) finally drops below half an ulp of
    // the largest double below 1.0. The routine used to saturate at
    // 7.9772948850000001, a single precision cutoff, where the true value is
    // still 1 - 2.36e-7. The branch below reaches 1.0 on its own from
    // 0.5 * ln(4 / eps) = 18.715 upwards, so past that this constant is really
    // a guard against Math.exp overflowing at 709.78, where the quotient of
    // two infinities would be NaN.
    private static final double TANH_CUTOFF = 19.061547465398498;

    /**
     * The square of the hyperbolic secant, {@code 1 / cosh(x)^2}. Computed
     * from {@link #cosh(double)} and accurate over the whole range; do not
     * derive it from {@code 1 - tanh(x)^2}, which cancels to exactly zero
     * from {@code |x| > 18.7} upwards, where this method still answers
     * 2.2e-16.
     *
     * @param x
     *            input value
     * @return {@code sech(x)^2} at {@code x}
     */
    public static double sech2(double x) {
        double y = cosh(x);
        return 1.0 / (y * y);
    }

    /**
     * The hyperbolic cosine of {@code x}.
     *
     * @param x
     *            input value
     * @return {@code cosh(x)}, {@code +Infinity} for an infinite argument
     *         and {@code NaN} for {@code NaN}
     */
    public static double cosh(double x) {
        if (Double.isNaN(x)) {
            return Double.NaN;
        }
        if (Double.isInfinite(x)) {
            return Double.POSITIVE_INFINITY;
        }
        double y = Math.exp(Math.abs(x));
        if (y < 94906265.62) {
            return 0.5 * (y + 1.0 / y);
        }
        return 0.5 * y;
    }

    /**
     * The hyperbolic tangent of {@code d}, accurate to about two ulp over
     * the whole range. It saturates to {@code +/-1} only where a
     * {@code double} does, so {@code 1 - tanh(d)^2} stays meaningful up to
     * {@code |d| = 18.7}; beyond that use {@link #sech2(double)}.
     *
     * @param d
     *            input value
     * @return {@code tanh(d)} in {@code [-1, 1]}, {@code NaN} for
     *         {@code NaN}
     */
    public static double tanh(double d) {
        double d2 = Math.abs(d);
        double d1;
        if (Double.isNaN(d)) {
            d1 = Double.NaN;
        } else if (d2 < 1.82501e-8) {
            d1 = d;
        } else if (d2 <= 1.0) {
            d1 = d * (1.0 + csevl(2.0 * d * d - 1.0, TANH_COEFF));
        } else if (d2 < TANH_CUTOFF) {
            d2 = Math.exp(d2);
            d1 = sign((d2 - 1.0 / d2) / (d2 + 1.0 / d2), d);
        } else {
            d1 = sign(1.0, d);
        }
        return d1;
    }

    private static double sign(double d, double d1) {
        double d2 = (d >= 0.0) ? d : -d;
        return (d1 >= 0.0) ? d2 : -d2;
    }

    private static double csevl(double d, double[] ad) {
        double d2 = 0.0;
        double d1 = 0.0;
        double d3 = 0.0;
        double d4 = 2.0 * d;
        for (int i = ad.length - 1; i >= 0; i--) {
            d3 = d2;
            d2 = d1;
            d1 = (d4 * d2 - d3) + ad[i];
        }
        return 0.5 * (d1 - d3);
    }

    private Trig() {
        throw new AssertionError();
    }
}
