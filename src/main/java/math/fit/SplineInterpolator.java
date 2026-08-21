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
package math.fit;

/**
 * Natural cubic spline interpolation from commons-math3.
 * <p>
 * The interpolant is twice continuously differentiable and its second
 * derivative vanishes at both ends, which is what "natural" names. That makes
 * it the smoother and, between the knots, the more accurate of the two schemes
 * in this package on smooth data, but it does not preserve the shape of the
 * data: over a monotone step it was measured to leave the range of the values
 * by 11% in both directions. Where that matters, use
 * {@link KrugerInterpolator}.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Spline_interpolation">Spline
 * interpolation</a>
 */
public final class SplineInterpolator {

    private final CubicSpline spline;

    /**
     * Interpolates the given points.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     */
    public SplineInterpolator(double[] points, double[] values) {
        spline = interpolate(points, values);
    }

    /**
     * The value of the interpolant at {@code point}.
     *
     * @param point
     *            the abscissa, within the range of the knots
     * @return the interpolated value
     */
    public double value(double point) {
        return spline.value(point);
    }

    /**
     * The interpolant itself, which also has the derivative and the integral.
     *
     * @return the piecewise cubic this interpolator computed
     * @since 1.5.2
     */
    public CubicSpline spline() {
        return spline;
    }

    /**
     * Interpolates the given points by a natural cubic spline.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     * @return the interpolating spline
     * @since 1.5.2
     */
    public static CubicSpline interpolate(double[] points, double[] values) {
        CubicSpline.checkNodes(points, values);
        return new CubicSpline(points, polynomials(points, values));
    }

    private static double[][] polynomials(double points[], double values[]) {
        // Number of intervals. The number of data points is n + 1
        final int n = points.length - 1;

        // Differences between knot points
        double h[] = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = points[i + 1] - points[i];
        }

        double mu[] = new double[n];
        double z[] = new double[n + 1];
        mu[0] = 0.0;
        z[0] = 0.0;
        double g = 0.0;
        //@formatter:off
        for (int i = 1; i < n; i++) {
            g = 2.0 * (points[i + 1]  - points[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / g;
            z[i] = (3.0 * (values[i + 1] * h[i - 1] - values[i] * (points[i + 1] - points[i - 1]) + values[i - 1] * h[i]) /
                    (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / g;
        }
        //@formatter:on

        // cubic spline coefficients -- b is linear, c quadratic, d is cubic
        // (original values are the constants)
        double b[] = new double[n];
        double c[] = new double[n + 1];
        double d[] = new double[n];

        z[n] = 0.0;
        c[n] = 0.0;

        //@formatter:off
        for (int j = n - 1; j >= 0; j--) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (values[j + 1] - values[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        //@formatter:on

        double[][] polynomials = new double[n][];
        for (int i = 0; i < n; i++) {
            double coefficients[] = new double[4];
            coefficients[0] = values[i];
            coefficients[1] = b[i];
            coefficients[2] = c[i];
            coefficients[3] = d[i];
            polynomials[i] = coefficients;
        }

        return polynomials;
    }
}
