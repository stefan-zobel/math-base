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
 * Kruger's constrained cubic spline interpolation
 * (http://www.korf.co.uk/spline.pdf).
 * <p>
 * The knot slopes are chosen so that the interpolant does not overshoot: it
 * stays within the range of the data and reproduces a monotone data set
 * monotonically, where the natural spline of {@link SplineInterpolator} was
 * measured to leave the range by 11% over a monotone step. The price is
 * smoothness -- the second derivative jumps at the knots, so the interpolant is
 * only once continuously differentiable -- and somewhat less accuracy between
 * the knots on smooth data.
 * <p>
 */
public final class KrugerInterpolator {

    private final CubicSpline spline;

    /**
     * Interpolates the given points.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     */
    public KrugerInterpolator(double[] points, double[] values) {
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
     * Interpolates the given points by a constrained cubic spline.
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

    /**
     * The slopes of the paper: the harmonic mean of the neighbouring secants
     * where they agree in sign, zero where they do not, which is what keeps the
     * interpolant inside the data.
     */
    private static double[] slopes(double[] dx, double[] dy) {
        final int n = dx.length;
        double f1[] = new double[n + 1];
        if (n == 1) {
            // the end formulas below need an interior slope, and with two
            // points there is none; the only sensible interpolant through two
            // points is the straight line through them
            f1[0] = dy[0] / dx[0];
            f1[1] = f1[0];
            return f1;
        }
        for (int i = 1; i < n; i++) {
            double slope = dy[i - 1] * dy[i];
            if (slope > 0.0) {
                // doesn't change sign
                f1[i] = 2.0 / (dx[i] / dy[i] + dx[i - 1] / dy[i - 1]);
            } else {
                // changes sign
                f1[i] = 0.0;
            }
        }
        f1[0] = 3.0 * dy[0] / (2.0 * dx[0]) - f1[1] / 2.0;
        f1[n] = 3.0 * dy[n - 1] / (2.0 * dx[n - 1]) - f1[n - 1] / 2.0;
        return f1;
    }

    private static double[][] polynomials(double[] points, double[] values) {
        // Number of intervals. The number of data points is n + 1
        final int n = points.length - 1;
        // Differences between knot points
        double dx[] = new double[n];
        double dy[] = new double[n];
        for (int i = 0; i < n; i++) {
            dx[i] = points[i + 1] - points[i];
            dy[i] = values[i + 1] - values[i];
        }

        double[] f1 = slopes(dx, dy);

        // The piece over [x_i, x_i+1] is the cubic that matches the value and
        // the slope at both of its ends, written in the local coordinate
        // t = x - x_i. Writing it in x instead, as the paper does, cancels the
        // abscissa away: on knots spread over six orders of magnitude that
        // costs every digit of the result.
        double[][] polynomials = new double[n][];
        for (int i = 0; i < n; i++) {
            double h = dx[i];
            double s = dy[i] / h;
            double coefficients[] = new double[4];
            coefficients[0] = values[i];
            coefficients[1] = f1[i];
            coefficients[2] = (3.0 * s - 2.0 * f1[i] - f1[i + 1]) / h;
            coefficients[3] = (f1[i] + f1[i + 1] - 2.0 * s) / (h * h);
            polynomials[i] = coefficients;
        }

        return polynomials;
    }
}
