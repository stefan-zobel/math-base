package math.fit;

/**
 * Akima's cubic interpolation, in the rule of the 1970 paper and in the
 * modified weighting of Matlab's {@code makima}.
 * <p>
 * The slope at a knot is a weighted average of the two secants meeting there,
 * weighted by how much the secants <em>beyond</em> them disagree, so a knot
 * whose neighbourhood is smooth gets a smooth slope and one next to an outlier
 * does not. That makes the scheme <em>local</em>: moving a single value changes
 * exactly the three segments on either side of it and leaves every other
 * segment bit for bit as it was, where the natural spline of
 * {@link SplineInterpolator} solves one system over all the knots and moves
 * every segment there is. It is the usual answer for data with a few outliers.
 * <p>
 * Where it sits between its two neighbours here, all measured: it is
 * {@code C^1} and not {@code C^2}; it overshoots the range of the data far
 * less than the natural spline but, unlike {@link KrugerInterpolator}, not
 * never -- 3.2% against 24.4% on an unevenly spaced monotone set where the
 * constrained spline gives 0.0%; and between the knots on smooth data it is
 * worse than the natural spline on coarse grids and better on fine ones
 * (5.1e-2 against 2.9e-2 at 11 knots, 9.8e-4 against 1.7e-3 at 41).
 * {@link Variant#CLASSIC} alone reproduces a quadratic exactly on a uniform
 * grid, which none of the other schemes in this package does.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Akima_spline">Akima spline</a>
 *
 * @since 1.5.2
 */
public final class AkimaInterpolator {

    /** Which weighting the slope rule uses. */
    public enum Variant {

        /**
         * The rule of the 1970 paper. Where all four secants around a knot
         * agree in pairs the weights both vanish and the slope is the average
         * of the two adjacent secants, which is right in the middle of smooth
         * data and wrong at a corner where a ramp meets a plateau: the average
         * of the incoming slope and zero tilts the curve into the flat part.
         * Measured on a ramp running into a plateau at 3, it rises to 3.074.
         */
        CLASSIC,

        /**
         * The weighting of Matlab's {@code makima} (R2019b), which adds
         * {@code |m_a + m_b| / 2} to each weight. That breaks the tie at such
         * a corner in favour of the flat side -- the same plateau stays exactly
         * flat -- at the price of no longer reproducing a quadratic exactly.
         * On data without flat stretches the two differ by a percent or so of
         * the curve, neither being the better one in general.
         */
        MODIFIED
    }

    private final CubicSpline spline;

    /**
     * Interpolates the given points by {@link Variant#CLASSIC}.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     */
    public AkimaInterpolator(double[] points, double[] values) {
        this(points, values, Variant.CLASSIC);
    }

    /**
     * Interpolates the given points.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     * @param variant
     *            which weighting the slope rule uses
     */
    public AkimaInterpolator(double[] points, double[] values, Variant variant) {
        spline = interpolate(points, values, variant);
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
     */
    public CubicSpline spline() {
        return spline;
    }

    /**
     * Interpolates the given points by {@link Variant#CLASSIC}.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     * @return the interpolating spline
     */
    public static CubicSpline interpolate(double[] points, double[] values) {
        return interpolate(points, values, Variant.CLASSIC);
    }

    /**
     * Interpolates the given points.
     *
     * @param points
     *            the abscissae, strictly increasing, at least two
     * @param values
     *            the values at those abscissae
     * @param variant
     *            which weighting the slope rule uses
     * @return the interpolating spline
     */
    public static CubicSpline interpolate(double[] points, double[] values, Variant variant) {
        if (variant == null) {
            throw new IllegalArgumentException("variant must not be null");
        }
        CubicSpline.checkNodes(points, values);
        return new CubicSpline(points, CubicSpline.hermite(points, values, slopes(points, values, variant)));
    }

    /**
     * The slope at knot {@code i} is
     * {@code (w1 m[i-1] + w2 m[i]) / (w1 + w2)} over the secants {@code m},
     * with {@code w1 = |m[i+1] - m[i]|} and {@code w2 = |m[i-1] - m[i-2]|}.
     * The two secants a knot does not have at either end are the artificial
     * ones the paper prescribes, extrapolated linearly from the two innermost.
     */
    private static double[] slopes(double[] points, double[] values, Variant variant) {
        final int n = points.length - 1;
        double[] slopes = new double[n + 1];
        if (n == 1) {
            // the rule needs two real secants and there is only one; through
            // two points the only sensible interpolant is the line
            slopes[0] = (values[1] - values[0]) / (points[1] - points[0]);
            slopes[1] = slopes[0];
            return slopes;
        }
        // m[k] lives at m[k + 2], so that the artificial secants m[-2], m[-1],
        // m[n] and m[n+1] have somewhere to go
        double[] m = new double[n + 4];
        for (int k = 0; k < n; k++) {
            m[k + 2] = (values[k + 1] - values[k]) / (points[k + 1] - points[k]);
        }
        m[1] = 2.0 * m[2] - m[3];
        m[0] = 2.0 * m[1] - m[2];
        m[n + 2] = 2.0 * m[n + 1] - m[n];
        m[n + 3] = 2.0 * m[n + 2] - m[n + 1];

        boolean modified = variant == Variant.MODIFIED;
        for (int i = 0; i <= n; i++) {
            double left = m[i + 1];
            double right = m[i + 2];
            double w1 = Math.abs(m[i + 3] - m[i + 2]);
            double w2 = Math.abs(m[i + 1] - m[i]);
            if (modified) {
                w1 += Math.abs(m[i + 3] + m[i + 2]) / 2.0;
                w2 += Math.abs(m[i + 1] + m[i]) / 2.0;
            }
            double sum = w1 + w2;
            if (sum > 0.0) {
                slopes[i] = (w1 * left + w2 * right) / sum;
            } else {
                // both neighbourhoods are straight; the rule has nothing to
                // prefer and the average is its continuous limit
                slopes[i] = (left + right) / 2.0;
            }
        }
        return slopes;
    }
}
