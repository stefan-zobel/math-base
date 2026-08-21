package math.fit;

import java.util.Arrays;

import math.fun.DiffDFunction;

/**
 * A cubic spline: a piecewise cubic polynomial over a strictly increasing set
 * of knots, evaluated in the local coordinate {@code t = x - knot(i)} of the
 * segment that contains {@code x}.
 * <p>
 * Instances come from {@link SplineInterpolator#interpolate(double[], double[])}
 * or {@link KrugerInterpolator#interpolate(double[], double[])}, which differ
 * in how they choose the coefficients. This class only evaluates them, but it
 * evaluates the derivative and the integral as well as the value, all three in
 * closed form. Since it is a {@link math.fun.DFunction}, it can be handed to
 * {@code math.solve.AdaptiveGaussKronrod} or {@code math.solve.RootFinder}
 * directly.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Spline_interpolation">Spline
 * interpolation</a>
 *
 * @since 1.5.2
 */
public final class CubicSpline implements DiffDFunction {

    private final double[] knots;
    private final double[][] coefficients;
    private final double min;
    private final double max;

    /**
     * @param knots
     *            strictly increasing knots, length {@code n + 1}
     * @param coefficients
     *            {@code n} arrays of 4, each {@code {a, b, c, d}} of the cubic
     *            {@code a + b t + c t^2 + d t^3} over
     *            {@code t = x - knots[i]}
     */
    CubicSpline(double[] knots, double[][] coefficients) {
        this.knots = knots.clone();
        this.coefficients = new double[coefficients.length][];
        for (int i = 0; i < coefficients.length; i++) {
            this.coefficients[i] = coefficients[i].clone();
        }
        this.min = this.knots[0];
        this.max = this.knots[this.knots.length - 1];
    }

    /**
     * The number of polynomial pieces, one less than the number of knots.
     *
     * @return the number of segments
     */
    public int segments() {
        return coefficients.length;
    }

    /**
     * The {@code i}-th knot.
     *
     * @param i
     *            the knot index, {@code 0} to {@link #segments()} inclusive
     * @return the abscissa of that knot
     */
    public double knot(int i) {
        if (i < 0 || i >= knots.length) {
            throw new IllegalArgumentException("knot index out of range [0, " + (knots.length - 1) + "] : " + i);
        }
        return knots[i];
    }

    /**
     * The first knot, below which this spline is not defined.
     *
     * @return the lower end of the interpolation range
     */
    public double lowerBound() {
        return min;
    }

    /**
     * The last knot, above which this spline is not defined.
     *
     * @return the upper end of the interpolation range
     */
    public double upperBound() {
        return max;
    }

    /**
     * A copy of the knots.
     *
     * @return the knots, strictly increasing
     */
    public double[] knots() {
        return knots.clone();
    }

    /**
     * A copy of the coefficients {@code {a, b, c, d}} of one segment, of the
     * cubic {@code a + b t + c t^2 + d t^3} in the local coordinate
     * {@code t = x - knot(segment)}.
     *
     * @param segment
     *            the segment index, {@code 0} to {@link #segments()} exclusive
     * @return the four coefficients of that piece
     */
    public double[] coefficients(int segment) {
        if (segment < 0 || segment >= coefficients.length) {
            throw new IllegalArgumentException(
                    "segment index out of range [0, " + (coefficients.length - 1) + "] : " + segment);
        }
        return coefficients[segment].clone();
    }

    /**
     * The value of the spline at {@code x}.
     *
     * @param x
     *            the abscissa, within {@link #lowerBound()} and
     *            {@link #upperBound()}
     * @return the interpolated value
     */
    public double value(double x) {
        int i = segmentOf(x);
        double[] k = coefficients[i];
        double t = x - knots[i];
        return t * (t * (t * k[3] + k[2]) + k[1]) + k[0];
    }

    @Override
    public double apply(double x) {
        return value(x);
    }

    @Override
    public double derivativeAt(double x) {
        int i = segmentOf(x);
        double[] k = coefficients[i];
        double t = x - knots[i];
        return t * (t * 3.0 * k[3] + 2.0 * k[2]) + k[1];
    }

    /**
     * The second derivative of the spline at {@code x}. It is continuous for a
     * {@link SplineInterpolator} spline and in general not for a
     * {@link KrugerInterpolator} one.
     *
     * @param x
     *            the abscissa, within {@link #lowerBound()} and
     *            {@link #upperBound()}
     * @return the second derivative at {@code x}
     */
    public double secondDerivativeAt(double x) {
        int i = segmentOf(x);
        double[] k = coefficients[i];
        double t = x - knots[i];
        return 6.0 * k[3] * t + 2.0 * k[2];
    }

    /**
     * The integral of the spline from {@code a} to {@code b}, in closed form.
     * It is antisymmetric: {@code integrate(b, a) == -integrate(a, b)}.
     *
     * @param a
     *            the lower limit, within {@link #lowerBound()} and
     *            {@link #upperBound()}
     * @param b
     *            the upper limit, within the same range
     * @return the definite integral
     */
    public double integrate(double a, double b) {
        checkInRange(a, "a");
        checkInRange(b, "b");
        if (a == b) {
            return 0.0;
        }
        if (a > b) {
            return -integrate(b, a);
        }
        int ia = segmentOf(a);
        int ib = segmentOf(b);
        if (ia == ib) {
            return antiderivative(ia, b - knots[ia]) - antiderivative(ia, a - knots[ia]);
        }
        double sum = antiderivative(ia, knots[ia + 1] - knots[ia]) - antiderivative(ia, a - knots[ia]);
        for (int i = ia + 1; i < ib; i++) {
            sum += antiderivative(i, knots[i + 1] - knots[i]);
        }
        return sum + antiderivative(ib, b - knots[ib]);
    }

    @Override
    public String toString() {
        return "CubicSpline[" + segments() + " segments over [" + min + ", " + max + "]]";
    }

    /** The antiderivative of piece {@code i}, vanishing at its left knot. */
    private double antiderivative(int i, double t) {
        double[] k = coefficients[i];
        return t * (k[0] + t * (k[1] / 2.0 + t * (k[2] / 3.0 + t * k[3] / 4.0)));
    }

    /** The index of the piece that owns {@code x}, after validating {@code x}. */
    private int segmentOf(double x) {
        checkInRange(x, "x");
        int i = Arrays.binarySearch(knots, x);
        if (i < 0) {
            i = -i - 2;
        }
        // there are only n - 1 pieces, so the last knot belongs to the last one
        if (i >= coefficients.length) {
            i = coefficients.length - 1;
        }
        return i;
    }

    /**
     * The precondition of the constructor, checked at the boundary of the
     * package: two arrays of equal length, at least two knots, everything
     * finite, and the knots strictly increasing.
     */
    static void checkNodes(double[] points, double[] values) {
        if (points == null || values == null) {
            throw new IllegalArgumentException("null argument");
        }
        if (points.length != values.length) {
            throw new IllegalArgumentException(
                    "points.length != values.length : " + points.length + " != " + values.length);
        }
        if (points.length < 2) {
            throw new IllegalArgumentException("at least 2 points are needed, got " + points.length);
        }
        for (int i = 0; i < points.length; i++) {
            if (Double.isNaN(points[i]) || Double.isInfinite(points[i])) {
                throw new IllegalArgumentException("points[" + i + "] is not finite : " + points[i]);
            }
            if (Double.isNaN(values[i]) || Double.isInfinite(values[i])) {
                throw new IllegalArgumentException("values[" + i + "] is not finite : " + values[i]);
            }
        }
        for (int i = 1; i < points.length; i++) {
            if (!(points[i - 1] < points[i])) {
                throw new IllegalArgumentException("points must be strictly increasing, but points[" + (i - 1)
                        + "] = " + points[i - 1] + " and points[" + i + "] = " + points[i]);
            }
        }
    }

    private void checkInRange(double x, String name) {
        // NaN fails both comparisons, so it is rejected explicitly
        if (Double.isNaN(x)) {
            throw new IllegalArgumentException(name + " is NaN");
        }
        if (x < min || x > max) {
            throw new IllegalArgumentException(name + " out of range [" + min + ", " + max + "] : " + x);
        }
    }
}
