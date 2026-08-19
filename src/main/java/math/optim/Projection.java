package math.optim;

import java.util.Arrays;

/**
 * The orthogonal projection onto a closed convex set: {@link #projectInto}
 * replaces its argument by the nearest point of the set. That is everything a
 * projected gradient method needs to know about a constraint, which is why
 * this is an interface with a single method rather than a description of a
 * set.
 * <p>
 * Implementations have to satisfy three properties, because
 * {@link SpectralProjectedGradient} relies on all of them: the set is convex
 * and closed, the result is the nearest point of it in the Euclidean norm, and
 * the operation is idempotent. Idempotence, and the related property that a
 * point already in the set comes back unchanged, hold exactly for
 * {@link #box(double[], double[])} and {@link #nonNegative(int)}; for
 * {@link #simplex(int, double)}, {@link #euclideanBall(int, double)} and
 * {@link #l1Ball(int, double)} they hold up to the rounding of the shift or
 * the scaling involved.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Projections_onto_convex_sets">
 * Wikipedia projections onto convex sets</a>.
 *
 * @since 1.5.2
 */
@FunctionalInterface
public interface Projection {

    /**
     * Replaces {@code x} by the point of the set nearest to it.
     *
     * @param x
     *            the point to project, overwritten by its projection
     */
    void projectInto(double[] x);

    /**
     * The projection onto the box {@code lower <= x <= upper}, which clamps
     * each coordinate into its own interval. An infinite bound leaves the
     * coordinate unconstrained on that side, and {@code lower[i] == upper[i]}
     * is the way to hold coordinate {@code i} fixed at that value.
     *
     * @param lower
     *            the lower bounds, none of them {@code NaN}, not modified
     * @param upper
     *            the upper bounds, none of them {@code NaN}, not modified,
     *            each {@code upper[i] >= lower[i]}
     * @return the projection onto the box
     */
    static Projection box(double[] lower, double[] upper) {
        if (lower == null) {
            throw new IllegalArgumentException("lower is null");
        }
        if (upper == null) {
            throw new IllegalArgumentException("upper is null");
        }
        if (lower.length != upper.length) {
            throw new IllegalArgumentException(
                    "lower has length " + lower.length + ", upper has length " + upper.length);
        }
        Projections.checkDimension(lower.length);
        for (int i = 0; i < lower.length; i++) {
            if (Double.isNaN(lower[i])) {
                throw new IllegalArgumentException("lower[" + i + "] is NaN");
            }
            if (Double.isNaN(upper[i])) {
                throw new IllegalArgumentException("upper[" + i + "] is NaN");
            }
            if (!(lower[i] <= upper[i])) {
                throw new IllegalArgumentException(
                        "lower[" + i + "] = " + lower[i] + " exceeds upper[" + i + "] = " + upper[i]);
            }
        }
        final double[] low = lower.clone();
        final double[] high = upper.clone();
        return x -> {
            Projections.checkLength(x, low.length);
            for (int i = 0; i < x.length; i++) {
                double v = x[i];
                if (v < low[i]) {
                    x[i] = low[i];
                } else if (v > high[i]) {
                    x[i] = high[i];
                }
            }
        };
    }

    /**
     * The projection onto the non-negative orthant {@code x >= 0}, the special
     * case of {@link #box(double[], double[])} that non-negative least squares
     * and non-negative matrix factorization ask for.
     *
     * @param dimension
     *            the length of the vectors to project, at least one
     * @return the projection onto the non-negative orthant
     */
    static Projection nonNegative(int dimension) {
        Projections.checkDimension(dimension);
        final int n = dimension;
        return x -> {
            Projections.checkLength(x, n);
            for (int i = 0; i < x.length; i++) {
                if (x[i] < 0.0) {
                    x[i] = 0.0;
                }
            }
        };
    }

    /**
     * The projection onto the simplex of all non-negative vectors whose
     * coordinates add up to {@code sum}, that is, onto the set of mixture
     * weights. This is not a bound constraint, and no purely bound-constrained
     * method can express it. Sorting dominates the cost, so this projection is
     * {@code O(n log n)} where the others are {@code O(n)}.
     *
     * @param dimension
     *            the length of the vectors to project, at least one
     * @param sum
     *            the value the coordinates have to add up to, finite and
     *            positive
     * @return the projection onto the simplex
     */
    static Projection simplex(int dimension, double sum) {
        Projections.checkDimension(dimension);
        Projections.checkPositive(sum, "sum");
        final int n = dimension;
        final double target = sum;
        return x -> {
            Projections.checkLength(x, n);
            Projections.shrinkToSum(x, target, x);
        };
    }

    /**
     * The projection onto the Euclidean ball of the given radius around the
     * origin, which leaves a point inside the ball alone and scales any other
     * one back onto the surface.
     *
     * @param dimension
     *            the length of the vectors to project, at least one
     * @param radius
     *            the radius of the ball, finite and positive
     * @return the projection onto the Euclidean ball
     */
    static Projection euclideanBall(int dimension, double radius) {
        Projections.checkDimension(dimension);
        Projections.checkPositive(radius, "radius");
        final int n = dimension;
        final double r = radius;
        return x -> {
            Projections.checkLength(x, n);
            double norm = Projections.scaledTwoNorm(x);
            if (norm > r) {
                double scale = r / norm;
                for (int i = 0; i < x.length; i++) {
                    x[i] *= scale;
                }
            }
        };
    }

    /**
     * The projection onto the {@code l1} ball of the given radius around the
     * origin. A point inside the ball is left alone; any other one has all of
     * its coordinates shrunk towards zero by a common amount, which sets the
     * small ones to exactly zero. That is the sparsity of the lasso imposed as
     * a hard restriction rather than as a penalty.
     *
     * @param dimension
     *            the length of the vectors to project, at least one
     * @param radius
     *            the radius of the ball, finite and positive
     * @return the projection onto the {@code l1} ball
     */
    static Projection l1Ball(int dimension, double radius) {
        Projections.checkDimension(dimension);
        Projections.checkPositive(radius, "radius");
        final int n = dimension;
        final double r = radius;
        return x -> {
            Projections.checkLength(x, n);
            double[] magnitude = new double[x.length];
            double oneNorm = 0.0;
            for (int i = 0; i < x.length; i++) {
                magnitude[i] = Math.abs(x[i]);
                oneNorm += magnitude[i];
            }
            if (oneNorm <= r) {
                return;
            }
            Projections.shrinkToSum(magnitude, r, magnitude);
            for (int i = 0; i < x.length; i++) {
                double m = magnitude[i];
                x[i] = (m > 0.0) ? Math.copySign(m, x[i]) : 0.0;
            }
        };
    }
}

/**
 * Shared machinery for the projections that {@link Projection} hands out. Kept
 * out of the interface because Java 8 has no private interface methods.
 */
final class Projections {

    private Projections() {
        throw new AssertionError();
    }

    static void checkDimension(int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be >= 1 : " + dimension);
        }
    }

    static void checkPositive(double value, String name) {
        if (!(value > 0.0) || Double.isInfinite(value)) {
            throw new IllegalArgumentException(name + " must be finite and > 0.0 : " + value);
        }
    }

    static void checkLength(double[] x, int expected) {
        if (x == null) {
            throw new IllegalArgumentException("x is null");
        }
        if (x.length != expected) {
            throw new IllegalArgumentException("x has length " + x.length + ", expected " + expected);
        }
    }

    /**
     * Shrinks every entry of {@code values} by the one common amount that
     * makes the non-negative part of the result add up to {@code target}, and
     * writes that result into {@code out}. {@code out} may be {@code values}.
     * This is the projection onto a simplex, and with the absolute values as
     * input it is also the radial part of the projection onto an {@code l1}
     * ball.
     * <p>
     * The classical form of Held, Wolfe and Crowder subtracts an absolute
     * threshold, which cancels: at {@code values = [3e200, 4e200]} and
     * {@code target = 1} the threshold is {@code 4e200 - 1}, which rounds back
     * to {@code 4e200}, and the point comes back unshrunk and far outside the
     * set. The shift is therefore carried as a pivot plus an offset. The pivot
     * is one of the sorted values, so subtracting it is exact where it
     * matters, and the offset is built from differences of <i>adjacent</i>
     * sorted values, which never cancel against {@code target}.
     *
     * @param values
     *            the coordinates, not modified unless {@code out} aliases them
     * @param target
     *            the sum the shrunk coordinates have to reach, positive
     * @param out
     *            receives the result, same length as {@code values}
     */
    static void shrinkToSum(double[] values, double target, double[] out) {
        int n = values.length;
        double[] u = values.clone();
        Arrays.sort(u);
        // u is ascending, so u[n - 1 - k] is the (k + 1)-st largest. The
        // prefix of the k + 1 largest stays positive exactly while
        // target > sum over that prefix of (u[j] - u[k]), a sum of
        // non-negative terms built up one adjacent difference at a time
        double pivot = u[n - 1];
        double offset = -target;
        double spread = 0.0;
        for (int k = 1; k < n; k++) {
            double uk = u[n - 1 - k];
            spread += k * (u[n - k] - uk);
            if (!(target > spread)) {
                break;
            }
            pivot = uk;
            offset = (spread - target) / (k + 1);
        }
        for (int i = 0; i < n; i++) {
            double v = (values[i] - pivot) - offset;
            out[i] = (v > 0.0) ? v : 0.0;
        }
    }

    /**
     * The Euclidean norm in the {@code dnrm2} manner, one pass with a running
     * scale, so that projecting a vector of magnitude {@code 1e200} back into
     * the unit ball does not overflow on the way.
     */
    static double scaledTwoNorm(double[] x) {
        double scale = 0.0;
        double sumOfSquares = 1.0;
        for (int i = 0; i < x.length; i++) {
            double xi = x[i];
            if (xi != 0.0) {
                double a = Math.abs(xi);
                if (scale < a) {
                    double ratio = scale / a;
                    sumOfSquares = 1.0 + sumOfSquares * ratio * ratio;
                    scale = a;
                } else {
                    double ratio = a / scale;
                    sumOfSquares += ratio * ratio;
                }
            }
        }
        return scale * Math.sqrt(sumOfSquares);
    }
}
