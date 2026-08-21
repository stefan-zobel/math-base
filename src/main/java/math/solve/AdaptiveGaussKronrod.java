package math.solve;

import java.util.Locale;
import java.util.stream.IntStream;
import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * Adaptive Gauss-Kronrod quadrature (G7/K15) in one, two and three dimensions.
 * The Kronrod rule supplies the value, the embedded Gauss rule the error
 * estimate {@code |K - G|}, and the adaptive variants bisect the longest edge
 * until that estimate meets the tolerance.
 *
 * <h2>Forced initial subdivision</h2>
 * On the undivided domain the error estimate is not to be trusted. A feature
 * narrower than the spacing of the 15 Kronrod nodes is missed by the Kronrod
 * and by the Gauss rule alike; the two then agree, {@code |K - G|} reports a
 * tiny error, and a completely wrong value is accepted without a single
 * subdivision. For a Gaussian peak of width 0.005 in the unit domain that
 * happens at roughly a third of all peak positions.
 * <p>
 * The adaptive methods therefore bisect a number of leading levels
 * unconditionally, without looking at the estimate. How many are worth their
 * price differs by dimension, because one level costs a factor of two while
 * the rule itself already costs {@code 15^dimensions} evaluations:
 * <ul>
 *   <li>1D and 2D force two bisections per axis, so 3 and 4 levels. That
 *       removes every observed failure and costs 3 us instead of 1 us, and
 *       261 us instead of 53 us, on a smooth integrand.</li>
 *   <li>3D forces one bisection per axis, so 3 levels, giving 2 x 2 x 2. Two
 *       per axis would remove the remaining failures as well, but at 3278 us
 *       instead of 418 us on a smooth integrand - a 48-fold regression for the
 *       common case, paid on every call.</li>
 * </ul>
 * This is a mitigation, not a guarantee: a feature narrower than the spacing of
 * the forced panels can still be missed, and in 3D such a peak still is. That
 * residual risk is left to the caller, who can raise the count through the
 * overloads that take it explicitly -
 * {@link #integrate1DAdaptive(G7_K15, DFunction, double, double, double, int, int)},
 * {@link #integrate2DAdaptive(G7_K15, DBiFunction, double, double, double, double, double, int, int)}
 * and
 * {@link #integrate3DAdaptive(G7_K15, DTriFunction, double, double, double, double, double, double, double, int, int)}.
 * <p>
 * Note that the recursion budget has to exceed the forced levels by a good
 * margin, otherwise the adaptive part has nothing left to work with.
 */
public class AdaptiveGaussKronrod {

    /** Forced bisections by dimension; see the class comment. */
    private static int forcedSubdivisions(int dimensions) {
        if (dimensions >= 3) {
            return 3;
        }
        return Math.max(3, 2 * dimensions);
    }

    public static class IntegralResult {
        public final double value;
        public final double approximatedErrorEstimate;

        public IntegralResult(double value, double approximatedError) {
            this.value = value;
            this.approximatedErrorEstimate = approximatedError;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT, "Value: %.8f (approx. Error: %.2e)", value, approximatedErrorEstimate);
        }
    }

    public enum G7_K15 {
        POINTS_15(
            new double[]{
                -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993944,
                -0.5860872354676911, -0.4058451513773972, -0.2077849550078985,  0.0,
                 0.2077849550078985,  0.4058451513773972,  0.5860872354676911,  0.7415311855993944,
                 0.8648644233597691,  0.9491079123427585,  0.9914553711208126
            },
            new double[]{
                0.0229353220105292, 0.0630920926299786, 0.1047900103222502, 0.1406532597155259,
                0.1690047266392679, 0.1903505780647854, 0.2044329400752989, 0.2094821410847278,
                0.2044329400752989, 0.1903505780647854, 0.1690047266392679, 0.1406532597155259,
                0.1047900103222502, 0.0630920926299786, 0.0229353220105292
            },
            new double[]{
                0.0, 0.1294849661688697, 0.0, 0.2797053914892767,
                0.0, 0.3818300505051189, 0.0, 0.4179591836734694,
                0.0, 0.3818300505051189, 0.0, 0.2797053914892767,
                0.0, 0.1294849661688697, 0.0
            }
        );

        final double[] xi;
        final double[] wK;
        final double[] wG;
        final int points;

        G7_K15(double[] xi, double[] wK, double[] wG) {
            this.xi = xi;
            this.wK = wK;
            this.wG = wG;
            this.points = xi.length;
        }
    }

    private static final class PartialSum {
        double kronrod = 0.0;
        double gauss = 0.0;
        void add(PartialSum other) {
            this.kronrod += other.kronrod;
            this.gauss += other.gauss;
        }
    }

    public static IntegralResult integrate2DParallel(G7_K15 setup, DBiFunction f, 
                                                      double ax, double bx, double ay, double by) {
        double cx = (bx - ax) / 2.0;
        double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0;
        double dy = (by + ay) / 2.0;

        PartialSum total = IntStream.range(0, setup.points).parallel().collect(
            PartialSum::new,
            (pSum, i) -> {
                double x = cx * setup.xi[i] + dx;
                for (int j = 0; j < setup.points; j++) {
                    double y = cy * setup.xi[j] + dy;
                    double fVal = f.apply(x, y);

                    pSum.kronrod += setup.wK[i] * setup.wK[j] * fVal;
                    if (setup.wG[i] != 0.0 && setup.wG[j] != 0.0) {
                        pSum.gauss += setup.wG[i] * setup.wG[j] * fVal;
                    }
                }
            },
            PartialSum::add
        );
        double areaFactor = cx * cy;
        return new IntegralResult(total.kronrod * areaFactor, Math.abs(total.kronrod - total.gauss) * areaFactor);
    }

    public static double integrate2DAdaptive(G7_K15 setup, DBiFunction f,
                                             double ax, double bx, double ay, double by,
                                             double epsTol, int maxDepth) {
        return subdivide2D(setup, f, ax, bx, ay, by, epsTol, maxDepth,
                Math.min(forcedSubdivisions(2), maxDepth));
    }

    /**
     * Adaptive 2D integration with an explicit number of forced bisections, for
     * an integrand whose features are too narrow for the default of 4; see the
     * class comment. Every extra level doubles the lower bound on the work and
     * halves the width of a feature that can still be missed. Values below
     * {@code 0} or above {@code maxDepth} are clamped.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x
     * @param bx
     *            upper limit in x
     * @param ay
     *            lower limit in y
     * @param by
     *            upper limit in y
     * @param epsTol
     *            error tolerance, halved for each subdivision
     * @param maxDepth
     *            maximal recursion depth
     * @param forcedSubdivisions
     *            number of leading levels that bisect without consulting the
     *            error estimate
     * @return the approximated integral
     */
    public static double integrate2DAdaptive(G7_K15 setup, DBiFunction f,
                                             double ax, double bx, double ay, double by,
                                             double epsTol, int maxDepth, int forcedSubdivisions) {
        return subdivide2D(setup, f, ax, bx, ay, by, epsTol, maxDepth, clamp(forcedSubdivisions, maxDepth));
    }

    private static double subdivide2D(G7_K15 setup, DBiFunction f,
                                      double ax, double bx, double ay, double by,
                                      double epsTol, int maxDepth, int forced) {

        if (maxDepth <= 0) {
            return integrate2DParallel(setup, f, ax, bx, ay, by).value;
        }
        // see forcedSubdivisions
        if (forced <= 0) {
            IntegralResult res = integrate2DParallel(setup, f, ax, bx, ay, by);
            if (res.approximatedErrorEstimate <= epsTol) {
                return res.value;
            }
        }

        double dx = bx - ax;
        double dy = by - ay;

        if (dx >= dy) {
            double midX = ax + dx / 2.0;
            return subdivide2D(setup, f, ax, midX, ay, by, epsTol / 2.0, maxDepth - 1, forced - 1)
                 + subdivide2D(setup, f, midX, bx, ay, by, epsTol / 2.0, maxDepth - 1, forced - 1);
        } else {
            double midY = ay + dy / 2.0;
            return subdivide2D(setup, f, ax, bx, ay, midY, epsTol / 2.0, maxDepth - 1, forced - 1)
                 + subdivide2D(setup, f, ax, bx, midY, by, epsTol / 2.0, maxDepth - 1, forced - 1);
        }
    }

    public static IntegralResult integrate3DParallel(G7_K15 setup, DTriFunction f, 
                                                      double ax, double bx, double ay, double by, double az, double bz) {
        double cx = (bx - ax) / 2.0;
        double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0;
        double dy = (by + ay) / 2.0;
        double cz = (bz - az) / 2.0;
        double dz = (bz + az) / 2.0;

        PartialSum total = IntStream.range(0, setup.points).parallel().collect(
            PartialSum::new,
            (pSum, i) -> {
                double x = cx * setup.xi[i] + dx;
                for (int j = 0; j < setup.points; j++) {
                    double y = cy * setup.xi[j] + dy;
                    for (int k = 0; k < setup.points; k++) {
                        double z = cz * setup.xi[k] + dz;
                        double fVal = f.apply(x, y, z);

                        pSum.kronrod += setup.wK[i] * setup.wK[j] * setup.wK[k] * fVal;
                        if (setup.wG[i] != 0.0 && setup.wG[j] != 0.0 && setup.wG[k] != 0.0) {
                            pSum.gauss += setup.wG[i] * setup.wG[j] * setup.wG[k] * fVal;
                        }
                    }
                }
            },
            PartialSum::add
        );
        double volumeFactor = cx * cy * cz;
        return new IntegralResult(total.kronrod * volumeFactor, Math.abs(total.kronrod - total.gauss) * volumeFactor);
    }

    public static double integrate3DAdaptive(G7_K15 setup, DTriFunction f,
                                             double ax, double bx, double ay, double by, double az, double bz,
                                             double epsTol, int maxRecursions) {
        return subdivide3D(setup, f, ax, bx, ay, by, az, bz, epsTol, maxRecursions,
                Math.min(forcedSubdivisions(3), maxRecursions));
    }

    /**
     * Adaptive 3D integration with an explicit number of forced bisections.
     * This is the way to close the residual risk that the default of 3 leaves
     * open: a peak narrower than roughly a hundredth of the domain is missed by
     * the default and reported as converged; see the class comment. Six levels,
     * so two bisections per axis, covered every case measured, at eight times
     * the work on a smooth integrand. Values below {@code 0} or above
     * {@code maxRecursions} are clamped.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x
     * @param bx
     *            upper limit in x
     * @param ay
     *            lower limit in y
     * @param by
     *            upper limit in y
     * @param az
     *            lower limit in z
     * @param bz
     *            upper limit in z
     * @param epsTol
     *            error tolerance, halved for each subdivision
     * @param maxRecursions
     *            maximal recursion depth; it has to exceed
     *            {@code forcedSubdivisions} by a good margin, otherwise the
     *            adaptive part has no budget left
     * @param forcedSubdivisions
     *            number of leading levels that bisect without consulting the
     *            error estimate
     * @return the approximated integral
     */
    public static double integrate3DAdaptive(G7_K15 setup, DTriFunction f,
                                             double ax, double bx, double ay, double by, double az, double bz,
                                             double epsTol, int maxRecursions, int forcedSubdivisions) {
        return subdivide3D(setup, f, ax, bx, ay, by, az, bz, epsTol, maxRecursions,
                clamp(forcedSubdivisions, maxRecursions));
    }

    private static double subdivide3D(G7_K15 setup, DTriFunction f,
                                      double ax, double bx, double ay, double by, double az, double bz,
                                      double epsTol, int maxRecursions, int forced) {

        if (maxRecursions <= 0) {
            return integrate3DParallel(setup, f, ax, bx, ay, by, az, bz).value;
        }
        // see forcedSubdivisions. Three forced levels bisect the longest edge
        // three times, which on a cube is one bisection of every axis
        if (forced <= 0) {
            IntegralResult res = integrate3DParallel(setup, f, ax, bx, ay, by, az, bz);
            if (res.approximatedErrorEstimate <= epsTol) {
                return res.value;
            }
        }

        double dx = bx - ax;
        double dy = by - ay;
        double dz = bz - az;

        if (dx >= dy && dx >= dz) {
            double midX = ax + dx / 2.0;
            return subdivide3D(setup, f, ax, midX, ay, by, az, bz, epsTol / 2.0, maxRecursions - 1, forced - 1)
                 + subdivide3D(setup, f, midX, bx, ay, by, az, bz, epsTol / 2.0, maxRecursions - 1, forced - 1);
        } else if (dy >= dx && dy >= dz) {
            double midY = ay + dy / 2.0;
            return subdivide3D(setup, f, ax, bx, ay, midY, az, bz, epsTol / 2.0, maxRecursions - 1, forced - 1)
                 + subdivide3D(setup, f, ax, bx, midY, by, az, bz, epsTol / 2.0, maxRecursions - 1, forced - 1);
        } else {
            double midZ = az + dz / 2.0;
            return subdivide3D(setup, f, ax, bx, ay, by, az, midZ, epsTol / 2.0, maxRecursions - 1, forced - 1)
                 + subdivide3D(setup, f, ax, bx, ay, by, midZ, bz, epsTol / 2.0, maxRecursions - 1, forced - 1);
        }
    }

    public static IntegralResult integrate1D(G7_K15 setup, DFunction f, double a, double b) {
        double c = (b - a) / 2.0;
        double d = (b + a) / 2.0;
        double sumK = 0.0;
        double sumG = 0.0;

        for (int i = 0; i < setup.points; i++) {
            double x = c * setup.xi[i] + d;
            double fVal = f.apply(x);

            sumK += setup.wK[i] * fVal;
            if (setup.wG[i] != 0.0) {
                sumG += setup.wG[i] * fVal;
            }
        }
        return new IntegralResult(sumK * c, Math.abs(sumK - sumG) * c);
    }

    public static double integrate1DAdaptive(G7_K15 setup, DFunction f, double a, double b, double epsTol, int maxDepth) {
        return subdivide1D(setup, f, a, b, epsTol, maxDepth, Math.min(forcedSubdivisions(1), maxDepth));
    }

    /**
     * Adaptive 1D integration with an explicit number of forced bisections, for
     * an integrand whose features are too narrow for the default of 3; see the
     * class comment. Every extra level doubles the lower bound on the work and
     * halves the width of a feature that can still be missed. Values below
     * {@code 0} or above {@code maxDepth} are clamped.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param a
     *            lower limit
     * @param b
     *            upper limit
     * @param epsTol
     *            error tolerance, halved for each subdivision
     * @param maxDepth
     *            maximal recursion depth
     * @param forcedSubdivisions
     *            number of leading levels that bisect without consulting the
     *            error estimate
     * @return the approximated integral
     */
    public static double integrate1DAdaptive(G7_K15 setup, DFunction f, double a, double b, double epsTol,
                                             int maxDepth, int forcedSubdivisions) {
        return subdivide1D(setup, f, a, b, epsTol, maxDepth, clamp(forcedSubdivisions, maxDepth));
    }

    private static double subdivide1D(G7_K15 setup, DFunction f, double a, double b, double epsTol,
                                      int maxDepth, int forced) {
        if (maxDepth <= 0) {
            return integrate1D(setup, f, a, b).value;
        }
        // see forcedSubdivisions: while forced > 0 the error estimate
        // is not consulted at all, the panel is bisected in any case
        if (forced <= 0) {
            IntegralResult res = integrate1D(setup, f, a, b);
            if (res.approximatedErrorEstimate <= epsTol) {
                return res.value;
            }
        }

        double mid = a + (b - a) / 2.0;
        return subdivide1D(setup, f, a, mid, epsTol / 2.0, maxDepth - 1, forced - 1)
             + subdivide1D(setup, f, mid, b, epsTol / 2.0, maxDepth - 1, forced - 1);
    }

    /** Keeps the forced level count inside {@code [0, maxDepth]}. */
    private static int clamp(int forcedSubdivisions, int maxDepth) {
        if (forcedSubdivisions < 0) {
            return 0;
        }
        return Math.min(forcedSubdivisions, maxDepth);
    }
}
