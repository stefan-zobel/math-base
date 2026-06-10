package math.solve;

import java.util.stream.IntStream;
import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

public class AdaptiveGaussKronrod {

    public static class IntegralResult {
        public final double value;
        public final double approximatedErrorEstimate;

        public IntegralResult(double value, double approximatedError) {
            this.value = value;
            this.approximatedErrorEstimate = approximatedError;
        }

        @Override
        public String toString() {
            return String.format("Value: %.8f (approx. Error: %.2e)", value, approximatedErrorEstimate);
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

        IntegralResult res = integrate2DParallel(setup, f, ax, bx, ay, by);
        
        if (res.approximatedErrorEstimate <= epsTol || maxDepth <= 0) {
            return res.value;
        }

        double dx = bx - ax; 
        double dy = by - ay;

        if (dx >= dy) {
            double midX = ax + dx / 2.0;
            return integrate2DAdaptive(setup, f, ax, midX, ay, by, epsTol / 2.0, maxDepth - 1)
                 + integrate2DAdaptive(setup, f, midX, bx, ay, by, epsTol / 2.0, maxDepth - 1);
        } else {
            double midY = ay + dy / 2.0;
            return integrate2DAdaptive(setup, f, ax, bx, ay, midY, epsTol / 2.0, maxDepth - 1)
                 + integrate2DAdaptive(setup, f, ax, bx, midY, by, epsTol / 2.0, maxDepth - 1);
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

        IntegralResult res = integrate3DParallel(setup, f, ax, bx, ay, by, az, bz);

        if (res.approximatedErrorEstimate <= epsTol || maxRecursions <= 0) {
            return res.value;
        }

        double dx = bx - ax;
        double dy = by - ay;
        double dz = bz - az;

        if (dx >= dy && dx >= dz) {
            double midX = ax + dx / 2.0;
            return integrate3DAdaptive(setup, f, ax, midX, ay, by, az, bz, epsTol / 2.0, maxRecursions - 1)
                 + integrate3DAdaptive(setup, f, midX, bx, ay, by, az, bz, epsTol / 2.0, maxRecursions - 1);
        } else if (dy >= dx && dy >= dz) {
            double midY = ay + dy / 2.0;
            return integrate3DAdaptive(setup, f, ax, bx, ay, midY, az, bz, epsTol / 2.0, maxRecursions - 1)
                 + integrate3DAdaptive(setup, f, ax, bx, midY, by, az, bz, epsTol / 2.0, maxRecursions - 1);
        } else {
            double midZ = az + dz / 2.0;
            return integrate3DAdaptive(setup, f, ax, bx, ay, by, az, midZ, epsTol / 2.0, maxRecursions - 1)
                 + integrate3DAdaptive(setup, f, ax, bx, ay, by, midZ, bz, epsTol / 2.0, maxRecursions - 1);
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
        IntegralResult res = integrate1D(setup, f, a, b);

        if (res.approximatedErrorEstimate <= epsTol || maxDepth <= 0) {
            return res.value;
        }

        double mid = a + (b - a) / 2.0;
        return integrate1DAdaptive(setup, f, a, mid, epsTol / 2.0, maxDepth - 1)
             + integrate1DAdaptive(setup, f, mid, b, epsTol / 2.0, maxDepth - 1);
    }
}
