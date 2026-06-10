package math.solve;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

public class InfiniteIntegrator {

    /**
     * Extended 1D integral that automatically handles infinite bounds.
     * Supports: [-inf, +inf], [a, +inf] and [-inf, b]
     */
    public static double integrate1DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DFunction f, 
                                               double a, double b, double epsTol) {

        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);

        // CASE A: Doubly infinite [-inf, +inf]
        if (aInf && bInf) {
            // Transformation: x = t / (1 - t^2) on the interval [-1, 1]
            DFunction transformed = t -> {
                double t2 = t * t;
                double divisor = 1.0 - t2;
                // Guard against exact division by zero at the boundaries -1 and 1
                if (Math.abs(divisor) < 1e-15) return 0.0; 

                double x = t / divisor;
                double derivative = (1.0 + t2) / (divisor * divisor);
                return f.apply(x) * derivative;
            };
            // We integrate the transformed function strictly from -1 to 1
            return MetaIntegrator.integrate1DSmart(setup, transformed, -1.0, 1.0, epsTol);
        }

        // CASE B: Semi-infinite upward [a, +inf]
        if (!aInf && bInf) {
            // Transformation: x = a + t / (1 - t) on the interval [0, 1]
            DFunction transformed = t -> {
                double divisor = 1.0 - t;
                if (Math.abs(divisor) < 1e-15) return 0.0;

                double x = a + t / divisor;
                double derivative = 1.0 / (divisor * divisor);
                return f.apply(x) * derivative;
            };
            return MetaIntegrator.integrate1DSmart(setup, transformed, 0.0, 1.0, epsTol);
        }

        // CASE C: Semi-infinite downward [-inf, b]
        if (aInf && !bInf) {
            // Transformation: x = b - t / (1 - t) on the interval [0, 1]
            DFunction transformed = t -> {
                double divisor = 1.0 - t;
                if (Math.abs(divisor) < 1e-15) return 0.0;

                double x = b - t / divisor;
                double derivative = 1.0 / (divisor * divisor);
                return f.apply(x) * derivative;
            };
            return MetaIntegrator.integrate1DSmart(setup, transformed, 0.0, 1.0, epsTol);
        }

        // CASE D: Ordinary finite integral [a, b]
        return MetaIntegrator.integrate1DSmart(setup, f, a, b, epsTol);
    }
    
    // =========================================================================
    // MULTIDIMENSIONAL INFINITE INTEGRATOR (2D)
    // =========================================================================
    public static double integrate2DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f, 
                                               double ax, double bx, double ay, double by, double epsTol) {

        // 1. Analyze bounds
        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);

        // 2. Set target integration bounds for the transformed variables
        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }

        // 3. Build dynamic wrapper function
        DBiFunction transformed = (tX, tY) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (axInf && bxInf) {
                double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = tX / div; jX = (1.0 + tX * tX) / (div * div);
            } else if (!axInf && bxInf) {
                double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = ax + tX / div; jX = 1.0 / (div * div);
            } else if (axInf && !bxInf) {
                double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = bx - tX / div; jX = 1.0 / (div * div);
            }

            // Transform Y-axis
            double y = tY, jY = 1.0;
            if (ayInf && byInf) {
                double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = tY / div; jY = (1.0 + tY * tY) / (div * div);
            } else if (!ayInf && byInf) {
                double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = ay + tY / div; jY = 1.0 / (div * div);
            } else if (ayInf && !byInf) {
                double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = by - tY / div; jY = 1.0 / (div * div);
            }

            return f.apply(x, y) * jX * jY;
        };

        // Delegate to the parallel 2D meta-integrator
        return MetaIntegrator.integrate2DSmart(setup, transformed, transAx, transBx, transAy, transBy, epsTol);
    }

    // =========================================================================
    // MULTIDIMENSIONAL INFINITE INTEGRATOR (3D)
    // =========================================================================
    public static double integrate3DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f, 
                                               double ax, double bx, double ay, double by, double az, double bz, double epsTol) {

        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);
        boolean azInf = (az == Double.NEGATIVE_INFINITY); boolean bzInf = (bz == Double.POSITIVE_INFINITY);

        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;
        double transAz = (azInf && bzInf) ? -1.0 : 0.0; double transBz = (azInf && bzInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }
        if (!azInf && !bzInf) { transAz = az; transBz = bz; }

        DTriFunction transformed = (tX, tY, tZ) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (axInf && bxInf) { double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0; x = tX / div; jX = (1.0 + tX * tX) / (div * div); }
            else if (!axInf && bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = ax + tX / div; jX = 1.0 / (div * div); }
            else if (axInf && !bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = bx - tX / div; jX = 1.0 / (div * div); }

            // Transform Y-axis
            double y = tY, jY = 1.0;
            if (ayInf && byInf) { double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0; y = tY / div; jY = (1.0 + tY * tY) / (div * div); }
            else if (!ayInf && byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = ay + tY / div; jY = 1.0 / (div * div); }
            else if (ayInf && !byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = by - tY / div; jY = 1.0 / (div * div); }

            // Transform Z-axis
            double z = tZ, jZ = 1.0;
            if (azInf && bzInf) { double div = 1.0 - tZ * tZ; if (Math.abs(div) < 1e-15) return 0.0; z = tZ / div; jZ = (1.0 + tZ * tZ) / (div * div); }
            else if (!azInf && bzInf) { double div = 1.0 - tZ; if (Math.abs(div) < 1e-15) return 0.0; z = az + tZ / div; jZ = 1.0 / (div * div); }
            else if (azInf && !bzInf) { double div = 1.0 - tZ; if (Math.abs(div) < 1e-15) return 0.0; z = bz - tZ / div; jZ = 1.0 / (div * div); }

            return f.apply(x, y, z) * jX * jY * jZ;
        };

        return MetaIntegrator.integrate3DSmart(setup, transformed, transAx, transBx, transAy, transBy, transAz, transBz, epsTol);
    }
}
