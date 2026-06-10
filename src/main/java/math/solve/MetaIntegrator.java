package math.solve;

import java.util.stream.IntStream;
import math.fun.DFunction;
import math.fun.DBiFunction;
import math.fun.DTriFunction;

public class MetaIntegrator {

    // 1D Configuration (N=64 -> 65 Points)
    private static final int CC_N = 64;
    private static final int CC_POINTS = CC_N + 1;
    private static final double[] CC_XI = new double[CC_POINTS];
    private static final double[] CC_W = new double[CC_POINTS];

    // 2D Configuration (N=32 -> 33 Points)
    private static final int CC_N2D = 32;
    private static final int CC_P_2D = CC_N2D + 1;
    private static final double[] CC_XI_2D = new double[CC_P_2D];
    private static final double[] CC_W_2D = new double[CC_P_2D];

    // 3D Configuration (N=16 -> 17 Points)
    private static final int CC_N3D = 16;
    private static final int CC_P_3D = CC_N3D + 1;
    private static final double[] CC_XI_3D = new double[CC_P_3D];
    private static final double[] CC_W_3D = new double[CC_P_3D];

    static {
        for (int i = 0; i <= CC_N; i++) {
            // Symmetric Tschebyscheff-Nodes on [-1, 1]
            CC_XI[i] = Math.cos(i * Math.PI / CC_N);

            // Canonical Clenshaw-Curtis formula
            double s = 1.0;
            for (int k = 1; k < CC_N / 2; k++) {
                s -= (2.0 / (4.0 * k * k - 1.0)) * Math.cos(2.0 * k * i * Math.PI / CC_N);
            }
            s -= (1.0 / (CC_N * CC_N - 1.0)) * Math.cos(i * Math.PI);

            double fac = (i == 0 || i == CC_N) ? 0.5 : 1.0;
            CC_W[i] = (2.0 / CC_N) * s * fac;
        }

        // Initialize 2D Matrix (N=32)
        for (int i = 0; i <= CC_N2D; i++) {
            CC_XI_2D[i] = Math.cos(i * Math.PI / CC_N2D);
            double s = 1.0;
            for (int k = 1; k < CC_N2D / 2; k++) {
                s -= (2.0 / (4.0 * k * k - 1.0)) * Math.cos(2.0 * k * i * Math.PI / CC_N2D);
            }
            s -= (1.0 / (CC_N2D * CC_N2D - 1.0)) * Math.cos(i * Math.PI);
            CC_W_2D[i] = (2.0 / CC_N2D) * s * ((i == 0 || i == CC_N2D) ? 0.5 : 1.0);
        }

        // Initialize 3D Matrix (N=16)
        for (int i = 0; i <= CC_N3D; i++) {
            CC_XI_3D[i] = Math.cos(i * Math.PI / CC_N3D);
            double s = 1.0;
            for (int k = 1; k < CC_N3D / 2; k++) {
                s -= (2.0 / (4.0 * k * k - 1.0)) * Math.cos(2.0 * k * i * Math.PI / CC_N3D);
            }
            s -= (1.0 / (CC_N3D * CC_N3D - 1.0)) * Math.cos(i * Math.PI);
            CC_W_3D[i] = (2.0 / CC_N3D) * s * ((i == 0 || i == CC_N3D) ? 0.5 : 1.0);
        }
    }

    // --- CORE CLENSHAW-CURTIS INTEGRATORS ---

    private static double performClenshawCurtis1D(DFunction f, double a, double b) {
        double c = (b - a) / 2.0; double d = (b + a) / 2.0;
        double sum = 0.0;
        for (int i = 0; i < CC_POINTS; i++) {
            sum += CC_W[i] * f.apply(c * CC_XI[i] + d);
        }
        return sum * c;
    }

    private static double performClenshawCurtis2D(DBiFunction f, double ax, double bx, double ay, double by) {
        double cx = (bx - ax) / 2.0; double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0; double dy = (by + ay) / 2.0;

        // Parallelized X-axis stream to fully engage multicore performance
        return cx * cy * IntStream.range(0, CC_P_2D).parallel().mapToDouble(i -> {
            double x = cx * CC_XI_2D[i] + dx;
            double sumY = 0.0;
            for (int j = 0; j < CC_P_2D; j++) {
                sumY += CC_W_2D[i] * CC_W_2D[j] * f.apply(x, cy * CC_XI_2D[j] + dy);
            }
            return sumY;
        }).sum();
    }

    private static double performClenshawCurtis3D(DTriFunction f, double ax, double bx, double ay, double by, double az, double bz) {
        double cx = (bx - ax) / 2.0; double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0; double dy = (by + ay) / 2.0;
        double cz = (bz - az) / 2.0; double dz = (bz + az) / 2.0;

        return cx * cy * cz * IntStream.range(0, CC_P_3D).parallel().mapToDouble(i -> {
            double x = cx * CC_XI_3D[i] + dx;
            double sumYZ = 0.0;
            for (int j = 0; j < CC_P_3D; j++) {
                double y = cy * CC_XI_3D[j] + dy;
                for (int k = 0; k < CC_P_3D; k++) {
                    sumYZ += CC_W_3D[i] * CC_W_3D[j] * CC_W_3D[k] * f.apply(x, y, cz * CC_XI_3D[k] + dz);
                }
            }
            return sumYZ;
        }).sum();
    }

    // ==========================================
    // SMART INTEGRATOR 1D
    // ==========================================
    public static double integrate1DSmart(AdaptiveGaussKronrod.G7_K15 setup, DFunction f, 
                                           double a, double b, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate1D(setup, f, a, b);

        if (firstStep.approximatedErrorEstimate < epsTol) return firstStep.value; // case 1: Smooth

        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate1D(setup, x -> Math.abs(f.apply(x)), a, b);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Aliasing triggers an oscillation switch ONLY if the values genuinely cancel each other out
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value)) && (oscillationIndex < 0.1);

        // case 2: Genuine Oscillation (Low oscillation index OR heavy aliasing noise)
        if ((oscillationIndex < 0.05 || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) { // case 2: Oscillation
            System.out.println("[Meta1D] Oscillation detected (Index: " + oscillationIndex + ", Noise: " + isAliasedOscillation + ") -> Switch to Clenshaw-Curtis.");
            return performClenshawCurtis1D(f, a, b);
        }

        return AdaptiveGaussKronrod.integrate1DAdaptive(setup, f, a, b, epsTol, 20); // case 3: Subdivision
    }

    // ==========================================
    // SMART INTEGRATOR 2D
    // ==========================================
    public static double _integrate2DSmart_(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f, 
                                           double ax, double bx, double ay, double by, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate2DParallel(setup, f, ax, bx, ay, by);

        if (firstStep.approximatedErrorEstimate < epsTol) return firstStep.value;

        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate2DParallel(setup, (x, y) -> Math.abs(f.apply(x, y)), ax, bx, ay, by);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Aliasing triggers an oscillation switch ONLY if the values genuinely cancel each other out
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value)) && (oscillationIndex < 0.1);

        // case 2: Genuine Oscillation (Low oscillation index OR heavy aliasing noise)
        if ((oscillationIndex < 0.05 || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) {
            System.out.println("[Meta2D] Oscillation detected (Index: " + oscillationIndex + ", Noise: " + isAliasedOscillation + ") -> Switch to Clenshaw-Curtis.");
            return performClenshawCurtis2D(f, ax, bx, ay, by);
        }

        return AdaptiveGaussKronrod.integrate2DAdaptive(setup, f, ax, bx, ay, by, epsTol, 10);
    }

    public static double integrate2DSmart(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f, 
            double ax, double bx, double ay, double by, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate2DParallel(setup, f, ax, bx, ay, by);

        if (firstStep.approximatedErrorEstimate < epsTol) return firstStep.value;

        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate2DParallel(setup, (x, y) -> Math.abs(f.apply(x, y)), ax, bx, ay, by);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Aliasing criteria: high error combined with clean cancellation signature
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value)) && (oscillationIndex < 0.1);

        if (((oscillationIndex < 0.05) || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) {
            System.out.println("[Meta2D] Oscillation detected -> Route to efficient Clenshaw-Curtis.");
            return performClenshawCurtis2D(f, ax, bx, ay, by);
        }

        return AdaptiveGaussKronrod.integrate2DAdaptive(setup, f, ax, bx, ay, by, epsTol, 10);
    }

    // ==========================================
    // SMART INTEGRATOR 3D
    // ==========================================

    public static double integrate3DSmart(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f, double ax, double bx,
            double ay, double by, double az, double bz, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate3DParallel(setup, f, ax, bx, ay,
                by, az, bz);

        if (firstStep.approximatedErrorEstimate < epsTol)
            return firstStep.value;

        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate3DParallel(setup,
                (x, y, z) -> Math.abs(f.apply(x, y, z)), ax, bx, ay, by, az, bz);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Block aliasing triggers on sharp 3D spikes by enforcing the cancellation index check (< 0.1)
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value))
                && (oscillationIndex < 0.1);

        if (((oscillationIndex < 0.05) || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) {
            System.out.println("[Meta3D] Oscillation detected -> Route to efficient Clenshaw-Curtis.");
            return performClenshawCurtis3D(f, ax, bx, ay, by, az, bz);
        }

        return AdaptiveGaussKronrod.integrate3DAdaptive(setup, f, ax, bx, ay, by, az, bz, epsTol, 8);
    }
}
