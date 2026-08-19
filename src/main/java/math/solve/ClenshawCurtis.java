package math.solve;

import java.util.stream.IntStream;
import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * Clenshaw-Curtis quadrature on Chebyshev-Lobatto points, refined by doubling
 * the grid until the change from one level to the next meets the tolerance. See
 * <a href="https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature">
 * Clenshaw-Curtis quadrature</a>.
 *
 * <h2>Why the doubling</h2>
 * A fixed grid of {@code N + 1} points resolves about {@code N / 2} half waves
 * and says nothing about what it missed: an integrand oscillating faster than
 * that is aliased and a wrong value is returned without any indication.
 * Chebyshev-Lobatto grids are nested, {@code cos(i*pi/N) == cos(2i*pi/(2N))}, so
 * point {@code i} of the {@code N} grid is point {@code 2i} of the {@code 2N}
 * grid. Evaluating at {@code N} and again at {@code 2N} therefore costs one
 * further set of {@code N} evaluations and yields {@code |I_2N - I_N|} as an
 * error estimate. Because of the reuse the whole ladder costs exactly as much
 * as its final level alone - {@code 65 + 64 + 128 + ... = 2049} points for a
 * final {@code N} of 2048 - and the same identity holds per axis in two and
 * three dimensions.
 * <p>
 * {@code |I_2N - I_N|} is an estimate, not a bound. It is usually pessimistic,
 * since Clenshaw-Curtis converges faster than any fixed order for a smooth
 * integrand, but two coarse grids can agree by accident. The ladder therefore
 * starts one level below the grid a caller would have chosen as a fixed one, so
 * that the first estimate compares that grid against its predecessor.
 *
 * @since 1.5.1
 */
public final class ClenshawCurtis {

    /** Default first grid in 1D, so that the first estimate is 33 against 65 points. */
    private static final int MIN_N_1D = 32;
    /** Default last grid in 1D, 2049 evaluations. */
    private static final int MAX_N_1D = 2048;

    /** Default first grid per axis in 2D, so that the first estimate is 17 against 33 points. */
    private static final int MIN_N_2D = 16;
    /** Default last grid per axis in 2D, 257^2 = 66049 evaluations. */
    private static final int MAX_N_2D = 256;

    /** Default first grid per axis in 3D, so that the first estimate is 9 against 17 points. */
    private static final int MIN_N_3D = 8;
    /** Default last grid per axis in 3D, 65^3 = 274625 evaluations and 2.2 MB of cached values. */
    private static final int MAX_N_3D = 64;

    private static final int MIN_TABLE_N = 8;
    private static final int MAX_TABLE_N = 2048;

    /** Node and weight tables by level, index {@code log2(N) - 3}. */
    private static final Table[] TABLES = new Table[10];

    /** The result of one integration, with the estimate that decided when to stop. */
    public static final class IntegralResult {

        /** The approximated integral, always taken from the finest grid computed. */
        public final double value;
        /** The change from the second finest to the finest grid. */
        public final double approximatedErrorEstimate;
        /** Points per axis of the finest grid, always one more than a power of two. */
        public final int points;
        /** Whether {@link #approximatedErrorEstimate} met the tolerance. */
        public final boolean converged;

        IntegralResult(double value, double approximatedErrorEstimate, int points, boolean converged) {
            this.value = value;
            this.approximatedErrorEstimate = approximatedErrorEstimate;
            this.points = points;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format("Value: %.8f (approx. Error: %.2e, %d points/axis, converged: %b)", value,
                    approximatedErrorEstimate, points, converged);
        }
    }

    private static final class Table {

        final double[] xi;
        final double[] w;

        Table(double[] xi, double[] w) {
            this.xi = xi;
            this.w = w;
        }
    }

    private ClenshawCurtis() {
        throw new AssertionError();
    }

    // =========================================================================
    // NODES AND WEIGHTS
    // =========================================================================

    private static synchronized Table table(int n) {
        int idx = Integer.numberOfTrailingZeros(n) - 3;
        Table t = TABLES[idx];
        if (t == null) {
            t = buildTable(n);
            TABLES[idx] = t;
        }
        return t;
    }

    /**
     * Nodes {@code cos(i*pi/n)} and the canonical Clenshaw-Curtis weights for an
     * even {@code n}. The cosines are taken from a table of {@code n + 1} values
     * rather than recomputed, which turns the {@code O(n^2)} construction into
     * multiply-adds; at {@code n = 2048} that saves 2.1 million transcendental
     * calls.
     */
    private static Table buildTable(int n) {
        double[] cosTab = new double[n + 1];
        for (int m = 0; m <= n; m++) {
            cosTab[m] = Math.cos(m * Math.PI / n);
        }
        double[] xi = new double[n + 1];
        double[] w = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            xi[i] = cosTab[i];
            double s = 1.0;
            for (int k = 1; k < n / 2; k++) {
                s -= (2.0 / (4.0 * k * k - 1.0)) * cosFold(cosTab, n, 2 * k * i);
            }
            // the k = n/2 term, where cos(2*k*i*pi/n) collapses to cos(i*pi)
            s -= (1.0 / ((double) n * n - 1.0)) * (((i & 1) == 0) ? 1.0 : -1.0);
            double c = (i == 0 || i == n) ? 0.5 : 1.0;
            w[i] = (2.0 / n) * s * c;
        }
        return new Table(xi, w);
    }

    /** {@code cos(m*pi/n)} for any non-negative {@code m}, folded into the table. */
    private static double cosFold(double[] cosTab, int n, int m) {
        int r = m % (2 * n);
        if (r > n) {
            r = 2 * n - r;
        }
        return cosTab[r];
    }

    /** Rounds down to a power of two inside the range of the node tables. */
    private static int normalize(int n) {
        if (n < MIN_TABLE_N) {
            return MIN_TABLE_N;
        }
        if (n > MAX_TABLE_N) {
            return MAX_TABLE_N;
        }
        return Integer.highestOneBit(n);
    }

    /** The ladder needs at least one doubling, otherwise there is no estimate. */
    private static int normalizeMax(int maxN, int minN) {
        int hi = normalize(maxN);
        if (hi < 2 * minN) {
            hi = Math.min(2 * minN, MAX_TABLE_N);
        }
        return hi;
    }

    // =========================================================================
    // ONE DIMENSION
    // =========================================================================

    /**
     * Integrates {@code f} over {@code [a, b]}, doubling the grid from 32 up to
     * at most 2048 subintervals.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit
     * @param b
     *            upper limit
     * @param epsTol
     *            error tolerance for {@code |I_2N - I_N|}
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate1D(DFunction f, double a, double b, double epsTol) {
        return integrate1D(f, a, b, epsTol, MIN_N_1D, MAX_N_1D);
    }

    /**
     * Integrates {@code f} over {@code [a, b]} with an explicit grid range, for
     * an integrand known to be harder or easier than the default assumes. Both
     * bounds are rounded down to a power of two in {@code [8, 2048]}, and
     * {@code maxN} is raised to {@code 2 * minN} if it is smaller, since a
     * single grid yields no error estimate.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit
     * @param b
     *            upper limit
     * @param epsTol
     *            error tolerance for {@code |I_2N - I_N|}
     * @param minN
     *            subintervals of the first grid, which has {@code minN + 1}
     *            points
     * @param maxN
     *            subintervals of the last grid
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate1D(DFunction f, double a, double b, double epsTol, int minN, int maxN) {
        int lo = normalize(minN);
        int hi = normalizeMax(maxN, lo);

        double c = (b - a) / 2.0;
        double d = (b + a) / 2.0;

        int n = lo;
        double[] fv = evaluate1D(f, n, c, d);
        double previous = sum1D(fv, n) * c;
        double err = Double.POSITIVE_INFINITY;

        while (n < hi) {
            fv = refine1D(f, fv, n, c, d);
            n *= 2;
            double value = sum1D(fv, n) * c;
            err = Math.abs(value - previous);
            previous = value;
            if (err <= epsTol) {
                return new IntegralResult(value, err, n + 1, true);
            }
        }
        return new IntegralResult(previous, err, n + 1, false);
    }

    private static double[] evaluate1D(DFunction f, int n, double c, double d) {
        double[] xi = table(n).xi;
        double[] fv = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            fv[i] = f.apply(c * xi[i] + d);
        }
        return fv;
    }

    /**
     * Point {@code i} of the coarse grid is point {@code 2i} of the fine one, so
     * only the odd indices have to be evaluated.
     */
    private static double[] refine1D(DFunction f, double[] coarse, int n, double c, double d) {
        int fine = 2 * n;
        double[] xi = table(fine).xi;
        double[] fv = new double[fine + 1];
        for (int i = 0; i <= n; i++) {
            fv[2 * i] = coarse[i];
        }
        for (int i = 1; i < fine; i += 2) {
            fv[i] = f.apply(c * xi[i] + d);
        }
        return fv;
    }

    private static double sum1D(double[] fv, int n) {
        double[] w = table(n).w;
        double sum = 0.0;
        for (int i = 0; i <= n; i++) {
            sum += w[i] * fv[i];
        }
        return sum;
    }

    // =========================================================================
    // TWO DIMENSIONS
    // =========================================================================

    /**
     * Integrates {@code f} over the rectangle, doubling both axes together from
     * 16 up to at most 256 subintervals per axis.
     *
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
     *            error tolerance for {@code |I_2N - I_N|}
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate2D(DBiFunction f, double ax, double bx, double ay, double by,
            double epsTol) {
        return integrate2D(f, ax, bx, ay, by, epsTol, MIN_N_2D, MAX_N_2D);
    }

    /**
     * Integrates {@code f} over the rectangle with an explicit grid range. Both
     * bounds are rounded down to a power of two in {@code [8, 2048]}, and
     * {@code maxN} is raised to {@code 2 * minN} if it is smaller, since a
     * single grid yields no error estimate. Note that one level costs four
     * times the previous one here.
     *
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
     *            error tolerance for {@code |I_2N - I_N|}
     * @param minN
     *            subintervals per axis of the first grid
     * @param maxN
     *            subintervals per axis of the last grid
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate2D(DBiFunction f, double ax, double bx, double ay, double by, double epsTol,
            int minN, int maxN) {
        int lo = normalize(minN);
        int hi = normalizeMax(maxN, lo);

        double cx = (bx - ax) / 2.0;
        double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0;
        double dy = (by + ay) / 2.0;
        double areaFactor = cx * cy;

        int n = lo;
        double[] fv = evaluate2D(f, n, cx, dx, cy, dy);
        double previous = sum2D(fv, n) * areaFactor;
        double err = Double.POSITIVE_INFINITY;

        while (n < hi) {
            fv = refine2D(f, fv, n, cx, dx, cy, dy);
            n *= 2;
            double value = sum2D(fv, n) * areaFactor;
            err = Math.abs(value - previous);
            previous = value;
            if (err <= epsTol) {
                return new IntegralResult(value, err, n + 1, true);
            }
        }
        return new IntegralResult(previous, err, n + 1, false);
    }

    private static double[] evaluate2D(DBiFunction f, int n, double cx, double dx, double cy, double dy) {
        double[] xi = table(n).xi;
        int m = n + 1;
        double[] fv = new double[m * m];
        IntStream.range(0, m).parallel().forEach(i -> {
            double x = cx * xi[i] + dx;
            int base = i * m;
            for (int j = 0; j < m; j++) {
                fv[base + j] = f.apply(x, cy * xi[j] + dy);
            }
        });
        return fv;
    }

    /**
     * An entry of the fine grid is inherited exactly when both of its indices
     * are even, so no marker array is needed to tell old from new.
     */
    private static double[] refine2D(DBiFunction f, double[] coarse, int n, double cx, double dx, double cy,
            double dy) {
        int fine = 2 * n;
        int mf = fine + 1;
        int mc = n + 1;
        double[] xi = table(fine).xi;
        double[] fv = new double[mf * mf];

        for (int i = 0; i < mc; i++) {
            int src = i * mc;
            int dst = 2 * i * mf;
            for (int j = 0; j < mc; j++) {
                fv[dst + 2 * j] = coarse[src + j];
            }
        }

        IntStream.range(0, mf).parallel().forEach(i -> {
            double x = cx * xi[i] + dx;
            int base = i * mf;
            boolean iEven = (i & 1) == 0;
            int start = iEven ? 1 : 0;
            int step = iEven ? 2 : 1;
            for (int j = start; j < mf; j += step) {
                fv[base + j] = f.apply(x, cy * xi[j] + dy);
            }
        });
        return fv;
    }

    private static double sum2D(double[] fv, int n) {
        double[] w = table(n).w;
        int m = n + 1;
        return IntStream.range(0, m).parallel().mapToDouble(i -> {
            double s = 0.0;
            int base = i * m;
            for (int j = 0; j < m; j++) {
                s += w[j] * fv[base + j];
            }
            return w[i] * s;
        }).sum();
    }

    // =========================================================================
    // THREE DIMENSIONS
    // =========================================================================

    /**
     * Integrates {@code f} over the box, doubling all three axes together from 8
     * up to at most 64 subintervals per axis.
     *
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
     *            error tolerance for {@code |I_2N - I_N|}
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate3D(DTriFunction f, double ax, double bx, double ay, double by, double az,
            double bz, double epsTol) {
        return integrate3D(f, ax, bx, ay, by, az, bz, epsTol, MIN_N_3D, MAX_N_3D);
    }

    /**
     * Integrates {@code f} over the box with an explicit grid range. Both bounds
     * are rounded down to a power of two in {@code [8, 2048]}, and {@code maxN}
     * is raised to {@code 2 * minN} if it is smaller, since a single grid yields
     * no error estimate. One level costs eight times the previous one here, and
     * the cached values need {@code 8 * (maxN + 1)^3} bytes, so the default
     * ceiling of 64 is not raised lightly.
     *
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
     *            error tolerance for {@code |I_2N - I_N|}
     * @param minN
     *            subintervals per axis of the first grid
     * @param maxN
     *            subintervals per axis of the last grid
     * @return the approximated integral, its error estimate, the grid used and
     *         whether the tolerance was met
     */
    public static IntegralResult integrate3D(DTriFunction f, double ax, double bx, double ay, double by, double az,
            double bz, double epsTol, int minN, int maxN) {
        int lo = normalize(minN);
        int hi = normalizeMax(maxN, lo);

        double cx = (bx - ax) / 2.0;
        double dx = (bx + ax) / 2.0;
        double cy = (by - ay) / 2.0;
        double dy = (by + ay) / 2.0;
        double cz = (bz - az) / 2.0;
        double dz = (bz + az) / 2.0;
        double volumeFactor = cx * cy * cz;

        int n = lo;
        double[] fv = evaluate3D(f, n, cx, dx, cy, dy, cz, dz);
        double previous = sum3D(fv, n) * volumeFactor;
        double err = Double.POSITIVE_INFINITY;

        while (n < hi) {
            fv = refine3D(f, fv, n, cx, dx, cy, dy, cz, dz);
            n *= 2;
            double value = sum3D(fv, n) * volumeFactor;
            err = Math.abs(value - previous);
            previous = value;
            if (err <= epsTol) {
                return new IntegralResult(value, err, n + 1, true);
            }
        }
        return new IntegralResult(previous, err, n + 1, false);
    }

    private static double[] evaluate3D(DTriFunction f, int n, double cx, double dx, double cy, double dy, double cz,
            double dz) {
        double[] xi = table(n).xi;
        int m = n + 1;
        double[] fv = new double[m * m * m];
        IntStream.range(0, m).parallel().forEach(i -> {
            double x = cx * xi[i] + dx;
            for (int j = 0; j < m; j++) {
                double y = cy * xi[j] + dy;
                int base = (i * m + j) * m;
                for (int k = 0; k < m; k++) {
                    fv[base + k] = f.apply(x, y, cz * xi[k] + dz);
                }
            }
        });
        return fv;
    }

    /** Inherited exactly when all three indices are even; see {@link #refine2D}. */
    private static double[] refine3D(DTriFunction f, double[] coarse, int n, double cx, double dx, double cy,
            double dy, double cz, double dz) {
        int fine = 2 * n;
        int mf = fine + 1;
        int mc = n + 1;
        double[] xi = table(fine).xi;
        double[] fv = new double[mf * mf * mf];

        for (int i = 0; i < mc; i++) {
            for (int j = 0; j < mc; j++) {
                int src = (i * mc + j) * mc;
                int dst = (2 * i * mf + 2 * j) * mf;
                for (int k = 0; k < mc; k++) {
                    fv[dst + 2 * k] = coarse[src + k];
                }
            }
        }

        IntStream.range(0, mf).parallel().forEach(i -> {
            double x = cx * xi[i] + dx;
            boolean iEven = (i & 1) == 0;
            for (int j = 0; j < mf; j++) {
                double y = cy * xi[j] + dy;
                int base = (i * mf + j) * mf;
                boolean bothEven = iEven && ((j & 1) == 0);
                int start = bothEven ? 1 : 0;
                int step = bothEven ? 2 : 1;
                for (int k = start; k < mf; k += step) {
                    fv[base + k] = f.apply(x, y, cz * xi[k] + dz);
                }
            }
        });
        return fv;
    }

    private static double sum3D(double[] fv, int n) {
        double[] w = table(n).w;
        int m = n + 1;
        return IntStream.range(0, m).parallel().mapToDouble(i -> {
            double si = 0.0;
            for (int j = 0; j < m; j++) {
                double sj = 0.0;
                int base = (i * m + j) * m;
                for (int k = 0; k < m; k++) {
                    sj += w[k] * fv[base + k];
                }
                si += w[j] * sj;
            }
            return w[i] * si;
        }).sum();
    }
}
