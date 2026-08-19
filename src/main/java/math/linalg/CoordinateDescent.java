package math.linalg;

/**
 * Cyclic coordinate descent for the elastic net, the engine behind
 * {@link Lasso}. Minimizes
 * {@code (1/(2n)) ||y - X b||^2 + lambda * (alpha ||b||_1 + ((1 - alpha)/2) ||b||^2)}
 * over a {@link Standardization}, where every column has
 * {@code sum_i x_ij^2 == n}. See
 * <a href="https://en.wikipedia.org/wiki/Coordinate_descent">coordinate
 * descent</a>.
 * <p>
 * With the other coordinates held fixed the minimum over {@code b_j} is
 * available in closed form, which is the whole reason the method is used here:
 *
 * <pre>
 * z     = b_j + (x_j . r) / n            r = y - X b, kept up to date
 * b_j'  = softThreshold(z, lambda*alpha) / (1 + lambda*(1 - alpha))
 * </pre>
 *
 * The soft threshold returns an exact {@code 0.0} rather than something near
 * it, so the sparsity the lasso is used for is a property of the arithmetic
 * and not of a cutoff applied afterwards.
 */
final class CoordinateDescent {

    /** Default convergence tolerance, relative to the size of the coefficients. */
    static final double DEFAULT_TOL = 1.0e-9;
    /** Default sweep budget. Well posed problems need a few dozen. */
    static final int DEFAULT_MAX_SWEEPS = 100_000;

    /** Full residual recomputation interval, guarding against drift. */
    private static final int REFRESH_EVERY = 100;

    /** The outcome of one fit on standardized data. */
    static final class Fit {

        /** Coefficients in the standardized scale. */
        final double[] beta;
        /** Number of sweeps performed. */
        final int sweeps;
        /** Whether a full sweep stopped moving anything before the budget ran out. */
        final boolean converged;
        /** Largest violation of the optimality conditions at the returned point. */
        final double maxKktViolation;

        Fit(double[] beta, int sweeps, boolean converged, double maxKktViolation) {
            this.beta = beta;
            this.sweeps = sweeps;
            this.converged = converged;
            this.maxKktViolation = maxKktViolation;
        }
    }

    /**
     * Fits the standardized problem.
     *
     * @param std
     *            the centered and scaled data
     * @param lambda
     *            the overall penalty, {@code 0} or greater
     * @param alpha
     *            the mixing parameter in {@code [0, 1]}, {@code 1} being the
     *            pure L1 penalty
     * @param warmStart
     *            coefficients to start from, or {@code null} to start at zero;
     *            copied, never modified
     * @param tol
     *            convergence tolerance; a sweep counts as final when it moves
     *            no coefficient by more than {@code tol * (1 + max |b_j|)}
     * @param maxSweeps
     *            the sweep budget
     * @return the fit, see {@link Fit}
     */
    static Fit solve(Standardization std, double lambda, double alpha, double[] warmStart, double tol, int maxSweeps) {
        int n = std.n;
        int p = std.p;
        double[] x = std.x;
        double[] y = std.y;

        double l1 = lambda * alpha;
        double denom = 1.0 + lambda * (1.0 - alpha);

        double[] b = new double[p];
        if (warmStart != null) {
            System.arraycopy(warmStart, 0, b, 0, p);
        }
        double[] r = new double[n];
        residuals(x, y, b, r, n, p);

        int[] active = new int[p];
        int sweeps = 0;
        int lastRefresh = 0;
        boolean converged = false;

        while (sweeps < maxSweeps) {
            if (sweeps - lastRefresh >= REFRESH_EVERY) {
                residuals(x, y, b, r, n, p);
                lastRefresh = sweeps;
            }
            double delta = sweep(x, r, b, n, null, p, l1, denom);
            ++sweeps;
            if (delta <= threshold(b, tol)) {
                converged = true;
                break;
            }

            // Cycle over the coefficients that are currently non-zero until
            // they settle. The full sweep above is what lets a new one in, and
            // the full sweep of the next round is what proves none is missing.
            int count = 0;
            for (int j = 0; j < p; j++) {
                if (b[j] != 0.0) {
                    active[count++] = j;
                }
            }
            while (sweeps < maxSweeps) {
                if (sweeps - lastRefresh >= REFRESH_EVERY) {
                    residuals(x, y, b, r, n, p);
                    lastRefresh = sweeps;
                }
                double d = sweep(x, r, b, n, active, count, l1, denom);
                ++sweeps;
                if (d <= threshold(b, tol)) {
                    break;
                }
            }
        }

        // the incremental updates have accumulated rounding by now, and the
        // optimality check below is only worth as much as the residual it uses
        residuals(x, y, b, r, n, p);
        return new Fit(b, sweeps, converged, kktViolation(x, r, b, n, p, l1, lambda * (1.0 - alpha)));
    }

    /** One pass over the given coordinates; returns the largest step taken. */
    private static double sweep(double[] x, double[] r, double[] b, int n, int[] idx, int count, double l1,
            double denom) {
        double maxDelta = 0.0;
        for (int k = 0; k < count; k++) {
            int j = (idx == null) ? k : idx[k];
            int col = j * n;
            double dot = 0.0;
            for (int i = 0; i < n; i++) {
                dot += x[col + i] * r[i];
            }
            double old = b[j];
            double updated = softThreshold(old + dot / n, l1) / denom;
            if (updated != old) {
                double delta = updated - old;
                b[j] = updated;
                for (int i = 0; i < n; i++) {
                    r[i] -= x[col + i] * delta;
                }
                double step = Math.abs(delta);
                if (step > maxDelta) {
                    maxDelta = step;
                }
            }
        }
        return maxDelta;
    }

    /**
     * {@code sign(z) * max(|z| - g, 0)}, written so that the zero it produces
     * is exact.
     */
    private static double softThreshold(double z, double g) {
        if (z > g) {
            return z - g;
        }
        if (z < -g) {
            return z + g;
        }
        return 0.0;
    }

    /** Convergence threshold, scaled by the size of the current solution. */
    private static double threshold(double[] b, double tol) {
        double maxAbs = 0.0;
        for (int j = 0; j < b.length; j++) {
            double a = Math.abs(b[j]);
            if (a > maxAbs) {
                maxAbs = a;
            }
        }
        return tol * (1.0 + maxAbs);
    }

    /** {@code r = y - X b}, from scratch. */
    private static void residuals(double[] x, double[] y, double[] b, double[] r, int n, int p) {
        System.arraycopy(y, 0, r, 0, n);
        for (int j = 0; j < p; j++) {
            double bj = b[j];
            if (bj != 0.0) {
                int col = j * n;
                for (int i = 0; i < n; i++) {
                    r[i] -= x[col + i] * bj;
                }
            }
        }
    }

    /**
     * Largest violation of the stationarity conditions:
     * {@code (x_j . r)/n == l1*sign(b_j) + l2*b_j} where {@code b_j != 0}, and
     * {@code |(x_j . r)/n| <= l1} where it is zero.
     */
    private static double kktViolation(double[] x, double[] r, double[] b, int n, int p, double l1, double l2) {
        double worst = 0.0;
        for (int j = 0; j < p; j++) {
            int col = j * n;
            double dot = 0.0;
            for (int i = 0; i < n; i++) {
                dot += x[col + i] * r[i];
            }
            double g = dot / n;
            double violation;
            if (b[j] != 0.0) {
                violation = Math.abs(g - (((b[j] > 0.0) ? l1 : -l1) + l2 * b[j]));
            } else {
                violation = Math.max(0.0, Math.abs(g) - l1);
            }
            if (violation > worst) {
                worst = violation;
            }
        }
        return worst;
    }

    private CoordinateDescent() {
        throw new AssertionError();
    }
}
