package math.optim;

import math.MathConsts;
import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;

/**
 * The Nelder-Mead downhill simplex method for the minimization of a
 * {@link DMultiFunction} without derivatives. A simplex of {@code n + 1}
 * vertices is moved through the search space by reflection, expansion,
 * contraction and shrinking until it collapses onto a minimum.
 * <p>
 * The method is a heuristic. For {@code n >= 2} there is no convergence proof,
 * and the simplex can degenerate and stall at a point that is not stationary.
 * Restarting with a fresh simplex around the best point found so far is the
 * usual remedy and is done once by default.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method">Wikipedia
 * Nelder-Mead method</a>.
 *
 * @since 1.5.1
 */
public final class NelderMead {

    /** Relative spread of the function values that stops the iteration. */
    private static final double DEFAULT_F_TOL = 1.0e-10;
    /** Relative diameter of the simplex that stops the iteration. */
    private static final double DEFAULT_X_TOL = 1.0e-8;
    /** Number of restarts if none is given. */
    private static final int DEFAULT_MAX_RESTARTS = 1;
    /** Factor of the default iteration budget, which grows as {@code n^2}. */
    private static final int ITERATION_FACTOR = 200;
    /** Initial step for a non-zero coordinate of the starting point. */
    private static final double RELATIVE_STEP = 0.05;
    /** Initial step for a coordinate of the starting point that is zero. */
    private static final double ABSOLUTE_STEP = 0.00025;

    private final double fTol;
    private final double xTol;
    private final int maxIterations;
    private final int maxRestarts;

    /**
     * Creates a Nelder-Mead minimizer with default settings.
     */
    public NelderMead() {
        this(DEFAULT_F_TOL, DEFAULT_X_TOL, 0, DEFAULT_MAX_RESTARTS);
    }

    /**
     * Creates a Nelder-Mead minimizer.
     *
     * @param fTol
     *            relative tolerance on the spread of the function values over
     *            the simplex
     * @param xTol
     *            relative tolerance on the diameter of the simplex
     * @param maxIterations
     *            maximal number of iterations, shared by all restarts;
     *            {@code 0} or less selects {@code 200 * n * (n + 1)}
     * @param maxRestarts
     *            how often the search may be restarted with a fresh simplex
     *            around the best point found so far
     */
    public NelderMead(double fTol, double xTol, int maxIterations, int maxRestarts) {
        if (!(fTol >= 0.0)) {
            throw new IllegalArgumentException("fTol must be >= 0.0 : " + fTol);
        }
        if (!(xTol >= 0.0)) {
            throw new IllegalArgumentException("xTol must be >= 0.0 : " + xTol);
        }
        if (maxRestarts < 0) {
            throw new IllegalArgumentException("maxRestarts must be >= 0 : " + maxRestarts);
        }
        this.fTol = fTol;
        this.xTol = xTol;
        this.maxIterations = maxIterations;
        this.maxRestarts = maxRestarts;
    }

    /**
     * Minimizes {@code f}, starting from {@code start}, with an initial simplex
     * derived from the starting point.
     *
     * @param f
     *            the function to minimize
     * @param start
     *            the starting point, not modified
     * @return the point the search has settled on, the function value there,
     *         the number of iterations and whether the tolerances were met
     */
    public DMultiFunctionEval minimize(DMultiFunction f, double[] start) {
        return minimize(f, start, null);
    }

    /**
     * Minimizes {@code f}, starting from {@code start}, with an initial simplex
     * spanned by the given per-coordinate steps. Supplying the steps is the way
     * to handle a badly scaled problem, where the default simplex is either far
     * too small or far too large in some of the coordinates.
     *
     * @param f
     *            the function to minimize
     * @param start
     *            the starting point, not modified
     * @param initialStep
     *            the edge of the initial simplex in each coordinate, all
     *            entries non-zero and finite, or {@code null} for the default
     * @return the point the search has settled on, the function value there,
     *         the number of iterations and whether the tolerances were met
     */
    public DMultiFunctionEval minimize(DMultiFunction f, double[] start, double[] initialStep) {
        if (f == null) {
            throw new IllegalArgumentException("f is null");
        }
        if (start == null) {
            throw new IllegalArgumentException("start is null");
        }
        int n = start.length;
        if (n < 1) {
            throw new IllegalArgumentException("start must have at least one component");
        }
        for (int i = 0; i < n; i++) {
            if (!isFinite(start[i])) {
                throw new IllegalArgumentException("start[" + i + "] is not finite : " + start[i]);
            }
        }
        double[] fixedStep = null;
        if (initialStep != null) {
            if (initialStep.length != n) {
                throw new IllegalArgumentException(
                        "initialStep has length " + initialStep.length + ", expected " + n);
            }
            for (int i = 0; i < n; i++) {
                if (initialStep[i] == 0.0 || !isFinite(initialStep[i])) {
                    throw new IllegalArgumentException(
                            "initialStep[" + i + "] must be non-zero and finite : " + initialStep[i]);
                }
            }
            fixedStep = initialStep.clone();
        }

        double[] best = start.clone();
        double bestValue = f.apply(best);
        if (!isFinite(bestValue)) {
            throw new IllegalArgumentException("f is not finite at the starting point : " + bestValue);
        }

        int budget = (maxIterations > 0) ? maxIterations : defaultBudget(n);
        int used = 0;
        boolean converged = false;
        for (int restart = 0; restart <= maxRestarts; restart++) {
            double[] step = (fixedStep != null) ? fixedStep : defaultStep(best);
            Simplex s = build(f, best, step);
            run(f, s, budget - used);
            used += s.iterations;
            converged = s.converged;
            if (s.fv[0] < bestValue) {
                bestValue = s.fv[0];
                System.arraycopy(s.x[0], 0, best, 0, n);
            }
            if (!converged || used >= budget) {
                break;
            }
        }
        return new DMultiFunctionEval(best, bestValue, used, converged);
    }

    /** The working simplex, kept sorted by ascending function value. */
    private static final class Simplex {
        final double[][] x;
        final double[] fv;
        int iterations;
        boolean converged;

        Simplex(int n) {
            x = new double[n + 1][n];
            fv = new double[n + 1];
        }
    }

    private void run(DMultiFunction f, Simplex s, int budget) {
        int n = s.fv.length - 1;
        // Nelder-Mead coefficients. Below n = 3 the adaptive choice of Gao and
        // Han either coincides with the classic one (n = 2) or degenerates
        // (sigma = 0 at n = 1, which would collapse the simplex to a point)
        double rho = 1.0;
        double chi;
        double gamma;
        double sigma;
        if (n <= 2) {
            chi = 2.0;
            gamma = 0.5;
            sigma = 0.5;
        } else {
            chi = 1.0 + 2.0 / n;
            gamma = 0.75 - 1.0 / (2.0 * n);
            sigma = 1.0 - 1.0 / n;
        }
        double[] centroid = new double[n];
        double[] trial = new double[n];
        double[] alt = new double[n];

        s.iterations = 0;
        s.converged = hasConverged(s);
        while (!s.converged && s.iterations < budget) {
            s.iterations++;
            // centroid of the best n vertices, i.e. of all but the worst
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += s.x[j][i];
                }
                centroid[i] = sum / n;
            }
            double[] worst = s.x[n];
            double fBest = s.fv[0];
            double fSecondWorst = s.fv[n - 1];
            double fWorst = s.fv[n];

            // A NaN never satisfies any of these comparisons and is therefore
            // never accepted; it also sorts as the worst vertex, so the search
            // walks away from a region where f is undefined
            double fr = evaluate(f, centroid, worst, rho, trial);
            if (fr < fBest) {
                double fe = evaluate(f, centroid, worst, rho * chi, alt);
                if (fe < fr) {
                    replaceWorst(s, alt, fe);
                } else {
                    replaceWorst(s, trial, fr);
                }
            } else if (fr < fSecondWorst) {
                replaceWorst(s, trial, fr);
            } else if (fr < fWorst) {
                double fc = evaluate(f, centroid, worst, rho * gamma, alt);
                if (fc <= fr) {
                    replaceWorst(s, alt, fc);
                } else {
                    shrink(f, s, sigma);
                }
            } else {
                double fc = evaluate(f, centroid, worst, -gamma, alt);
                if (fc < fWorst) {
                    replaceWorst(s, alt, fc);
                } else {
                    shrink(f, s, sigma);
                }
            }
            s.converged = hasConverged(s);
        }
    }

    private boolean hasConverged(Simplex s) {
        int n = s.fv.length - 1;
        double fLow = s.fv[0];
        double fHigh = s.fv[n];
        // negated so that a NaN spread counts as "not converged"
        if (!(fHigh - fLow <= fTol * 0.5 * (Math.abs(fLow) + Math.abs(fHigh)) + MathConsts.MIN_TOL)) {
            return false;
        }
        double[] x0 = s.x[0];
        double scale = 0.0;
        for (int i = 0; i < n; i++) {
            double a = Math.abs(x0[i]);
            if (a > scale) {
                scale = a;
            }
        }
        double diameter = 0.0;
        for (int j = 1; j <= n; j++) {
            double[] xj = s.x[j];
            for (int i = 0; i < n; i++) {
                double d = Math.abs(xj[i] - x0[i]);
                if (d > diameter) {
                    diameter = d;
                }
            }
        }
        return diameter <= xTol * (1.0 + scale);
    }

    /** Reflects {@code away} through {@code centroid} and evaluates {@code f}. */
    private static double evaluate(DMultiFunction f, double[] centroid, double[] away, double coefficient,
            double[] out) {
        for (int i = 0; i < out.length; i++) {
            out[i] = centroid[i] + coefficient * (centroid[i] - away[i]);
        }
        return f.apply(out);
    }

    private static Simplex build(DMultiFunction f, double[] base, double[] step) {
        int n = base.length;
        Simplex s = new Simplex(n);
        System.arraycopy(base, 0, s.x[0], 0, n);
        for (int i = 0; i < n; i++) {
            double[] v = s.x[i + 1];
            System.arraycopy(base, 0, v, 0, n);
            v[i] = base[i] + step[i];
        }
        for (int j = 0; j <= n; j++) {
            s.fv[j] = f.apply(s.x[j]);
        }
        sort(s);
        return s;
    }

    /**
     * The number of iterations to reach a given accuracy grows roughly like
     * {@code n^2}, so a budget linear in {@code n} would cut off well behaved
     * problems in higher dimensions. This is a cap, not a cost.
     */
    private static int defaultBudget(int n) {
        long budget = (long) ITERATION_FACTOR * n * (n + 1L);
        return (budget > Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) budget;
    }

    private static double[] defaultStep(double[] x) {
        double[] step = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            step[i] = (x[i] != 0.0) ? RELATIVE_STEP * Math.abs(x[i]) : ABSOLUTE_STEP;
        }
        return step;
    }

    private static void shrink(DMultiFunction f, Simplex s, double sigma) {
        int n = s.fv.length - 1;
        double[] x0 = s.x[0];
        for (int j = 1; j <= n; j++) {
            double[] xj = s.x[j];
            for (int i = 0; i < n; i++) {
                xj[i] = x0[i] + sigma * (xj[i] - x0[i]);
            }
            s.fv[j] = f.apply(xj);
        }
        sort(s);
    }

    /**
     * Overwrites the worst vertex and moves it to its place. Everything below
     * the last index is still sorted, so one insertion step is enough.
     */
    private static void replaceWorst(Simplex s, double[] point, double value) {
        int n = s.fv.length - 1;
        System.arraycopy(point, 0, s.x[n], 0, point.length);
        double[] xn = s.x[n];
        s.fv[n] = value;
        int i = n - 1;
        while (i >= 0 && Double.compare(s.fv[i], value) > 0) {
            s.fv[i + 1] = s.fv[i];
            s.x[i + 1] = s.x[i];
            i--;
        }
        s.fv[i + 1] = value;
        s.x[i + 1] = xn;
    }

    /** Insertion sort by ascending value; NaN sorts last, i.e. as the worst. */
    private static void sort(Simplex s) {
        int m = s.fv.length;
        for (int j = 1; j < m; j++) {
            double fj = s.fv[j];
            double[] xj = s.x[j];
            int i = j - 1;
            while (i >= 0 && Double.compare(s.fv[i], fj) > 0) {
                s.fv[i + 1] = s.fv[i];
                s.x[i + 1] = s.x[i];
                i--;
            }
            s.fv[i + 1] = fj;
            s.x[i + 1] = xj;
        }
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
