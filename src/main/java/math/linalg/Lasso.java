package math.linalg;

import java.util.Locale;

import math.rng.SplitMix64;

/**
 * Lasso and elastic net regression, the least squares fit under an L1 penalty
 * or under a mixture of an L1 and an L2 penalty:
 * {@code min (1/(2n)) ||y - X b||^2 + lambda * (alpha ||b||_1 + ((1 - alpha)/2) ||b||^2)}.
 * See <a href="https://en.wikipedia.org/wiki/Lasso_(statistics)">lasso</a> and
 * <a href="https://en.wikipedia.org/wiki/Elastic_net_regularization">elastic
 * net regularization</a>.
 * <p>
 * {@code alpha == 1} is the lasso, {@code alpha == 0} the ridge, anything
 * between the two the elastic net. Unlike {@link Ridge}, which shrinks every
 * coefficient towards zero without ever reaching it, the L1 penalty sets
 * coefficients <em>exactly</em> to zero and so selects variables. That is what
 * the method is used for, and it is a property of the soft threshold in
 * {@link CoordinateDescent} rather than of a cutoff applied afterwards. For a
 * group of strongly correlated predictors the pure lasso keeps one member of
 * the group and discards the rest more or less arbitrarily; a mixture around
 * {@code alpha == 0.5} keeps the group together, which is why the parameter
 * exists. For {@code alpha == 0} prefer {@link Ridge#estimate}, which solves
 * that case directly instead of iteratively.
 * <p>
 * The penalty is not scale invariant, so the columns are centered and scaled to
 * unit root mean square before the fit and the coefficients are returned in the
 * original scale. The design matrix passed in must therefore carry <em>no</em>
 * intercept column: the intercept is estimated outside the penalty as
 * {@code ybar - sum(beta_j * xbar_j)}. The factor {@code 1/(2n)} in the
 * objective is the usual convention and makes a given {@code lambda} mean the
 * same thing at different sample sizes.
 * <p>
 * No t values, p values or confidence intervals are reported. Beyond the bias
 * that already rules them out for {@link Ridge}, these coefficients come from a
 * model that was <em>selected</em> using the same data, so the classical
 * sampling distributions do not apply to them at all.
 * <p>
 * Choosing {@code lambda} is the real problem in practice, and
 * {@link #path(DMatrix, DMatrix, double)} and
 * {@link #cv(DMatrix, DMatrix, double, int)} exist for it: the former walks a
 * whole grid of penalties for a fraction of what the same number of independent
 * fits would cost, because each point starts from the previous solution, and the
 * latter picks a point on that grid by k-fold cross-validation.
 *
 * @since 1.5.2
 */
public final class Lasso {

    /** Default number of points on a regularization path. */
    private static final int DEFAULT_PATH_LENGTH = 100;
    /** Fixed seed of the fold assignment, so a run without a seed is reproducible. */
    private static final long DEFAULT_CV_SEED = 0x9E3779B97F4A7C15L;
    /** Below this the grid would run into a problem that is not identified. */
    private static final double MIN_ALPHA_FOR_PATH = 1.0e-6;

    /** The result of one elastic net fit at a single penalty. */
    public static final class Result {

        /** Coefficients in the original scale of the columns, without the intercept. */
        public final double[] beta;
        /** Intercept, estimated outside the penalty. */
        public final double intercept;
        /** The penalty this fit used. */
        public final double lambda;
        /** The mixing parameter this fit used. */
        public final double alpha;
        /** How many entries of {@link #beta} are not exactly zero. */
        public final int nonZeros;
        /** Fitted values, {@code intercept + X beta}. */
        public final double[] fitted;
        /** Residuals, {@code y - fitted}. */
        public final double[] residuals;
        /** Coefficient of determination, {@code 1 - RSS/SST}. */
        public final double rSquared;
        /** Number of coordinate descent sweeps performed. */
        public final int sweeps;
        /** Whether the sweeps settled before the budget ran out. */
        public final boolean converged;

        Result(double[] beta, double intercept, double lambda, double alpha, int nonZeros, double[] fitted,
                double[] residuals, double rSquared, int sweeps, boolean converged) {
            this.beta = beta;
            this.intercept = intercept;
            this.lambda = lambda;
            this.alpha = alpha;
            this.nonZeros = nonZeros;
            this.fitted = fitted;
            this.residuals = residuals;
            this.rSquared = rSquared;
            this.sweeps = sweeps;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT,
                    "Lasso(lambda=%.6g, alpha=%.3f): R^2 = %.6f, %d of %d coefficients non-zero", lambda, alpha,
                    rSquared, nonZeros, beta.length);
        }
    }

    /** A whole regularization path, one fit per penalty. */
    public static final class Path {

        /** The penalties, descending; {@code lambdas[0]} zeroes every coefficient. */
        public final double[] lambdas;
        /** Coefficients per penalty, {@code [lambdas.length][p]}, original scale. */
        public final double[][] beta;
        /** Intercepts per penalty. */
        public final double[] intercepts;
        /** Number of non-zero coefficients per penalty. */
        public final int[] nonZeros;
        /** The mixing parameter this path used. */
        public final double alpha;
        /** Whether every fit on the path converged. */
        public final boolean converged;

        Path(double[] lambdas, double[][] beta, double[] intercepts, int[] nonZeros, double alpha, boolean converged) {
            this.lambdas = lambdas;
            this.beta = beta;
            this.intercepts = intercepts;
            this.nonZeros = nonZeros;
            this.alpha = alpha;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT,
                    "Lasso.Path(alpha=%.3f): %d penalties from %.6g down to %.6g, up to %d non-zeros", alpha,
                    lambdas.length, lambdas[0], lambdas[lambdas.length - 1], nonZeros[nonZeros.length - 1]);
        }
    }

    /** The result of k-fold cross-validation over a regularization path. */
    public static final class CvResult {

        /** The path, fitted on all rows. */
        public final Path path;
        /** Mean squared prediction error per penalty, averaged over the folds. */
        public final double[] meanError;
        /** Standard error of {@link #meanError} across the folds. */
        public final double[] standardError;
        /** The penalty with the smallest {@link #meanError}. */
        public final double lambdaMin;
        /** Largest penalty whose error is within one standard error of the minimum. */
        public final double lambda1se;
        /** The fit on all rows at {@link #lambdaMin}. */
        public final Result best;
        /** How many folds were used. */
        public final int folds;

        CvResult(Path path, double[] meanError, double[] standardError, double lambdaMin, double lambda1se,
                Result best, int folds) {
            this.path = path;
            this.meanError = meanError;
            this.standardError = standardError;
            this.lambdaMin = lambdaMin;
            this.lambda1se = lambda1se;
            this.best = best;
            this.folds = folds;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT,
                    "Lasso.CvResult(%d folds): lambdaMin = %.6g (%d non-zeros), lambda1se = %.6g", folds, lambdaMin,
                    best.nonZeros, lambda1se);
        }
    }

    /**
     * Fits {@code y} on {@code X} under the pure lasso penalty, that is with
     * {@code alpha == 1}.
     *
     * @param X
     *            the design matrix, {@code n x p} and <em>without</em> an
     *            intercept column
     * @param y
     *            the regressand, {@code n x 1}
     * @param lambda
     *            the penalty, {@code 0} or greater
     * @return the fit, see {@link Result}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code lambda} is negative
     *             or not finite, or if a column of {@code X} is constant
     */
    public static Result estimate(DMatrix X, DMatrix y, double lambda) {
        return estimate(X, y, lambda, 1.0);
    }

    /**
     * Fits {@code y} on {@code X} under the elastic net penalty.
     *
     * @param X
     *            the design matrix, {@code n x p} and <em>without</em> an
     *            intercept column
     * @param y
     *            the regressand, {@code n x 1}
     * @param lambda
     *            the penalty, {@code 0} or greater
     * @param alpha
     *            the mixing parameter in {@code [0, 1]}; {@code 1} is the
     *            lasso, {@code 0} the ridge, for which {@link Ridge#estimate}
     *            is the more direct route
     * @return the fit, see {@link Result}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code lambda} is negative
     *             or not finite, if {@code alpha} lies outside {@code [0, 1]},
     *             or if a column of {@code X} is constant
     */
    public static Result estimate(DMatrix X, DMatrix y, double lambda, double alpha) {
        checkDimensions(X, y);
        checkLambda(lambda);
        checkAlpha(alpha, false);

        int n = X.numRows();
        int p = X.numColumns();
        Standardization std = Standardization.of(X.getArrayUnsafe(), y.getArrayUnsafe(), n, p, null);
        CoordinateDescent.Fit fit = CoordinateDescent.solve(std, lambda, alpha, null, CoordinateDescent.DEFAULT_TOL,
                CoordinateDescent.DEFAULT_MAX_SWEEPS);
        return result(X.getArrayUnsafe(), y.getArrayUnsafe(), n, p, std, fit.beta, lambda, alpha, fit.sweeps,
                fit.converged);
    }

    /**
     * The smallest penalty at which every coefficient of the elastic net fit is
     * zero, {@code max_j |x_j . y| / (n * alpha)} on the standardized data,
     * rounded up by one unit in the last place so that the guarantee holds in
     * floating point as well.
     *
     * @param X
     *            the design matrix, {@code n x p} and without an intercept
     *            column
     * @param y
     *            the regressand, {@code n x 1}
     * @param alpha
     *            the mixing parameter, greater than zero; the quantity is
     *            unbounded at {@code alpha == 0}
     * @return the penalty above which the fit is the empty model
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} lies outside
     *             {@code (0, 1]}, or if a column of {@code X} is constant
     */
    public static double lambdaMax(DMatrix X, DMatrix y, double alpha) {
        checkDimensions(X, y);
        checkAlpha(alpha, true);
        Standardization std = Standardization.of(X.getArrayUnsafe(), y.getArrayUnsafe(), X.numRows(), X.numColumns(),
                null);
        return lambdaMax(std, alpha);
    }

    /**
     * A regularization path of 100 penalties, log-spaced from
     * {@link #lambdaMax} down to a small fraction of it.
     *
     * @param X
     *            the design matrix, {@code n x p} and without an intercept
     *            column
     * @param y
     *            the regressand, {@code n x 1}
     * @param alpha
     *            the mixing parameter, greater than zero
     * @return the path, see {@link Path}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} lies outside
     *             {@code (0, 1]}, or if a column of {@code X} is constant
     */
    public static Path path(DMatrix X, DMatrix y, double alpha) {
        return path(X, y, alpha, DEFAULT_PATH_LENGTH, defaultMinRatio(X.numRows(), X.numColumns()));
    }

    /**
     * A regularization path over a grid of the given size. Each fit starts from
     * the solution of the previous, larger penalty, so a point on the path is
     * several times cheaper than the same fit started from zero would be.
     *
     * @param X
     *            the design matrix, {@code n x p} and without an intercept
     *            column
     * @param y
     *            the regressand, {@code n x 1}
     * @param alpha
     *            the mixing parameter, greater than zero
     * @param nLambda
     *            the number of penalties, at least {@code 1}
     * @param lambdaMinRatio
     *            the smallest penalty as a fraction of {@link #lambdaMax}, in
     *            {@code (0, 1)}
     * @return the path, see {@link Path}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} lies outside
     *             {@code (0, 1]}, if {@code nLambda} or {@code lambdaMinRatio}
     *             is out of range, or if a column of {@code X} is constant
     */
    public static Path path(DMatrix X, DMatrix y, double alpha, int nLambda, double lambdaMinRatio) {
        checkDimensions(X, y);
        checkAlpha(alpha, true);
        checkGrid(nLambda, lambdaMinRatio);

        int n = X.numRows();
        int p = X.numColumns();
        double[] xs = X.getArrayUnsafe();
        double[] ys = y.getArrayUnsafe();
        Standardization std = Standardization.of(xs, ys, n, p, null);
        double[] lambdas = grid(lambdaMax(std, alpha), nLambda, lambdaMinRatio);
        return path(std, lambdas, alpha);
    }

    /**
     * Picks a penalty by k-fold cross-validation, with a fixed default seed for
     * the fold assignment.
     *
     * @param X
     *            the design matrix, {@code n x p} and without an intercept
     *            column
     * @param y
     *            the regressand, {@code n x 1}
     * @param alpha
     *            the mixing parameter, greater than zero
     * @param folds
     *            the number of folds, at least {@code 2} and at most {@code n}
     * @return the cross-validation, see {@link CvResult}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} lies outside
     *             {@code (0, 1]}, if {@code folds} is out of range, or if a
     *             column of {@code X} is constant within a training fold
     */
    public static CvResult cv(DMatrix X, DMatrix y, double alpha, int folds) {
        return cv(X, y, alpha, folds, DEFAULT_CV_SEED);
    }

    /**
     * Picks a penalty by k-fold cross-validation.
     * <p>
     * The grid is computed once on all rows and reused in every fold, since
     * errors from different grids cannot be averaged. Each fold standardizes on
     * its <em>training</em> rows alone: standardizing on all rows first would
     * let the held out response into the fit and bend the error curve
     * downwards. Rows are assigned to folds by a shuffle driven by
     * {@code seed}, so the same seed reproduces the same partition and thus the
     * same answer.
     * <p>
     * Two penalties are reported. {@code lambdaMin} minimizes the mean error;
     * {@code lambda1se} is the largest penalty whose error is still within one
     * standard error of that minimum, the more conservative and usually sparser
     * choice.
     *
     * @param X
     *            the design matrix, {@code n x p} and without an intercept
     *            column
     * @param y
     *            the regressand, {@code n x 1}
     * @param alpha
     *            the mixing parameter, greater than zero
     * @param folds
     *            the number of folds, at least {@code 2} and at most {@code n}
     * @param seed
     *            seed of the fold assignment
     * @return the cross-validation, see {@link CvResult}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code alpha} lies outside
     *             {@code (0, 1]}, if {@code folds} is out of range, or if a
     *             column of {@code X} is constant within a training fold
     */
    public static CvResult cv(DMatrix X, DMatrix y, double alpha, int folds, long seed) {
        checkDimensions(X, y);
        checkAlpha(alpha, true);
        int n = X.numRows();
        int p = X.numColumns();
        if (folds < 2 || folds > n) {
            throw new IllegalArgumentException("folds must lie between 2 and the row count " + n + " : " + folds);
        }

        double[] xs = X.getArrayUnsafe();
        double[] ys = y.getArrayUnsafe();
        Standardization all = Standardization.of(xs, ys, n, p, null);
        double[] lambdas = grid(lambdaMax(all, alpha), DEFAULT_PATH_LENGTH, defaultMinRatio(n, p));
        int nLambda = lambdas.length;

        int[] foldOf = assignFolds(n, folds, seed);
        // squared prediction error, per fold and penalty, on the held out rows
        double[][] foldError = new double[folds][nLambda];
        for (int k = 0; k < folds; k++) {
            int testCount = 0;
            for (int i = 0; i < n; i++) {
                if (foldOf[i] == k) {
                    ++testCount;
                }
            }
            int[] train = new int[n - testCount];
            int[] test = new int[testCount];
            int ti = 0;
            int si = 0;
            for (int i = 0; i < n; i++) {
                if (foldOf[i] == k) {
                    test[si++] = i;
                } else {
                    train[ti++] = i;
                }
            }

            Standardization std;
            try {
                std = Standardization.of(xs, ys, n, p, train);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("fold " + k + " of " + folds + ": " + e.getMessage(), e);
            }
            Path foldPath = path(std, lambdas, alpha);
            for (int l = 0; l < nLambda; l++) {
                double[] beta = foldPath.beta[l];
                double intercept = foldPath.intercepts[l];
                double sum = 0.0;
                for (int t = 0; t < test.length; t++) {
                    int i = test[t];
                    double v = intercept;
                    for (int j = 0; j < p; j++) {
                        v += beta[j] * xs[j * n + i];
                    }
                    double e = ys[i] - v;
                    sum += e * e;
                }
                foldError[k][l] = sum / test.length;
            }
        }

        double[] meanError = new double[nLambda];
        double[] standardError = new double[nLambda];
        for (int l = 0; l < nLambda; l++) {
            double mean = 0.0;
            for (int k = 0; k < folds; k++) {
                mean += foldError[k][l];
            }
            mean /= folds;
            double ss = 0.0;
            for (int k = 0; k < folds; k++) {
                double d = foldError[k][l] - mean;
                ss += d * d;
            }
            meanError[l] = mean;
            standardError[l] = Math.sqrt(ss / (folds - 1.0)) / Math.sqrt(folds);
        }

        int best = 0;
        for (int l = 1; l < nLambda; l++) {
            if (meanError[l] < meanError[best]) {
                best = l;
            }
        }
        double lambdaMin = lambdas[best];
        // the grid descends, so the first point inside the band is the largest
        double bound = meanError[best] + standardError[best];
        double lambda1se = lambdaMin;
        for (int l = 0; l <= best; l++) {
            if (meanError[l] <= bound) {
                lambda1se = lambdas[l];
                break;
            }
        }

        Path fullPath = path(all, lambdas, alpha);
        Result bestFit = estimate(X, y, lambdaMin, alpha);
        return new CvResult(fullPath, meanError, standardError, lambdaMin, lambda1se, bestFit, folds);
    }

    /** Fits the given grid on already standardized data, warm starting down the path. */
    private static Path path(Standardization std, double[] lambdas, double alpha) {
        int nLambda = lambdas.length;
        double[][] beta = new double[nLambda][];
        double[] intercepts = new double[nLambda];
        int[] nonZeros = new int[nLambda];
        boolean converged = true;

        double[] warm = null;
        for (int l = 0; l < nLambda; l++) {
            CoordinateDescent.Fit fit = CoordinateDescent.solve(std, lambdas[l], alpha, warm,
                    CoordinateDescent.DEFAULT_TOL, CoordinateDescent.DEFAULT_MAX_SWEEPS);
            warm = fit.beta;
            converged &= fit.converged;
            double[] b = std.unscale(fit.beta);
            beta[l] = b;
            intercepts[l] = std.intercept(b);
            nonZeros[l] = countNonZeros(fit.beta);
        }
        return new Path(lambdas, beta, intercepts, nonZeros, alpha, converged);
    }

    /** Builds the reported result from a fit on standardized data. */
    private static Result result(double[] xs, double[] ys, int n, int p, Standardization std, double[] betaScaled,
            double lambda, double alpha, int sweeps, boolean converged) {
        double[] beta = std.unscale(betaScaled);
        double intercept = std.intercept(beta);

        double[] fitted = new double[n];
        double[] residuals = new double[n];
        double rss = 0.0;
        double sst = 0.0;
        for (int i = 0; i < n; i++) {
            double v = intercept;
            for (int j = 0; j < p; j++) {
                v += beta[j] * xs[j * n + i];
            }
            fitted[i] = v;
            double r = ys[i] - v;
            residuals[i] = r;
            rss += r * r;
            double c = ys[i] - std.yBar;
            sst += c * c;
        }
        double rSquared = (sst == 0.0) ? 1.0 : (1.0 - rss / sst);
        return new Result(beta, intercept, lambda, alpha, countNonZeros(betaScaled), fitted, residuals, rSquared,
                sweeps, converged);
    }

    /** Counts on the standardized coefficients, which share their zeros with the unscaled ones. */
    private static int countNonZeros(double[] beta) {
        int count = 0;
        for (int j = 0; j < beta.length; j++) {
            if (beta[j] != 0.0) {
                ++count;
            }
        }
        return count;
    }

    private static double lambdaMax(Standardization std, double alpha) {
        int n = std.n;
        double maxDot = 0.0;
        for (int j = 0; j < std.p; j++) {
            int col = j * n;
            double dot = 0.0;
            for (int i = 0; i < n; i++) {
                dot += std.x[col + i] * std.y[i];
            }
            double a = Math.abs(dot);
            if (a > maxDot) {
                maxDot = a;
            }
        }
        // one ulp up, so that the soft threshold really does return zero for
        // the column that attains the maximum rather than a rounding remnant
        return Math.nextUp(maxDot / (n * alpha));
    }

    private static double[] grid(double lambdaMax, int nLambda, double lambdaMinRatio) {
        double[] lambdas = new double[nLambda];
        lambdas[0] = lambdaMax;
        if (nLambda == 1) {
            return lambdas;
        }
        double logMax = Math.log(lambdaMax);
        double logMin = Math.log(lambdaMinRatio * lambdaMax);
        for (int l = 1; l < nLambda; l++) {
            lambdas[l] = Math.exp(logMax + (logMin - logMax) * (l / (nLambda - 1.0)));
        }
        return lambdas;
    }

    /**
     * glmnet's rule: with at least as many columns as rows the unpenalized
     * problem is not identified, so the grid stops well short of it.
     */
    private static double defaultMinRatio(int n, int p) {
        return (n > p) ? 1.0e-4 : 1.0e-2;
    }

    /** Balanced fold labels from a Fisher-Yates shuffle of the row indices. */
    private static int[] assignFolds(int n, int folds, long seed) {
        int[] order = new int[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        SplitMix64 rng = new SplitMix64(seed);
        for (int i = n - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            int tmp = order[i];
            order[i] = order[j];
            order[j] = tmp;
        }
        int[] foldOf = new int[n];
        for (int i = 0; i < n; i++) {
            foldOf[order[i]] = i % folds;
        }
        return foldOf;
    }

    private static void checkDimensions(DMatrix X, DMatrix y) {
        if (X.numRows() != y.numRows()) {
            throw new IllegalArgumentException("X.numRows != y.numRows : " + X.numRows() + " != " + y.numRows());
        }
        if (y.numColumns() != 1) {
            throw new IllegalArgumentException("y must have exactly one column, has " + y.numColumns());
        }
    }

    private static void checkLambda(double lambda) {
        if (lambda < 0.0 || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalArgumentException("lambda must be finite and non-negative : " + lambda);
        }
    }

    private static void checkAlpha(double alpha, boolean positive) {
        if (Double.isNaN(alpha) || alpha < 0.0 || alpha > 1.0) {
            throw new IllegalArgumentException("alpha must lie in [0, 1] : " + alpha);
        }
        if (positive && alpha < MIN_ALPHA_FOR_PATH) {
            throw new IllegalArgumentException("alpha must be greater than " + MIN_ALPHA_FOR_PATH
                    + " here, since the largest useful penalty grows without bound as alpha approaches 0"
                    + "; for a pure L2 penalty use Ridge.estimate : " + alpha);
        }
    }

    private static void checkGrid(int nLambda, double lambdaMinRatio) {
        if (nLambda < 1) {
            throw new IllegalArgumentException("nLambda must be at least 1 : " + nLambda);
        }
        if (!(lambdaMinRatio > 0.0) || lambdaMinRatio >= 1.0) {
            throw new IllegalArgumentException("lambdaMinRatio must lie in (0, 1) : " + lambdaMinRatio);
        }
    }

    private Lasso() {
        throw new AssertionError();
    }
}
