package math.linalg;

import math.MathConsts;

/**
 * Shared core of the least squares estimators in this package. Everything here
 * works on a thin singular value decomposition {@code X = U D V'} and covers
 * ordinary and ridge regression at once, the former being the latter at
 * {@code lambda == 0}.
 * <p>
 * Going through the singular values rather than through the normal equations
 * costs more arithmetic and saves the squaring of the condition number:
 * {@code cond(X'X) == cond(X)^2}, so a design that is merely awkward becomes
 * unsolvable as soon as it is squared.
 * <p>
 * The decomposition is the caller's: hand in a
 * {@link FlatParallelJacobiSVD.Result} of the design matrix and the methods
 * here read it. That is what makes them reusable -- one decomposition serves a
 * fit, its variances, its effective degrees of freedom and its rank check
 * without being formed four times.
 *
 * @since 1.5.2
 */
public final class SvdLeastSquares {

    /**
     * {@code beta = V diag(d_i / (d_i^2 + lambda)) U' y}. At {@code lambda == 0}
     * this is {@code V D^-1 U' y}, the ordinary least squares solution.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param y
     *            the regressand, length {@code svd.m}
     * @param lambda
     *            ridge penalty, {@code 0} for ordinary least squares
     * @return the coefficients, length {@code svd.n}
     * @throws IllegalArgumentException
     *             if {@code svd} or {@code y} is null, if {@code y} is not of
     *             length {@code svd.m}, or if {@code lambda} is negative or not
     *             finite
     */
    public static double[] solve(FlatParallelJacobiSVD.Result svd, double[] y, double lambda) {
        checkRegressand(svd, y);
        checkPenalty(lambda, "lambda");
        int m = svd.m;
        int n = svd.n;
        // z_i = (d_i / (d_i^2 + lambda)) * (U' y)_i
        double[] z = new double[n];
        for (int i = 0; i < n; i++) {
            double d = svd.sigma[i];
            double uty = 0.0;
            int col = i * m;
            for (int k = 0; k < m; k++) {
                uty += svd.U[col + k] * y[k];
            }
            z[i] = filter(d, lambda) * uty;
        }
        double[] beta = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += svd.V[j * n + i] * z[j];
            }
            beta[i] = sum;
        }
        return beta;
    }

    /**
     * {@code beta = V diag(1/d_i) U' y} over the singular values above
     * {@code tol}, the rest discarded. Where {@link #solve} shrinks a small
     * singular value by a penalty, this one drops it: at {@code lambda == 0}
     * the filter is {@code 1/d}, which is unbounded, so a caller that cannot
     * refuse a rank deficient matrix needs the truncated pseudo-inverse
     * instead. The result is the minimum norm solution of the reduced problem.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param y
     *            the regressand, length {@code svd.m}
     * @param tol
     *            singular values at or below this are discarded, customarily
     *            {@code sigma[0] * max(m, n) * eps}
     * @return the coefficients, length {@code svd.n}
     * @throws IllegalArgumentException
     *             if {@code svd} or {@code y} is null, if {@code y} is not of
     *             length {@code svd.m}, or if {@code tol} is negative or not
     *             finite
     */
    public static double[] solveTruncated(FlatParallelJacobiSVD.Result svd, double[] y, double tol) {
        checkRegressand(svd, y);
        checkPenalty(tol, "tol");
        int m = svd.m;
        int n = svd.n;
        double[] z = new double[n];
        for (int i = 0; i < n; i++) {
            double d = svd.sigma[i];
            if (d <= tol) {
                continue;
            }
            double uty = 0.0;
            int col = i * m;
            for (int k = 0; k < m; k++) {
                uty += svd.U[col + k] * y[k];
            }
            z[i] = uty / d;
        }
        double[] beta = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += svd.V[j * n + i] * z[j];
            }
            beta[i] = sum;
        }
        return beta;
    }

    /**
     * Diagonal of {@code V diag(d_i^2 / (d_i^2 + lambda)^2) V'}, which becomes
     * the variances of the coefficients once multiplied by the residual
     * variance. Every entry is a sum of squares over positive denominators and
     * therefore cannot come out negative, which is what the normal equations
     * could not guarantee.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param lambda
     *            ridge penalty, {@code 0} for ordinary least squares
     * @return the diagonal, length {@code svd.n}
     * @throws IllegalArgumentException
     *             if {@code svd} is null or {@code lambda} is negative or not
     *             finite
     */
    public static double[] varianceDiagonal(FlatParallelJacobiSVD.Result svd, double lambda) {
        checkDecomposition(svd);
        checkPenalty(lambda, "lambda");
        int n = svd.n;
        double[] f = new double[n];
        for (int i = 0; i < n; i++) {
            double fi = filter(svd.sigma[i], lambda);
            f[i] = fi * fi;
        }
        double[] diag = new double[n];
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                double v = svd.V[i * n + j];
                sum += v * v * f[i];
            }
            diag[j] = sum;
        }
        return diag;
    }

    /**
     * Full variance-covariance matrix {@code V diag(d_i^2 / (d_i^2 + lambda)^2)
     * V'} in the flat column-major layout of {@link DMatrix}, again still to be
     * scaled by the residual variance.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param lambda
     *            ridge penalty, {@code 0} for ordinary least squares
     * @return an {@code n x n} matrix, column-major
     * @throws IllegalArgumentException
     *             if {@code svd} is null or {@code lambda} is negative or not
     *             finite
     */
    public static double[] varianceMatrix(FlatParallelJacobiSVD.Result svd, double lambda) {
        checkDecomposition(svd);
        checkPenalty(lambda, "lambda");
        int n = svd.n;
        double[] f = new double[n];
        for (int i = 0; i < n; i++) {
            double fi = filter(svd.sigma[i], lambda);
            f[i] = fi * fi;
        }
        double[] cov = new double[n * n];
        for (int c = 0; c < n; c++) {
            for (int r = 0; r <= c; r++) {
                double sum = 0.0;
                for (int i = 0; i < n; i++) {
                    sum += svd.V[i * n + r] * svd.V[i * n + c] * f[i];
                }
                cov[c * n + r] = sum;
                cov[r * n + c] = sum;
            }
        }
        return cov;
    }

    /**
     * Effective degrees of freedom {@code sum d_i^2 / (d_i^2 + lambda)}, equal
     * to the number of columns at {@code lambda == 0} and falling towards zero
     * as the penalty grows.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param lambda
     *            ridge penalty
     * @return the effective degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code svd} is null or {@code lambda} is negative or not
     *             finite
     */
    public static double effectiveDf(FlatParallelJacobiSVD.Result svd, double lambda) {
        checkDecomposition(svd);
        checkPenalty(lambda, "lambda");
        double df = 0.0;
        for (int i = 0; i < svd.n; i++) {
            double dd = svd.sigma[i] * svd.sigma[i];
            df += dd / (dd + lambda);
        }
        return df;
    }

    /**
     * Index of the first singular value that falls below
     * {@code sigma[0] * max(m, n) * eps}, the customary numerical rank
     * criterion. Equivalent to {@link #rankDeficientAt(FlatParallelJacobiSVD.Result, double)}
     * at {@link #defaultRankTolerance(FlatParallelJacobiSVD.Result)}.
     *
     * @param svd
     *            decomposition of the design matrix
     * @return the index of the first negligible singular value, or {@code -1}
     *         if the matrix has full numerical rank
     * @throws IllegalArgumentException
     *             if {@code svd} is null
     */
    public static int rankDeficientAt(FlatParallelJacobiSVD.Result svd) {
        checkDecomposition(svd);
        return rankDeficientAt(svd, defaultRankTolerance(svd));
    }

    /**
     * Index of the first singular value that falls at or below
     * {@code relativeTolerance * sigma[0]}.
     * <p>
     * The tolerance is relative to the largest singular value and therefore
     * dimensionless: scaling the whole design matrix by a constant leaves the
     * verdict where it was. {@code 0.0} accepts every design that is not
     * exactly singular, which is the loosest defensible setting -- a singular
     * value of exactly zero has to be refused here, because the ordinary least
     * squares filter of {@link #solve} would silently hand back a truncated
     * fit for it rather than fail.
     * <p>
     * The customary criterion, {@link #defaultRankTolerance}, calls a design
     * rank deficient at the level where rounding alone could have produced the
     * singular value. Between that and zero lie the designs that are ill
     * conditioned rather than rank deficient, whose answers the singular value
     * route still reaches.
     *
     * @param svd
     *            decomposition of the design matrix
     * @param relativeTolerance
     *            a singular value is negligible when it is at or below this
     *            multiple of {@code sigma[0]}; finite and in {@code [0, 1)}
     * @return the index of the first negligible singular value, or {@code -1}
     *         if the matrix has full numerical rank at this tolerance
     * @throws IllegalArgumentException
     *             if {@code svd} is null or {@code relativeTolerance} is not
     *             finite or outside {@code [0, 1)}
     * @since 1.5.2
     */
    public static int rankDeficientAt(FlatParallelJacobiSVD.Result svd, double relativeTolerance) {
        checkDecomposition(svd);
        checkRankTolerance(relativeTolerance);
        if (svd.n == 0) {
            return -1;
        }
        double tol = svd.sigma[0] * relativeTolerance;
        for (int i = 0; i < svd.n; i++) {
            if (svd.sigma[i] <= tol) {
                return i;
            }
        }
        return -1;
    }

    /**
     * The customary numerical rank tolerance {@code max(m, n) * eps}, relative
     * to the largest singular value.
     *
     * @param svd
     *            decomposition of the design matrix
     * @return the tolerance {@link #rankDeficientAt(FlatParallelJacobiSVD.Result)} uses
     * @throws IllegalArgumentException
     *             if {@code svd} is null
     * @since 1.5.2
     */
    public static double defaultRankTolerance(FlatParallelJacobiSVD.Result svd) {
        checkDecomposition(svd);
        return Math.max(svd.m, svd.n) * MathConsts.MACH_EPS_DBL;
    }

    /**
     * Condition number {@code sigma[0] / sigma[n-1]} of the decomposed matrix,
     * {@code Infinity} if the smallest singular value is zero.
     *
     * @param svd
     *            decomposition of the design matrix
     * @return the condition number, or {@code NaN} if there are no columns
     * @throws IllegalArgumentException
     *             if {@code svd} is null
     * @since 1.5.2
     */
    public static double conditionNumber(FlatParallelJacobiSVD.Result svd) {
        checkDecomposition(svd);
        if (svd.n == 0) {
            return Double.NaN;
        }
        return svd.sigma[0] / svd.sigma[svd.n - 1];
    }

    private static void checkRankTolerance(double relativeTolerance) {
        if (!(relativeTolerance >= 0.0) || relativeTolerance >= 1.0) {
            throw new IllegalArgumentException(
                    "relativeTolerance must be in [0, 1) : " + relativeTolerance);
        }
    }

    private static void checkDecomposition(FlatParallelJacobiSVD.Result svd) {
        if (svd == null) {
            throw new IllegalArgumentException("svd is null");
        }
    }

    private static void checkRegressand(FlatParallelJacobiSVD.Result svd, double[] y) {
        checkDecomposition(svd);
        if (y == null) {
            throw new IllegalArgumentException("y is null");
        }
        if (y.length != svd.m) {
            throw new IllegalArgumentException("y.length != svd.m : " + y.length + " != " + svd.m);
        }
    }

    private static void checkPenalty(double value, String name) {
        if (value < 0.0 || Double.isNaN(value) || Double.isInfinite(value)) {
            throw new IllegalArgumentException(name + " must be finite and non-negative : " + value);
        }
    }

    /** {@code d / (d^2 + lambda)}, the shrinkage applied to one singular value. */
    private static double filter(double d, double lambda) {
        double denom = d * d + lambda;
        if (denom == 0.0) {
            return 0.0;
        }
        return d / denom;
    }

    private SvdLeastSquares() {
        throw new AssertionError();
    }
}
