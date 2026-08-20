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
 */
final class SvdLeastSquares {

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
     */
    static double[] solve(FlatParallelJacobiSVD.Result svd, double[] y, double lambda) {
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
     */
    static double[] solveTruncated(FlatParallelJacobiSVD.Result svd, double[] y, double tol) {
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
     */
    static double[] varianceDiagonal(FlatParallelJacobiSVD.Result svd, double lambda) {
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
     */
    static double[] varianceMatrix(FlatParallelJacobiSVD.Result svd, double lambda) {
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
     */
    static double effectiveDf(FlatParallelJacobiSVD.Result svd, double lambda) {
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
     * criterion.
     *
     * @param svd
     *            decomposition of the design matrix
     * @return the index of the first negligible singular value, or {@code -1}
     *         if the matrix has full numerical rank
     */
    static int rankDeficientAt(FlatParallelJacobiSVD.Result svd) {
        if (svd.n == 0) {
            return -1;
        }
        double tol = svd.sigma[0] * Math.max(svd.m, svd.n) * MathConsts.MACH_EPS_DBL;
        for (int i = 0; i < svd.n; i++) {
            if (svd.sigma[i] <= tol) {
                return i;
            }
        }
        return -1;
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
