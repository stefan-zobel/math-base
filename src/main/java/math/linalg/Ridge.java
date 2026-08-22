package math.linalg;

import java.util.Locale;

/**
 * Ridge regression, the least squares fit under an L2 penalty on the
 * coefficients:
 * {@code min ||y - X b||^2 + lambda * ||b||^2}. See
 * <a href="https://en.wikipedia.org/wiki/Ridge_regression">ridge
 * regression</a>.
 * <p>
 * Solved through the singular values of the design matrix, which shrinks the
 * contribution of the direction belonging to {@code d_i} by
 * {@code d_i^2 / (d_i^2 + lambda)}. The small singular values that make an
 * ill-conditioned least squares problem explode are damped the hardest, which
 * is the whole point of the method.
 * <p>
 * The penalty is not scale invariant, so the columns are centred and scaled to
 * unit root mean square before the fit and the coefficients are returned in the
 * original scale. Consequently the design matrix passed in must carry
 * <em>no</em> intercept column: the intercept is estimated outside the penalty
 * as {@code ybar - sum(beta_j * xbar_j)}.
 * <p>
 * No t values, p values or confidence intervals are reported. The ridge
 * estimator is biased, so an interval built from these standard errors would
 * cover {@code E[beta(lambda)]} rather than the true coefficient, and reporting
 * it in the shape of {@link LSSummary} would invite exactly that misreading.
 * The standard errors and the effective degrees of freedom are given because
 * both are well defined and both are needed to choose {@code lambda}.
 *
 * @since 1.5.2
 */
public final class Ridge {

    /** The result of one ridge fit. */
    public static final class Result {

        /** Coefficients in the original scale of the columns, without the intercept. */
        public final double[] beta;
        /** Intercept, estimated outside the penalty. */
        public final double intercept;
        /** The penalty this fit used. */
        public final double lambda;
        /** Fitted values, {@code intercept + X beta}. */
        public final double[] fitted;
        /** Residuals, {@code y - fitted}. */
        public final double[] residuals;
        /** Coefficient of determination, {@code 1 - RSS/SST}. */
        public final double rSquared;
        /** Residual variance, {@code RSS / (n - effectiveDf)}. */
        public final double sigmaHatSquared;
        /** Standard errors of {@link #beta}, in the same scale. */
        public final double[] standardErrors;
        /** {@code sum d_i^2 / (d_i^2 + lambda)}, equal to the column count at {@code lambda == 0}. */
        public final double effectiveDf;
        /** Whether the underlying decomposition converged. */
        public final boolean converged;

        Result(double[] beta, double intercept, double lambda, double[] fitted, double[] residuals, double rSquared,
                double sigmaHatSquared, double[] standardErrors, double effectiveDf, boolean converged) {
            this.beta = beta;
            this.intercept = intercept;
            this.lambda = lambda;
            this.fitted = fitted;
            this.residuals = residuals;
            this.rSquared = rSquared;
            this.sigmaHatSquared = sigmaHatSquared;
            this.standardErrors = standardErrors;
            this.effectiveDf = effectiveDf;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT, "Ridge(lambda=%.6g): R^2 = %.6f, effective df = %.3f, %d coefficients",
                    lambda, rSquared, effectiveDf, beta.length);
        }
    }

    /**
     * Fits {@code y} on {@code X} under the penalty {@code lambda}.
     *
     * @param X
     *            the design matrix, {@code n x p} with {@code n >= p} and
     *            <em>without</em> an intercept column
     * @param y
     *            the regressand, {@code n x 1}
     * @param lambda
     *            the penalty, {@code 0} or greater; {@code 0} gives the
     *            ordinary least squares fit on centred data
     * @return the fit, see {@link Result}
     * @throws IllegalArgumentException
     *             if the dimensions do not match, if {@code X} has more columns
     *             than rows, if {@code lambda} is negative or not finite, or if
     *             a column of {@code X} is constant
     */
    public static Result estimate(DMatrix X, DMatrix y, double lambda) {
        if (X.numRows() != y.numRows()) {
            throw new IllegalArgumentException("X.numRows != y.numRows : " + X.numRows() + " != " + y.numRows());
        }
        if (y.numColumns() != 1) {
            throw new IllegalArgumentException("y must have exactly one column, has " + y.numColumns());
        }
        if (X.numRows() < X.numColumns()) {
            throw new IllegalArgumentException(
                    "X must not have more columns than rows : " + X.numColumns() + " > " + X.numRows());
        }
        if (lambda < 0.0 || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalArgumentException("lambda must be finite and non-negative : " + lambda);
        }

        int n = X.numRows();
        int p = X.numColumns();
        double[] xs = X.getArrayUnsafe();
        double[] ys = y.getArrayUnsafe();

        ScaledDesign std = ScaledDesign.of(xs, ys, n, p, null);
        double yBar = std.yBar;

        // decomposeInPlace consumes std.x, which is ours and not needed afterwards
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decomposeInPlace(std.x, n, p);
        if (lambda == 0.0 && SvdLeastSquares.rankDeficientAt(svd) >= 0) {
            throw new IllegalArgumentException(
                    "X is rank deficient and lambda is 0, which is the case ridge regression exists to avoid; "
                            + "use a positive lambda");
        }

        double[] beta = std.unscale(SvdLeastSquares.solve(svd, std.y, lambda));
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
            double c = ys[i] - yBar;
            sst += c * c;
        }

        double effectiveDf = SvdLeastSquares.effectiveDf(svd, lambda);
        double sigmaHatSquared = (n > effectiveDf) ? (rss / (n - effectiveDf)) : Double.NaN;
        double[] varDiag = SvdLeastSquares.varianceDiagonal(svd, lambda);
        double[] standardErrors = new double[p];
        for (int j = 0; j < p; j++) {
            // beta_j = betaScaled_j / scale_j, so the variance divides by scale_j^2
            standardErrors[j] = Math.sqrt(sigmaHatSquared * varDiag[j]) / std.scale[j];
        }
        double rSquared = (sst == 0.0) ? 1.0 : (1.0 - rss / sst);

        return new Result(beta, intercept, lambda, fitted, residuals, rSquared, sigmaHatSquared, standardErrors,
                effectiveDf, svd.converged);
    }

    private Ridge() {
        throw new AssertionError();
    }
}
