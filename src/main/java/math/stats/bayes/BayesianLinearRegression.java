package math.stats.bayes;

import math.cern.FastGamma;
import math.distribution.InverseGamma;
import math.distribution.MultivariateNormal;
import math.distribution.StudentT;
import math.fun.DFunction;
import math.linalg.DMatrix;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.ScaledDesign;
import math.linalg.SvdLeastSquares;
import math.optim.BrentMinimizer;

/**
 * Linear regression with a conjugate prior: ridge regression with the shrinkage
 * stated as a prior instead of as a penalty, which turns it into a posterior
 * one can ask questions of.
 * <p>
 * The model is {@code y = alpha + X beta + e} with {@code e ~ N(0, sigma^2)},
 * a normal prior {@code beta | sigma^2 ~ N(0, sigma^2 / lambda)} on the
 * standardized coefficients, and an inverse gamma prior on {@code sigma^2}. The
 * intercept is given a flat prior and integrated out, which is what centering
 * the response amounts to.
 * <p>
 * <b>The posterior mean is the ridge estimate.</b> Not approximately -- it is
 * the same vector from the same solve, so
 * {@code math.linalg.Ridge.estimate(X, y, lambda).beta} and {@link #mean} agree
 * to round-off. What the prior adds is everything around it:
 * <ul>
 * <li><b>Credible intervals for a shrunk estimate.</b> {@code Ridge} reports no
 * intervals on purpose -- its javadoc says an interval from its standard errors
 * would cover {@code E[beta(lambda)]} rather than the true coefficient, because
 * the estimator is biased. A posterior interval has no such problem: the
 * shrinkage is a prior that was stated, not a bias that was not.</li>
 * <li><b>A way to choose {@code lambda}.</b> {@code Ridge} takes it as an
 * argument and nothing in this library could produce one.
 * {@link #selectPenalty(DMatrix, DMatrix)} maximizes the marginal likelihood,
 * which needs no cross-validation and no held-out data.</li>
 * <li><b>A predictive distribution</b> that carries the uncertainty of the
 * coefficients and of {@code sigma^2}, as a Student t.</li>
 * </ul>
 * <p>
 * <b>Conventions are {@code Ridge}'s throughout</b>, and a design matrix that
 * carries an intercept column will be refused by the standardization: the
 * columns are centered and scaled to unit root mean square before the fit, the
 * coefficients come back in the original scale, and the intercept is estimated
 * outside the prior. A prior on the coefficients is no more scale invariant
 * than a penalty is.
 * <p>
 * Instances are immutable and can be shared between threads.
 *
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Bayesian_linear_regression">Wikipedia
 *      Bayesian linear regression</a>
 * @since 1.5.3
 */
public final class BayesianLinearRegression {

    /** An interval and the posterior probability it was built to hold. */
    public static final class Interval {

        /** The lower end. */
        public final double lower;
        /** The upper end. */
        public final double upper;
        /** The posterior mass between them, as asked for. */
        public final double level;

        Interval(double lower, double upper, double level) {
            this.lower = lower;
            this.upper = upper;
            this.level = level;
        }

        /**
         * The width of the interval.
         *
         * @return {@code upper - lower}
         */
        public double width() {
            return upper - lower;
        }

        @Override
        public String toString() {
            return "[" + lower + ", " + upper + "] at " + level;
        }
    }

    private final int n;
    private final int p;
    private final double lambda;
    private final double[] xBar;
    private final double[] scale;
    @SuppressWarnings("unused")
    private final double yBar;

    /** Posterior mean in the original scale, without the intercept. */
    private final double[] beta;
    private final double intercept;
    /** {@code A^-1} in the original scale, column-major, to be scaled by sigma^2. */
    private final double[] aInverse;
    private final double shape;
    private final double rate;
    private final double effectiveDf;
    private final double logEvidence;
    private final boolean converged;

    /**
     * The posterior of a linear regression under a ridge prior of strength
     * {@code lambda} and the reference prior {@code p(sigma^2) ~ 1/sigma^2}.
     *
     * @param X
     *            the design matrix, {@code n x p}, <b>without</b> an intercept
     *            column and with no constant column. Not modified
     * @param y
     *            the response, {@code n x 1}. Not modified
     * @param lambda
     *            the prior precision of the standardized coefficients, the same
     *            number {@code math.linalg.Ridge} takes as a penalty. Not
     *            negative; see {@link #logEvidence()} for what {@code 0} means
     * @return the posterior
     * @throws IllegalArgumentException
     *             if either matrix is {@code null}, if their row counts differ,
     *             if {@code y} is not a single column, if {@code X} has more
     *             columns than rows or holds a constant column, if
     *             {@code lambda} is negative or not finite, or if
     *             {@code lambda} is {@code 0} and the design is rank deficient
     */
    public static BayesianLinearRegression of(DMatrix X, DMatrix y, double lambda) {
        return new BayesianLinearRegression(X, y, lambda, 0.0, 0.0);
    }

    /**
     * The posterior under a ridge prior of strength {@code lambda} and a proper
     * inverse gamma prior on {@code sigma^2}.
     *
     * @param X
     *            the design matrix, {@code n x p}, without an intercept column.
     *            Not modified
     * @param y
     *            the response, {@code n x 1}. Not modified
     * @param lambda
     *            the prior precision of the standardized coefficients
     * @param a0
     *            the shape of the inverse gamma prior on {@code sigma^2},
     *            strictly positive
     * @param s0
     *            its scale, strictly positive
     * @return the posterior
     * @throws IllegalArgumentException
     *             under the conditions {@link #of(DMatrix, DMatrix, double)}
     *             states, and additionally if {@code a0} or {@code s0} is not
     *             finite and strictly positive
     */
    public static BayesianLinearRegression of(DMatrix X, DMatrix y, double lambda, double a0, double s0) {
        if (!(a0 > 0.0) || Double.isInfinite(a0)) {
            throw new IllegalArgumentException("a0 must be finite and strictly positive : " + a0);
        }
        if (!(s0 > 0.0) || Double.isInfinite(s0)) {
            throw new IllegalArgumentException("s0 must be finite and strictly positive : " + s0);
        }
        return new BayesianLinearRegression(X, y, lambda, a0, s0);
    }

    private BayesianLinearRegression(DMatrix X, DMatrix y, double lambda, double a0, double s0) {
        checkArguments(X, y, lambda);
        int rows = X.numRows();
        int cols = X.numColumns();
        double[] xs = X.getArrayUnsafe();
        double[] ys = y.getArrayUnsafe();

        ScaledDesign std = ScaledDesign.of(xs, ys, rows, cols, null);
        // decomposeInPlace consumes what it is given, and the standardized
        // design is wanted afterwards for nothing, so it goes straight in --
        // this is what Ridge does
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decomposeInPlace(std.x, rows, cols);
        if (lambda == 0.0 && SvdLeastSquares.rankDeficientAt(svd) >= 0) {
            throw new IllegalArgumentException("X is rank deficient and lambda is 0, so the flat prior on the "
                    + "coefficients is improper in a direction the data does not identify; use a positive lambda");
        }

        double[] scaled = SvdLeastSquares.solve(svd, std.y, lambda);
        this.n = rows;
        this.p = cols;
        this.lambda = lambda;
        this.xBar = std.xBar.clone();
        this.scale = std.scale.clone();
        this.yBar = std.yBar;
        this.beta = std.unscale(scaled);
        this.intercept = std.intercept(beta);
        this.converged = svd.converged;
        this.effectiveDf = SvdLeastSquares.effectiveDf(svd, lambda);

        double rss = 0.0;
        for (int i = 0; i < rows; i++) {
            double fit = intercept;
            for (int j = 0; j < cols; j++) {
                fit += beta[j] * xs[j * rows + i];
            }
            double r = ys[i] - fit;
            rss += r * r;
        }
        double penaltyTerm = 0.0;
        for (int j = 0; j < cols; j++) {
            penaltyTerm += scaled[j] * scaled[j];
        }

        // The degrees of freedom differ between the two priors, and not by a
        // rounding: a proper prior on beta contributes its own power of sigma^2
        // which cancels the one the integration over beta produces, while a
        // flat prior does not. So lambda == 0 is the *reference analysis* and
        // not the limit of lambda -> 0, and it is the one that reproduces
        // ordinary least squares exactly
        this.shape = a0 + (lambda == 0.0 ? (rows - 1 - cols) : (rows - 1)) / 2.0;
        this.rate = s0 + 0.5 * (rss + lambda * penaltyTerm);
        if (!(shape > 0.0)) {
            throw new IllegalArgumentException("there are not enough observations for " + cols
                    + " coefficients and an intercept: the posterior of sigma squared has shape " + shape);
        }

        this.aInverse = posteriorPrecisionInverse(svd, lambda, scale, cols);
        this.logEvidence = logEvidence(svd, rows, cols, lambda, a0, s0, rate, shape);
    }

    // ------------------------------------------------------------------
    // the posterior
    // ------------------------------------------------------------------

    /**
     * The number of coefficients, not counting the intercept.
     *
     * @return {@code p}
     */
    public int coefficients() {
        return p;
    }

    /**
     * The prior precision this posterior was formed under.
     *
     * @return {@code lambda}
     */
    public double penalty() {
        return lambda;
    }

    /**
     * The intercept, estimated outside the prior as
     * {@code yBar - sum(beta_j xBar_j)}.
     *
     * @return the intercept
     */
    public double intercept() {
        return intercept;
    }

    /**
     * Writes the posterior mean of the coefficients, in the original scale of
     * the columns, into {@code out}.
     * <p>
     * This is the ridge estimate at the same {@code lambda}, to round-off.
     *
     * @param out
     *            where the result is written, one entry per coefficient. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or is not as long as there are
     *             coefficients
     */
    public void mean(double[] out) {
        if (out == null || out.length != p) {
            throw new IllegalArgumentException("out must hold " + p + " coefficients");
        }
        System.arraycopy(beta, 0, out, 0, p);
    }

    /**
     * The posterior covariance of the coefficients <b>given
     * {@code sigma^2 == 1}</b>, that is {@code (Z'Z + lambda I)^-1} carried back
     * into the original scale.
     * <p>
     * Multiply by a variance to get a covariance, or use
     * {@link #given(double)}, which does. It is <em>not</em> the sampling
     * variance of the ridge estimator that
     * {@code SvdLeastSquares.varianceMatrix} returns: that is
     * {@code V diag(d^2/(d^2+lambda)^2) V'} and this is
     * {@code V diag(1/(d^2+lambda)) V'}. The two coincide only at
     * {@code lambda == 0}.
     *
     * @return a fresh {@code p x p} matrix, exactly symmetric
     */
    public DMatrix covariance() {
        DMatrix m = new DMatrix(p, p);
        double[] a = m.getArrayUnsafe();
        System.arraycopy(aInverse, 0, a, 0, aInverse.length);
        return m;
    }

    /**
     * The posterior of the error variance, an inverse gamma.
     * <p>
     * Its mean is {@code rate / (shape - 1)} where a least squares fit would
     * report {@code rss / (n - p - 1)}; the two agree at {@code lambda == 0} up
     * to exactly that factor, which is what an inverse gamma mean is.
     *
     * @return the posterior of {@code sigma^2}
     */
    public InverseGamma variance() {
        return new InverseGamma(shape, rate);
    }

    /**
     * The posterior of the coefficients for a known error variance, which is a
     * multivariate normal.
     * <p>
     * With {@code sigma^2} unknown the joint posterior is a multivariate
     * Student t, which this library has no class for; what it does have are the
     * one dimensional marginals, {@link #marginal(int)}, and this conditional.
     *
     * @param sigmaSquared
     *            the error variance to condition on, strictly positive
     * @return the conditional posterior of the coefficients
     * @throws IllegalArgumentException
     *             if {@code sigmaSquared} is not finite and strictly positive
     */
    public MultivariateNormal given(double sigmaSquared) {
        if (!(sigmaSquared > 0.0) || Double.isInfinite(sigmaSquared)) {
            throw new IllegalArgumentException("sigmaSquared must be finite and strictly positive : "
                    + sigmaSquared);
        }
        DMatrix cov = covariance();
        cov.scaleInplace(sigmaSquared);
        return new MultivariateNormal(beta, cov);
    }

    /**
     * The effective number of parameters the fit used,
     * {@code sum d_i^2 / (d_i^2 + lambda)}, which falls from {@code p} at no
     * shrinkage towards zero as the prior tightens.
     *
     * @return the effective degrees of freedom
     */
    public double effectiveDegreesOfFreedom() {
        return effectiveDf;
    }

    /**
     * Whether the underlying decomposition converged.
     *
     * @return {@code false} if the Jacobi sweeps were exhausted, in which case
     *         everything here may be inaccurate
     */
    public boolean converged() {
        return converged;
    }

    // ------------------------------------------------------------------
    // marginals, intervals, prediction
    // ------------------------------------------------------------------

    /**
     * The degrees of freedom of every coefficient marginal and of the
     * predictive, {@code 2 * shape}.
     *
     * @return the degrees of freedom
     */
    public double marginalDegreesOfFreedom() {
        return 2.0 * shape;
    }

    /**
     * The marginal posterior of one coefficient, as the <b>standard</b> Student
     * t it is once centered and scaled.
     * <p>
     * {@link StudentT} carries no location and no scale, so the marginal of
     * coefficient {@code j} is
     * {@code mean[j] + coefficientScale(j) * marginal(j)}. For an interval,
     * {@link #credibleInterval(int, double)} does that arithmetic.
     *
     * @param j
     *            the coefficient
     * @return the standard Student t on {@link #marginalDegreesOfFreedom()}
     * @throws IndexOutOfBoundsException
     *             if {@code j} is not a coefficient
     */
    public StudentT marginal(int j) {
        checkCoefficient(j);
        return new StudentT(marginalDegreesOfFreedom());
    }

    /**
     * The scale of the marginal posterior of one coefficient,
     * {@code sqrt(rate / shape * (A^-1)_jj)}.
     *
     * @param j
     *            the coefficient
     * @return the scale, strictly positive
     * @throws IndexOutOfBoundsException
     *             if {@code j} is not a coefficient
     */
    public double coefficientScale(int j) {
        checkCoefficient(j);
        return Math.sqrt(rate / shape * aInverse[j * p + j]);
    }

    /**
     * The equal-tailed credible interval for one coefficient.
     * <p>
     * Unlike a confidence interval built from a ridge standard error, this one
     * is honest about a shrunk estimate: it is a statement about the
     * coefficient given the prior, not about the expectation of an estimator.
     *
     * @param j
     *            the coefficient
     * @param level
     *            the posterior mass the interval is to hold, in {@code (0, 1)}
     * @return the interval
     * @throws IndexOutOfBoundsException
     *             if {@code j} is not a coefficient
     * @throws IllegalArgumentException
     *             if {@code level} is not in {@code (0, 1)}
     */
    public Interval credibleInterval(int j, double level) {
        checkCoefficient(j);
        checkLevel(level);
        double half = new StudentT(marginalDegreesOfFreedom()).inverseCdf(0.5 * (1.0 + level))
                * coefficientScale(j);
        return new Interval(beta[j] - half, beta[j] + half, level);
    }

    /**
     * The mean of the posterior predictive at a new point,
     * {@code intercept + x . beta}.
     *
     * @param x
     *            the new point, one entry per column of the design. Not
     *            modified
     * @return the predictive mean
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is not as long as there are
     *             coefficients
     */
    public double predictiveMean(double[] x) {
        checkPoint(x);
        double sum = intercept;
        for (int j = 0; j < p; j++) {
            sum += beta[j] * x[j];
        }
        return sum;
    }

    /**
     * The scale of the posterior predictive at a new point.
     * <p>
     * It carries three things: the error variance, the uncertainty of the
     * coefficients through {@code z' A^-1 z}, and the uncertainty of the
     * intercept through {@code 1/n}.
     *
     * @param x
     *            the new point, one entry per column of the design. Not
     *            modified
     * @return the predictive scale, strictly positive
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is not as long as there are
     *             coefficients
     */
    public double predictiveScale(double[] x) {
        checkPoint(x);
        double[] z = new double[p];
        for (int j = 0; j < p; j++) {
            z[j] = (x[j] - xBar[j]) / scale[j];
        }
        // z' A^-1 z, with A^-1 held in the original scale, so the scaling is
        // undone on the way in
        double quadratic = 0.0;
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                quadratic += z[i] * (aInverse[j * p + i] * scale[i] * scale[j]) * z[j];
            }
        }
        return Math.sqrt(rate / shape * (1.0 + 1.0 / n + quadratic));
    }

    /**
     * The posterior predictive at a new point, as the standard Student t it is
     * once centered on {@link #predictiveMean(double[])} and scaled by
     * {@link #predictiveScale(double[])}.
     *
     * @param x
     *            the new point. Not modified
     * @return the standard Student t on {@link #marginalDegreesOfFreedom()}
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or is of the wrong length
     */
    public StudentT predictive(double[] x) {
        checkPoint(x);
        return new StudentT(marginalDegreesOfFreedom());
    }

    /**
     * The equal-tailed prediction interval for a new observation at {@code x}.
     *
     * @param x
     *            the new point. Not modified
     * @param level
     *            the predictive mass the interval is to hold, in {@code (0, 1)}
     * @return the interval
     * @throws IllegalArgumentException
     *             if {@code x} is of the wrong length or {@code level} is not
     *             in {@code (0, 1)}
     */
    public Interval predictionInterval(double[] x, double level) {
        checkLevel(level);
        double center = predictiveMean(x);
        double half = new StudentT(marginalDegreesOfFreedom()).inverseCdf(0.5 * (1.0 + level))
                * predictiveScale(x);
        return new Interval(center - half, center + half, level);
    }

    // ------------------------------------------------------------------
    // the evidence
    // ------------------------------------------------------------------

    /**
     * The log marginal likelihood of the response under this prior.
     * <p>
     * <b>It is relative, and by construction.</b> The intercept carries a flat
     * prior, which has no normalizing constant, and under
     * {@link #of(DMatrix, DMatrix, double)} so does {@code sigma^2}. What the
     * missing constant depends on is the number of observations and nothing
     * else -- not on {@code lambda}, not on the design -- so it cancels in
     * <b>any comparison of models fitted to the same response</b>, which is the
     * case a Bayes factor is for, and in the maximization
     * {@link #selectPenalty(DMatrix, DMatrix)} performs. It does not cancel
     * between different responses, and no number here should be read as an
     * absolute probability.
     * <p>
     * At {@code lambda == 0} it is {@link Double#NaN}: a flat prior on the
     * coefficients leaves a constant that depends on how many of them there
     * are, so the one thing a Bayes factor would be used for -- comparing
     * designs of different width -- is exactly what it cannot do.
     *
     * @return the log marginal likelihood up to a constant depending only on
     *         the number of observations, or {@link Double#NaN} at
     *         {@code lambda == 0}
     */
    public double logEvidence() {
        return logEvidence;
    }

    /**
     * The log marginal likelihood at one penalty, under the reference prior on
     * {@code sigma^2}, without building a posterior.
     *
     * @param X
     *            the design matrix, without an intercept column. Not modified
     * @param y
     *            the response. Not modified
     * @param lambda
     *            the prior precision, strictly positive
     * @return the log marginal likelihood, relative in the sense
     *         {@link #logEvidence()} describes
     * @throws IllegalArgumentException
     *             under the conditions {@link #of(DMatrix, DMatrix, double)}
     *             states
     */
    public static double logEvidence(DMatrix X, DMatrix y, double lambda) {
        return of(X, y, lambda).logEvidence();
    }

    /**
     * The log marginal likelihood at one penalty under a proper prior on
     * {@code sigma^2}.
     *
     * @param X
     *            the design matrix, without an intercept column. Not modified
     * @param y
     *            the response. Not modified
     * @param lambda
     *            the prior precision, strictly positive
     * @param a0
     *            the shape of the inverse gamma prior, strictly positive
     * @param s0
     *            its scale, strictly positive
     * @return the log marginal likelihood
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #of(DMatrix, DMatrix, double, double, double)} states
     */
    public static double logEvidence(DMatrix X, DMatrix y, double lambda, double a0, double s0) {
        return of(X, y, lambda, a0, s0).logEvidence();
    }

    /**
     * The penalty that maximizes the marginal likelihood -- the answer
     * {@code math.linalg.Ridge} needs and cannot produce, arrived at without
     * cross-validation and without held-out data.
     * <p>
     * The search runs over {@code log(lambda)} and not over {@code lambda},
     * because the evidence is flat in the one and curved in the other, and
     * because the squared singular values of the design -- which is what
     * {@code lambda} is compared against -- span orders of magnitude. The
     * bracket starts from the mean of those squared singular values, since
     * below the smallest of them the prior touches nothing and above the
     * largest it annihilates everything.
     *
     * @param X
     *            the design matrix, without an intercept column. Not modified
     * @param y
     *            the response. Not modified
     * @return the maximizing penalty, strictly positive
     * @throws IllegalArgumentException
     *             under the conditions {@link #of(DMatrix, DMatrix, double)}
     *             states, or if no maximum can be bracketed
     */
    public static double selectPenalty(DMatrix X, DMatrix y) {
        return selectPenalty(X, y, 0.0, 0.0);
    }

    /**
     * The penalty that maximizes the marginal likelihood under a proper prior
     * on {@code sigma^2}.
     *
     * @param X
     *            the design matrix, without an intercept column. Not modified
     * @param y
     *            the response. Not modified
     * @param a0
     *            the shape of the inverse gamma prior, strictly positive
     * @param s0
     *            its scale, strictly positive
     * @return the maximizing penalty, strictly positive
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #of(DMatrix, DMatrix, double, double, double)} states,
     *             or if no maximum can be bracketed
     */
    public static double selectPenalty(DMatrix X, DMatrix y, double a0, double s0) {
        boolean proper = a0 > 0.0 || s0 > 0.0;
        if (proper) {
            if (!(a0 > 0.0) || Double.isInfinite(a0)) {
                throw new IllegalArgumentException("a0 must be finite and strictly positive : " + a0);
            }
            if (!(s0 > 0.0) || Double.isInfinite(s0)) {
                throw new IllegalArgumentException("s0 must be finite and strictly positive : " + s0);
            }
        }
        checkArguments(X, y, 1.0);
        int rows = X.numRows();
        int cols = X.numColumns();
        // one decomposition for the whole search: the evidence depends on the
        // design only through the singular values and the fitted residual, and
        // both come from it
        ScaledDesign std = ScaledDesign.of(X.getArrayUnsafe(), y.getArrayUnsafe(), rows, cols, null);
        double[] centered = std.y.clone();
        final FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decomposeInPlace(std.x, rows, cols);

        double meanSquare = 0.0;
        for (int i = 0; i < cols; i++) {
            meanSquare += svd.sigma[i] * svd.sigma[i];
        }
        meanSquare /= cols;

        final double[] yc = centered;
        final int nn = rows;
        final int pp = cols;
        final double aa = a0;
        final double ss = s0;
        DFunction negative = new DFunction() {
            @Override
            public double apply(double logLambda) {
                double value = evidenceAt(svd, yc, nn, pp, Math.exp(logLambda), aa, ss);
                return Double.isNaN(value) ? Double.MAX_VALUE / 8.0 : -value;
            }
        };
        BrentMinimizer minimizer = new BrentMinimizer();
        double start = Math.log(meanSquare);
        BrentMinimizer.Bracket bracket = minimizer.bracket(negative, start, start + 1.0);
        if (!bracket.bracketed) {
            throw new IllegalArgumentException(
                    "the marginal likelihood has no interior maximum in lambda for this design; "
                            + "it is monotone, which means the data prefers no shrinkage or total shrinkage");
        }
        return Math.exp(minimizer.minimize(negative, bracket).x);
    }

    // ------------------------------------------------------------------
    // the arithmetic
    // ------------------------------------------------------------------

    /**
     * {@code (Z'Z + lambda I)^-1 = V diag(1/(d_i^2 + lambda)) V'}, carried into
     * the original scale by dividing row {@code i} and column {@code j} by
     * their scales.
     */
    private static double[] posteriorPrecisionInverse(FlatParallelJacobiSVD.Result svd, double lambda,
            double[] scale, int p) {
        double[] a = new double[p * p];
        for (int k = 0; k < p; k++) {
            double d = svd.sigma[k];
            double f = 1.0 / (d * d + lambda);
            for (int r = 0; r < p; r++) {
                double vr = svd.V[k * p + r] * f;
                for (int c = 0; c < p; c++) {
                    a[c * p + r] += vr * svd.V[k * p + c];
                }
            }
        }
        for (int c = 0; c < p; c++) {
            for (int r = 0; r < p; r++) {
                a[c * p + r] /= scale[r] * scale[c];
            }
        }
        // exactly symmetric, so that MultivariateNormal and CholeskyDecomp are
        // handed what they expect rather than something symmetric to round-off
        for (int c = 0; c < p; c++) {
            for (int r = c + 1; r < p; r++) {
                double m = 0.5 * (a[c * p + r] + a[r * p + c]);
                a[c * p + r] = m;
                a[r * p + c] = m;
            }
        }
        return a;
    }

    /** The evidence of a posterior that has already been formed. */
    private static double logEvidence(FlatParallelJacobiSVD.Result svd, int n, int p, double lambda, double a0,
            double s0, double rate, double shape) {
        if (lambda == 0.0) {
            return Double.NaN;
        }
        double logDet = 0.0;
        for (int i = 0; i < p; i++) {
            double d = svd.sigma[i];
            logDet += Math.log(d * d + lambda);
        }
        double value = 0.5 * p * Math.log(lambda) - 0.5 * logDet - 0.5 * (n - 1) * Math.log(2.0 * Math.PI)
                - 0.5 * Math.log(n) - shape * Math.log(rate) + FastGamma.logGamma(shape);
        if (a0 > 0.0) {
            value += a0 * Math.log(s0) - FastGamma.logGamma(a0);
        }
        return value;
    }

    /** The evidence at one penalty, from a decomposition computed once. */
    private static double evidenceAt(FlatParallelJacobiSVD.Result svd, double[] centeredY, int n, int p,
            double lambda, double a0, double s0) {
        if (!(lambda > 0.0) || Double.isInfinite(lambda)) {
            return Double.NaN;
        }
        double[] scaled = SvdLeastSquares.solve(svd, centeredY, lambda);
        // the residual of the standardized fit, which is the residual of the
        // original one because the transform is a change of coordinates
        double rss = residualSumOfSquares(svd, centeredY, scaled, n, p);
        double penalty = 0.0;
        for (int j = 0; j < p; j++) {
            penalty += scaled[j] * scaled[j];
        }
        double shape = a0 + (n - 1) / 2.0;
        double rate = s0 + 0.5 * (rss + lambda * penalty);
        return logEvidence(svd, n, p, lambda, a0, s0, rate, shape);
    }

    /**
     * {@code ||yc - Z b||^2} where {@code Z = U D V'}, formed from the
     * decomposition because the standardized design was consumed by it.
     */
    private static double residualSumOfSquares(FlatParallelJacobiSVD.Result svd, double[] centeredY,
            double[] scaled, int n, int p) {
        // Z b = U (D (V' b)), and U is n x p
        double[] vtb = new double[p];
        for (int k = 0; k < p; k++) {
            double sum = 0.0;
            for (int j = 0; j < p; j++) {
                sum += svd.V[k * p + j] * scaled[j];
            }
            vtb[k] = sum * svd.sigma[k];
        }
        double rss = 0.0;
        for (int i = 0; i < n; i++) {
            double fit = 0.0;
            for (int k = 0; k < p; k++) {
                fit += svd.U[k * n + i] * vtb[k];
            }
            double r = centeredY[i] - fit;
            rss += r * r;
        }
        return rss;
    }

    // ------------------------------------------------------------------
    // validation
    // ------------------------------------------------------------------

    private static void checkArguments(DMatrix X, DMatrix y, double lambda) {
        if (X == null || y == null) {
            throw new IllegalArgumentException("X and y must not be null");
        }
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
        if (!(lambda >= 0.0) || Double.isInfinite(lambda)) {
            throw new IllegalArgumentException("lambda must be finite and non-negative : " + lambda);
        }
    }

    private void checkCoefficient(int j) {
        if (j < 0 || j >= p) {
            throw new IndexOutOfBoundsException("no such coefficient : " + j);
        }
    }

    private void checkPoint(double[] x) {
        if (x == null || x.length != p) {
            throw new IllegalArgumentException("x must hold " + p + " coordinates");
        }
    }

    private static void checkLevel(double level) {
        if (Double.isNaN(level) || !(level > 0.0) || !(level < 1.0)) {
            throw new IllegalArgumentException("level must be in (0, 1) : " + level);
        }
    }
}
