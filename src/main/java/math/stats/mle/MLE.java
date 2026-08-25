/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * Derived from SSJ.
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 */
package math.stats.mle;

import java.util.Arrays;

import math.MathConsts;
import math.cern.FastGamma;
import math.cern.GammaFun;
import math.fun.DFunction;
import math.fun.DMultiFunctionEval;
import math.fun.NumericalDiffDMultiFunction;
import math.optim.BrentMinimizer;
import math.optim.CGOptimizer;
import math.solve.RootFinder;

/**
 * Provides methods for maximum likelihood estimation of distribution
 * parameters.
 *
 * @since 1.5.3
 */
public final class MLE {

    private static final double HUGE = 1.0e200;
    /** Iterations the Student t fit may take. */
    private static final int T_MAX_ITER = 500;
    /**
     * Relative log-likelihood gain below which the Student t fit has settled.
     * ECME never lowers the likelihood, so the gain is what is left to win.
     */
    private static final double T_TOL = 1.0e-12;
    /** Degrees of freedom the Student t fit starts from. */
    private static final double T_DF_START = 4.0;
    /** Floor of the search for the degrees of freedom. */
    private static final double T_DF_MIN = 1.0e-2;
    /**
     * Above this a t and its normal limit are the same distribution. Telling
     * them apart means resolving an excess kurtosis of {@code 6 / (nu - 4)}
     * against a standard error of {@code sqrt(24 / n)}, which at
     * {@code nu = 1e4} would take more than {@code 1e8} observations. A
     * maximizer beyond this point carries no information, so it is reported as
     * {@link Double#POSITIVE_INFINITY} instead of as a number.
     */
    private static final double T_DF_MAX = 1.0e4;
    /** Scales the median absolute deviation into a standard deviation. */
    private static final double MAD_TO_SIGMA = 1.4826022185056018;
    /**
     * Maximizes the profile likelihood in {@code log(nu)}. The bracketing is
     * kept short on purpose: the caller has already established that a
     * finite maximum exists, so a long walk means something is wrong.
     */
    private static final BrentMinimizer DF_MINIMIZER = new BrentMinimizer(1.0e-8, 200, 60);
    private static final String NO_OBS_MSG = "No observations (x[].length = 0)";
    /** How far {@link #bracketRoot} may double its half width. */
    private static final int MAX_EXPANSIONS = 60;
    /** Shape parameters are strictly positive; this is the floor of the search. */
    private static final double MIN_POSITIVE_SHAPE = 1.0e-10;
    /**
     * Largest score residual the Beta likelihood solve accepts as a solution.
     * The residual is a difference of digamma values and therefore
     * dimensionless; measured over 5400 samples the solver reaches 1.3e-14.
     */
    private static final double BETA_TOL = 1.0e-12;
    /** Newton iterations the Beta likelihood solve may take. */
    private static final int BETA_MAX_ITER = 200;
    /** Step halvings per iteration before the Beta solve gives up. */
    private static final int BETA_MAX_HALVINGS = 60;

    private static final class GammaMLE implements DFunction {
        private final int n;
        private final double empiricalMean;
        private final double sumLn;

        GammaMLE(int n, double empiricalMean, double sumLn) {
            this.n = n;
            this.empiricalMean = empiricalMean;
            this.sumLn = sumLn;
        }

        @Override
        public double apply(double x) {
            if (x <= 0.0) {
                return HUGE;
            }
            return (n * Math.log(empiricalMean / x) + n * GammaFun.digamma(x) - sumLn);
        }
    }

    private static final class WeibullMLE implements DFunction {
        private final double xi[];
        private final double lnXi[];
        private double sumLnXi = 0.0;

        WeibullMLE(double x[]) {
            xi = x.clone();
            lnXi = new double[x.length];

            for (int i = 0; i < x.length; i++) {
                // the public entry points reject anything outside (0, infinity)
                double lnx = Math.log(x[i]);
                lnXi[i] = lnx;
                sumLnXi += lnx;
            }
        }

        @Override
        public double apply(double x) {
            if (x <= 0.0) {
                return HUGE;
            }
            double sumXiLnXi = 0.0;
            double sumXi = 0.0;
            for (int i = 0; i < xi.length; i++) {
                double xalpha = Math.pow(xi[i], x);
                sumXiLnXi += xalpha * lnXi[i];
                sumXi += xalpha;
            }
            return x * (xi.length * sumXiLnXi - sumLnXi * sumXi) - (xi.length * sumXi);
        }
    }

    private static final class ChiSquareMLE extends NumericalDiffDMultiFunction {
        private static final double TERM = MathConsts.LN_2 / 2.0;
        private final double sumLnHalfth;
        private final int n;

        ChiSquareMLE(double sumLn, int n) {
            this.sumLnHalfth = sumLn / 2.0;
            this.n = n;
        }

        @Override
        public double apply(double[] point) {
            double x = point[0];
            if (x <= 0.0) {
                return -HUGE;
            }
            return x * sumLnHalfth - n * FastGamma.logGamma(x / 2.0) - (n * x) * TERM;
        }
    }

    /**
     * Estimates the parameters {@code shape} ({@code k}) and {@code scale}
     * (&theta;) of the Gamma distribution from the observations {@code x} using
     * the maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters {@code k} and &theta;
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value that is not finite and strictly positive
     */
    public static ParGamma getGammaMLE(double[] x) {
        int n = requirePositive(x);
        requireAtLeast(n, 2);
        double sum = 0.0;
        double sumLn = 0.0;

        for (int i = 0; i < n; i++) {
            sum += x[i];
            sumLn += Math.log(x[i]);
        }
        double empiricalMean = sum / (double) n;

        sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += (x[i] - empiricalMean) * (x[i] - empiricalMean);
        }

        double alphaMME = (empiricalMean * empiricalMean * (double) n) / sum;
        GammaMLE equation = new GammaMLE(n, empiricalMean, sumLn);
        // the moment estimate is only a guess at where the root is; a fixed
        // interval around it misses once the shape outgrows the interval
        double[] bracket = bracketRoot(equation, alphaMME);
        if (bracket == null) {
            throw new IllegalArgumentException(
                    "no Gamma shape solves the likelihood equation for this sample (moment estimate "
                            + alphaMME + ")");
        }

        // an absolute tolerance is the wrong shape of request: 1e-7 at a shape
        // of 0.1 is a relative accuracy of 1e-6. Ask for everything the solver
        // has and let its own relative term set the floor
        double shape = RootFinder.brentDekker(bracket[0], bracket[1], equation, MathConsts.MIN_TOL);
        return new ParGamma(shape, empiricalMean / shape);
    }

    /**
     * Estimates the parameters &mu; and &sigma; of the LogNormal distribution
     * from the observations {@code x} using the maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters &mu; and &sigma;
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value that is not finite and strictly positive
     */
    public static ParLogNormal getLogNormalMLE(double[] x) {
        int n = requirePositive(x);
        requireAtLeast(n, 2);
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += Math.log(x[i]);
        }

        double mu_hat = sum / n;
        double tmp;
        sum = 0.0;

        for (int i = 0; i < n; i++) {
            tmp = Math.log(x[i]) - mu_hat;
            sum += (tmp * tmp);
        }

        return new ParLogNormal(mu_hat, Math.sqrt(sum / n));
    }

    /**
     * Estimates the parameters {@code scale} (&lambda;) and {@code shape}
     * ({@code k}) of the Weibull distribution from the observations {@code x}
     * using the maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters &lambda; and {@code k}
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value that is not finite and strictly positive
     */
    public static ParWeibull getWeibullMLE(double[] x) {
        int n = requirePositive(x);
        requireAtLeast(n, 2);
        double sumLn = 0.0;
        double sumLn2 = 0.0;

        for (int i = 0; i < x.length; i++) {
            double lnxi = Math.log(x[i]);
            sumLn += lnxi;
            sumLn2 += (lnxi * lnxi);
        }

        double alpha0 = Math.sqrt((double) n / ((6.0 / MathConsts.PI_SQUARED) * (sumLn2 - sumLn * sumLn / (double) n)));
        WeibullMLE equation = new WeibullMLE(x);
        // as for the Gamma shape: search for the bracket, do not assume one
        double[] bracket = bracketRoot(equation, alpha0);
        if (bracket == null) {
            throw new IllegalArgumentException(
                    "no Weibull shape solves the likelihood equation for this sample (initial guess "
                            + alpha0 + ")");
        }

        double k = RootFinder.brentDekker(bracket[0], bracket[1], equation, MathConsts.MIN_TOL);

        double sumXalpha = 0.0;
        for (int i = 0; i < x.length; i++) {
            sumXalpha += Math.pow(x[i], k);
        }
        double scale = 1.0 / (Math.pow((double) n / sumXalpha, 1.0 / k));

        return new ParWeibull(scale, k);
    }

    /**
     * Estimates the location &mu;, the scale &sigma; and the degrees of freedom
     * &nu; of a StudentT distribution from the observations {@code x} using the
     * maximum likelihood method. The degrees of freedom are double-valued.
     * <p>
     * The fit is the ECME of Liu and Rubin, which treats the t as a normal
     * whose variance is itself random: given &nu;, the weights
     * {@code w_i = (nu + 1) / (nu + d_i^2)} turn the estimate of &mu; and
     * &sigma; into a weighted mean and a weighted variance, and given those,
     * &nu; maximizes the observed-data likelihood along one line. Working in
     * {@code log(nu)} keeps the parameter positive and makes the step scale
     * with the value.
     * <p>
     * The search starts from the median and the median absolute deviation
     * rather than from the mean and the second moment. Inverting
     * {@code E[X^2] = nu / (nu - 2)} looks like a cheap way to a starting
     * value, but the map flattens as &nu; grows -- {@code dnu/dm = -nu^2/2},
     * so a one per cent error in the second moment is sixteen per cent in
     * &nu; at &nu; = 30 -- and once the sample second moment drops below one,
     * which is routine above &nu; = 25, the starting value turns negative.
     * <p>
     * The returned {@link ParStudentT#df} is
     * {@link Double#POSITIVE_INFINITY} when the likelihood is still rising
     * where the t and the normal have become the same distribution. That is
     * not a failure of the search: the degrees of freedom of a t are weakly
     * identified once the tails are light, and for such a sample there is no
     * finite maximizer. {@link ParStudentT#location} and
     * {@link ParStudentT#scale} are then the normal estimates.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return the location, the scale and the degrees of freedom
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value that is not finite
     */
    public static ParStudentT getStudentTMLE(double[] x) {
        int n = requireFinite(x);
        requireAtLeast(n, 2);
        return fitT(x, n, Double.NaN);
    }

    /**
     * Estimates the location {@code x0} and the scale &gamma; of the Cauchy
     * distribution from the observations {@code x} using the maximum
     * likelihood method.
     * <p>
     * The Cauchy distribution is the StudentT at one degree of freedom, so
     * this is {@link #getStudentTMLE(double[])} with &nu; held at {@code 1}
     * instead of estimated -- a different estimator, and the right one for a
     * caller who knows the tail. Measured against the free fit on Cauchy
     * data, it costs about 4 per cent less error in the location and 13 to 25
     * per cent less in the scale, and runs about thirty times faster, because
     * a fixed &nu; removes the one-dimensional search that otherwise runs on
     * every iteration.
     * <p>
     * Holding &nu; at some other value is deliberately not offered. The
     * location would survive it -- the weights are even, so the location
     * equation is unbiased at the centre of symmetry whatever &nu; is -- but
     * the scale would not: on Cauchy data an assumed &nu; drives the scale to
     * {@code gamma * sqrt(nu)}, an error that does not shrink with the sample.
     * <p>
     * The location likelihood of a Cauchy is the textbook multimodal one, but
     * only with the scale held fixed. The joint problem solved here has a
     * single stationary point, which is why the fit does not depend on where
     * it starts. That fails below three observations, and three is therefore
     * the minimum.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return the location and the scale
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than three
             observations, or holds a value that is not finite
     */
    public static ParCauchy getCauchyMLE(double[] x) {
        int n = requireFinite(x);
        // at two observations the starts land on different optima, and at one
        // the scale collapses; both were measured
        requireAtLeast(n, 3);
        ParStudentT fit = fitT(x, n, 1.0);
        return new ParCauchy(fit.location, fit.scale, fit.converged);
    }

    /**
     * The shared t fit. {@code fixedDf} of {@link Double#NaN} means the
     * degrees of freedom are estimated along with the rest, which is ECME;
     * any other value holds them there, which is plain EM and monotone in the
     * observed-data likelihood just the same.
     */
    private static ParStudentT fitT(double[] x, int n, double fixedDf) {
        double[] sorted = x.clone();
        Arrays.sort(sorted);
        double location = median(sorted);
        double scale = MAD_TO_SIGMA * medianAbsoluteDeviation(sorted, location);
        if (!(scale > 0.0)) {
            // more than half of the sample sits on the median; fall back to
            // the second moment, which only fails if the sample is constant
            scale = standardDeviation(x, location);
        }
        if (!(scale > 0.0)) {
            throw new IllegalArgumentException("the sample has no scale: every observation is " + x[0]);
        }

        double df = Double.isNaN(fixedDf) ? T_DF_START : fixedDf;
        double[] weight = new double[n];
        boolean converged = false;
        double previous = Double.NEGATIVE_INFINITY;

        for (int it = 0; it < T_MAX_ITER; it++) {
            // E-step: how much each observation counts under the current fit
            for (int i = 0; i < n; i++) {
                double d = (x[i] - location) / scale;
                weight[i] = (df == Double.POSITIVE_INFINITY) ? 1.0 : (df + 1.0) / (df + d * d);
            }

            // CM-steps 1 and 2: the weighted mean, then the weighted spread
            // around it
            double sumW = 0.0;
            double sumWX = 0.0;
            for (int i = 0; i < n; i++) {
                sumW += weight[i];
                sumWX += weight[i] * x[i];
            }
            double nextLocation = sumWX / sumW;
            double sumWD = 0.0;
            for (int i = 0; i < n; i++) {
                double d = x[i] - nextLocation;
                sumWD += weight[i] * d * d;
            }
            double nextScale = Math.sqrt(sumWD / n);
            if (!(nextScale > 0.0)) {
                // the fit has collapsed onto a point mass, where the
                // likelihood is unbounded. That is a degeneracy of the model,
                // not an estimate, and it leaves the loop unconverged
                location = nextLocation;
                scale = nextScale;
                break;
            }

            // CM-step 3: nu against the observed-data likelihood, not against
            // the completed one, which is what makes this ECME rather than EM.
            // A caller that fixed nu skips it, and the loop is then the EM for
            // a t of known degrees of freedom
            double nextDf = Double.isNaN(fixedDf) ? maximizeDf(x, nextLocation, nextScale, df) : fixedDf;

            location = nextLocation;
            scale = nextScale;
            df = nextDf;

            // the parameters cannot serve as the test: nu sits on a flat
            // likelihood and is found by a search of finite accuracy, so it
            // never repeats itself to the precision the others reach
            double current = tLogLikelihood(x, location, scale, df);
            if (current - previous <= T_TOL * (1.0 + Math.abs(current))) {
                converged = true;
                break;
            }
            previous = current;
        }

        if (converged && df == Double.POSITIVE_INFINITY) {
            // in the limit every weight is one, so the fixed point is the
            // sample mean and the sample standard deviation, exactly. The
            // stopping rule fires in the iteration nu goes infinite, one step
            // before those are reached
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += x[i];
            }
            location = sum / n;
            scale = standardDeviation(x, location);
        }

        return new ParStudentT(location, scale, df, converged);
    }

    /**
     * Maximizes the observed-data log-likelihood of a t over the degrees of
     * freedom, holding the location and the scale fixed, and returns
     * {@link Double#POSITIVE_INFINITY} if there is no finite maximizer.
     */
    private static double maximizeDf(final double[] x, final double location, final double scale,
            double start) {
        DFunction negative = new DFunction() {
            @Override
            public double apply(double logDf) {
                return -tLogLikelihood(x, location, scale, Math.exp(logDf));
            }
        };
        // still rising where the two distributions have already merged means
        // it rises all the way, and no number is the answer
        if (negative.apply(Math.log(T_DF_MAX)) <= negative.apply(Math.log(T_DF_MAX / 10.0))) {
            return Double.POSITIVE_INFINITY;
        }
        double from = Math.log(Math.min(Math.max(start, T_DF_MIN), T_DF_MAX));
        BrentMinimizer.Bracket bracket = DF_MINIMIZER.bracket(negative, from, from + 0.5);
        if (!bracket.bracketed) {
            return Double.POSITIVE_INFINITY;
        }
        return Math.exp(DF_MINIMIZER.minimize(negative, bracket).x);
    }

    /**
     * Log-likelihood of a t with the given location, scale and degrees of
     * freedom, and of its normal limit for {@code df} infinite.
     */
    private static double tLogLikelihood(double[] x, double location, double scale, double df) {
        int n = x.length;
        double sum = 0.0;
        if (df == Double.POSITIVE_INFINITY) {
            for (int i = 0; i < n; i++) {
                double d = (x[i] - location) / scale;
                sum += d * d;
            }
            return -n * (Math.log(scale) + 0.5 * Math.log(MathConsts.TWO_PI)) - 0.5 * sum;
        }
        double constant = FastGamma.logGamma((df + 1.0) / 2.0) - FastGamma.logGamma(df / 2.0)
                - 0.5 * Math.log(Math.PI * df) - Math.log(scale);
        for (int i = 0; i < n; i++) {
            double d = (x[i] - location) / scale;
            sum += Math.log1p(d * d / df);
        }
        return n * constant - (df + 1.0) * 0.5 * sum;
    }

    /** Median of an already sorted sample. */
    private static double median(double[] sorted) {
        int n = sorted.length;
        if ((n & 1) == 1) {
            return sorted[n / 2];
        }
        return 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    }

    /** Median of the absolute deviations from {@code center}. */
    private static double medianAbsoluteDeviation(double[] sorted, double center) {
        double[] deviations = new double[sorted.length];
        for (int i = 0; i < sorted.length; i++) {
            deviations[i] = Math.abs(sorted[i] - center);
        }
        Arrays.sort(deviations);
        return median(deviations);
    }

    /** Root mean square deviation from {@code center}. */
    private static double standardDeviation(double[] x, double center) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            double d = x[i] - center;
            sum += d * d;
        }
        return Math.sqrt(sum / x.length);
    }

    /**
     * Estimates the parameters {@code alpha} (&alpha;) and {@code beta}
     * (&beta;) of the Beta distribution from the observations {@code x} using
     * the maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters {@code alpha} (&alpha;) and {@code beta}
     *         (&beta;)
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value outside {@code (0, 1)}
     */
    public static ParBeta getBetaMLE(double[] x) {
        int n = requireUnitInterval(x);
        requireAtLeast(n, 2);
        double sum = 0.0;
        double a = 0.0;
        double b = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
            a += Math.log(x[i]);
            b += Math.log1p(-x[i]);
        }
        double mean = sum / n;

        sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += (x[i] - mean) * (x[i] - mean);
        }
        double var = sum / (n - 1);

        // the score equations are psi(alpha) - psi(alpha + beta) = mean(ln x)
        // and its counterpart in 1 - x. Against the raw sums instead of the
        // means the right hand side grows like -n, and the only root left is
        // at alpha, beta ~ 1/n, unreachable from any sane start
        a /= n;
        b /= n;

        double alpha0 = mean * ((mean * (1.0 - mean) / var) - 1.0);
        double beta0 = (1.0 - mean) * ((mean * (1.0 - mean) / var) - 1.0);
        if (!(alpha0 > 0.0) || !(beta0 > 0.0)) {
            // the moment estimate is positive only for var < mean(1 - mean);
            // start from the uniform rather than outside the parameter space
            alpha0 = 1.0;
            beta0 = 1.0;
        }

        return solveBetaScore(a, b, alpha0, beta0);
    }

    /**
     * Estimates the parameter {@code k} (degrees of freedom) of the ChiSquare
     * distribution from the observations {@code x} using the maximum likelihood
     * method. Note that this implementation allows for double-valued
     * estimators.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameter {@code k} (degrees of freedom)
     * @throws IllegalArgumentException
             if {@code x} is {@code null} or empty, or holds a value that is
             not finite and strictly positive
     */
    public static ParChiSquare getChiSquareMLE(double[] x) {
        int n = requirePositive(x);
        double sumLn = 0.0;
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
            sumLn += Math.log(x[i]);
        }

        // the likelihood is a function of the one parameter k, so the start is
        // one number: the moment estimate, since E[X] = k
        double mean = sum / n;
        double[] start = new double[] { mean > 0.0 ? mean : 1.0 };
        DMultiFunctionEval res = CGOptimizer.maximize(new ChiSquareMLE(sumLn, n), start);
        return new ParChiSquare(res.point[0]);
    }

    /**
     * Estimates the parameter &lambda; (rate) of the Exponential distribution
     * from the observations {@code x} using the maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the exponential rate parameter &lambda;
     * @throws IllegalArgumentException
             if {@code x} is {@code null} or empty, or holds a value that is
             not finite and strictly positive
     */
    public static ParExponential getExponentialMLE(double[] x) {
        int n = requirePositive(x);
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        return new ParExponential((double) n / sum);
    }

    /**
     * Estimates the parameters {@code mean} (&mu;) and {@code stdDev} (&sigma;)
     * of the Normal distribution from the observations {@code x} using the
     * maximum likelihood method.
     *
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters {@code mean} (&mu;) and {@code stdDev}
     *         (&sigma;)
     * @throws IllegalArgumentException
             if {@code x} is {@code null}, holds fewer than two observations,
             or holds a value that is not finite
     */
    public static ParNormal getNormalMLE(double[] x) {
        int n = requireFinite(x);
        requireAtLeast(n, 2);
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        double sigma = sum / n;

        sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            double dev = x[i] - sigma;
            sum = sum + (dev * dev);
        }
        return new ParNormal(sigma, Math.sqrt(sum / n));
    }

    /**
     * Solves the two Beta score equations
     * {@code psi(a) - psi(a + b) = meanLogX} and
     * {@code psi(b) - psi(a + b) = meanLog1mX} by Newton's method on the 2x2
     * system, halving the step whenever it would leave the positive quadrant.
     * <p>
     * Levenberg-Marquardt was measured here and stalls: a trial step that
     * overshoots into a non-positive parameter has no finite residual, so the
     * boundary acts as a wall, the trust region collapses and the search
     * returns its starting point -- which is the moment estimate, a plausible
     * number that is not a maximum-likelihood estimate. Over 5400 samples the
     * worst score residual was 6.4 for {@code lmder1} against 1.3e-14 here.
     */
    private static ParBeta solveBetaScore(double meanLogX, double meanLog1mX, double a, double b) {
        boolean converged = false;
        for (int it = 0; it < BETA_MAX_ITER; it++) {
            double psiSum = GammaFun.digamma(a + b);
            double f1 = GammaFun.digamma(a) - psiSum - meanLogX;
            double f2 = GammaFun.digamma(b) - psiSum - meanLog1mX;
            if (Math.abs(f1) <= BETA_TOL && Math.abs(f2) <= BETA_TOL) {
                // the equations are what has to hold, so they are what is
                // tested: the step length is not reachable at this precision
                converged = true;
                break;
            }
            double triSum = GammaFun.trigamma(a + b);
            double j11 = GammaFun.trigamma(a) - triSum;
            double j12 = -triSum;
            double j22 = GammaFun.trigamma(b) - triSum;
            double det = j11 * j22 - j12 * j12;
            double da = (j22 * f1 - j12 * f2) / det;
            double db = (j11 * f2 - j12 * f1) / det;
            int halvings = 0;
            // a NaN step fails this test too, so a singular Jacobian ends the
            // iteration rather than poisoning the estimate
            while (!(a - da > 0.0 && b - db > 0.0) && halvings < BETA_MAX_HALVINGS) {
                da *= 0.5;
                db *= 0.5;
                halvings++;
            }
            if (halvings >= BETA_MAX_HALVINGS) {
                break;
            }
            a -= da;
            b -= db;
        }
        return new ParBeta(a, b, converged);
    }

    /**
     * Steps outwards from {@code guess} until {@code f} changes sign, doubling
     * the half width every time, and returns the bracketing pair or
     * {@code null} if there is none within {@link #MAX_EXPANSIONS} steps.
     * <p>
     * {@link RootFinder#brentDekker} does not verify that its interval
     * brackets a root -- handed one that does not, it iterates to its
     * maximum and returns a number that looks like an answer. Every caller
     * here therefore finds the bracket first and treats a missing one as bad
     * input rather than as a result.
     */
    private static double[] bracketRoot(DFunction f, double guess) {
        double center = (guess > 0.0 && guess < Double.POSITIVE_INFINITY) ? guess : 1.0;
        double width = 0.5 * center;
        for (int i = 0; i < MAX_EXPANSIONS; i++) {
            double lo = Math.max(center - width, MIN_POSITIVE_SHAPE);
            double hi = center + width;
            double flo = f.apply(lo);
            double fhi = f.apply(hi);
            if ((flo <= 0.0 && fhi >= 0.0) || (flo >= 0.0 && fhi <= 0.0)) {
                return new double[] { lo, hi };
            }
            width *= 2.0;
        }
        return null;
    }

    /**
     * Rejects a sample that is absent, empty, or not made of numbers, and
     * returns its length.
     */
    private static int requireFinite(double[] x) {
        if (x == null) {
            throw new IllegalArgumentException("x must not be null");
        }
        if (x.length == 0) {
            throw new IllegalArgumentException(NO_OBS_MSG);
        }
        for (int i = 0; i < x.length; i++) {
            if (!isFinite(x[i])) {
                throw new IllegalArgumentException("x[" + i + "] is not a finite number : " + x[i]);
            }
        }
        return x.length;
    }

    /**
     * Rejects a sample that does not lie in {@code (0, infinity)}, which is
     * where the Gamma, Weibull, LogNormal, Exponential and ChiSquare
     * likelihoods are defined.
     */
    private static int requirePositive(double[] x) {
        int n = requireFinite(x);
        for (int i = 0; i < n; i++) {
            if (x[i] <= 0.0) {
                throw new IllegalArgumentException("x[" + i + "] must be strictly positive : " + x[i]);
            }
        }
        return n;
    }

    /** Rejects a sample that does not lie strictly inside {@code (0, 1)}. */
    private static int requireUnitInterval(double[] x) {
        int n = requireFinite(x);
        for (int i = 0; i < n; i++) {
            if (!(x[i] > 0.0 && x[i] < 1.0)) {
                throw new IllegalArgumentException("x[" + i + "] must lie in (0, 1) : " + x[i]);
            }
        }
        return n;
    }

    /**
     * Every estimator that has a scale or a shape to find needs a second
     * observation to see it in; the Cauchy fit needs a third.
     */
    private static void requireAtLeast(int n, int minimum) {
        if (n < minimum) {
            throw new IllegalArgumentException(
                    "at least " + minimum + " observations are needed, got " + n);
        }
    }

    /** {@code Double.isFinite} of Java 8, spelled out for the release 8 target. */
    private static boolean isFinite(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v);
    }

    private MLE() {
        throw new AssertionError();
    }
}
