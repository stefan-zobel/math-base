package math.stats.mle;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.cern.FastGamma;
import math.cern.GammaFun;
import math.distribution.Beta;
import math.distribution.Cauchy;
import math.distribution.ChiSquare;
import math.distribution.ContinuousDistribution;
import math.distribution.Exponential;
import math.distribution.Gamma;
import math.distribution.LogNormal;
import math.distribution.Normal;
import math.distribution.StudentT;
import math.distribution.Weibull;
import math.fun.DFunction;
import math.solve.RootFinder;

/**
 * The estimators are checked against the equations they claim to solve rather
 * than against remembered numbers: a maximum-likelihood estimate is the point
 * where the score vanishes, and that is a property the test can evaluate
 * itself. Where a closed form exists the test recomputes it.
 */
public final class MLETest {

    /** One call to an estimator, so a whole row of them can be swept. */
    private interface Estimator {
        void fit(double[] x);
    }

    private static void rejects(String name, Estimator estimator, double[] x, String what) {
        try {
            estimator.fit(x);
            fail(name + " accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(name + " threw without saying what was wrong", expected.getMessage() != null);
        }
    }

    /** Deterministic uniforms in {@code (0,1)}. */
    private static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) + 0.5) * 0x1.0p-53;
        }
    }

    private static double[] sample(ContinuousDistribution d, int n, long seed) {
        Lcg lcg = new Lcg(seed);
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = d.inverseCdf(lcg.next());
        }
        return x;
    }

    // ---------------------------------------------------------- the scores --

    /**
     * {@code |psi(a) - psi(a+b) - mean(ln x)| + |psi(b) - psi(a+b) - mean(ln(1-x))|},
     * which vanishes exactly at the Beta maximum-likelihood estimate.
     */
    private static double betaScore(double[] x, double alpha, double beta) {
        int n = x.length;
        double meanLogX = 0.0;
        double meanLog1mX = 0.0;
        for (int i = 0; i < n; i++) {
            meanLogX += Math.log(x[i]);
            meanLog1mX += Math.log1p(-x[i]);
        }
        double psiSum = GammaFun.digamma(alpha + beta);
        return Math.abs(GammaFun.digamma(alpha) - psiSum - meanLogX / n)
                + Math.abs(GammaFun.digamma(beta) - psiSum - meanLog1mX / n);
    }

    /** {@code |ln k - psi(k) - (ln mean(x) - mean(ln x))|}, zero at the Gamma shape. */
    private static double gammaScore(double[] x, double shape) {
        int n = x.length;
        double sum = 0.0;
        double sumLn = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i];
            sumLn += Math.log(x[i]);
        }
        return Math.abs(Math.log(shape) - GammaFun.digamma(shape) - (Math.log(sum / n) - sumLn / n));
    }

    /** {@code |sum x^k ln x / sum x^k - mean(ln x) - 1/k|}, zero at the Weibull shape. */
    private static double weibullScore(double[] x, double shape) {
        double sumXkLnX = 0.0;
        double sumXk = 0.0;
        double sumLn = 0.0;
        for (int i = 0; i < x.length; i++) {
            double lnx = Math.log(x[i]);
            double xk = Math.exp(shape * lnx);
            sumXkLnX += xk * lnx;
            sumXk += xk;
            sumLn += lnx;
        }
        return Math.abs(sumXkLnX / sumXk - sumLn / x.length - 1.0 / shape);
    }

    /** Gamma log-likelihood with the scale profiled out at {@code mean(x)/k}. */
    private static double gammaProfile(double[] x, double shape) {
        int n = x.length;
        double sum = 0.0;
        double sumLn = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i];
            sumLn += Math.log(x[i]);
        }
        double scale = sum / (n * shape);
        return (shape - 1.0) * sumLn - n * shape - n * FastGamma.logGamma(shape) - n * shape * Math.log(scale);
    }

    /** Weibull log-likelihood with the scale profiled out. */
    private static double weibullProfile(double[] x, double shape) {
        int n = x.length;
        double sumLn = 0.0;
        double sumXk = 0.0;
        for (int i = 0; i < n; i++) {
            sumLn += Math.log(x[i]);
            sumXk += Math.pow(x[i], shape);
        }
        return n * Math.log(shape) - n * Math.log(sumXk / n) + (shape - 1.0) * sumLn - n;
    }

    // ------------------------------------------------------------------ Beta --

    @Test
    public void testTheBetaEstimateSolvesTheScoreEquations() {
        // the score equations carry a factor 1/n. Against the raw sums the
        // right hand side grows like -n, the only root left sits at
        // alpha, beta ~ 1/n, and the search stalls on its starting point
        double[][] shapes = { { 2.0, 5.0 }, { 0.5, 0.5 }, { 5.0, 1.0 }, { 3.0, 3.0 }, { 0.2, 8.0 }, { 20.0, 20.0 } };
        double worst = 0.0;
        String at = "";
        for (int s = 0; s < shapes.length; s++) {
            for (long seed = 1; seed <= 40; seed++) {
                for (int n : new int[] { 10, 50, 500 }) {
                    double[] x = sample(new Beta(shapes[s][0], shapes[s][1]), n, seed * 7919L + 1);
                    ParBeta p = MLE.getBetaMLE(x);
                    assertTrue("Beta(" + shapes[s][0] + "," + shapes[s][1] + ") n=" + n + " seed=" + seed
                            + " did not converge", p.converged);
                    assertTrue("the estimate is not a Beta parameter pair : " + p.alpha + ", " + p.beta, p.isValid());
                    double score = betaScore(x, p.alpha, p.beta);
                    if (score > worst) {
                        worst = score;
                        at = "Beta(" + shapes[s][0] + "," + shapes[s][1] + ") n=" + n + " seed=" + seed;
                    }
                }
            }
        }
        // measured over 5400 samples: 1.9e-12, bounded by the solver's own
        // acceptance threshold
        assertTrue("worst Beta score residual " + worst + " at " + at, worst < 1.0e-10);
    }

    @Test
    public void testTheMomentEstimateDoesNotSolveTheBetaScoreEquations() {
        // pins what the defect was: the method used to return its start, and
        // the start is the moment estimate, which is a different estimator
        double[] x = sample(new Beta(2.0, 5.0), 1000, 53L);
        int n = x.length;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i];
        }
        double mean = sum / n;
        double centered = 0.0;
        for (int i = 0; i < n; i++) {
            centered += (x[i] - mean) * (x[i] - mean);
        }
        double var = centered / (n - 1);
        double momentAlpha = mean * ((mean * (1.0 - mean) / var) - 1.0);
        double momentBeta = (1.0 - mean) * ((mean * (1.0 - mean) / var) - 1.0);

        ParBeta p = MLE.getBetaMLE(x);
        assertTrue("the moment estimate already solves the score equations, so this sample proves nothing",
                betaScore(x, momentAlpha, momentBeta) > 1.0e-6);
        assertTrue("the fit came back as the moment estimate", p.alpha != momentAlpha || p.beta != momentBeta);
        assertTrue("the fit does not solve the score equations", betaScore(x, p.alpha, p.beta) < 1.0e-10);
    }

    // ----------------------------------------------------------------- Gamma --

    @Test
    public void testTheGammaShapeSolvesItsScoreEquation() {
        // the shape used to be sought in a fixed interval of half width 10
        // around the moment estimate. Brent does not check that an interval
        // brackets a root, so a shape outside it came back silently wrong
        double worst = 0.0;
        String at = "";
        for (double k : new double[] { 0.1, 0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0 }) {
            for (long seed = 1; seed <= 40; seed++) {
                for (int n : new int[] { 8, 30, 200 }) {
                    double[] x = sample(new Gamma(k, 1.0), n, seed * 7919L + 1);
                    ParGamma p = MLE.getGammaMLE(x);
                    assertTrue("k=" + k + " n=" + n + " seed=" + seed + " : " + p.shape + ", " + p.scale, p.isValid());
                    assertTrue("the shape must be positive : " + p.shape, p.shape > 0.0);
                    double score = gammaScore(x, p.shape);
                    if (score > worst) {
                        worst = score;
                        at = "k=" + k + " n=" + n + " seed=" + seed;
                    }
                }
            }
        }
        // measured worst over 7200 samples: 1.3e-12
        assertTrue("worst Gamma score residual " + worst + " at " + at, worst < 1.0e-9);
    }

    @Test
    public void testTheGammaShapeMaximizesTheProfileLikelihood() {
        // k = 50 at n = 8, seed 160 of this generator used to come back as
        // 99.197 where the maximum is at 90.151
        double worst = Double.MAX_VALUE;
        String at = "";
        for (double k : new double[] { 0.5, 3.0, 50.0, 1000.0 }) {
            for (long seed = 150; seed <= 190; seed++) {
                for (int n : new int[] { 8, 30 }) {
                    double[] x = sample(new Gamma(k, 1.0), n, seed * 7919L + 1);
                    double shape = MLE.getGammaMLE(x).shape;
                    double margin = gammaProfile(x, shape)
                            - Math.max(gammaProfile(x, shape * 0.999), gammaProfile(x, shape * 1.001));
                    if (margin < worst) {
                        worst = margin;
                        at = "k=" + k + " n=" + n + " seed=" + seed + " shape=" + shape;
                    }
                }
            }
        }
        assertTrue("a neighbour of the estimate has the higher likelihood, by " + (-worst) + ", at " + at, worst > 0.0);
    }

    // --------------------------------------------------------------- Weibull --

    @Test
    public void testTheWeibullShapeSolvesItsScoreEquation() {
        double worst = 0.0;
        String at = "";
        for (double k : new double[] { 0.1, 0.5, 2.0, 10.0, 60.0, 200.0 }) {
            for (long seed = 1; seed <= 40; seed++) {
                for (int n : new int[] { 8, 30, 200 }) {
                    double[] x = sample(new Weibull(1.0, k), n, seed * 104729L + 3);
                    ParWeibull p = MLE.getWeibullMLE(x);
                    assertTrue("k=" + k + " n=" + n + " seed=" + seed + " : " + p.shape + ", " + p.scale, p.isValid());
                    assertTrue("the shape must be positive : " + p.shape, p.shape > 0.0);
                    double score = weibullScore(x, p.shape);
                    if (score > worst) {
                        worst = score;
                        at = "k=" + k + " n=" + n + " seed=" + seed;
                    }
                }
            }
        }
        // measured worst over 5400 samples: 2.0e-12
        assertTrue("worst Weibull score residual " + worst + " at " + at, worst < 1.0e-9);
    }

    @Test
    public void testTheWeibullShapeMaximizesTheProfileLikelihood() {
        // k = 60 at n = 8, seed 117 of this generator used to come back as
        // 108.476 with a log-likelihood of 21.941, where the maximum is 77.408
        // with 23.069
        double worst = Double.MAX_VALUE;
        String at = "";
        for (double k : new double[] { 0.5, 2.0, 60.0, 200.0 }) {
            for (long seed = 100; seed <= 140; seed++) {
                for (int n : new int[] { 8, 30 }) {
                    double[] x = sample(new Weibull(1.0, k), n, seed * 104729L + 3);
                    double shape = MLE.getWeibullMLE(x).shape;
                    double margin = weibullProfile(x, shape)
                            - Math.max(weibullProfile(x, shape * 0.999), weibullProfile(x, shape * 1.001));
                    if (margin < worst) {
                        worst = margin;
                        at = "k=" + k + " n=" + n + " seed=" + seed + " shape=" + shape;
                    }
                }
            }
        }
        assertTrue("a neighbour of the estimate has the higher likelihood, by " + (-worst) + ", at " + at, worst > 0.0);
    }

    @Test
    public void testTheWeibullScaleFollowsFromTheShape() {
        // the scale is not searched for, it is the k-th root of the k-th
        // sample moment, so it has to reproduce exactly
        double[] x = sample(new Weibull(2.0, 1.7), 500, 41L);
        ParWeibull p = MLE.getWeibullMLE(x);
        double sumXk = 0.0;
        for (int i = 0; i < x.length; i++) {
            sumXk += Math.pow(x[i], p.shape);
        }
        assertEquals("scale", Math.pow(sumXk / x.length, 1.0 / p.shape), p.scale, 1.0e-12);
    }

    // ------------------------------------------------------------- ChiSquare --

    @Test
    public void testTheChiSquareEstimateMatchesTheExactRoot() {
        // the search used to start at 0.001 with a stray second coordinate the
        // likelihood does not depend on
        double worst = 0.0;
        String at = "";
        for (double k : new double[] { 0.5, 1.0, 7.0, 25.0, 100.0 }) {
            for (long seed = 1; seed <= 8; seed++) {
                final double[] x = sample(new ChiSquare(k), 500, seed * 7919L + 1);
                final int n = x.length;
                double s = 0.0;
                for (int i = 0; i < n; i++) {
                    s += Math.log(x[i]);
                }
                final double sumLn = s;
                // dl/dk = sumLn/2 - (n/2) psi(k/2) - (n/2) ln 2
                DFunction derivative = new DFunction() {
                    @Override
                    public double apply(double df) {
                        return sumLn / 2.0 - (n / 2.0) * GammaFun.digamma(df / 2.0) - n * Math.log(2.0) / 2.0;
                    }
                };
                double exact = RootFinder.brentDekker(1.0e-4, 5000.0, derivative, 1.0e-13);
                double got = MLE.getChiSquareMLE(x).degreesOfFreedom;
                double rel = Math.abs(got - exact) / exact;
                if (rel > worst) {
                    worst = rel;
                    at = "k=" + k + " seed=" + seed + " exact=" + exact + " got=" + got;
                }
            }
        }
        // measured worst: 1.2e-6 from the moment start, 1.2e-5 from the old one
        assertTrue("worst relative deviation from the exact root " + worst + " at " + at, worst < 5.0e-6);
    }

    // ------------------------------------------------------- the closed forms --

    @Test
    public void testTheNormalEstimateIsTheSampleMeanAndDeviation() {
        double[] x = sample(new Normal(3.0, 2.0), 500, 99L);
        ParNormal p = MLE.getNormalMLE(x);
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        double mean = sum / x.length;
        double centered = 0.0;
        for (int i = 0; i < x.length; i++) {
            centered += (x[i] - mean) * (x[i] - mean);
        }
        assertEquals("mean", mean, p.mean, 0.0);
        assertEquals("stdDev", Math.sqrt(centered / x.length), p.stdDev, 1.0e-15);
        assertTrue(p.isValid());
    }

    @Test
    public void testTheLogNormalEstimateIsTheNormalEstimateOfTheLogs() {
        double[] x = sample(new LogNormal(0.5, 1.2), 500, 7L);
        double[] logs = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            logs[i] = Math.log(x[i]);
        }
        ParLogNormal p = MLE.getLogNormalMLE(x);
        ParNormal onLogs = MLE.getNormalMLE(logs);
        assertEquals("mu", onLogs.mean, p.mu, 1.0e-14);
        assertEquals("sigma", onLogs.stdDev, p.sigma, 1.0e-14);
        assertTrue(p.isValid());
    }

    @Test
    public void testTheExponentialEstimateIsTheReciprocalMean() {
        double[] x = sample(new Exponential(2.5), 500, 11L);
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        ParExponential p = MLE.getExponentialMLE(x);
        assertEquals("lambda", x.length / sum, p.lambda, 0.0);
        assertTrue(p.isValid());
    }

    // -------------------------------------------------------------- StudentT --

    /** Log-likelihood of a t with the given parameters, and of its normal limit. */
    private static double tLogLikelihood(double[] x, double location, double scale, double df) {
        int n = x.length;
        double sum = 0.0;
        if (df == Double.POSITIVE_INFINITY) {
            for (int i = 0; i < n; i++) {
                double d = (x[i] - location) / scale;
                sum += d * d;
            }
            return -n * (Math.log(scale) + 0.5 * Math.log(2.0 * Math.PI)) - 0.5 * sum;
        }
        double constant = FastGamma.logGamma((df + 1.0) / 2.0) - FastGamma.logGamma(df / 2.0)
                - 0.5 * Math.log(Math.PI * df) - Math.log(scale);
        for (int i = 0; i < n; i++) {
            double d = (x[i] - location) / scale;
            sum += Math.log1p(d * d / df);
        }
        return n * constant - (df + 1.0) * 0.5 * sum;
    }

    @Test
    public void testTheDegreesOfFreedomAreAlwaysPositive() {
        // the start used to invert E[X^2] = nu/(nu-2). That map flattens as nu
        // grows, and once the sample second moment slips below one -- routine
        // above nu = 25 -- the start, and with it the answer, turns negative
        for (double nu : new double[] { 1.0, 2.5, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0 }) {
            for (long seed = 1; seed <= 12; seed++) {
                double[] x = sample(new StudentT(nu), 500, seed * 7919L + 1);
                ParStudentT p = MLE.getStudentTMLE(x);
                String at = "nu=" + nu + " seed=" + seed;
                assertTrue(at + " did not converge", p.converged);
                assertTrue(at + " gave df = " + p.df, p.df > 0.0);
                assertTrue(at + " gave scale = " + p.scale, p.scale > 0.0);
            }
        }
    }

    @Test
    public void testALightTailedSampleNoLongerComesBackNegative() {
        // StudentT(30) at n = 1000 from this seed has a sample second moment of
        // 0.9889, so the old moment start was 2m/(m-1) = -178.99, and the grid
        // walk never left it
        double[] x = sample(new StudentT(30.0), 1000, 89L);
        ParStudentT p = MLE.getStudentTMLE(x);
        assertTrue("df = " + p.df, p.df > 0.0);
        assertTrue("did not converge", p.converged);
    }

    @Test
    public void testTheFitIsEquivariantUnderLocationAndScale() {
        // fitting a*x + b has to give a*location + b, a*scale and the same df.
        // The estimator is exactly equivariant; only the stopping rule is not,
        // so the two runs may part company in the last digits
        double worstLocation = 0.0;
        double worstScale = 0.0;
        double worstDf = 0.0;
        for (double nu : new double[] { 2.0, 5.0 }) {
            for (double[] ab : new double[][] { { 8.0, 0.0 }, { 0.25, 0.0 }, { 8.0, 100.0 } }) {
                for (long seed = 1; seed <= 5; seed++) {
                    double[] x = sample(new StudentT(nu), 400, seed * 7919L + 1);
                    double[] y = new double[x.length];
                    for (int i = 0; i < x.length; i++) {
                        y[i] = ab[0] * x[i] + ab[1];
                    }
                    ParStudentT px = MLE.getStudentTMLE(x);
                    ParStudentT py = MLE.getStudentTMLE(y);
                    double wantLocation = ab[0] * px.location + ab[1];
                    double wantScale = ab[0] * px.scale;
                    worstLocation = Math.max(worstLocation,
                            Math.abs(py.location - wantLocation) / Math.max(1.0, Math.abs(wantLocation)));
                    worstScale = Math.max(worstScale, Math.abs(py.scale - wantScale) / wantScale);
                    assertTrue("df is finite on one side only : " + px.df + " and " + py.df,
                            (px.df == Double.POSITIVE_INFINITY) == (py.df == Double.POSITIVE_INFINITY));
                    if (px.df < Double.POSITIVE_INFINITY) {
                        worstDf = Math.max(worstDf, Math.abs(py.df - px.df) / px.df);
                    }
                }
            }
        }
        // measured: 1.4e-6, 2.0e-6 and 2.0e-4
        assertTrue("location is not equivariant, off by " + worstLocation, worstLocation < 1.0e-4);
        assertTrue("scale is not equivariant, off by " + worstScale, worstScale < 1.0e-4);
        assertTrue("df is not invariant, off by " + worstDf, worstDf < 1.0e-2);
    }

    @Test
    public void testTheStudentTEstimateMaximizesTheLikelihood() {
        double worst = Double.MAX_VALUE;
        String at = "";
        for (double nu : new double[] { 1.5, 3.0, 8.0 }) {
            for (long seed = 1; seed <= 8; seed++) {
                double[] x = sample(new StudentT(nu), 400, seed * 7919L + 1);
                ParStudentT p = MLE.getStudentTMLE(x);
                double best = tLogLikelihood(x, p.location, p.scale, p.df);
                double e = 1.0e-3;
                double[][] around = { { p.location * (1.0 + e) + e, p.scale, p.df },
                        { p.location * (1.0 - e) - e, p.scale, p.df },
                        { p.location, p.scale * (1.0 + e), p.df },
                        { p.location, p.scale * (1.0 - e), p.df },
                        { p.location, p.scale, p.df * (1.0 + e) },
                        { p.location, p.scale, p.df * (1.0 - e) } };
                for (int k = 0; k < around.length; k++) {
                    double margin = best - tLogLikelihood(x, around[k][0], around[k][1], around[k][2]);
                    if (margin < worst) {
                        worst = margin;
                        at = "nu=" + nu + " seed=" + seed + " neighbour " + k;
                    }
                }
            }
        }
        assertTrue("a neighbour has the higher likelihood, by " + (-worst) + ", at " + at, worst > 0.0);
    }

    @Test
    public void testASampleThatDoesNotIdentifyTheTailReportsInfinity() {
        // normal data cannot say where the tail of a t would be, so the
        // likelihood rises all the way and there is no finite maximizer. The
        // location and the scale are the normal estimates then, exactly: every
        // weight is one in the limit
        int found = 0;
        for (long seed = 1; seed <= 20 && found == 0; seed++) {
            double[] x = sample(new Normal(3.0, 2.0), 1000, seed * 7919L + 1);
            ParStudentT p = MLE.getStudentTMLE(x);
            if (p.df != Double.POSITIVE_INFINITY) {
                continue;
            }
            found++;
            assertTrue("an infinite df still has to converge", p.converged);
            assertTrue("infinite degrees of freedom are not a StudentT parameter", !p.isValid());
            double sum = 0.0;
            for (int i = 0; i < x.length; i++) {
                sum += x[i];
            }
            double mean = sum / x.length;
            double centered = 0.0;
            for (int i = 0; i < x.length; i++) {
                double d = x[i] - mean;
                centered += d * d;
            }
            assertEquals("location", mean, p.location, 0.0);
            assertEquals("scale", Math.sqrt(centered / x.length), p.scale, 0.0);
        }
        assertTrue("no normal sample in 20 was reported as unidentified", found > 0);
    }

    @Test
    public void testAHeavyTailIsRecovered() {
        // measured over 40 seeds: [2.42, 4.97] at n = 500 and [2.75, 3.33] at 5000
        for (long seed = 1; seed <= 10; seed++) {
            double[] x = sample(new StudentT(3.0), 2000, seed * 7919L + 1);
            ParStudentT p = MLE.getStudentTMLE(x);
            assertTrue("seed=" + seed + " gave df = " + p.df, p.df > 1.5 && p.df < 8.0);
            assertTrue("seed=" + seed + " gave scale = " + p.scale, p.scale > 0.5 && p.scale < 2.0);
            assertTrue("seed=" + seed + " gave location = " + p.location, Math.abs(p.location) < 0.5);
            assertTrue("seed=" + seed + " is not a valid fit", p.isValid());
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void testASingleObservationHasNoScale() {
        MLE.getStudentTMLE(new double[] { 1.0 });
    }

    @Test(expected = IllegalArgumentException.class)
    public void testAConstantSampleHasNoScale() {
        MLE.getStudentTMLE(new double[] { 2.0, 2.0, 2.0, 2.0 });
    }

    // ------------------------------------------------------------- the support --

    @Test
    public void testEveryEstimatorRejectsDataOutsideItsSupport() {
        double[] empty = new double[0];
        double[] withNaN = { 0.1, Double.NaN, 0.2 };
        double[] withInfinity = { 0.1, Double.POSITIVE_INFINITY, 0.2 };
        double[] one = { 0.5 };

        String[] positiveNames = { "getGammaMLE", "getWeibullMLE", "getLogNormalMLE", "getExponentialMLE",
                "getChiSquareMLE" };
        Estimator[] positive = { MLE::getGammaMLE, MLE::getWeibullMLE, MLE::getLogNormalMLE,
                MLE::getExponentialMLE, MLE::getChiSquareMLE };
        for (int i = 0; i < positive.length; i++) {
            rejects(positiveNames[i], positive[i], null, "null");
            rejects(positiveNames[i], positive[i], empty, "an empty sample");
            rejects(positiveNames[i], positive[i], withNaN, "a NaN");
            rejects(positiveNames[i], positive[i], withInfinity, "an infinity");
            rejects(positiveNames[i], positive[i], new double[] { 0.1, 0.0, 0.2 }, "a zero");
            rejects(positiveNames[i], positive[i], new double[] { 0.1, -3.0, 0.2 }, "a negative value");
        }

        String[] anywhereNames = { "getNormalMLE", "getStudentTMLE" };
        Estimator[] anywhere = { MLE::getNormalMLE, MLE::getStudentTMLE };
        for (int i = 0; i < anywhere.length; i++) {
            rejects(anywhereNames[i], anywhere[i], null, "null");
            rejects(anywhereNames[i], anywhere[i], empty, "an empty sample");
            rejects(anywhereNames[i], anywhere[i], withNaN, "a NaN");
            rejects(anywhereNames[i], anywhere[i], withInfinity, "an infinity");
            rejects(anywhereNames[i], anywhere[i], one, "a single observation");
        }

        rejects("getBetaMLE", MLE::getBetaMLE, null, "null");
        rejects("getBetaMLE", MLE::getBetaMLE, empty, "an empty sample");
        rejects("getBetaMLE", MLE::getBetaMLE, withNaN, "a NaN");
        rejects("getBetaMLE", MLE::getBetaMLE, new double[] { 0.2, 0.0, 0.4 }, "a zero");
        rejects("getBetaMLE", MLE::getBetaMLE, new double[] { 0.2, 1.0, 0.4 }, "a one");
        rejects("getBetaMLE", MLE::getBetaMLE, new double[] { 0.2, 1.5, 0.4 }, "a value above one");
        rejects("getBetaMLE", MLE::getBetaMLE, new double[] { 0.2, -0.5, 0.4 }, "a negative value");
        rejects("getBetaMLE", MLE::getBetaMLE, one, "a single observation");

        rejects("getCauchyMLE", MLE::getCauchyMLE, null, "null");
        rejects("getCauchyMLE", MLE::getCauchyMLE, empty, "an empty sample");
        rejects("getCauchyMLE", MLE::getCauchyMLE, withNaN, "a NaN");
        rejects("getCauchyMLE", MLE::getCauchyMLE, withInfinity, "an infinity");
        rejects("getCauchyMLE", MLE::getCauchyMLE, one, "a single observation");
        // two observations converge, but onto whichever optimum the start is
        // nearest: measured over 300 samples the starts spread by 117
        rejects("getCauchyMLE", MLE::getCauchyMLE, new double[] { 0.1, 0.2 }, "two observations");

        // the two estimators that need only one observation take one
        assertTrue(MLE.getExponentialMLE(one).isValid());
        assertTrue(MLE.getChiSquareMLE(one).isValid());
    }

    @Test
    public void testASingleOutOfSupportPointNoLongerCorruptsTheFit() {
        // ln(x) of a value the distribution cannot take used to be replaced by
        // ln(MIN_NORMAL / 2), about -709. One 0.0 among 200 Gamma observations
        // turned a shape of 3.0 into 0.196, and isValid() said true
        double[] x = sample(new Gamma(3.0, 2.0), 200, 5L);
        double clean = MLE.getGammaMLE(x).shape;
        assertTrue("the untouched sample already fits badly : " + clean, clean > 1.5 && clean < 6.0);
        x[0] = 0.0;
        try {
            MLE.getGammaMLE(x);
            fail("a zero observation was accepted, shape came back as " + MLE.getGammaMLE(x).shape);
        } catch (IllegalArgumentException expected) {
            assertTrue("the message does not name the offender : " + expected.getMessage(),
                    expected.getMessage().contains("x[0]"));
        }
    }

    // ---------------------------------------------------------------- Cauchy --

    /** Log-likelihood of a Cauchy, which is the t at one degree of freedom. */
    private static double cauchyLogLikelihood(double[] x, double location, double scale) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            double d = (x[i] - location) / scale;
            sum += Math.log1p(d * d);
        }
        return -x.length * Math.log(Math.PI * scale) - sum;
    }

    @Test
    public void testTheCauchyEstimateMaximizesTheLikelihood() {
        double worst = Double.MAX_VALUE;
        String at = "";
        for (int n : new int[] { 10, 100, 1000 }) {
            for (long seed = 1; seed <= 10; seed++) {
                double[] x = sample(new Cauchy(2.0, 3.0), n, seed * 7919L + 1);
                ParCauchy p = MLE.getCauchyMLE(x);
                assertTrue("n=" + n + " seed=" + seed + " is not a valid fit", p.isValid());
                double best = cauchyLogLikelihood(x, p.location, p.scale);
                double e = 1.0e-3;
                double[][] around = { { p.location * (1.0 + e) + e, p.scale },
                        { p.location * (1.0 - e) - e, p.scale }, { p.location, p.scale * (1.0 + e) },
                        { p.location, p.scale * (1.0 - e) } };
                for (int k = 0; k < around.length; k++) {
                    double margin = best - cauchyLogLikelihood(x, around[k][0], around[k][1]);
                    if (margin < worst) {
                        worst = margin;
                        at = "n=" + n + " seed=" + seed + " neighbour " + k;
                    }
                }
            }
        }
        assertTrue("a neighbour has the higher likelihood, by " + (-worst) + ", at " + at, worst > 0.0);
    }

    @Test
    public void testTheCauchyEstimateIsTheGlobalMaximum() {
        // the location likelihood of a Cauchy is the textbook multimodal one,
        // but only with the scale held fixed. The joint problem has a single
        // stationary point, so no point of a coarse grid may beat the fit --
        // which is also why the search does not need to be started carefully
        for (int n : new int[] { 5, 20, 60 }) {
            for (long seed = 1; seed <= 6; seed++) {
                double[] x = sample(new Cauchy(2.0, 3.0), n, seed * 7919L + 1);
                ParCauchy p = MLE.getCauchyMLE(x);
                double best = cauchyLogLikelihood(x, p.location, p.scale);
                double lo = x[0];
                double hi = x[0];
                for (int i = 1; i < n; i++) {
                    lo = Math.min(lo, x[i]);
                    hi = Math.max(hi, x[i]);
                }
                for (int a = 0; a <= 80; a++) {
                    double location = lo + (hi - lo) * a / 80.0;
                    for (int b = -30; b <= 30; b++) {
                        double scale = p.scale * Math.pow(10.0, b / 10.0);
                        double here = cauchyLogLikelihood(x, location, scale);
                        assertTrue("n=" + n + " seed=" + seed + " : the grid point (" + location + ", " + scale
                                + ") beats the fit (" + p.location + ", " + p.scale + ") by " + (here - best),
                                here <= best);
                    }
                }
            }
        }
    }

    @Test
    public void testTheCauchyFitIsEquivariantUnderLocationAndScale() {
        double worstLocation = 0.0;
        double worstScale = 0.0;
        for (double[] ab : new double[][] { { 8.0, 0.0 }, { 0.25, 0.0 }, { 8.0, 100.0 } }) {
            for (long seed = 1; seed <= 8; seed++) {
                double[] x = sample(new Cauchy(2.0, 3.0), 200, seed * 7919L + 1);
                double[] y = new double[x.length];
                for (int i = 0; i < x.length; i++) {
                    y[i] = ab[0] * x[i] + ab[1];
                }
                ParCauchy px = MLE.getCauchyMLE(x);
                ParCauchy py = MLE.getCauchyMLE(y);
                double wantLocation = ab[0] * px.location + ab[1];
                double wantScale = ab[0] * px.scale;
                worstLocation = Math.max(worstLocation,
                        Math.abs(py.location - wantLocation) / Math.max(1.0, Math.abs(wantLocation)));
                worstScale = Math.max(worstScale, Math.abs(py.scale - wantScale) / wantScale);
            }
        }
        // measured: 3.1e-5 and 3.7e-6
        assertTrue("location is not equivariant, off by " + worstLocation, worstLocation < 1.0e-3);
        assertTrue("scale is not equivariant, off by " + worstScale, worstScale < 1.0e-4);
    }

    @Test
    public void testTheCauchyFitBeatsTheFreeStudentTFitOnCauchyData() {
        // this is what makes the estimator worth having: the t fit spends a
        // parameter on degrees of freedom that are known to be 1, and the pair
        // it comes back with is therefore not the Cauchy maximum
        double worstGap = Double.MAX_VALUE;
        String at = "";
        for (int n : new int[] { 20, 100, 1000 }) {
            for (long seed = 1; seed <= 15; seed++) {
                double[] x = sample(new Cauchy(2.0, 3.0), n, seed * 7919L + 1);
                ParCauchy c = MLE.getCauchyMLE(x);
                ParStudentT t = MLE.getStudentTMLE(x);
                double gap = cauchyLogLikelihood(x, c.location, c.scale)
                        - cauchyLogLikelihood(x, t.location, t.scale);
                if (gap < worstGap) {
                    worstGap = gap;
                    at = "n=" + n + " seed=" + seed + " df=" + t.df;
                }
            }
        }
        assertTrue("the free t fit is the better Cauchy, by " + (-worstGap) + ", at " + at, worstGap >= 0.0);
    }

    @Test
    public void testTheCauchyEstimateRecoversTheParameters() {
        // measured over 60 seeds: location in [1.798, 2.214] and scale in
        // [2.780, 3.237] at n = 2000, none of them failing to converge
        for (long seed = 1; seed <= 12; seed++) {
            ParCauchy p = MLE.getCauchyMLE(sample(new Cauchy(2.0, 3.0), 2000, seed * 7919L + 1));
            assertTrue("seed=" + seed + " did not converge", p.converged);
            assertTrue("seed=" + seed + " gave location = " + p.location, Math.abs(p.location - 2.0) < 0.5);
            assertTrue("seed=" + seed + " gave scale = " + p.scale, Math.abs(p.scale - 3.0) < 0.5);
            assertTrue("seed=" + seed + " is not a valid fit", p.isValid());
        }
    }

    @Test
    public void testTheSmallestCauchySampleIsThree() {
        ParCauchy p = MLE.getCauchyMLE(sample(new Cauchy(2.0, 3.0), 3, 5L));
        assertTrue("three observations must be enough", p.isValid());
    }
}
