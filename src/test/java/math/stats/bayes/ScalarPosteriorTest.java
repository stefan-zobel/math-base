package math.stats.bayes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.cern.Arithmetic;
import math.cern.FastGamma;
import math.distribution.Beta;
import math.distribution.Binomial;
import math.distribution.ContinuousDistribution;
import math.distribution.Gamma;
import math.distribution.LogNormal;
import math.distribution.Normal;
import math.distribution.Poisson;
import math.distribution.StudentT;
import math.fun.DFunction;

/**
 * {@link ScalarPosterior}, which normalizes a log posterior by quadrature.
 * <p>
 * The whole of it is checked the same way: hand it the log density of a law the
 * library already knows and demand the law back. Every member has an
 * independent closed form to be wrong against, so nothing here rests on a
 * golden value. The conjugate cases at the end go further and hand it a
 * posterior it has to <em>derive</em> -- prior times likelihood -- whose answer
 * is a distribution the library can also write down.
 * <p>
 * Tolerances were measured before they were written. Over the six laws below
 * the worst figures seen were {@code 3.4e-16} on the log evidence,
 * {@code 2.3e-16} relative on the density, {@code 3.8e-15} on the distribution
 * function and {@code 6.7e-15} relative on the quantile; what is asserted is
 * roughly two orders of magnitude above each.
 */
public final class ScalarPosteriorTest {

    private static final double INF = Double.POSITIVE_INFINITY;

    private static final Normal NORMAL = new Normal(3.0, 2.0);
    private static final Beta BETA = new Beta(2.0, 5.0);
    private static final Gamma GAMMA = new Gamma(3.0, 2.0);
    private static final LogNormal LOG_NORMAL = new LogNormal(0.0, 1.0);
    private static final StudentT T5 = new StudentT(5.0);

    /** The five laws, with the range each of them lives on. */
    private static final ContinuousDistribution[] LAWS = { NORMAL, BETA, GAMMA, LOG_NORMAL, T5 };
    private static final double[] LOWER = { -INF, 0.0, 0.0, 0.0, -INF };
    private static final double[] UPPER = { INF, 1.0, INF, INF, INF };
    private static final String[] NAMES = { "Normal(3, 2)", "Beta(2, 5)", "Gamma(3, 2)", "LogNormal(0, 1)",
            "StudentT(5)" };

    /** Built once: a construction is a few dozen integrations. */
    private static final ScalarPosterior[] POSTERIORS = new ScalarPosterior[LAWS.length];

    static {
        for (int i = 0; i < LAWS.length; i++) {
            POSTERIORS[i] = posteriorOf(LAWS[i], LOWER[i], UPPER[i]);
        }
    }

    private static ScalarPosterior posteriorOf(final ContinuousDistribution law, double lo, double hi) {
        return ScalarPosterior.of(new DFunction() {
            @Override
            public double apply(double x) {
                return law.logPdf(x);
            }
        }, lo, hi);
    }

    // ------------------------------------------------------------------
    // the evidence and the density
    // ------------------------------------------------------------------

    @Test
    public void anAlreadyNormalizedDensityHasNoEvidenceToFind() {
        // a density that already integrates to one leaves nothing for the
        // normalizer to do, so the log evidence is the one number that is known
        // in advance for every law at once. Worst figure measured: 3.4e-16
        for (int i = 0; i < LAWS.length; i++) {
            assertEquals(NAMES[i] + " is already normalized", 0.0, POSTERIORS[i].logEvidence(), 1.0e-13);
        }
    }

    @Test
    public void theEvidenceRecoversAConstantThatWasAddedToIt() {
        // what makes the evidence usable for model comparison: an unnormalized
        // posterior differs from a normalized one by a constant factor, and
        // that factor has to come back exactly as it went in -- including one
        // large enough that the density it produces is not representable
        for (final double offset : new double[] { 12.345, -12.345, 700.0, -5000.0 }) {
            final Normal law = NORMAL;
            ScalarPosterior p = ScalarPosterior.of(new DFunction() {
                @Override
                public double apply(double x) {
                    return law.logPdf(x) + offset;
                }
            }, -INF, INF);
            assertEquals("an offset of " + offset + " was not recovered", offset, p.logEvidence(),
                    1.0e-12 * Math.max(1.0, Math.abs(offset)));
        }
    }

    @Test
    public void theDensityIsTheDensityItWasGiven() {
        for (int i = 0; i < LAWS.length; i++) {
            ContinuousDistribution law = LAWS[i];
            ScalarPosterior p = POSTERIORS[i];
            for (int q = 1; q < 100; q++) {
                double x = law.inverseCdf(q / 100.0);
                if (!(x > LOWER[i]) || !(x < UPPER[i])) {
                    continue;
                }
                assertEquals(NAMES[i] + ": the density differs at " + x, law.pdf(x), p.pdf(x),
                        1.0e-13 * Math.max(1.0, law.pdf(x)));
                assertEquals(NAMES[i] + ": the log density differs at " + x, law.logPdf(x), p.logPdf(x), 1.0e-12);
            }
        }
    }

    @Test
    public void outsideTheRangeThereIsNoDensity() {
        ScalarPosterior p = POSTERIORS[1];
        for (double x : new double[] { -1.0, -1.0e-9, 1.0 + 1.0e-9, 2.0, Double.NaN }) {
            assertEquals("density at " + x, 0.0, p.pdf(x), 0.0);
            assertEquals("log density at " + x, Double.NEGATIVE_INFINITY, p.logPdf(x), 0.0);
        }
    }

    @Test
    public void theWholeThingWouldBeZeroWithoutTheShift() {
        // the test the design exists for. Five thousand observations put the
        // unnormalized density below the smallest double at every point a rule
        // will evaluate, so without the mode shift the normalizer is zero and
        // every member here is 0/0
        double[] data = normalSample(5000, 7.0, 1.5);
        DFunction logPosterior = normalLogLikelihood(data, 1.5);

        assertEquals("the unshifted density has not underflowed, so this test proves nothing", 0.0,
                Math.exp(logPosterior.apply(7.0)), 0.0);

        ScalarPosterior p = ScalarPosterior.of(logPosterior, -INF, INF, 0.0, 1.0e-10);
        assertTrue("the log evidence is not finite : " + p.logEvidence(), isFinite(p.logEvidence()));
        assertTrue("the density at the mode is not positive : " + p.pdf(p.mode()), p.pdf(p.mode()) > 0.0);

        // under a flat prior the posterior of the mean is exactly normal, which
        // is what all of it is checked against
        Normal exact = new Normal(average(data), 1.5 / Math.sqrt(data.length));
        assertEquals("the mode is not the sample average", exact.mean(), p.mode(), 1.0e-6);
        assertEquals("the mean is wrong", exact.mean(), p.mean(), 1.0e-10);
        assertEquals("the variance is wrong", exact.variance(), p.variance(), 1.0e-10 * exact.variance());
        assertEquals("the density is wrong at the peak", exact.pdf(exact.mean()), p.pdf(exact.mean()),
                1.0e-10 * exact.pdf(exact.mean()));
        assertTrue("the two routes to the normalizer disagree : " + p.errorEstimate(), p.errorEstimate() < 1.0e-11);
    }

    @Test
    public void theModeIsWhereTheMassIs() {
        // BrentMinimizer resolves an abscissa to about the square root of the
        // machine epsilon, and says so: a function is quadratically flat at its
        // extremum. So this is a 1e-6 statement, not a 1e-15 one
        assertEquals("Normal", 3.0, POSTERIORS[0].mode(), 1.0e-6);
        assertEquals("Beta(2, 5) peaks at (a-1)/(a+b-2)", 0.2, POSTERIORS[1].mode(), 1.0e-6);
        assertEquals("Gamma(3, 2) peaks at (k-1) theta", 4.0, POSTERIORS[2].mode(), 1.0e-6);
        assertEquals("LogNormal(0, 1) peaks at exp(mu - sigma^2)", Math.exp(-1.0), POSTERIORS[3].mode(), 1.0e-6);
        assertEquals("StudentT", 0.0, POSTERIORS[4].mode(), 1.0e-6);
    }

    @Test
    public void theScalesDifferWhereTheLawIsSkewed() {
        // half a nat below the peak is exactly one standard deviation for a
        // normal, which is the calibration; for a skewed law the two sides
        // differ, and building the panel ladder out of one number would put
        // half of it in the wrong place
        assertEquals("a normal falls half a nat at one standard deviation", 2.0, POSTERIORS[0].leftScale(), 1.0e-9);
        assertEquals("a normal falls half a nat at one standard deviation", 2.0, POSTERIORS[0].rightScale(), 1.0e-9);
        // the two half-widths are found by separate root searches at the square
        // root of the machine epsilon, so symmetry is a 1e-6 statement and not
        // a 1e-15 one. Measured gap: 1.9e-8
        assertEquals("a Student t is symmetric", POSTERIORS[4].leftScale(), POSTERIORS[4].rightScale(), 1.0e-6);

        for (int i : new int[] { 1, 2, 3 }) {
            assertTrue(NAMES[i] + " is skewed right, so its right half-width must be the larger",
                    POSTERIORS[i].rightScale() > POSTERIORS[i].leftScale());
        }
    }

    @Test
    public void theTwoRoutesToTheNormalizerAgree() {
        // the panel ladder and one adaptive integration over the whole range
        // subdivide the same integral differently, so this owes nothing to the
        // tolerance either of them was asked for. Worst figure measured: 5.1e-16
        for (int i = 0; i < LAWS.length; i++) {
            double gap = POSTERIORS[i].errorEstimate();
            assertTrue(NAMES[i] + ": the two routes disagree by " + gap, gap < 1.0e-12);
        }
    }

    @Test
    public void aPosteriorFarFromTheOriginIsStillFound() {
        // InfiniteIntegrator's own javadoc names this as the case its probe
        // ladder cannot see -- a unit-width feature at 1e4 underflows to zero
        // at every rung -- and the cure it names is a center, which this class
        // always has because it locates the mode first
        Normal far = new Normal(1.0e4, 1.0);
        ScalarPosterior p = posteriorOf(far, -INF, INF);
        assertEquals("the evidence is wrong", 0.0, p.logEvidence(), 1.0e-12);
        assertEquals("the mode is wrong", 1.0e4, p.mode(), 1.0e-4);
        assertEquals("the mean is wrong", 1.0e4, p.mean(), 1.0e-8);
        assertEquals("the variance is wrong", 1.0, p.variance(), 1.0e-9);
        assertTrue("the two routes disagree : " + p.errorEstimate(), p.errorEstimate() < 1.0e-11);
        for (double q : new double[] { 0.01, 0.25, 0.5, 0.75, 0.99 }) {
            assertEquals("the quantile is wrong at " + q, far.inverseCdf(q), p.inverseCdf(q), 1.0e-9);
        }
    }

    @Test
    public void aPeakAtTheEndOfTheRangeIsStillAPeak() {
        // an exponential posterior peaks exactly on the boundary, where there
        // is no room to the left at all and the width has to come from the
        // other side
        ScalarPosterior p = ScalarPosterior.of(new DFunction() {
            @Override
            public double apply(double x) {
                return -x;
            }
        }, 0.0, INF);
        assertEquals("the peak is not at the end", 0.0, p.mode(), 1.0e-6);
        assertEquals("the evidence of exp(-x) over (0, inf) is log(1)", 0.0, p.logEvidence(), 1.0e-12);
        assertEquals("the mean of a unit exponential", 1.0, p.mean(), 1.0e-10);
        assertEquals("the variance of a unit exponential", 1.0, p.variance(), 1.0e-9);
        assertEquals("the median of a unit exponential", Math.log(2.0), p.inverseCdf(0.5), 1.0e-10);
    }

    @Test
    public void thePosteriorRejectsWhatItCannotIntegrate() {
        final Normal law = NORMAL;
        DFunction ok = new DFunction() {
            @Override
            public double apply(double x) {
                return law.logPdf(x);
            }
        };
        refuses("a null integrand", null, 0.0, 1.0, 0.5, 1.0e-10);
        refuses("an empty range", ok, 1.0, 1.0, 1.0, 1.0e-10);
        refuses("a reversed range", ok, 1.0, 0.0, 0.5, 1.0e-10);
        refuses("a NaN lower end", ok, Double.NaN, 1.0, 0.5, 1.0e-10);
        refuses("a NaN upper end", ok, 0.0, Double.NaN, 0.5, 1.0e-10);
        refuses("a start outside the range", ok, 0.0, 1.0, 2.0, 1.0e-10);
        refuses("an infinite start", ok, -INF, INF, INF, 1.0e-10);
        refuses("a zero tolerance", ok, -INF, INF, 0.0, 0.0);
        refuses("a negative tolerance", ok, -INF, INF, 0.0, -1.0e-10);

        // and the one that is about the posterior rather than the arguments: a
        // flat log posterior over the whole line cannot be normalized, and the
        // width search is what notices
        refuses("an improper posterior", new DFunction() {
            @Override
            public double apply(double x) {
                return 0.0;
            }
        }, -INF, INF, 0.0, 1.0e-10);
    }

    // ------------------------------------------------------------------
    // the distribution function
    // ------------------------------------------------------------------

    @Test
    public void theCdfIsMonotoneAcrossEveryPanel() {
        // the panel sums are formed once and the partial panel is clamped into
        // the gap it belongs to, so this is exact rather than true to within a
        // quadrature tolerance -- which is what makes the root search in
        // inverseCdf a well posed question out in the tails
        for (int i = 0; i < LAWS.length; i++) {
            ScalarPosterior p = POSTERIORS[i];
            double from = isFinite(LOWER[i]) ? LOWER[i] : p.mode() - 30.0 * p.leftScale();
            double to = isFinite(UPPER[i]) ? UPPER[i] : p.mode() + 30.0 * p.rightScale();
            double previous = 0.0;
            for (int k = 0; k <= 2000; k++) {
                double x = from + (to - from) * k / 2000.0;
                double c = p.cdf(x);
                assertTrue(NAMES[i] + ": the cdf fell at " + x + ", " + c + " after " + previous, c >= previous);
                assertTrue(NAMES[i] + ": the cdf left [0, 1] at " + x + " : " + c, c >= 0.0 && c <= 1.0);
                previous = c;
            }
            assertEquals(NAMES[i] + ": the cdf does not start at zero", 0.0, p.cdf(LOWER[i]), 0.0);
            assertEquals(NAMES[i] + ": the cdf does not close at one", 1.0, p.cdf(UPPER[i]), 0.0);
        }
    }

    @Test
    public void theCdfIsTheCdfItWasGiven() {
        for (int i = 0; i < LAWS.length; i++) {
            ContinuousDistribution law = LAWS[i];
            ScalarPosterior p = POSTERIORS[i];
            for (int q = 1; q < 100; q++) {
                double x = law.inverseCdf(q / 100.0);
                if (!(x > LOWER[i]) || !(x < UPPER[i])) {
                    continue;
                }
                assertEquals(NAMES[i] + ": the cdf differs at " + x, law.cdf(x), p.cdf(x), 1.0e-12);
            }
        }
    }

    @Test
    public void theQuantileInvertsTheCdf() {
        for (int i = 0; i < LAWS.length; i++) {
            ScalarPosterior p = POSTERIORS[i];
            for (int q = 1; q < 200; q++) {
                double target = q / 200.0;
                double x = p.inverseCdf(target);
                assertEquals(NAMES[i] + ": cdf(inverseCdf(" + target + ")) missed", target, p.cdf(x), 1.0e-11);
            }
            assertEquals(NAMES[i] + ": inverseCdf(0)", LOWER[i], p.inverseCdf(0.0), 0.0);
            assertEquals(NAMES[i] + ": inverseCdf(1)", UPPER[i], p.inverseCdf(1.0), 0.0);
        }
    }

    @Test
    public void theQuantileIsTheQuantileItWasGiven() {
        for (int i = 0; i < LAWS.length; i++) {
            ContinuousDistribution law = LAWS[i];
            ScalarPosterior p = POSTERIORS[i];
            for (int q = 1; q < 100; q++) {
                double target = q / 100.0;
                double expected = law.inverseCdf(target);
                assertEquals(NAMES[i] + ": the quantile differs at " + target, expected, p.inverseCdf(target),
                        1.0e-12 * Math.max(1.0, Math.abs(expected)));
            }
        }
    }

    @Test
    public void theInheritedDefaultsWork() {
        // probability(x0, x1) and supportLowerBound come from
        // ContinuousDistribution, which is the whole point of implementing it
        for (int i = 0; i < LAWS.length; i++) {
            ScalarPosterior p = POSTERIORS[i];
            assertEquals(NAMES[i] + ": the support starts elsewhere", LOWER[i], p.supportLowerBound(), 0.0);
            assertEquals(NAMES[i] + ": the support ends elsewhere", UPPER[i], p.supportUpperBound(), 0.0);
            double a = p.inverseCdf(0.1);
            double b = p.inverseCdf(0.9);
            assertEquals(NAMES[i] + ": probability(a, b) is not the cdf difference", 0.8, p.probability(a, b),
                    1.0e-11);
        }
    }

    @Test
    public void theQuantileRejectsWhatItCannotAnswer() {
        ScalarPosterior p = POSTERIORS[0];
        for (double q : new double[] { -1.0e-15, -1.0, 1.0 + 1.0e-15, 2.0, Double.NaN, INF, -INF }) {
            try {
                p.inverseCdf(q);
                fail("the quantile accepted p = " + q);
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------
    // the moments
    // ------------------------------------------------------------------

    @Test
    public void theMomentsAreTheMomentsItWasGiven() {
        for (int i = 0; i < LAWS.length; i++) {
            ContinuousDistribution law = LAWS[i];
            ScalarPosterior p = POSTERIORS[i];
            assertEquals(NAMES[i] + ": the mean differs", law.mean(), p.mean(),
                    1.0e-10 * Math.max(1.0, Math.abs(law.mean())));
            assertEquals(NAMES[i] + ": the variance differs", law.variance(), p.variance(),
                    1.0e-10 * Math.max(1.0, law.variance()));
        }
    }

    @Test
    public void theVarianceOfAHeavyTailIsATruncation() {
        // the one place where a number comes back that should not, and it is
        // documented rather than detected. A Student t with two degrees of
        // freedom has no variance, and what is reported over [-R, R] grows by
        // exactly log(100) for every factor of a hundred in R -- the signature
        // of a logarithmic divergence being cut off rather than integrated
        DFunction t2 = new DFunction() {
            private final StudentT law = new StudentT(2.0);

            @Override
            public double apply(double x) {
                return law.logPdf(x);
            }
        };
        double[] variances = new double[3];
        double[] ranges = { 1.0e2, 1.0e4, 1.0e6 };
        for (int i = 0; i < ranges.length; i++) {
            variances[i] = ScalarPosterior.of(t2, -ranges[i], ranges[i], 0.0, 1.0e-10).variance();
        }
        // log(100) per decade, so twice that for each factor of a hundred here.
        // Measured: 9.2093 against the 9.2103 the divergence predicts
        double firstStep = variances[1] - variances[0];
        double secondStep = variances[2] - variances[1];
        double perHundred = 2.0 * Math.log(100.0);
        assertEquals("the growth per factor of a hundred is not 2 log(100)", perHundred, firstStep, 0.01);
        assertEquals("the growth per factor of a hundred is not 2 log(100)", perHundred, secondStep, 0.01);

        // where the variance does exist the range makes no difference
        DFunction t5 = new DFunction() {
            @Override
            public double apply(double x) {
                return T5.logPdf(x);
            }
        };
        // where the variance does exist, widening the range settles it instead
        // of growing it: the tail beyond a hundred is worth 3.2e-5 and the tail
        // beyond ten thousand nothing at all. The two behaviours are five
        // orders of magnitude apart, which is what makes them distinguishable
        double narrow = ScalarPosterior.of(t5, -1.0e2, 1.0e2, 0.0, 1.0e-10).variance();
        double wide = ScalarPosterior.of(t5, -1.0e4, 1.0e4, 0.0, 1.0e-10).variance();
        assertEquals("StudentT(5) has a variance and it is 5/3", 5.0 / 3.0, wide, 1.0e-9);
        double settled = Math.abs(wide - narrow);
        assertTrue("a variance that exists moved by " + settled + " over the same widening", settled < 1.0e-4);
        assertTrue("the two behaviours are not far enough apart to tell apart", firstStep > 1000.0 * settled);
    }

    @Test
    public void theMeanIsNotTheMode() {
        // a skewed posterior pulls its mean towards the long tail, and all
        // three of these are skewed to the right
        for (int i : new int[] { 1, 2, 3 }) {
            assertTrue(NAMES[i] + ": the mean should sit above the mode",
                    POSTERIORS[i].mean() > POSTERIORS[i].mode());
        }
        assertEquals("a symmetric posterior has them in the same place", POSTERIORS[0].mean(),
                POSTERIORS[0].mode(), 1.0e-6);
    }

    // ------------------------------------------------------------------
    // credible intervals
    // ------------------------------------------------------------------

    @Test
    public void theIntervalHoldsTheMassItPromises() {
        for (double level : new double[] { 0.5, 0.8, 0.95, 0.99 }) {
            for (int i = 0; i < LAWS.length; i++) {
                ScalarPosterior p = POSTERIORS[i];
                ScalarPosterior.Interval equal = p.credibleInterval(level);
                assertEquals(NAMES[i] + ": the equal-tailed interval at " + level + " holds the wrong mass", level,
                        equal.mass, 1.0e-9);
                assertEquals(NAMES[i] + ": mass and cdf disagree", p.cdf(equal.upper) - p.cdf(equal.lower),
                        equal.mass, 0.0);
                assertTrue(NAMES[i] + ": the interval is not an interval", equal.upper > equal.lower);
                assertEquals(NAMES[i] + ": width", equal.upper - equal.lower, equal.width(), 0.0);

                ScalarPosterior.Interval dense = p.highestDensityInterval(level);
                assertEquals(NAMES[i] + ": the density interval at " + level + " holds the wrong mass", level,
                        dense.mass, 1.0e-9);
                assertTrue(NAMES[i] + ": the density interval is not an interval", dense.upper > dense.lower);
            }
        }
    }

    @Test
    public void theEqualTailedIntervalIsTheQuantilePair() {
        ScalarPosterior p = POSTERIORS[3];
        ScalarPosterior.Interval interval = p.credibleInterval(0.95);
        // not bit for bit: the interval forms its tail as (1 - level) / 2,
        // which is 0.025000000000000022 and not 0.025, and a root search
        // answers the question it was actually asked
        assertEquals("the lower end is not the 2.5% quantile", p.inverseCdf(0.025), interval.lower, 1.0e-12);
        assertEquals("the upper end is not the 97.5% quantile", p.inverseCdf(0.975), interval.upper, 1.0e-12);
    }

    @Test
    public void theDensityIntervalIsNarrowerWhereTheLawIsSkewed() {
        // the point of the thing: for a skewed posterior the equal-tailed
        // interval is not the shortest one holding that mass, and it puts its
        // lower end where the density is higher than at points it excludes
        for (int i : new int[] { 1, 2, 3 }) {
            ScalarPosterior p = POSTERIORS[i];
            ScalarPosterior.Interval equal = p.credibleInterval(0.95);
            ScalarPosterior.Interval dense = p.highestDensityInterval(0.95);
            assertTrue(NAMES[i] + ": the density interval is not shorter : " + dense.width() + " vs "
                    + equal.width(), dense.width() < equal.width());
        }
        // and for a symmetric one the two coincide
        ScalarPosterior symmetric = POSTERIORS[0];
        ScalarPosterior.Interval equal = symmetric.credibleInterval(0.95);
        ScalarPosterior.Interval dense = symmetric.highestDensityInterval(0.95);
        assertEquals("a symmetric posterior should give the same two intervals", equal.lower, dense.lower, 1.0e-6);
        assertEquals("a symmetric posterior should give the same two intervals", equal.upper, dense.upper, 1.0e-6);
    }

    @Test
    public void theDensityIntervalHasEqualDensityAtItsEnds() {
        // which is what defines it, and is the condition it is found through
        for (int i : new int[] { 1, 2, 3 }) {
            ScalarPosterior p = POSTERIORS[i];
            ScalarPosterior.Interval dense = p.highestDensityInterval(0.9);
            double atLower = p.pdf(dense.lower);
            double atUpper = p.pdf(dense.upper);
            assertEquals(NAMES[i] + ": the density is not level at the ends", atLower, atUpper,
                    1.0e-5 * Math.max(atLower, atUpper));
        }
    }

    @Test
    public void bothKindsOfIntervalRejectALevelOutsideZeroAndOne() {
        ScalarPosterior p = POSTERIORS[0];
        for (double level : new double[] { 0.0, 1.0, -0.5, 1.5, Double.NaN, INF }) {
            try {
                p.credibleInterval(level);
                fail("the equal-tailed interval accepted a level of " + level);
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
            try {
                p.highestDensityInterval(level);
                fail("the density interval accepted a level of " + level);
            } catch (IllegalArgumentException expected) {
                assertTrue("threw without a message", expected.getMessage() != null);
            }
        }
    }

    // ------------------------------------------------------------------
    // the conjugate families, where the answer is arithmetic
    // ------------------------------------------------------------------

    @Test
    public void theBetaBinomialPosteriorIsTheOneTheArithmeticGives() {
        // here the posterior is derived rather than handed over: prior times
        // likelihood, integrated. That it comes back as the Beta the conjugate
        // update names is a statement about the whole apparatus at once
        final double a = 2.0;
        final double b = 3.0;
        final int trials = 20;
        final int successes = 7;
        final Beta prior = new Beta(a, b);
        ScalarPosterior p = ScalarPosterior.of(new DFunction() {
            @Override
            public double apply(double x) {
                return prior.logPdf(x) + new Binomial(trials, x).logPmf(successes);
            }
        }, 0.0, 1.0);

        Beta exact = new Beta(a + successes, b + trials - successes);
        for (int q = 1; q < 100; q++) {
            double x = q / 100.0;
            assertEquals("the density differs at " + x, exact.pdf(x), p.pdf(x), 1.0e-11);
            assertEquals("the cdf differs at " + x, exact.cdf(x), p.cdf(x), 1.0e-11);
            assertEquals("the quantile differs at " + x, exact.inverseCdf(x), p.inverseCdf(x), 1.0e-10);
        }
        assertEquals("the mean differs", exact.mean(), p.mean(), 1.0e-11);
        assertEquals("the variance differs", exact.variance(), p.variance(), 1.0e-12);
    }

    @Test
    public void theBetaBinomialEvidenceIsTheClosedFormOne() {
        // the marginal likelihood of the data, which is what a Bayes factor is
        // a ratio of, against the beta-binomial mass in closed form
        final double a = 2.0;
        final double b = 3.0;
        final Beta prior = new Beta(a, b);
        for (final int trials : new int[] { 5, 20, 60 }) {
            for (final int successes : new int[] { 0, 1, trials / 3, trials }) {
                ScalarPosterior p = ScalarPosterior.of(new DFunction() {
                    @Override
                    public double apply(double x) {
                        return prior.logPdf(x) + new Binomial(trials, x).logPmf(successes);
                    }
                }, 0.0, 1.0);
                double logChoose = Arithmetic.logFactorial(trials) - Arithmetic.logFactorial(successes)
                        - Arithmetic.logFactorial(trials - successes);
                double expected = logChoose + logBeta(a + successes, b + trials - successes) - logBeta(a, b);
                assertEquals("the evidence at " + successes + " of " + trials + " is wrong", expected,
                        p.logEvidence(), 1.0e-11 * Math.max(1.0, Math.abs(expected)));
            }
        }
    }

    @Test
    public void theGammaPoissonPosteriorIsTheOneTheArithmeticGives() {
        final double shape = 2.5;
        final double scale = 1.5;
        final int[] counts = { 3, 7, 2, 5, 4, 9, 1, 6 };
        final Gamma prior = new Gamma(shape, scale);
        ScalarPosterior p = ScalarPosterior.of(new DFunction() {
            @Override
            public double apply(double lambda) {
                double value = prior.logPdf(lambda);
                Poisson likelihood = new Poisson(lambda);
                for (int i = 0; i < counts.length; i++) {
                    value += likelihood.logPmf(counts[i]);
                }
                return value;
            }
        }, 0.0, INF);

        int total = 0;
        for (int i = 0; i < counts.length; i++) {
            total += counts[i];
        }
        Gamma exact = new Gamma(shape + total, scale / (1.0 + counts.length * scale));
        assertEquals("the mean differs", exact.mean(), p.mean(), 1.0e-10);
        assertEquals("the variance differs", exact.variance(), p.variance(), 1.0e-10);
        for (double q : new double[] { 0.025, 0.25, 0.5, 0.75, 0.975 }) {
            assertEquals("the quantile differs at " + q, exact.inverseCdf(q), p.inverseCdf(q), 1.0e-9);
        }
    }

    @Test
    public void theNormalNormalPosteriorIsTheOneTheArithmeticGives() {
        final double priorMean = 0.0;
        final double priorSd = 3.0;
        final double sigma = 1.5;
        final double[] data = normalSample(40, 7.0, sigma);
        final Normal prior = new Normal(priorMean, priorSd);
        ScalarPosterior p = ScalarPosterior.of(new DFunction() {
            @Override
            public double apply(double mu) {
                return prior.logPdf(mu) + normalLogLikelihood(data, sigma).apply(mu);
            }
        }, -INF, INF);

        double precision = 1.0 / (priorSd * priorSd) + data.length / (sigma * sigma);
        double posteriorMean = (priorMean / (priorSd * priorSd) + data.length * average(data) / (sigma * sigma))
                / precision;
        Normal exact = new Normal(posteriorMean, Math.sqrt(1.0 / precision));
        assertEquals("the mean differs", exact.mean(), p.mean(), 1.0e-10);
        assertEquals("the variance differs", exact.variance(), p.variance(), 1.0e-11);
        for (double q : new double[] { 0.025, 0.5, 0.975 }) {
            assertEquals("the quantile differs at " + q, exact.inverseCdf(q), p.inverseCdf(q), 1.0e-9);
        }
    }

    @Test
    public void andItStillWorksWithFiveThousandObservations() {
        // the conjugate answer does not care how many observations there are,
        // so this is an exact statement about the shift at a sample size where
        // the unshifted integrand is identically zero
        final double shape = 2.0;
        final double scale = 1.0;
        final int[] counts = poissonSample(5000, 4.0);
        final Gamma prior = new Gamma(shape, scale);
        DFunction logPosterior = new DFunction() {
            @Override
            public double apply(double lambda) {
                double value = prior.logPdf(lambda);
                Poisson likelihood = new Poisson(lambda);
                for (int i = 0; i < counts.length; i++) {
                    value += likelihood.logPmf(counts[i]);
                }
                return value;
            }
        };
        assertEquals("the unshifted density has not underflowed, so this test proves nothing", 0.0,
                Math.exp(logPosterior.apply(4.0)), 0.0);

        ScalarPosterior p = ScalarPosterior.of(logPosterior, 0.0, INF, 4.0, 1.0e-10);
        long total = 0L;
        for (int i = 0; i < counts.length; i++) {
            total += counts[i];
        }
        Gamma exact = new Gamma(shape + total, scale / (1.0 + counts.length * scale));
        assertEquals("the mean differs", exact.mean(), p.mean(), 1.0e-9 * exact.mean());
        assertEquals("the variance differs", exact.variance(), p.variance(), 1.0e-8 * exact.variance());
        for (double q : new double[] { 0.025, 0.5, 0.975 }) {
            assertEquals("the quantile differs at " + q, exact.inverseCdf(q), p.inverseCdf(q),
                    1.0e-9 * exact.mean());
        }
        assertTrue("the two routes disagree : " + p.errorEstimate(), p.errorEstimate() < 1.0e-10);
    }

    // ------------------------------------------------------------------
    // helpers
    // ------------------------------------------------------------------

    private static double logBeta(double x, double y) {
        return FastGamma.logGamma(x) + FastGamma.logGamma(y) - FastGamma.logGamma(x + y);
    }

    /** The log likelihood of {@code data} for a normal mean, up to a constant. */
    private static DFunction normalLogLikelihood(final double[] data, final double sigma) {
        return new DFunction() {
            @Override
            public double apply(double mu) {
                double sum = 0.0;
                for (int i = 0; i < data.length; i++) {
                    double d = data[i] - mu;
                    sum += d * d;
                }
                return -sum / (2.0 * sigma * sigma);
            }
        };
    }

    private static double average(double[] data) {
        double sum = 0.0;
        for (int i = 0; i < data.length; i++) {
            sum += data[i];
        }
        return sum / data.length;
    }

    /** Box-Muller over the LCG the tests in this repository use. */
    private static double[] normalSample(int n, double mean, double sd) {
        long lcg = 20260827L;
        double[] data = new double[n];
        for (int i = 0; i < n; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double u1 = ((lcg >>> 11) + 0.5) * 0x1.0p-53;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double u2 = ((lcg >>> 11) + 0.5) * 0x1.0p-53;
            data[i] = mean + sd * Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        }
        return data;
    }

    /** Knuth's multiplication method, which is enough for a mean of four. */
    private static int[] poissonSample(int n, double mean) {
        long lcg = 987654321L;
        double limit = Math.exp(-mean);
        int[] counts = new int[n];
        for (int i = 0; i < n; i++) {
            int k = 0;
            double product = 1.0;
            while (true) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                product *= ((lcg >>> 11) + 0.5) * 0x1.0p-53;
                if (product <= limit) {
                    break;
                }
                k++;
            }
            counts[i] = k;
        }
        return counts;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static void refuses(String what, DFunction f, double lo, double hi, double start, double epsTol) {
        try {
            ScalarPosterior.of(f, lo, hi, start, epsTol);
            fail("the factory accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue(what + ": threw without a message", expected.getMessage() != null);
        }
    }
}
