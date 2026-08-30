package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.solve.Quadrature;

/**
 * The inverse gamma against the gamma it is the reciprocal of, which is an
 * implementation already in the package and shares nothing with it but the
 * incomplete gamma function underneath.
 * <p>
 * The other cross-check lives in {@code math.rng.SingleDrawTest}, where
 * {@code PseudoRandom.nextInverseGamma} -- written years before this class and
 * in another package -- is fitted against it.
 */
public class InverseGammaTest {

    /** The law of {@code 1/X} for {@code X ~ InverseGamma(alpha, beta)}. */
    private static Gamma reciprocalOf(double alpha, double beta) {
        // Gamma takes a scale and the reciprocal has rate beta
        return new Gamma(alpha, 1.0 / beta);
    }

    private static final double[] SHAPES = { 0.25, 0.5, 1.0, 1.5, 2.5, 7.0, 40.0 };
    private static final double[] SCALES = { 0.25, 1.0, 3.0, 40.0 };

    @Test
    public void theDensityIsTheGammaDensityThroughTheJacobian() {
        // if X ~ InverseGamma(a, b) then Y = 1/X ~ Gamma(a, rate b), so
        // f_X(x) = f_Y(1/x) / x^2. The right-hand side is computed by a class
        // that knows nothing about this one
        double worst = 0.0;
        String at = "";
        for (double alpha : SHAPES) {
            for (double beta : SCALES) {
                InverseGamma d = new InverseGamma(alpha, beta);
                Gamma y = reciprocalOf(alpha, beta);
                for (double x = 1.0e-3; x < 1.0e3; x *= 1.4) {
                    double direct = d.pdf(x);
                    double jacobian = y.pdf(1.0 / x) / (x * x);
                    if (!(direct > 1.0e-280) || !(jacobian > 1.0e-280)) {
                        continue;
                    }
                    double relative = Math.abs(direct - jacobian) / direct;
                    if (relative > worst) {
                        worst = relative;
                        at = "alpha=" + alpha + " beta=" + beta + " x=" + x;
                    }
                }
            }
        }
        // measured over 112 settings: 1.3e-12, which is the gamma density's own
        // accuracy far out and not a disagreement about the law
        assertTrue("the two differ by " + worst + " relative, " + at, worst < 1.0e-10);
    }

    @Test
    public void theDistributionFunctionKeepsTheLowerTailASubtractionLoses() {
        // cdf(x) is Q(alpha, beta/x) taken directly. The obvious route, one
        // minus the reciprocal gamma's cdf, agrees wherever the subtraction can
        // afford it and returns exactly zero where it cannot -- which is most
        // of the lower tail
        InverseGamma d = new InverseGamma(2.5, 3.0);
        Gamma y = reciprocalOf(2.5, 3.0);

        // measured: 1.00e-62 and 3.14e-24 against an answer of exactly zero
        assertEquals("the subtraction should have nothing left here", 0.0, 1.0 - y.cdf(1.0 / 0.02), 0.0);
        assertTrue("cdf(0.02) is " + d.cdf(0.02), d.cdf(0.02) > 1.0e-63 && d.cdf(0.02) < 1.0e-61);
        assertEquals("nor here", 0.0, 1.0 - y.cdf(1.0 / 0.05), 0.0);
        assertTrue("cdf(0.05) is " + d.cdf(0.05), d.cdf(0.05) > 1.0e-24 && d.cdf(0.05) < 1.0e-23);

        // and where the subtraction does survive, the two agree
        double worst = 0.0;
        for (double alpha : SHAPES) {
            for (double beta : SCALES) {
                InverseGamma g = new InverseGamma(alpha, beta);
                Gamma r = reciprocalOf(alpha, beta);
                for (double x = 1.0e-3; x < 1.0e4; x *= 1.4) {
                    double viaSubtraction = 1.0 - r.cdf(1.0 / x);
                    if (viaSubtraction > 1.0e-8 && viaSubtraction < 1.0 - 1.0e-8) {
                        worst = Math.max(worst, Math.abs(g.cdf(x) - viaSubtraction) / viaSubtraction);
                    }
                }
            }
        }
        // measured: 1.7e-9, and the direct route is the accurate one of the two
        assertTrue("the two differ by " + worst + " relative", worst < 1.0e-7);
    }

    @Test
    public void theQuantileInvertsTheDistributionFunction() {
        double worst = 0.0;
        String at = "";
        for (double alpha : SHAPES) {
            for (double beta : SCALES) {
                InverseGamma d = new InverseGamma(alpha, beta);
                for (double p : new double[] { 1.0e-6, 1.0e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999,
                        1.0 - 1.0e-6 }) {
                    double x = d.inverseCdf(p);
                    assertTrue("alpha=" + alpha + " beta=" + beta + " p=" + p + ": the quantile is " + x,
                            x > 0.0 && x < Double.POSITIVE_INFINITY);
                    double relative = Math.abs(d.cdf(x) - p) / p;
                    if (relative > worst) {
                        worst = relative;
                        at = "alpha=" + alpha + " beta=" + beta + " p=" + p;
                    }
                }
            }
        }
        // measured: 9.2e-15, worst at alpha = 0.5 where no mean exists to start
        // the search from and it comes from the reciprocal gamma instead
        assertTrue("the round trip closes to only " + worst + " relative, " + at, worst < 1.0e-12);
    }

    @Test
    public void theScaleIsAScale() {
        // beta multiplies the variable, so cdf(x) at (alpha, beta) is cdf(c x)
        // at (alpha, c beta). With c a power of two the ratio beta/x inside is
        // formed from exactly scaled operands and comes out identical, so this
        // holds bit for bit rather than nearly
        for (double alpha : SHAPES) {
            for (double beta : SCALES) {
                for (int exponent : new int[] { -60, -8, -1, 1, 8, 60 }) {
                    double c = Math.scalb(1.0, exponent);
                    InverseGamma plain = new InverseGamma(alpha, beta);
                    InverseGamma scaled = new InverseGamma(alpha, c * beta);
                    for (double x = 1.0e-2; x < 1.0e2; x *= 3.0) {
                        String at = "alpha=" + alpha + " beta=" + beta + " 2^" + exponent + " x=" + x;
                        assertEquals(at, plain.cdf(x), scaled.cdf(c * x), 0.0);
                        // the density carries the Jacobian, so its logarithm
                        // differs by exactly log(c) -- and that one addition is
                        // why this is not bit for bit
                        assertEquals(at + ": the log density", plain.logPdf(x) - Math.log(c),
                                scaled.logPdf(c * x), 1.0e-12);
                    }
                }
            }
        }
    }

    @Test
    public void theMomentsAreTheIntegralsTheyClaimToBe() {
        // only where the tail is light enough for the quadrature to reach them:
        // an inverse gamma with a shape near one has a tail the rule cannot
        // integrate to anything like full precision, which is a fact about the
        // distribution rather than about this class
        for (double alpha : new double[] { 7.0, 40.0 }) {
            for (double beta : new double[] { 1.0, 3.0 }) {
                final InverseGamma d = new InverseGamma(alpha, beta);
                // the window has to hold the mass rather than merely contain
                // it: at alpha = 40 the density is a needle around 0.026, and
                // integrating it over (0, 4000] returned 1.1e-59 -- the rule
                // subdivided a domain a hundred thousand times too wide and
                // never found the peak
                double top = d.mean() + 1000.0 * Math.sqrt(d.variance());
                double mean = Quadrature.integrate(t -> t * d.pdf(t), 0.0, top, 1.0e-12);
                double variance = Quadrature.integrate(t -> (t - d.mean()) * (t - d.mean()) * d.pdf(t), 0.0, top,
                        1.0e-12);
                String at = "alpha=" + alpha + " beta=" + beta;
                // measured at 1000 standard deviations: 1.4e-14 on the mean
                // and 2.2e-11 on the variance. Sixty is not enough for the
                // lighter shape -- 1.2e-07 and 1.8e-05 -- since an inverse
                // gamma carries more of its second moment far out than a
                // multiple of its own spread suggests
                assertEquals(at + ": mean", d.mean(), mean, 1.0e-12 * Math.abs(d.mean()));
                assertEquals(at + ": variance", d.variance(), variance, 1.0e-9 * Math.abs(d.variance()));
            }
        }
    }

    @Test
    public void theMomentsThatDoNotExistAreNotReported() {
        // the FisherF convention: a moment that does not exist is NaN, not a
        // number that happens to come out of the formula
        for (double alpha : new double[] { 0.25, 0.5, 1.0 }) {
            InverseGamma d = new InverseGamma(alpha, 2.0);
            assertTrue("alpha=" + alpha + ": a mean was reported", Double.isNaN(d.mean()));
            assertTrue("alpha=" + alpha + ": a variance was reported", Double.isNaN(d.variance()));
        }
        for (double alpha : new double[] { 1.5, 2.0 }) {
            InverseGamma d = new InverseGamma(alpha, 2.0);
            assertTrue("alpha=" + alpha + ": no mean was reported", !Double.isNaN(d.mean()));
            assertTrue("alpha=" + alpha + ": a variance was reported", Double.isNaN(d.variance()));
        }
        InverseGamma both = new InverseGamma(3.0, 2.0);
        assertEquals("mean", 1.0, both.mean(), 1.0e-15);
        assertEquals("variance", 1.0, both.variance(), 1.0e-15);
    }

    @Test
    public void theVarianceSurvivesAScaleTheSquareWouldNot() {
        // beta^2 leaves the double range from about 1.3e154 upwards. That is
        // not always a problem the class can solve -- InverseGamma(3, 1e200)
        // has a variance of 2.5e399, which is no double whatever order the
        // arithmetic is done in -- but where the answer is representable and
        // the square of the scale alone is not, forming the ratio first is what
        // reaches it
        InverseGamma d = new InverseGamma(1.0e50, 1.0e200);
        double variance = d.variance();
        assertTrue("the variance came back as " + variance, variance > 0.0 && !Double.isInfinite(variance));
        // beta / (alpha - 1) is about 1e150, its square about 1e300, and
        // dividing by alpha - 2 lands at about 1e250
        assertTrue("the variance is " + variance, variance > 1.0e249 && variance < 1.0e251);
        assertEquals("the square of the scale is representable after all", Double.POSITIVE_INFINITY,
                1.0e200 * 1.0e200, 0.0);

        // and where it genuinely is not a double, an infinity is the honest
        // answer rather than a wrong finite one
        assertEquals("a variance of 2.5e399", Double.POSITIVE_INFINITY, new InverseGamma(3.0, 1.0e200).variance(),
                0.0);
    }

    @Test
    public void theParametersComeBackOut() {
        InverseGamma d = new InverseGamma(2.5, 7.0);
        assertEquals("shape", 2.5, d.getShape(), 0.0);
        assertEquals("scale", 7.0, d.getScale(), 0.0);
        assertEquals("the support starts at zero", 0.0, d.supportLowerBound(), 0.0);
        assertEquals("and does not end", Double.POSITIVE_INFINITY, d.supportUpperBound(), 0.0);
    }

    @Test
    public void theInverseGammaRejectsWhatItCannotBe() {
        refuses("a zero shape", 0.0, 1.0);
        refuses("a negative shape", -1.0, 1.0);
        refuses("an infinite shape", Double.POSITIVE_INFINITY, 1.0);
        refuses("a NaN shape", Double.NaN, 1.0);
        refuses("a zero scale", 1.0, 0.0);
        refuses("a negative scale", 1.0, -1.0);
        refuses("an infinite scale", 1.0, Double.POSITIVE_INFINITY);
        refuses("a NaN scale", 1.0, Double.NaN);
    }

    private static void refuses(String what, double shape, double scale) {
        try {
            new InverseGamma(shape, scale);
            fail("InverseGamma accepted " + what);
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without saying what was wrong", expected.getMessage() != null);
        }
    }
}
