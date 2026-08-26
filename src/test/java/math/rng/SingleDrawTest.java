package math.rng;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

import org.junit.Test;

import math.distribution.Beta;
import math.distribution.ContinuousDistribution;
import math.distribution.FisherF;
import math.distribution.InverseGamma;
import math.distribution.Uniform;
import math.stats.Alternative;
import math.stats.HypothesisTests;

/**
 * The single-draw entry points on {@link PseudoRandom}.
 * <p>
 * Fifteen of the seventeen families are held to a far stronger statement than
 * any distributional test could make: from the same seed the single draw
 * must be the <b>identical</b> {@code double} the stream would produce. That
 * makes the evidence the stream tests already carry transfer whole, and no
 * family needs its distribution established a second time.
 * <p>
 * Beta and F cannot be held to it, because they deliberately take both of their
 * underlying gamma variates from one generator where the stream splits a second
 * one off. Those two are checked against {@link Beta} and {@link FisherF}
 * instead, and not by one fit but by the uniformity of many.
 */
public final class SingleDrawTest {

    private static final long[] SEEDS = { 1L, 42L, -7L, 0x2545F4914F6CDD1DL, 987654321L };
    private static final int DRAWS = 500;

    // ------------------------------------------------------------------------
    // ten families, bit for bit
    // ------------------------------------------------------------------------

    private static void identical(String what, ToDoubleFunction<PseudoRandom> single,
            Function<PseudoRandom, DoubleStream> stream) {
        for (int s = 0; s < SEEDS.length; s++) {
            PseudoRandom one = DefaultRng.newPseudoRandom(SEEDS[s]);
            PseudoRandom two = DefaultRng.newPseudoRandom(SEEDS[s]);
            for (int i = 0; i < DRAWS; i++) {
                double direct = single.applyAsDouble(one);
                double viaStream = stream.apply(two).toArray()[0];
                assertEquals(what + ", seed " + SEEDS[s] + ", draw " + i, viaStream, direct, 0.0);
            }
        }
    }

    @Test
    public void testCauchyIsBitForBitTheStream() {
        identical("cauchy(0, 1)", p -> p.nextCauchy(0.0, 1.0), p -> p.cauchy(1L, 0.0, 1.0));
        identical("cauchy(-3.5, 0.25)", p -> p.nextCauchy(-3.5, 0.25), p -> p.cauchy(1L, -3.5, 0.25));
        identical("cauchy(1e6, 1e-6)", p -> p.nextCauchy(1.0e6, 1.0e-6), p -> p.cauchy(1L, 1.0e6, 1.0e-6));
    }

    @Test
    public void testExponentialIsBitForBitTheStream() {
        identical("exponential(1)", p -> p.nextExponential(1.0), p -> p.exponential(1L, 1.0));
        identical("exponential(0.001)", p -> p.nextExponential(0.001), p -> p.exponential(1L, 0.001));
        identical("exponential(1000)", p -> p.nextExponential(1000.0), p -> p.exponential(1L, 1000.0));
    }

    @Test
    public void testGammaIsBitForBitTheStream() {
        // both sides of the shape-one boundary, where the algorithm changes
        identical("gamma(0.01, 1)", p -> p.nextGamma(0.01, 1.0), p -> p.gamma(1L, 0.01, 1.0));
        identical("gamma(0.9, 2)", p -> p.nextGamma(0.9, 2.0), p -> p.gamma(1L, 0.9, 2.0));
        identical("gamma(1, 1)", p -> p.nextGamma(1.0, 1.0), p -> p.gamma(1L, 1.0, 1.0));
        identical("gamma(7.5, 0.25)", p -> p.nextGamma(7.5, 0.25), p -> p.gamma(1L, 7.5, 0.25));
        identical("gamma(500, 3)", p -> p.nextGamma(500.0, 3.0), p -> p.gamma(1L, 500.0, 3.0));
    }

    @Test
    public void testChiSquareIsBitForBitTheStream() {
        identical("chiSquare(1)", p -> p.nextChiSquare(1.0), p -> p.chiSquare(1L, 1.0));
        identical("chiSquare(2.5)", p -> p.nextChiSquare(2.5), p -> p.chiSquare(1L, 2.5));
        identical("chiSquare(1000)", p -> p.nextChiSquare(1000.0), p -> p.chiSquare(1L, 1000.0));
    }

    @Test
    public void testLogNormalIsBitForBitTheStream() {
        identical("logNormal(0, 1)", p -> p.nextLogNormal(0.0, 1.0), p -> p.logNormal(1L, 0.0, 1.0));
        identical("logNormal(-2, 0.5)", p -> p.nextLogNormal(-2.0, 0.5), p -> p.logNormal(1L, -2.0, 0.5));
    }

    @Test
    public void testStudentTIsBitForBitTheStream() {
        identical("studentT(1)", p -> p.nextStudentT(1.0), p -> p.studentT(1L, 1.0));
        identical("studentT(4.5)", p -> p.nextStudentT(4.5), p -> p.studentT(1L, 4.5));
        identical("studentT(200)", p -> p.nextStudentT(200.0), p -> p.studentT(1L, 200.0));
    }

    @Test
    public void testWeibullIsBitForBitTheStream() {
        identical("weibull(1, 1)", p -> p.nextWeibull(1.0, 1.0), p -> p.weibull(1L, 1.0, 1.0));
        identical("weibull(2.5, 0.3)", p -> p.nextWeibull(2.5, 0.3), p -> p.weibull(1L, 2.5, 0.3));
    }

    @Test
    public void testTruncatedStandardNormalIsBitForBitTheStream() {
        // the symmetric case takes the reflected branch in both paths, and the
        // wide one is where the error function has saturated and the clamp acts
        identical("truncated(-2, 2)", p -> p.nextTruncatedStandardNormal(-2.0, 2.0),
                p -> p.truncatedStandardNormal(1L, -2.0, 2.0));
        identical("truncated(0.5, 3)", p -> p.nextTruncatedStandardNormal(0.5, 3.0),
                p -> p.truncatedStandardNormal(1L, 0.5, 3.0));
        identical("truncated(-9, 9)", p -> p.nextTruncatedStandardNormal(-9.0, 9.0),
                p -> p.truncatedStandardNormal(1L, -9.0, 9.0));
    }

    @Test
    public void testLeCunNormalIsBitForBitTheStream() {
        identical("leCun(1)", p -> p.nextLeCunNormal(1.0), p -> p.leCunNormal(1L, 1.0));
        identical("leCun(0.04)", p -> p.nextLeCunNormal(0.04), p -> p.leCunNormal(1L, 0.04));
    }

    @Test
    public void testInverseGammaIsBitForBitTheStream() {
        // the single draw takes the scale as stated and the stream inverts it in
        // its constructor; the two have to arrive at the same reciprocal
        identical("inverseGamma(3, 2)", p -> p.nextInverseGamma(3.0, 2.0), p -> p.inverseGamma(1L, 3.0, 2.0));
        identical("inverseGamma(0.5, 7)", p -> p.nextInverseGamma(0.5, 7.0), p -> p.inverseGamma(1L, 0.5, 7.0));
    }

    // ------------------------------------------------------------------------
    // the five discrete families, bit for bit as well
    // ------------------------------------------------------------------------

    private static void identicalCount(String what, ToIntFunction<PseudoRandom> single,
            Function<PseudoRandom, IntStream> stream) {
        for (int s = 0; s < SEEDS.length; s++) {
            PseudoRandom one = DefaultRng.newPseudoRandom(SEEDS[s]);
            PseudoRandom two = DefaultRng.newPseudoRandom(SEEDS[s]);
            for (int i = 0; i < DRAWS; i++) {
                int direct = single.applyAsInt(one);
                int viaStream = stream.apply(two).toArray()[0];
                assertEquals(what + ", seed " + SEEDS[s] + ", draw " + i, viaStream, direct);
            }
        }
    }

    @Test
    public void testPoissonIsBitForBitTheStream() {
        // the stream keeps its derived constants where the single draw builds
        // them per call, so this is a statement about those constants as much
        // as about the draw. Both sides of the crossover from multiplication to
        // rejection are covered
        for (double lambda : new double[] { 0.5, 3.0, 29.0, 30.0, 45.0, 500.0, 700.0 }) {
            identicalCount("poisson(" + lambda + ")", p -> p.nextPoisson(lambda),
                    p -> p.poisson(1L, lambda));
        }
    }

    @Test
    public void testBinomialIsBitForBitTheStream() {
        // n p on both sides of the crossover, and p on both sides of one half,
        // where the draw reflects
        int[][] settings = { { 1, 5 }, { 10, 30 }, { 20, 50 }, { 50, 50 }, { 300, 40 }, { 1000, 90 } };
        for (int s = 0; s < settings.length; s++) {
            int n = settings[s][0];
            double p = settings[s][1] / 100.0;
            identicalCount("binomial(" + n + ", " + p + ")", g -> g.nextBinomial(n, p),
                    g -> g.binomial(1L, n, p));
        }
    }

    @Test
    public void testGeometricIsBitForBitTheStream() {
        for (double p : new double[] { 0.01, 0.2, 0.5, 0.9, 1.0 }) {
            identicalCount("geometric(" + p + ")", g -> g.nextGeometric(p), g -> g.geometric(1L, p));
        }
    }

    @Test
    public void testNegativeBinomialIsBitForBitTheStream() {
        int[][] settings = { { 1, 30 }, { 5, 50 }, { 20, 10 }, { 100, 80 } };
        for (int s = 0; s < settings.length; s++) {
            int r = settings[s][0];
            double p = settings[s][1] / 100.0;
            identicalCount("negativeBinomial(" + r + ", " + p + ")", g -> g.nextNegativeBinomial(r, p),
                    g -> g.negativeBinomial(1L, r, p));
        }
    }

    @Test
    public void testHypergeometricIsBitForBitTheStream() {
        int[][] settings = { { 50, 20, 10 }, { 1000, 300, 100 }, { 10, 5, 5 }, { 100000, 50000, 2 } };
        for (int s = 0; s < settings.length; s++) {
            int population = settings[s][0];
            int successes = settings[s][1];
            int draws = settings[s][2];
            identicalCount("hypergeometric(" + population + ", " + successes + ", " + draws + ")",
                    g -> g.nextHypergeometric(population, successes, draws),
                    g -> g.hypergeometric(1L, population, successes, draws));
        }
    }

    // ------------------------------------------------------------------------
    // the two that cannot be, and are checked against their distribution
    // ------------------------------------------------------------------------

    private static final int FITS = 60;
    private static final int FIT_DRAWS = 2000;

    private static void fitsUniformly(String what, ContinuousDistribution law, long start,
            ToDoubleFunction<PseudoRandom> draw) {
        double[] p = new double[FITS];
        long lcg = start;
        for (int f = 0; f < FITS; f++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            PseudoRandom prng = DefaultRng.newPseudoRandom(lcg);
            double[] x = new double[FIT_DRAWS];
            for (int i = 0; i < FIT_DRAWS; i++) {
                x[i] = draw.applyAsDouble(prng);
            }
            p[f] = HypothesisTests.kolmogorovSmirnov(x, law, Alternative.TWO_SIDED).pValue;
        }
        double uniformity = HypothesisTests.kolmogorovSmirnov(p, new Uniform(0.0, 1.0), Alternative.TWO_SIDED).pValue;
        // measured over these exact settings: 0.14 at worst, 0.99 at best. A
        // single draw that took its two gamma variates from a state the other
        // had already touched would not land in that range
        assertTrue(what + ": the p-values of " + FITS + " fits are not uniform, p = " + uniformity,
                uniformity > 0.01);
        double[] sorted = p.clone();
        Arrays.sort(sorted);
        assertTrue(what + ": the median fit is " + sorted[FITS / 2], sorted[FITS / 2] > 0.2);
    }

    @Test
    public void testBetaFollowsTheBetaDistribution() {
        double[][] shapes = { { 0.5, 0.5 }, { 1.0, 1.0 }, { 2.0, 5.0 }, { 0.1, 3.0 }, { 30.0, 30.0 } };
        for (int s = 0; s < shapes.length; s++) {
            double alpha = shapes[s][0];
            double beta = shapes[s][1];
            fitsUniformly("Beta(" + alpha + ", " + beta + ")", new Beta(alpha, beta), 0x2545F4914F6CDD1DL,
                    p -> p.nextBeta(alpha, beta));
        }
    }

    @Test
    public void testFisherFFollowsTheFDistribution() {
        int[][] dfs = { { 1, 1 }, { 2, 5 }, { 10, 10 }, { 3, 100 }, { 100, 3 } };
        for (int s = 0; s < dfs.length; s++) {
            int d1 = dfs[s][0];
            int d2 = dfs[s][1];
            fitsUniformly("F(" + d1 + ", " + d2 + ")", new FisherF(d1, d2), 0x9E3779B97F4A7C15L,
                    p -> p.nextFisherF(d1, d2));
        }
    }

    @Test
    public void testInverseGammaFollowsTheInverseGammaDistribution() {
        // the stream and the distribution class were written years apart and
        // in different packages: the sampler inverts a gamma variate, the
        // class integrates an incomplete gamma. Nothing but the parameter
        // convention is shared, and that is exactly what this pins
        double[][] shapes = { { 0.5, 1.0 }, { 1.0, 1.0 }, { 2.5, 3.0 }, { 7.0, 0.25 }, { 30.0, 30.0 } };
        for (int s = 0; s < shapes.length; s++) {
            double alpha = shapes[s][0];
            double beta = shapes[s][1];
            fitsUniformly("InverseGamma(" + alpha + ", " + beta + ")", new InverseGamma(alpha, beta),
                    0x14057B7EF767814FL, p -> p.nextInverseGamma(alpha, beta));
        }
    }

    @Test
    public void testConsecutiveBetaDrawsAreNotSeriallyCorrelated() {
        // the one thing a fit cannot see: with a single generator the first
        // variate of draw i + 1 comes from the state the second variate of draw
        // i left behind. Measured over 400000 draws the worst lag-one
        // correlation was 2.2e-3 against a standard error of 1.6e-3
        int n = 200000;
        double[][] shapes = { { 0.5, 0.5 }, { 2.0, 5.0 }, { 30.0, 30.0 } };
        for (int s = 0; s < shapes.length; s++) {
            PseudoRandom prng = DefaultRng.newPseudoRandom(31L + s);
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = prng.nextBeta(shapes[s][0], shapes[s][1]);
            }
            double r = lagOne(x);
            // 2.5 standard errors at this sample size
            assertTrue("Beta(" + shapes[s][0] + ", " + shapes[s][1] + ") lag-one correlation is " + r,
                    Math.abs(r) < 2.5 / Math.sqrt(n));
        }
    }

    private static double lagOne(double[] x) {
        int n = x.length;
        double mean = 0.0;
        for (int i = 0; i < n; i++) {
            mean += x[i];
        }
        mean /= n;
        double numerator = 0.0;
        double denominator = 0.0;
        for (int i = 0; i < n; i++) {
            double d = x[i] - mean;
            denominator += d * d;
            if (i + 1 < n) {
                numerator += d * (x[i + 1] - mean);
            }
        }
        return numerator / denominator;
    }

    // ------------------------------------------------------------------------
    // whatever the stream rejects, the single draw rejects
    // ------------------------------------------------------------------------

    private static void bothReject(String what, Runnable single, Runnable stream) {
        String fromSingle = throwsWith(single);
        String fromStream = throwsWith(stream);
        if (fromSingle == null) {
            fail(what + ": the single draw accepted it, the stream said \"" + fromStream + "\"");
        }
        if (fromStream == null) {
            fail(what + ": the stream accepted it, the single draw said \"" + fromSingle + "\"");
        }
        // Both refuse, which is the parity that matters. The wording is not
        // asserted to agree: AbstractRng64 range-checks some arguments itself
        // before it ever builds a spliterator, so for those the stream answers
        // with its own message and the single draw with the spliterator's
        assertTrue(what + ": the single draw threw without a message", fromSingle.length() > 0);
        assertTrue(what + ": the stream threw without a message", fromStream.length() > 0);
    }

    private static String throwsWith(Runnable r) {
        try {
            r.run();
            return null;
        } catch (IllegalArgumentException e) {
            return e.getMessage();
        }
    }

    @Test
    public void testTheSingleDrawRejectsExactlyWhatTheStreamRejects() {
        PseudoRandom p = DefaultRng.newPseudoRandom(5L);
        for (double bad : new double[] { 0.0, -1.0, -1.0e-300, Double.NEGATIVE_INFINITY,
                Double.NaN }) {
            bothReject("cauchy scale " + bad, () -> p.nextCauchy(0.0, bad), () -> p.cauchy(1L, 0.0, bad));
            bothReject("exponential rate " + bad, () -> p.nextExponential(bad), () -> p.exponential(1L, bad));
            bothReject("gamma shape " + bad, () -> p.nextGamma(bad, 1.0), () -> p.gamma(1L, bad, 1.0));
            bothReject("gamma scale " + bad, () -> p.nextGamma(1.0, bad), () -> p.gamma(1L, 1.0, bad));
            bothReject("beta alpha " + bad, () -> p.nextBeta(bad, 1.0), () -> p.beta(1L, bad, 1.0));
            bothReject("beta beta " + bad, () -> p.nextBeta(1.0, bad), () -> p.beta(1L, 1.0, bad));
            bothReject("chiSquare df " + bad, () -> p.nextChiSquare(bad), () -> p.chiSquare(1L, bad));
            bothReject("logNormal sigma " + bad, () -> p.nextLogNormal(0.0, bad), () -> p.logNormal(1L, 0.0, bad));
            bothReject("studentT df " + bad, () -> p.nextStudentT(bad), () -> p.studentT(1L, bad));
            bothReject("weibull scale " + bad, () -> p.nextWeibull(bad, 1.0), () -> p.weibull(1L, bad, 1.0));
            bothReject("weibull shape " + bad, () -> p.nextWeibull(1.0, bad), () -> p.weibull(1L, 1.0, bad));
            bothReject("leCun sigma " + bad, () -> p.nextLeCunNormal(bad), () -> p.leCunNormal(1L, bad));
            bothReject("inverseGamma alpha " + bad, () -> p.nextInverseGamma(bad, 1.0),
                    () -> p.inverseGamma(1L, bad, 1.0));
            bothReject("inverseGamma beta " + bad, () -> p.nextInverseGamma(1.0, bad),
                    () -> p.inverseGamma(1L, 1.0, bad));
        }
        for (int bad : new int[] { 0, -1, Integer.MIN_VALUE }) {
            bothReject("fisherF numerator " + bad, () -> p.nextFisherF(bad, 1), () -> p.fisherF(1L, bad, 1));
            bothReject("fisherF denominator " + bad, () -> p.nextFisherF(1, bad), () -> p.fisherF(1L, 1, bad));
        }
        for (double[] bounds : new double[][] { { 1.0, 1.0 }, { 2.0, -2.0 }, { 0.0, Double.NaN },
                { Double.NaN, 0.0 } }) {
            double lower = bounds[0];
            double upper = bounds[1];
            bothReject("truncated (" + lower + ", " + upper + ")",
                    () -> p.nextTruncatedStandardNormal(lower, upper),
                    () -> p.truncatedStandardNormal(1L, lower, upper));
        }

        // the discrete five, whose checks were written later and do name NaN
        for (double bad : new double[] { -1.0, Double.NaN, 1.0e30 }) {
            bothReject("poisson mean " + bad, () -> p.nextPoisson(bad), () -> p.poisson(1L, bad));
        }
        for (double bad : new double[] { -0.5, 1.5, Double.NaN }) {
            bothReject("binomial p " + bad, () -> p.nextBinomial(10, bad), () -> p.binomial(1L, 10, bad));
            bothReject("geometric p " + bad, () -> p.nextGeometric(bad), () -> p.geometric(1L, bad));
            bothReject("negativeBinomial p " + bad, () -> p.nextNegativeBinomial(3, bad),
                    () -> p.negativeBinomial(1L, 3, bad));
        }
        bothReject("binomial n -1", () -> p.nextBinomial(-1, 0.5), () -> p.binomial(1L, -1, 0.5));
        bothReject("geometric p 0", () -> p.nextGeometric(0.0), () -> p.geometric(1L, 0.0));
        bothReject("negativeBinomial r 0", () -> p.nextNegativeBinomial(0, 0.5),
                () -> p.negativeBinomial(1L, 0, 0.5));
        bothReject("hypergeometric population -1", () -> p.nextHypergeometric(-1, 0, 0),
                () -> p.hypergeometric(1L, -1, 0, 0));
        bothReject("hypergeometric successes above the population", () -> p.nextHypergeometric(10, 11, 5),
                () -> p.hypergeometric(1L, 10, 11, 5));
        bothReject("hypergeometric draws above the population", () -> p.nextHypergeometric(10, 5, 11),
                () -> p.hypergeometric(1L, 10, 5, 11));
    }

    @Test
    public void testEveryFamilyRefusesANaNParameter() {
        // A "<= 0.0" test is false for NaN, so for a long time a NaN shape or
        // scale travelled through this package and came back as a NaN variate
        // with no diagnostic at all. Every such check is now written
        // "!(x > 0.0)", which is false for NaN as well, and this is the test
        // that says so for both paths at once -- there is no family left that
        // answers the question differently
        PseudoRandom p = DefaultRng.newPseudoRandom(11L);
        double n = Double.NaN;

        rejectsNaN("nextGaussian", () -> p.nextGaussian(0.0, n));
        rejectsNaN("normal stream", () -> p.normal(1L, 0.0, n).toArray());
        rejectsNaN("nextCauchy", () -> p.nextCauchy(0.0, n));
        rejectsNaN("cauchy stream", () -> p.cauchy(1L, 0.0, n).toArray());
        rejectsNaN("nextExponential", () -> p.nextExponential(n));
        rejectsNaN("exponential stream", () -> p.exponential(1L, n).toArray());
        rejectsNaN("nextGamma shape", () -> p.nextGamma(n, 1.0));
        rejectsNaN("nextGamma scale", () -> p.nextGamma(1.0, n));
        rejectsNaN("gamma stream", () -> p.gamma(1L, n, 1.0).toArray());
        rejectsNaN("nextBeta alpha", () -> p.nextBeta(n, 1.0));
        rejectsNaN("nextBeta beta", () -> p.nextBeta(1.0, n));
        rejectsNaN("beta stream", () -> p.beta(1L, n, 1.0).toArray());
        rejectsNaN("nextChiSquare", () -> p.nextChiSquare(n));
        rejectsNaN("chiSquare stream", () -> p.chiSquare(1L, n).toArray());
        rejectsNaN("nextLogNormal", () -> p.nextLogNormal(0.0, n));
        rejectsNaN("logNormal stream", () -> p.logNormal(1L, 0.0, n).toArray());
        rejectsNaN("nextStudentT", () -> p.nextStudentT(n));
        rejectsNaN("studentT stream", () -> p.studentT(1L, n).toArray());
        rejectsNaN("nextWeibull scale", () -> p.nextWeibull(n, 1.0));
        rejectsNaN("nextWeibull shape", () -> p.nextWeibull(1.0, n));
        rejectsNaN("weibull stream", () -> p.weibull(1L, n, 1.0).toArray());
        rejectsNaN("nextLeCunNormal", () -> p.nextLeCunNormal(n));
        rejectsNaN("leCunNormal stream", () -> p.leCunNormal(1L, n).toArray());
        rejectsNaN("nextInverseGamma alpha", () -> p.nextInverseGamma(n, 1.0));
        rejectsNaN("nextInverseGamma beta", () -> p.nextInverseGamma(1.0, n));
        rejectsNaN("inverseGamma stream", () -> p.inverseGamma(1L, n, 1.0).toArray());
        rejectsNaN("nextTruncatedStandardNormal", () -> p.nextTruncatedStandardNormal(n, 1.0));
        rejectsNaN("nextPoisson", () -> p.nextPoisson(n));
        rejectsNaN("nextBinomial", () -> p.nextBinomial(10, n));
        rejectsNaN("nextGeometric", () -> p.nextGeometric(n));
        rejectsNaN("nextNegativeBinomial", () -> p.nextNegativeBinomial(3, n));

        // and the draw still works with the same parameter made finite, so the
        // tightening did not simply refuse everything
        finite("nextCauchy", p.nextCauchy(0.0, 1.0));
        finite("nextGamma", p.nextGamma(1.0, 1.0));
        finite("nextWeibull", p.nextWeibull(1.0, 1.0));
    }

    private static void rejectsNaN(String what, Runnable r) {
        try {
            r.run();
            fail(what + " accepted a NaN parameter");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + " threw without a message", expected.getMessage() != null);
        }
    }

    // ------------------------------------------------------------------------
    // the bodies really are inherited, and beta no longer needs two generators
    // ------------------------------------------------------------------------

    @Test
    public void testAGeneratorFromOutsideTheHierarchyGetsAllTwelve() {
        // Delegate implements PseudoRandom directly: it is neither a
        // SplittablePseudoRandom nor an AbstractRng64, and it refuses to hand
        // out a stream at all. If the twelve draws work on it, they are
        // inherited from the interface and nothing had to be added to any
        // implementation
        PseudoRandom foreign = new Delegate(DefaultRng.newPseudoRandom(99L));
        assertTrue("the test generator can be detached after all, which would make this test vacuous",
                PseudoRandomSpliterator.detach(foreign) == null);

        finite("nextCauchy", foreign.nextCauchy(0.0, 1.0));
        finite("nextExponential", foreign.nextExponential(2.0));
        finite("nextGamma", foreign.nextGamma(3.0, 1.5));
        finite("nextBeta", foreign.nextBeta(2.0, 3.0));
        finite("nextChiSquare", foreign.nextChiSquare(4.0));
        finite("nextFisherF", foreign.nextFisherF(3, 7));
        finite("nextLogNormal", foreign.nextLogNormal(0.0, 1.0));
        finite("nextStudentT", foreign.nextStudentT(5.0));
        finite("nextWeibull", foreign.nextWeibull(1.0, 2.0));
        finite("nextTruncatedStandardNormal", foreign.nextTruncatedStandardNormal(-1.0, 1.0));
        finite("nextLeCunNormal", foreign.nextLeCunNormal(0.5));
        finite("nextInverseGamma", foreign.nextInverseGamma(3.0, 2.0));

        assertTrue("nextPoisson", foreign.nextPoisson(4.0) >= 0);
        assertTrue("nextBinomial", foreign.nextBinomial(20, 0.3) >= 0);
        assertTrue("nextGeometric", foreign.nextGeometric(0.25) >= 0);
        assertTrue("nextNegativeBinomial", foreign.nextNegativeBinomial(4, 0.4) >= 0);
        assertTrue("nextHypergeometric", foreign.nextHypergeometric(50, 20, 10) >= 0);
    }

    @Test
    public void testBetaAndFisherFWorkWhereTheirStreamCannot() {
        // the stream needs a second generator and says so; the single draw does
        // not, which is the whole reason it takes both variates from one
        PseudoRandom foreign = new Delegate(DefaultRng.newPseudoRandom(4L));
        for (int i = 0; i < 100; i++) {
            double b = foreign.nextBeta(2.0, 5.0);
            assertTrue("a beta draw left [0, 1]: " + b, b >= 0.0 && b <= 1.0);
            double f = foreign.nextFisherF(4, 9);
            assertTrue("an F draw is not positive: " + f, f >= 0.0 && !Double.isNaN(f));
        }
        try {
            new BetaSpliterator(foreign, 0L, 10L, 2.0, 5.0);
            fail("the beta stream accepted a generator it cannot split");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without naming the generator",
                    expected.getMessage().contains("neither splittable nor re-creatable"));
        }
    }

    @Test
    public void testTheSingleDrawsAreReproducibleFromTheSeed() {
        PseudoRandom one = DefaultRng.newPseudoRandom(2024L);
        PseudoRandom two = DefaultRng.newPseudoRandom(2024L);
        for (int i = 0; i < 300; i++) {
            assertEquals("beta at " + i, one.nextBeta(1.5, 2.5), two.nextBeta(1.5, 2.5), 0.0);
            assertEquals("fisherF at " + i, one.nextFisherF(3, 4), two.nextFisherF(3, 4), 0.0);
            assertEquals("gamma at " + i, one.nextGamma(0.7, 2.0), two.nextGamma(0.7, 2.0), 0.0);
        }
    }

    private static void finite(String what, double value) {
        assertTrue(what + " returned " + value, !Double.isNaN(value) && !Double.isInfinite(value));
    }

    /**
     * A generator that implements {@link PseudoRandom} and nothing else: it
     * forwards the primitives and has no streams at all. Standing outside both
     * {@code SplittablePseudoRandom} and {@code AbstractRng64}, it is the case
     * the single draws had to be made to work for.
     */
    private static final class Delegate implements PseudoRandom {

        private final PseudoRandom d;

        Delegate(PseudoRandom d) {
            this.d = d;
        }

        private static UnsupportedOperationException noStreams() {
            return new UnsupportedOperationException("this generator has no streams");
        }

        @Override
        public IntStream ints() {
            throw noStreams();
        }

        @Override
        public IntStream ints(long streamSize) {
            throw noStreams();
        }

        @Override
        public IntStream ints(long streamSize, int min, int max) {
            throw noStreams();
        }

        @Override
        public IntStream ints(int min, int max) {
            throw noStreams();
        }

        @Override
        public LongStream longs() {
            throw noStreams();
        }

        @Override
        public LongStream longs(long streamSize) {
            throw noStreams();
        }

        @Override
        public LongStream longs(long streamSize, long min, long max) {
            throw noStreams();
        }

        @Override
        public LongStream longs(long min, long max) {
            throw noStreams();
        }

        @Override
        public DoubleStream doubles() {
            throw noStreams();
        }

        @Override
        public DoubleStream doubles(long streamSize) {
            throw noStreams();
        }

        @Override
        public DoubleStream doubles(long streamSize, double min, double max) {
            throw noStreams();
        }

        @Override
        public DoubleStream doubles(double min, double max) {
            throw noStreams();
        }

        @Override
        public DoubleStream normal(double mu, double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream normal(long streamSize, double mu, double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream cauchy(double location, double scale) {
            throw noStreams();
        }

        @Override
        public DoubleStream cauchy(long streamSize, double location, double scale) {
            throw noStreams();
        }

        @Override
        public DoubleStream exponential(double lambda) {
            throw noStreams();
        }

        @Override
        public DoubleStream exponential(long streamSize, double lambda) {
            throw noStreams();
        }

        @Override
        public DoubleStream gamma(double k, double theta) {
            throw noStreams();
        }

        @Override
        public DoubleStream gamma(long streamSize, double k, double theta) {
            throw noStreams();
        }

        @Override
        public DoubleStream beta(double alpha, double beta) {
            throw noStreams();
        }

        @Override
        public DoubleStream beta(long streamSize, double alpha, double beta) {
            throw noStreams();
        }

        @Override
        public DoubleStream chiSquare(double k) {
            throw noStreams();
        }

        @Override
        public DoubleStream chiSquare(long streamSize, double k) {
            throw noStreams();
        }

        @Override
        public DoubleStream fisherF(int numeratorDF, int denominatorDF) {
            throw noStreams();
        }

        @Override
        public DoubleStream fisherF(long streamSize, int numeratorDF, int denominatorDF) {
            throw noStreams();
        }

        @Override
        public DoubleStream logNormal(double mu, double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream logNormal(long streamSize, double mu, double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream studentT(double df) {
            throw noStreams();
        }

        @Override
        public DoubleStream studentT(long streamSize, double df) {
            throw noStreams();
        }

        @Override
        public DoubleStream weibull(double scale, double shape) {
            throw noStreams();
        }

        @Override
        public DoubleStream weibull(long streamSize, double scale, double shape) {
            throw noStreams();
        }

        @Override
        public DoubleStream truncatedStandardNormal(double min, double max) {
            throw noStreams();
        }

        @Override
        public DoubleStream truncatedStandardNormal(long streamSize, double min, double max) {
            throw noStreams();
        }

        @Override
        public DoubleStream leCunNormal(double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream leCunNormal(long streamSize, double sigma) {
            throw noStreams();
        }

        @Override
        public DoubleStream inverseGamma(double alpha, double beta) {
            throw noStreams();
        }

        @Override
        public DoubleStream inverseGamma(long streamSize, double alpha, double beta) {
            throw noStreams();
        }

        @Override
        public IntStream categorical(double[] weights) {
            throw noStreams();
        }

        @Override
        public IntStream categorical(long streamSize, double[] weights) {
            throw noStreams();
        }

        @Override
        public IntStream categorical(AliasTable table) {
            throw noStreams();
        }

        @Override
        public IntStream categorical(long streamSize, AliasTable table) {
            throw noStreams();
        }

        @Override
        public IntStream poisson(double lambda) {
            throw noStreams();
        }

        @Override
        public IntStream poisson(long streamSize, double lambda) {
            throw noStreams();
        }

        @Override
        public IntStream binomial(int n, double p) {
            throw noStreams();
        }

        @Override
        public IntStream binomial(long streamSize, int n, double p) {
            throw noStreams();
        }

        @Override
        public IntStream geometric(double p) {
            throw noStreams();
        }

        @Override
        public IntStream geometric(long streamSize, double p) {
            throw noStreams();
        }

        @Override
        public IntStream negativeBinomial(int r, double p) {
            throw noStreams();
        }

        @Override
        public IntStream negativeBinomial(long streamSize, int r, double p) {
            throw noStreams();
        }

        @Override
        public IntStream hypergeometric(int population, int successes, int draws) {
            throw noStreams();
        }

        @Override
        public IntStream hypergeometric(long streamSize, int population, int successes, int draws) {
            throw noStreams();
        }

        @Override
        public long nextLong() {
            return d.nextLong();
        }

        @Override
        public double nextDouble() {
            return d.nextDouble();
        }

        @Override
        public double nextDouble(double min, double max) {
            return d.nextDouble(min, max);
        }

        @Override
        public double nextGaussian() {
            return d.nextGaussian();
        }

        @Override
        public double nextGaussian(double mean, double stdDeviation) {
            return d.nextGaussian(mean, stdDeviation);
        }

        @Override
        public float nextFloat() {
            return d.nextFloat();
        }

        @Override
        public float nextFloat(float min, float max) {
            return d.nextFloat(min, max);
        }

        @Override
        public int nextInt() {
            return d.nextInt();
        }

        @Override
        public void nextBytes(byte[] bytes) {
            d.nextBytes(bytes);
        }

        @Override
        public boolean nextBoolean() {
            return d.nextBoolean();
        }

        @Override
        public long nextLong(long n) {
            return d.nextLong(n);
        }

        @Override
        public int nextInt(int n) {
            return d.nextInt(n);
        }

        @Override
        public int nextInt(int min, int max) {
            return d.nextInt(min, max);
        }

        @Override
        public long nextLong(long min, long max) {
            return d.nextLong(min, max);
        }

        @Override
        public int next(int bits) {
            return d.next(bits);
        }

        @Override
        public void nextLongs(long[] longs) {
            d.nextLongs(longs);
        }

        @Override
        public void nextDoubles(double[] doubles) {
            d.nextDoubles(doubles);
        }

        @Override
        public int[] intsSampledWithoutReplacement(int min, int max, int count) {
            return d.intsSampledWithoutReplacement(min, max, count);
        }

        @Override
        public String getAlgorithm() {
            return d.getAlgorithm();
        }

        @Override
        public long[] getSeed() {
            return d.getSeed();
        }
    }
}
