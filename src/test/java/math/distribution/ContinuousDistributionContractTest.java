package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * The contract of {@link ContinuousDistribution} at the edges of the double
 * range, checked over every implementation at once rather than one class at a
 * time, so that a distribution added later without the rule fails here.
 * <p>
 * {@link ContinuousDistributionTest} is the numerical battery over windows in
 * the interior; this is only about totality: what happens outside the support,
 * at either infinity, at the top of the double range, and at
 * {@code Double.NaN}.
 */
public class ContinuousDistributionContractTest {

    private static final double NEG = Double.NEGATIVE_INFINITY;
    private static final double POS = Double.POSITIVE_INFINITY;
    private static final double MAX = Double.MAX_VALUE;

    /**
     * Every implementation, at the parameter settings that separate the
     * branches inside them: a shape below, at and above one, one and two
     * degrees of freedom in either position of {@code FisherF}, and a rate
     * large enough that {@code rate * x} overflows before {@code x} does.
     */
    private static ContinuousDistribution[] all() {
        return new ContinuousDistribution[] { new Beta(0.5, 0.5), new Beta(2.0, 3.0), new Beta(0.3, 4.0),
                new Cauchy(0.0, 1.0), new Cauchy(1000.0, 1.0), new ChiSquare(1.0), new ChiSquare(2.0),
                new ChiSquare(7.0), new Exponential(1.0), new Exponential(1.0e-9), new FisherF(1, 1),
                new FisherF(2, 2), new FisherF(5, 2), new FisherF(10, 20), new Gamma(0.5, 1.0), new Gamma(1.0, 1.0),
                new Gamma(2.0, 1.0), new Gamma(9.0, 0.5), new Gamma(2.0, 1.0e-9), new LogNormal(0.0, 1.0),
                new Normal(0.0, 1.0), new Normal(100.0, 5.0), new StudentT(1.0), new StudentT(30.0),
                new Uniform(0.0, 1.0), new Uniform(-3.0, 7.0), new Weibull(1.0, 0.5), new Weibull(1.0, 1.0),
                new Weibull(1.0, 2.0) };
    }

    private static String name(ContinuousDistribution d) {
        return d.getClass().getSimpleName() + "[" + d.supportLowerBound() + ", " + d.supportUpperBound() + "]";
    }

    /**
     * A density is zero beyond either end of the support. The point just past
     * the bound is the interesting one: it is as close to the mass as a
     * {@code double} can get without being in it.
     */
    @Test
    public void theDensityVanishesOutsideTheSupport() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            double lo = d.supportLowerBound();
            double hi = d.supportUpperBound();
            if (lo != NEG) {
                assertEquals(name(d) + ": the density does not vanish just below the support",
                        0.0, d.pdf(Math.nextAfter(lo, NEG)), 0.0);
                assertEquals(name(d) + ": the density does not vanish far below the support",
                        0.0, d.pdf(-MAX), 0.0);
            }
            if (hi != POS) {
                assertEquals(name(d) + ": the density does not vanish just above the support",
                        0.0, d.pdf(Math.nextAfter(hi, POS)), 0.0);
                assertEquals(name(d) + ": the density does not vanish far above the support",
                        0.0, d.pdf(MAX), 0.0);
            }
        }
    }

    /**
     * A density vanishes at both infinities, and at the top of the double range
     * where {@code rate * x} or {@code x^k} overflows on the way there. Before
     * the repair {@code Gamma} with a shape above one, {@code ChiSquare} above
     * two degrees of freedom, {@code Weibull} with a shape other than one and
     * {@code FisherF(1, 1)} all returned {@code NaN} from one of these.
     */
    @Test
    public void theDensityVanishesAtInfinity() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            assertEquals(name(d) + ": the density at minus infinity", 0.0, d.pdf(NEG), 0.0);
            assertEquals(name(d) + ": the density at plus infinity", 0.0, d.pdf(POS), 0.0);
            assertEquals(name(d) + ": the density at 1e308", 0.0, d.pdf(1.0e308), 0.0);
            assertEquals(name(d) + ": the density at MAX_VALUE", 0.0, d.pdf(MAX), 0.0);
        }
    }

    /**
     * The distribution function runs from zero to one and reaches both. The
     * lower end is asked for negligibility rather than for zero, because the
     * Cauchy tail is genuinely {@code 1.8e-309} at {@code -Double.MAX_VALUE}
     * and only the limit itself is zero.
     */
    @Test
    public void theDistributionFunctionRunsFromZeroToOne() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            assertEquals(name(d) + ": the distribution function at minus infinity", 0.0, d.cdf(NEG), 0.0);
            assertEquals(name(d) + ": the distribution function at plus infinity", 1.0, d.cdf(POS), 0.0);
            assertEquals(name(d) + ": the distribution function at 1e308", 1.0, d.cdf(1.0e308), 0.0);
            assertEquals(name(d) + ": the distribution function at MAX_VALUE", 1.0, d.cdf(MAX), 0.0);
            double bottom = d.cdf(-MAX);
            assertTrue(name(d) + ": the distribution function is " + bottom + " at -MAX_VALUE, not negligible",
                    bottom >= 0.0 && bottom <= 1.0e-300);
            double lo = d.supportLowerBound();
            double hi = d.supportUpperBound();
            if (lo != NEG) {
                assertEquals(name(d) + ": the distribution function at the lower bound", 0.0, d.cdf(lo), 0.0);
            }
            if (hi != POS) {
                assertEquals(name(d) + ": the distribution function at the upper bound", 1.0, d.cdf(hi), 0.0);
            }
        }
    }

    /**
     * A point that is not a number has neither a density nor a probability.
     * {@code Uniform} used to answer {@code 1 / (b - a)} here, because neither
     * {@code NaN < a} nor {@code NaN > b} holds and the guard let it through --
     * a plausible number for a question that was never asked.
     */
    @Test
    public void aQuestionThatIsNotANumberIsNotAnswered() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            double density = d.pdf(Double.NaN);
            double probability = d.cdf(Double.NaN);
            assertTrue(name(d) + ": the density at NaN is " + density, Double.isNaN(density));
            assertTrue(name(d) + ": the distribution function at NaN is " + probability, Double.isNaN(probability));
        }
    }

    /**
     * Nothing but {@code NaN} itself may produce {@code NaN}. This is the rule
     * the four density faults above and both distribution function faults break
     * at once, stated once for every point the other tests visit.
     */
    @Test
    public void onlyNotANumberProducesNotANumber() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            double lo = d.supportLowerBound();
            double hi = d.supportUpperBound();
            double[] xs = new double[] { NEG, -MAX, -1.0, 0.0, lo, hi, 1.0, 1.0e15, 1.0e308, MAX, POS,
                    Math.nextAfter(lo, NEG), Math.nextAfter(hi, POS) };
            for (int k = 0; k < xs.length; ++k) {
                assertTrue(name(d) + ": the density at " + xs[k] + " is not a number",
                        !Double.isNaN(d.pdf(xs[k])));
                assertTrue(name(d) + ": the distribution function at " + xs[k] + " is not a number",
                        !Double.isNaN(d.cdf(xs[k])));
            }
        }
    }

    /**
     * The two bounds are not decoration: they bracket the whole mass, they are
     * ordered, and the mean lies between them wherever one exists.
     */
    @Test
    public void theSupportBoundsBracketTheMass() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            double lo = d.supportLowerBound();
            double hi = d.supportUpperBound();
            assertTrue(name(d) + ": the bounds are not ordered", lo < hi);
            assertEquals(name(d) + ": the support does not hold all of the mass", 1.0, d.cdf(hi) - d.cdf(lo), 0.0);
            double mu = d.mean();
            if (!Double.isNaN(mu)) {
                assertTrue(name(d) + ": the mean " + mu + " lies outside the support", mu >= lo && mu <= hi);
            }
        }
    }

    /**
     * The quantile agrees with the bounds at both ends. {@code Gamma},
     * {@code ChiSquare}, {@code Exponential} and {@code Weibull} used to answer
     * {@code Double.MAX_VALUE} at a probability of one where {@code FisherF},
     * {@code LogNormal}, {@code Normal}, {@code StudentT} and {@code Cauchy}
     * answered infinity -- the same unbounded support, two different answers.
     */
    @Test
    public void theQuantileEndsWhereTheSupportEnds() {
        ContinuousDistribution[] ds = all();
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            assertEquals(name(d) + ": the quantile at zero", d.supportLowerBound(), d.inverseCdf(0.0), 0.0);
            assertEquals(name(d) + ": the quantile below zero", d.supportLowerBound(), d.inverseCdf(-1.0), 0.0);
            assertEquals(name(d) + ": the quantile at one", d.supportUpperBound(), d.inverseCdf(1.0), 0.0);
            assertEquals(name(d) + ": the quantile above one", d.supportUpperBound(), d.inverseCdf(2.0), 0.0);
        }
    }

    /**
     * None of the new guards has swallowed a real value. The density is finite
     * and non-negative at interior points, the distribution function stays in
     * [0, 1], and neither is flat at zero across the whole support -- which is
     * what a guard that returned too early would look like.
     */
    @Test
    public void theInteriorIsUndisturbed() {
        ContinuousDistribution[] ds = all();
        long lcg = 20260822L;
        for (int i = 0; i < ds.length; ++i) {
            ContinuousDistribution d = ds[i];
            double lo = d.supportLowerBound();
            double hi = d.supportUpperBound();
            // a window that is inside the support and holds the bulk
            double left = (lo == NEG) ? d.inverseCdf(0.01) : lo;
            double right = (hi == POS) ? d.inverseCdf(0.99) : hi;
            double reached = 0.0;
            double largest = 0.0;
            for (int k = 0; k <= 20; ++k) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double u = ((lcg >>> 11) * 0x1.0p-53) * 0.9 + 0.05;
                double x = left + (right - left) * ((k / 20.0) * 0.5 + u * 0.5);
                double density = d.pdf(x);
                assertTrue(name(d) + " at x = " + x + ": the density is " + density,
                        !Double.isNaN(density) && density >= 0.0);
                largest = Math.max(largest, density);
                double probability = d.cdf(x);
                assertTrue(name(d) + " at x = " + x + ": the distribution function is " + probability,
                        probability >= 0.0 && probability <= 1.0);
                reached = Math.max(reached, probability);
            }
            assertTrue(name(d) + ": the density is zero across its own support", largest > 0.0);
            assertTrue(name(d) + ": the distribution function never left zero", reached > 0.0);
        }
    }

    /**
     * A negative argument is answered rather than refused. {@code Gamma.pdf}
     * threw {@code IllegalArgumentException} below zero where its six siblings
     * with the same half line returned {@code 0.0}; the three shapes here are
     * the three branches the method takes at the origin.
     */
    @Test
    public void aNegativeArgumentNoLongerThrows() {
        double[] shapes = new double[] { 0.5, 1.0, 2.0 };
        for (int i = 0; i < shapes.length; ++i) {
            Gamma g = new Gamma(shapes[i], 1.0);
            assertEquals("Gamma(" + shapes[i] + "): the density at minus one", 0.0, g.pdf(-1.0), 0.0);
            assertEquals("Gamma(" + shapes[i] + "): the density at minus infinity", 0.0, g.pdf(NEG), 0.0);
            assertEquals("Gamma(" + shapes[i] + "): the density just below zero", 0.0,
                    g.pdf(Math.nextAfter(0.0, NEG)), 0.0);
        }
    }

    /**
     * The density at the lower bound follows the shape, not the arithmetic.
     * {@code Weibull} computed {@code x^(k-1) * x} at the origin, which reads
     * {@code infinity * zero} for a shape below one and gave {@code NaN} where
     * {@code Gamma} and {@code Beta} both report the pole they have.
     */
    @Test
    public void aPoleAtTheLowerBoundIsReportedAsOne() {
        assertEquals("Weibull with a shape below one has a pole at zero", POS, new Weibull(1.0, 0.5).pdf(0.0), 0.0);
        assertEquals("Weibull with a shape of one has its rate at zero", 1.0, new Weibull(1.0, 1.0).pdf(0.0), 0.0);
        assertEquals("Weibull with a shape above one vanishes at zero", 0.0, new Weibull(1.0, 2.0).pdf(0.0), 0.0);
        assertEquals("Gamma with a shape below one has a pole at zero", POS, new Gamma(0.5, 1.0).pdf(0.0), 0.0);
        assertEquals("Beta with both shapes below one has a pole at zero", POS, new Beta(0.5, 0.5).pdf(0.0), 0.0);
        assertEquals("Beta with both shapes below one has a pole at one", POS, new Beta(0.5, 0.5).pdf(1.0), 0.0);
    }

    /**
     * The overflow that produced the two named faults, at the argument that
     * produces it rather than at infinity. {@code Gamma(9, 0.5)} has a rate of
     * two, so {@code rate * x} leaves the double range at half of
     * {@code Double.MAX_VALUE}; {@code FisherF} with two or more numerator
     * degrees of freedom does the same with {@code d1 * x}.
     */
    @Test
    public void theProductInsideTheDistributionFunctionNoLongerOverflows() {
        assertEquals("Gamma(9, 0.5) at 1e308", 1.0, new Gamma(9.0, 0.5).cdf(1.0e308), 0.0);
        assertEquals("Gamma(2, 1e-9) at 1e308", 1.0, new Gamma(2.0, 1.0e-9).cdf(1.0e308), 0.0);
        assertEquals("FisherF(2, 2) at 1e308", 1.0, new FisherF(2, 2).cdf(1.0e308), 0.0);
        assertEquals("FisherF(10, 20) at 1e308", 1.0, new FisherF(10, 20).cdf(1.0e308), 0.0);
        // and the same product in the density
        assertEquals("Gamma(9, 0.5) density at 1e308", 0.0, new Gamma(9.0, 0.5).pdf(1.0e308), 0.0);
        assertEquals("FisherF(1, 1) density at 1e308", 0.0, new FisherF(1, 1).pdf(1.0e308), 0.0);
        assertEquals("Weibull(1, 2) density at 1e308", 0.0, new Weibull(1.0, 2.0).pdf(1.0e308), 0.0);
    }

    /**
     * The bounds themselves, class by class, so that an override that is
     * quietly dropped or given the wrong constant is caught here rather than
     * through the behavior that depends on it.
     */
    @Test
    public void theBoundsAreTheOnesTheDistributionsHave() {
        assertBounds(new Beta(2.0, 3.0), 0.0, 1.0);
        assertBounds(new Cauchy(0.0, 1.0), NEG, POS);
        assertBounds(new ChiSquare(3.0), 0.0, POS);
        assertBounds(new Exponential(1.0), 0.0, POS);
        assertBounds(new FisherF(3, 4), 0.0, POS);
        assertBounds(new Gamma(2.0, 1.0), 0.0, POS);
        assertBounds(new LogNormal(0.0, 1.0), 0.0, POS);
        assertBounds(new Normal(0.0, 1.0), NEG, POS);
        assertBounds(new StudentT(5.0), NEG, POS);
        assertBounds(new Uniform(-3.0, 7.0), -3.0, 7.0);
        assertBounds(new Weibull(1.0, 2.0), 0.0, POS);
    }

    private static void assertBounds(ContinuousDistribution d, double lo, double hi) {
        assertEquals(d.getClass().getSimpleName() + ": the lower bound", lo, d.supportLowerBound(), 0.0);
        assertEquals(d.getClass().getSimpleName() + ": the upper bound", hi, d.supportUpperBound(), 0.0);
    }
}
