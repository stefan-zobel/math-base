package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.fun.DFunction;
import math.solve.RootFinder;

/**
 * {@link StudentT#inverseCdf(double)} is the only caller of
 * {@link ContinuousDistribution#findRoot(double, double, double, double)} in
 * the library, and it used to start the iteration at {@link StudentT#mean()},
 * which is {@code NaN} for one degree of freedom or fewer because the mean does
 * not exist there. Every quantile below that threshold came back {@code NaN},
 * silently, although the distribution is the Cauchy and its quantile is
 * {@code tan(pi (p - 1/2))} -- a value the class could confirm through its own
 * {@code cdf}. The search now starts at zero, the median of every symmetric t.
 */
public class StudentTTest {

    private static final double[] PS = { 1.0e-6, 0.001, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.999,
            1.0 - 1.0e-6 };

    /**
     * An independent root of {@code cdf(x) = p}. The bracket widens until it
     * holds the root, because {@code brentDekker} documents that as the
     * caller's job and the tails here are heavy enough to need it: with half a
     * degree of freedom the tail falls off like {@code |x|^-0.5}, so the
     * millionth quantile sits near {@code -1e11}.
     */
    private static double invert(final StudentT t, final double p) {
        final DFunction equation = new DFunction() {
            @Override
            public double apply(double x) {
                return t.cdf(x) - p;
            }
        };
        double bound = 10.0;
        while (equation.apply(-bound) * equation.apply(bound) > 0.0 && bound < 1.0e300) {
            bound *= 100.0;
        }
        return RootFinder.brentDekker(-bound, bound, equation, 1.0e-14);
    }

    /** With one degree of freedom the quantile is the Cauchy quantile, in closed form. */
    @Test
    public void quantilesWithOneDegreeOfFreedomAreCauchy() {
        StudentT t = new StudentT(1.0);
        double[] ps = { 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975 };
        for (int k = 0; k < ps.length; ++k) {
            double expected = Math.tan(Math.PI * (ps[k] - 0.5));
            double actual = t.inverseCdf(ps[k]);
            assertTrue("p = " + ps[k] + " gave " + actual, !Double.isNaN(actual));
            assertEquals("p = " + ps[k], expected, actual, 1.0e-9 * Math.max(1.0, Math.abs(expected)));
        }
    }

    /** Below one degree of freedom the mean does not exist either, and it still has quantiles. */
    @Test
    public void quantilesExistWhereTheMeanDoesNot() {
        StudentT t = new StudentT(0.5);
        assertTrue("the mean is not supposed to exist here", Double.isNaN(t.mean()));

        double q = t.inverseCdf(0.75);
        assertTrue("inverseCdf returned " + q, !Double.isNaN(q));
        assertEquals("the quantile must satisfy its own equation", 0.75, t.cdf(q), 1.0e-9);
        assertEquals(invert(t, 0.75), q, 1.0e-8 * Math.max(1.0, Math.abs(q)));
    }

    /** The quantile against an independent root of the same equation, over the whole range. */
    @Test
    public void inverseCdfAgreesWithAnIndependentRoot() {
        double[] dfs = { 0.5, 1.0, 1.5, 2.0, 5.0, 99.0, 200.0 };
        for (int d = 0; d < dfs.length; ++d) {
            StudentT t = new StudentT(dfs[d]);
            for (int k = 0; k < PS.length; ++k) {
                double own = t.inverseCdf(PS[k]);
                double byRoot = invert(t, PS[k]);
                assertTrue("df " + dfs[d] + ", p " + PS[k] + " gave " + own, !Double.isNaN(own));
                assertEquals("df " + dfs[d] + ", p " + PS[k], byRoot, own,
                        1.0e-8 * Math.max(1.0, Math.abs(byRoot)));
            }
        }
    }

    /** Whatever the quantile is, feeding it back through the cdf has to return the probability. */
    @Test
    public void theQuantileSatisfiesItsOwnEquation() {
        double[] dfs = { 1.0, 2.0, 5.0, 99.0 };
        for (int d = 0; d < dfs.length; ++d) {
            StudentT t = new StudentT(dfs[d]);
            for (int k = 0; k < PS.length; ++k) {
                assertEquals("df " + dfs[d] + ", p " + PS[k], PS[k], t.cdf(t.inverseCdf(PS[k])), 1.0e-10);
            }
        }
    }

    /** The distribution is symmetric, so its quantiles are too. */
    @Test
    public void quantilesAreSymmetric() {
        double[] dfs = { 1.0, 3.0, 30.0 };
        for (int d = 0; d < dfs.length; ++d) {
            StudentT t = new StudentT(dfs[d]);
            double[] ps = { 0.01, 0.05, 0.25, 0.4 };
            for (int k = 0; k < ps.length; ++k) {
                double low = t.inverseCdf(ps[k]);
                double high = t.inverseCdf(1.0 - ps[k]);
                assertEquals("df " + dfs[d] + ", p " + ps[k], -low, high, 1.0e-9 * Math.max(1.0, Math.abs(high)));
            }
        }
    }

    /** The value every table of the t distribution prints. */
    @Test
    public void theTabulatedQuantileIsReproduced() {
        assertEquals(1.98422, new StudentT(99.0).inverseCdf(0.975), 1.0e-5);
        assertEquals(2.05954, new StudentT(25.0).inverseCdf(0.975), 1.0e-5);
        assertEquals(2.05553, new StudentT(26.0).inverseCdf(0.975), 1.0e-5);
        assertEquals(12.70620, new StudentT(1.0).inverseCdf(0.975), 1.0e-5);
    }

    /** As the degrees of freedom grow the t becomes the normal. */
    @Test
    public void largeDegreesOfFreedomApproachTheNormal() {
        double normalQuantile = new Normal().inverseCdf(0.975);
        double small = Math.abs(new StudentT(10.0).inverseCdf(0.975) - normalQuantile);
        double large = Math.abs(new StudentT(10000.0).inverseCdf(0.975) - normalQuantile);
        assertTrue("t(10) is " + small + " away, t(10000) is " + large, large < small);
        assertEquals(normalQuantile, new StudentT(1.0e7).inverseCdf(0.975), 1.0e-4);
    }
}
