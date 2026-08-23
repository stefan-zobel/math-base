package math.distribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.ClenshawCurtis;
import math.solve.Quadrature;
import math.solve.RootFinder;

/**
 * The eleven continuous distributions against the quadrature and root finding
 * of {@code math.solve}. Every distribution states a density, a distribution
 * function, a quantile and two moments, and those are bound together by
 * identities that nothing in this repository checked before: the two packages
 * had no relationship at all.
 * <p>
 * <b>Why this is not a distribution checking itself.</b> For seven of the
 * eleven the distribution function routes through the Colt ports in
 * {@code math.cern} -- {@code Normal}, {@code Beta}, {@code Gamma},
 * {@code LogNormal}, {@code StudentT}, and {@code ChiSquare} and
 * {@code FisherF} by delegation -- while the density is an elementary closed
 * form written in the class itself. Integrating the one and comparing it
 * against the other puts two different pieces of code on the two sides. The
 * other four -- {@code Weibull}, {@code Cauchy}, {@code Uniform},
 * {@code Exponential} -- are elementary on both sides and in the same file, so
 * there this checks the relationship between density and distribution function
 * rather than two independent implementations of it.
 * <p>
 * <b>The floors were measured, not chosen.</b> Each row carries what this
 * library achieved when the sweep was written, rounded down by roughly half a
 * digit to a digit. They guard against regression; they are not claims that the
 * library cannot do better.
 */
public class ContinuousDistributionTest {

    /** The rule is exact well before this; see the depth study in the report. */
    private static final int DEPTH = 22;

    /** Halved at every level, so it is the depth that governs, not this. */
    private static final double TOL = 1.0e-13;

    static final double[] PS = { 0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999 };

    /** The density is unbounded at an endpoint of the support. */
    static final int SINGULAR = 1;
    /** No mean exists, so no variance does either. */
    static final int NO_MEAN = 2;
    /** A mean exists but no finite variance does. */
    static final int NO_VARIANCE = 4;

    // -----------------------------------------------------------------
    // digits
    // -----------------------------------------------------------------

    /** Digits of absolute agreement, which is the scale a probability has. */
    static double absDigits(double got, double want) {
        if (got == want) {
            return 17.0;
        }
        return clip(-Math.log10(Math.abs(got - want)));
    }

    /** Digits of relative agreement, with an absolute floor of one. */
    static double relDigits(double got, double want) {
        if (got == want) {
            return 17.0;
        }
        return clip(-Math.log10(Math.abs(got - want) / Math.max(1.0, Math.abs(want))));
    }

    private static double clip(double d) {
        if (Double.isNaN(d)) {
            return 0.0;
        }
        return Math.max(0.0, Math.min(17.0, d));
    }

    // -----------------------------------------------------------------
    // the table
    // -----------------------------------------------------------------

    static final class Row {
        final String name;
        final String params;
        final ContinuousDistribution d;
        final double lo;
        final double hi;
        double mlo;
        double mhi;
        int flags;
        double tail;
        double aFloor;
        double cFloor;
        double dxFloor;
        double dpFloor;
        double muFloor;
        double varFloor;

        Row(String name, String params, ContinuousDistribution d, double lo, double hi) {
            this.name = name;
            this.params = params;
            this.d = d;
            this.lo = lo;
            this.hi = hi;
            this.mlo = lo;
            this.mhi = hi;
        }

        Row moments(double a, double b) {
            this.mlo = a;
            this.mhi = b;
            return this;
        }

        Row flags(int f) {
            this.flags = f;
            return this;
        }

        /** Mass the window is allowed to leave outside itself. */
        Row tail(double t) {
            this.tail = t;
            return this;
        }

        Row floors(double a, double c, double dx, double dp, double mu, double var) {
            this.aFloor = a;
            this.cFloor = c;
            this.dxFloor = dx;
            this.dpFloor = dp;
            this.muFloor = mu;
            this.varFloor = var;
            return this;
        }

        boolean is(int f) {
            return (flags & f) != 0;
        }

        boolean hasMean() {
            return !is(NO_MEAN);
        }

        boolean hasVariance() {
            return !is(NO_MEAN) && !is(NO_VARIANCE);
        }

        boolean onHalfLine() {
            return lo == 0.0;
        }

        @Override
        public String toString() {
            return name + "(" + params + ")";
        }
    }

    private static Row row(String name, String params, ContinuousDistribution d, double lo, double hi) {
        return new Row(name, params, d, lo, hi);
    }

    /**
     * The rows. Windows are data, deliberately: deriving them from
     * {@code inverseCdf} would make the sweep circular with one of the things
     * it is testing.
     */
    static List<Row> rows() {
        List<Row> r = new ArrayList<Row>();
        r.add(row("Normal", "0, 1", new Normal(0.0, 1.0), -10.0, 10.0).moments(-12.0, 12.0)
                .floors(15.0, 15.0, 13.0, 15.0, 15.0, 15.0));
        r.add(row("Normal", "3, 2", new Normal(3.0, 2.0), -17.0, 23.0).moments(-21.0, 27.0)
                .floors(15.0, 15.0, 13.0, 15.0, 14.0, 15.0));
        r.add(row("Normal", "299.85, 0.079", new Normal(299.8524, 0.079011), 298.9, 300.8).moments(298.75, 300.96)
                .floors(12.0, 13.0, 14.0, 12.0, 13.0, 14.0));
        r.add(row("Uniform", "0, 1", new Uniform(0.0, 1.0), 0.0, 1.0)
                .floors(16.0, 16.0, 13.0, 16.0, 16.0, 16.0));
        r.add(row("Uniform", "-3, 7", new Uniform(-3.0, 7.0), -3.0, 7.0)
                .floors(16.0, 15.0, 13.0, 16.0, 16.0, 15.0));
        r.add(row("Exponential", "lambda=1", new Exponential(1.0), 0.0, 40.0).moments(0.0, 60.0)
                .floors(16.0, 15.0, 13.0, 16.0, 16.0, 16.0));
        r.add(row("Exponential", "lambda=0.1", new Exponential(0.1), 0.0, 400.0).moments(0.0, 500.0)
                .floors(15.0, 15.0, 13.0, 16.0, 15.0, 15.0));
        r.add(row("Exponential", "lambda=7", new Exponential(7.0), 0.0, 6.0).moments(0.0, 9.0)
                .floors(15.0, 16.0, 13.0, 16.0, 16.0, 16.0));
        r.add(row("Cauchy", "0, 1", new Cauchy(0.0, 1.0), -1000.0, 1000.0).flags(NO_MEAN).tail(7.0e-4)
                .floors(16.0, 15.0, 13.0, 15.0, 0.0, 0.0));
        r.add(row("Cauchy", "2, 0.5", new Cauchy(2.0, 0.5), -498.0, 502.0).flags(NO_MEAN).tail(7.0e-4)
                .floors(16.0, 15.0, 13.0, 15.0, 0.0, 0.0));
        r.add(row("ChiSquare", "df=1", new ChiSquare(1.0), 0.0, 80.0).moments(0.0, 100.0).flags(SINGULAR)
                .floors(3.5, 4.0, 13.0, 15.0, 11.0, 4.0));
        r.add(row("ChiSquare", "df=2", new ChiSquare(2.0), 0.0, 80.0).moments(0.0, 100.0)
                .floors(16.0, 15.0, 13.0, 16.0, 16.0, 16.0));
        r.add(row("ChiSquare", "df=3", new ChiSquare(3.0), 0.0, 80.0).moments(0.0, 100.0)
                .floors(11.0, 13.0, 13.0, 15.0, 14.0, 11.0));
        r.add(row("ChiSquare", "df=30", new ChiSquare(30.0), 0.0, 150.0).moments(0.0, 150.0)
                .floors(14.0, 14.0, 14.0, 14.0, 14.0, 14.0));
        // the round-trip floor here is 15 and not 16 for the reason spelled
        // out at df=8 below: the quantile of p = 0.9 comes back exact or one
        // ulp away, nothing in between, and which of the two is decided by
        // the last bit of the density
        r.add(row("StudentT", "df=1", new StudentT(1.0), -1000.0, 1000.0).flags(NO_MEAN).tail(7.0e-4)
                .floors(14.0, 14.0, 12.0, 15.0, 0.0, 0.0));
        r.add(row("StudentT", "df=2", new StudentT(2.0), -500.0, 500.0).flags(NO_VARIANCE).tail(5.0e-6)
                .floors(15.0, 15.0, 13.0, 15.0, 15.0, 0.0));
        // mean floor 15, same reason as df=8 below
        r.add(row("StudentT", "df=3", new StudentT(3.0), -300.0, 300.0).moments(-3000.0, 3000.0).tail(1.0e-7)
                .floors(15.0, 15.0, 12.0, 15.0, 15.0, 3.0));
        // the mean floor here is 15 and not 16 because the mean is exactly
        // zero: relDigits then reports either 17, when the quadrature happens
        // to cancel to the last bit, or 15.9546, when it is left holding one
        // ulp of one. Which of the two a row draws is decided by the last bit
        // of the density, and it moved from df=100 to df=8 when StudentT
        // stopped forming its normalizing constant as a difference of two
        // logarithms -- with the density itself bit for bit unchanged there.
        // Normal(0, 1), df=2 and df=100 already carry 15 for the same reason
        r.add(row("StudentT", "df=8", new StudentT(8.0), -60.0, 60.0).moments(-600.0, 600.0).tail(1.0e-11)
                .floors(14.0, 14.0, 13.0, 15.0, 15.0, 13.0));
        r.add(row("StudentT", "df=100", new StudentT(100.0), -30.0, 30.0).moments(-60.0, 60.0)
                .floors(13.0, 13.0, 13.0, 14.0, 15.0, 13.0));
        r.add(row("Beta", "2, 5", new Beta(2.0, 5.0), 0.0, 1.0)
                .floors(14.0, 14.0, 13.0, 15.0, 14.0, 15.0));
        r.add(row("Beta", "2, 2", new Beta(2.0, 2.0), 0.0, 1.0)
                .floors(14.0, 15.0, 13.0, 16.0, 15.0, 16.0));
        r.add(row("Beta", "0.5, 0.5", new Beta(0.5, 0.5), 0.0, 1.0).flags(SINGULAR)
                .floors(4.0, 5.0, 13.0, 13.0, 4.5, 5.0));
        r.add(row("Beta", "5, 1", new Beta(5.0, 1.0), 0.0, 1.0)
                .floors(14.0, 14.0, 13.0, 15.0, 14.0, 16.0));
        r.add(row("Beta", "0.3, 4", new Beta(0.3, 4.0), 0.0, 1.0).flags(SINGULAR)
                .floors(2.0, 2.5, 13.0, 13.0, 12.0, 4.5));
        r.add(row("Gamma", "k=0.5, theta=1", new Gamma(0.5, 1.0), 0.0, 80.0).moments(0.0, 100.0).flags(SINGULAR)
                .floors(3.5, 4.0, 13.0, 15.0, 11.0, 4.0));
        r.add(row("Gamma", "k=1, theta=1", new Gamma(1.0, 1.0), 0.0, 45.0).moments(0.0, 60.0)
                .floors(16.0, 15.0, 13.0, 16.0, 16.0, 16.0));
        r.add(row("Gamma", "k=2, theta=1", new Gamma(2.0, 1.0), 0.0, 50.0).moments(0.0, 60.0)
                .floors(14.0, 14.0, 13.0, 15.0, 14.0, 14.0));
        r.add(row("Gamma", "k=9, theta=0.5", new Gamma(9.0, 0.5), 0.0, 30.0).moments(0.0, 40.0)
                .floors(14.0, 14.0, 13.0, 15.0, 14.0, 14.0));
        r.add(row("LogNormal", "0, 1", new LogNormal(0.0, 1.0), 0.0, 1000.0).moments(0.0, 1.0e4).tail(3.0e-12)
                .floors(15.0, 15.0, 13.0, 15.0, 15.0, 12.0));
        r.add(row("LogNormal", "0, 2", new LogNormal(0.0, 2.0), 0.0, 1.0e5).moments(0.0, 1.0e7).tail(5.0e-9)
                .floors(6.0, 9.0, 13.0, 15.0, 5.5, 4.5));
        r.add(row("LogNormal", "3, 0.25", new LogNormal(3.0, 0.25), 0.0, 200.0).moments(0.0, 400.0)
                .floors(15.0, 15.0, 14.0, 15.0, 16.0, 14.0));
        r.add(row("FisherF", "1, 1", new FisherF(1, 1), 0.0, 1.0e6).flags(SINGULAR | NO_MEAN).tail(7.0e-4)
                .floors(2.0, 2.5, 13.0, 14.0, 0.0, 0.0));
        r.add(row("FisherF", "5, 2", new FisherF(5, 2), 0.0, 1.0e5).flags(NO_MEAN).tail(2.0e-5)
                .floors(10.0, 12.0, 13.0, 15.0, 0.0, 0.0));
        r.add(row("FisherF", "10, 20", new FisherF(10, 20), 0.0, 500.0).moments(0.0, 2000.0).tail(1.0e-15)
                .floors(14.0, 14.0, 13.0, 15.0, 14.0, 14.0));
        r.add(row("Weibull", "lambda=1, k=1", new Weibull(1.0, 1.0), 0.0, 45.0).moments(0.0, 60.0)
                .floors(16.0, 15.0, 13.0, 16.0, 14.0, 15.0));
        r.add(row("Weibull", "lambda=1, k=0.5", new Weibull(1.0, 0.5), 0.0, 900.0).moments(0.0, 1200.0)
                .flags(SINGULAR).tail(1.0e-13).floors(3.0, 3.5, 13.0, 16.0, 10.0, 4.0));
        r.add(row("Weibull", "lambda=2, k=5", new Weibull(2.0, 5.0), 0.0, 12.0).moments(0.0, 16.0)
                .floors(16.0, 16.0, 13.0, 15.0, 14.0, 13.0));
        r.add(row("Weibull", "lambda=0.3, k=3", new Weibull(0.3, 3.0), 0.0, 3.0).moments(0.0, 5.0)
                .floors(15.0, 15.0, 13.0, 16.0, 14.0, 15.0));
        return r;
    }

    // -----------------------------------------------------------------
    // quadrature helpers
    // -----------------------------------------------------------------

    static DFunction density(final ContinuousDistribution d) {
        return new DFunction() {
            @Override
            public double apply(double x) {
                return d.pdf(x);
            }
        };
    }

    static double integral(ContinuousDistribution d, double a, double b) {
        return AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, density(d), a, b, TOL, DEPTH);
    }

    static double moment(final ContinuousDistribution d, final int power, final double centre, double a,
            double b) {
        return AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                double w = (power == 1) ? x : (x - centre) * (x - centre);
                return w * d.pdf(x);
            }
        }, a, b, TOL, DEPTH);
    }

    /**
     * A bracket that provably contains the root of {@code cdf(x) - p}, widened
     * outwards from the window until the distribution function straddles
     * {@code p}. {@code RootFinder.brentDekker} states that the interval must
     * contain the root and does not check it: handed one that does not, it
     * returns a confident wrong answer.
     */
    static double[] bracket(ContinuousDistribution d, double lo, double hi, double p) {
        double a = lo;
        double width = hi - lo;
        for (int i = 0; i < 200 && d.cdf(a) >= p; ++i) {
            a -= width;
            width += width;
        }
        double b = hi;
        width = hi - lo;
        for (int i = 0; i < 200 && d.cdf(b) <= p; ++i) {
            b += width;
            width += width;
        }
        return new double[] { a, b };
    }

    // -----------------------------------------------------------------
    // the identities
    // -----------------------------------------------------------------

    /** The table has to describe what it claims to describe. */
    @Test
    public void theTableIsIntact() {
        List<Row> rows = rows();
        assertEquals(38, rows.size());
        int distinct = 0;
        String previous = "";
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            assertTrue(row + ": empty window", row.hi > row.lo);
            assertTrue(row + ": moment window does not contain the window", row.mlo <= row.lo && row.mhi >= row.hi);
            assertTrue(row + ": tail mass", row.tail >= 0.0 && row.tail < 1.0e-3);
            assertTrue(row + ": no floor set", row.aFloor > 0.0);
            assertTrue(row + ": the density must be finite inside the window",
                    !Double.isNaN(row.d.pdf(0.5 * (row.lo + row.hi))));
            if (!row.name.equals(previous)) {
                distinct++;
                previous = row.name;
            }
        }
        assertEquals("every continuous distribution in the package is covered", 11, distinct);
    }

    /**
     * <b>A.</b> The integral of the density over the window is the difference
     * of the distribution function across it. This holds on any window
     * whatsoever, so unlike the others it cannot be tuned into passing by
     * choosing where to look.
     */
    @Test
    public void theIntegralOfTheDensityIsTheDistributionFunction() {
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            double got = integral(row.d, row.lo, row.hi);
            double want = row.d.cdf(row.hi) - row.d.cdf(row.lo);
            double digits = absDigits(got, want);
            assertTrue(row + ": density and distribution function agree to only " + digits + " digits",
                    digits >= row.aFloor);
        }
    }

    /**
     * <b>B.</b> The window holds the mass it claims to, and the distribution
     * function runs from zero to one over the whole support.
     */
    @Test
    public void theWindowHoldsTheMassItClaims() {
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            double mass = row.d.cdf(row.hi) - row.d.cdf(row.lo);
            assertTrue(row + ": the window holds " + mass + ", short by more than the declared tail",
                    mass >= 1.0 - row.tail);
            assertTrue(row + ": the window holds more than all of the mass", mass <= 1.0 + 1.0e-12);
            assertEquals(row + ": the distribution function must vanish at the bottom of the range", 0.0,
                    row.d.cdf(-Double.MAX_VALUE), 1.0e-300);
            assertEquals(row + ": the distribution function must reach one at the top of the range", 1.0,
                    row.d.cdf(Double.MAX_VALUE), 0.0);

            double previous = -1.0;
            for (int k = 0; k <= 40; ++k) {
                double x = row.lo + (row.hi - row.lo) * k / 40.0;
                double c = row.d.cdf(x);
                assertTrue(row + " at x = " + x + ": the distribution function is not a number",
                        !Double.isNaN(c));
                assertTrue(row + " at x = " + x + ": the distribution function decreased", c >= previous);
                previous = c;
            }
        }
    }

    /**
     * The usable range of the distribution function, which is now the whole
     * one. Every row is correct and monotone from {@code 1e3} up to
     * {@code Double.POSITIVE_INFINITY} inclusive.
     * <p>
     * This test used to stop at {@code 1e15} on purpose, because {@code Gamma}
     * and {@code FisherF} returned {@code NaN} rather than one at infinity and
     * within a factor of two of {@code Double.MAX_VALUE} -- the rate or the
     * numerator degrees of freedom made the product inside them overflow before
     * the argument itself did. {@code ContinuousDistributionContractTest} states
     * that rule over every implementation; what is left here is that the table
     * of this file obeys it too, monotonically, all the way out.
     */
    @Test
    public void theDistributionFunctionIsUsableAcrossTheWholeOrdinaryRange() {
        double[] xs = new double[] { 1.0e3, 1.0e15, 1.0e308, Double.MAX_VALUE, Double.POSITIVE_INFINITY };
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            double previous = -1.0;
            for (int k = 0; k < xs.length; ++k) {
                double c = row.d.cdf(xs[k]);
                assertTrue(row + ": the distribution function is not a number at " + xs[k], !Double.isNaN(c));
                assertTrue(row + ": the distribution function decreased below " + xs[k], c >= previous);
                assertTrue(row + ": the distribution function exceeds one at " + xs[k], c <= 1.0);
                previous = c;
            }
            assertEquals(row + ": the distribution function must reach one at infinity", 1.0,
                    row.d.cdf(Double.POSITIVE_INFINITY), 0.0);
        }
    }

    /**
     * <b>C.</b> The distribution function is the running integral of its own
     * density, panel by panel across the window, and not merely in total.
     */
    @Test
    public void theDistributionFunctionIsTheRunningIntegral() {
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            double running = row.d.cdf(row.lo);
            double previous = row.lo;
            for (int k = 1; k <= 9; ++k) {
                double x = row.lo + (row.hi - row.lo) * k / 10.0;
                running += integral(row.d, previous, x);
                previous = x;
                double digits = absDigits(running, row.d.cdf(x));
                assertTrue(row + " at x = " + x + ": running integral and cdf agree to only " + digits + " digits",
                        digits >= row.cFloor);
            }
        }
    }

    /**
     * <b>D.</b> The quantile against an independent root of the same equation.
     * {@code RootFinder.brentDekker} is not the Newton and bisection hybrid of
     * {@code ContinuousDistribution.findRoot} that most of these delegate to,
     * so this is two root finders and not one twice.
     * <p>
     * The position is checked with a floor of one and the probability
     * absolutely, because where the density is flat many positions carry the
     * same probability: there the round trip is the meaningful statement and
     * the position is not.
     */
    @Test
    public void theQuantileAgreesWithAnIndependentRoot() {
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            final Row row = rows.get(i);
            for (int k = 0; k < PS.length; ++k) {
                final double p = PS[k];
                double[] br = bracket(row.d, row.lo, row.hi, p);
                double byRoot = RootFinder.brentDekker(br[0], br[1], new DFunction() {
                    @Override
                    public double apply(double x) {
                        return row.d.cdf(x) - p;
                    }
                }, 1.0e-14);
                double inverse = row.d.inverseCdf(p);

                assertTrue(row + " at p = " + p + ": the quantile is not a number", !Double.isNaN(inverse));
                double dx = relDigits(inverse, byRoot);
                assertTrue(row + " at p = " + p + ": the two roots agree to only " + dx + " digits",
                        dx >= row.dxFloor);
                double dp = absDigits(row.d.cdf(inverse), p);
                assertTrue(row + " at p = " + p + ": the round trip closes to only " + dp + " digits",
                        dp >= row.dpFloor);
            }
        }
    }

    /**
     * <b>E.</b> The mean and the variance are the first two moments of the
     * density, integrated. Where a moment does not exist the distribution has
     * to say so rather than return a number, which is the harder half of the
     * claim.
     */
    @Test
    public void theMomentsAreTheIntegralsTheyClaimToBe() {
        List<Row> rows = rows();
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            if (!row.hasMean()) {
                assertTrue(row + ": there is no mean, so it must not report one", Double.isNaN(row.d.mean()));
                assertTrue(row + ": there is no variance either", Double.isNaN(row.d.variance()));
                continue;
            }
            double mu = row.d.mean();
            assertTrue(row + ": the mean is not a number", !Double.isNaN(mu));
            double first = moment(row.d, 1, 0.0, row.mlo, row.mhi);
            double digits = relDigits(first, mu);
            assertTrue(row + ": mean and its integral agree to only " + digits + " digits", digits >= row.muFloor);

            if (!row.hasVariance()) {
                double v = row.d.variance();
                assertTrue(row + ": there is no finite variance, so it must not report one",
                        Double.isNaN(v) || Double.isInfinite(v));
                continue;
            }
            double second = moment(row.d, 2, mu, row.mlo, row.mhi);
            double vDigits = relDigits(second, row.d.variance());
            assertTrue(row + ": variance and its integral agree to only " + vDigits + " digits",
                    vDigits >= row.varFloor);
        }
    }

    // -----------------------------------------------------------------
    // the awkward half
    // -----------------------------------------------------------------

    /**
     * A density that is unbounded at an endpoint is where the bisecting rule is
     * the wrong tool: the panel next to the singularity never becomes small
     * enough, and every extra level of recursion buys a fixed fraction of a
     * digit at twice the work. Measured over depths 8 to 30 the gain is about
     * 0.15 digits per level and it does not converge.
     * <p>
     * The point of the test is that this is a property of the quadrature and
     * not of the distributions: the same densities are reproduced to full
     * precision once the singular endpoint is excluded.
     */
    @Test
    public void singularEndpointsCostDigitsThatAnInnerWindowGivesBack() {
        List<Row> rows = rows();
        int singular = 0;
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            if (!row.is(SINGULAR)) {
                continue;
            }
            singular++;
            double outer = absDigits(integral(row.d, row.lo, row.hi),
                    row.d.cdf(row.hi) - row.d.cdf(row.lo));
            assertTrue(row + ": the singular endpoint should cost digits, but " + outer + " were reached",
                    outer < 8.0);

            // Beta(0.5, 0.5) is unbounded at both ends, so stepping in at one
            // of them recovers nothing: 5.1 digits against 15.0 for both.
            double step = 1.0e-6 * (row.hi - row.lo);
            double a = row.lo + step;
            double b = row.hi - step;
            double got = integral(row.d, a, b);
            double want = row.d.cdf(b) - row.d.cdf(a);
            double digits = absDigits(got, want);
            double surrendered = 1.0 - want;
            assertTrue(row + ": stepping in past " + surrendered + " of the mass still reaches only " + digits
                    + " digits", digits >= 13.0);
            assertTrue(row + ": the inner window must recover what the singular one lost", digits > outer + 3.0);
        }
        assertEquals("six rows carry an unbounded density", 6, singular);
    }

    /**
     * Cauchy and Student t with one or two degrees of freedom have no moments,
     * and the integrals that would define them have to diverge rather than
     * settle on a number. They diverge logarithmically, which is a sharper
     * statement than growth: over a symmetric window of half width {@code w}
     * the absolute first moment of the Cauchy grows by
     * {@code (2/pi) * ln(10)} per decade of {@code w}, and the second moment of
     * a Student t on two degrees of freedom by {@code 2 * ln(10)}.
     */
    @Test
    public void theMomentsThatDoNotExistDivergeAtTheRateTheyShould() {
        double[] widths = { 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6 };
        double[] cauchy = new double[widths.length];
        double[] student = new double[widths.length];
        Cauchy c = new Cauchy(0.0, 1.0);
        StudentT t2 = new StudentT(2.0);
        for (int i = 0; i < widths.length; ++i) {
            cauchy[i] = absoluteMoment(c, 1, widths[i]);
            student[i] = absoluteMoment(t2, 2, widths[i]);
        }

        double perDecadeCauchy = 2.0 * Math.log(10.0) / Math.PI;
        double perDecadeStudent = 2.0 * Math.log(10.0);
        for (int i = 1; i < widths.length; ++i) {
            assertEquals("Cauchy: the absolute first moment must grow by (2/pi) ln 10 per decade", perDecadeCauchy,
                    cauchy[i] - cauchy[i - 1], 1.0e-3);
            assertEquals("StudentT(2): the second moment must grow by 2 ln 10 per decade", perDecadeStudent,
                    student[i] - student[i - 1], 1.0e-3);
        }
        assertTrue("the divergence has to be visible over four decades", cauchy[4] > 2.9 * cauchy[0] / 1.0);
    }

    private static double absoluteMoment(final ContinuousDistribution d, final int power, double w) {
        return AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                double a = Math.abs(x);
                return (power == 1 ? a : a * a) * d.pdf(x);
            }
        }, -w, w, TOL, DEPTH);
    }

    // -----------------------------------------------------------------
    // the integrators against each other
    // -----------------------------------------------------------------

    /**
     * The one-sided infinite route against the complement of the distribution
     * function. This is that method's first exercise outside its own unit test,
     * and the half line is where it is sound: the doubly infinite form of the
     * same substitution is not, and has its own test below.
     */
    @Test
    public void theOneSidedInfiniteRouteReproducesTheComplement() {
        List<Row> rows = rows();
        double[] cuts = { 1.0, 5.0 };
        int checked = 0;
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            if (!row.onHalfLine() || row.is(SINGULAR) || row.hi <= 5.0) {
                continue;
            }
            for (int k = 0; k < cuts.length; ++k) {
                double got = Quadrature.integrate(density(row.d), cuts[k],
                        Double.POSITIVE_INFINITY, TOL);
                double want = 1.0 - row.d.cdf(cuts[k]);
                double digits = absDigits(got, want);
                assertTrue(row + " from " + cuts[k] + ": the infinite route reaches only " + digits + " digits",
                        digits >= 14.0);
                checked++;
            }
        }
        assertTrue("the half line should carry most of the table: " + checked, checked >= 20);
    }

    /**
     * The doubly infinite substitution {@code x = t / (1 - t^2)} maps the line
     * onto {@code [-1, 1]} and pushes mass far from the origin into a spike
     * against the endpoints. Whether the subdivision finds that spike turns on
     * how heavy the tail is and on how far the mass sits measured in its own
     * width, not on the distance alone: a Cauchy centered at a thousand is
     * integrated exactly, a normal of unit width at a hundred is not.
     * <p>
     * What it must never do is return the near-zero quietly. It no longer does:
     * the probe names the place it should have looked, the line is split there,
     * and the double-exponential rule integrates both halves. The refusal is
     * kept for the cases that rule cannot resolve either.
     */
    @Test
    public void theDoublyInfiniteRouteHoldsForHeavyTailsAndNowRecoversThinOnesFarOut() {
        assertEquals("a Cauchy is integrated wherever it sits", 1.0, wholeLine(new Cauchy(1000.0, 1.0)), 1.0e-8);
        assertEquals("so is a Student t", 1.0, wholeLine(new StudentT(3.0)), 1.0e-9);
        assertEquals("a normal at the origin is fine", 1.0, wholeLine(new Normal(0.0, 1.0)), 1.0e-9);
        assertEquals("and still fine at thirty", 1.0, wholeLine(new Normal(30.0, 1.0)), 1.0e-9);
        assertEquals("scale delays it rather than preventing it", 1.0, wholeLine(new Normal(100.0, 10.0)), 1.0e-9);

        // These four used to be refused. The algebraic substitution misses them,
        // and the refusal named the point where the mass sits and told the
        // caller to split the range there. That split is now made for the
        // caller, both halves going to the double-exponential rule, whose nodes
        // crowd towards the finite end of each half - which is the peak itself.
        recovered(new Normal(100.0, 1.0));
        recovered(new Normal(137.0, 1.0));
        recovered(new Normal(300.0, 10.0));
        recovered(new Normal(-300.0, 1.0));
    }

    @Test
    public void theRefusalSurvivesWhereTheFallbackCannotAnswerEither() {
        // The fallback is trusted only when it reports convergence, so an
        // integrand it cannot resolve leaves the refusal standing. A fast
        // oscillation riding on mass that sits far out is such a case: the
        // probe finds the mass, the substitution never sampled it, and no
        // double-exponential node set resolves the oscillation.
        final DFunction bell = density(new Normal(200.0, 1.0));
        DFunction f = x -> bell.apply(x) * (1.0 + 0.9 * Math.sin(1.0e4 * x));
        try {
            wholeLine(f);
            org.junit.Assert.fail("an unresolvable oscillation far from the origin must still be refused");
        } catch (ArithmeticException expected) {
            assertTrue("the refusal must still say what went wrong: " + expected.getMessage(),
                    expected.getMessage().contains("never sampled the integrand where its mass lies"));
        }
    }

    /**
     * A density the substitution alone misses, which has to come back as one
     * anyway. The remedy the refusal used to propose - split the line where the
     * mass is - is now applied for the caller, so what was the message is now
     * the answer.
     */
    private static void recovered(ContinuousDistribution d) {
        assertEquals("the substitution misses " + d + " and the fallback has to catch it", 1.0, wholeLine(d),
                1.0e-9);
    }

    /**
     * The other half of a check like this, and the half that decides whether it
     * is worth having: it must not cry wolf. An integral that is genuinely zero
     * by cancellation, an integrand that is zero everywhere, one that is merely
     * very small, and heavy tails all have to pass untouched.
     */
    @Test
    public void theRefusalDoesNotFireOnIntegralsThatAreHonestlySmall() {
        final Normal standard = new Normal(0.0, 1.0);
        assertEquals("an odd integrand cancels to zero and is sampled where it is large", 0.0,
                wholeLine(new DFunction() {
                    @Override
                    public double apply(double x) {
                        return x * standard.pdf(x);
                    }
                }), 1.0e-12);
        assertEquals("the zero function integrates to zero", 0.0, wholeLine(new DFunction() {
            @Override
            public double apply(double x) {
                return 0.0;
            }
        }), 0.0);
        assertEquals("a genuinely tiny integral is not a missed one", 1.0e-200, wholeLine(new DFunction() {
            @Override
            public double apply(double x) {
                return 1.0e-200 * standard.pdf(x);
            }
        }), 1.0e-208);
        assertEquals("heavy tails are fine", Math.PI, wholeLine(new DFunction() {
            @Override
            public double apply(double x) {
                return 1.0 / (1.0 + x * x);
            }
        }), 1.0e-9);
        assertEquals("and so is a slowly decaying one", 2.0, wholeLine(new DFunction() {
            @Override
            public double apply(double x) {
                return Math.exp(-Math.abs(x));
            }
        }), 1.0e-9);
    }

    /**
     * The one-sided forms share the fault and share the remedy: mass far from
     * the finite end is missed the same way, and recovered the same way.
     */
    @Test
    public void theOneSidedFormsRecoverTheSameWay() {
        Normal far = new Normal(300.0, 1.0);
        assertEquals("the one-sided substitution misses the mass and the fallback catches it", 1.0,
                Quadrature.integrate(density(far), 0.0,
                        Double.POSITIVE_INFINITY, TOL),
                1.0e-9);
        assertEquals("and cutting close to the mass by hand is all it ever took", 1.0,
                Quadrature.integrate(density(far), 290.0,
                        Double.POSITIVE_INFINITY, TOL),
                1.0e-9);
    }

    private static double wholeLine(ContinuousDistribution d) {
        return wholeLine(density(d));
    }

    private static double wholeLine(DFunction f) {
        return Quadrature.integrate(f, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, TOL);
    }

    /**
     * Clenshaw-Curtis clusters its nodes at the ends of the interval and the
     * bisecting rule does not, so the two are genuinely different quadratures.
     * Where Clenshaw-Curtis reports convergence it has to agree with adaptive
     * Gauss-Kronrod, and where it cannot it has to say so: its grid is
     * Chebyshev-Lobatto and therefore includes the endpoints, so an unbounded
     * density is evaluated at the singularity itself. It reports that honestly
     * rather than returning a value it does not have.
     */
    @Test
    public void clenshawCurtisAgreesWhereverItClaimsToHaveConverged() {
        List<Row> rows = rows();
        int converged = 0;
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            ClenshawCurtis.IntegralResult cc = ClenshawCurtis.integrate1D(density(row.d), row.lo, row.hi, TOL);
            if (row.is(SINGULAR)) {
                assertTrue(row + ": an unbounded endpoint must not be reported as converged", !cc.converged);
                continue;
            }
            if (!cc.converged) {
                continue;
            }
            converged++;
            double agk = integral(row.d, row.lo, row.hi);
            assertEquals(row + ": the two quadratures disagree although Clenshaw-Curtis reports convergence", agk,
                    cc.value, 1.0e-11);
        }
        assertTrue("Clenshaw-Curtis should converge on most of the smooth rows: " + converged, converged >= 20);
    }
}
