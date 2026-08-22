package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.concurrent.atomic.AtomicLong;

import org.junit.Test;

import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.DoubleExponential.IntegralResult;

/**
 * Tests for the double-exponential rule. What is asserted here are the
 * properties the rule is supposed to have - machine precision at an integrable
 * endpoint singularity, a node set that doubles, and above all a convergence
 * flag that can be believed - rather than golden values.
 */
public class DoubleExponentialTest {

    private static final double TOL = 1.0e-13;
    private static final double INF = Double.POSITIVE_INFINITY;
    private static final double N_INF = Double.NEGATIVE_INFINITY;

    private static double relative(double got, double want) {
        return Math.abs(got - want) / Math.abs(want);
    }

    @Test
    public void testEndpointSingularitiesReachMachinePrecision() {
        // the case the rule exists for: an integrable pole at a limit, where
        // adaptive bisection gains about a seventh of a digit per level
        DFunction[] f = { x -> 1.0 / Math.sqrt(x), x -> Math.log(x), x -> Math.pow(x, -0.9),
                x -> Math.sqrt(x) * Math.log(1.0 / x) };
        double[] exact = { 2.0, -1.0, 10.0, 4.0 / 9.0 };
        String[] name = { "1/sqrt(x)", "log(x)", "x^-0.9", "sqrt(x)*log(1/x)" };
        for (int i = 0; i < f.length; i++) {
            IntegralResult r = DoubleExponential.integrate1D(f[i], 0.0, 1.0, TOL);
            assertTrue(name[i] + " must converge", r.converged);
            assertTrue(name[i] + " reached " + relative(r.value, exact[i]),
                    relative(r.value, exact[i]) <= 1.0e-14);
            assertTrue(name[i] + " took " + r.evaluations + " evaluations", r.evaluations <= 400L);
        }

        // a density with a pole at the origin, over a window that is 900 wide
        IntegralResult w = DoubleExponential.integrate1D(x -> 0.5 * Math.pow(x, -0.5) * Math.exp(-Math.sqrt(x)),
                0.0, 900.0, TOL);
        assertTrue("the Weibull density must converge", w.converged);
        assertEquals(1.0 - Math.exp(-30.0), w.value, 1.0e-14);
    }

    @Test
    public void testTheHalfLineAndTheWholeLine() {
        // a half line and a pole at the finite end at the same time, which the
        // algebraic substitutions have no route for
        IntegralResult h = DoubleExponential.integrate1D(x -> Math.exp(-x) / Math.sqrt(x), 0.0, INF, TOL);
        assertTrue("exp(-x)/sqrt(x) must converge", h.converged);
        assertTrue(relative(h.value, Math.sqrt(Math.PI)) <= 1.0e-14);

        IntegralResult heavy = DoubleExponential.integrate1D(x -> 1.0 / (1.0 + x * x), 0.0, INF, TOL);
        assertTrue("a heavy tail must converge", heavy.converged);
        assertTrue(relative(heavy.value, 0.5 * Math.PI) <= 1.0e-14);

        IntegralResult g = DoubleExponential.integrate1D(x -> Math.exp(-x * x), N_INF, INF, TOL);
        assertTrue("exp(-x*x) over the line must converge", g.converged);
        assertTrue(relative(g.value, Math.sqrt(Math.PI)) <= 1.0e-14);

        IntegralResult c = DoubleExponential.integrate1D(x -> 1.0 / (Math.PI * (1.0 + x * x)), N_INF, INF, TOL);
        assertTrue("the Cauchy density over the line must converge", c.converged);
        assertTrue(relative(c.value, 1.0) <= 1.0e-14);
    }

    @Test
    public void testTheReflectedFormIsBitIdenticalToTheUpperOne() {
        // (-inf, b] is [-b, inf) of the reflected integrand and nothing else,
        // so the two must not differ in the last place either
        DFunction f = x -> Math.exp(-x * x) * (1.0 + 0.5 * x);
        IntegralResult up = DoubleExponential.integrate1D(f, 2.0, INF, TOL);
        IntegralResult down = DoubleExponential.integrate1D(x -> f.apply(-x), N_INF, -2.0, TOL);
        assertTrue(up.converged);
        assertEquals("the reflection must agree bit for bit", up.value, down.value, 0.0);
        assertEquals(up.evaluations, down.evaluations);
    }

    @Test
    public void testAnIntegrandThatIsZeroWhereverItLooksIsNotCertified() {
        // the rule cannot tell a function that is zero everywhere from one that
        // is zero everywhere it looked, and must not pretend otherwise
        IntegralResult r = DoubleExponential.integrate1D(x -> 0.0, 0.0, 1.0, TOL);
        assertFalse("a zero integrand must not be certified", r.converged);
        assertEquals(0.0, r.value, 0.0);
    }

    @Test
    public void testAnIntegralThatIsGenuinelyZeroOrGenuinelyTinyStillConverges() {
        // the tolerance is scaled by the integral of |f|, not by the result, so
        // neither of these is confused with having missed the integrand
        IntegralResult zero = DoubleExponential.integrate1D(x -> x - 0.5, 0.0, 1.0, TOL);
        assertTrue("an integral that is genuinely zero must converge", zero.converged);
        assertEquals(0.0, zero.value, 1.0e-15);

        IntegralResult tiny = DoubleExponential.integrate1D(x -> 1.0e-200 * x, 0.0, 1.0, TOL);
        assertTrue("an integral that is genuinely 1e-200 must converge", tiny.converged);
        assertTrue(relative(tiny.value, 0.5e-200) <= 1.0e-14);
    }

    @Test
    public void testMassTooFarOutIsRefusedRatherThanReturned() {
        // when a single node lands on the mass, halving the step halves the
        // value: 9.8e-154, 4.9e-154, 2.4e-154 is a miss, not a convergence
        IntegralResult far = DoubleExponential.integrate1D(normal(3000.0), N_INF, INF, TOL);
        assertFalse("mass at 3000 widths out must not be certified", far.converged);

        IntegralResult further = DoubleExponential.integrate1D(normal(100000.0), N_INF, INF, TOL);
        assertFalse("mass at 100000 widths out must not be certified", further.converged);

        // and the same integrand centered where the rule can see it is fine
        IntegralResult here = DoubleExponential.integrate1D(normal(0.0), N_INF, INF, TOL);
        assertTrue("mass at the origin must converge", here.converged);
        assertTrue(relative(here.value, 1.0) <= 1.0e-14);
    }

    private static DFunction normal(final double mu) {
        return x -> {
            double z = x - mu;
            return Math.exp(-0.5 * z * z) / Math.sqrt(2.0 * Math.PI);
        };
    }

    @Test
    public void testHalvingTheStepDoublesTheNodeSet() {
        // the nodes are t = j*h, so the previous level is the even j: one level
        // costs as much as everything before it put together, which is what
        // makes a generous ceiling free
        double[][] limits = { { 0.0, 1.0 }, { 0.0, INF }, { N_INF, INF } };
        String[] member = { "tanh-sinh", "exp-sinh", "sinh-sinh" };
        DFunction f = x -> Math.exp(-x * x);
        for (int m = 0; m < limits.length; m++) {
            long previous = 0L;
            for (int level = 2; level <= 9; level++) {
                // with minLevel == maxLevel the ladder computes exactly this
                // level and stops there whether or not it converges, so the
                // count is the cost of the level and nothing else
                IntegralResult r = DoubleExponential.integrate1D(f, limits[m][0], limits[m][1], 0.0, level, level);
                assertEquals(member[m] + " level", level, r.level);
                if (previous > 0L) {
                    double ratio = (double) r.evaluations / previous;
                    assertTrue(member[m] + " level " + level + " cost ratio " + ratio,
                            ratio > 1.97 && ratio < 2.03);
                }
                previous = r.evaluations;
            }
        }
    }

    @Test
    public void testAPoleAtANonZeroEndpointIsReportedAsUnresolved() {
        // b - delta collapses onto b once delta falls under one ulp of b, while
        // a + delta stays distinct down to 1e-307 when a is zero. The same
        // integrand therefore succeeds on one half of its interval and not on
        // the other, and the rule has to say which.
        DFunction f = x -> 1.0 / Math.sqrt(x * (1.0 - x));

        IntegralResult left = DoubleExponential.integrate1D(f, 0.0, 0.5, TOL);
        assertTrue("the pole at zero is resolved", left.converged);
        assertTrue(relative(left.value, 0.5 * Math.PI) <= 1.0e-14);

        IntegralResult right = DoubleExponential.integrate1D(f, 0.5, 1.0, TOL);
        assertFalse("the pole at one cannot be resolved and must not be certified", right.converged);
        // it is still worth several digits more than bisection manages
        assertTrue("but it is not nonsense either: " + right.value,
                relative(right.value, 0.5 * Math.PI) < 1.0e-7);
    }

    @Test
    public void testGaussKronrodRemainsTheBetterRuleAwayFromTheEnds() {
        // this rule joins the package beside the adaptive one rather than
        // replacing it, and the test is where that decision is recorded
        DFunction pole = x -> Math.pow(x, -0.9);
        double gkPole = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, pole, 0.0, 1.0, TOL, 22);
        IntegralResult dePole = DoubleExponential.integrate1D(pole, 0.0, 1.0, TOL);
        assertTrue("bisection must lose badly at a pole, got " + gkPole, relative(gkPole, 10.0) > 1.0e-3);
        assertTrue("and this rule must win it", relative(dePole.value, 10.0) <= 1.0e-14);

        DFunction poly = x -> x * x * x * x * x;
        double gkPoly = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, poly, 0.0, 1.0, TOL, 22);
        IntegralResult dePoly = DoubleExponential.integrate1D(poly, 0.0, 1.0, TOL);
        assertEquals("Gauss-Kronrod is exact on a polynomial", 1.0 / 6.0, gkPoly, 0.0);
        assertTrue("a fixed node set is not", dePoly.value != 1.0 / 6.0);
    }

    @Test
    public void testTheValueDoesNotDependOnTheToleranceAsked() {
        // asking for more than double arithmetic can give must not change the
        // answer, only what the flag says about it
        DFunction[] f = { x -> 1.0 / Math.sqrt(x), x -> Math.pow(x, -0.9), x -> Math.exp(x),
                x -> Math.sin(50.0 * x) };
        for (int i = 0; i < f.length; i++) {
            IntegralResult loose = DoubleExponential.integrate1D(f[i], 0.0, 1.0, 1.0e-9);
            IntegralResult tight = DoubleExponential.integrate1D(f[i], 0.0, 1.0, 1.0e-14);
            assertEquals("case " + i, loose.value, tight.value, 1.0e-12 * Math.abs(1.0 + loose.value));
        }
    }

    @Test
    public void testTheErrorEstimateCoversTheTrueErrorToWithinAFewUlp() {
        // not the plain statement: two levels can come out bit identical, which
        // makes the estimate exactly zero while the value is still an ulp or
        // two away from the exact answer. Measured on sqrt(x)*log(1/x) and on
        // exp(-x)/sqrt(x), both of which report an estimate of zero.
        DFunction[] f = { x -> 1.0 / Math.sqrt(x), x -> Math.log(x), x -> Math.sqrt(x) * Math.log(1.0 / x),
                x -> Math.exp(x), x -> Math.sin(50.0 * x) };
        double[] exact = { 2.0, -1.0, 4.0 / 9.0, Math.E - 1.0, (1.0 - Math.cos(50.0)) / 50.0 };
        for (int i = 0; i < f.length; i++) {
            IntegralResult r = DoubleExponential.integrate1D(f[i], 0.0, 1.0, TOL);
            double trueError = Math.abs(r.value - exact[i]);
            double slack = 4.0 * Math.ulp(Math.abs(r.value));
            assertTrue("case " + i + ": estimate " + r.approximatedErrorEstimate + " below true error "
                    + trueError, r.approximatedErrorEstimate + slack >= trueError);
        }
    }

    @Test
    public void testTheLimitsAreChecked() {
        DFunction f = x -> 1.0;
        double[][] bad = { { Double.NaN, 1.0 }, { 0.0, Double.NaN }, { INF, 1.0 }, { 0.0, N_INF } };
        for (double[] limits : bad) {
            try {
                DoubleExponential.integrate1D(f, limits[0], limits[1], TOL);
                fail("[" + limits[0] + ", " + limits[1] + "] must be rejected");
            } catch (IllegalArgumentException expected) {
                // this is the point
            }
        }
    }

    @Test
    public void testDegenerateAndClampedArguments() {
        IntegralResult empty = DoubleExponential.integrate1D(x -> 1.0 / Math.sqrt(x), 1.0, 1.0, TOL);
        assertTrue("an empty interval is zero and known to be", empty.converged);
        assertEquals(0.0, empty.value, 0.0);
        assertEquals(0L, empty.evaluations);

        // maxLevel below minLevel is raised to it rather than skipping the ladder
        IntegralResult clamped = DoubleExponential.integrate1D(x -> Math.exp(x), 0.0, 1.0, 0.0, 5, 1);
        assertEquals(5, clamped.level);

        // and minLevel below the floor is raised to the floor
        IntegralResult floored = DoubleExponential.integrate1D(x -> Math.exp(x), 0.0, 1.0, 0.0, -3, 2);
        assertEquals(2, floored.level);
    }

    @Test
    public void testTheLadderStopsAsSoonAsItCanAndCountsHonestly() {
        // the evaluation count has to be the real one, since the case for this
        // rule rests on it being cheaper and not only more accurate
        final AtomicLong calls = new AtomicLong();
        DFunction counted = x -> {
            calls.incrementAndGet();
            return 1.0 / Math.sqrt(x);
        };
        IntegralResult r = DoubleExponential.integrate1D(counted, 0.0, 1.0, TOL);
        assertTrue(r.converged);
        assertEquals("the reported count must be the real one", calls.get(), r.evaluations);
        assertTrue("17 digits at a pole for " + r.evaluations + " evaluations", r.evaluations < 200L);
    }
}
