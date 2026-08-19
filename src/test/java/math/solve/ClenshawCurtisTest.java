package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.concurrent.atomic.AtomicLong;

import org.junit.Test;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * Tests for the refining Clenshaw-Curtis rule. The properties asserted here are
 * the ones the previous fixed-grid implementation inside MetaIntegrator could
 * not offer: an error estimate, an honest non-convergence signal, and reuse of
 * every evaluation across the grid levels.
 */
public class ClenshawCurtisTest {

    private static final double TOL = 1.0e-10;

    /** Closed form of the integral of {@code sin(k*x)} over {@code [0, 1]}. */
    private static double sineExact(double k) {
        return (1.0 - Math.cos(k)) / k;
    }

    @Test
    public void testPolynomialExactnessAtTheFirstGrid() {
        for (int p = 0; p <= 8; p++) {
            final int degree = p;
            DFunction f = x -> Math.pow(x, degree);
            double exact = Math.pow(2.0, degree + 1) / (degree + 1);

            ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate1D(f, 0.0, 2.0, TOL);

            assertTrue("x^" + degree + " must converge", r.converged);
            assertEquals("x^" + degree, exact, r.value, 1.0e-13 * exact);
            assertEquals("a polynomial must not need more than the first two grids", 65, r.points);
        }
    }

    @Test
    public void testSmoothIntegrandCostsNoMoreThanTheOldFixedGrid() {
        // The fixed grid this rule replaced used 65 / 33 / 17 points per axis.
        // A smooth integrand has to stop exactly there, not one level later.
        ClenshawCurtis.IntegralResult r1 = ClenshawCurtis.integrate1D(x -> Math.exp(x), 0.0, 1.0, TOL);
        assertTrue(r1.converged);
        assertEquals(65, r1.points);
        assertEquals(Math.E - 1.0, r1.value, 1.0e-13);

        ClenshawCurtis.IntegralResult r2 = ClenshawCurtis.integrate2D((x, y) -> x * x + y * y, 0.0, 2.0, 0.0, 2.0, TOL);
        assertTrue(r2.converged);
        assertEquals(33, r2.points);
        assertEquals(32.0 / 3.0, r2.value, 1.0e-12);

        ClenshawCurtis.IntegralResult r3 = ClenshawCurtis.integrate3D((x, y, z) -> x * y * z, 0.0, 1.0, 0.0, 1.0, 0.0,
                1.0, TOL);
        assertTrue(r3.converged);
        assertEquals(17, r3.points);
        assertEquals(0.125, r3.value, 1.0e-14);
    }

    @Test
    public void testFastOscillationIsResolvedInsteadOfAliased() {
        // On the old fixed grid of 65 points these came back with relative
        // errors of 3.6e-04 (k = 200) and 2.2e+01 with the wrong sign
        // (k = 1000), reported as if they were converged.
        for (int k : new int[] { 30, 200, 1000 }) {
            final double freq = k;
            DFunction f = x -> Math.sin(freq * x);
            double exact = sineExact(freq);

            ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate1D(f, 0.0, 1.0, TOL);

            assertTrue("sin(" + k + "x) must converge", r.converged);
            assertEquals("sin(" + k + "x)", exact, r.value, 1.0e-10 * Math.abs(exact));
        }
    }

    @Test
    public void testOscillationBeyondTheFinestGridIsReportedAsNotConverged() {
        // sin(5000x) has 1592 half waves on [0, 1] and the finest grid has 2049
        // points. The point of the whole class is that this is signalled rather
        // than returned as a converged value.
        ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate1D(x -> Math.sin(5000.0 * x), 0.0, 1.0, TOL);

        assertFalse("an unresolvable oscillation must not claim convergence", r.converged);
        assertEquals("the finest grid must have been reached", 2049, r.points);
        assertTrue("the estimate has to be above the tolerance", r.approximatedErrorEstimate > TOL);
    }

    @Test
    public void testErrorEstimateCoversTheTrueError() {
        DFunction[] fs = { x -> Math.exp(x), x -> 1.0 / (1.0 + x * x), x -> Math.sin(3.0 * x),
                x -> Math.pow(x, 9) };
        double[] exact = { Math.E - 1.0, Math.PI / 4.0, (1.0 - Math.cos(3.0)) / 3.0, 0.1 };

        for (int i = 0; i < fs.length; i++) {
            ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate1D(fs[i], 0.0, 1.0, TOL);
            double trueError = Math.abs(r.value - exact[i]);

            assertTrue("case " + i + ": estimate " + r.approximatedErrorEstimate + " below true error " + trueError,
                    trueError <= Math.max(r.approximatedErrorEstimate, 1.0e-13));
        }
    }

    @Test
    public void testNestedGridsReuseEveryEvaluation() {
        // Grid N is a subset of grid 2N, so the whole ladder must cost exactly
        // as much as its finest level - no point twice, none missed.
        AtomicLong calls = new AtomicLong();

        DFunction f1 = x -> {
            calls.incrementAndGet();
            return Math.exp(x);
        };
        ClenshawCurtis.IntegralResult r1 = ClenshawCurtis.integrate1D(f1, 0.0, 1.0, TOL);
        assertEquals("1D", r1.points, calls.get());

        calls.set(0L);
        DBiFunction f2 = (x, y) -> {
            calls.incrementAndGet();
            return Math.exp(x + y);
        };
        ClenshawCurtis.IntegralResult r2 = ClenshawCurtis.integrate2D(f2, 0.0, 1.0, 0.0, 1.0, TOL);
        assertEquals("2D", (long) r2.points * r2.points, calls.get());

        calls.set(0L);
        DTriFunction f3 = (x, y, z) -> {
            calls.incrementAndGet();
            return Math.exp(x + y + z);
        };
        ClenshawCurtis.IntegralResult r3 = ClenshawCurtis.integrate3D(f3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOL);
        assertEquals("3D", (long) r3.points * r3.points * r3.points, calls.get());
    }

    @Test
    public void testSeparableIntegrandsMatchTheProductOfTheirFactors() {
        double oneD = Math.E - 1.0;

        ClenshawCurtis.IntegralResult r2 = ClenshawCurtis.integrate2D((x, y) -> Math.exp(x) * Math.exp(y), 0.0, 1.0,
                0.0, 1.0, TOL);
        assertEquals(oneD * oneD, r2.value, 1.0e-12);

        ClenshawCurtis.IntegralResult r3 = ClenshawCurtis.integrate3D(
                (x, y, z) -> Math.exp(x) * Math.exp(y) * Math.exp(z), 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOL);
        assertEquals(oneD * oneD * oneD, r3.value, 1.0e-12);
    }

    @Test
    public void testOscillationInTwoAndThreeDimensions() {
        for (int k : new int[] { 20, 60, 200 }) {
            final double freq = k;
            DBiFunction f = (x, y) -> Math.sin(freq * x) * Math.cos(freq * y);
            double exact = sineExact(freq) * (Math.sin(freq) / freq);

            ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate2D(f, 0.0, 1.0, 0.0, 1.0, TOL);

            assertTrue("2D k = " + k, r.converged);
            assertEquals("2D k = " + k, exact, r.value, 1.0e-10 * Math.abs(exact));
        }

        for (int k : new int[] { 10, 35 }) {
            final double freq = k;
            DTriFunction f = (x, y, z) -> Math.sin(freq * x) * Math.cos(freq * y) * Math.sin(freq * z);
            double exact = sineExact(freq) * (Math.sin(freq) / freq) * sineExact(freq);

            ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate3D(f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOL);

            assertTrue("3D k = " + k, r.converged);
            assertEquals("3D k = " + k, exact, r.value, 1.0e-9 * Math.abs(exact));
        }
    }

    @Test
    public void testExplicitGridRangeIsRoundedAndAlwaysAllowsOneDoubling() {
        DFunction f = x -> Math.exp(x);

        // 100 rounds down to 64, and a maximum below twice the minimum is
        // raised, since a single grid cannot produce an estimate
        ClenshawCurtis.IntegralResult r = ClenshawCurtis.integrate1D(f, 0.0, 1.0, TOL, 100, 100);
        assertEquals(129, r.points);
        assertTrue(r.converged);
        assertEquals(Math.E - 1.0, r.value, 1.0e-13);

        // out of range on both ends: clamped to 8 and to 2048
        ClenshawCurtis.IntegralResult low = ClenshawCurtis.integrate1D(f, 0.0, 1.0, TOL, 1, 3);
        assertEquals(17, low.points);

        ClenshawCurtis.IntegralResult high = ClenshawCurtis.integrate1D(x -> Math.sin(5000.0 * x), 0.0, 1.0, TOL, 4096,
                999999);
        assertEquals(2049, high.points);
    }

    @Test
    public void testRepeatedCallsAgreeBitForBit() {
        DBiFunction f2 = (x, y) -> Math.sin(20.0 * x) * Math.cos(20.0 * y);
        assertEquals(ClenshawCurtis.integrate2D(f2, 0.0, 1.0, 0.0, 1.0, TOL).value,
                ClenshawCurtis.integrate2D(f2, 0.0, 1.0, 0.0, 1.0, TOL).value, 0.0);

        DTriFunction f3 = (x, y, z) -> Math.exp(x + y) * (z * z + x * y);
        assertEquals(ClenshawCurtis.integrate3D(f3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOL).value,
                ClenshawCurtis.integrate3D(f3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TOL).value, 0.0);
    }
}
