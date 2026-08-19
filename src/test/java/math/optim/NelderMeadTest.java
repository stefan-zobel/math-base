package math.optim;

import static org.junit.Assert.*;

import org.junit.Test;

import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;

public class NelderMeadTest {

    /** Extended Rosenbrock, minimum 0 at (1, ..., 1). */
    private static final DMultiFunction ROSENBROCK = x -> {
        double s = 0.0;
        for (int i = 0; i + 1 < x.length; i++) {
            double a = x[i + 1] - x[i] * x[i];
            double b = 1.0 - x[i];
            s += 100.0 * a * a + b * b;
        }
        return s;
    };

    /** Minimum 0 at (1, 2, ..., n). */
    private static final DMultiFunction QUADRATIC = x -> {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            double d = x[i] - (i + 1);
            s += d * d;
        }
        return s;
    };

    /** Powell's singular quartic, minimum 0 at the origin with a singular Hessian. */
    private static final DMultiFunction POWELL = x -> {
        double a = x[0] + 10.0 * x[1];
        double b = x[2] - x[3];
        double c = x[1] - 2.0 * x[2];
        double d = x[0] - x[3];
        return a * a + 5.0 * b * b + c * c * c * c + 10.0 * d * d * d * d;
    };

    /** Himmelblau's function, four minima of value 0. */
    private static final DMultiFunction HIMMELBLAU = x -> {
        double a = x[0] * x[0] + x[1] - 11.0;
        double b = x[0] + x[1] * x[1] - 7.0;
        return a * a + b * b;
    };

    private static final double[][] HIMMELBLAU_MINIMA = {
        { 3.0, 2.0 },
        { -2.805118086952745, 3.131312518250573 },
        { -3.779310253377747, -3.283185991286170 },
        { 3.584428340330492, -1.848126526964403 }
    };

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testRosenbrock2D() {
        DMultiFunctionEval r = new NelderMead().minimize(ROSENBROCK, new double[] { -1.2, 1.0 });
        assertTrue(r.converged);
        assertTrue("f = " + r.value, r.value < 1.0e-10);
        assertArrayEquals(new double[] { 1.0, 1.0 }, r.point, 1.0e-6);
    }

    @Test
    public void testRosenbrock5D() {
        // n >= 3, so this runs on the adaptive coefficients
        double[] start = { -1.2, 1.0, -1.2, 1.0, -1.2 };
        DMultiFunctionEval r = new NelderMead().minimize(ROSENBROCK, start);
        assertTrue(r.converged);
        assertTrue("f = " + r.value, r.value < 1.0e-10);
        assertArrayEquals(new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 }, r.point, 1.0e-6);
    }

    @Test
    public void testQuadraticInSeveralDimensions() {
        // n == 1 uses the classic coefficients, n == 8 the adaptive ones
        int[] dims = { 1, 2, 8 };
        for (int k = 0; k < dims.length; k++) {
            int n = dims[k];
            DMultiFunctionEval r = new NelderMead().minimize(QUADRATIC, new double[n]);
            assertTrue("n = " + n, r.converged);
            assertTrue("n = " + n + ", f = " + r.value, r.value < 1.0e-12);
            double[] expected = new double[n];
            for (int i = 0; i < n; i++) {
                expected[i] = i + 1;
            }
            assertArrayEquals("n = " + n, expected, r.point, 1.0e-6);
        }
    }

    @Test
    public void testPowellSingularQuartic() {
        // the quartic valley makes the point far less accurate than the value
        DMultiFunctionEval r = new NelderMead().minimize(POWELL, new double[] { 3.0, -1.0, 0.0, 1.0 });
        assertTrue(r.converged);
        assertTrue("f = " + r.value, r.value < 1.0e-16);
        assertArrayEquals(new double[4], r.point, 1.0e-4);
    }

    @Test
    public void testHimmelblauFindsTheNearestMinimum() {
        double[][] starts = { { 0.0, 0.0 }, { -4.0, 4.0 }, { 4.0, -4.0 }, { -4.0, -4.0 } };
        NelderMead nm = new NelderMead();
        for (int k = 0; k < starts.length; k++) {
            DMultiFunctionEval r = nm.minimize(HIMMELBLAU, starts[k]);
            assertTrue(r.converged);
            assertTrue("f = " + r.value, r.value < 1.0e-10);
            boolean matched = false;
            for (int m = 0; m < HIMMELBLAU_MINIMA.length; m++) {
                double dx = Math.abs(r.point[0] - HIMMELBLAU_MINIMA[m][0]);
                double dy = Math.abs(r.point[1] - HIMMELBLAU_MINIMA[m][1]);
                if (dx < 1.0e-6 && dy < 1.0e-6) {
                    matched = true;
                }
            }
            assertTrue("start " + k + " settled at (" + r.point[0] + ", " + r.point[1] + ")", matched);
        }
    }

    @Test
    public void testTranslationInvariance() {
        double[] shift = { 2.5, -7.0, 0.75 };
        DMultiFunction shifted = x -> {
            double[] y = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                y[i] = x[i] - shift[i];
            }
            return QUADRATIC.apply(y);
        };
        double[] step = { 0.1, 0.1, 0.1 };
        NelderMead nm = new NelderMead();
        DMultiFunctionEval plain = nm.minimize(QUADRATIC, new double[3], step);
        DMultiFunctionEval moved = nm.minimize(shifted, shift.clone(), step);
        assertTrue(plain.converged && moved.converged);
        for (int i = 0; i < 3; i++) {
            assertEquals(plain.point[i] + shift[i], moved.point[i], 1.0e-7);
        }
    }

    @Test
    public void testExplicitStepOnABadlyScaledProblem() {
        // the optimum is (1, 1e6), so a simplex of the same edge in both
        // coordinates is a poor fit
        DMultiFunction f = x -> {
            double a = x[0] - 1.0;
            double b = x[1] * 1.0e-6 - 1.0;
            return a * a + b * b;
        };
        NelderMead nm = new NelderMead();
        DMultiFunctionEval auto = nm.minimize(f, new double[2]);
        DMultiFunctionEval scaled = nm.minimize(f, new double[2], new double[] { 0.1, 1.0e5 });
        assertTrue(auto.converged && scaled.converged);
        assertEquals(1.0, scaled.point[0], 1.0e-6);
        assertEquals(1.0e6, scaled.point[1], 1.0);
        assertTrue("steps matched to the problem should not need more iterations",
                scaled.iterations < auto.iterations);
    }

    @Test
    public void testWalksAwayFromARegionWhereTheFunctionIsUndefined() {
        // NaN sorts as the worst vertex and is never accepted as an improvement
        DMultiFunction f = x -> {
            if (x[0] < -1.0) {
                return Double.NaN;
            }
            return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] + 3.0) * (x[1] + 3.0);
        };
        DMultiFunctionEval r = new NelderMead().minimize(f, new double[2]);
        assertTrue(r.converged);
        assertArrayEquals(new double[] { 2.0, -3.0 }, r.point, 1.0e-6);
    }

    @Test
    public void testIsDeterministic() {
        NelderMead nm = new NelderMead();
        DMultiFunctionEval a = nm.minimize(ROSENBROCK, new double[] { -1.2, 1.0, 0.5 });
        DMultiFunctionEval b = nm.minimize(ROSENBROCK, new double[] { -1.2, 1.0, 0.5 });
        assertArrayEquals(a.point, b.point, 0.0);
        assertEquals(a.value, b.value, 0.0);
        assertEquals(a.iterations, b.iterations);
    }

    @Test
    public void testExhaustedBudgetReportsFailureButKeepsTheBestPoint() {
        double[] start = { -1.2, 1.0 };
        DMultiFunctionEval r = new NelderMead(1.0e-10, 1.0e-8, 5, 0).minimize(ROSENBROCK, start);
        assertFalse(r.converged);
        assertEquals(5, r.iterations);
        assertTrue(r.value < ROSENBROCK.apply(start));
        assertFalse(Double.isNaN(r.point[0]) || Double.isNaN(r.point[1]));
    }

    @Test
    public void testArgumentValidation() {
        NelderMead nm = new NelderMead();
        expectIae(() -> nm.minimize(null, new double[2]));
        expectIae(() -> nm.minimize(QUADRATIC, null));
        expectIae(() -> nm.minimize(QUADRATIC, new double[0]));
        expectIae(() -> nm.minimize(QUADRATIC, new double[] { Double.NaN, 0.0 }));
        expectIae(() -> nm.minimize(QUADRATIC, new double[] { Double.POSITIVE_INFINITY, 0.0 }));
        expectIae(() -> nm.minimize(QUADRATIC, new double[2], new double[3]));
        expectIae(() -> nm.minimize(QUADRATIC, new double[2], new double[] { 0.1, 0.0 }));
        expectIae(() -> nm.minimize(x -> Double.NaN, new double[2]));
        expectIae(() -> new NelderMead(-1.0, 1.0e-8, 0, 1));
        expectIae(() -> new NelderMead(1.0e-10, Double.NaN, 0, 1));
        expectIae(() -> new NelderMead(1.0e-10, 1.0e-8, 0, -1));
    }

    private static void expectIae(Runnable r) {
        try {
            r.run();
            fail("expected an IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // that is the point of the call
        }
    }
}
