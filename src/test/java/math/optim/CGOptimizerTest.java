package math.optim;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

import math.fun.DMultiFunctionEval;
import math.fun.DiffDMultiFunction;

/**
 * Characterization tests for {@link CGOptimizer}, which had no coverage at all.
 * They record what the class does rather than what it ought to do, evaluation
 * counts included, so that lifting the one-dimensional minimizer out of it can
 * be shown to change nothing. A count that moves is either explained or a bug.
 */
public class CGOptimizerTest {

    private int values;
    private int derivatives;

    @Before
    public void resetCounters() {
        values = 0;
        derivatives = 0;
    }

    /** Counts every evaluation so the cost of a change is visible. */
    private abstract class Counting implements DiffDMultiFunction {

        @Override
        public final double apply(double[] x) {
            values++;
            return f(x);
        }

        @Override
        public final void derivativeAt(double[] x, double[] grad) {
            derivatives++;
            d(x, grad);
        }

        abstract double f(double[] x);

        abstract void d(double[] x, double[] grad);
    }

    /** {@code 3(x-2)^2 + 5(y+1)^2}, minimal at {@code (2, -1)}. */
    private final class Quadratic extends Counting {

        @Override
        double f(double[] x) {
            return 3.0 * (x[0] - 2.0) * (x[0] - 2.0) + 5.0 * (x[1] + 1.0) * (x[1] + 1.0);
        }

        @Override
        void d(double[] x, double[] grad) {
            grad[0] = 6.0 * (x[0] - 2.0);
            grad[1] = 10.0 * (x[1] + 1.0);
        }
    }

    /** The same quadratic negated, so its maximum is at {@code (2, -1)}. */
    private final class NegatedQuadratic extends Counting {

        @Override
        double f(double[] x) {
            return -(3.0 * (x[0] - 2.0) * (x[0] - 2.0) + 5.0 * (x[1] + 1.0) * (x[1] + 1.0));
        }

        @Override
        void d(double[] x, double[] grad) {
            grad[0] = -6.0 * (x[0] - 2.0);
            grad[1] = -10.0 * (x[1] + 1.0);
        }
    }

    /** Rosenbrock in any dimension, minimal at the all-ones point. */
    private final class Rosenbrock extends Counting {

        @Override
        double f(double[] x) {
            double s = 0.0;
            for (int i = 0; i + 1 < x.length; i++) {
                double a = x[i + 1] - x[i] * x[i];
                double b = 1.0 - x[i];
                s += 100.0 * a * a + b * b;
            }
            return s;
        }

        @Override
        void d(double[] x, double[] grad) {
            for (int i = 0; i < grad.length; i++) {
                grad[i] = 0.0;
            }
            for (int i = 0; i + 1 < x.length; i++) {
                double a = x[i + 1] - x[i] * x[i];
                grad[i] += -400.0 * x[i] * a - 2.0 * (1.0 - x[i]);
                grad[i + 1] += 200.0 * a;
            }
        }
    }

    /** {@code s*x^2 + y^2}: raising {@code s} makes the problem badly scaled. */
    private final class Scaled extends Counting {

        private final double s;

        Scaled(double s) {
            this.s = s;
        }

        @Override
        double f(double[] x) {
            return s * x[0] * x[0] + x[1] * x[1];
        }

        @Override
        void d(double[] x, double[] grad) {
            grad[0] = 2.0 * s * x[0];
            grad[1] = 2.0 * x[1];
        }
    }

    /** Unbounded below, so no minimum exists along any descent direction. */
    private final class Unbounded extends Counting {

        @Override
        double f(double[] x) {
            return -(x[0] * x[0] + x[1] * x[1]);
        }

        @Override
        void d(double[] x, double[] grad) {
            grad[0] = -2.0 * x[0];
            grad[1] = -2.0 * x[1];
        }
    }

    /** The well behaved case, pinned exactly. */
    @Test
    public void testQuadratic() {
        DMultiFunctionEval r = CGOptimizer.minimize(new Quadratic(), new double[] { 5.0, 5.0 });

        assertEquals("x", 2.0, r.point[0], 0.0);
        assertEquals("y", -1.0, r.point[1], 1.0e-15);
        assertEquals("value", 0.0, r.value, 1.0e-30);
        assertTrue(r.converged);
        assertEquals("iterations", 8, r.iterations);
        assertEquals("function evaluations", 127, values);
        assertEquals("derivative evaluations", 89, derivatives);
    }

    /**
     * {@code maximize} is {@code minimize} on the negated function, and the sign
     * is put back on the way out. The only coverage the {@code Minimand} wrapper
     * has ever had.
     */
    @Test
    public void testMaximizeMirrorsMinimize() {
        DMultiFunctionEval r = CGOptimizer.maximize(new NegatedQuadratic(), new double[] { 5.0, 5.0 });

        assertEquals("x", 2.0, r.point[0], 0.0);
        assertEquals("y", -1.0, r.point[1], 1.0e-15);
        assertEquals("value", 0.0, r.value, 1.0e-30);
        assertTrue("the reported value must be the maximum, not its negation", r.value <= 0.0);
        assertTrue(r.converged);
        assertEquals("iterations", 8, r.iterations);
        assertEquals("function evaluations", 127, values);
        assertEquals("derivative evaluations", 89, derivatives);
    }

    /**
     * Started on the minimum, the line minimum sits at {@code x == 0}. That used
     * to be the worst case: the purely relative tolerance collapsed to zero
     * there, the convergence test could never fire, and confirming the answer
     * cost 108 function and 103 derivative evaluations. With the absolute floor
     * under the tolerance restored it costs 32 and 29.
     */
    @Test
    public void testStartingOnTheMinimumIsCheapToConfirm() {
        DMultiFunctionEval r = CGOptimizer.minimize(new Quadratic(), new double[] { 2.0, -1.0 });

        assertEquals("x", 2.0, r.point[0], 0.0);
        assertEquals("y", -1.0, r.point[1], 0.0);
        assertEquals("value", 0.0, r.value, 0.0);
        assertTrue(r.converged);
        assertEquals("iterations", 1, r.iterations);
        assertEquals("function evaluations", 32, values);
        assertEquals("derivative evaluations", 29, derivatives);
    }

    /** Non-convex, and the class handles it. */
    @Test
    public void testRosenbrock2d() {
        DMultiFunctionEval r = CGOptimizer.minimize(new Rosenbrock(), new double[] { -1.2, 1.0 });

        assertEquals("x", 1.0, r.point[0], 1.0e-14);
        assertEquals("y", 1.0, r.point[1], 1.0e-14);
        assertEquals("value", 0.0, r.value, 1.0e-30);
        assertTrue(r.converged);
        assertEquals("iterations", 31, r.iterations);
        assertEquals("function evaluations", 320, values);
        assertEquals("derivative evaluations", 225, derivatives);
    }

    /**
     * The same in four dimensions. This is the one case whose trajectory the
     * absolute tolerance floor visibly changes: two more iterations, 28 fewer
     * function and 33 more derivative evaluations, and a final value of 1.6e-27
     * where it used to be 1.5e-28. Both are the minimum to thirteen digits --
     * a different path down a non-convex valley, not a loss of accuracy.
     */
    @Test
    public void testRosenbrock4d() {
        DMultiFunctionEval r = CGOptimizer.minimize(new Rosenbrock(), new double[] { 3.0, -1.0, 0.0, 1.0 });

        for (int i = 0; i < 4; i++) {
            assertEquals("point[" + i + "]", 1.0, r.point[i], 1.0e-12);
        }
        assertEquals("value", 0.0, r.value, 1.0e-26);
        assertTrue(r.converged);
        assertEquals("iterations", 72, r.iterations);
        assertEquals("function evaluations", 744, values);
        assertEquals("derivative evaluations", 514, derivatives);
    }

    /**
     * Bad scaling costs the method nothing -- conditioning up to {@code 1e8}
     * converges in fewer iterations than the well scaled case, not more.
     */
    @Test
    public void testBadScalingIsHandled() {
        int[] expectedIterations = { 3, 9, 7, 6, 6 };
        int[] expectedValues = { 125, 85, 69, 62, 62 };
        int[] expectedDerivatives = { 11, 32, 27, 23, 23 };
        double[] scales = { 1.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8 };

        for (int k = 0; k < scales.length; k++) {
            resetCounters();
            DMultiFunctionEval r = CGOptimizer.minimize(new Scaled(scales[k]), new double[] { 1.0, 1.0 });

            assertEquals("s=" + scales[k] + " x", 0.0, r.point[0], 1.0e-20);
            assertEquals("s=" + scales[k] + " y", 0.0, r.point[1], 1.0e-20);
            assertEquals("s=" + scales[k] + " value", 0.0, r.value, 1.0e-60);
            assertTrue("s=" + scales[k], r.converged);
            assertEquals("s=" + scales[k] + " iterations", expectedIterations[k], r.iterations);
            assertEquals("s=" + scales[k] + " function evaluations", expectedValues[k], values);
            assertEquals("s=" + scales[k] + " derivative evaluations", expectedDerivatives[k], derivatives);
        }
    }

    /**
     * No minimum exists along any descent direction. The bracketing used to walk
     * out of the finite range unnoticed -- its loop ends when {@code fb > fc} is
     * false, which a NaN satisfies just as well as a real bracket -- and the run
     * then spent the full budget of 10001 iterations and 1 060 121 function
     * evaluations arriving at a NaN. It now stops at the first direction that
     * cannot be bracketed, keeps the point it started from, and says it did not
     * converge.
     */
    @Test
    public void testUnboundedFunctionStopsAtOnceInsteadOfReturningNaN() {
        DMultiFunctionEval r = CGOptimizer.minimize(new Unbounded(), new double[] { 1.0, 1.0 });

        assertEquals("x", 1.0, r.point[0], 0.0);
        assertEquals("y", 1.0, r.point[1], 0.0);
        assertEquals("value", -2.0, r.value, 0.0);
        assertFalse("nothing was minimized, and that is reported", r.converged);
        assertEquals("iterations", 1, r.iterations);
        assertEquals("function evaluations", 205, values);
        assertEquals("derivative evaluations", 1, derivatives);
    }
}
