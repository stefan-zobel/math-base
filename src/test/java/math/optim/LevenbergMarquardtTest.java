package math.optim;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.MathConsts;
import math.fun.DiffDVectorFunction;
import math.minpack.Lmder_fcn;
import math.minpack.MghProblems;
import math.minpack.MghProblems.Problem;
import math.minpack.Minpack_f77;

/**
 * Tests for {@link LevenbergMarquardt}. The facade adds no arithmetic, so the
 * strongest thing that can be said about it is that it adds no arithmetic: the
 * same problems of {@link MghProblems} are put through both doors -- the raw
 * FORTRAN calling convention below and the facade -- and every number that
 * comes back has to agree <em>exactly</em>. Anything else would be a defect in
 * the index arithmetic, which is all the facade really does.
 * <p>
 * The published minima are asserted separately, so that the public API is
 * stated against the source rather than only against the port.
 */
public class LevenbergMarquardtTest {

    /** The tolerance MINPACK recommends as its lower limit. */
    private static final double SQRT_EPS = Math.sqrt(Math.ulp(1.0));

    /** What the raw port returns, for comparison against the facade. */
    private static final class Raw {

        double[] parameters;
        double[] residuals;
        double sumOfSquares;
        int info;
        int functionEvaluations;
        int jacobianEvaluations;
    }

    /**
     * Calls {@code lmder_f77} the way a caller had to before this class
     * existed: one-based arrays, a row-major Jacobian, results through
     * {@code int[1]} boxes, nineteen positional arguments.
     */
    private static Raw runRaw(Problem p, double tol) {
        int m = p.m;
        int n = p.n;
        double[] x = new double[n + 1];
        System.arraycopy(p.start, 0, x, 1, n);
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        Minpack_f77.lmder_f77(new Bridge(p), m, n, x, fvec, fjac, tol, tol, 0.0, p.maxEvaluations, diag, 1, 100.0, 0,
                info, nfev, njev, ipvt, qtf);

        Raw raw = new Raw();
        raw.parameters = new double[n];
        System.arraycopy(x, 1, raw.parameters, 0, n);
        raw.residuals = new double[m];
        System.arraycopy(fvec, 1, raw.residuals, 0, m);
        double ssq = 0.0;
        for (int i = 0; i < m; i++) {
            ssq += raw.residuals[i] * raw.residuals[i];
        }
        raw.sumOfSquares = ssq;
        raw.info = info[1];
        raw.functionEvaluations = nfev[1];
        raw.jacobianEvaluations = njev[1];
        return raw;
    }

    private static final class Bridge implements Lmder_fcn {

        private final Problem p;
        private final double[] x;
        private final double[] r;
        private final double[] jac;

        Bridge(Problem p) {
            this.p = p;
            this.x = new double[p.n];
            this.r = new double[p.m];
            this.jac = new double[p.m * p.n];
        }

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, double[][] fjac, int[] iflag) {
            System.arraycopy(xOneBased, 1, x, 0, n);
            if (iflag[1] == 1) {
                p.valueAt(x, r);
                System.arraycopy(r, 0, fvec, 1, m);
            } else {
                p.jacobianAt(x, jac);
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < m; i++) {
                        fjac[i + 1][j + 1] = jac[j * m + i];
                    }
                }
            }
        }
    }

    /** The load-bearing test: two doors, one engine, identical numbers. */
    @Test
    public void testTheFacadeReproducesTheRawPortExactly() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Raw raw = runRaw(p, SQRT_EPS);

            LevenbergMarquardt lm = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations, 100.0);
            LevenbergMarquardt.Result r = lm.solve(p, p.start, p.m);

            assertArrayEquals(p.name + ": parameters", raw.parameters, r.parameters, 0.0);
            assertArrayEquals(p.name + ": residuals", raw.residuals, r.residuals, 0.0);
            assertEquals(p.name + ": sum of squares", raw.sumOfSquares, r.sumOfSquares, 0.0);
            assertEquals(p.name + ": function evaluations", raw.functionEvaluations, r.functionEvaluations);
            assertEquals(p.name + ": Jacobian evaluations", raw.jacobianEvaluations, r.jacobianEvaluations);
            assertEquals(p.name + ": status", raw.info, r.status.code());
        }
    }

    /** And the same statement made against the source rather than the port. */
    @Test
    public void testThePublishedMinimaAreReached() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            LevenbergMarquardt lm = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations, 100.0);
            LevenbergMarquardt.Result r = lm.solve(p, p.start, p.m);

            assertTrue(p.name + ": " + r.status, r.converged);
            assertEquals(p.name + ": sum of squares", p.minimum, r.sumOfSquares, p.minimumTolerance);
            if (p.solution != null) {
                for (int j = 0; j < p.n; j++) {
                    double scale = Math.max(Math.abs(p.solution[j]), 1.0);
                    assertEquals(p.name + ": x[" + j + "]", p.solution[j], r.parameters[j],
                            p.solutionTolerance * scale);
                }
            }
        }
    }

    /** The residuals reported belong to the parameters reported. */
    @Test
    public void testTheReportedResidualsBelongToTheReportedParameters() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(p, p.start, p.m);

            double[] expected = new double[p.m];
            p.valueAt(r.parameters, expected);
            assertArrayEquals(p.name, expected, r.residuals, 0.0);
        }
    }

    /**
     * The distinction the integer loses. MINPACK returns 8 for Powell singular,
     * which is larger than every success code and means the scaled gradient is
     * down at machine precision -- a stationary point, and a better outcome
     * than the tolerance that was asked for. Read as a number it sorts after
     * "budget exhausted"; read as a status it does not.
     */
    @Test
    public void testTheStatusSeparatesSuccessFromMagnitude() {
        Problem powellSingular = byName("Powell singular");
        LevenbergMarquardt.Result good = new LevenbergMarquardt().solve(powellSingular, powellSingular.start,
                powellSingular.m);

        assertEquals(LevenbergMarquardt.Status.GRADIENT_AT_MACHINE_PRECISION, good.status);
        assertEquals("the raw code really is the larger one", 8, good.status.code());
        assertTrue("and it is a success", good.converged);
        assertTrue("with 5 reserved for the budget",
                LevenbergMarquardt.Status.TOO_MANY_EVALUATIONS.code() < good.status.code());
        assertFalse(LevenbergMarquardt.Status.TOO_MANY_EVALUATIONS.isSuccess());
    }

    /**
     * An exhausted budget is not reported as success, and the parameters still
     * come back -- on this problem three evaluations already improve on the
     * starting point by four orders of magnitude, and discarding that because
     * the status is not green would be throwing the answer away.
     */
    @Test
    public void testAnExhaustedBudgetKeepsTheParametersItReached() {
        Problem meyer = byName("Meyer");
        LevenbergMarquardt lm = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, 3, 100.0);
        LevenbergMarquardt.Result r = lm.solve(meyer, meyer.start, meyer.m);

        assertEquals(LevenbergMarquardt.Status.TOO_MANY_EVALUATIONS, r.status);
        assertFalse(r.converged);
        assertEquals("parameters", meyer.n, r.parameters.length);
        assertTrue("the point improved on the start: " + r.sumOfSquares + " against "
                + meyer.sumOfSquaresAt(meyer.start), r.sumOfSquares < meyer.sumOfSquaresAt(meyer.start));
    }

    /**
     * All five defaults at once: stating them explicitly must change nothing,
     * on every problem, bit for bit. That pins the tolerances, the step bound
     * and the budget together without needing a problem on which each of them
     * happens to bind.
     */
    @Test
    public void testStatingTheDefaultsExplicitlyChangesNothing() {
        double tol = Math.sqrt(MathConsts.MACH_EPS_DBL);
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            LevenbergMarquardt.Result byDefault = new LevenbergMarquardt().solve(p, p.start, p.m);
            LevenbergMarquardt.Result stated = new LevenbergMarquardt(tol, tol, 0.0, 100 * (p.n + 1), 100.0)
                    .solve(p, p.start, p.m);

            assertArrayEquals(p.name + ": parameters", byDefault.parameters, stated.parameters, 0.0);
            assertEquals(p.name + ": sum of squares", byDefault.sumOfSquares, stated.sumOfSquares, 0.0);
            assertEquals(p.name + ": evaluations", byDefault.functionEvaluations, stated.functionEvaluations);
            assertEquals(p.name + ": status", byDefault.status, stated.status);
        }
    }

    /**
     * And the budget on its own, which the test above cannot separate because
     * none of the eleven problems exhausts it. {@code exp(-x)} falls towards
     * zero without ever reaching it, so no tolerance can fire and the only
     * thing that stops the search is the budget -- which must be MINPACK's own
     * {@code 100 * (n + 1)}.
     */
    @Test
    public void testTheDefaultBudgetIsHundredTimesNPlusOne() {
        DiffDVectorFunction downhillForever = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] values) {
                values[0] = Math.exp(-x[0]);
            }

            @Override
            public void jacobianAt(double[] x, double[] jacobian) {
                jacobian[0] = -Math.exp(-x[0]);
            }
        };

        LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(downhillForever, new double[] { 0.0 }, 1);

        assertEquals(LevenbergMarquardt.Status.TOO_MANY_EVALUATIONS, r.status);
        assertFalse(r.converged);
        assertEquals("the budget for n = 1", 200, r.functionEvaluations);
        assertTrue("and the search really was making progress: " + r.sumOfSquares, r.sumOfSquares < 1.0e-12);
    }

    /** The starting point belongs to the caller and comes back untouched. */
    @Test
    public void testTheStartingPointIsNotModified() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            double[] start = p.start.clone();
            LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(p, start, p.m);

            assertArrayEquals(p.name + ": the caller's array was written to", p.start, start, 0.0);
            assertNotSame(p.name, start, r.parameters);
        }
    }

    /**
     * A tolerance below what double precision can deliver is reported as such
     * rather than as a failure to converge, and the two are different statuses.
     */
    @Test
    public void testATooStringentToleranceIsReportedAsTooStringent() {
        Problem bard = byName("Bard");
        LevenbergMarquardt lm = new LevenbergMarquardt(1.0e-30, 1.0e-30, 0.0, 1000, 100.0);
        LevenbergMarquardt.Result r = lm.solve(bard, bard.start, bard.m);

        assertTrue("status = " + r.status,
                r.status == LevenbergMarquardt.Status.VALUE_TOLERANCE_TOO_SMALL
                        || r.status == LevenbergMarquardt.Status.STEP_TOLERANCE_TOO_SMALL
                        || r.status == LevenbergMarquardt.Status.GRADIENT_AT_MACHINE_PRECISION);
        assertEquals("the minimum was reached all the same", bard.minimum, r.sumOfSquares, bard.minimumTolerance);
    }

    /** A gradient tolerance that is asked for is honoured. */
    @Test
    public void testAGradientToleranceStopsTheSearch() {
        Problem bard = byName("Bard");
        LevenbergMarquardt lm = new LevenbergMarquardt(0.0, 0.0, 1.0e-2, 1000, 100.0);
        LevenbergMarquardt.Result r = lm.solve(bard, bard.start, bard.m);

        assertEquals(LevenbergMarquardt.Status.GRADIENT_TOLERANCE_REACHED, r.status);
        assertTrue(r.converged);
    }

    /** Arguments are checked, not passed on to be misread. */
    @Test
    public void testArgumentValidation() {
        expectIae("valueTolerance", -1.0e-8, 1.0e-8, 0.0, 100, 100.0);
        expectIae("valueTolerance", Double.NaN, 1.0e-8, 0.0, 100, 100.0);
        expectIae("stepTolerance", 1.0e-8, -1.0e-8, 0.0, 100, 100.0);
        expectIae("gradientTolerance", 1.0e-8, 1.0e-8, -1.0, 100, 100.0);
        expectIae("maxEvaluations", 1.0e-8, 1.0e-8, 0.0, -1, 100.0);
        expectIae("initialStepBound", 1.0e-8, 1.0e-8, 0.0, 100, 0.0);
        expectIae("initialStepBound", 1.0e-8, 1.0e-8, 0.0, 100, Double.POSITIVE_INFINITY);

        final Problem p = byName("Rosenbrock");
        LevenbergMarquardt lm = new LevenbergMarquardt();
        expectIae("null function", lm, null, new double[] { 1.0, 1.0 }, 2);
        expectIae("null start", lm, p, null, 2);
        expectIae("empty start", lm, p, new double[0], 2);
        expectIae("m < n", lm, p, new double[] { 1.0, 1.0 }, 1);
        expectIae("NaN in start", lm, p, new double[] { Double.NaN, 1.0 }, 2);
        expectIae("infinity in start", lm, p, new double[] { Double.POSITIVE_INFINITY, 1.0 }, 2);
    }

    /** A residual count above {@code m} is a caller error the facade cannot see. */
    @Test
    public void testAResidualCountLargerThanTheFunctionProducesIsTheCallersProblem() {
        // The function writes m entries; claiming more only means the extra
        // residuals stay zero, which is a valid, if pointless, problem. The
        // point of asserting it is that nothing blows up in the index
        // arithmetic, which is where a facade over one-based FORTRAN would.
        DiffDVectorFunction f = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] values) {
                values[0] = x[0] - 3.0;
                values[1] = x[1] + 2.0;
            }

            @Override
            public void jacobianAt(double[] x, double[] jacobian) {
                for (int k = 0; k < jacobian.length; k++) {
                    jacobian[k] = 0.0;
                }
                jacobian[0 * 4 + 0] = 1.0;
                jacobian[1 * 4 + 1] = 1.0;
            }
        };

        LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(f, new double[] { 0.0, 0.0 }, 4);

        assertEquals("x[0]", 3.0, r.parameters[0], 1.0e-12);
        assertEquals("x[1]", -2.0, r.parameters[1], 1.0e-12);
        assertEquals("residuals", 4, r.residuals.length);
        assertEquals("the two that were never written stay zero", 0.0, r.residuals[2], 0.0);
        assertEquals("the two that were never written stay zero", 0.0, r.residuals[3], 0.0);
    }

    private static void expectIae(String what, double valueTolerance, double stepTolerance, double gradientTolerance,
            int maxEvaluations, double initialStepBound) {
        try {
            new LevenbergMarquardt(valueTolerance, stepTolerance, gradientTolerance, maxEvaluations,
                    initialStepBound);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    private static void expectIae(String what, LevenbergMarquardt lm, DiffDVectorFunction f, double[] start, int m) {
        try {
            lm.solve(f, start, m);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    private static Problem byName(String name) {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            if (name.equals(all[k].name)) {
                return all[k];
            }
        }
        throw new AssertionError(name + " is not in the collection");
    }
}
