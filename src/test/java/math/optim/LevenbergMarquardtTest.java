package math.optim;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.MathConsts;
import math.fun.DVectorFunction;
import math.fun.DiffDVectorFunction;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.minpack.Lmder_fcn;
import math.minpack.Lmdif_fcn;
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

    /** The same for {@code lmdif_f77}, which asks only for values. */
    private static Raw runRawWithoutJacobian(Problem p, double tol) {
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

        Minpack_f77.lmdif_f77(new Bridge(p), m, n, x, fvec, tol, tol, 0.0, p.maxEvaluations, 0.0, diag, 1, 100.0, 0,
                info, nfev, fjac, ipvt, qtf);

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
        raw.jacobianEvaluations = 0;
        return raw;
    }

    private static final class Bridge implements Lmder_fcn, Lmdif_fcn {

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

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, int[] iflag) {
            System.arraycopy(xOneBased, 1, x, 0, n);
            p.valueAt(x, r);
            System.arraycopy(r, 0, fvec, 1, m);
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

    /** The same for the derivative-free path, which goes through {@code lmdif}. */
    @Test
    public void testTheDerivativeFreePathReproducesTheRawPortExactly() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Raw raw = runRawWithoutJacobian(p, SQRT_EPS);

            LevenbergMarquardt lm = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations, 100.0);
            LevenbergMarquardt.Result r = lm.solve((DVectorFunction) p, p.start, p.m);

            assertArrayEquals(p.name + ": parameters", raw.parameters, r.parameters, 0.0);
            assertArrayEquals(p.name + ": residuals", raw.residuals, r.residuals, 0.0);
            assertEquals(p.name + ": sum of squares", raw.sumOfSquares, r.sumOfSquares, 0.0);
            assertEquals(p.name + ": function evaluations", raw.functionEvaluations, r.functionEvaluations);
            assertEquals(p.name + ": status", raw.info, r.status.code());
            assertEquals(p.name + ": the Jacobian is approximated, never evaluated", 0, r.jacobianEvaluations);
        }
    }

    /**
     * Not supplying a derivative reaches the same minima and never costs less.
     * The agreement asserted is on the minimum rather than on the point, and
     * loosely: a forward difference is a different algorithm, not a slower way
     * of running the same one.
     */
    @Test
    public void testTheDerivativeFreePathCostsMoreAndFindsTheSame() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            LevenbergMarquardt lm = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations, 100.0);
            LevenbergMarquardt.Result withJacobian = lm.solve(p, p.start, p.m);
            LevenbergMarquardt.Result withoutJacobian = lm.solve((DVectorFunction) p, p.start, p.m);

            double scale = Math.max(Math.abs(withJacobian.sumOfSquares), 1.0e-8);
            assertEquals(p.name + ": sum of squares", withJacobian.sumOfSquares, withoutJacobian.sumOfSquares,
                    1.0e-4 * scale);
            assertTrue(p.name + ": " + withoutJacobian.functionEvaluations + " against "
                    + withJacobian.functionEvaluations,
                    withoutJacobian.functionEvaluations >= withJacobian.functionEvaluations);
            assertTrue(p.name + ": a Jacobian was evaluated on the analytic path",
                    withJacobian.jacobianEvaluations > 0);
        }
    }

    /**
     * Which of the two overloads runs is decided by the <em>static</em> type of
     * the argument, as Java always decides it. Worth an executable statement
     * rather than only a sentence of documentation, because the same object
     * takes both paths here and the cheaper one is not the default it looks
     * like.
     */
    @Test
    public void testTheOverloadIsChosenByTheStaticType() {
        Problem rosenbrock = byName("Rosenbrock");
        DVectorFunction sameObjectSeenAsValuesOnly = rosenbrock;
        LevenbergMarquardt lm = new LevenbergMarquardt();

        LevenbergMarquardt.Result declared = lm.solve(rosenbrock, rosenbrock.start, rosenbrock.m);
        LevenbergMarquardt.Result widened = lm.solve(sameObjectSeenAsValuesOnly, rosenbrock.start, rosenbrock.m);

        assertTrue("the declared type carries the Jacobian", declared.jacobianEvaluations > 0);
        assertEquals("the widened one does not", 0, widened.jacobianEvaluations);
        assertTrue("and pays for it: " + widened.functionEvaluations + " against " + declared.functionEvaluations,
                widened.functionEvaluations > declared.functionEvaluations);
        for (int j = 0; j < rosenbrock.n; j++) {
            assertEquals("both still land on the minimum", 1.0, declared.parameters[j], 1.0e-8);
            assertEquals("both still land on the minimum", 1.0, widened.parameters[j], 1.0e-8);
        }
    }

    /**
     * The sharpest edge on the derivative-free path, and the reason
     * {@code functionAccuracy} is a constructor argument rather than something
     * left at MINPACK's default.
     * <p>
     * The residuals below are a straight line whose values are rounded to 1e-6,
     * which is what a residual computed by a simulation, a quadrature or a
     * lookup table looks like. Left to assume machine precision, the difference
     * quotient is taken over a step of about 1.5e-8, every rounded value comes
     * out unchanged, and the approximated Jacobian is <em>identically zero</em>.
     * A zero matrix is orthogonal to every residual, so MINPACK's gradient test
     * fires, and the run reports {@code GRADIENT_TOLERANCE_REACHED} -- a
     * success -- after three evaluations, on the starting point, having moved
     * nothing.
     * <p>
     * Told that the values are only good to 1e-6, the same solver picks a step
     * of 1e-3 instead and lands on the exact answer in seven evaluations.
     */
    @Test
    public void testANoisyResidualNeedsToBeDeclaredOrTheFitSilentlySucceedsAtTheStart() {
        DVectorFunction quantizedLine = new DVectorFunction() {
            @Override
            public void valueAt(double[] p, double[] r) {
                for (int i = 0; i < 20; i++) {
                    double t = i / 10.0;
                    double coarse = Math.rint((p[0] * t + p[1]) * 1.0e6) / 1.0e6;
                    r[i] = coarse - (2.0 * t - 1.0);
                }
            }
        };
        double[] start = { 0.0, 0.0 };

        LevenbergMarquardt assumingFullPrecision = new LevenbergMarquardt(1.0e-10, 1.0e-10, 0.0, 2000, 100.0);
        LevenbergMarquardt.Result trap = assumingFullPrecision.solve(quantizedLine, start, 20);

        assertEquals("a zero Jacobian is orthogonal to everything, so this reports success",
                LevenbergMarquardt.Status.GRADIENT_TOLERANCE_REACHED, trap.status);
        assertTrue("and it is a success by every measure the library has", trap.converged);
        assertArrayEquals("on the starting point, unmoved", start, trap.parameters, 0.0);
        assertTrue("with the sum of squares it started at: " + trap.sumOfSquares, trap.sumOfSquares > 1.0);

        LevenbergMarquardt told = new LevenbergMarquardt(1.0e-10, 1.0e-10, 0.0, 2000, 100.0, 1.0e-6);
        LevenbergMarquardt.Result fit = told.solve(quantizedLine, start, 20);

        assertTrue(fit.converged);
        assertEquals("slope", 2.0, fit.parameters[0], 1.0e-9);
        assertEquals("intercept", -1.0, fit.parameters[1], 1.0e-9);
        assertTrue("and cheaply: " + fit.functionEvaluations, fit.functionEvaluations < 20);
    }

    /**
     * {@code functionAccuracy} is meaningless where the derivative is supplied,
     * and changing it must therefore change nothing at all on that path.
     */
    @Test
    public void testFunctionAccuracyDoesNotTouchTheAnalyticPath() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            LevenbergMarquardt.Result byDefault = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations,
                    100.0).solve(p, p.start, p.m);
            LevenbergMarquardt.Result declared = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations,
                    100.0, 1.0e-3).solve(p, p.start, p.m);

            assertArrayEquals(p.name, byDefault.parameters, declared.parameters, 0.0);
            assertEquals(p.name, byDefault.functionEvaluations, declared.functionEvaluations);
            assertEquals(p.name, byDefault.status, declared.status);
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

    /** The inline generator the test conventions ask for, seeded per test. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            this.state = seed;
        }

        /** Uniform on {@code [-0.5, 0.5)}. */
        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) * 0x1.0p-53) - 0.5;
        }
    }

    /** A linear model {@code X b}, which is the one case with an outside answer. */
    private static final class LinearModel implements DiffDVectorFunction {

        final DMatrix x;
        final DMatrix y;
        private final double[] design;
        private final double[] target;
        private final int m;
        private final int n;

        LinearModel(int m, int n, long seed, double noise) {
            this.m = m;
            this.n = n;
            Lcg lcg = new Lcg(seed);
            double[] beta = { 2.0, -1.5, 0.5, 3.0 };
            x = new DMatrix(m, n);
            for (int i = 0; i < m; i++) {
                x.set(i, 0, 1.0);
                for (int j = 1; j < n; j++) {
                    x.set(i, j, 4.0 * lcg.next());
                }
            }
            y = new DMatrix(m, 1);
            for (int i = 0; i < m; i++) {
                double v = 0.0;
                for (int j = 0; j < n; j++) {
                    v += x.get(i, j) * beta[j];
                }
                y.set(i, 0, v + noise * lcg.next());
            }
            // DMatrix stores column-major with idx(i, j) == j * rows + i, which
            // is the layout DJacobian asks for, so the design matrix is the
            // Jacobian -- the same array, not a rearrangement of it
            design = x.getArrayUnsafe();
            target = y.getArrayUnsafe();
        }

        @Override
        public void valueAt(double[] b, double[] r) {
            for (int i = 0; i < m; i++) {
                double v = 0.0;
                for (int j = 0; j < n; j++) {
                    v += design[j * m + i] * b[j];
                }
                r[i] = v - target[i];
            }
        }

        @Override
        public void jacobianAt(double[] b, double[] jac) {
            System.arraycopy(design, 0, jac, 0, m * n);
        }
    }

    /**
     * The check that is not a comparison against MINPACK. A linear least
     * squares problem is a nonlinear one whose Jacobian happens to be constant,
     * so {@link OLS} -- which answers it through the singular values, sharing
     * no code and no idea with a trust region method -- must give the same
     * coefficients.
     * <p>
     * It agrees to 5.6e-16 relative, and the reason it agrees that closely is
     * worth knowing: on a linear model the Gauss Newton step is exact, so the
     * first step lands on the answer and the run is over in three evaluations
     * from any starting point, including one a million away. The tolerance
     * never binds, which is why the loose default reaches the same place the
     * tight setting does.
     */
    @Test
    public void testALinearProblemReproducesOrdinaryLeastSquares() {
        LinearModel model = new LinearModel(40, 4, 20260820L, 0.05);
        LSSummary ols = OLS.estimate(0.05, model.x, model.y);

        double[][] starts = { { 0.0, 0.0, 0.0, 0.0 }, { 100.0, -100.0, 100.0, -100.0 },
                { 1.0e6, 1.0e6, 1.0e6, 1.0e6 } };
        for (int s = 0; s < starts.length; s++) {
            LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(model, starts[s], 40);

            assertTrue("start " + s + ": " + r.status, r.converged);
            assertTrue("start " + s + ": an exact step needs no iterating, " + r.functionEvaluations,
                    r.functionEvaluations <= 5);
            for (int j = 0; j < 4; j++) {
                double reference = ols.getBeta().get(j);
                double scale = Math.max(Math.abs(reference), 1.0);
                assertEquals("start " + s + ": beta[" + j + "]", reference, r.parameters[j], 1.0e-12 * scale);
            }
        }
    }

    /** And without a derivative, where the difference of a linear model is exact. */
    @Test
    public void testTheDerivativeFreePathAlsoReproducesOrdinaryLeastSquares() {
        LinearModel model = new LinearModel(40, 4, 20260820L, 0.05);
        LSSummary ols = OLS.estimate(0.05, model.x, model.y);

        LevenbergMarquardt.Result r = new LevenbergMarquardt().solve((DVectorFunction) model,
                new double[] { 0.0, 0.0, 0.0, 0.0 }, 40);

        assertTrue(r.status.toString(), r.converged);
        for (int j = 0; j < 4; j++) {
            double reference = ols.getBeta().get(j);
            double scale = Math.max(Math.abs(reference), 1.0);
            assertEquals("beta[" + j + "]", reference, r.parameters[j], 1.0e-7 * scale);
        }
    }

    /**
     * The case the class actually exists for, which no linear method answers:
     * three parameters entering an exponential decay {@code A exp(-k t) + C}.
     * <p>
     * The parameters are recovered only as well as the noise allows, so the
     * assertion that carries the weight is the one that does not depend on the
     * noise at all -- the fit must be at least as good as the truth it was
     * generated from. A fit that is worse than the generating parameters has
     * not converged, whatever it reports.
     */
    @Test
    public void testAnExponentialDecayIsRecovered() {
        final int m = 30;
        final double[] truth = { 5.0, 0.7, 1.2 };
        final double[] t = new double[m];
        final double[] y = new double[m];
        Lcg lcg = new Lcg(987654321L);
        for (int i = 0; i < m; i++) {
            t[i] = 5.0 * i / (m - 1);
            y[i] = truth[0] * Math.exp(-truth[1] * t[i]) + truth[2] + 1.0e-3 * lcg.next();
        }

        DiffDVectorFunction decay = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] p, double[] r) {
                for (int i = 0; i < m; i++) {
                    r[i] = p[0] * Math.exp(-p[1] * t[i]) + p[2] - y[i];
                }
            }

            @Override
            public void jacobianAt(double[] p, double[] jac) {
                for (int i = 0; i < m; i++) {
                    double e = Math.exp(-p[1] * t[i]);
                    jac[0 * m + i] = e;
                    jac[1 * m + i] = -p[0] * t[i] * e;
                    jac[2 * m + i] = 1.0;
                }
            }
        };

        LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(decay, new double[] { 1.0, 1.0, 0.0 }, m);

        assertTrue(r.status.toString(), r.converged);
        double[] atTruth = new double[m];
        decay.valueAt(truth, atTruth);
        double truthSumOfSquares = 0.0;
        for (int i = 0; i < m; i++) {
            truthSumOfSquares += atTruth[i] * atTruth[i];
        }
        assertTrue("the fit must beat the parameters the data came from: " + r.sumOfSquares + " against "
                + truthSumOfSquares, r.sumOfSquares <= truthSumOfSquares);

        assertEquals("amplitude", truth[0], r.parameters[0], 1.0e-3);
        assertEquals("rate", truth[1], r.parameters[1], 1.0e-3);
        assertEquals("offset", truth[2], r.parameters[2], 1.0e-3);
    }

    /**
     * Multiplying every residual by a constant multiplies the sum of squares by
     * its square and leaves the minimizer where it was. MINPACK scales the
     * parameters from the columns of the Jacobian, which scale with the
     * residuals, so the scaled problem is the same problem.
     * <p>
     * The same problem, but not the same arithmetic: over sixteen orders of
     * magnitude of scaling the answer moves by at most three units in the last
     * place, never more. Floating-point multiplication is not quite
     * distributive, so the column norms, the factorization and the trust region
     * radius all round a little differently. The bound below is stated at 1e-14
     * relative, which is some forty ulps of room, and still four orders of
     * magnitude tighter than any agreement that could happen by accident.
     */
    @Test
    public void testScalingTheResidualsLeavesTheMinimizerWhereItWas() {
        final LinearModel model = new LinearModel(40, 4, 4242L, 0.05);
        double[] start = { 0.0, 0.0, 0.0, 0.0 };
        LevenbergMarquardt.Result plain = new LevenbergMarquardt().solve(model, start, 40);

        double[] factors = { 1.0e-8, 1.0e-6, 1.0e-3, 1.0e3, 1.0e6, 1.0e8 };
        for (int k = 0; k < factors.length; k++) {
            final double c = factors[k];
            DiffDVectorFunction scaled = new DiffDVectorFunction() {
                @Override
                public void valueAt(double[] b, double[] r) {
                    model.valueAt(b, r);
                    for (int i = 0; i < r.length; i++) {
                        r[i] *= c;
                    }
                }

                @Override
                public void jacobianAt(double[] b, double[] jac) {
                    model.jacobianAt(b, jac);
                    for (int i = 0; i < jac.length; i++) {
                        jac[i] *= c;
                    }
                }
            };

            LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(scaled, start, 40);

            for (int j = 0; j < plain.parameters.length; j++) {
                double scale = Math.max(Math.abs(plain.parameters[j]), 1.0);
                assertEquals("factor " + c + ": the minimizer moved at parameter " + j, plain.parameters[j],
                        r.parameters[j], 1.0e-14 * scale);
            }
            assertEquals("factor " + c + ": the sum of squares scales by the square",
                    plain.sumOfSquares * c * c, r.sumOfSquares, 1.0e-12 * plain.sumOfSquares * c * c);
        }
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
