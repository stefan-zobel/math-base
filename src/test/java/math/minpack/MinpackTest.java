package math.minpack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.linalg.VectorOps;
import math.minpack.MghProblems.Problem;

/**
 * Tests for the MINPACK port, which had none at all: 3590 lines, no test and no
 * caller. The oracle is external -- {@link MghProblems}, the collection of More,
 * Garbow and Hillstrom -- and the port is exercised through its own FORTRAN
 * calling convention here, which is what {@code math.optim.LevenbergMarquardt}
 * exists to hide.
 * <p>
 * The published numbers are not the load-bearing assertion, though. Two checks
 * come first and depend on nothing but the problem definitions: that their
 * analytic Jacobians agree with finite differences, so that a wrong derivation
 * cannot be mistaken for a defect in the port, and that the gradient of the sum
 * of squares vanishes at every reported solution, which is the defining
 * equation of a least squares minimum. The published minima then catch what
 * those two cannot -- a mistyped data constant, which is self-consistent and
 * would pass both.
 */
public class MinpackTest {

    /** Recommended lower limit for the MINPACK tolerances. */
    private static final double SQRT_EPS = Math.sqrt(Math.ulp(1.0));

    /** What a run of the port returns, unpacked back to zero-based arrays. */
    private static final class Outcome {

        double[] x;
        double[] residuals;
        double sumOfSquares;
        int info;
        int functionEvaluations;
        int jacobianEvaluations;
    }

    /**
     * Bridges a {@link Problem} to the one-based FORTRAN calling convention,
     * including the turn from the column-major Jacobian the library uses to the
     * row-major two-dimensional one MINPACK wants.
     */
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

    /**
     * Runs {@code lmder_f77}, the full driver -- {@code lmder1_f77} keeps its
     * evaluation counts to itself.
     */
    private static Outcome runLmder(Problem p, double tol) {
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

        return unpack(p, x, fvec, info, nfev, njev);
    }

    /** Runs {@code lmdif_f77}, which builds the Jacobian by forward differences. */
    private static Outcome runLmdif(Problem p, double tol) {
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

        return unpack(p, x, fvec, info, nfev, new int[2]);
    }

    private static Outcome unpack(Problem p, double[] x, double[] fvec, int[] info, int[] nfev, int[] njev) {
        Outcome o = new Outcome();
        o.x = new double[p.n];
        System.arraycopy(x, 1, o.x, 0, p.n);
        o.residuals = new double[p.m];
        System.arraycopy(fvec, 1, o.residuals, 0, p.m);
        double ssq = 0.0;
        for (int i = 0; i < p.m; i++) {
            ssq += o.residuals[i] * o.residuals[i];
        }
        o.sumOfSquares = ssq;
        o.info = info[1];
        o.functionEvaluations = nfev[1];
        o.jacobianEvaluations = njev[1];
        return o;
    }

    /**
     * The success codes. One and two are the two tolerances met, three is both,
     * and four is a residual orthogonal to the Jacobian -- a stationary point.
     * <p>
     * Eight belongs here too, and only looks like a failure: MINPACK sets it
     * when the scaled gradient norm is down at machine precision, which is a
     * stationary point by any measure, and files it under "gtol is too small"
     * merely because a smaller one was asked for. {@code lmder1_f77} makes the
     * same judgement in the port itself, remapping eight to four before it
     * returns. The distinction matters for the facade -- read as a plain
     * number, eight is larger than every success code and smaller than nothing.
     */
    private static boolean isSuccess(int info) {
        return (info >= 1 && info <= 4) || info == 8;
    }

    /**
     * Every analytic Jacobian of the collection, against central differences.
     * This runs first on purpose: a wrong derivation there would send the solver
     * somewhere odd and look exactly like a defect in the port. It pins the
     * column-major layout at the same time, since a transposed Jacobian would
     * disagree with the differences everywhere off the diagonal.
     */
    @Test
    public void testAnalyticJacobiansAgreeWithFiniteDifferences() {
        double h = Math.cbrt(Math.ulp(1.0));
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            double[] analytic = new double[p.m * p.n];
            p.jacobianAt(p.start, analytic);

            double[] plus = new double[p.m];
            double[] minus = new double[p.m];
            double[] x = p.start.clone();
            for (int j = 0; j < p.n; j++) {
                double step = h * Math.max(Math.abs(p.start[j]), 1.0);
                x[j] = p.start[j] + step;
                p.valueAt(x, plus);
                x[j] = p.start[j] - step;
                p.valueAt(x, minus);
                x[j] = p.start[j];

                for (int i = 0; i < p.m; i++) {
                    double numeric = (plus[i] - minus[i]) / (2.0 * step);
                    double scale = Math.max(Math.abs(numeric), 1.0);
                    assertEquals(p.name + ": d f[" + i + "] / d x[" + j + "]", numeric, analytic[j * p.m + i],
                            1.0e-5 * scale);
                }
            }
        }
    }

    /**
     * The problems whose residuals vanish at the solution are the ones where
     * both the minimizer and the minimum are known exactly.
     */
    @Test
    public void testZeroResidualProblemsReachTheirKnownSolution() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            if (!p.hasZeroResidual()) {
                continue;
            }
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": info = " + o.info, isSuccess(o.info));
            assertEquals(p.name + ": sum of squares", 0.0, o.sumOfSquares, p.minimumTolerance);
            for (int j = 0; j < p.n; j++) {
                double scale = Math.max(Math.abs(p.solution[j]), 1.0);
                assertEquals(p.name + ": x[" + j + "]", p.solution[j], o.x[j], p.solutionTolerance * scale);
            }
        }
    }

    /**
     * The rest are data fitting problems with a residual that does not vanish.
     * Only the minimal sum of squares is published for them, to six digits, and
     * the tolerance here is that of the published number rather than of the
     * algorithm.
     */
    @Test
    public void testNonZeroResidualProblemsReachThePublishedMinimum() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            if (p.hasZeroResidual()) {
                continue;
            }
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": info = " + o.info, isSuccess(o.info));
            assertEquals(p.name + ": sum of squares", p.minimum, o.sumOfSquares, p.minimumTolerance);
        }
    }

    /**
     * The check that needs no published number: at a least squares minimum the
     * residual is orthogonal to every column of the Jacobian, which is the
     * defining equation. Stated as a cosine so that it is scale free and one
     * bound fits problems whose residuals differ by ten orders of magnitude.
     * <p>
     * Only the problems with a residual that does not vanish are asked. Where
     * it does, the cosine is the angle of a vector of length 1e-16 and carries
     * no information; for those the sum of squares itself is the stronger
     * statement, and the test above makes it.
     * <p>
     * The bound has to be stated relative to the tolerance that was asked for,
     * because MINPACK stops on the value and the step, never on the gradient
     * unless {@code gtol} says so, and these runs pass {@code gtol == 0}. The
     * gradient is therefore only as small as the other two force it to be, and
     * the rate differs by regime:
     * <ul>
     * <li>Kowalik and Osborne has a minimal sum of squares of 3e-4. The Gauss
     * Newton model is good there, convergence is quadratic, and the cosine
     * falls <em>linearly</em> with the tolerance -- 5.8e-6 at 1.5e-8, 5.5e-8 at
     * 1e-12, a factor of about 400 throughout.</li>
     * <li>Brown and Dennis has a minimal sum of squares of 85822. Only the
     * value test binds there, and the value is quadratic at the minimum, so a
     * relative tolerance of {@code tol} in it pins the point to
     * {@code sqrt(tol)} and the gradient with it -- 1.6e-4 at 1.5e-8, 9.3e-7 at
     * 1e-12, both within a factor of two of {@code sqrt(tol)}.</li>
     * </ul>
     * The second is the weaker of the two and sets the bound. Running at two
     * tolerances three orders of magnitude apart asserts the coupling rather
     * than a single number, which is the part that would be hard to reproduce
     * by accident.
     */
    @Test
    public void testTheGradientVanishesAtEverySolution() {
        double[] tolerances = { SQRT_EPS, 1.0e-12 };
        Problem[] all = MghProblems.all();
        for (int t = 0; t < tolerances.length; t++) {
            double tol = tolerances[t];
            for (int k = 0; k < all.length; k++) {
                Problem p = all[k];
                if (p.hasZeroResidual()) {
                    continue;
                }
                Outcome o = runLmder(p, tol);

                double[] jac = new double[p.m * p.n];
                p.jacobianAt(o.x, jac);
                double residualNorm = VectorOps.twoNorm(o.residuals);
                for (int j = 0; j < p.n; j++) {
                    double dot = 0.0;
                    double columnNorm = 0.0;
                    for (int i = 0; i < p.m; i++) {
                        dot += jac[j * p.m + i] * o.residuals[i];
                        columnNorm += jac[j * p.m + i] * jac[j * p.m + i];
                    }
                    columnNorm = Math.sqrt(columnNorm);
                    if (columnNorm == 0.0) {
                        continue;
                    }
                    double cosine = Math.abs(dot) / (columnNorm * residualNorm);
                    assertTrue(p.name + " at tol " + tol + ": column " + j
                            + " is not orthogonal to the residual, cosine = " + cosine,
                            cosine < 10.0 * Math.sqrt(tol));
                }
            }
        }
    }

    /**
     * {@code fvec} comes back as an output argument, so it can disagree with
     * the returned point without anything complaining. It must not.
     */
    @Test
    public void testReturnedResidualsBelongToTheReturnedPoint() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome o = runLmder(p, SQRT_EPS);

            double[] expected = new double[p.m];
            p.valueAt(o.x, expected);
            for (int i = 0; i < p.m; i++) {
                assertEquals(p.name + ": residual " + i, expected[i], o.residuals[i], 0.0);
            }
        }
    }

    /**
     * The solver never returns a point worse than the one it was given. Trivial
     * to state, and the kind of thing a mistranslated index breaks first.
     */
    @Test
    public void testNoProblemEndsWorseThanItStarted() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": " + o.sumOfSquares + " > " + p.sumOfSquaresAt(p.start),
                    o.sumOfSquares <= p.sumOfSquaresAt(p.start));
        }
    }

    /**
     * {@code lmdif} builds the Jacobian by forward differences instead of
     * asking for it. That is a different algorithm, so what is asserted here is
     * the minimum reached and the price paid for it, never the path or the
     * status code -- see the test below for the problem where those two part
     * company. The budgets each problem carries are sized for this path rather
     * than for the analytic one, which is the cheaper of the two throughout:
     * Brown and Dennis, for instance, converges in 254 evaluations with a
     * Jacobian and 1229 without.
     */
    @Test
    public void testLmdifReachesTheSameMinimaWithoutAnAnalyticJacobian() {
        Problem[] all = MghProblems.all();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome withJacobian = runLmder(p, SQRT_EPS);
            Outcome withoutJacobian = runLmdif(p, SQRT_EPS);

            double scale = Math.max(Math.abs(withJacobian.sumOfSquares), 1.0e-8);
            assertEquals(p.name + ": sum of squares", withJacobian.sumOfSquares, withoutJacobian.sumOfSquares,
                    1.0e-4 * scale);
            assertTrue(p.name + ": not asking for the Jacobian cannot be cheaper, "
                    + withoutJacobian.functionEvaluations + " against " + withJacobian.functionEvaluations,
                    withoutJacobian.functionEvaluations >= withJacobian.functionEvaluations);
            assertTrue(withJacobian.jacobianEvaluations >= 0);
        }
    }

    /**
     * Powell singular has a rank deficient Jacobian at its minimum, and it is
     * the one problem here that separates <em>having</em> the answer from
     * <em>being able to say so</em>. With an analytic Jacobian the orthogonality
     * test fires and the run ends at once. Without one, the noise of the
     * forward difference keeps that test from firing, and MINPACK falls back on
     * the relative step test -- which converges linearly on a singular problem
     * and needs about 98000 evaluations to be satisfied, having reached a sum
     * of squares of 2e-41 within the first thousand.
     * <p>
     * So {@code info = 5} here reports an exhausted budget on a point that is
     * exact to eleven digits. Any facade over this port that maps five to
     * "failed" without also handing back the point would be throwing away a
     * correct answer.
     */
    @Test
    public void testTheDerivativeFreePathPaysForARankDeficientJacobian() {
        Problem powellSingular = byName("Powell singular");

        Outcome withJacobian = runLmder(powellSingular, SQRT_EPS);
        Outcome withoutJacobian = runLmdif(powellSingular, SQRT_EPS);

        assertEquals("with an analytic Jacobian the orthogonality test fires", 8, withJacobian.info);
        assertTrue("and it is cheap: " + withJacobian.functionEvaluations, withJacobian.functionEvaluations < 100);

        assertEquals("without one the budget runs out", 5, withoutJacobian.info);
        assertTrue("yet the point is already exact, sum of squares = " + withoutJacobian.sumOfSquares,
                withoutJacobian.sumOfSquares < 1.0e-30);
        for (int j = 0; j < powellSingular.n; j++) {
            assertEquals("x[" + j + "]", 0.0, withoutJacobian.x[j], 1.0e-8);
        }
    }

    /**
     * {@code enorm_f77} carries its own three-accumulator scaling, and the
     * expected values here are stated in closed form rather than taken from
     * another implementation.
     * <p>
     * This test used to assert the opposite of what it asserts now, and it is
     * the reason the change was noticed. {@link VectorOps#twoNorm} accumulated
     * the squares directly, so it overflowed above about {@code 1e154} and
     * underflowed below about {@code 1e-162}; three of the six vectors below
     * were listed as out of its range and it was asserted to fail on them,
     * "so the day {@code twoNorm} is rewritten in the BLAS {@code dnrm2} style
     * this test notices". That day came. The two are now independent
     * implementations of the same quantity -- three accumulators split by
     * magnitude here, one scaling pass by a power of two there -- and they are
     * held to agreeing on every vector.
     */
    @Test
    public void testEnormAgreesWithTwoNormAtEveryMagnitude() {
        double root3 = Math.sqrt(3.0);
        double[][] vectors = { { 3.0, 4.0 }, { 1.0e300, 1.0e300, 1.0e300 }, { 1.0e-300, 1.0e-300, 1.0e-300 },
                { 0.0, 0.0, 0.0 }, { 1.0e300, 1.0e-300, 1.0 }, { -3.0, 0.0, 4.0 } };
        double[] expected = { 5.0, root3 * 1.0e300, root3 * 1.0e-300, 0.0, 1.0e300, 5.0 };

        for (int k = 0; k < vectors.length; k++) {
            double[] v = vectors[k];
            double[] oneBased = new double[v.length + 1];
            System.arraycopy(v, 0, oneBased, 1, v.length);

            double enorm = Minpack_f77.enorm_f77(v.length, oneBased);
            assertEquals("vector " + k, expected[k], enorm, 1.0e-14 * Math.max(expected[k], Double.MIN_NORMAL));

            double twoNorm = VectorOps.twoNorm(v);
            assertEquals("vector " + k + ": the two implementations must agree", enorm, twoNorm,
                    1.0e-14 * Math.max(enorm, Double.MIN_NORMAL));
        }
    }

    /**
     * The trap this port sets for a caller who does not read the FORTRAN
     * comments: {@code info = 5} means the evaluation budget ran out. It is a
     * larger number than the four success codes and looks like a better one.
     */
    @Test
    public void testBudgetExhaustionIsReportedAsFiveRatherThanSuccess() {
        Problem meyer = byName("Meyer");

        int m = meyer.m;
        int n = meyer.n;
        double[] x = new double[n + 1];
        System.arraycopy(meyer.start, 0, x, 1, n);
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        Minpack_f77.lmder_f77(new Bridge(meyer), m, n, x, fvec, fjac, SQRT_EPS, SQRT_EPS, 0.0, 3, diag, 1, 100.0, 0,
                info, nfev, njev, ipvt, qtf);

        assertEquals("info", 5, info[1]);
        assertFalse("5 is not one of the success codes", isSuccess(info[1]));
    }

    /** Improper input is reported as {@code info = 0}, not thrown. */
    @Test
    public void testImproperInputIsReportedAsZero() {
        Problem p = MghProblems.all()[0];
        int m = p.m;
        int n = p.n;
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        // a negative tolerance
        int[] info = new int[2];
        Minpack_f77.lmder_f77(new Bridge(p), m, n, new double[n + 1], fvec, fjac, -1.0, SQRT_EPS, 0.0, 100, diag, 1,
                100.0, 0, info, nfev, njev, ipvt, qtf);
        assertEquals("negative ftol", 0, info[1]);

        // more variables than functions
        info = new int[2];
        Minpack_f77.lmder_f77(new Bridge(p), 1, n, new double[n + 1], fvec, fjac, SQRT_EPS, SQRT_EPS, 0.0, 100, diag,
                1, 100.0, 0, info, nfev, njev, ipvt, qtf);
        assertEquals("m < n", 0, info[1]);
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
