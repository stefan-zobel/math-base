package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorFunction;
import math.fun.DiffDVectorFunction;
import math.minpack.MghProblems;
import math.optim.LevenbergMarquardt;
import math.solve.NonlinearEquationsSolver.Result;
import math.solve.NonlinearEquationsSolver.Status;

/**
 * The systems are the four of {@link MghProblems} that are square and have a
 * vanishing residual -- so each is a nonlinear system with a published root --
 * plus the Broyden tridiagonal and banded functions, which are scalable and
 * make the {@code O(n^2)} half of Broyden's method visible.
 */
public class NonlinearEquationsSolverTest {

    /** One square system, with what is known about its root. */
    private static final class Sys {

        final String name;
        final int n;
        final double[] start;
        final DiffDVectorFunction f;
        /** The published root, or {@code null} where none is published. */
        final double[] root;
        /**
         * How far two methods may differ on this system. It is a property of
         * the system, not of the solver: Powell singular has a Jacobian that is
         * rank deficient exactly at its root, so the root is only pinned down
         * to about six digits by any method at all, while Rosenbrock is exact.
         */
        final double tolerance;

        Sys(String name, int n, double[] start, DiffDVectorFunction f, double[] root, double tolerance) {
            this.name = name;
            this.n = n;
            this.start = start;
            this.f = f;
            this.root = root;
            this.tolerance = tolerance;
        }

        DVectorFunction withoutTheDerivative() {
            final DiffDVectorFunction inner = f;
            return new DVectorFunction() {
                @Override
                public void valueAt(double[] x, double[] values) {
                    inner.valueAt(x, values);
                }
            };
        }
    }

    /** The largest residual any of these systems needs to be called solved. */
    private static final double RESIDUAL = NonlinearEquationsSolver.DEFAULT_RESIDUAL_TOLERANCE;

    private static NonlinearEquationsSolver solver() {
        return new NonlinearEquationsSolver();
    }

    /** MGH 30, the Broyden tridiagonal function, as a system of {@code n}. */
    private static Sys broydenTridiagonal(final int n) {
        double[] start = new double[n];
        for (int i = 0; i < n; ++i) {
            start[i] = -1.0;
        }
        return new Sys("Broyden tridiagonal n=" + n, n, start, new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                for (int i = 0; i < n; ++i) {
                    double left = (i == 0) ? 0.0 : x[i - 1];
                    double right = (i == n - 1) ? 0.0 : x[i + 1];
                    f[i] = (3.0 - 2.0 * x[i]) * x[i] - left - 2.0 * right + 1.0;
                }
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                for (int c = 0; c < n * n; ++c) {
                    j[c] = 0.0;
                }
                for (int i = 0; i < n; ++i) {
                    j[i * n + i] = 3.0 - 4.0 * x[i];
                    if (i > 0) {
                        j[(i - 1) * n + i] = -1.0;
                    }
                    if (i < n - 1) {
                        j[(i + 1) * n + i] = -2.0;
                    }
                }
            }
        }, null, 1.0e-10);
    }

    /** MGH 31, the Broyden banded function, with {@code ml = 5} and {@code mu = 1}. */
    private static Sys broydenBanded(final int n) {
        double[] start = new double[n];
        for (int i = 0; i < n; ++i) {
            start[i] = -1.0;
        }
        return new Sys("Broyden banded n=" + n, n, start, new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                for (int i = 0; i < n; ++i) {
                    int lo = Math.max(0, i - 5);
                    int hi = Math.min(n - 1, i + 1);
                    double sum = 0.0;
                    for (int j = lo; j <= hi; ++j) {
                        if (j != i) {
                            sum += x[j] * (1.0 + x[j]);
                        }
                    }
                    f[i] = x[i] * (2.0 + 5.0 * x[i] * x[i]) + 1.0 - sum;
                }
            }

            @Override
            public void jacobianAt(double[] x, double[] jac) {
                for (int c = 0; c < n * n; ++c) {
                    jac[c] = 0.0;
                }
                for (int i = 0; i < n; ++i) {
                    int lo = Math.max(0, i - 5);
                    int hi = Math.min(n - 1, i + 1);
                    jac[i * n + i] = 2.0 + 15.0 * x[i] * x[i];
                    for (int j = lo; j <= hi; ++j) {
                        if (j != i) {
                            jac[j * n + i] = -(1.0 + 2.0 * x[j]);
                        }
                    }
                }
            }
        }, null, 1.0e-10);
    }

    /**
     * The systems, freshly built each time because
     * {@link MghProblems#all()} is documented to be worth taking fresh.
     */
    private static Sys[] systems() {
        MghProblems.Problem[] all = MghProblems.all();
        Sys[] out = new Sys[7];
        int k = 0;
        for (int i = 0; i < all.length; ++i) {
            MghProblems.Problem p = all[i];
            if (p.m == p.n && p.hasZeroResidual()) {
                out[k++] = new Sys(p.name, p.n, p.start, p, p.solution, p.solutionTolerance);
            }
        }
        assertEquals("the square zero residual problems of the collection", 4, k);
        out[k++] = broydenTridiagonal(10);
        out[k++] = broydenTridiagonal(50);
        out[k++] = broydenBanded(10);
        return out;
    }

    private long lcg = 987654321987654321L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) / (double) (1L << 53)) * 2.0 - 1.0;
    }

    /**
     * Twelve starting points per system: the published one, six drawn within a
     * quarter of it and five within its own magnitude.
     */
    private double[] startingPoint(Sys sys, int index) {
        if (index == 0) {
            return sys.start.clone();
        }
        double spread = (index < 6) ? 0.25 : 1.0;
        double[] out = new double[sys.n];
        for (int i = 0; i < sys.n; ++i) {
            out[i] = sys.start[i] + spread * next() * Math.max(Math.abs(sys.start[i]), 1.0);
        }
        return out;
    }

    private void seedFor(int system) {
        lcg = 987654321987654321L + 7919L * system;
    }

    private static double worstRelative(double[] a, double[] b) {
        double worst = 0.0;
        for (int i = 0; i < a.length; ++i) {
            worst = Math.max(worst, Math.abs(a[i] - b[i]) / Math.max(Math.abs(b[i]), 1.0));
        }
        return worst;
    }

    private static double infinityNorm(double[] a) {
        double worst = 0.0;
        for (int i = 0; i < a.length; ++i) {
            worst = Math.max(worst, Math.abs(a[i]));
        }
        return worst;
    }

    /**
     * The residual of the defining equation, which is the only thing a root
     * finder can be judged on without a second algorithm: whenever the solver
     * says it solved the system, {@code F} really is below the tolerance there,
     * and the residuals it reports are the ones the function produces at the
     * point it reports.
     */
    @Test
    public void testTheDefiningEquationIsSatisfied() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        double[] values = new double[64];
        for (int s = 0; s < systems.length; ++s) {
            Sys sys = systems[s];
            for (int method = 0; method < 4; ++method) {
                seedFor(s);
                for (int t = 0; t < 12; ++t) {
                    double[] start = startingPoint(sys, t);
                    Result r = run(solver, sys, method, start);
                    String where = sys.name + ", " + nameOf(method) + ", start " + t;
                    if (!r.converged) {
                        assertFalse(where, r.status == Status.SOLVED);
                        continue;
                    }
                    assertSame(where, Status.SOLVED, r.status);
                    assertTrue(where + " : ||F|| = " + r.residualNorm, r.residualNorm <= RESIDUAL);
                    if (values.length < sys.n) {
                        values = new double[sys.n];
                    }
                    sys.f.valueAt(r.solution.clone(), values);
                    for (int i = 0; i < sys.n; ++i) {
                        assertEquals(where + " residual " + i, values[i], r.residuals[i], 0.0);
                    }
                    assertEquals(where, infinityNorm(r.residuals), r.residualNorm, 0.0);
                }
            }
        }
    }

    /**
     * Newton reaches a root from every one of the eighty four starting points,
     * with the derivative and without it. That is the claim the class is worth
     * having at all, so it is asserted as a count rather than left implicit in
     * the test above.
     */
    @Test
    public void testNewtonSolvesEverySystemFromEveryStartingPoint() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        for (int method = 0; method <= 1; ++method) {
            int solved = 0;
            int attempted = 0;
            for (int s = 0; s < systems.length; ++s) {
                seedFor(s);
                for (int t = 0; t < 12; ++t) {
                    ++attempted;
                    if (run(solver, systems[s], method, startingPoint(systems[s], t)).converged) {
                        ++solved;
                    }
                }
            }
            assertEquals(nameOf(method), attempted, solved);
        }
    }

    /**
     * Broyden reaches the same roots, except on Powell badly scaled, whose two
     * equations differ by nine orders of magnitude and whose rank one update
     * therefore ages badly. Asserted as a floor rather than exactly, so that
     * the test says what is expected without pinning a count that a harmless
     * change in the arithmetic could move.
     */
    @Test
    public void testBroydenSolvesEverySystemThatIsNotBadlyScaled() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        for (int method = 2; method <= 3; ++method) {
            for (int s = 0; s < systems.length; ++s) {
                seedFor(s);
                int solved = 0;
                for (int t = 0; t < 12; ++t) {
                    if (run(solver, systems[s], method, startingPoint(systems[s], t)).converged) {
                        ++solved;
                    }
                }
                int expected = "Powell badly scaled".equals(systems[s].name) ? 8 : 12;
                assertTrue(systems[s].name + ", " + nameOf(method) + " : solved " + solved + " of 12",
                        solved >= expected);
            }
        }
    }

    /**
     * Against the roots More, Garbow and Hillstrom publish, to the accuracy
     * they publish them with.
     */
    @Test
    public void testThePublishedRootsAreReached() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        int checked = 0;
        for (int s = 0; s < systems.length; ++s) {
            Sys sys = systems[s];
            if (sys.root == null) {
                continue;
            }
            ++checked;
            for (int method = 0; method < 4; ++method) {
                Result r = run(solver, sys, method, sys.start.clone());
                String where = sys.name + ", " + nameOf(method);
                assertTrue(where, r.converged);
                assertTrue(where + " : off by " + worstRelative(r.solution, sys.root),
                        worstRelative(r.solution, sys.root) <= sys.tolerance);
            }
        }
        assertEquals(4, checked);
    }

    /**
     * The independent algorithm the house rules ask for. A square system with a
     * vanishing residual is a least squares problem whose minimizer is the
     * root, so the Levenberg-Marquardt facade over MINPACK answers the same
     * question by a completely different route -- a trust region on the normal
     * equations rather than a line search on the merit function.
     */
    @Test
    public void testNewtonAgreesWithLevenbergMarquardt() {
        NonlinearEquationsSolver solver = solver();
        LevenbergMarquardt reference = new LevenbergMarquardt(1.0e-14, 1.0e-14, 0.0, 4000, 100.0);
        Sys[] systems = systems();
        int compared = 0;
        for (int s = 0; s < systems.length; ++s) {
            Sys sys = systems[s];
            seedFor(s);
            for (int t = 0; t < 12; ++t) {
                double[] start = startingPoint(sys, t);
                Result mine = solver.newton(sys.f, start.clone());
                LevenbergMarquardt.Result theirs = reference.solve(sys.f, start.clone(), sys.n);
                if (!mine.converged || !theirs.converged || theirs.sumOfSquares > 1.0e-20) {
                    // the reference did not reach a root of its own, so there
                    // is nothing to compare against
                    continue;
                }
                ++compared;
                assertTrue(sys.name + ", start " + t + " : off by " + worstRelative(mine.solution, theirs.parameters),
                        worstRelative(mine.solution, theirs.parameters) <= sys.tolerance);
            }
        }
        assertTrue("the reference reached a root often enough to be one : " + compared, compared >= 70);
    }

    /**
     * All four entry points land on the same root. The approximated Jacobian
     * costs accuracy in the <em>steps</em> and none in the answer, because the
     * fixed point of the iteration is where {@code F} vanishes and does not
     * depend on how the derivative was obtained.
     */
    @Test
    public void testTheFourEntryPointsFindTheSameRoot() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        for (int s = 0; s < systems.length; ++s) {
            Sys sys = systems[s];
            seedFor(s);
            for (int t = 0; t < 12; ++t) {
                double[] start = startingPoint(sys, t);
                Result exact = solver.newton(sys.f, start.clone());
                if (!exact.converged) {
                    continue;
                }
                for (int method = 1; method < 4; ++method) {
                    Result other = run(solver, sys, method, start.clone());
                    if (!other.converged) {
                        continue;
                    }
                    assertTrue(sys.name + ", " + nameOf(method) + ", start " + t + " : off by "
                            + worstRelative(other.solution, exact.solution),
                            worstRelative(other.solution, exact.solution) <= sys.tolerance);
                }
            }
        }
    }

    /**
     * Broyden pays for its rank one update with steps and saves on
     * factorizations: over a whole run it asks for one Jacobian where Newton
     * asks for one per step. On the differenced paths, where a Jacobian is
     * {@code n} further evaluations of {@code F}, that is the difference
     * between the two methods stated in the only currency that matters.
     */
    @Test
    public void testBroydenEvaluatesOneJacobianAndNewtonOnePerStep() {
        NonlinearEquationsSolver solver = solver();
        Sys big = broydenTridiagonal(50);
        Result newton = solver.newton(big.f, big.start.clone());
        Result broyden = solver.broyden(big.f, big.start.clone());
        assertTrue(newton.converged && broyden.converged);
        assertEquals("one Jacobian per step", newton.iterations, newton.jacobianEvaluations);
        assertEquals("one Jacobian for the whole run", 1, broyden.jacobianEvaluations);
        assertEquals("and no restart was needed", 0, broyden.restarts);
        assertTrue("Broyden takes more steps : " + broyden.iterations + " against " + newton.iterations,
                broyden.iterations > newton.iterations);

        Result newtonDifferenced = solver.newton(big.withoutTheDerivative(), big.start.clone());
        Result broydenDifferenced = solver.broyden(big.withoutTheDerivative(), big.start.clone());
        assertTrue(newtonDifferenced.converged && broydenDifferenced.converged);
        assertTrue("and far fewer evaluations of F : " + broydenDifferenced.functionEvaluations + " against "
                + newtonDifferenced.functionEvaluations,
                broydenDifferenced.functionEvaluations < newtonDifferenced.functionEvaluations / 3);
    }

    /**
     * A system with no root is not reported as solved from any starting point,
     * and the residual that comes back with it says so rather than looking
     * small.
     */
    @Test
    public void testASystemWithoutARootIsNotReportedAsSolved() {
        NonlinearEquationsSolver solver = solver();
        DiffDVectorFunction noRoot = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = x[0] * x[0] + 1.0;
                f[1] = x[1];
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                j[0] = 2.0 * x[0];
                j[1] = 0.0;
                j[2] = 0.0;
                j[3] = 1.0;
            }
        };
        double[][] starts = { { 3.0, 4.0 }, { -5.0, 2.0 }, { 0.0, 0.0 }, { 1.0e3, -1.0e3 }, { 0.25, -0.25 } };
        for (int i = 0; i < starts.length; ++i) {
            Result newton = solver.newton(noRoot, starts[i].clone());
            Result broyden = solver.broyden(noRoot, starts[i].clone());
            for (int k = 0; k < 2; ++k) {
                Result r = (k == 0) ? newton : broyden;
                String where = "start (" + starts[i][0] + ", " + starts[i][1] + "), " + ((k == 0) ? "newton"
                        : "broyden");
                assertFalse(where, r.converged);
                assertFalse(where, r.status == Status.SOLVED);
                assertFalse(where, r.status.isSuccess());
                assertTrue(where + " : ||F|| = " + r.residualNorm, r.residualNorm >= 1.0);
            }
        }
    }

    /**
     * At a local minimum of {@code 0.5 F'F} where {@code F} is not zero the
     * Jacobian is necessarily singular -- otherwise {@code J'F = 0} would force
     * {@code F = 0} -- so the two outcomes arrive together and the naming has
     * to separate them. It is the vanishing gradient that is reported, not the
     * singular Jacobian, because that is what tells the caller that no further
     * iteration will help.
     */
    @Test
    public void testAStationaryPointOfTheMeritFunctionIsNamedAsOne() {
        NonlinearEquationsSolver solver = solver();
        // two parallel lines: the merit function has a minimum along the whole
        // line x0 + x1 = 2, where the residuals are 1 and -1
        DiffDVectorFunction parallel = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = x[0] + x[1] - 1.0;
                f[1] = x[0] + x[1] - 3.0;
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                j[0] = 1.0;
                j[1] = 1.0;
                j[2] = 1.0;
                j[3] = 1.0;
            }
        };
        Result newton = solver.newton(parallel, new double[] { 1.0, 1.0 });
        Result broyden = solver.broyden(parallel, new double[] { 1.0, 1.0 });
        assertSame(Status.SPURIOUS_MINIMUM, newton.status);
        assertSame(Status.SPURIOUS_MINIMUM, broyden.status);
        assertFalse(newton.converged);
        assertFalse(broyden.converged);
        assertEquals(1.0, newton.residualNorm, 1.0e-12);
        assertEquals("the point is on the line where the merit function is least",
                2.0, newton.solution[0] + newton.solution[1], 1.0e-9);
    }

    /**
     * A Jacobian that loses rank does not end a Newton run: Powell singular is
     * rank deficient exactly at its root and is solved anyway, and a system
     * whose Jacobian has a zero column throughout is reported honestly instead
     * of throwing. {@link LinearEquationsSolver} rejects such a matrix outright,
     * which is right for a linear system and wrong for a step of a nonlinear
     * one.
     */
    @Test
    public void testASingularJacobianIsNotFatal() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        Sys powell = null;
        for (int s = 0; s < systems.length; ++s) {
            if ("Powell singular".equals(systems[s].name)) {
                powell = systems[s];
            }
        }
        assertNotNull(powell);
        for (int method = 0; method < 4; ++method) {
            Result r = run(solver, powell, method, powell.start.clone());
            assertTrue("Powell singular, " + nameOf(method), r.converged);
        }

        // x1 does not appear in either equation, so the second column of the
        // Jacobian is zero everywhere and no Newton step exists
        DiffDVectorFunction rankOne = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = x[0] - 2.0;
                f[1] = 2.0 * x[0] - 4.0;
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                j[0] = 1.0;
                j[1] = 2.0;
                j[2] = 0.0;
                j[3] = 0.0;
            }
        };
        Result r = solver.newton(rankOne, new double[] { 5.0, 7.0 });
        assertTrue("a root exists along x0 = 2 and the gradient direction finds it", r.converged);
        assertEquals(2.0, r.solution[0], 1.0e-9);
        assertEquals("nothing moved the coordinate the equations do not mention", 7.0, r.solution[1], 0.0);
    }

    /**
     * Multiplying every equation by a constant leaves the root where it was
     * over twelve orders of magnitude. The tolerance is scaled along with the
     * equations, which is the whole content of the statement: the iteration is
     * invariant under the scaling and the stopping rule is not, since it is an
     * absolute bound on a quantity that just moved. Asking for the unscaled
     * tolerance instead is the case
     * {@link #testATooTightToleranceIsReportedRatherThanRoundedAway} covers.
     */
    @Test
    public void testScalingTheEquationsLeavesTheRootWhereItWas() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        double[] factors = { 1.0e-6, 1.0e-3, 3.0, 1.0e3, 1.0e6 };
        for (int s = 0; s < systems.length; ++s) {
            final Sys sys = systems[s];
            Result base = solver.newton(sys.f, sys.start.clone());
            assertTrue(sys.name, base.converged);
            for (int k = 0; k < factors.length; ++k) {
                final double factor = factors[k];
                NonlinearEquationsSolver scaledSolver = new NonlinearEquationsSolver(
                        factor * NonlinearEquationsSolver.DEFAULT_RESIDUAL_TOLERANCE,
                        NonlinearEquationsSolver.DEFAULT_STEP_TOLERANCE,
                        NonlinearEquationsSolver.DEFAULT_MAX_ITERATIONS);
                DiffDVectorFunction scaled = new DiffDVectorFunction() {
                    @Override
                    public void valueAt(double[] x, double[] values) {
                        sys.f.valueAt(x, values);
                        for (int i = 0; i < values.length; ++i) {
                            values[i] *= factor;
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jacobian) {
                        sys.f.jacobianAt(x, jacobian);
                        for (int i = 0; i < jacobian.length; ++i) {
                            jacobian[i] *= factor;
                        }
                    }
                };
                Result r = scaledSolver.newton(scaled, sys.start.clone());
                String where = sys.name + " scaled by " + factor;
                assertTrue(where, r.converged);
                assertTrue(where + " : moved by " + worstRelative(r.solution, base.solution),
                        worstRelative(r.solution, base.solution) <= sys.tolerance);
            }
        }
    }

    /**
     * An absolute tolerance that the scale of the equations puts out of reach
     * is reported as unmet rather than met.
     */
    @Test
    public void testATooTightToleranceIsReportedRatherThanRoundedAway() {
        Sys sys = broydenTridiagonal(10);
        Result reachable = new NonlinearEquationsSolver(1.0e-12, NonlinearEquationsSolver.DEFAULT_STEP_TOLERANCE,
                200).newton(sys.f, sys.start.clone());
        assertTrue(reachable.converged);
        Result unreachable = new NonlinearEquationsSolver(1.0e-20,
                NonlinearEquationsSolver.DEFAULT_STEP_TOLERANCE, 200).newton(sys.f, sys.start.clone());
        assertFalse(unreachable.converged);
        assertFalse(unreachable.status == Status.SOLVED);
        assertTrue("but it is the same point : " + worstRelative(unreachable.solution, reachable.solution),
                worstRelative(unreachable.solution, reachable.solution) <= 1.0e-12);
    }

    /** An exhausted budget is reported as an exhausted budget. */
    @Test
    public void testAnExhaustedBudgetKeepsThePointItReached() {
        Sys sys = broydenTridiagonal(10);
        NonlinearEquationsSolver stingy = new NonlinearEquationsSolver(
                NonlinearEquationsSolver.DEFAULT_RESIDUAL_TOLERANCE,
                NonlinearEquationsSolver.DEFAULT_STEP_TOLERANCE, 2);
        Result r = stingy.newton(sys.f, sys.start.clone());
        assertSame(Status.ITERATION_LIMIT, r.status);
        assertFalse(r.converged);
        assertEquals(2, r.iterations);
        assertTrue("and it got closer than it started : " + r.residualNorm, r.residualNorm < 1.0);
    }

    /** A starting point that is already a root is recognized without a step. */
    @Test
    public void testAStartingPointThatIsAlreadyARootCostsNoIteration() {
        NonlinearEquationsSolver solver = solver();
        Sys sys = broydenTridiagonal(10);
        double[] root = solver.newton(sys.f, sys.start.clone()).solution;
        for (int method = 0; method < 4; ++method) {
            Result r = run(solver, sys, method, root.clone());
            assertSame(nameOf(method), Status.SOLVED, r.status);
            assertEquals(nameOf(method), 0, r.iterations);
            assertEquals(nameOf(method), 0, r.jacobianEvaluations);
            assertEquals(nameOf(method), 1, r.functionEvaluations);
        }
    }

    /** Residuals that are not finite where the search starts end it at once. */
    @Test
    public void testResidualsThatAreNotFiniteAtTheStartAreReported() {
        NonlinearEquationsSolver solver = solver();
        DiffDVectorFunction unbounded = new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] x, double[] f) {
                f[0] = 1.0 / (x[0] - 1.0);
                f[1] = x[1];
            }

            @Override
            public void jacobianAt(double[] x, double[] j) {
                j[0] = -1.0 / ((x[0] - 1.0) * (x[0] - 1.0));
                j[1] = 0.0;
                j[2] = 0.0;
                j[3] = 1.0;
            }
        };
        Result r = solver.newton(unbounded, new double[] { 1.0, 0.0 });
        assertSame(Status.NOT_FINITE, r.status);
        assertFalse(r.converged);
        assertEquals(0, r.iterations);
    }

    /** The array the caller hands in is never written to. */
    @Test
    public void testTheStartingPointIsNotModified() {
        NonlinearEquationsSolver solver = solver();
        Sys[] systems = systems();
        for (int s = 0; s < systems.length; ++s) {
            Sys sys = systems[s];
            for (int method = 0; method < 4; ++method) {
                double[] start = sys.start.clone();
                double[] copy = start.clone();
                run(solver, sys, method, start);
                for (int i = 0; i < start.length; ++i) {
                    assertEquals(sys.name + ", " + nameOf(method) + ", component " + i, copy[i], start[i], 0.0);
                }
            }
        }
    }

    @Test
    public void testArgumentValidation() {
        NonlinearEquationsSolver solver = solver();
        DiffDVectorFunction f = broydenTridiagonal(2).f;
        rejects(solver, null, new double[] { 1.0, 1.0 });
        rejects(solver, f, null);
        rejects(solver, f, new double[0]);
        rejects(solver, f, new double[] { 1.0, Double.NaN });
        rejects(solver, f, new double[] { Double.POSITIVE_INFINITY, 1.0 });

        double[][] badConfigurations = { { 0.0, 1.0e-14, 200.0 }, { -1.0, 1.0e-14, 200.0 },
                { Double.POSITIVE_INFINITY, 1.0e-14, 200.0 }, { Double.NaN, 1.0e-14, 200.0 },
                { 1.0e-12, 0.0, 200.0 }, { 1.0e-12, -1.0, 200.0 }, { 1.0e-12, Double.NaN, 200.0 },
                { 1.0e-12, 1.0e-14, 0.0 }, { 1.0e-12, 1.0e-14, -3.0 } };
        for (int i = 0; i < badConfigurations.length; ++i) {
            double[] c = badConfigurations[i];
            try {
                new NonlinearEquationsSolver(c[0], c[1], (int) c[2]);
                fail("configuration " + i + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertNotNull(expected.getMessage());
            }
        }
    }

    private static void rejects(NonlinearEquationsSolver solver, DiffDVectorFunction f, double[] start) {
        try {
            solver.newton(f, start);
            fail("newton accepted it");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            solver.broyden(f, start);
            fail("broyden accepted it");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    private static Result run(NonlinearEquationsSolver solver, Sys sys, int method, double[] start) {
        switch (method) {
        case 0:
            return solver.newton(sys.f, start);
        case 1:
            return solver.newton(sys.withoutTheDerivative(), start);
        case 2:
            return solver.broyden(sys.f, start);
        default:
            return solver.broyden(sys.withoutTheDerivative(), start);
        }
    }

    private static String nameOf(int method) {
        return new String[] { "newton", "newton with a differenced Jacobian", "broyden",
                "broyden with a differenced Jacobian" }[method];
    }
}
