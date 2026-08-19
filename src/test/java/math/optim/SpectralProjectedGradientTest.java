package math.optim;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.MathConsts;
import math.fun.DMultiFunctionEval;
import math.fun.DiffDMultiFunction;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;

public class SpectralProjectedGradientTest {

    /** No constraint at all, for comparison against the unconstrained methods. */
    private static final Projection UNCONSTRAINED = x -> {
        // the whole space, so there is nothing to do
    };

    /**
     * {@code 0.5 * ||x - a||^2}. Its minimum over any convex set is the
     * projection of {@code a} onto that set, so this one problem turns every
     * projection into an oracle for the minimizer.
     */
    private static final class Distance implements DiffDMultiFunction {

        private final double[] a;

        Distance(double[] a) {
            this.a = a;
        }

        @Override
        public double apply(double[] x) {
            double s = 0.0;
            for (int i = 0; i < x.length; i++) {
                double d = x[i] - a[i];
                s += d * d;
            }
            return 0.5 * s;
        }

        @Override
        public void derivativeAt(double[] x, double[] grad) {
            for (int i = 0; i < x.length; i++) {
                grad[i] = x[i] - a[i];
            }
        }
    }

    /** {@code 0.5 * sum scale[i] * (x[i] - centre[i])^2}, condition max/min scale. */
    private static final class Conditioned implements DiffDMultiFunction {

        private final double[] scale;
        private final double[] centre;

        Conditioned(double[] scale, double[] centre) {
            this.scale = scale;
            this.centre = centre;
        }

        @Override
        public double apply(double[] x) {
            double s = 0.0;
            for (int i = 0; i < x.length; i++) {
                double d = x[i] - centre[i];
                s += scale[i] * d * d;
            }
            return 0.5 * s;
        }

        @Override
        public void derivativeAt(double[] x, double[] grad) {
            for (int i = 0; i < x.length; i++) {
                grad[i] = scale[i] * (x[i] - centre[i]);
            }
        }
    }

    /** Extended Rosenbrock, minimum 0 at (1, ..., 1). */
    private static final class Rosenbrock implements DiffDMultiFunction {

        @Override
        public double apply(double[] x) {
            double s = 0.0;
            for (int i = 0; i + 1 < x.length; i++) {
                double a = x[i + 1] - x[i] * x[i];
                double b = 1.0 - x[i];
                s += 100.0 * a * a + b * b;
            }
            return s;
        }

        @Override
        public void derivativeAt(double[] x, double[] grad) {
            Arrays.fill(grad, 0.0);
            for (int i = 0; i + 1 < x.length; i++) {
                double a = x[i + 1] - x[i] * x[i];
                grad[i] += -400.0 * x[i] * a - 2.0 * (1.0 - x[i]);
                grad[i + 1] += 200.0 * a;
            }
        }
    }

    // -----------------------------------------------------------------
    // the projection as an oracle for the minimizer
    // -----------------------------------------------------------------

    /**
     * Minimizing the distance to {@code a} over a set has the projection of
     * {@code a} as its exact answer, so the minimizer has to reproduce the
     * projection it is driven by, for every one of the five sets.
     */
    @Test
    public void testTheProjectionProblemReproducesTheProjection() {
        double[] a = { -7.0, 2.5, 9.0, 0.4 };
        Projection[] sets = { Projection.box(fill(4, -1.0), fill(4, 2.0)), Projection.nonNegative(4),
                Projection.simplex(4, 1.0), Projection.euclideanBall(4, 2.0), Projection.l1Ball(4, 2.0) };
        String[] names = { "box", "nonNegative", "simplex", "euclideanBall", "l1Ball" };
        SpectralProjectedGradient spg = new SpectralProjectedGradient();
        for (int k = 0; k < sets.length; k++) {
            double[] expected = a.clone();
            sets[k].projectInto(expected);
            SpectralProjectedGradient.Result r = spg.minimize(new Distance(a), fill(4, 5.0), sets[k]);
            assertTrue(names[k] + " did not converge : " + r.status, r.converged);
            // the Hessian is the identity, so one step of the right length
            // solves this exactly; the first spectral step is 1/||pg|| rather
            // than 1, which costs one more
            assertTrue(names[k] + " took " + r.iterations + " iterations", r.iterations <= 2);
            assertArrayEquals(names[k], expected, r.point, 1.0e-15);
        }
    }

    /**
     * The property the whole class exists for: a coordinate whose optimum lies
     * outside its bound comes back sitting <i>on</i> the bound, bit for bit,
     * where a reparametrization with {@code exp} or a logistic could only
     * approach it.
     */
    @Test
    public void testActiveBoxCoordinatesLandExactlyOnTheBound() {
        int n = 5;
        double[] lower = new double[n];
        double[] upper = fill(n, 1.0);
        double[] centre = new double[n];
        double[] scale = new double[n];
        for (int i = 0; i < n; i++) {
            // every unconstrained optimum is far outside the box
            centre[i] = (i % 2 == 0) ? 5.0 : -5.0;
            scale[i] = Math.pow(10.0, i);
        }
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                .minimize(new Conditioned(scale, centre), new double[n], lower, upper);
        assertTrue(r.status.toString(), r.converged);
        for (int i = 0; i < n; i++) {
            double bound = (centre[i] > 0.0) ? upper[i] : lower[i];
            assertEquals("coordinate " + i + " is not exactly on its bound", bound, r.point[i], 0.0);
        }
    }

    // -----------------------------------------------------------------
    // agreement with the unconstrained optimizers
    // -----------------------------------------------------------------

    @Test
    public void testUnconstrainedMatchesTheClosedFormMinimum() {
        double[] scale = { 1.0, 0.8, 0.6, 0.4 };
        double[] centre = { 1.0, -2.0, 3.0, -4.0 };
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                .minimize(new Conditioned(scale, centre), new double[] { 9.0, 9.0, 9.0, 9.0 }, UNCONSTRAINED);
        assertTrue(r.status.toString(), r.converged);
        // the tolerance is one on the gradient, which is
        // scale[i] * (x[i] - centre[i]) here, so the accuracy of a coordinate
        // is that tolerance divided by its own scale and nothing better
        for (int i = 0; i < centre.length; i++) {
            assertEquals("coordinate " + i, centre[i], r.point[i],
                    r.projectedGradientNorm / scale[i] + 1.0e-15);
        }
        assertTrue("value " + r.value, r.value <= 1.0e-14);
    }

    /**
     * A conjugate gradient method answering the same unconstrained question by
     * an unrelated route.
     */
    @Test
    public void testUnconstrainedAgreesWithCGOptimizer() {
        double[] scale = { 1.0, 0.8, 0.6, 0.4 };
        double[] centre = { 1.0, -2.0, 3.0, -4.0 };
        Conditioned f = new Conditioned(scale, centre);
        double[] start = { 9.0, 9.0, 9.0, 9.0 };
        SpectralProjectedGradient.Result mine = new SpectralProjectedGradient().minimize(f, start, UNCONSTRAINED);
        DMultiFunctionEval theirs = CGOptimizer.minimize(f, start.clone());
        assertTrue(theirs.converged);
        assertArrayEquals(theirs.point, mine.point, 1.0e-6);
    }

    // -----------------------------------------------------------------
    // invariants that hold whatever the outcome
    // -----------------------------------------------------------------

    @Test
    public void testEveryReturnedPointIsFeasible() {
        Projection[] sets = { Projection.box(fill(4, -1.0), fill(4, 2.0)), Projection.nonNegative(4),
                Projection.simplex(4, 1.0), Projection.euclideanBall(4, 2.0), Projection.l1Ball(4, 2.0) };
        String[] names = { "box", "nonNegative", "simplex", "euclideanBall", "l1Ball" };
        DiffDMultiFunction[] problems = { new Rosenbrock(), new Distance(new double[] { 8.0, -8.0, 0.5, 3.0 }),
                new Conditioned(new double[] { 1.0, 1.0e2, 1.0e4, 1.0e6 }, new double[] { 0.3, -9.0, 4.0, 0.0 }) };
        SpectralProjectedGradient spg = new SpectralProjectedGradient(1.0e-8, 1.0e-12, 60, 10);
        for (int k = 0; k < sets.length; k++) {
            for (DiffDMultiFunction f : problems) {
                SpectralProjectedGradient.Result r = spg.minimize(f, new double[] { 0.4, 0.1, 0.2, 0.3 }, sets[k]);
                assertNotNull(r.status);
                assertFeasible(names[k], sets[k], r.point);
            }
        }
    }

    @Test
    public void testTheStartIsNotModifiedAndAnInfeasibleStartIsProjected() {
        double[] start = { 5.0, -5.0, 5.0 };
        double[] untouched = start.clone();
        Projection p = Projection.box(fill(3, -1.0), fill(3, 1.0));
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                .minimize(new Distance(new double[] { 0.25, 0.5, 0.75 }), start, p);
        assertArrayEquals("the starting point was modified", untouched, start, 0.0);
        assertFeasible("box", p, r.point);
        assertArrayEquals(new double[] { 0.25, 0.5, 0.75 }, r.point, 1.0e-15);
    }

    @Test
    public void testStartingAtTheConstrainedSolutionConvergesWithoutAnIteration() {
        double[] a = { 9.0, -9.0 };
        Projection p = Projection.box(fill(2, -1.0), fill(2, 1.0));
        double[] solution = a.clone();
        p.projectInto(solution);
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient().minimize(new Distance(a), solution, p);
        assertEquals(SpectralProjectedGradient.Status.KKT_TOLERANCE_REACHED, r.status);
        assertEquals(0, r.iterations);
        assertEquals(1, r.functionEvaluations);
        assertEquals(1, r.gradientEvaluations);
        assertEquals(0.0, r.projectedGradientNorm, 0.0);
        assertArrayEquals(solution, r.point, 0.0);
    }

    @Test
    public void testAnExhaustedBudgetIsReportedAndThePointIsStillFeasible() {
        Projection p = Projection.box(fill(2, -5.0), fill(2, 5.0));
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient(1.0e-10, 1.0e-14, 3, 10)
                .minimize(new Rosenbrock(), new double[] { -1.2, 1.0 }, p);
        assertEquals(SpectralProjectedGradient.Status.TOO_MANY_ITERATIONS, r.status);
        assertTrue(!r.converged);
        assertEquals(3, r.iterations);
        assertFeasible("box", p, r.point);
    }

    /**
     * The whole point of the window: with {@code memory == 1} the search is
     * monotone, so a larger budget can never give a worse value, and with
     * {@code memory == 10} the value is free to rise for a while. Both halves
     * are observable from outside by running the same problem with a growing
     * cap, which is the only handle a caller has on the inner iteration.
     */
    @Test
    public void testTheWindowLetsTheValueRiseAndAWindowOfOneDoesNot() {
        Projection p = Projection.box(fill(2, -5.0), fill(2, 5.0));
        double[] start = { -1.2, 1.0 };

        double previous = Double.POSITIVE_INFINITY;
        for (int budget = 1; budget <= 40; budget++) {
            SpectralProjectedGradient.Result r = new SpectralProjectedGradient(1.0e-10, 1.0e-14, budget, 1)
                    .minimize(new Rosenbrock(), start, p);
            assertTrue("memory 1, budget " + budget + " gave " + r.value + " after " + previous,
                    r.value <= previous);
            previous = r.value;
        }

        boolean rose = false;
        previous = Double.POSITIVE_INFINITY;
        for (int budget = 1; budget <= 40; budget++) {
            SpectralProjectedGradient.Result r = new SpectralProjectedGradient(1.0e-10, 1.0e-14, budget, 10)
                    .minimize(new Rosenbrock(), start, p);
            if (r.value > previous) {
                rose = true;
            }
            previous = r.value;
        }
        assertTrue("the window of ten never accepted an increase, so it is not in effect", rose);
    }

    /**
     * What the window buys, measured rather than cited. The same code with
     * {@code memory == 1} needs several times as many iterations on the same
     * problem, and on Rosenbrock in ten dimensions it is the difference
     * between reaching the tolerance and not.
     */
    @Test
    public void testTheNonmonotoneWindowMeetsTheCostItWasAddedFor() {
        Projection p = Projection.box(fill(2, -5.0), fill(2, 5.0));
        double[] start = { -1.2, 1.0 };
        SpectralProjectedGradient.Result monotone = new SpectralProjectedGradient(
                Math.sqrt(MathConsts.MACH_EPS_DBL), 1.0e-12, 0, 1).minimize(new Rosenbrock(), start, p);
        SpectralProjectedGradient.Result windowed = new SpectralProjectedGradient(
                Math.sqrt(MathConsts.MACH_EPS_DBL), 1.0e-12, 0, 10).minimize(new Rosenbrock(), start, p);
        assertTrue("monotone : " + monotone.status, monotone.converged);
        assertTrue("windowed : " + windowed.status, windowed.converged);
        assertTrue("the window did not pay: " + windowed.iterations + " against " + monotone.iterations,
                2 * windowed.iterations < monotone.iterations);
    }

    @Test
    public void testANonFiniteGradientIsReported() {
        DiffDMultiFunction broken = new DiffDMultiFunction() {
            @Override
            public double apply(double[] x) {
                return x[0] * x[0];
            }

            @Override
            public void derivativeAt(double[] x, double[] grad) {
                grad[0] = Double.POSITIVE_INFINITY;
            }
        };
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient().minimize(broken, new double[] { 2.0 },
                UNCONSTRAINED);
        assertEquals(SpectralProjectedGradient.Status.NOT_FINITE, r.status);
        assertTrue(!r.converged);
        assertEquals(0, r.iterations);
    }

    // -----------------------------------------------------------------
    // what the spectral step buys
    // -----------------------------------------------------------------

    /**
     * The Barzilai-Borwein step is the whole reason this method is usable, and
     * the numbers here are the measurement rather than a citation. With a
     * fixed step of one, the very same code reached <b>none</b> of these six
     * problems: the three quadratics ran out after 1600 iterations with the
     * projected gradient still at {@code 4.1e-7}, {@code 6.9e-5} and
     * {@code 5.7e-3}, and the two Rosenbrocks ran out after 1200 at
     * {@code 7.7e-3} and {@code 1.1e-3}. With the spectral step and the
     * window the measured counts are 10, 10, 8, 69 and 57 iterations; the caps
     * below leave room for two to three times that.
     */
    @Test
    public void testTheSpectralStepMeetsTheCostItWasAddedFor() {
        SpectralProjectedGradient spg = new SpectralProjectedGradient();
        int[] caps = { 40, 40, 40 };
        int c = 0;
        for (int exponent = 2; exponent <= 6; exponent += 2) {
            int m = 6;
            double[] scale = new double[m];
            double[] centre = new double[m];
            for (int i = 0; i < m; i++) {
                scale[i] = Math.pow(10.0, (exponent * i) / (double) (m - 1));
                centre[i] = (i % 2 == 0) ? 0.5 * (i + 1) / (double) m : 4.0;
            }
            SpectralProjectedGradient.Result r = spg.minimize(new Conditioned(scale, centre), new double[m],
                    fill(m, -1.0), fill(m, 1.0));
            assertTrue("condition 1e" + exponent + " : " + r.status, r.converged);
            assertTrue("condition 1e" + exponent + " needed " + r.iterations + " iterations",
                    r.iterations <= caps[c]);
            c++;
        }

        SpectralProjectedGradient.Result interior = spg.minimize(new Rosenbrock(), new double[] { -1.2, 1.0 },
                fill(2, -5.0), fill(2, 5.0));
        assertTrue("interior : " + interior.status, interior.converged);
        assertTrue("interior needed " + interior.iterations, interior.iterations <= 200);

        SpectralProjectedGradient.Result cutOff = spg.minimize(new Rosenbrock(), new double[] { -1.2, 1.0 },
                new double[] { -5.0, -5.0 }, new double[] { 0.5, 5.0 });
        assertTrue("cut off : " + cutOff.status, cutOff.converged);
        assertTrue("cut off needed " + cutOff.iterations, cutOff.iterations <= 150);
        // the minimum of Rosenbrock over x0 <= 0.5 sits on that bound, and
        // there it is 0.25 at (0.5, 0.25)
        assertEquals(0.5, cutOff.point[0], 0.0);
        assertEquals(0.25, cutOff.point[1], 1.0e-9);
        assertEquals(0.25, cutOff.value, 1.0e-12);
    }

    // -----------------------------------------------------------------
    // independent oracles: everything above says the class agrees with
    // itself, these say it is right
    // -----------------------------------------------------------------

    /**
     * For a strictly convex quadratic over a box in four variables the exact
     * answer can be had by brute force: every coordinate is at its lower
     * bound, at its upper bound or free, so enumerating all {@code 3^n}
     * patterns, solving the reduced system for the free block and keeping the
     * one that is feasible with the right multiplier signs gives the solution
     * by a route that has nothing in common with a projected gradient method.
     */
    @Test
    public void testAgreesWithTheActiveSetSolutionOfABoxConstrainedQuadratic() {
        Lcg rng = new Lcg(20260820L);
        int n = 4;
        double[] lower = fill(n, -1.0);
        double[] upper = fill(n, 1.0);
        for (int trial = 0; trial < 25; trial++) {
            double[][] q = randomPositiveDefinite(n, rng);
            double[] b = new double[n];
            for (int i = 0; i < n; i++) {
                // large enough that the solution is on the boundary often
                b[i] = 4.0 * rng.centred();
            }
            double[] expected = activeSetSolution(q, b, lower, upper);
            assertNotNull("trial " + trial + ": no active set pattern was admissible", expected);
            SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                    .minimize(new Quadratic(q, b), new double[n], lower, upper);
            assertTrue("trial " + trial + " : " + r.status, r.converged);
            // Q = M'M + I has its smallest eigenvalue at or above one, so the
            // distance to the exact minimizer is bounded by the projected
            // gradient itself, and by sqrt(n) times it per coordinate once the
            // infinity norm is turned into a Euclidean one
            double bound = Math.sqrt(n) * r.projectedGradientNorm + 1.0e-12;
            assertArrayEquals("trial " + trial, expected, r.point, bound);
        }
    }

    /**
     * Non-negative least squares on a design whose ordinary solution is
     * already positive: the constraint is inactive at the optimum, so the two
     * have to agree, and {@link math.linalg.OLS} gets there through a singular
     * value decomposition that shares no line of code with this class.
     */
    @Test
    public void testNonNegativeLeastSquaresAgreesWithOlsWhereTheConstraintIsInactive() {
        Lcg rng = new Lcg(4711L);
        int rows = 40;
        int columns = 4;
        double[] truth = { 2.0, 1.5, 0.75, 3.0 };
        double[][] x = new double[rows][columns];
        double[] y = new double[rows];
        DMatrix design = new DMatrix(rows, columns);
        DMatrix response = new DMatrix(rows, 1);
        for (int i = 0; i < rows; i++) {
            double value = 0.0;
            for (int j = 0; j < columns; j++) {
                x[i][j] = (j == 0) ? 1.0 : rng.centred();
                design.set(i, j, x[i][j]);
                value += truth[j] * x[i][j];
            }
            y[i] = value + 0.01 * rng.centred();
            response.set(i, 0, y[i]);
        }
        LSSummary ols = OLS.estimate(0.05, design, response);
        double[] reference = new double[columns];
        for (int j = 0; j < columns; j++) {
            reference[j] = ols.getBeta().get(j);
            assertTrue("the ordinary solution has to be positive for this test to mean anything : "
                    + reference[j], reference[j] > 0.1);
        }
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                .minimize(new SumOfSquares(x, y), new double[columns], Projection.nonNegative(columns));
        assertTrue(r.status.toString(), r.converged);
        assertArrayEquals(reference, r.point, 1.0e-7);
    }

    /**
     * The other half: a coefficient whose ordinary estimate is negative is
     * driven to exactly zero rather than to something small, which is the
     * point of solving the constrained problem instead of clipping afterwards.
     */
    @Test
    public void testNonNegativeLeastSquaresPinsANegativeCoefficientToExactlyZero() {
        Lcg rng = new Lcg(881L);
        int rows = 40;
        int columns = 3;
        double[] truth = { 2.0, -3.0, 1.0 };
        double[][] x = new double[rows][columns];
        double[] y = new double[rows];
        DMatrix design = new DMatrix(rows, columns);
        DMatrix response = new DMatrix(rows, 1);
        for (int i = 0; i < rows; i++) {
            double value = 0.0;
            for (int j = 0; j < columns; j++) {
                x[i][j] = (j == 0) ? 1.0 : rng.centred();
                design.set(i, j, x[i][j]);
                value += truth[j] * x[i][j];
            }
            y[i] = value + 0.01 * rng.centred();
            response.set(i, 0, y[i]);
        }
        assertTrue("the ordinary estimate has to be negative for this test to mean anything",
                OLS.estimate(0.05, design, response).getBeta().get(1) < -1.0);
        SpectralProjectedGradient.Result r = new SpectralProjectedGradient()
                .minimize(new SumOfSquares(x, y), new double[columns], Projection.nonNegative(columns));
        assertTrue(r.status.toString(), r.converged);
        assertEquals("the constrained coefficient is not exactly zero", 0.0, r.point[1], 0.0);
    }

    /**
     * On the simplex the optimality condition is that the gradient is the same
     * on the support and no smaller off it. Checked directly, without a
     * reference implementation, because there is none in this repository.
     */
    @Test
    public void testTheSimplexSolutionSatisfiesTheKktConditions() {
        Lcg rng = new Lcg(1234L);
        int n = 6;
        for (int trial = 0; trial < 10; trial++) {
            double[][] q = randomPositiveDefinite(n, rng);
            double[] b = new double[n];
            for (int i = 0; i < n; i++) {
                b[i] = 3.0 * rng.centred();
            }
            Quadratic f = new Quadratic(q, b);
            SpectralProjectedGradient.Result r = new SpectralProjectedGradient().minimize(f,
                    fill(n, 1.0 / n), Projection.simplex(n, 1.0));
            assertTrue("trial " + trial + " : " + r.status, r.converged);

            double total = 0.0;
            for (int i = 0; i < n; i++) {
                assertTrue("negative weight " + r.point[i], r.point[i] >= 0.0);
                total += r.point[i];
            }
            assertEquals(1.0, total, 1.0e-12);

            double[] gradient = new double[n];
            f.derivativeAt(r.point, gradient);
            double multiplier = Double.POSITIVE_INFINITY;
            for (int i = 0; i < n; i++) {
                if (r.point[i] > 0.0 && gradient[i] < multiplier) {
                    multiplier = gradient[i];
                }
            }
            for (int i = 0; i < n; i++) {
                if (r.point[i] > 1.0e-10) {
                    assertEquals("the gradient is not level on the support", multiplier, gradient[i], 1.0e-7);
                } else {
                    assertTrue("a zero weight has a gradient below the multiplier",
                            gradient[i] >= multiplier - 1.0e-7);
                }
            }
        }
    }

    /**
     * The badly scaled family {@code s * x0^2 + x1^2}, as
     * {@code CGOptimizerTest} runs it, but in a box that cuts the first
     * coordinate off from its optimum. The bound is met exactly no matter how
     * violent the scaling.
     */
    @Test
    public void testTheBadlyScaledFamilyInABox() {
        for (double s : new double[] { 1.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8 }) {
            DiffDMultiFunction f = new Conditioned(new double[] { 2.0 * s, 2.0 }, new double[] { 0.0, 0.0 });
            SpectralProjectedGradient.Result r = new SpectralProjectedGradient().minimize(f,
                    new double[] { 2.0, 0.9 }, new double[] { 0.5, -1.0 }, new double[] { 2.0, 1.0 });
            assertTrue("s = " + s + " : " + r.status, r.converged);
            assertEquals("s = " + s, 0.5, r.point[0], 0.0);
            assertEquals("s = " + s, 0.0, r.point[1], 1.0e-8);
        }
    }

    // -----------------------------------------------------------------
    // validation
    // -----------------------------------------------------------------

    @Test
    public void testConstructorValidation() {
        expectIae("negative kktTolerance", () -> new SpectralProjectedGradient(-1.0e-8, 1.0e-12, 0, 10));
        expectIae("NaN kktTolerance", () -> new SpectralProjectedGradient(Double.NaN, 1.0e-12, 0, 10));
        expectIae("negative stepTolerance", () -> new SpectralProjectedGradient(1.0e-8, -1.0e-12, 0, 10));
        expectIae("NaN stepTolerance", () -> new SpectralProjectedGradient(1.0e-8, Double.NaN, 0, 10));
        expectIae("memory 0", () -> new SpectralProjectedGradient(1.0e-8, 1.0e-12, 0, 0));
        expectIae("negative memory", () -> new SpectralProjectedGradient(1.0e-8, 1.0e-12, 0, -1));
    }

    @Test
    public void testMinimizeValidation() {
        SpectralProjectedGradient spg = new SpectralProjectedGradient();
        Distance f = new Distance(new double[] { 0.0, 0.0 });
        double[] start = { 1.0, 1.0 };
        expectIae("null f", () -> spg.minimize(null, start, UNCONSTRAINED));
        expectIae("null start", () -> spg.minimize(f, null, UNCONSTRAINED));
        expectIae("null projection", () -> spg.minimize(f, start, (Projection) null));
        expectIae("empty start", () -> spg.minimize(f, new double[0], UNCONSTRAINED));
        expectIae("NaN start", () -> spg.minimize(f, new double[] { Double.NaN, 0.0 }, UNCONSTRAINED));
        expectIae("infinite start", () -> spg.minimize(f, new double[] { Double.POSITIVE_INFINITY, 0.0 },
                UNCONSTRAINED));
    }

    @Test
    public void testANonFiniteValueAtTheStartIsRejected() {
        DiffDMultiFunction broken = new DiffDMultiFunction() {
            @Override
            public double apply(double[] x) {
                return Double.NaN;
            }

            @Override
            public void derivativeAt(double[] x, double[] grad) {
                grad[0] = 0.0;
            }
        };
        expectIae("NaN value at the start",
                () -> new SpectralProjectedGradient().minimize(broken, new double[] { 1.0 }, UNCONSTRAINED));
    }

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /** Deterministic LCG, so a failure is always reproducible. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            state = seed;
        }

        /** Uniform in {@code [0, 1)}. */
        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }

        /** Uniform in {@code [-1, 1)}. */
        double centred() {
            return 2.0 * next() - 1.0;
        }
    }

    /** {@code 0.5 * x'Qx + b'x} with {@code Q} symmetric positive definite. */
    private static final class Quadratic implements DiffDMultiFunction {

        private final double[][] q;
        private final double[] b;

        Quadratic(double[][] q, double[] b) {
            this.q = q;
            this.b = b;
        }

        @Override
        public double apply(double[] x) {
            double s = 0.0;
            for (int i = 0; i < x.length; i++) {
                double t = 0.0;
                for (int j = 0; j < x.length; j++) {
                    t += q[i][j] * x[j];
                }
                s += x[i] * (0.5 * t + b[i]);
            }
            return s;
        }

        @Override
        public void derivativeAt(double[] x, double[] grad) {
            for (int i = 0; i < x.length; i++) {
                double t = b[i];
                for (int j = 0; j < x.length; j++) {
                    t += q[i][j] * x[j];
                }
                grad[i] = t;
            }
        }
    }

    /** {@code ||X beta - y||^2}, the least squares objective in its own right. */
    private static final class SumOfSquares implements DiffDMultiFunction {

        private final double[][] x;
        private final double[] y;

        SumOfSquares(double[][] x, double[] y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public double apply(double[] beta) {
            double s = 0.0;
            for (int i = 0; i < x.length; i++) {
                double r = -y[i];
                for (int j = 0; j < beta.length; j++) {
                    r += x[i][j] * beta[j];
                }
                s += r * r;
            }
            return s;
        }

        @Override
        public void derivativeAt(double[] beta, double[] grad) {
            Arrays.fill(grad, 0.0);
            for (int i = 0; i < x.length; i++) {
                double r = -y[i];
                for (int j = 0; j < beta.length; j++) {
                    r += x[i][j] * beta[j];
                }
                for (int j = 0; j < beta.length; j++) {
                    grad[j] += 2.0 * r * x[i][j];
                }
            }
        }
    }

    private static double[][] randomPositiveDefinite(int n, Lcg rng) {
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i][j] = rng.centred();
            }
        }
        double[][] q = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = (i == j) ? 1.0 : 0.0;
                for (int k = 0; k < n; k++) {
                    s += m[k][i] * m[k][j];
                }
                q[i][j] = s;
            }
        }
        return q;
    }

    /**
     * The exact minimizer of {@code 0.5 x'Qx + b'x} over a box, by trying every
     * assignment of every coordinate to lower bound, free or upper bound. The
     * one that solves its reduced system inside the box and has the right
     * multiplier signs is the solution, and for a strictly convex problem it is
     * the only one.
     */
    private static double[] activeSetSolution(double[][] q, double[] b, double[] lower, double[] upper) {
        int n = b.length;
        int patterns = 1;
        for (int i = 0; i < n; i++) {
            patterns *= 3;
        }
        for (int pattern = 0; pattern < patterns; pattern++) {
            int[] state = new int[n];
            int rest = pattern;
            for (int i = 0; i < n; i++) {
                state[i] = rest % 3;
                rest /= 3;
            }
            double[] x = new double[n];
            int freeCount = 0;
            for (int i = 0; i < n; i++) {
                if (state[i] == 0) {
                    x[i] = lower[i];
                } else if (state[i] == 2) {
                    x[i] = upper[i];
                } else {
                    freeCount++;
                }
            }
            if (freeCount > 0) {
                int[] index = new int[freeCount];
                int k = 0;
                for (int i = 0; i < n; i++) {
                    if (state[i] == 1) {
                        index[k++] = i;
                    }
                }
                double[][] a = new double[freeCount][freeCount];
                double[] rhs = new double[freeCount];
                for (int r = 0; r < freeCount; r++) {
                    for (int c = 0; c < freeCount; c++) {
                        a[r][c] = q[index[r]][index[c]];
                    }
                    double s = b[index[r]];
                    for (int i = 0; i < n; i++) {
                        if (state[i] != 1) {
                            s += q[index[r]][i] * x[i];
                        }
                    }
                    rhs[r] = -s;
                }
                double[] solution = solveLinear(a, rhs);
                for (int r = 0; r < freeCount; r++) {
                    x[index[r]] = solution[r];
                }
            }
            boolean admissible = true;
            for (int i = 0; i < n && admissible; i++) {
                if (state[i] == 1 && (x[i] < lower[i] - 1.0e-10 || x[i] > upper[i] + 1.0e-10)) {
                    admissible = false;
                }
            }
            if (!admissible) {
                continue;
            }
            for (int i = 0; i < n && admissible; i++) {
                double g = b[i];
                for (int j = 0; j < n; j++) {
                    g += q[i][j] * x[j];
                }
                // at a lower bound the gradient has to push out of the box,
                // at an upper bound the other way, and a free coordinate is
                // stationary
                if (state[i] == 0 && g < -1.0e-9) {
                    admissible = false;
                } else if (state[i] == 2 && g > 1.0e-9) {
                    admissible = false;
                } else if (state[i] == 1 && Math.abs(g) > 1.0e-9) {
                    admissible = false;
                }
            }
            if (admissible) {
                return x;
            }
        }
        return null;
    }

    /** Gaussian elimination with partial pivoting, for the small blocks above. */
    private static double[] solveLinear(double[][] matrix, double[] rightHandSide) {
        int m = rightHandSide.length;
        double[][] a = new double[m][];
        for (int i = 0; i < m; i++) {
            a[i] = matrix[i].clone();
        }
        double[] b = rightHandSide.clone();
        for (int k = 0; k < m; k++) {
            int pivot = k;
            for (int i = k + 1; i < m; i++) {
                if (Math.abs(a[i][k]) > Math.abs(a[pivot][k])) {
                    pivot = i;
                }
            }
            double[] rowSwap = a[k];
            a[k] = a[pivot];
            a[pivot] = rowSwap;
            double valueSwap = b[k];
            b[k] = b[pivot];
            b[pivot] = valueSwap;
            for (int i = k + 1; i < m; i++) {
                double factor = a[i][k] / a[k][k];
                for (int j = k; j < m; j++) {
                    a[i][j] -= factor * a[k][j];
                }
                b[i] -= factor * b[k];
            }
        }
        double[] x = new double[m];
        for (int i = m - 1; i >= 0; i--) {
            double s = b[i];
            for (int j = i + 1; j < m; j++) {
                s -= a[i][j] * x[j];
            }
            x[i] = s / a[i][i];
        }
        return x;
    }

    /**
     * A point is in the set exactly when the projection leaves it where it is,
     * which is one check for all five sets and needs no knowledge of any of
     * them.
     */
    private static void assertFeasible(String name, Projection p, double[] point) {
        double[] projected = point.clone();
        p.projectInto(projected);
        for (int i = 0; i < point.length; i++) {
            double scale = 1.0 + Math.abs(point[i]);
            assertEquals(name + ": coordinate " + i + " is outside the set", projected[i], point[i],
                    1.0e-12 * scale);
        }
    }

    private static double[] fill(int n, double value) {
        double[] a = new double[n];
        Arrays.fill(a, value);
        return a;
    }

    private static void expectIae(String what, Runnable r) {
        try {
            r.run();
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
