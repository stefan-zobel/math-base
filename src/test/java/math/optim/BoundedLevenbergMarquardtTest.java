package math.optim;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.fun.DVectorFunction;
import math.fun.DiffDMultiFunction;
import math.fun.DiffDVectorFunction;
import math.linalg.BoundedLeastSquares;
import math.linalg.DMatrix;
import math.linalg.OLS;
import math.minpack.MghProblems;
import math.minpack.MghProblems.Problem;

/**
 * Tests for {@link BoundedLevenbergMarquardt}. The class exists so that the
 * model is never evaluated outside its own box, so the load-bearing assertion
 * is exactly that: every problem is wrapped in a function that records where it
 * was called, and no call may land outside. The rest rests on three oracles
 * that share no code with it -- {@link BoundedLeastSquares} where the model is
 * linear, {@link SpectralProjectedGradient} where it is not, and the exhaustive
 * enumeration of active sets for the small cases.
 */
public class BoundedLevenbergMarquardtTest {

    /** The tolerance MINPACK recommends as its lower limit. */
    private static final double SQRT_EPS = Math.sqrt(Math.ulp(1.0));

    // ----------------------------------------------------------------- box --

    @Test
    public void testEveryPointTheModelSeesIsInsideTheBox() {
        int calls = 0;
        for (Problem p : MghProblems.all()) {
            for (int variant = 0; variant < 4; variant++) {
                double[] lower = lowerFor(p, variant);
                double[] upper = upperFor(p, variant);
                Caged caged = new Caged(p, lower, upper);
                solver(5 * p.maxEvaluations).solve(caged, p.start, p.m, lower, upper);
                assertEquals(p.name + ", variant " + variant + ": the model was evaluated outside its box",
                        0, caged.violations);
                calls += caged.calls;
            }
        }
        assertTrue("the wrapper should have seen a substantial number of calls", calls > 1000);
    }

    @Test
    public void testAParameterAtABoundEqualsThatBoundExactly() {
        int held = 0;
        for (Problem p : MghProblems.all()) {
            for (int variant = 0; variant < 4; variant++) {
                double[] lower = lowerFor(p, variant);
                double[] upper = upperFor(p, variant);
                BoundedLevenbergMarquardt.Result r = solver(5 * p.maxEvaluations).solve(p, p.start, p.m, lower,
                        upper);
                for (int j = 0; j < p.n; j++) {
                    if (r.activeSet[j] == BoundedLeastSquares.Bound.AT_LOWER) {
                        assertEquals(p.name + ", parameter " + j, lower[j], r.parameters[j], 0.0);
                        held++;
                    } else if (r.activeSet[j] == BoundedLeastSquares.Bound.AT_UPPER) {
                        assertEquals(p.name + ", parameter " + j, upper[j], r.parameters[j], 0.0);
                        held++;
                    }
                }
            }
        }
        assertTrue("no parameter ever ended up at a bound, so the test proved nothing", held > 20);
    }

    @Test
    public void testTheActiveSetAndTheParametersAgreeInBothDirections() {
        for (Problem p : MghProblems.all()) {
            for (int variant = 0; variant < 4; variant++) {
                double[] lower = lowerFor(p, variant);
                double[] upper = upperFor(p, variant);
                BoundedLevenbergMarquardt.Result r = solver(5 * p.maxEvaluations).solve(p, p.start, p.m, lower,
                        upper);
                int bound = 0;
                for (int j = 0; j < p.n; j++) {
                    boolean atLower = r.parameters[j] == lower[j];
                    boolean atUpper = r.parameters[j] == upper[j];
                    if (r.activeSet[j] == BoundedLeastSquares.Bound.FREE) {
                        assertFalse(p.name + ": free but at its lower bound", atLower);
                        assertFalse(p.name + ": free but at its upper bound", atUpper);
                        assertTrue(p.name + ": free but outside the box",
                                r.parameters[j] > lower[j] && r.parameters[j] < upper[j]);
                    } else {
                        assertTrue(p.name + ": held but not at a bound", atLower || atUpper);
                        bound++;
                    }
                }
                assertEquals(p.name + ": atBound disagrees with the active set", bound, r.atBound);
            }
        }
    }

    @Test
    public void testAStartingPointOutsideTheBoxIsMovedOntoIt() {
        double[] t = grid(40, 0.0, 5.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        Decay model = new Decay(t, y);
        double[] lower = { 0.5, 0.25, -1.0 };
        double[] upper = { 10.0, 10.0, 1.0 };
        // every component of the start violates the box
        double[] start = { -100.0, 1.0e6, 50.0 };
        Caged caged = new Caged(model, lower, upper);
        BoundedLevenbergMarquardt.Result r = solver(2000).solve(caged, start, t.length, lower, upper);
        assertEquals("clamping the start is not allowed to escape the box", 0, caged.violations);
        assertTrue("the fit should have converged from the clamped start", r.converged);
        assertEquals("amplitude", 3.0, r.parameters[0], 1.0e-6);
        assertEquals("rate", 1.5, r.parameters[1], 1.0e-6);
        assertEquals("offset", 0.5, r.parameters[2], 1.0e-6);
    }

    @Test
    public void testAPinnedParameterKeepsItsValueExactly() {
        double[] t = grid(40, 0.0, 5.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        Decay model = new Decay(t, y);
        double pinned = 0.75;
        double[] lower = { 0.0, 0.01, pinned };
        double[] upper = { 100.0, 100.0, pinned };
        BoundedLevenbergMarquardt.Result r = solver(2000).solve(model, new double[] { 1.0, 0.5, 0.0 }, t.length,
                lower, upper);
        assertEquals("a pinned parameter must come back untouched", pinned, r.parameters[2], 0.0);
        assertTrue("a pinned parameter is not free",
                r.activeSet[2] != BoundedLeastSquares.Bound.FREE);
        assertEquals("it consumed no degree of freedom", t.length - 2, r.degreesOfFreedom);
    }

    // ------------------------------------------------------------- oracles --

    @Test
    public void testALinearModelReproducesBoundedLeastSquares() {
        // a least squares problem whose Jacobian happens to be constant IS a
        // bounded linear one, so the exact and finite solver has to agree
        double worstParameter = 0.0;
        double worstSumOfSquares = 0.0;
        int compared = 0;
        for (int trial = 0; trial < 24; trial++) {
            int m = 30 + 7 * trial;
            int n = 3 + trial % 8;
            Lcg lcg = new Lcg(20260820L + trial);
            double[] a = lcg.array(m * n);
            double[] b = lcg.array(m);
            double[] lower = new double[n];
            double[] upper = new double[n];
            for (int j = 0; j < n; j++) {
                int variant = trial % 3;
                if (variant == 0) {
                    lower[j] = 0.0;
                    upper[j] = Double.POSITIVE_INFINITY;
                } else if (variant == 1) {
                    lower[j] = -0.2;
                    upper[j] = 0.2;
                } else {
                    lower[j] = (j % 3 == 0) ? 0.0 : Double.NEGATIVE_INFINITY;
                    upper[j] = (j % 3 == 1) ? 0.1 : Double.POSITIVE_INFINITY;
                }
            }
            BoundedLeastSquares.Result exact = BoundedLeastSquares.bounded(matrix(m, n, a), matrix(m, 1, b), lower,
                    upper);
            BoundedLevenbergMarquardt tight = new BoundedLevenbergMarquardt(1.0e-12, 1.0e-12, 0.0, 500, 1.0e-3);
            BoundedLevenbergMarquardt.Result r = tight.solve(new Linear(a, b, m, n), new double[n], m, lower,
                    upper);

            double exactSumOfSquares = 0.0;
            for (int i = 0; i < m; i++) {
                exactSumOfSquares += exact.residuals[i] * exact.residuals[i];
            }
            for (int j = 0; j < n; j++) {
                assertEquals("trial " + trial + ", parameter " + j + ": the active sets must agree exactly",
                        exact.activeSet[j], r.activeSet[j]);
                worstParameter = Math.max(worstParameter, Math.abs(r.parameters[j] - exact.solution[j]));
                compared++;
            }
            worstSumOfSquares = Math.max(worstSumOfSquares,
                    Math.abs(r.sumOfSquares - exactSumOfSquares) / Math.max(1.0, exactSumOfSquares));
        }
        assertTrue("nothing was compared", compared > 100);
        // measured: 1.6e-9 and 8.1e-16 over these designs
        assertTrue("worst parameter difference was " + worstParameter, worstParameter < 1.0e-7);
        assertTrue("worst relative difference in the sum of squares was " + worstSumOfSquares,
                worstSumOfSquares < 1.0e-11);
    }

    @Test
    public void testTheCollectionReachesThePublishedMinima() {
        for (Problem p : MghProblems.all()) {
            BoundedLevenbergMarquardt.Result r = solver(p.maxEvaluations).solve(p, p.start, p.m,
                    fill(p.n, Double.NEGATIVE_INFINITY), fill(p.n, Double.POSITIVE_INFINITY));
            assertTrue(p.name + " did not converge : " + r.status, r.converged);
            assertEquals(p.name, p.minimum, r.sumOfSquares,
                    p.minimumTolerance * Math.max(1.0, Math.abs(p.minimum)));
            assertEquals(p.name + ": no bound may bind when there is none", 0, r.atBound);
        }
    }

    @Test
    public void testItAgreesWithLevenbergMarquardtWhereNoBoundBinds() {
        // not bit for bit: this is a different trust region rule on the same
        // problem, and both stop at a tolerance rather than at a point.
        // Measured worst disagreement over the collection: 3.4e-5 relative, on
        // Brown and Dennis, whose minimum is flat enough that MINPACK's own
        // documentation warns about it
        double worst = 0.0;
        for (Problem p : MghProblems.all()) {
            LevenbergMarquardt reference = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, p.maxEvaluations, 100.0);
            LevenbergMarquardt.Result expected = reference.solve(p, p.start, p.m);
            BoundedLevenbergMarquardt.Result actual = solver(p.maxEvaluations).solve(p, p.start, p.m,
                    fill(p.n, Double.NEGATIVE_INFINITY), fill(p.n, Double.POSITIVE_INFINITY));
            for (int j = 0; j < p.n; j++) {
                double relative = Math.abs(actual.parameters[j] - expected.parameters[j])
                        / Math.max(1.0, Math.abs(expected.parameters[j]));
                worst = Math.max(worst, relative);
                assertEquals(p.name + ", parameter " + j, expected.parameters[j], actual.parameters[j],
                        1.0e-3 * Math.max(1.0, Math.abs(expected.parameters[j])));
            }
        }
        assertTrue("the agreement got worse than it was measured to be : " + worst, worst < 1.0e-3);
    }

    @Test
    public void testAHeldParameterIsNeverPulledBackIntoTheBox() {
        // the multiplier sign, which is what a wrong active set would break.
        // It has to hold whatever tolerance stopped the search, unlike
        // stationarity of the free parameters
        int held = 0;
        for (Problem p : MghProblems.all()) {
            for (int variant = 0; variant < 4; variant++) {
                double[] lower = lowerFor(p, variant);
                double[] upper = upperFor(p, variant);
                BoundedLevenbergMarquardt.Result r = solver(5 * p.maxEvaluations).solve(p, p.start, p.m, lower,
                        upper);
                if (!isFinite(r.residuals)) {
                    continue;
                }
                double[] jacobian = new double[p.m * p.n];
                p.jacobianAt(r.parameters, jacobian);
                double residualNorm = norm(r.residuals);
                for (int j = 0; j < p.n; j++) {
                    if (r.activeSet[j] == BoundedLeastSquares.Bound.FREE) {
                        continue;
                    }
                    held++;
                    double gradient = 0.0;
                    double columnNorm = 0.0;
                    for (int i = 0; i < p.m; i++) {
                        gradient += jacobian[j * p.m + i] * r.residuals[i];
                        columnNorm += jacobian[j * p.m + i] * jacobian[j * p.m + i];
                    }
                    columnNorm = Math.sqrt(columnNorm);
                    double inward = (r.activeSet[j] == BoundedLeastSquares.Bound.AT_LOWER)
                            ? Math.max(0.0, -gradient) : Math.max(0.0, gradient);
                    double floor = columnNorm * residualNorm;
                    // measured worst over these runs: 9.3e-45 relative
                    assertTrue(p.name + ", variant " + variant + ", parameter " + j
                            + " is held at a bound its gradient wants to leave : " + inward,
                            inward <= 1.0e-8 * floor);
                }
            }
        }
        assertTrue("no parameter was ever held, so the test proved nothing", held > 20);
    }

    @Test
    public void testItAgreesWithSpectralProjectedGradient() {
        // the same bounded problem by a first order method over a projection,
        // sharing no line of code with the active set inside the step
        double worst = 0.0;
        int compared = 0;
        for (Problem p : MghProblems.all()) {
            double[] lower = fill(p.n, -1.0);
            double[] upper = fill(p.n, 1.0);
            BoundedLevenbergMarquardt.Result mine = solver(5 * p.maxEvaluations).solve(p, p.start, p.m, lower,
                    upper);
            SpectralProjectedGradient spg = new SpectralProjectedGradient(1.0e-10, 1.0e-14, 20000, 10);
            SpectralProjectedGradient.Result reference = spg.minimize(new Squares(p, p.m, p.n), p.start.clone(),
                    lower, upper);
            if (!reference.converged) {
                continue;
            }
            compared++;
            double relative = Math.abs(mine.sumOfSquares - reference.value) / Math.max(1.0, Math.abs(reference.value));
            worst = Math.max(worst, relative);
            assertEquals(p.name, reference.value, mine.sumOfSquares,
                    1.0e-5 * Math.max(1.0, Math.abs(reference.value)));
        }
        assertTrue("SPG converged on too few problems to prove anything", compared >= 8);
        assertTrue("the agreement got worse than it was measured to be : " + worst, worst < 1.0e-5);
    }

    @Test
    public void testAgreesWithTheExhaustiveActiveSetSolution() {
        // every parameter is at its lower bound, at its upper bound or free;
        // enumerating all 3^n patterns and keeping the feasible one with the
        // right multiplier signs is the answer by pure combinatorics
        int checked = 0;
        for (int trial = 0; trial < 20; trial++) {
            int n = 4;
            int m = 12;
            Lcg lcg = new Lcg(9876543L + trial);
            double[] a = lcg.array(m * n);
            double[] b = lcg.array(m);
            double[] lower = new double[n];
            double[] upper = new double[n];
            for (int j = 0; j < n; j++) {
                if (trial % 3 == 0) {
                    lower[j] = 0.0;
                    upper[j] = Double.POSITIVE_INFINITY;
                } else if (trial % 3 == 1) {
                    lower[j] = -0.15;
                    upper[j] = 0.15;
                } else {
                    lower[j] = (j == 0) ? 0.05 : -0.3;
                    upper[j] = (j == 1) ? 0.05 : 0.3;
                }
            }
            double[][] q = new double[n][n];
            double[] linear = new double[n];
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    double s = 0.0;
                    for (int i = 0; i < m; i++) {
                        s += a[r * m + i] * a[c * m + i];
                    }
                    q[r][c] = s;
                }
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    s += a[r * m + i] * b[i];
                }
                linear[r] = -s;
            }
            double[] expected = exhaustive(q, linear, lower, upper);
            if (expected == null) {
                continue;
            }
            BoundedLevenbergMarquardt tight = new BoundedLevenbergMarquardt(1.0e-12, 1.0e-12, 0.0, 500, 1.0e-3);
            BoundedLevenbergMarquardt.Result r = tight.solve(new Linear(a, b, m, n), new double[n], m, lower, upper);
            assertArrayEquals("trial " + trial, expected, r.parameters, 1.0e-8);
            checked++;
        }
        assertTrue("the enumeration answered too few problems : " + checked, checked >= 15);
    }

    @Test
    public void testTheCovarianceOfTheFreeBlockMatchesOrdinaryLeastSquares() {
        // on a linear model with a bound that binds, the free parameters are an
        // ordinary least squares fit over the columns that are left, which OLS
        // answers through the singular values
        int checked = 0;
        for (int trial = 0; trial < 10; trial++) {
            int m = 40;
            int n = 5;
            Lcg lcg = new Lcg(555L + trial);
            double[] a = lcg.array(m * n);
            double[] b = lcg.array(m);
            double[] lower = fill(n, Double.NEGATIVE_INFINITY);
            double[] upper = fill(n, Double.POSITIVE_INFINITY);
            lower[0] = 0.5;

            BoundedLevenbergMarquardt tight = new BoundedLevenbergMarquardt(1.0e-13, 1.0e-13, 0.0, 500, 1.0e-3);
            BoundedLevenbergMarquardt.Result r = tight.solve(new Linear(a, b, m, n), new double[] { 1.0, 0, 0, 0, 0 },
                    m, lower, upper);
            if (r.activeSet[0] != BoundedLeastSquares.Bound.AT_LOWER) {
                continue;
            }
            assertNotNull("the covariance of the free block exists here", r.covariance);
            assertEquals("a held parameter has no variance", 0.0, r.covariance[0], 0.0);
            assertEquals("a held parameter has no standard error", 0.0, r.standardErrors[0], 0.0);
            assertEquals("a held parameter consumed no degree of freedom", m - (n - 1), r.degreesOfFreedom);

            // the reduced problem: the held column moved to the right hand side
            int free = n - 1;
            double[] reduced = new double[m * free];
            System.arraycopy(a, m, reduced, 0, m * free);
            double[] target = new double[m];
            for (int i = 0; i < m; i++) {
                target[i] = b[i] - a[i] * r.parameters[0];
            }
            math.linalg.LSSummary ols = OLS.estimate(0.05, matrix(m, free, reduced), matrix(m, 1, target));
            math.list.DoubleList expected = ols.getCoefficientStandardErrors();
            for (int j = 0; j < free; j++) {
                assertEquals("trial " + trial + ", standard error " + j, expected.get(j), r.standardErrors[j + 1],
                        1.0e-9 * Math.max(1.0, Math.abs(expected.get(j))));
            }
            checked++;
        }
        assertTrue("the bound never bound, so the test proved nothing : " + checked, checked >= 5);
    }

    // --------------------------------------------------- without a Jacobian --

    @Test
    public void testTheDerivativeFreePathReachesThePublishedMinima() {
        for (Problem p : MghProblems.all()) {
            BoundedLevenbergMarquardt.Result r = solver(20 * p.maxEvaluations).solve((DVectorFunction) p, p.start,
                    p.m, fill(p.n, Double.NEGATIVE_INFINITY), fill(p.n, Double.POSITIVE_INFINITY));
            // measured worst over the collection: 6.3e-7 relative, on Meyer
            assertEquals(p.name + " (" + r.status + ")", p.minimum, r.sumOfSquares,
                    1.0e-4 * Math.max(1.0, Math.abs(p.minimum)));
            assertTrue(p.name + ": the approximation should cost more than the derivative",
                    r.functionEvaluations > r.jacobianEvaluations);
        }
    }

    @Test
    public void testAModelThatIsUndefinedOutsideItsBoxIsFittedAnyway() {
        // this is the whole reason the box is carried inside the step. A
        // forward difference that probed outward, or an inner solve that walked
        // out of the box, would see NaN here and MINPACK's answer to a Jacobian
        // it cannot form is silent success
        for (Problem p : MghProblems.all()) {
            if (p.n > 4) {
                continue;
            }
            double[] lower = new double[p.n];
            double[] upper = new double[p.n];
            for (int j = 0; j < p.n; j++) {
                lower[j] = Math.min(-1.0, p.start[j]);
                upper[j] = Math.max(1.0, p.start[j]);
            }
            Undefined model = new Undefined(p, lower, upper);
            BoundedLevenbergMarquardt.Result r = solver(20 * p.maxEvaluations).solve(model, p.start, p.m, lower,
                    upper);
            assertEquals(p.name + ": the model was asked for a value outside its domain", 0, model.outside);
            assertTrue(p.name + ": the residuals must be finite", isFinite(r.residuals));
            assertFalse(p.name + ": the sum of squares must be a number", Double.isNaN(r.sumOfSquares));

            // it also has to be a genuine minimum of the restricted problem,
            // not merely a finite number: the analytic path over the same box
            // is the reference
            BoundedLevenbergMarquardt.Result reference = solver(20 * p.maxEvaluations).solve(p, p.start, p.m, lower,
                    upper);
            assertTrue(p.name + ": differenced " + r.sumOfSquares + " against analytic "
                    + reference.sumOfSquares,
                    r.sumOfSquares <= reference.sumOfSquares + 1.0e-6 * Math.max(1.0, reference.sumOfSquares));
        }
    }

    @Test
    public void testTheDifferenceIsTakenInwardsAtABound() {
        // a parameter sitting on its upper bound may only be probed downwards
        double bound = 2.0;
        double[] t = grid(30, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        double[] lower = { 0.0, 0.01, -bound };
        double[] upper = { 100.0, 100.0, bound };
        Recorder recorder = new Recorder(new Decay(t, y));
        solver(2000).solve(recorder, new double[] { 1.0, 0.5, bound }, t.length, lower, upper);
        assertTrue("the recorder saw nothing", recorder.calls > 5);
        for (int k = 0; k < recorder.calls; k++) {
            assertTrue("the third parameter was probed above its bound : " + recorder.thirdParameter[k],
                    recorder.thirdParameter[k] <= bound);
        }
    }

    @Test
    public void testABoxNarrowerThanTheDifferenceStepStillProducesADerivative() {
        // the step has to shrink to fit rather than cross, and the fit must
        // still move: the parameter is pinned to a hair's width, the other two
        // are not
        double[] t = grid(30, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        double[] lower = { 0.0, 0.01, 0.5 - 1.0e-13 };
        double[] upper = { 100.0, 100.0, 0.5 + 1.0e-13 };
        Undefined model = new Undefined(new Decay(t, y), lower, upper);
        BoundedLevenbergMarquardt.Result r = solver(4000).solve(model, new double[] { 1.0, 0.5, 0.5 }, t.length,
                lower, upper);
        assertEquals("the difference step left the box", 0, model.outside);
        assertTrue("the fit did not move : " + r.status, r.iterations > 0);
        assertEquals("amplitude", 3.0, r.parameters[0], 1.0e-5);
        assertEquals("rate", 1.5, r.parameters[1], 1.0e-5);
    }

    @Test
    public void testAParameterNearZeroIsStillEstimated() {
        // the fit must not depend on where a parameter started. Without the
        // floor the step for an offset at 1e-10 is 1e-18, the residuals come
        // back unchanged in every bit, the column is zero and the parameter
        // is never moved again -- and the search reports success
        double[] t = grid(30, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        double[] lower = { 0.0, 0.01, -10.0 };
        double[] upper = { 100.0, 100.0, 10.0 };
        double[] starts = { 0.5, 1.0e-6, 1.0e-10, 1.0e-12, 0.0 };
        for (int k = 0; k < starts.length; k++) {
            Undefined model = new Undefined(new Decay(t, y), lower, upper);
            BoundedLevenbergMarquardt.Result r = solver(2000).solve(model,
                    new double[] { 1.0, 0.5, starts[k] }, t.length, lower, upper);
            String at = "from " + starts[k] + " (" + r.status + ")";
            assertEquals(at + ": the difference step left the box", 0, model.outside);
            assertTrue(at + ": the search did not converge", r.converged);
            assertEquals(at + ": amplitude", 3.0, r.parameters[0], 1.0e-6);
            assertEquals(at + ": rate", 1.5, r.parameters[1], 1.0e-6);
            assertEquals(at + ": offset", 0.5, r.parameters[2], 1.0e-6);
        }
    }

    @Test
    public void testTheStepDoesNotShrinkWithTheParameterItIsTakenFrom() {
        // the first Jacobian is one evaluation at the point plus one probe
        // per parameter, so the first four calls hold it. The probe of the
        // third parameter has to sit a step away from 1e-10 rather than
        // 1e-10 times a step
        double start = 1.0e-10;
        double[] t = grid(30, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        double[] lower = { 0.0, 0.01, -10.0 };
        double[] upper = { 100.0, 100.0, 10.0 };
        Recorder recorder = new Recorder(new Decay(t, y));
        solver(2000).solve(recorder, new double[] { 1.0, 0.5, start }, t.length, lower, upper);
        assertTrue("the first Jacobian was never formed : " + recorder.calls, recorder.calls >= 4);
        double largest = 0.0;
        for (int k = 0; k < 4; k++) {
            largest = Math.max(largest, Math.abs(recorder.thirdParameter[k] - start));
        }
        assertTrue("the third parameter was probed with a step of " + largest, largest > 1.0e-9);
    }

    // ------------------------------------------------------ degenerate ends --

    @Test
    public void testEveryBoundInfiniteReproducesTheUnconstrainedFit() {
        double[] t = grid(50, 0.0, 5.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.01);
        Decay model = new Decay(t, y);
        double[] start = { 1.0, 0.5, 0.0 };
        BoundedLevenbergMarquardt.Result r = solver(2000).solve(model, start, t.length,
                fill(3, Double.NEGATIVE_INFINITY), fill(3, Double.POSITIVE_INFINITY));
        LevenbergMarquardt reference = new LevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, 2000, 100.0);
        LevenbergMarquardt.Result expected = reference.solve(model, start, t.length);
        assertEquals("no parameter can be at a bound when there is none", 0, r.atBound);
        for (int j = 0; j < 3; j++) {
            assertEquals("parameter " + j, expected.parameters[j], r.parameters[j], 1.0e-8);
        }
    }

    @Test
    public void testASingleParameter() {
        Lcg lcg = new Lcg(4242L);
        int m = 20;
        double[] a = lcg.array(m);
        double[] b = lcg.array(m);
        double[] lower = { 0.25 };
        double[] upper = { Double.POSITIVE_INFINITY };
        BoundedLevenbergMarquardt.Result r = new BoundedLevenbergMarquardt(1.0e-12, 1.0e-12, 0.0, 500, 1.0e-3)
                .solve(new Linear(a, b, m, 1), new double[] { 1.0 }, m, lower, upper);
        BoundedLeastSquares.Result exact = BoundedLeastSquares.bounded(matrix(m, 1, a), matrix(m, 1, b), lower,
                upper);
        assertEquals(exact.activeSet[0], r.activeSet[0]);
        assertEquals(exact.solution[0], r.parameters[0], 1.0e-9);
    }

    @Test
    public void testEveryParameterCanEndUpAtABound() {
        // a box placed entirely on one side of the unconstrained solution
        double[] t = grid(30, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        double[] lower = { 10.0, 20.0, 5.0 };
        double[] upper = { 11.0, 21.0, 6.0 };
        BoundedLevenbergMarquardt.Result r = solver(2000).solve(new Decay(t, y), new double[] { 10.5, 20.5, 5.5 },
                t.length, lower, upper);
        assertEquals("every parameter should have been pushed onto a bound", 3, r.atBound);
        for (int j = 0; j < 3; j++) {
            assertTrue("parameter " + j + " is not exactly on a bound",
                    r.parameters[j] == lower[j] || r.parameters[j] == upper[j]);
        }
        assertEquals("nothing is free, so nothing is estimated", t.length, r.degreesOfFreedom);
        assertNull("no free parameter means no covariance", r.covariance);
        assertNull("and no standard errors", r.standardErrors);
    }

    @Test
    public void testAJacobianWithAZeroColumn() {
        // a parameter the residuals do not depend on: the scaling must not
        // divide by its column norm, and the search must still finish
        double[] t = grid(25, 0.0, 4.0);
        double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
        Ignored model = new Ignored(new Decay(t, y));
        double[] lower = { 0.0, 0.01, -10.0, -1.0 };
        double[] upper = { 100.0, 100.0, 10.0, 1.0 };
        BoundedLevenbergMarquardt.Result r = solver(2000).solve(model, new double[] { 1.0, 0.5, 0.0, 0.3 },
                t.length, lower, upper);
        assertTrue("the search should still finish : " + r.status, r.converged);
        assertEquals("amplitude", 3.0, r.parameters[0], 1.0e-6);
        assertEquals("the ignored parameter cannot move", 0.3, r.parameters[3], 1.0e-12);
        assertNull("a rank deficient free block has no covariance", r.covariance);
    }

    @Test
    public void testAnExhaustedBudgetIsNotReportedAsConvergence() {
        Problem hard = null;
        for (Problem p : MghProblems.all()) {
            if ("Meyer".equals(p.name)) {
                hard = p;
            }
        }
        assertNotNull(hard);
        BoundedLevenbergMarquardt starved = new BoundedLevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, 4, 1.0e-3);
        BoundedLevenbergMarquardt.Result r = starved.solve(hard, hard.start, hard.m,
                fill(hard.n, Double.NEGATIVE_INFINITY), fill(hard.n, Double.POSITIVE_INFINITY));
        assertEquals(BoundedLevenbergMarquardt.Status.TOO_MANY_EVALUATIONS, r.status);
        assertFalse("an exhausted budget is not convergence", r.converged);
        assertFalse("but the point it reached still comes back", r.status.isSuccess());
        assertTrue("and it is a real point", isFinite(r.parameters));
    }

    @Test
    public void testAModelThatIsNotFiniteAtTheStartIsReportedRatherThanThrown() {
        // Helical valley on x >= 0 clamps its start onto the origin, where the
        // model is undefined. That is a fact about the caller's problem, and it
        // has to be said rather than thrown or silently called converged
        Problem helical = null;
        for (Problem p : MghProblems.all()) {
            if ("Helical valley".equals(p.name)) {
                helical = p;
            }
        }
        assertNotNull(helical);
        BoundedLevenbergMarquardt.Result r = solver(500).solve(helical, helical.start, helical.m,
                fill(helical.n, 0.0), fill(helical.n, Double.POSITIVE_INFINITY));
        assertEquals(BoundedLevenbergMarquardt.Status.NOT_FINITE, r.status);
        assertFalse(r.converged);
        assertEquals("it should have stopped at the first evaluation", 1, r.functionEvaluations);
        assertArrayEquals("the clamped starting point comes back", new double[] { 0.0, 0.0, 0.0 }, r.parameters,
                0.0);
    }

    @Test
    public void testTheCovarianceIsNullWhereItCannotExist() {
        // as many residuals as parameters leaves no spread to estimate from
        Lcg lcg = new Lcg(31415L);
        int n = 4;
        double[] a = lcg.array(n * n);
        double[] b = lcg.array(n);
        BoundedLevenbergMarquardt.Result r = solver(500).solve(new Linear(a, b, n, n), new double[n], n,
                fill(n, -10.0), fill(n, 10.0));
        assertEquals("m == n leaves no degrees of freedom", 0, r.degreesOfFreedom);
        assertNull(r.covariance);
        assertNull(r.standardErrors);
    }

    @Test
    public void testDetectionIsInvariantUnderARescalingOfTheParameters() {
        // the damping multiplies the squared column norms of the Jacobian, so
        // stretching a parameter must not change which bound binds
        double reference = Double.NaN;
        for (double factor : new double[] { 1.0, 1.0e3, 1.0e6 }) {
            double[] t = grid(40, 0.0, 5.0);
            double[] y = decayData(t, 3.0, 1.5, 0.5, 0.0);
            Stretched model = new Stretched(new Decay(t, y), factor);
            // the true amplitude is 3, so the upper bound of 2 binds and the
            // rate that comes out is the constrained one rather than 1.5
            double[] lower = { 0.0, 0.01 / factor, -10.0 };
            double[] upper = { 2.0, 100.0 / factor, 10.0 };
            BoundedLevenbergMarquardt.Result r = solver(4000).solve(model,
                    new double[] { 1.0, 0.5 / factor, 0.0 }, t.length, lower, upper);
            assertEquals("factor " + factor + ": the amplitude must be pushed onto its upper bound",
                    BoundedLeastSquares.Bound.AT_UPPER, r.activeSet[0]);
            assertEquals("factor " + factor + ": the amplitude", 2.0, r.parameters[0], 0.0);
            assertEquals("factor " + factor + ": the rate is free", BoundedLeastSquares.Bound.FREE,
                    r.activeSet[1]);
            double rate = r.parameters[1] * factor;
            if (Double.isNaN(reference)) {
                reference = rate;
            } else {
                assertEquals("factor " + factor + ": stretching the rate changed the fit", reference, rate,
                        1.0e-6 * reference);
            }
        }
    }

    // ----------------------------------------------------------- rejections --

    @Test
    public void testTheConstructorRejectsNonsense() {
        expectIae("negative valueTolerance", -1.0, SQRT_EPS, 0.0, 100, 1.0);
        expectIae("negative stepTolerance", SQRT_EPS, -1.0, 0.0, 100, 1.0);
        expectIae("negative kktTolerance", SQRT_EPS, SQRT_EPS, -1.0, 100, 1.0);
        expectIae("negative budget", SQRT_EPS, SQRT_EPS, 0.0, -3, 1.0);
        expectIae("zero damping", SQRT_EPS, SQRT_EPS, 0.0, 100, 0.0);
        expectIae("infinite damping", SQRT_EPS, SQRT_EPS, 0.0, 100, Double.POSITIVE_INFINITY);
        try {
            new BoundedLevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, 100, 1.0, -1.0);
            fail("a negative functionAccuracy should be rejected");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    @Test
    public void testTheArgumentsAreChecked() {
        DiffDVectorFunction f = new Linear(new double[] { 1.0, 1.0 }, new double[] { 1.0, 1.0 }, 2, 1);
        double[] start = { 0.0 };
        double[] lower = { 0.0 };
        double[] upper = { 1.0 };
        expectIae("null function", null, start, 2, lower, upper);
        expectIae("null start", f, null, 2, lower, upper);
        expectIae("null lower", f, start, 2, null, upper);
        expectIae("null upper", f, start, 2, lower, null);
        expectIae("empty start", f, new double[0], 2, new double[0], new double[0]);
        expectIae("too few residuals", f, start, 0, lower, upper);
        expectIae("bounds of the wrong length", f, start, 2, new double[] { 0.0, 0.0 }, upper);
        expectIae("a NaN bound", f, start, 2, new double[] { Double.NaN }, upper);
        expectIae("lower above upper", f, start, 2, new double[] { 2.0 }, upper);
        expectIae("an empty box", f, start, 2, new double[] { Double.POSITIVE_INFINITY }, upper);
        expectIae("a non-finite start", f, new double[] { Double.NaN }, 2, lower, upper);
    }

    // ------------------------------------------------------------- fixtures --

    private static BoundedLevenbergMarquardt solver(int budget) {
        return new BoundedLevenbergMarquardt(SQRT_EPS, SQRT_EPS, 0.0, budget, 1.0e-3);
    }

    /** A linear model stated as a nonlinear one: the Jacobian is constant. */
    private static final class Linear implements DiffDVectorFunction {

        private final double[] a;
        private final double[] b;
        private final int m;
        private final int n;

        Linear(double[] a, double[] b, int m, int n) {
            this.a = a;
            this.b = b;
            this.m = m;
            this.n = n;
        }

        @Override
        public void valueAt(double[] x, double[] r) {
            for (int i = 0; i < m; i++) {
                r[i] = -b[i];
            }
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    r[i] += a[j * m + i] * x[j];
                }
            }
        }

        @Override
        public void jacobianAt(double[] x, double[] jacobian) {
            System.arraycopy(a, 0, jacobian, 0, m * n);
        }
    }

    /** {@code A exp(-k t) + C}, the shape a curve fit actually has. */
    private static final class Decay implements DiffDVectorFunction {

        private final double[] t;
        private final double[] y;

        Decay(double[] t, double[] y) {
            this.t = t;
            this.y = y;
        }

        @Override
        public void valueAt(double[] p, double[] r) {
            for (int i = 0; i < t.length; i++) {
                r[i] = p[0] * Math.exp(-p[1] * t[i]) + p[2] - y[i];
            }
        }

        @Override
        public void jacobianAt(double[] p, double[] jacobian) {
            int m = t.length;
            for (int i = 0; i < m; i++) {
                double e = Math.exp(-p[1] * t[i]);
                jacobian[i] = e;
                jacobian[m + i] = -p[0] * t[i] * e;
                jacobian[2 * m + i] = 1.0;
            }
        }
    }

    /** The same model in coordinates where the rate is stretched by a factor. */
    private static final class Stretched implements DiffDVectorFunction {

        private final DiffDVectorFunction inner;
        private final double factor;
        private final double[] scratch = new double[3];

        Stretched(DiffDVectorFunction inner, double factor) {
            this.inner = inner;
            this.factor = factor;
        }

        @Override
        public void valueAt(double[] p, double[] r) {
            scratch[0] = p[0];
            scratch[1] = p[1] * factor;
            scratch[2] = p[2];
            inner.valueAt(scratch, r);
        }

        @Override
        public void jacobianAt(double[] p, double[] jacobian) {
            scratch[0] = p[0];
            scratch[1] = p[1] * factor;
            scratch[2] = p[2];
            inner.jacobianAt(scratch, jacobian);
            int m = jacobian.length / 3;
            for (int i = 0; i < m; i++) {
                jacobian[m + i] *= factor;
            }
        }
    }

    /** The same model with a fourth parameter nothing depends on. */
    private static final class Ignored implements DiffDVectorFunction {

        private final DiffDVectorFunction inner;

        Ignored(DiffDVectorFunction inner) {
            this.inner = inner;
        }

        @Override
        public void valueAt(double[] p, double[] r) {
            inner.valueAt(new double[] { p[0], p[1], p[2] }, r);
        }

        @Override
        public void jacobianAt(double[] p, double[] jacobian) {
            int m = jacobian.length / 4;
            double[] inner3 = new double[3 * m];
            inner.jacobianAt(new double[] { p[0], p[1], p[2] }, inner3);
            System.arraycopy(inner3, 0, jacobian, 0, 3 * m);
            Arrays.fill(jacobian, 3 * m, 4 * m, 0.0);
        }
    }

    /** Records where the model was asked for a value, and refuses nothing. */
    private static final class Caged implements DiffDVectorFunction {

        private final DiffDVectorFunction inner;
        private final double[] lower;
        private final double[] upper;
        int violations;
        int calls;

        Caged(DiffDVectorFunction inner, double[] lower, double[] upper) {
            this.inner = inner;
            this.lower = lower;
            this.upper = upper;
        }

        private void check(double[] x) {
            calls++;
            for (int j = 0; j < x.length; j++) {
                if (x[j] < lower[j] || x[j] > upper[j]) {
                    violations++;
                }
            }
        }

        @Override
        public void valueAt(double[] x, double[] r) {
            check(x);
            inner.valueAt(x, r);
        }

        @Override
        public void jacobianAt(double[] x, double[] jacobian) {
            check(x);
            inner.jacobianAt(x, jacobian);
        }
    }

    /** A model with no derivative that is genuinely undefined outside its box. */
    private static final class Undefined implements DVectorFunction {

        private final DiffDVectorFunction inner;
        private final double[] lower;
        private final double[] upper;
        int outside;

        Undefined(DiffDVectorFunction inner, double[] lower, double[] upper) {
            this.inner = inner;
            this.lower = lower;
            this.upper = upper;
        }

        @Override
        public void valueAt(double[] x, double[] r) {
            for (int j = 0; j < x.length; j++) {
                if (x[j] < lower[j] || x[j] > upper[j]) {
                    outside++;
                    Arrays.fill(r, Double.NaN);
                    return;
                }
            }
            inner.valueAt(x, r);
        }
    }

    /** Remembers the third parameter of every point it was called at. */
    private static final class Recorder implements DVectorFunction {

        private final DiffDVectorFunction inner;
        final double[] thirdParameter = new double[100000];
        int calls;

        Recorder(DiffDVectorFunction inner) {
            this.inner = inner;
        }

        @Override
        public void valueAt(double[] x, double[] r) {
            if (calls < thirdParameter.length) {
                thirdParameter[calls] = x[2];
            }
            calls++;
            inner.valueAt(x, r);
        }
    }

    /** The sum of squares of a residual function, as a scalar for SPG. */
    private static final class Squares implements DiffDMultiFunction {

        private final DiffDVectorFunction f;
        private final int m;
        private final int n;

        Squares(DiffDVectorFunction f, int m, int n) {
            this.f = f;
            this.m = m;
            this.n = n;
        }

        @Override
        public double apply(double[] x) {
            double[] r = new double[m];
            f.valueAt(x, r);
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                sum += r[i] * r[i];
            }
            return sum;
        }

        @Override
        public void derivativeAt(double[] x, double[] gradient) {
            double[] r = new double[m];
            double[] jacobian = new double[m * n];
            f.valueAt(x, r);
            f.jacobianAt(x, jacobian);
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int i = 0; i < m; i++) {
                    sum += jacobian[j * m + i] * r[i];
                }
                gradient[j] = 2.0 * sum;
            }
        }
    }

    /** The deterministic generator this repository uses in tests. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            this.state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
        }

        double[] array(int length) {
            double[] a = new double[length];
            for (int i = 0; i < length; i++) {
                a[i] = next();
            }
            return a;
        }
    }

    private static double[] lowerFor(Problem p, int variant) {
        double[] lower = fill(p.n, Double.NEGATIVE_INFINITY);
        for (int j = 0; j < p.n; j++) {
            if (variant == 0) {
                lower[j] = 0.0;
            } else if (variant == 1) {
                lower[j] = -1.0;
            } else if (variant == 3 && j % 2 == 0) {
                lower[j] = 0.1;
            }
        }
        return lower;
    }

    private static double[] upperFor(Problem p, int variant) {
        double[] upper = fill(p.n, Double.POSITIVE_INFINITY);
        for (int j = 0; j < p.n; j++) {
            if (variant == 1) {
                upper[j] = 1.0;
            } else if (variant == 2) {
                upper[j] = 0.5;
            } else if (variant == 3 && j % 2 == 1) {
                upper[j] = 2.0;
            }
        }
        return upper;
    }

    private static double[] fill(int n, double value) {
        double[] a = new double[n];
        Arrays.fill(a, value);
        return a;
    }

    private static double[] grid(int m, double from, double to) {
        double[] t = new double[m];
        for (int i = 0; i < m; i++) {
            t[i] = from + (to - from) * i / (m - 1.0);
        }
        return t;
    }

    private static double[] decayData(double[] t, double a, double k, double c, double noise) {
        Lcg lcg = new Lcg(20260821L);
        double[] y = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            y[i] = a * Math.exp(-k * t[i]) + c + noise * lcg.next();
        }
        return y;
    }

    private static DMatrix matrix(int rows, int cols, double[] values) {
        DMatrix d = new DMatrix(rows, cols);
        System.arraycopy(values, 0, d.getArrayUnsafe(), 0, rows * cols);
        return d;
    }

    private static boolean isFinite(double[] values) {
        for (int i = 0; i < values.length; i++) {
            if (Double.isNaN(values[i]) || Double.isInfinite(values[i])) {
                return false;
            }
        }
        return true;
    }

    private static double norm(double[] values) {
        double sum = 0.0;
        for (int i = 0; i < values.length; i++) {
            sum += values[i] * values[i];
        }
        return Math.sqrt(sum);
    }

    /**
     * The exact minimizer of {@code 0.5 x'Qx + b'x} over a box, by trying every
     * assignment of every coordinate to lower bound, free or upper bound. The
     * one that solves its reduced system inside the box and has the right
     * multiplier signs is the solution.
     */
    private static double[] exhaustive(double[][] q, double[] b, double[] lower, double[] upper) {
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
            boolean usable = true;
            for (int i = 0; i < n; i++) {
                if (state[i] != 1 && Double.isInfinite(x[i])) {
                    usable = false;
                }
            }
            if (!usable) {
                continue;
            }
            if (freeCount > 0) {
                int[] index = new int[freeCount];
                int k = 0;
                for (int i = 0; i < n; i++) {
                    if (state[i] == 1) {
                        index[k] = i;
                        k++;
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
            double[] swap = a[k];
            a[k] = a[pivot];
            a[pivot] = swap;
            double t = b[k];
            b[k] = b[pivot];
            b[pivot] = t;
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

    private static void expectIae(String what, double valueTolerance, double stepTolerance, double kktTolerance,
            int maxEvaluations, double initialDamping) {
        try {
            new BoundedLevenbergMarquardt(valueTolerance, stepTolerance, kktTolerance, maxEvaluations,
                    initialDamping);
            fail(what + " should have been rejected");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    private static void expectIae(String what, DiffDVectorFunction f, double[] start, int residualCount,
            double[] lower, double[] upper) {
        try {
            new BoundedLevenbergMarquardt().solve(f, start, residualCount, lower, upper);
            fail(what + " should have been rejected");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }
}
