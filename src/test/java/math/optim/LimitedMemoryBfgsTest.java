package math.optim;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.optim.OwlQnTest.Rosenbrock;

/**
 * Tests for {@link LimitedMemoryBFGS} and the {@link BackTrackLineSearch} it
 * drives, concentrated on the edges an objective outside the smooth strictly
 * concave case reaches: a violated curvature condition, a line search that
 * cannot step, and an aborted run.
 */
public class LimitedMemoryBfgsTest {

    /** {@code -(a - x)^2 - b (y - x^2)^2}, maximized at {@code (a, a^2)}. */
    private static final class Concave implements Optimizable.ByGradientValue {

        private final double[] w = new double[2];
        private final boolean poisonGradient;

        Concave() {
            this(false);
        }

        Concave(boolean poisonGradient) {
            this.poisonGradient = poisonGradient;
        }

        @Override
        public int getNumParameters() {
            return 2;
        }

        @Override
        public void getParameters(double[] buffer) {
            System.arraycopy(w, 0, buffer, 0, 2);
        }

        @Override
        public double getParameter(int index) {
            return w[index];
        }

        @Override
        public void setParameters(double[] params) {
            System.arraycopy(params, 0, w, 0, 2);
        }

        @Override
        public void setParameter(int index, double value) {
            w[index] = value;
        }

        @Override
        public double getValue() {
            return -3.0 * (w[0] - 2.0) * (w[0] - 2.0) - 5.0 * (w[1] + 1.0) * (w[1] + 1.0);
        }

        @Override
        public void getValueGradient(double[] buffer) {
            buffer[0] = poisonGradient ? Double.NaN : (-6.0 * (w[0] - 2.0));
            buffer[1] = -10.0 * (w[1] + 1.0);
        }
    }

    /**
     * The ordinary case still lands on the maximum. The accuracy is what the
     * default gradient tolerance of 1e-3 buys; the tightened one below is what
     * was unreachable before.
     */
    @Test
    public void testConcaveQuadraticIsMaximized() {
        Concave f = new Concave();
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f);

        assertTrue(bfgs.optimize());
        assertTrue(bfgs.isConverged());

        double[] w = new double[2];
        f.getParameters(w);
        assertEquals("w[0]", 2.0, w[0], 1.0e-4);
        assertEquals("w[1]", -1.0, w[1], 1.0e-4);
    }

    /**
     * Tightening the stopping rules of the search alone does not make this
     * class more accurate, and the constructor says so. What binds is the step
     * tolerance of the internal {@code BackTrackLineSearch}: the search stops
     * at the same point no matter how much precision is asked for -- it merely
     * stops by running out of steps instead of by meeting a tolerance. The
     * ceiling is real, which is what makes the test below worth having.
     */
    @Test
    public void testTighteningTheStoppingRulesAloneIsBoundedByTheLineSearch() {
        Concave loose = new Concave();
        assertTrue(new LimitedMemoryBFGS(loose, 1000, 1.0e-4, 1.0e-3, 4).optimize());
        double[] wLoose = new double[2];
        loose.getParameters(wLoose);

        Concave tight = new Concave();
        assertFalse("cannot meet a tolerance the line search cannot resolve",
                new LimitedMemoryBFGS(tight, 1000, 1.0e-30, 1.0e-30, 4).optimize());
        double[] wTight = new double[2];
        tight.getParameters(wTight);

        assertEquals("w[0] unchanged by the tighter request", wLoose[0], wTight[0], 0.0);
        assertEquals("w[1] unchanged by the tighter request", wLoose[1], wTight[1], 0.0);
        // and both are only as good as that step tolerance allows
        assertEquals("w[0]", 2.0, wTight[0], 1.0e-4);
        assertTrue("the default line search resolves no better than 1e-8",
                Math.abs(wTight[0] - 2.0) > 1.0e-8);
    }

    /**
     * Tightening the step tolerances along with them does. They were reachable
     * only through two setters on a package-private class whose only instance
     * was private to this one, so the ceiling above could not be lifted from
     * outside; they are constructor arguments now. On this quadratic the
     * default settles 2.8e-6 away from the maximum and the opened-up line
     * search lands on it exactly, which is five orders of magnitude of
     * accuracy that no caller could ask for before.
     */
    @Test
    public void testOpeningTheLineSearchLiftsTheCeiling() {
        Concave f = new Concave();
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f, 1000, 1.0e-30, 1.0e-30, 4, 1.0e-20, 1.0e-20);

        assertTrue(bfgs.optimize());
        assertTrue(bfgs.isConverged());

        double[] w = new double[2];
        f.getParameters(w);
        assertEquals("w[0]", 2.0, w[0], 0.0);
        assertEquals("w[1]", -1.0, w[1], 0.0);
    }

    /**
     * The step tolerances trade accuracy against work monotonically rather than
     * being an all-or-nothing switch, so a caller can ask for what it needs.
     */
    @Test
    public void testAccuracyImprovesMonotonicallyWithTheStepTolerances() {
        double[] steps = { 1.0e-4, 1.0e-8, 1.0e-12, 1.0e-16 };
        double previousError = Double.MAX_VALUE;
        for (int k = 0; k < steps.length; k++) {
            Concave f = new Concave();
            new LimitedMemoryBFGS(f, 1000, 1.0e-30, 1.0e-30, 4, steps[k], steps[k]).optimize();
            double[] w = new double[2];
            f.getParameters(w);

            double error = Math.abs(w[0] - 2.0) + Math.abs(w[1] + 1.0);
            assertTrue("step tolerance " + steps[k] + ": error " + error + " is not below " + previousError,
                    error < previousError);
            previousError = error;
        }
    }

    /**
     * The defaults of the new constructor are the values the line search always
     * had, so nothing moves for a caller that does not ask for them.
     */
    @Test
    public void testTheStepToleranceDefaultsReproduceTheOldBehaviour() {
        Concave shorthand = new Concave();
        new LimitedMemoryBFGS(shorthand, 1000, 1.0e-4, 1.0e-3, 4).optimize();
        double[] wShorthand = new double[2];
        shorthand.getParameters(wShorthand);

        Concave explicit = new Concave();
        new LimitedMemoryBFGS(explicit, 1000, 1.0e-4, 1.0e-3, 4, 1.0e-7, 1.0e-4).optimize();
        double[] wExplicit = new double[2];
        explicit.getParameters(wExplicit);

        assertEquals("w[0]", wShorthand[0], wExplicit[0], 0.0);
        assertEquals("w[1]", wShorthand[1], wExplicit[1], 0.0);
    }

    /**
     * The curvature condition {@code s . y < 0} this class needs holds for a
     * strictly concave objective and fails on one that is not. The Rosenbrock
     * function violates it from its classic starting point; skipping the update
     * keeps the search alive where treating it as a caller error ended it.
     */
    @Test
    public void testNonConcaveObjectiveSurvivesAViolatedCurvatureCondition() {
        Rosenbrock f = new Rosenbrock(new double[] { -1.2, 1.0 });
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f);

        bfgs.optimize();

        double[] w = new double[2];
        f.getParameters(w);
        assertEquals("w[0]", 1.0, w[0], 1.0e-3);
        assertEquals("w[1]", 1.0, w[1], 1.0e-3);
        assertTrue("value = " + (-f.getValue()), -f.getValue() < 1.0e-7);
    }

    /**
     * A line search that cannot step is the flat region near the maximum as
     * often as it is a failure. The parameters it settled on are the result and
     * must survive; only the convergence flag may say that the tolerances were
     * never met.
     */
    @Test
    public void testLineSearchThatCannotStepKeepsTheBestPointFound() {
        Rosenbrock f = new Rosenbrock(new double[] { 3.0, -1.0, 0.0, 1.0 });
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f);

        boolean reported = bfgs.optimize();

        double[] w = new double[4];
        f.getParameters(w);
        for (int i = 0; i < w.length; i++) {
            assertEquals("w[" + i + "]", 1.0, w[i], 1.0e-3);
        }
        assertTrue("value = " + (-f.getValue()), -f.getValue() < 1.0e-8);
        if (!reported) {
            assertFalse("a run that did not meet its tolerances is not converged", bfgs.isConverged());
        }
    }

    /**
     * Exhausting the iteration budget is not convergence. Before the budget was
     * reachable this could not be asserted at all.
     */
    @Test
    public void testBudgetExhaustionIsNotConvergence() {
        Rosenbrock f = new Rosenbrock(new double[] { -2.0, 3.0, -2.0, 3.0, -2.0, 3.0 });
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f, 3, 1.0e-14, 1.0e-14, 4);

        assertFalse("optimize() must not report success", bfgs.optimize());
        assertFalse("isConverged() must agree", bfgs.isConverged());
        assertTrue("iterations = " + bfgs.getIteration(), bfgs.getIteration() >= 3);
    }

    /** A tightened gradient tolerance reaches the maximum more closely. */
    @Test
    public void testTightenedGradientToleranceReachesTheMaximumMoreClosely() {
        Rosenbrock loose = new Rosenbrock(new double[] { -1.2, 1.0 });
        new LimitedMemoryBFGS(loose, 1000, 1.0e-4, 1.0e-1, 4).optimize();

        Rosenbrock tight = new Rosenbrock(new double[] { -1.2, 1.0 });
        new LimitedMemoryBFGS(tight, 1000, 1.0e-4, 1.0e-8, 4).optimize();

        assertTrue("loose = " + (-loose.getValue()) + ", tight = " + (-tight.getValue()),
                -tight.getValue() <= -loose.getValue());
    }

    /** An evaluator that stops the run has not made it converge. */
    @Test
    public void testEvaluatorAbortIsNotConvergence() {
        Rosenbrock f = new Rosenbrock(new double[] { -1.2, 1.0 });
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(f);
        bfgs.setEvaluator(new OptimizerEvaluator.ByGradient() {
            @Override
            public boolean evaluate(Optimizable.ByGradientValue maxable, int iter) {
                return false;
            }
        });

        assertFalse("optimize() must not report success", bfgs.optimize());
        assertFalse("isConverged() must agree", bfgs.isConverged());
    }

    /** The stopping rules are checked, not silently accepted. */
    @Test
    public void testArgumentValidation() {
        expectIae(null, 1000, 1.0e-4, 1.0e-3, 4);
        expectIae("maxIterations", 0, 1.0e-4, 1.0e-3, 4);
        expectIae("tolerance", 1000, 0.0, 1.0e-3, 4);
        expectIae("tolerance", 1000, Double.NaN, 1.0e-3, 4);
        expectIae("gradientTolerance", 1000, 1.0e-4, -1.0e-9, 4);
        expectIae("m", 1000, 1.0e-4, 1.0e-3, 0);
        expectIae("m", 1000, 1.0e-4, 1.0e-3, 101);
        expectIae("stepRelTolerance", 0.0, 1.0e-4);
        expectIae("stepRelTolerance", Double.NaN, 1.0e-4);
        expectIae("stepRelTolerance", -1.0e-7, 1.0e-4);
        expectIae("stepAbsTolerance", 1.0e-7, 0.0);
        expectIae("stepAbsTolerance", 1.0e-7, Double.POSITIVE_INFINITY);
    }

    private static void expectIae(String what, int maxIterations, double tolerance, double gradientTolerance,
            int m) {
        Optimizable.ByGradientValue f = (what == null) ? null : new Concave();
        try {
            new LimitedMemoryBFGS(f, maxIterations, tolerance, gradientTolerance, m);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    private static void expectIae(String what, double stepRelTolerance, double stepAbsTolerance) {
        try {
            new LimitedMemoryBFGS(new Concave(), 1000, 1.0e-4, 1.0e-3, 4, stepRelTolerance, stepAbsTolerance);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    /**
     * A NaN in the gradient used to be guarded by an assertion, which is
     * disabled in the shipped artifact, so it propagated into the parameters
     * instead.
     */
    @Test
    public void testLineSearchRejectsANaNGradient() {
        Concave f = new Concave(true);
        BackTrackLineSearch search = new BackTrackLineSearch(f);
        try {
            search.optimize(new double[] { 1.0, 1.0 }, 1.0);
            fail("expected InvalidOptimizableException");
        } catch (InvalidOptimizableException expected) {
            // as it should be
        }
    }

    /** Nothing has happened yet, and the class says so rather than saying false. */
    @Test
    public void testTerminationBeforeTheFirstRun() {
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(new Concave());

        assertEquals(Termination.NOT_STARTED, bfgs.getTermination());
        assertFalse(bfgs.isConverged());
        assertTrue("gradient norm " + bfgs.getGradientNorm(), Double.isNaN(bfgs.getGradientNorm()));
    }

    /**
     * Five different outcomes that a single boolean had to stand for, each one
     * naming itself now.
     */
    @Test
    public void testEachStoppingRuleNamesItself() {
        Rosenbrock tightened = new Rosenbrock(new double[] { -1.2, 1.0 });
        LimitedMemoryBFGS a = new LimitedMemoryBFGS(tightened, 1000, 1.0e-12, 1.0e-9, 4, 1.0e-14, 1.0e-14);
        assertTrue(a.optimize());
        assertEquals(Termination.VALUE_TOLERANCE, a.getTermination());
        assertTrue(a.getTermination().isConvergence());

        Rosenbrock budgeted = new Rosenbrock(new double[] { -2.0, 3.0, -2.0, 3.0, -2.0, 3.0 });
        LimitedMemoryBFGS b = new LimitedMemoryBFGS(budgeted, 3, 1.0e-14, 1.0e-14, 4);
        assertFalse(b.optimize());
        assertEquals(Termination.ITERATION_LIMIT, b.getTermination());

        Rosenbrock aborted = new Rosenbrock(new double[] { -1.2, 1.0 });
        LimitedMemoryBFGS c = new LimitedMemoryBFGS(aborted);
        c.setEvaluator(new OptimizerEvaluator.ByGradient() {
            @Override
            public boolean evaluate(Optimizable.ByGradientValue maxable, int iter) {
                return false;
            }
        });
        assertFalse(c.optimize());
        assertEquals(Termination.EVALUATOR_ABORTED, c.getTermination());

        Rosenbrock partial = new Rosenbrock(new double[] { -2.0, 3.0, -2.0, 3.0, -2.0, 3.0 });
        LimitedMemoryBFGS d = new LimitedMemoryBFGS(partial);
        assertFalse("two iterations of many are not a result", d.optimize(2));
        assertEquals(Termination.PARTIAL_RUN, d.getTermination());
        assertFalse(d.isConverged());

        Rosenbrock stalling = new Rosenbrock(new double[] { 3.0, -1.0, 0.0, 1.0 });
        LimitedMemoryBFGS e = new LimitedMemoryBFGS(stalling);
        assertFalse(e.optimize());
        assertEquals(Termination.LINE_SEARCH_STALLED, e.getTermination());

    }

    /** A gradient rule that fired is always one of the convergences. */
    @Test
    public void testStationaryImpliesConvergence() {
        for (Termination t : Termination.values()) {
            assertFalse(t + " is stationary but not a convergence", t.isStationary() && !t.isConvergence());
        }
        assertTrue(Termination.GRADIENT_TOLERANCE.isStationary());
        assertTrue(Termination.ZERO_GRADIENT.isStationary());
        assertFalse("a value rule says nothing about the gradient", Termination.VALUE_TOLERANCE.isStationary());
        assertTrue("but it is still one of the convergences", Termination.VALUE_TOLERANCE.isConvergence());
        assertFalse(Termination.NOT_STARTED.isConvergence());
    }

    /**
     * The reason the reason is worth having. Both of these runs report
     * {@code false}, and they mean opposite things: one stopped a hair from the
     * maximum, the other never got near it. Only the termination and the
     * gradient norm can tell them apart.
     */
    @Test
    public void testTheGradientNormSeparatesAStallFromAnExhaustedBudget() {
        Rosenbrock nearlyThere = new Rosenbrock(new double[] { -1.2, 1.0 });
        LimitedMemoryBFGS stalled = new LimitedMemoryBFGS(nearlyThere);
        boolean stallReported = stalled.optimize();

        Rosenbrock nowhereNear = new Rosenbrock(new double[] { -2.0, 3.0, -2.0, 3.0, -2.0, 3.0 });
        LimitedMemoryBFGS budget = new LimitedMemoryBFGS(nowhereNear, 3, 1.0e-14, 1.0e-14, 4);

        assertFalse("the budget run must not report success", budget.optimize());
        assertEquals(Termination.ITERATION_LIMIT, budget.getTermination());
        assertTrue("budget gradient norm " + budget.getGradientNorm(), budget.getGradientNorm() > 1.0);

        if (!stallReported) {
            assertEquals(Termination.LINE_SEARCH_STALLED, stalled.getTermination());
            assertTrue("stalled gradient norm " + stalled.getGradientNorm(), stalled.getGradientNorm() < 1.0e-2);
            assertTrue("the two must be orders of magnitude apart",
                    budget.getGradientNorm() > 1000.0 * stalled.getGradientNorm());
        }
    }

}
