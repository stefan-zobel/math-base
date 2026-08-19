package math.optim;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DFunction;
import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;
import math.fun.DiffDFunction;
import math.solve.RootFinder;

/**
 * Tests for {@link BrentMinimizer}. The strongest of these does not compare the
 * two variants against each other but against an entirely different question
 * answered by an entirely different algorithm: the abscissa of a minimum is a
 * root of the derivative, and {@link RootFinder#brentDekker} finds roots.
 */
public class BrentMinimizerTest {

    /** A function with a known minimum, its derivative, and a search range. */
    private static final class Case implements DiffDFunction {

        final String name;
        final double xMin;
        final double from;
        final double to;
        final double rootFrom;
        final double rootTo;
        private final DFunction f;
        private final DFunction df;
        int evaluations;

        Case(String name, double xMin, double from, double to, double rootFrom, double rootTo, DFunction f,
                DFunction df) {
            this.name = name;
            this.xMin = xMin;
            this.from = from;
            this.to = to;
            this.rootFrom = rootFrom;
            this.rootTo = rootTo;
            this.f = f;
            this.df = df;
        }

        @Override
        public double apply(double x) {
            evaluations++;
            return f.apply(x);
        }

        @Override
        public double derivativeAt(double x) {
            return df.apply(x);
        }
    }

    private static Case[] cases() {
        return new Case[] {
                new Case("(x-3)^2", 3.0, 0.0, 1.0, 0.0, 10.0, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return (x - 3.0) * (x - 3.0);
                    }
                }, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return 2.0 * (x - 3.0);
                    }
                }),
                // f' = x^2 (4x - 9), so the minimum is at 9/4
                new Case("x^4-3x^3+2", 2.25, 1.0, 1.5, 1.0, 5.0, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return x * x * x * x - 3.0 * x * x * x + 2.0;
                    }
                }, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return 4.0 * x * x * x - 9.0 * x * x;
                    }
                }),
                new Case("cos(x)", Math.PI, 1.0, 2.0, 1.0, 5.0, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return Math.cos(x);
                    }
                }, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return -Math.sin(x);
                    }
                }),
                // minimum exactly at zero, which is what ZEPS protects
                new Case("|x|^1.5", 0.0, -1.0, -0.5, -1.0, 1.0, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return Math.pow(Math.abs(x), 1.5);
                    }
                }, new DFunction() {
                    @Override
                    public double apply(double x) {
                        return 1.5 * Math.signum(x) * Math.sqrt(Math.abs(x));
                    }
                }) };
    }

    /** The bracket must actually straddle a minimum. */
    @Test
    public void testBracketStraddlesTheMinimum() {
        BrentMinimizer m = new BrentMinimizer();
        Case[] all = cases();
        for (int k = 0; k < all.length; k++) {
            Case c = all[k];
            BrentMinimizer.Bracket br = m.bracket(c, c.from, c.to);

            assertTrue(c.name, br.bracketed);
            assertTrue(c.name + ": f(b) <= f(a)", br.fb <= br.fa);
            assertTrue(c.name + ": f(b) <= f(c)", br.fb <= br.fc);
            double lo = Math.min(br.a, br.c);
            double hi = Math.max(br.a, br.c);
            assertTrue(c.name + ": b inside", lo <= br.b && br.b <= hi);
            assertTrue(c.name + ": the minimum is inside " + br, lo <= c.xMin && c.xMin <= hi);
            assertTrue(c.evaluations + " evaluations for " + c.name, c.evaluations >= 0);
        }
    }

    /** Bracketing works whichever way the initial interval points. */
    @Test
    public void testBracketFindsTheMinimumFromEitherSide() {
        BrentMinimizer m = new BrentMinimizer();
        BrentMinimizer.Bracket right = m.bracket(cases()[0], 0.0, 1.0);
        BrentMinimizer.Bracket left = m.bracket(cases()[0], 9.0, 8.0);
        BrentMinimizer.Bracket over = m.bracket(cases()[0], 0.0, 100.0);

        assertTrue("from the left", right.bracketed);
        assertTrue("from the right", left.bracketed);
        assertTrue("overshooting", over.bracketed);
        for (int k = 0; k < 3; k++) {
            BrentMinimizer.Bracket br = (k == 0) ? right : (k == 1) ? left : over;
            assertTrue(br.toString(), Math.min(br.a, br.c) <= 3.0 && 3.0 <= Math.max(br.a, br.c));
        }
    }

    /**
     * A function without a minimum in the search direction must be reported as
     * such, cheaply. Left to run on infinities, the loop ends only because a
     * comparison against a NaN is false, which is indistinguishable from
     * success.
     */
    @Test
    public void testFunctionWithoutAMinimumIsReportedAndCostsLittle() {
        final int[] evaluations = new int[1];
        DFunction downhillForever = new DFunction() {
            @Override
            public double apply(double x) {
                evaluations[0]++;
                return -(x * x);
            }
        };

        BrentMinimizer.Bracket br = new BrentMinimizer().bracket(downhillForever, 0.0, 0.01);

        assertFalse("no minimum exists, so none can be bracketed", br.bracketed);
        assertTrue("evaluations = " + evaluations[0], evaluations[0] < 2000);
    }

    /** Both variants find the known minima. */
    @Test
    public void testKnownMinima() {
        BrentMinimizer m = new BrentMinimizer();
        Case[] withDerivative = cases();
        Case[] withoutDerivative = cases();
        for (int k = 0; k < withDerivative.length; k++) {
            Case cd = withDerivative[k];
            BrentMinimizer.Result rd = m.minimize(cd, m.bracket(cd, cd.from, cd.to));
            assertTrue(cd.name + " (derivative)", rd.converged);
            assertEquals(cd.name + " (derivative)", cd.xMin, rd.x, 1.0e-7);

            Case cf = withoutDerivative[k];
            BrentMinimizer.Result rf = m.minimize((DFunction) cf, m.bracket(cf, cf.from, cf.to));
            assertTrue(cf.name + " (derivative free)", rf.converged);
            assertEquals(cf.name + " (derivative free)", cf.xMin, rf.x, 1.0e-7);

            assertEquals(cd.name + ": the two variants agree", rd.x, rf.x, 1.0e-7);
            assertEquals(cd.name + ": and on the value", rd.value, rf.value, 1.0e-12);
        }
    }

    /**
     * The abscissa of a minimum is a root of the derivative. Brent-Dekker finds
     * that root by an unrelated method, so the two must meet.
     */
    @Test
    public void testAgreesWithTheRootOfTheDerivative() {
        BrentMinimizer m = new BrentMinimizer();
        Case[] all = cases();
        for (int k = 0; k < all.length; k++) {
            final Case c = all[k];
            double root = RootFinder.brentDekker(c.rootFrom, c.rootTo, new DFunction() {
                @Override
                public double apply(double x) {
                    return c.derivativeAt(x);
                }
            }, 1.0e-12);

            BrentMinimizer.Result r = m.minimize(c, m.bracket(c, c.from, c.to));
            assertEquals(c.name + ": minimizer vs root of f'", root, r.x, 1.0e-7);
        }
    }

    /**
     * A minimum sitting exactly at zero is what the absolute floor under the
     * tolerance is for. Purely relative, the tolerance vanishes there and the
     * convergence test can never fire, so the search runs to its budget.
     */
    @Test
    public void testMinimumAtZeroConvergesQuickly() {
        Case c = cases()[3];
        BrentMinimizer m = new BrentMinimizer();
        BrentMinimizer.Result r = m.minimize(c, m.bracket(c, c.from, c.to));

        assertTrue("must not run to the budget", r.converged);
        assertEquals("x", 0.0, r.x, 1.0e-7);
        assertTrue("iterations = " + r.iterations, r.iterations < 60);
    }

    /** Running out of iterations is not convergence. */
    @Test
    public void testBudgetExhaustionIsNotConvergence() {
        Case c = cases()[2];
        BrentMinimizer m = new BrentMinimizer(1.0e-14, 2, 200);
        BrentMinimizer.Result r = m.minimize(c, m.bracket(c, c.from, c.to));

        assertFalse("optimize must not claim success", r.converged);
        assertEquals("iterations", 2, r.iterations);
    }

    /**
     * Cross-check against a different family of method entirely: the simplex
     * search, run in one dimension.
     */
    @Test
    public void testAgreesWithNelderMeadInOneDimension() {
        BrentMinimizer m = new BrentMinimizer();
        Case[] all = cases();
        for (int k = 0; k < all.length; k++) {
            final Case c = all[k];
            BrentMinimizer.Result r = m.minimize(c, m.bracket(c, c.from, c.to));

            DMultiFunctionEval nm = new NelderMead().minimize(new DMultiFunction() {
                @Override
                public double apply(double[] x) {
                    return c.apply(x[0]);
                }
            }, new double[] { c.from });

            assertEquals(c.name + ": Brent vs Nelder-Mead", nm.point[0], r.x, 1.0e-6);
        }
    }

    /** Arguments are checked, not silently accepted. */
    @Test
    public void testArgumentValidation() {
        expectIae("tol", 0.0, 100, 200);
        expectIae("tol", Double.NaN, 100, 200);
        expectIae("maxIterations", 1.0e-8, 0, 200);
        expectIae("maxBracketSteps", 1.0e-8, 100, 0);

        BrentMinimizer m = new BrentMinimizer();
        Case c = cases()[0];
        try {
            m.bracket(null, 0.0, 1.0);
            fail("expected IllegalArgumentException for a null function");
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
        try {
            m.bracket(c, Double.POSITIVE_INFINITY, 1.0);
            fail("expected IllegalArgumentException for a non-finite endpoint");
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
        try {
            m.bracket(c, 1.0, 1.0);
            fail("expected IllegalArgumentException for a degenerate interval");
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
        try {
            m.minimize(c, new BrentMinimizer().bracket(new DFunction() {
                @Override
                public double apply(double x) {
                    return -(x * x);
                }
            }, 0.0, 0.01));
            fail("expected IllegalArgumentException for an unbracketed input");
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    private static void expectIae(String what, double tol, int maxIterations, int maxBracketSteps) {
        try {
            new BrentMinimizer(tol, maxIterations, maxBracketSteps);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }
}
