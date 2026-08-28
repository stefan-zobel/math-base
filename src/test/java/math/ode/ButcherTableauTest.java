package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The coefficients of the shipped methods, checked against the order conditions
 * rather than against a table they were copied from.
 * <p>
 * A mistyped coefficient is the failure this class exists to catch, and it is a
 * silent one: nothing throws, the method simply advances at a lower order than
 * it claims, and the symptom turns up much later as an accuracy that does not
 * improve the way halving the step size says it should. The order conditions
 * decide the question outright and without solving anything -- for order
 * {@code p} there is one equation per rooted tree with at most {@code p} nodes,
 * 8 of them through order four and 17 through order five, and a tableau either
 * satisfies them or is not of that order.
 * <p>
 * The measurements the tolerances below rest on: over both tableaux the row
 * sums deviate from the stage times by at most {@code 2.2e-16}, every condition
 * a method is supposed to satisfy comes out below {@code 2.3e-16}, and every
 * condition it is supposed to <em>fail</em> is missed by at least
 * {@code 6.6e-05}. Fourteen orders of magnitude separate the two, so nothing
 * here is close to its threshold.
 */
public final class ButcherTableauTest {

    /** Measured worst case 2.3e-16 over both tableaux and all conditions. */
    private static final double SATISFIED = 1.0e-14;
    /** Measured smallest violation 6.6e-05, on the embedded Dormand-Prince pair. */
    private static final double VIOLATED = 1.0e-06;

    @Test
    public void testTheRowsOfTheStageMatrixSumToTheStageTimes() {
        assertEquals(0.0, worstRowSum(ButcherTableau.CLASSIC_RK4), 0.0);
        assertTrue(worstRowSum(ButcherTableau.DORMAND_PRINCE_45) <= SATISFIED);
    }

    @Test
    public void testTheClassicalMethodSatisfiesTheConditionsThroughOrderFour() {
        double[] r = conditions(ButcherTableau.CLASSIC_RK4, ButcherTableau.CLASSIC_RK4.b());
        assertSatisfied(r, 0, 8);
    }

    @Test
    public void testTheClassicalMethodFailsEveryConditionOfOrderFive() {
        double[] r = conditions(ButcherTableau.CLASSIC_RK4, ButcherTableau.CLASSIC_RK4.b());
        assertViolated(r, 8, 17);
    }

    @Test
    public void testDormandPrinceSatisfiesAllSeventeenConditionsThroughOrderFive() {
        double[] r = conditions(ButcherTableau.DORMAND_PRINCE_45, ButcherTableau.DORMAND_PRINCE_45.b());
        assertSatisfied(r, 0, 17);
    }

    @Test
    public void testTheEmbeddedSolutionIsOfOrderFour() {
        double[] r = conditions(ButcherTableau.DORMAND_PRINCE_45, ButcherTableau.DORMAND_PRINCE_45.bStar());
        assertSatisfied(r, 0, 8);
    }

    /**
     * The half of the embedded pair that makes it an error estimate rather than
     * a second answer to the same question: if the lower order weights happened
     * to satisfy the order five conditions too, the difference of the two
     * solutions would be zero to leading order and would estimate nothing.
     */
    @Test
    public void testTheEmbeddedSolutionFailsEveryConditionOfOrderFive() {
        double[] r = conditions(ButcherTableau.DORMAND_PRINCE_45, ButcherTableau.DORMAND_PRINCE_45.bStar());
        assertViolated(r, 8, 17);
    }

    @Test
    public void testTheLastStageOfDormandPrinceIsTheFirstStageOfTheNextStep() {
        ButcherTableau t = ButcherTableau.DORMAND_PRINCE_45;
        double[] b = t.b();
        double[][] a = t.a();
        int last = b.length - 1;
        assertEquals(1.0, t.c()[last], 0.0);
        assertEquals(0.0, b[last], 0.0);
        assertEquals(last, a[last].length);
        for (int j = 0; j < last; ++j) {
            // bit for bit: the two arrays are written from the same literals
            assertEquals("a[" + last + "][" + j + "]", b[j], a[last][j], 0.0);
        }
        assertTrue(t.isFsal());
        assertFalse(ButcherTableau.CLASSIC_RK4.isFsal());
    }

    @Test
    public void testTheContinuousExtensionReproducesTheEndsOfTheStep() {
        ButcherTableau t = ButcherTableau.DORMAND_PRINCE_45;
        double[][] d = t.dense();
        double[] b = t.b();
        for (int i = 0; i < d.length; ++i) {
            double atEnd = 0.0;
            for (int j = 0; j < d[i].length; ++j) {
                atEnd += d[i][j];
            }
            assertEquals("stage " + i + " at the end of the step", b[i], atEnd, SATISFIED);
        }
    }

    /**
     * The other two of the four conditions that make the continuous extension
     * an interpolant: it leaves the start of the step along the first stage and
     * arrives at the end along the last one, which is the derivative of the
     * solution at each of the two points.
     */
    @Test
    public void testTheContinuousExtensionLeavesAndArrivesAlongTheSolution() {
        ButcherTableau t = ButcherTableau.DORMAND_PRINCE_45;
        double[][] d = t.dense();
        int last = d.length - 1;
        for (int i = 0; i < d.length; ++i) {
            double slopeAtStart = d[i][0];
            double slopeAtEnd = 0.0;
            for (int j = 0; j < d[i].length; ++j) {
                slopeAtEnd += (j + 1) * d[i][j];
            }
            // exact, because the linear coefficient is written as one or zero
            assertEquals("stage " + i + " at the start", (i == 0) ? 1.0 : 0.0, slopeAtStart, 0.0);
            assertEquals("stage " + i + " at the end", (i == last) ? 1.0 : 0.0, slopeAtEnd, SATISFIED);
        }
    }

    @Test
    public void testWhatTheTableauxSayAboutThemselves() {
        ButcherTableau rk4 = ButcherTableau.CLASSIC_RK4;
        assertEquals(4, rk4.order());
        assertEquals(4, rk4.stages());
        assertEquals(4, rk4.denseStages());
        assertFalse(rk4.hasEmbedded());
        assertEquals(0, rk4.embeddedOrder());
        assertNull(rk4.bStar());
        assertFalse(rk4.hasDenseOutput());
        assertEquals(0, rk4.denseDegree());
        assertNull(rk4.dense());

        ButcherTableau dp = ButcherTableau.DORMAND_PRINCE_45;
        assertEquals(5, dp.order());
        assertEquals(7, dp.stages());
        assertEquals(7, dp.denseStages());
        assertTrue(dp.hasEmbedded());
        assertEquals(4, dp.embeddedOrder());
        assertNotNull(dp.bStar());
        assertTrue(dp.hasDenseOutput());
        assertEquals(4, dp.denseDegree());
        assertNotNull(dp.dense());
    }

    @Test
    public void testTheAccessorsHandOutCopies() {
        ButcherTableau t = ButcherTableau.DORMAND_PRINCE_45;
        double[] b = t.b();
        double keep = b[0];
        b[0] = 42.0;
        assertEquals(keep, t.b()[0], 0.0);

        double[][] a = t.a();
        double keptA = a[1][0];
        a[1][0] = 42.0;
        assertEquals(keptA, t.a()[1][0], 0.0);

        double[][] d = t.dense();
        double keptD = d[0][0];
        d[0][0] = 42.0;
        assertEquals(keptD, t.dense()[0][0], 0.0);
    }

    @Test
    public void testARowThatDoesNotSumToItsStageTimeIsRefused() {
        double[] c = { 0.0, 0.5 };
        double[][] a = { {}, { 0.4 } };
        double[] b = { 0.0, 1.0 };
        try {
            new ButcherTableau("wrong", 2, c, a, b);
            fail("a row summing to 0.4 against a stage time of 0.5 must not be accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("row 1"));
        }
    }

    @Test
    public void testAContinuousExtensionThatMissesTheEndOfTheStepIsRefused() {
        double[] c = { 0.0, 0.5 };
        double[][] a = { {}, { 0.5 } };
        double[] b = { 0.0, 1.0 };
        double[][] dense = { { 0.0, 0.0 }, { 1.0, 0.5 } };
        try {
            new ButcherTableau("wrong", 2, c, a, b, null, 0, dense);
            fail("an extension reaching 1.5 where b is 1.0 must not be accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("continuous extension"));
        }
    }

    @Test
    public void testTheShapesAreChecked() {
        double[] c = { 0.0, 0.5 };
        double[][] a = { {}, { 0.5 } };
        double[] b = { 0.0, 1.0 };
        refuse(null, 2, c, a, b, "name");
        refuse("x", 0, c, a, b, "order");
        refuse("x", 2, c, new double[][] { {} }, b, "rows");
        refuse("x", 2, c, new double[][] { {}, {} }, b, "length 1");
        refuse("x", 2, c, a, new double[] { 1.0, 0.0, 0.0 }, "length 3");
    }

    @Test
    public void testTwoTableauxWithTheSameCoefficientsAreEqual() {
        ButcherTableau t = ButcherTableau.DORMAND_PRINCE_45;
        ButcherTableau same = new ButcherTableau("under another name", t.order(), t.c(), t.a(), t.b(),
                t.bStar(), t.embeddedOrder(), t.dense());
        assertEquals(t, same);
        assertEquals(t.hashCode(), same.hashCode());
        assertFalse(t.equals(ButcherTableau.CLASSIC_RK4));
        assertTrue(t.toString().contains("FSAL"));
    }

    private static void refuse(String name, int order, double[] c, double[][] a, double[] b, String hint) {
        try {
            new ButcherTableau(name, order, c, a, b);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }

    private static void assertSatisfied(double[] residuals, int from, int to) {
        for (int i = from; i < to; ++i) {
            assertEquals("order condition " + (i + 1), 0.0, residuals[i], SATISFIED);
        }
    }

    private static void assertViolated(double[] residuals, int from, int to) {
        for (int i = from; i < to; ++i) {
            assertTrue("order condition " + (i + 1) + " is satisfied at " + residuals[i]
                    + " where it must not be", Math.abs(residuals[i]) > VIOLATED);
        }
    }

    private static double worstRowSum(ButcherTableau t) {
        double[] c = t.c();
        double[][] a = t.a();
        double worst = 0.0;
        for (int i = 0; i < a.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += a[i][j];
            }
            worst = Math.max(worst, Math.abs(sum - c[i]));
        }
        return worst;
    }

    /**
     * The seventeen order conditions in the standard ordering, as the amount by
     * which each is missed: the first satisfies order one, the first two order
     * two, the first four order three, the first eight order four and all
     * seventeen order five.
     */
    private static double[] conditions(ButcherTableau t, double[] w) {
        double[] c = t.c();
        double[][] a = t.a();
        double[] ac = apply(a, c);
        double[] ac2 = apply(a, power(c, 2));
        double[] ac3 = apply(a, power(c, 3));
        double[] aac = apply(a, ac);
        double[] acac = apply(a, product(c, ac));
        double[] aac2 = apply(a, ac2);
        double[] aaac = apply(a, aac);
        double[] ones = new double[c.length];
        for (int i = 0; i < ones.length; ++i) {
            ones[i] = 1.0;
        }
        return new double[] {
                dot(w, ones) - 1.0,
                dot(w, c) - 1.0 / 2.0,
                dot(w, power(c, 2)) - 1.0 / 3.0,
                dot(w, ac) - 1.0 / 6.0,
                dot(w, power(c, 3)) - 1.0 / 4.0,
                dot(w, product(c, ac)) - 1.0 / 8.0,
                dot(w, ac2) - 1.0 / 12.0,
                dot(w, aac) - 1.0 / 24.0,
                dot(w, power(c, 4)) - 1.0 / 5.0,
                dot(w, product(power(c, 2), ac)) - 1.0 / 10.0,
                dot(w, product(ac, ac)) - 1.0 / 20.0,
                dot(w, product(c, ac2)) - 1.0 / 15.0,
                dot(w, ac3) - 1.0 / 20.0,
                dot(w, product(c, aac)) - 1.0 / 30.0,
                dot(w, acac) - 1.0 / 40.0,
                dot(w, aac2) - 1.0 / 60.0,
                dot(w, aaac) - 1.0 / 120.0 };
    }

    private static double[] apply(double[][] a, double[] v) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += a[i][j] * v[j];
            }
            out[i] = sum;
        }
        return out;
    }

    private static double[] power(double[] v, int p) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; ++i) {
            out[i] = Math.pow(v[i], p);
        }
        return out;
    }

    private static double[] product(double[] x, double[] y) {
        double[] out = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            out[i] = x[i] * y[i];
        }
        return out;
    }

    private static double dot(double[] w, double[] x) {
        double sum = 0.0;
        for (int i = 0; i < w.length; ++i) {
            sum += w[i] * x[i];
        }
        return sum;
    }
}
