package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The coefficients, checked against the conditions every Nystroem method
 * satisfies whatever its order.
 * <p>
 * This is where a transcription error gets caught, and it is not a hypothetical
 * concern: the first transcription of {@link NystromTableau#RKN6_4} read a digit
 * group wrong and produced a method of order two that ran, converged and
 * returned plausible numbers. The row-sum condition below points straight at the
 * row that is wrong.
 */
public class NystromTableauTest {

    private static final double EPS = 1.0e-12;

    @Test
    public void testTheRowsOfTheStageMatrixSumToHalfTheSquaredStageTimes() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[] c = t.c();
        double[][] a = t.a();
        for (int i = 0; i < c.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < a[i].length; ++j) {
                sum += a[i][j];
            }
            assertEquals("row " + i, 0.5 * c[i] * c[i], sum, EPS);
        }
    }

    @Test
    public void testTheWeightsSumToWhatTheyMust() {
        NystromTableau t = NystromTableau.RKN6_4;
        assertEquals("sum of bbar", 0.5, sum(t.bbar()), EPS);
        assertEquals("sum of b", 1.0, sum(t.b()), EPS);
        assertEquals("sum of bbarStar", 0.5, sum(t.bbarStar()), EPS);
        assertEquals("sum of bStar", 1.0, sum(t.bStar()), EPS);
    }

    /**
     * The first stage sits at the start of the step and the last at its end,
     * which is what lets the step after it reuse the evaluation.
     */
    @Test
    public void testTheFirstSameAsLastPropertyIsReadOffTheCoefficients() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[] c = t.c();
        double[][] a = t.a();
        double[] bbar = t.bbar();
        assertEquals("the first stage time", 0.0, c[0], 0.0);
        assertEquals("the last stage time", 1.0, c[c.length - 1], 0.0);
        for (int j = 0; j < c.length - 1; ++j) {
            assertEquals("a[last][" + j + "] against bbar[" + j + "]", bbar[j],
                    a[c.length - 1][j], 0.0);
        }
        assertEquals("the last position weight", 0.0, bbar[bbar.length - 1], 0.0);
        assertTrue("so the tableau says it is FSAL", t.isFsal());
    }

    /**
     * A continuous extension has to meet the step it extends: at the end of the
     * step each stage's polynomial must come to the weight that stage carries.
     */
    @Test
    public void testTheContinuousExtensionMeetsTheEndOfTheStep() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[][] p = t.densePosition();
        double[][] q = t.denseVelocity();
        double[] bbar = t.bbar();
        double[] b = t.b();
        for (int i = 0; i < bbar.length; ++i) {
            assertEquals("position polynomial " + i + " at theta = 1", bbar[i], sum(p[i]), EPS);
            assertEquals("velocity polynomial " + i + " at theta = 1", b[i], sum(q[i]), EPS);
        }
    }

    @Test
    public void testWhatTheTableauSaysAboutItself() {
        NystromTableau t = NystromTableau.RKN6_4;
        assertEquals("name", "RKN6(4)6FM", t.name());
        assertEquals("order", 6, t.order());
        assertEquals("embedded order", 4, t.embeddedOrder());
        assertEquals("stages", 6, t.stages());
        assertTrue("has an error estimate", t.hasErrorEstimate());
        assertTrue("has dense output", t.hasDenseOutput());
        assertTrue("says so", t.toString().contains("RKN6(4)6FM"));
        assertTrue("and names both orders", t.toString().contains("order 6(4)"));
    }

    @Test
    public void testEveryAccessorHandsOutACopy() {
        NystromTableau t = NystromTableau.RKN6_4;
        assertNotSame("c", t.c(), t.c());
        assertNotSame("a", t.a(), t.a());
        assertNotSame("bbar", t.bbar(), t.bbar());
        assertNotSame("b", t.b(), t.b());
        assertNotSame("bbarStar", t.bbarStar(), t.bbarStar());
        assertNotSame("bStar", t.bStar(), t.bStar());
        assertNotSame("densePosition", t.densePosition(), t.densePosition());
        assertNotSame("denseVelocity", t.denseVelocity(), t.denseVelocity());
        double[] c = t.c();
        c[0] = 17.0;
        assertEquals("the tableau is unmoved", 0.0, t.c()[0], 0.0);
    }

    /**
     * The condition that caught the real transcription error, which was a
     * denominator with two zeros too many and so a row out by 37 %.
     */
    @Test
    public void testARowThatDoesNotSumToItsStageTimeIsRefused() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[][] a = t.a();
        a[3][0] *= 100.0;
        try {
            new NystromTableau("broken", 6, 4, t.c(), a, t.bbar(), t.b(), t.bbarStar(), t.bStar(),
                    null, null);
            fail("a row that does not sum to half its squared stage time was accepted");
        } catch (IllegalArgumentException e) {
            assertTrue("names the row, not " + e.getMessage(), e.getMessage().contains("row 3"));
        }
    }

    @Test
    public void testWeightsThatDoNotSumToTheirTotalsAreRefused() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[] bbar = t.bbar();
        bbar[0] += 0.25;
        try {
            new NystromTableau("broken", 6, 4, t.c(), t.a(), bbar, t.b(), null, null, null, null);
            fail("weights that do not sum to one half were accepted");
        } catch (IllegalArgumentException e) {
            assertTrue("says which sum, not " + e.getMessage(),
                    e.getMessage().contains("the sum of bbar"));
        }
    }

    @Test
    public void testAContinuousExtensionThatMissesTheEndOfTheStepIsRefused() {
        NystromTableau t = NystromTableau.RKN6_4;
        double[][] dense = t.densePosition();
        dense[0][0] += 0.5;
        try {
            new NystromTableau("broken", 6, 4, t.c(), t.a(), t.bbar(), t.b(), t.bbarStar(),
                    t.bStar(), dense, t.denseVelocity());
            fail("a continuous extension that misses the end of the step was accepted");
        } catch (IllegalArgumentException e) {
            assertTrue("says so, not " + e.getMessage(),
                    e.getMessage().contains("does not meet the step it extends"));
        }
    }

    @Test
    public void testTheShapesAreChecked() {
        NystromTableau t = NystromTableau.RKN6_4;
        try {
            new NystromTableau("short", 6, 4, t.c(), t.a(), new double[3], t.b(), null, null, null,
                    null);
            fail("a weight vector of the wrong length was accepted");
        } catch (IllegalArgumentException expected) {
            // the message names the stage count
        }
        try {
            new NystromTableau("half embedded", 6, 4, t.c(), t.a(), t.bbar(), t.b(), t.bbarStar(),
                    null, null, null);
            fail("one embedded weight vector without the other was accepted");
        } catch (IllegalArgumentException expected) {
            // an embedded solution needs both
        }
        try {
            double[] c = t.c();
            c[0] = 0.5;
            new NystromTableau("shifted", 6, 4, c, t.a(), t.bbar(), t.b(), null, null, null, null);
            fail("a first stage away from the start of the step was accepted");
        } catch (IllegalArgumentException expected) {
            // c[0] has to be zero
        }
    }

    @Test
    public void testATableauWithoutAnErrorEstimateIsAllowedButSaysSo() {
        NystromTableau t = NystromTableau.RKN6_4;
        NystromTableau plain = new NystromTableau("no estimate", 6, 0, t.c(), t.a(), t.bbar(),
                t.b(), null, null, null, null);
        assertFalse("has no error estimate", plain.hasErrorEstimate());
        assertFalse("and no dense output", plain.hasDenseOutput());
        assertEquals("its embedded order is zero", 0, plain.embeddedOrder());
        assertEquals("bbarStar is null", null, plain.bbarStar());
    }

    @Test
    public void testTwoTableauxWithTheSameCoefficientsAreEqual() {
        NystromTableau t = NystromTableau.RKN6_4;
        NystromTableau same = new NystromTableau("RKN6(4)6FM", 6, 4, t.c(), t.a(), t.bbar(), t.b(),
                t.bbarStar(), t.bStar(), t.densePosition(), t.denseVelocity());
        assertEquals("equal", t, same);
        assertEquals("and hash alike", t.hashCode(), same.hashCode());
    }

    private static double sum(double[] x) {
        double s = 0.0;
        for (int i = 0; i < x.length; ++i) {
            s += x[i];
        }
        return s;
    }
}
