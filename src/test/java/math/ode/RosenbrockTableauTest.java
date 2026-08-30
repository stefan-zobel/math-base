package math.ode;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * The coefficients, and the shapes they have to have.
 * <p>
 * Nothing here checks an order condition. Those are awkward to write down for a
 * Rosenbrock method and the check that matters is available elsewhere and for
 * free: {@link RosenbrockTest} runs the method on {@code y' = lambda y}, where
 * one step is exactly {@code R(z)} with {@code R} rational, and reads the
 * order, the error constant and L-stability off that. What is asserted here is
 * that {@link RosenbrockTableau#ROS2} is the closed form it claims to be rather
 * than a decimal someone typed, and that the constructor refuses a table that
 * does not fit together.
 */
public final class RosenbrockTableauTest {

    private static final double ROOT = Math.sqrt(2.0);
    private static final double STEP = 2.0 - ROOT;

    @Test
    public void testRos2IsTheClosedFormAndNotATranscription() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        assertEquals("gamma is one plus one over the root of two", 1.0 + 1.0 / ROOT, t.gamma(), 0.0);
        assertEquals(STEP, t.a()[1][0], 0.0);
        assertEquals(-2.0 * STEP, t.c()[1][0], 0.0);
        assertArrayEquals(new double[] { t.gamma(), -t.gamma() }, t.d(), 0.0);
        assertArrayEquals(new double[] { 1.5 * STEP, 0.5 * STEP }, t.b(), 0.0);
        assertArrayEquals(new double[] { 0.5 * STEP, 0.5 * STEP }, t.bError(), 0.0);
        assertArrayEquals(new double[] { 0.0, 1.0 }, t.alpha(), 0.0);
    }

    /**
     * The same numbers as they are published, to fifteen places. This is the
     * one place the decimals appear, and it is a test rather than the code.
     */
    @Test
    public void testTheClosedFormAgreesWithThePublishedDecimals() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        assertEquals(1.7071067811865475, t.gamma(), 1.0e-15);
        assertEquals(0.585786437626905, t.a()[1][0], 1.0e-15);
        assertEquals(-1.17157287525381, t.c()[1][0], 1.0e-14);
        assertEquals(0.8786796564403574, t.b()[0], 1.0e-15);
        assertEquals(0.2928932188134525, t.b()[1], 1.0e-15);
        // 1 - 1/sqrt(2) and (2 - sqrt(2))/2 are the same number and not the
        // same double: they part company two bits from the end
        assertEquals("and the second weight is also one minus one over the root of two",
                1.0 - 1.0 / ROOT, t.b()[1], 5.0e-16);
    }

    @Test
    public void testTheShapeIsWhatTheStepperExpects() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        assertEquals("ROS2", t.name());
        assertEquals(2, t.stages());
        assertEquals(2, t.order());
        assertEquals(1, t.embeddedOrder());
        assertTrue(t.hasErrorEstimate());
        assertFalse("two stages carry no continuous extension", t.hasDenseOutput());
        assertEquals(0, t.denseRows());
        assertNull(t.dense());
        assertFalse("the last stage is not the solution", t.isStifflyAccurate());
        assertEquals("the stage matrix is ragged", 0, t.a()[0].length);
        assertEquals(1, t.a()[1].length);
        assertTrue(t.toString(), t.toString().contains("ROS2"));
        assertTrue(t.toString(), t.toString().contains("2 stages"));
    }

    @Test
    public void testEveryAccessorHandsOutACopy() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        t.alpha()[1] = 99.0;
        t.a()[1][0] = 99.0;
        t.c()[1][0] = 99.0;
        t.d()[0] = 99.0;
        t.b()[0] = 99.0;
        t.bError()[0] = 99.0;
        assertEquals(1.0, t.alpha()[1], 0.0);
        assertEquals(STEP, t.a()[1][0], 0.0);
        assertEquals(-2.0 * STEP, t.c()[1][0], 0.0);
        assertEquals(t.gamma(), t.d()[0], 0.0);
        assertEquals(1.5 * STEP, t.b()[0], 0.0);
        assertEquals(0.5 * STEP, t.bError()[0], 0.0);
    }

    @Test
    public void testATableauEqualsOneBuiltFromTheSameNumbers() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        RosenbrockTableau same = new RosenbrockTableau(t.name(), t.order(), t.embeddedOrder(), t.gamma(),
                t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError(), t.dense());
        assertEquals(t, same);
        assertEquals(t.hashCode(), same.hashCode());
        RosenbrockTableau other = new RosenbrockTableau("other", t.order(), t.embeddedOrder(), t.gamma(),
                t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError(), t.dense());
        assertNotEquals(t, other);
        assertNotEquals(t, null);
        assertNotEquals(t, "not a tableau");
    }

    @Test
    public void testADenseTableauReportsItsRows() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        RosenbrockTableau dense = new RosenbrockTableau("with dense", t.order(), t.embeddedOrder(),
                t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError(),
                new double[][] { { 1.0, 2.0 }, { 3.0, 4.0 } });
        assertTrue(dense.hasDenseOutput());
        assertEquals(2, dense.denseRows());
        assertEquals(4.0, dense.dense()[1][1], 0.0);
        dense.dense()[1][1] = 99.0;
        assertEquals("and that is a copy too", 4.0, dense.dense()[1][1], 0.0);
    }

    @Test
    public void testATableauWithoutAnErrorEstimateIsAllowed() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        RosenbrockTableau plain = new RosenbrockTableau("fixed step only", 2, 0, t.gamma(), t.alpha(),
                t.a(), t.c(), t.d(), t.b(), null, null);
        assertFalse(plain.hasErrorEstimate());
        assertNull(plain.bError());
        assertEquals(0, plain.embeddedOrder());
        assertTrue(plain.toString(), plain.toString().contains("order 2,"));
    }

    @Test
    public void testTheConstructorRefusesATableThatDoesNotFitTogether() {
        RosenbrockTableau t = RosenbrockTableau.ROS2;
        refuses("name", null, 2, 1, t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError());
        refuses("order", "x", 0, 1, t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError());
        refuses("gamma", "x", 2, 1, 0.0, t.alpha(), t.a(), t.c(), t.d(), t.b(), t.bError());
        refuses("alpha[0] must be zero", "x", 2, 1, t.gamma(), new double[] { 0.5, 1.0 }, t.a(), t.c(),
                t.d(), t.b(), t.bError());
        refuses("a[1] must hold 1", "x", 2, 1, t.gamma(), t.alpha(), new double[][] { {}, {} }, t.c(),
                t.d(), t.b(), t.bError());
        refuses("c must have 2 rows", "x", 2, 1, t.gamma(), t.alpha(), t.a(), new double[][] { {} },
                t.d(), t.b(), t.bError());
        refuses("d must hold 2", "x", 2, 1, t.gamma(), t.alpha(), t.a(), t.c(), new double[] { 1.0 },
                t.b(), t.bError());
        refuses("not finite", "x", 2, 1, t.gamma(), t.alpha(), t.a(), t.c(), t.d(),
                new double[] { Double.NaN, 1.0 }, t.bError());
        refuses("embeddedOrder is 1 without", "x", 2, 1, t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(),
                null);
        refuses("without an embeddedOrder", "x", 2, 0, t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(),
                t.bError());
        try {
            new RosenbrockTableau("x", 2, 1, t.gamma(), t.alpha(), t.a(), t.c(), t.d(), t.b(),
                    t.bError(), new double[][] { { 1.0 } });
            fail("expected a refusal about the dense row");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("dense[0]"));
        }
    }

    private static void refuses(String fragment, String name, int order, int embeddedOrder, double gamma,
            double[] alpha, double[][] a, double[][] c, double[] d, double[] b, double[] bError) {
        try {
            new RosenbrockTableau(name, order, embeddedOrder, gamma, alpha, a, c, d, b, bError, null);
            fail("expected a refusal mentioning " + fragment);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(fragment));
        }
    }
}
