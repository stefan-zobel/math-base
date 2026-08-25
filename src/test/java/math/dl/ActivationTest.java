package math.dl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.function.DoubleUnaryOperator;

import org.junit.Test;

import math.MathConsts;

/**
 * The four activation functions of the package against references built from
 * {@link Math} alone, plus the invariants that hold whatever the reference
 * says: continuity, agreement between the analytic derivative and a central
 * difference, agreement between the {@code float} entry points and the
 * {@code double} ones, and the behavior at the edges of the range.
 * <p>
 * Every tolerance here was measured before it was written down; the value each
 * one was measured at is in the comment beside it.
 */
public class ActivationTest {

    private static final double K = MathConsts.SQRT_TWO_BY_PI;

    /**
     * Relative part of the comparison. Worst measured over 2 000 001 samples
     * per band: values 2.280e-13, derivatives 9.420e-11.
     */
    private static final double REL_TOL = 1e-12;

    /**
     * Absolute part, without which the three derivatives cannot be compared at
     * all: each of them crosses zero -- GELU at -0.7525, Mish at -1.1924,
     * Swish at -1.2785 -- and there the two terms of the sum cancel, so a
     * relative measure reports 1e-10 for a difference of 3e-15. Worst measured
     * absolute difference over the same samples: 8.438e-15.
     */
    private static final double ABS_TOL = 1e-13;

    /** Worst measured against a central difference: 4.001e-11. */
    private static final double NUMERIC_TOL = 1e-9;

    /** Worst measured: 5.960e-08, which is half an ulp of a float. */
    private static final double FLOAT_TOL = 1e-7;

    /**
     * Where each value becomes exactly {@code x}, and each derivative exactly
     * 1.0, found by bisection. {@code MISH_SATURATES} is {@code 27.5 * ln(2)},
     * the same constant the tangent saturates at, because the softplus Mish
     * feeds it approaches its own argument.
     */
    private static final double MISH_SATURATES = 19.061547465398494;
    private static final double MISH_D_SATURATES = 20.573551;
    private static final double SWISH_SATURATES = 36.736801;
    private static final double SWISH_D_SATURATES = 40.436534;

    /** Where each one reaches zero: Math.exp underflows resp. overflows there. */
    private static final double MISH_VANISHES = -745.13;
    private static final double SWISH_VANISHES = -709.78;

    // ----- the four, each with a reference and the band it is valid on -------

    /** What a float entry point looks like; there is no primitive interface for it. */
    private interface FloatOp {
        float apply(float x);
    }

    private static final class Act {
        final String name;
        final DoubleUnaryOperator value;
        final DoubleUnaryOperator deriv;
        final DoubleUnaryOperator valueRef;
        final DoubleUnaryOperator derivRef;
        final FloatOp valueF;
        final FloatOp derivF;
        /** The band on which the implementation claims to be accurate. */
        final double lo;
        final double hi;

        Act(String name, DoubleUnaryOperator value, DoubleUnaryOperator deriv, DoubleUnaryOperator valueRef,
                DoubleUnaryOperator derivRef, FloatOp valueF, FloatOp derivF, double lo, double hi) {
            this.name = name;
            this.value = value;
            this.deriv = deriv;
            this.valueRef = valueRef;
            this.derivRef = derivRef;
            this.valueF = valueF;
            this.derivF = derivF;
            this.lo = lo;
            this.hi = hi;
        }
    }

    private static double sigma(double x) {
        return 1.0 / (1.0 + Math.exp(-x));
    }

    private static Act[] all() {
        return new Act[] {
                new Act("GELU", GELU::gelu, GELU::dgelu_dx, ActivationTest::geluRef, ActivationTest::dgeluRef,
                        GELU::geluF, GELU::dgeluF_dx, -25.0, 10.0),
                // the upper bound of each band is the point at which the value
                // becomes exactly x, so the whole accurate range is covered and
                // nothing beyond it
                new Act("Mish", Mish::mish, Mish::dmish_dx, ActivationTest::mishRef, ActivationTest::dmishRef,
                        Mish::mishF, Mish::dmishF_dx, -40.0, MISH_SATURATES),
                new Act("Swish", Swish::swish, Swish::dswish_dx, ActivationTest::swishRef,
                        ActivationTest::dswishRef, Swish::swishF, Swish::dswishF_dx, -40.0, SWISH_SATURATES),
                new Act("RELU", RELU::relu, RELU::drelu_dx, ActivationTest::reluRef, ActivationTest::dreluRef,
                        RELU::reluF, RELU::dreluF_dx, -20.0, 20.0) };
    }

    private static double geluRef(double x) {
        // sigma(2y) is 0.5 * (1 + tanh(y)) without the cancellation
        return x * sigma(2.0 * K * (x + 0.044715 * x * x * x));
    }

    private static double dgeluRef(double x) {
        double x2 = x * x;
        double y = K * (x + 0.044715 * x2 * x);
        double c = Math.cosh(y);
        return sigma(2.0 * y) + 0.5 * x / (c * c) * K * (1.0 + 0.134145 * x2);
    }

    private static double mishRef(double x) {
        return x * Math.tanh(Math.log1p(Math.exp(x)));
    }

    private static double dmishRef(double x) {
        double sp = Math.log1p(Math.exp(x));
        double c = Math.cosh(sp);
        return Math.tanh(sp) + x * sigma(x) / (c * c);
    }

    private static double swishRef(double x) {
        return x * sigma(x);
    }

    private static double dswishRef(double x) {
        // 1 - sigma(x) is sigma(-x) exactly, so this does not cancel
        return sigma(x) + x * sigma(x) * sigma(-x);
    }

    private static double reluRef(double x) {
        return x > 0.0 ? x : 0.0;
    }

    private static double dreluRef(double x) {
        return x > 0.0 ? 1.0 : 0.0;
    }

    // ----- against the references --------------------------------------------

    @Test
    public void theValuesAgreeWithAReferenceBuiltFromMath() {
        Act[] acts = all();
        for (int a = 0; a < acts.length; ++a) {
            Act act = acts[a];
            for (int i = 0; i <= 400000; ++i) {
                double x = act.lo + (act.hi - act.lo) * i / 400000.0;
                double got = act.value.applyAsDouble(x);
                double ref = act.valueRef.applyAsDouble(x);
                assertClose(act.name + ".value at " + x, ref, got);
            }
        }
    }

    @Test
    public void theDerivativesAgreeWithAReferenceBuiltFromMath() {
        Act[] acts = all();
        for (int a = 0; a < acts.length; ++a) {
            Act act = acts[a];
            for (int i = 0; i <= 400000; ++i) {
                double x = act.lo + (act.hi - act.lo) * i / 400000.0;
                double got = act.deriv.applyAsDouble(x);
                double ref = act.derivRef.applyAsDouble(x);
                assertClose(act.name + ".derivative at " + x, ref, got);
            }
        }
    }

    @Test
    public void theDerivativesAgreeWithACentralDifference() {
        // an oracle that knows nothing about how either side is computed
        Act[] acts = all();
        for (int a = 0; a < acts.length; ++a) {
            Act act = acts[a];
            if ("RELU".equals(act.name)) {
                continue; // its derivative is a step, a difference quotient says nothing
            }
            for (int i = 0; i <= 20000; ++i) {
                double x = act.lo + 1.0 + (act.hi - act.lo - 2.0) * i / 20000.0;
                double h = 1e-5 * Math.max(1.0, Math.abs(x));
                double num = (act.value.applyAsDouble(x + h) - act.value.applyAsDouble(x - h)) / (2.0 * h);
                double ana = act.deriv.applyAsDouble(x);
                assertTrue(act.name + " derivative at " + x + ": analytic " + ana + ", difference quotient "
                        + num, Math.abs(ana - num) <= NUMERIC_TOL * Math.max(Math.abs(num), 1.0));
            }
        }
    }

    // ----- invariants that hold whatever the reference says ------------------

    @Test
    public void theValuesAreContinuous() {
        // sampled finely, the largest step of a smooth function is about its
        // sample spacing times the steepest slope, so the ratio stays small. A
        // cutoff that truncates a value the function has not reached yet shows
        // up as a ratio in the hundreds -- Mish measured 119.26 here before
        // its cutoff was moved, and 1.09 after
        assertStepRatio("GELU.gelu", GELU::gelu, -35.0, 10.0, 5.0);
        assertStepRatio("GELU.dgelu_dx", GELU::dgelu_dx, -25.0, 10.0, 5.0);
        assertStepRatio("Mish.mish", Mish::mish, -40.0, 25.0, 5.0);
        assertStepRatio("Mish.dmish_dx", Mish::dmish_dx, -40.0, 25.0, 5.0);
        assertStepRatio("Swish.swish", Swish::swish, -40.0, 45.0, 5.0);
        assertStepRatio("Swish.dswish_dx", Swish::dswish_dx, -40.0, 45.0, 5.0);
        assertStepRatio("RELU.relu", RELU::relu, -20.0, 20.0, 5.0);

        // a sweep cannot see a step smaller than its own sample spacing --
        // GELU used to jump by 5.8e-07 at -4.861 and none of the sweeps above
        // noticed. So look at the cutoffs themselves, which is where a step
        // would sit
        assertNoStepAt("GELU.gelu at -30", GELU::gelu, -30.0);
        assertNoStepAt("GELU.gelu at 7.09", GELU::gelu, 7.09);
        assertNoStepAt("GELU.dgelu_dx at -21.5", GELU::dgelu_dx, -21.5);
        assertNoStepAt("GELU.dgelu_dx at 7.45", GELU::dgelu_dx, 7.45);
        assertNoStepAt("Mish.mish, lower", Mish::mish, MISH_VANISHES);
        assertNoStepAt("Mish.mish, upper", Mish::mish, MISH_SATURATES);
        assertNoStepAt("Mish.dmish_dx, lower", Mish::dmish_dx, MISH_VANISHES);
        assertNoStepAt("Mish.dmish_dx, upper", Mish::dmish_dx, MISH_D_SATURATES);
        assertNoStepAt("Swish.swish, lower", Swish::swish, SWISH_VANISHES);
        assertNoStepAt("Swish.swish, upper", Swish::swish, SWISH_SATURATES);
        assertNoStepAt("Swish.dswish_dx, lower", Swish::dswish_dx, SWISH_VANISHES);
        assertNoStepAt("Swish.dswish_dx, upper", Swish::dswish_dx, SWISH_D_SATURATES);
    }

    @Test
    public void theFloatVariantsAreTheDoubleOnesRounded() {
        Act[] acts = all();
        long lcg = 31L;
        for (int i = 0; i < 200000; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            float x = (float) (40.0 * (((lcg >>> 11) * 0x1.0p-53) - 0.5));
            for (int a = 0; a < acts.length; ++a) {
                Act act = acts[a];
                assertFloatMatches(act.name + ".valueF at " + x, act.value.applyAsDouble(x),
                        act.valueF.apply(x));
                assertFloatMatches(act.name + ".derivF at " + x, act.deriv.applyAsDouble(x),
                        act.derivF.apply(x));
            }
        }
    }

    @Test
    public void theEdgesBehave() {
        Act[] acts = all();
        for (int a = 0; a < acts.length; ++a) {
            Act act = acts[a];
            String n = act.name;
            assertEquals(n + " at +Inf", Double.POSITIVE_INFINITY,
                    act.value.applyAsDouble(Double.POSITIVE_INFINITY), 0.0);
            assertEquals(n + " derivative at +Inf", 1.0, act.deriv.applyAsDouble(Double.POSITIVE_INFINITY), 0.0);
            assertEquals(n + " at -Inf", 0.0, Math.abs(act.value.applyAsDouble(Double.NEGATIVE_INFINITY)), 0.0);
            assertEquals(n + " derivative at -Inf", 0.0,
                    Math.abs(act.deriv.applyAsDouble(Double.NEGATIVE_INFINITY)), 0.0);
            assertEquals(n + " at 0", 0.0, Math.abs(act.value.applyAsDouble(0.0)), 0.0);
            assertEquals(n + " at -0.0", 0.0, Math.abs(act.value.applyAsDouble(-0.0)), 0.0);
            assertEquals(n + " at MAX_VALUE", Double.MAX_VALUE, act.value.applyAsDouble(Double.MAX_VALUE), 0.0);
            assertEquals(n + " at -MAX_VALUE", 0.0,
                    Math.abs(act.value.applyAsDouble(-Double.MAX_VALUE)), 0.0);
            assertEquals(n + " derivative at -MAX_VALUE", 0.0,
                    Math.abs(act.deriv.applyAsDouble(-Double.MAX_VALUE)), 0.0);
        }

        // the derivative at zero is the one value each of them has in closed form
        assertEquals("GELU'(0)", 0.5, GELU.dgelu_dx(0.0), 0.0);
        assertEquals("Swish'(0)", 0.5, Swish.dswish_dx(0.0), 0.0);
        assertEquals("Mish'(0)", Math.tanh(Math.log(2.0)), Mish.dmish_dx(0.0), 1e-15);
        assertEquals("RELU'(0)", 0.0, RELU.drelu_dx(0.0), 0.0);
    }

    @Test
    public void nanPropagates() {
        assertTrue("GELU", Double.isNaN(GELU.gelu(Double.NaN)));
        assertTrue("GELU derivative", Double.isNaN(GELU.dgelu_dx(Double.NaN)));
        assertTrue("Mish", Double.isNaN(Mish.mish(Double.NaN)));
        assertTrue("Mish derivative", Double.isNaN(Mish.dmish_dx(Double.NaN)));
        assertTrue("Swish", Double.isNaN(Swish.swish(Double.NaN)));
        assertTrue("Swish derivative", Double.isNaN(Swish.dswish_dx(Double.NaN)));

        // RELU used to answer zero for all four of these, because x > 0.0 is
        // false for a NaN and the else branch caught it. A NaN turned into a
        // plausible zero is how a defect upstream stops being visible
        assertTrue("RELU", Double.isNaN(RELU.relu(Double.NaN)));
        assertTrue("RELU derivative", Double.isNaN(RELU.drelu_dx(Double.NaN)));
        assertTrue("RELU float", Float.isNaN(RELU.reluF(Float.NaN)));
        assertTrue("RELU float derivative", Float.isNaN(RELU.dreluF_dx(Float.NaN)));

        // and it still clamps everything else the way it did
        assertEquals("RELU(-0.0)", 0.0, RELU.relu(-0.0), 0.0);
        assertEquals("RELU(MIN_VALUE)", Double.MIN_VALUE, RELU.relu(Double.MIN_VALUE), 0.0);
        assertEquals("RELU(-Inf)", 0.0, RELU.relu(Double.NEGATIVE_INFINITY), 0.0);
    }

    @Test
    public void theCutoffsSitWhereADoubleGivesUp() {
        // GELU
        assertTrue("GELU still resolves -20", GELU.gelu(-20.0) < 0.0);
        assertClose("GELU at -20", geluRef(-20.0), GELU.gelu(-20.0));
        assertEquals("GELU value is exactly x above 7.09", 8.0, GELU.gelu(8.0), 0.0);
        assertEquals("GELU derivative is exactly 1 above 7.45", 1.0, GELU.dgelu_dx(8.0), 0.0);
        assertEquals("GELU value underflows below -30", 0.0, Math.abs(GELU.gelu(-31.0)), 0.0);

        // Mish used to give up at -8, where the true value is still
        // -2.683e-03, and at 16, where mish(16) is 15.999999999999595 rather
        // than 16
        assertClose("Mish at -8.0000001", mishRef(-8.0000001), Mish.mish(-8.0000001));
        assertClose("Mish at -40", mishRef(-40.0), Mish.mish(-40.0));
        assertClose("Mish at -300", mishRef(-300.0), Mish.mish(-300.0));
        assertTrue("Mish still resolves -700", Mish.mish(-700.0) < 0.0);
        assertClose("Mish at 16", mishRef(16.0), Mish.mish(16.0));
        assertTrue("and 16 is not yet exactly x", Mish.mish(16.0) != 16.0);
        assertEquals("but 20 is", 20.0, Mish.mish(20.0), 0.0);
        assertEquals("below the vanishing point", 0.0, Math.abs(Mish.mish(-746.0)), 0.0);

        // Swish used to give up at -17, where the true value is still
        // -7.038e-07, and at 17, where swish(17) is 16.999999296210614
        assertClose("Swish at -17.00001", swishRef(-17.00001), Swish.swish(-17.00001));
        assertClose("Swish at -40", swishRef(-40.0), Swish.swish(-40.0));
        assertTrue("Swish still resolves -700", Swish.swish(-700.0) < 0.0);
        assertClose("Swish at 17", swishRef(17.0), Swish.swish(17.0));
        assertTrue("and 17 is not yet exactly x", Swish.swish(17.0) != 17.0);
        assertEquals("but 37 is", 37.0, Swish.swish(37.0), 0.0);
        assertEquals("below the vanishing point", 0.0, Math.abs(Swish.swish(-711.0)), 0.0);
    }

    // ----- helpers -----------------------------------------------------------

    private static void assertClose(String what, double expected, double actual) {
        assertTrue(what + " : expected " + expected + " but was " + actual,
                Math.abs(actual - expected) <= REL_TOL * Math.abs(expected) + ABS_TOL);
    }

    private static void assertFloatMatches(String what, double exact, float rounded) {
        double mag = Math.abs(exact);
        if (mag <= 1e-30 || mag >= 1e30) {
            return; // outside float range the cast can only underflow or overflow
        }
        assertTrue(what + " : double " + exact + ", float " + rounded,
                Math.abs(rounded - exact) <= FLOAT_TOL * mag);
    }

    /** No step at a cutoff means the two sides differ by about the slope times 2 eps. */
    private static void assertNoStepAt(String what, DoubleUnaryOperator f, double c) {
        double eps = 1e-9 * Math.max(1.0, Math.abs(c));
        double below = f.applyAsDouble(c - eps);
        double above = f.applyAsDouble(c + eps);
        double scale = Math.max(Math.abs(below), Math.abs(above));
        assertTrue(what + " : " + below + " on one side, " + above + " on the other",
                Math.abs(above - below) <= 1e-6 * scale + ABS_TOL);
    }

    private static void assertStepRatio(String what, DoubleUnaryOperator f, double lo, double hi, double bound) {
        int n = 400000;
        double spacing = (hi - lo) / n;
        double worst = 0.0;
        double at = 0.0;
        double prev = f.applyAsDouble(lo);
        for (int i = 1; i <= n; ++i) {
            double x = lo + (hi - lo) * i / (double) n;
            double v = f.applyAsDouble(x);
            double step = Math.abs(v - prev);
            if (step > worst) {
                worst = step;
                at = x;
            }
            prev = v;
        }
        double ratio = worst / spacing;
        assertTrue(what + " : largest step " + worst + " at x = " + at + " is " + ratio
                + " times the sample spacing, bound " + bound, ratio < bound);
    }
}
