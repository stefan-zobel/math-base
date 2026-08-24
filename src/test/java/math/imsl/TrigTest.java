package math.imsl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * {@link Trig} against {@link Math}, which is the independent implementation of
 * the same three functions. The routines are kept because they are two to three
 * times faster, so they have to be checked against the accurate ones rather
 * than replaced by them.
 */
public class TrigTest {

    /** Measured worst relative error of {@code tanh}: 4.405e-16, three ulp, on (0, 1]. */
    private static final double TANH_TOL = 1e-15;

    /** Measured worst relative error of {@code cosh}: 3.140e-16. */
    private static final double COSH_TOL = 1e-15;

    /** Measured worst relative error of {@code sech2}: 8.349e-16, it squares a cosh. */
    private static final double SECH2_TOL = 2e-15;

    /** The point where a double first rounds tanh to 1.0, {@code 0.5 * ln(8 / eps)}. */
    private static final double SATURATION = 19.061547465398498;

    // ----- against the JDK ---------------------------------------------------

    @Test
    public void tanhAgreesWithTheJdkOverItsWholeRange() {
        // the bands are the branches: the linear one, the Chebyshev one, the
        // exponential one, and the saturated one
        double[][] bands = { { 0.0, 1e-8 }, { 0.0, 1.0 }, { 1.0, 8.0 }, { 8.0, SATURATION },
                { SATURATION, 100.0 }, { 100.0, 800.0 }, { -800.0, -1e-12 }, { -20.0, 20.0 } };
        double worst = 0.0;
        double worstUlp = 0.0;
        for (int b = 0; b < bands.length; ++b) {
            for (int i = 0; i <= 200000; ++i) {
                double x = bands[b][0] + (bands[b][1] - bands[b][0]) * i / 200000.0;
                double ref = Math.tanh(x);
                double got = Trig.tanh(x);
                double diff = Math.abs(got - ref);
                if (ref != 0.0) {
                    double rel = diff / Math.abs(ref);
                    if (rel > worst) {
                        worst = rel;
                    }
                    double ulps = diff / Math.ulp(ref);
                    if (ulps > worstUlp) {
                        worstUlp = ulps;
                    }
                }
                assertTrue("tanh(" + x + ") = " + got + ", JDK says " + ref,
                        diff <= TANH_TOL * Math.max(Math.abs(ref), Double.MIN_NORMAL));
            }
        }
        assertTrue("worst relative error " + worst, worst < TANH_TOL);
        assertTrue("worst error in ulp " + worstUlp, worstUlp <= 4.0);
    }

    @Test
    public void coshAgreesWithTheJdk() {
        for (int i = 0; i <= 400000; ++i) {
            double x = -40.0 + 80.0 * i / 400000.0;
            double ref = Math.cosh(x);
            double got = Trig.cosh(x);
            assertTrue("cosh(" + x + ") = " + got + ", JDK says " + ref,
                    Math.abs(got - ref) <= COSH_TOL * ref);
        }
        assertEquals("an infinite argument", Double.POSITIVE_INFINITY, Trig.cosh(800.0), 0.0);
        assertEquals("an infinite argument", Double.POSITIVE_INFINITY, Trig.cosh(Double.NEGATIVE_INFINITY), 0.0);
        assertTrue("NaN", Double.isNaN(Trig.cosh(Double.NaN)));
    }

    @Test
    public void sech2AgreesWithTheJdk() {
        for (int i = 0; i <= 400000; ++i) {
            double x = -30.0 + 60.0 * i / 400000.0;
            double ref = 1.0 / (Math.cosh(x) * Math.cosh(x));
            double got = Trig.sech2(x);
            assertTrue("sech2(" + x + ") = " + got + ", JDK says " + ref,
                    Math.abs(got - ref) <= SECH2_TOL * Math.max(ref, Double.MIN_NORMAL));
        }
    }

    // ----- the defect this class was written for -----------------------------

    @Test
    public void tanhSaturatesOnlyWhereADoubleDoes() {
        // it used to saturate at 7.9772948850000001, a single precision cutoff,
        // where the true value is still 1 - 2.36e-7
        assertTrue("tanh(8) must not be 1", Trig.tanh(8.0) < 1.0);
        assertEquals("tanh(8) against the JDK", Math.tanh(8.0), Trig.tanh(8.0), 1e-15);
        assertEquals("the deficit at the old cutoff", 2.3552642847146643e-07,
                1.0 - Trig.tanh(7.9772948850000001), 1e-15);

        // the exponential form itself reaches exactly 1.0 at
        // 0.5 * ln(4 / eps) = 18.714973875118524, one ulp before a double
        // rounds tanh to 1.0, so the constant above that is an overflow guard
        // for Math.exp rather than the place the value saturates
        assertTrue("not saturated below 18.7", Trig.tanh(18.5) < 1.0);
        assertEquals("and exact against the JDK there", Math.tanh(18.5), Trig.tanh(18.5), 0.0);
        assertEquals("at saturation", 1.0, Trig.tanh(SATURATION), 0.0);
        assertEquals("beyond saturation", 1.0, Trig.tanh(20.0), 0.0);

        // Math.exp overflows past 709.78 and (Inf - 0) / (Inf + 0) is NaN, so
        // the saturating branch has to keep catching these
        assertEquals("past the exp overflow", 1.0, Trig.tanh(710.0), 0.0);
        assertEquals("past the exp overflow", 1.0, Trig.tanh(800.0), 0.0);
        assertEquals("past the exp overflow", -1.0, Trig.tanh(-800.0), 0.0);
        assertEquals("Double.MAX_VALUE", 1.0, Trig.tanh(Double.MAX_VALUE), 0.0);
    }

    @Test
    public void theSquaredSecantOutlivesTheIdentity() {
        // 1 - tanh(x)^2 does not fail at a point, it erodes: the subtraction
        // cancels harder the closer tanh gets to 1. Measured against sech2 it
        // is 3.7e-15 off at x = 2, 2.2e-10 at x = 8, 1.8e-06 at x = 12,
        // 5.0e-03 at x = 16, and exactly zero from x = 18.715 on. sech2 goes
        // through cosh and stays accurate throughout. Before the cutoff was
        // fixed the identity died at x = 7.978 already, where the true value
        // is 4.7e-07
        assertTrue("the identity survives x = 8 now", 1.0 - Trig.tanh(8.0) * Trig.tanh(8.0) > 0.0);
        assertEquals("and is still roughly right there", Trig.sech2(8.0),
                1.0 - Trig.tanh(8.0) * Trig.tanh(8.0), 1e-9 * Trig.sech2(8.0));

        // where the identity has lost most of its digits, sech2 has lost none
        assertEquals("sech2 at 16", 1.0 / (Math.cosh(16.0) * Math.cosh(16.0)), Trig.sech2(16.0),
                SECH2_TOL * Trig.sech2(16.0));

        // and eventually the identity is gone altogether, which is why sech2
        // exists at all
        assertEquals("the identity is gone at 20", 0.0, 1.0 - Trig.tanh(20.0) * Trig.tanh(20.0), 0.0);
        assertTrue("sech2 is not", Trig.sech2(20.0) > 0.0);
        assertEquals("and is still right", 1.0 / (Math.cosh(20.0) * Math.cosh(20.0)), Trig.sech2(20.0),
                SECH2_TOL * Trig.sech2(20.0));
    }

    // ----- shape and edges ---------------------------------------------------

    @Test
    public void tanhIsOddBitForBit() {
        long lcg = 7L;
        for (int i = 0; i < 500000; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double x = 60.0 * (((lcg >>> 11) * 0x1.0p-53) - 0.5);
            assertSameBits("tanh(-x) at x = " + x, -Trig.tanh(x), Trig.tanh(-x));
        }
    }

    @Test
    public void theEdgeCasesAreWhatTheyShouldBe() {
        assertSameBits("tanh(0)", 0.0, Trig.tanh(0.0));
        assertSameBits("tanh(-0.0) keeps its sign", -0.0, Trig.tanh(-0.0));
        assertSameBits("the smallest subnormal passes through", Double.MIN_VALUE,
                Trig.tanh(Double.MIN_VALUE));
        assertEquals("tanh(+Inf)", 1.0, Trig.tanh(Double.POSITIVE_INFINITY), 0.0);
        assertEquals("tanh(-Inf)", -1.0, Trig.tanh(Double.NEGATIVE_INFINITY), 0.0);
        assertTrue("tanh(NaN)", Double.isNaN(Trig.tanh(Double.NaN)));
        assertTrue("sech2(NaN)", Double.isNaN(Trig.sech2(Double.NaN)));
        assertEquals("sech2(+Inf)", 0.0, Trig.sech2(Double.POSITIVE_INFINITY), 0.0);
        assertEquals("cosh(0)", 1.0, Trig.cosh(0.0), 0.0);
        assertEquals("sech2(0)", 1.0, Trig.sech2(0.0), 0.0);

        // tanh maps into [-1, 1] everywhere, which the old cutoff also did but
        // by giving up too early
        long lcg = 13L;
        for (int i = 0; i < 200000; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double x = 1600.0 * (((lcg >>> 11) * 0x1.0p-53) - 0.5);
            double t = Trig.tanh(x);
            assertFalse("tanh(" + x + ") is NaN", Double.isNaN(t));
            assertTrue("tanh(" + x + ") = " + t + " is out of range", t >= -1.0 && t <= 1.0);
        }
    }

    // ----- helpers -----------------------------------------------------------

    private static void assertSameBits(String what, double expected, double actual) {
        assertEquals(what + " : expected " + expected + " but was " + actual,
                Double.doubleToRawLongBits(expected), Double.doubleToRawLongBits(actual));
    }
}
