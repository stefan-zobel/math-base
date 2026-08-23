package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * The facade: that it is pure delegation, and that the one route it adds --
 * the double-exponential rule by name -- reaches what that rule reaches and
 * refuses what it cannot do.
 */
public class QuadratureTest {

    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;
    private static final double INF = Double.POSITIVE_INFINITY;
    private static final double N_INF = Double.NEGATIVE_INFINITY;
    private static final double TOL = 1.0e-10;

    // =========================================================================
    // PURE DELEGATION
    // =========================================================================

    /**
     * Every {@code integrate} is the call underneath with the one rule the
     * layer below accepts filled in, so the two agree to the bit. A tolerance
     * here would hide a facade that quietly does something else.
     */
    @Test
    public void theFacadeIsTheCallUnderneath() {
        DFunction gauss = x -> Math.exp(-x * x);
        assertEquals("whole line", InfiniteIntegrator.integrate1DInfinite(RULE, gauss, N_INF, INF, TOL),
                Quadrature.integrate(gauss, N_INF, INF, TOL), 0.0);
        assertEquals("[0, +inf)", InfiniteIntegrator.integrate1DInfinite(RULE, gauss, 0.0, INF, TOL),
                Quadrature.integrate(gauss, 0.0, INF, TOL), 0.0);
        assertEquals("(-inf, 2]", InfiniteIntegrator.integrate1DInfinite(RULE, gauss, N_INF, 2.0, TOL),
                Quadrature.integrate(gauss, N_INF, 2.0, TOL), 0.0);
        assertEquals("finite", InfiniteIntegrator.integrate1DInfinite(RULE, gauss, -1.0, 3.0, TOL),
                Quadrature.integrate(gauss, -1.0, 3.0, TOL), 0.0);

        DFunction far = x -> Math.exp(-0.5 * (x - 1.0e4) * (x - 1.0e4)) / Math.sqrt(2.0 * Math.PI);
        assertEquals("centered", InfiniteIntegrator.integrate1DInfinite(RULE, far, N_INF, INF, TOL, 1.0e4),
                Quadrature.integrate(far, N_INF, INF, TOL, 1.0e4), 0.0);

        DBiFunction plane = (x, y) -> Math.exp(-(x * x + y * y));
        assertEquals("the plane",
                InfiniteIntegrator.integrate2DInfinite(RULE, plane, N_INF, INF, N_INF, INF, 1.0e-9),
                Quadrature.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9), 0.0);
        assertEquals("the plane, centered",
                InfiniteIntegrator.integrate2DInfinite(RULE, plane, N_INF, INF, N_INF, INF, 1.0e-9, 0.0, 0.0),
                Quadrature.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9, 0.0, 0.0), 0.0);

        DTriFunction space = (x, y, z) -> Math.exp(-(x * x + y * y + z * z));
        assertEquals("space",
                InfiniteIntegrator.integrate3DInfinite(RULE, space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6),
                Quadrature.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6), 0.0);
        assertEquals("space, centered", InfiniteIntegrator.integrate3DInfinite(RULE, space, N_INF, INF, N_INF, INF,
                N_INF, INF, 1.0e-6, 0.0, 0.0, 0.0),
                Quadrature.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6, 0.0, 0.0, 0.0), 0.0);
    }

    /** The dimension comes from the integrand, and the values are right. */
    @Test
    public void oneVerbForOneTwoAndThreeDimensions() {
        assertEquals("sqrt(pi)", Math.sqrt(Math.PI), Quadrature.integrate((DFunction) x -> Math.exp(-x * x), N_INF,
                INF, TOL), 1.0e-9);
        assertEquals("pi", Math.PI, Quadrature.integrate((DBiFunction) (x, y) -> Math.exp(-(x * x + y * y)), N_INF,
                INF, N_INF, INF, 1.0e-9), 1.0e-7);
        assertEquals("pi^(3/2)", Math.pow(Math.PI, 1.5), Quadrature.integrate(
                (DTriFunction) (x, y, z) -> Math.exp(-(x * x + y * y + z * z)), N_INF, INF, N_INF, INF, N_INF, INF,
                1.0e-6), 1.0e-5);
        assertEquals("x^2 over [0, 3]", 9.0, Quadrature.integrate((DFunction) x -> x * x, 0.0, 3.0, TOL), 1.0e-9);
    }

    // =========================================================================
    // THE ROUTE THE HEURISTIC CANNOT TAKE
    // =========================================================================

    /**
     * Nothing routes to {@link DoubleExponential} by itself, because deciding
     * whether an integrand has an endpoint singularity is guesswork. Said by
     * name it reaches what that rule reaches -- the same integrals
     * {@code DoubleExponentialTest} asserts, through the facade.
     */
    @Test
    public void theSingularRouteIsThereWhenTheCallerNamesIt() {
        assertEquals("1/sqrt(x) over [0, 1]", 2.0,
                Quadrature.integrateSingular(x -> 1.0 / Math.sqrt(x), 0.0, 1.0, 1.0e-12), 1.0e-12);
        assertEquals("-log(x) over [0, 1]", 1.0, Quadrature.integrateSingular(x -> -Math.log(x), 0.0, 1.0, 1.0e-12),
                1.0e-12);
        assertEquals("e^-x/sqrt(x) over [0, +inf)", Math.sqrt(Math.PI),
                Quadrature.integrateSingular(x -> Math.exp(-x) / Math.sqrt(x), 0.0, INF, 1.0e-12), 1.0e-12);
        assertEquals("e^(-x^2) over the whole line", Math.sqrt(Math.PI),
                Quadrature.integrateSingular(x -> Math.exp(-x * x), N_INF, INF, 1.0e-12), 1.0e-12);
    }

    /**
     * The rule reports whether it reached the tolerance and this facade returns
     * a bare number, so a run that did not converge has to be refused rather
     * than handed out without its flag.
     */
    @Test
    public void aRunThatDidNotConvergeIsRefusedRatherThanReturned() {
        // an oscillation no double-exponential node set resolves, asked for a
        // tolerance it cannot reach
        DFunction wild = x -> Math.sin(1.0e6 / (x + 1.0e-12)) / Math.sqrt(x);
        try {
            double v = Quadrature.integrateSingular(wild, 0.0, 1.0, 1.0e-14);
            fail("a run that did not converge should have been refused, but it returned " + v);
        } catch (ArithmeticException expected) {
            String m = expected.getMessage();
            assertTrue("the refusal has to say what it stopped at, was: " + m, m.contains("level"));
            assertTrue("the refusal has to name the tolerance it missed, was: " + m, m.contains("1.000e-14"));
        }
    }

    /**
     * The refusal formats numbers, so it has to read the same wherever the
     * default locale points; see {@code math.OutputLocaleTest} for the printers
     * this mirrors.
     */
    @Test
    public void theRefusalDoesNotDependOnTheDefaultLocale() {
        DFunction wild = x -> Math.sin(1.0e6 / (x + 1.0e-12)) / Math.sqrt(x);
        java.util.Locale previous = java.util.Locale.getDefault();
        String german;
        String neutral;
        try {
            java.util.Locale.setDefault(java.util.Locale.GERMANY);
            german = refusalOf(wild);
            java.util.Locale.setDefault(java.util.Locale.ROOT);
            neutral = refusalOf(wild);
        } finally {
            java.util.Locale.setDefault(previous);
        }
        assertEquals("the refusal depends on the default locale", neutral, german);
        assertTrue("the refusal printed a decimal comma: " + german, !german.matches("(?s).*\\d,\\d.*"));
    }

    private static String refusalOf(DFunction f) {
        try {
            Quadrature.integrateSingular(f, 0.0, 1.0, 1.0e-14);
            fail("expected a refusal");
            return null;
        } catch (ArithmeticException expected) {
            return expected.getMessage();
        }
    }
}
