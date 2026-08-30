package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
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


    // =========================================================================
    // THE ERROR-RETURNING FORM
    // =========================================================================

    /**
     * {@code WithError} is the same integral, not a second one. The bare form
     * delegates to it, so the two agree to the bit; a tolerance here would hide
     * a second route quietly doing something else.
     */
    @Test
    public void theErrorFormIsTheSameIntegral() {
        DFunction gauss = x -> Math.exp(-x * x);
        assertEquals("whole line", Quadrature.integrate(gauss, N_INF, INF, TOL),
                Quadrature.WithError.integrate(gauss, N_INF, INF, TOL).value, 0.0);
        assertEquals("[0, +inf)", Quadrature.integrate(gauss, 0.0, INF, TOL),
                Quadrature.WithError.integrate(gauss, 0.0, INF, TOL).value, 0.0);
        assertEquals("(-inf, 2]", Quadrature.integrate(gauss, N_INF, 2.0, TOL),
                Quadrature.WithError.integrate(gauss, N_INF, 2.0, TOL).value, 0.0);
        assertEquals("finite", Quadrature.integrate(gauss, -1.0, 3.0, TOL),
                Quadrature.WithError.integrate(gauss, -1.0, 3.0, TOL).value, 0.0);

        DFunction far = x -> Math.exp(-0.5 * (x - 1.0e4) * (x - 1.0e4)) / Math.sqrt(2.0 * Math.PI);
        assertEquals("centered", Quadrature.integrate(far, N_INF, INF, TOL, 1.0e4),
                Quadrature.WithError.integrate(far, N_INF, INF, TOL, 1.0e4).value, 0.0);

        DBiFunction plane = (x, y) -> Math.exp(-(x * x + y * y));
        assertEquals("the plane", Quadrature.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9),
                Quadrature.WithError.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9).value, 0.0);
        assertEquals("the plane, centered", Quadrature.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9, 0.0, 0.0),
                Quadrature.WithError.integrate(plane, N_INF, INF, N_INF, INF, 1.0e-9, 0.0, 0.0).value, 0.0);

        DTriFunction space = (x, y, z) -> Math.exp(-(x * x + y * y + z * z));
        assertEquals("space", Quadrature.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6),
                Quadrature.WithError.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6).value, 0.0);
        assertEquals("space, centered",
                Quadrature.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6, 0.0, 0.0, 0.0),
                Quadrature.WithError.integrate(space, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6, 0.0, 0.0,
                        0.0).value, 0.0);

        assertEquals("1/sqrt(x), by name", Quadrature.integrateSingular(x -> 1.0 / Math.sqrt(x), 0.0, 1.0, 1.0e-12),
                Quadrature.WithError.integrateSingular(x -> 1.0 / Math.sqrt(x), 0.0, 1.0, 1.0e-12).value, 0.0);
        assertEquals("-log(x), by name", Quadrature.integrateSingular(x -> -Math.log(x), 0.0, 1.0, 1.0e-12),
                Quadrature.WithError.integrateSingular(x -> -Math.log(x), 0.0, 1.0, 1.0e-12).value, 0.0);
    }

    /**
     * Where the rule reaches the tolerance it says so, and the summed panel
     * estimate is inside the tolerance the caller asked for -- which it has to
     * be, since the tolerance is halved at every level and the panel
     * tolerances of a full tree sum back to it.
     */
    @Test
    public void aConvergedRunReportsAnEstimateInsideItsTolerance() {
        DFunction gauss = x -> Math.exp(-x * x);
        assertConverged("gauss, whole line", Quadrature.WithError.integrate(gauss, N_INF, INF, TOL), TOL);
        assertConverged("gauss, [0, +inf)", Quadrature.WithError.integrate(gauss, 0.0, INF, TOL), TOL);
        assertConverged("x^2 over [0, 3]",
                Quadrature.WithError.integrate((DFunction) x -> x * x, 0.0, 3.0, TOL), TOL);
        assertConverged("the plane",
                Quadrature.WithError.integrate((DBiFunction) (x, y) -> Math.exp(-(x * x + y * y)),
                        N_INF, INF, N_INF, INF, 1.0e-9), 1.0e-9);
        assertConverged("space",
                Quadrature.WithError.integrate((DTriFunction) (x, y, z) -> Math.exp(-(x * x + y * y + z * z)),
                        N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6), 1.0e-6);
    }

    private static void assertConverged(String what, Quadrature.Result r, double epsTol) {
        assertTrue(what + " should have met its tolerance, was " + r, r.toleranceMet);
        assertTrue(what + " estimate " + r.approximatedErrorEstimate + " exceeds the tolerance " + epsTol,
                r.approximatedErrorEstimate <= epsTol);
    }

    /**
     * An endpoint singularity on the ordinary route is exactly what the flag is
     * for: the value really is wrong, and nothing in a bare {@code double}
     * could have said so. {@code integrateSingular} is the answer, and the
     * point here is that the caller is told to reach for it.
     */
    @Test
    public void anEndpointSingularityOnTheOrdinaryRouteReportsItself() {
        Quadrature.Result root = Quadrature.WithError.integrate((DFunction) x -> 1.0 / Math.sqrt(x), 0.0, 1.0, TOL);
        assertFalse("1/sqrt(x) on the ordinary route should not claim convergence, was " + root, root.toleranceMet);
        assertTrue("and it should be wrong by more than the tolerance, was " + root.value,
                Math.abs(root.value - 2.0) > TOL);

        Quadrature.Result log = Quadrature.WithError.integrate((DFunction) x -> Math.log(x), 0.0, 1.0, TOL);
        assertFalse("log(x) on the ordinary route should not claim convergence, was " + log, log.toleranceMet);

        // the named route gets both of them right, and says so
        assertConverged("1/sqrt(x) by name",
                Quadrature.WithError.integrateSingular(x -> 1.0 / Math.sqrt(x), 0.0, 1.0, 1.0e-12), 1.0e-12);
    }

    /**
     * The counterpart of
     * {@link #aRunThatDidNotConvergeIsRefusedRatherThanReturned()}: the same
     * integrand that the bare form has to refuse comes back from
     * {@code WithError} with its flag, and the caller decides.
     */
    @Test
    public void aRunThatDidNotConvergeComesBackWithItsFlag() {
        DFunction wild = x -> Math.sin(1.0e6 / (x + 1.0e-12)) / Math.sqrt(x);
        Quadrature.Result r = Quadrature.WithError.integrateSingular(wild, 0.0, 1.0, 1.0e-14);
        assertFalse("the run did not converge and should say so, was " + r, r.toleranceMet);
        assertTrue("a value is still handed out, was " + r.value,
                !Double.isNaN(r.value) && !Double.isInfinite(r.value));
        assertTrue("and the estimate says how badly, was " + r.approximatedErrorEstimate,
                r.approximatedErrorEstimate > 1.0e-14);
    }

    /**
     * The caveat the class comment of {@code WithError} carries, pinned so it
     * cannot quietly stop being true: {@code toleranceMet} rules out an exhausted
     * budget, not a missed feature. A Gaussian of width 0.005 in the unit cube
     * is narrower than the forced panels of the 3D subdivision, so the Kronrod
     * and the Gauss rule miss it alike, agree, and report a tiny error on an
     * answer that is wrong by 99.6 percent.
     * <p>
     * The same peak in one dimension is resolved to machine precision, which is
     * what makes this a property of the 3D panel spacing and not of the
     * integrand.
     */
    @Test
    public void toleranceMetRulesOutAnExhaustedBudgetAndNotAMissedPeak() {
        final double w = 0.005;
        DTriFunction peak = (x, y, z) -> Math.exp(-(((x - 0.37) * (x - 0.37) + (y - 0.61) * (y - 0.61)
                + (z - 0.44) * (z - 0.44)) / (w * w)));
        double exact = Math.pow(Math.PI, 1.5) * w * w * w;
        Quadrature.Result r = Quadrature.WithError.integrate(peak, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0e-6);
        assertTrue("the 3D peak is missed and reported as within tolerance, which is the point, was " + r,
                r.toleranceMet);
        assertTrue("and the answer is wrong by most of itself, was " + r.value,
                Math.abs(r.value - exact) / exact > 0.9);
        assertTrue("while the estimate stays inside the tolerance, was " + r.approximatedErrorEstimate,
                r.approximatedErrorEstimate <= 1.0e-6);

        DFunction same = x -> Math.exp(-((x - 0.37) / w) * ((x - 0.37) / w));
        Quadrature.Result oneD = Quadrature.WithError.integrate(same, 0.0, 1.0, TOL);
        assertEquals("the same peak in one dimension is resolved", w * Math.sqrt(Math.PI), oneD.value, 1.0e-15);
        assertTrue("and reports it, was " + oneD, oneD.toleranceMet);
    }

    /** The result prints the same wherever the default locale points. */
    @Test
    public void theResultDoesNotDependOnTheDefaultLocale() {
        java.util.Locale previous = java.util.Locale.getDefault();
        String german;
        String neutral;
        try {
            java.util.Locale.setDefault(java.util.Locale.GERMANY);
            german = Quadrature.WithError.integrate((DFunction) x -> x * x, 0.0, 3.0, TOL).toString();
            java.util.Locale.setDefault(java.util.Locale.ROOT);
            neutral = Quadrature.WithError.integrate((DFunction) x -> x * x, 0.0, 3.0, TOL).toString();
        } finally {
            java.util.Locale.setDefault(previous);
        }
        assertEquals("the result depends on the default locale", neutral, german);
        assertTrue("the result printed a decimal comma: " + german, !german.matches("(?s).*\\d,\\d.*"));
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
