package math.solve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * The caller-supplied center of the substitution: what it buys, what it costs
 * the calls that do not use it, and that naming the wrong point does not buy a
 * wrong answer.
 * <p>
 * The blind spot these tests close is the one the class comment of
 * {@link InfiniteIntegrator} describes: mass so far from the anchor that it
 * underflows to zero at every rung of the probing ladder, so that no amount of
 * sampling can find it and the substitution returns a silent zero. Every
 * tolerance below is the measured error rounded up, not a guess.
 */
public class InfiniteIntegratorCenterTest {

    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;
    private static final double INF = Double.POSITIVE_INFINITY;
    private static final double N_INF = Double.NEGATIVE_INFINITY;
    private static final double NONE = Double.NaN;
    private static final double TOL = 1.0e-10;

    /** A normal density, whose integral over its whole support is one. */
    private static DFunction normal(final double mu, final double sigma) {
        return x -> {
            double d = (x - mu) / sigma;
            return Math.exp(-0.5 * d * d) / (sigma * Math.sqrt(2.0 * Math.PI));
        };
    }

    /** A unit-mass Gaussian blob centered at {@code (cx, cy)}. */
    private static DBiFunction blob2D(final double cx, final double cy, final double sigma) {
        return (x, y) -> {
            double dx = (x - cx) / sigma;
            double dy = (y - cy) / sigma;
            return Math.exp(-0.5 * (dx * dx + dy * dy)) / (2.0 * Math.PI * sigma * sigma);
        };
    }

    /** The same in three dimensions, on the diagonal. */
    private static DTriFunction blob3D(final double c, final double sigma) {
        return (x, y, z) -> {
            double dx = (x - c) / sigma;
            double dy = (y - c) / sigma;
            double dz = (z - c) / sigma;
            double norm = Math.pow(2.0 * Math.PI, 1.5) * sigma * sigma * sigma;
            return Math.exp(-0.5 * (dx * dx + dy * dy + dz * dz)) / norm;
        };
    }

    private static double whole1D(DFunction f, double center) {
        return InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, INF, TOL, center);
    }

    // =========================================================================
    // WHAT THE CENTER IS FOR
    // =========================================================================

    /**
     * The case the class comment ends on. A unit-width normal at {@code 1e4}
     * underflows to exactly zero at every rung of the ladder, so the guard has
     * nothing to refine around and the substitution's silent zero stands.
     * Named, the same integrand comes back to thirteen digits.
     */
    @Test
    public void massTooFarOutForAnyProbeIsReachedByNamingIt() {
        double[] far = { 1.0e4, 1.0e5, 1.0e6 };
        for (int i = 0; i < far.length; ++i) {
            DFunction f = normal(far[i], 1.0);
            double blind = InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, INF, TOL);
            assertTrue("without a center the mass at " + far[i] + " is still supposed to be out of reach, but the "
                    + "call returned " + blind + " -- if the probe can see it now, this test has become the wrong "
                    + "witness for the center", Math.abs(blind - 1.0) > 0.5);
            // measured error 7.1e-14, 5.1e-13 and 8.9e-13 for the three
            assertEquals("centered on " + far[i], 1.0, whole1D(f, far[i]), 1.0e-11);
        }
    }

    /** Semi-infinite in both directions, where the finite end pins the anchor. */
    @Test
    public void aSemiInfiniteIntervalIsCenteredToo() {
        // measured error 3.1e-14 upward and 3.1e-14 downward
        assertEquals("[0, +inf) with the mass at 1e4", 1.0,
                InfiniteIntegrator.integrate1DInfinite(RULE, normal(1.0e4, 1.0), 0.0, INF, TOL, 1.0e4), 1.0e-11);
        assertEquals("(-inf, 0] with the mass at -1e4", 1.0,
                InfiniteIntegrator.integrate1DInfinite(RULE, normal(-1.0e4, 1.0), N_INF, 0.0, TOL, -1.0e4), 1.0e-11);
    }

    /**
     * Centering a semi-infinite interval cannot be done by scaling its
     * substitution: {@code x = a + s*t/(1-t)} does carry {@code t = 1/2} onto
     * the center, but its derivative there is {@code 4s}, so a feature of unit
     * width arrives {@code 4s} times narrower and the rule misses it exactly as
     * before. Measured with that scaling in place: {@code 4.2e-98} for the case
     * above. The substitution used instead has derivative one at the center,
     * and this is the witness for that difference -- the mass keeps its width,
     * so a rule that resolves it at the origin resolves it at {@code 1e4}.
     */
    @Test
    public void theCenteredSubstitutionKeepsTheWidthOfTheMass() {
        for (double sigma : new double[] { 0.5, 1.0, 10.0, 100.0 }) {
            assertEquals("[0, +inf), sigma " + sigma, 1.0,
                    InfiniteIntegrator.integrate1DInfinite(RULE, normal(1.0e4, sigma), 0.0, INF, TOL, 1.0e4), 1.0e-9);
            assertEquals("whole line, sigma " + sigma, 1.0, whole1D(normal(1.0e4, sigma), 1.0e4), 1.0e-9);
        }
    }

    /** Two dimensions, where the class comment says the failure bites sooner. */
    @Test
    public void thePlaneIsCenteredPerAxis() {
        // uncentered this one is refused, because the probe does find the mass
        // and the substitution demonstrably did not
        try {
            InfiniteIntegrator.integrate2DInfinite(RULE, blob2D(100.0, 100.0, 1.0), N_INF, INF, N_INF, INF, 1.0e-9);
            fail("the uncentered call is supposed to refuse mass at (100, 100)");
        } catch (ArithmeticException expected) {
            // that is the contract
        }
        // measured error 8.9e-16
        assertEquals(1.0, InfiniteIntegrator.integrate2DInfinite(RULE, blob2D(100.0, 100.0, 1.0), N_INF, INF, N_INF,
                INF, 1.0e-9, 100.0, 100.0), 1.0e-11);

        // and at (1000, 1000) the probe cannot even see it, so the uncentered
        // call returns zero rather than refusing
        double blind = InfiniteIntegrator.integrate2DInfinite(RULE, blob2D(1000.0, 1000.0, 1.0), N_INF, INF, N_INF,
                INF, 1.0e-9);
        assertTrue("mass at (1000, 1000) is supposed to be out of the probe's reach, but the uncentered call "
                + "returned " + blind, Math.abs(blind - 1.0) > 0.5);
        // measured error 7.8e-15
        assertEquals(1.0, InfiniteIntegrator.integrate2DInfinite(RULE, blob2D(1000.0, 1000.0, 1.0), N_INF, INF,
                N_INF, INF, 1.0e-9, 1000.0, 1000.0), 1.0e-11);
    }

    /**
     * Three dimensions. The tolerance is loose on purpose: at {@code 1e-6} this
     * costs ten million evaluations of the integrand and a quarter of a second,
     * where {@code 1e-9} costs fifty million for the same twelve digits.
     */
    @Test
    public void spaceIsCenteredPerAxis() {
        double blind = InfiniteIntegrator.integrate3DInfinite(RULE, blob3D(500.0, 1.0), N_INF, INF, N_INF, INF, N_INF,
                INF, 1.0e-6);
        assertTrue("mass at (500, 500, 500) is supposed to be out of the probe's reach, but the uncentered call "
                + "returned " + blind, Math.abs(blind - 1.0) > 0.5);
        // measured error 1.5e-11
        assertEquals(1.0, InfiniteIntegrator.integrate3DInfinite(RULE, blob3D(500.0, 1.0), N_INF, INF, N_INF, INF,
                N_INF, INF, 1.0e-6, 500.0, 500.0, 500.0), 1.0e-8);
    }

    // =========================================================================
    // A CENTER IS A PROMISE, AND IT IS STILL CHECKED
    // =========================================================================

    /**
     * The probe samples the center in addition to everything it sampled before,
     * so a center that is nowhere near the mass leaves the guard intact. In one
     * dimension the guard does better than refuse: it splits where the mass
     * really is. Both of these came back wrong before that was added -- the
     * second returned {@code 3.1e-72}, because the substitution did sample the
     * mass, just badly, which without a center cannot happen.
     */
    @Test
    public void aCenterInTheWrongPlaceDoesNotBuyAWrongAnswer() {
        assertEquals("mass at 0, center at 1e4", 1.0, whole1D(normal(0.0, 1.0), 1.0e4), 1.0e-9);
        assertEquals("mass at 30, center at -30", 1.0, whole1D(normal(30.0, 1.0), -30.0), 1.0e-9);
        assertEquals("mass at -7.5, center at 250", 1.0, whole1D(normal(-7.5, 0.25), 250.0), 1.0e-9);
    }

    /**
     * What a center cannot do, stated rather than hidden: one that is no better
     * than the anchor gives back the anchor's blind spot. The center is a way
     * of telling the library something it cannot find out for itself, so a
     * caller who tells it nothing is where they were.
     */
    @Test
    public void aCenterNoBetterThanTheAnchorGivesBackTheAnchorsBlindSpot() {
        double atOrigin = whole1D(normal(1.0e4, 1.0), 0.0);
        assertTrue("centering on the origin is what the substitution already did, so this cannot improve on it, "
                + "but it returned " + atOrigin, Math.abs(atOrigin - 1.0) > 0.5);
    }

    // =========================================================================
    // THAT NOTHING ELSE MOVED
    // =========================================================================

    /**
     * The overload without a center and the overload given {@code NaN} are one
     * code path, so they agree to the bit and not merely to a tolerance.
     */
    @Test
    public void theUncenteredPathIsBitForBitWhatItWas() {
        double[][] cases = { { 0.0, 1.0 }, { 30.0, 1.0 }, { 100.0, 10.0 }, { -7.5, 0.25 }, { 3.0, 2.0 } };
        for (int i = 0; i < cases.length; ++i) {
            DFunction f = normal(cases[i][0], cases[i][1]);
            String what = "normal(" + cases[i][0] + ", " + cases[i][1] + ")";
            assertEquals(what + " over the whole line",
                    InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, INF, TOL),
                    InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, INF, TOL, NONE), 0.0);
            assertEquals(what + " over [0, +inf)", InfiniteIntegrator.integrate1DInfinite(RULE, f, 0.0, INF, TOL),
                    InfiniteIntegrator.integrate1DInfinite(RULE, f, 0.0, INF, TOL, NONE), 0.0);
            assertEquals(what + " over (-inf, 0]", InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, 0.0, TOL),
                    InfiniteIntegrator.integrate1DInfinite(RULE, f, N_INF, 0.0, TOL, NONE), 0.0);
        }

        DBiFunction g = blob2D(0.0, 0.0, 1.0);
        assertEquals("the plane", InfiniteIntegrator.integrate2DInfinite(RULE, g, N_INF, INF, N_INF, INF, 1.0e-9),
                InfiniteIntegrator.integrate2DInfinite(RULE, g, N_INF, INF, N_INF, INF, 1.0e-9, NONE, NONE), 0.0);

        DTriFunction h = blob3D(0.0, 1.0);
        assertEquals("space",
                InfiniteIntegrator.integrate3DInfinite(RULE, h, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6),
                InfiniteIntegrator.integrate3DInfinite(RULE, h, N_INF, INF, N_INF, INF, N_INF, INF, 1.0e-6, NONE,
                        NONE, NONE),
                0.0);
    }

    /**
     * An interval with two finite ends is not substituted at all, so there is
     * nothing to center and a center given for one is dropped rather than
     * silently changing the rule.
     */
    @Test
    public void aCenterOnAFiniteIntervalIsIgnored() {
        DFunction square = x -> x * x;
        double plain = InfiniteIntegrator.integrate1DInfinite(RULE, square, 0.0, 3.0, TOL);
        assertEquals(9.0, plain, 1.0e-9);
        for (double c : new double[] { 0.0, 1.5, 2.0, 3.0, -100.0, 1.0e6 }) {
            assertEquals("center " + c + " has nothing to act on", plain,
                    InfiniteIntegrator.integrate1DInfinite(RULE, square, 0.0, 3.0, TOL, c), 0.0);
        }
    }

    /** One axis centered, the other left alone. */
    @Test
    public void oneAxisMayBeCenteredAndTheOtherNot() {
        // mass far out in x, at the origin in y
        DBiFunction f = blob2D(1000.0, 0.0, 1.0);
        assertEquals(1.0, InfiniteIntegrator.integrate2DInfinite(RULE, f, N_INF, INF, N_INF, INF, 1.0e-9, 1000.0,
                NONE), 1.0e-9);
    }

    // =========================================================================
    // VALIDATION
    // =========================================================================

    @Test
    public void anImpossibleCenterIsRejected() {
        DFunction f = normal(0.0, 1.0);
        // infinite, in either direction, on any interval
        expectRejected(f, N_INF, INF, INF);
        expectRejected(f, N_INF, INF, N_INF);
        expectRejected(f, 0.0, INF, INF);
        // not strictly inside an interval that has a finite end
        expectRejected(f, 0.0, INF, 0.0);
        expectRejected(f, 0.0, INF, -1.0);
        expectRejected(f, N_INF, 5.0, 5.0);
        expectRejected(f, N_INF, 5.0, 6.0);
    }

    private static void expectRejected(DFunction f, double a, double b, double center) {
        try {
            InfiniteIntegrator.integrate1DInfinite(RULE, f, a, b, TOL, center);
            fail("center " + center + " on [" + a + ", " + b + "] should have been rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message has to name the parameter, was: " + expected.getMessage(),
                    expected.getMessage().contains("center"));
        }
    }

    @Test
    public void aRejectedCenterNamesItsAxisInTwoAndThreeDimensions() {
        DBiFunction f = blob2D(0.0, 0.0, 1.0);
        try {
            InfiniteIntegrator.integrate2DInfinite(RULE, f, 0.0, INF, N_INF, INF, 1.0e-9, 0.0, 0.0);
            fail("centerX = 0 on [0, +inf) should have been rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("was: " + expected.getMessage(), expected.getMessage().contains("centerX"));
        }
        DTriFunction g = blob3D(0.0, 1.0);
        try {
            InfiniteIntegrator.integrate3DInfinite(RULE, g, N_INF, INF, N_INF, INF, N_INF, 2.0, 1.0e-6, 0.0, 0.0,
                    2.0);
            fail("centerZ = 2 on (-inf, 2] should have been rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("was: " + expected.getMessage(), expected.getMessage().contains("centerZ"));
        }
    }
}
