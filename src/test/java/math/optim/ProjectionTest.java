package math.optim;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

public class ProjectionTest {

    private static final double[] LOWER = { -1.0, 0.0, -3.5, 2.0 };
    private static final double[] UPPER = { 1.0, 4.0, -0.5, 2.0 };

    /** Deterministic LCG, so a failure is always reproducible. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            state = seed;
        }

        /** Uniform in {@code [0, 1)}. */
        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }

        /** Uniform in {@code [lo, hi)}. */
        double next(double lo, double hi) {
            return lo + next() * (hi - lo);
        }
    }

    /** Draws a point that is known to be inside the set. */
    private interface Sampler {
        double[] sample(Lcg rng);
    }

    // -----------------------------------------------------------------
    // box
    // -----------------------------------------------------------------

    @Test
    public void testBoxClampsEachCoordinateIntoItsInterval() {
        Projection p = Projection.box(LOWER, UPPER);
        double[] x = { -7.0, 2.0, 9.0, -100.0 };
        p.projectInto(x);
        assertArrayEquals(new double[] { -1.0, 2.0, -0.5, 2.0 }, x, 0.0);
    }

    @Test
    public void testBoxHandlesInfiniteBounds() {
        Projection p = Projection.box(new double[] { Double.NEGATIVE_INFINITY, 0.0 },
                new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY });
        double[] x = { -1.0e300, -5.0 };
        p.projectInto(x);
        assertArrayEquals(new double[] { -1.0e300, 0.0 }, x, 0.0);
    }

    @Test
    public void testBoxWithEqualBoundsPinsTheCoordinateExactly() {
        // LOWER[3] == UPPER[3] == 2.0, which is how a parameter is held fixed
        Projection p = Projection.box(LOWER, UPPER);
        Lcg rng = new Lcg(20260820L);
        for (int trial = 0; trial < 100; trial++) {
            double[] x = { rng.next(-9.0, 9.0), rng.next(-9.0, 9.0), rng.next(-9.0, 9.0), rng.next(-9.0, 9.0) };
            p.projectInto(x);
            assertEquals(2.0, x[3], 0.0);
        }
    }

    @Test
    public void testBoxDoesNotDependOnTheBoundArraysAfterwards() {
        double[] lower = { -1.0, -1.0 };
        double[] upper = { 1.0, 1.0 };
        Projection p = Projection.box(lower, upper);
        lower[0] = -100.0;
        upper[1] = 100.0;
        double[] x = { -50.0, 50.0 };
        p.projectInto(x);
        assertArrayEquals(new double[] { -1.0, 1.0 }, x, 0.0);
    }

    // -----------------------------------------------------------------
    // non-negative orthant
    // -----------------------------------------------------------------

    @Test
    public void testNonNegativeIsTheMaximumWithZero() {
        Projection p = Projection.nonNegative(5);
        double[] x = { -1.0, 0.0, 3.0, -1.0e-300, 1.0e300 };
        p.projectInto(x);
        assertArrayEquals(new double[] { 0.0, 0.0, 3.0, 0.0, 1.0e300 }, x, 0.0);
    }

    // -----------------------------------------------------------------
    // simplex
    // -----------------------------------------------------------------

    @Test
    public void testSimplexImageIsFeasible() {
        Lcg rng = new Lcg(7L);
        for (int n = 1; n <= 12; n++) {
            for (double sum : new double[] { 1.0, 0.25, 17.0 }) {
                Projection p = Projection.simplex(n, sum);
                for (int trial = 0; trial < 50; trial++) {
                    double[] x = new double[n];
                    for (int i = 0; i < n; i++) {
                        x[i] = rng.next(-6.0, 6.0);
                    }
                    p.projectInto(x);
                    double total = 0.0;
                    for (int i = 0; i < n; i++) {
                        assertTrue("negative coordinate " + x[i], x[i] >= 0.0);
                        total += x[i];
                    }
                    assertEquals(sum, total, 1.0e-14 * sum);
                }
            }
        }
    }

    @Test
    public void testSimplexOfAnAllEqualVectorIsTheBarycentre() {
        for (int n = 1; n <= 8; n++) {
            double[] x = new double[n];
            java.util.Arrays.fill(x, -4.25);
            Projection.simplex(n, 3.0).projectInto(x);
            for (int i = 0; i < n; i++) {
                assertEquals(3.0 / n, x[i], 1.0e-15);
            }
        }
    }

    @Test
    public void testSimplexOfAFeasiblePointIsThatPoint() {
        Lcg rng = new Lcg(99L);
        int n = 6;
        Projection p = Projection.simplex(n, 1.0);
        for (int trial = 0; trial < 100; trial++) {
            double[] x = new double[n];
            double total = 0.0;
            for (int i = 0; i < n; i++) {
                x[i] = rng.next();
                total += x[i];
            }
            for (int i = 0; i < n; i++) {
                x[i] /= total;
            }
            double[] projected = x.clone();
            p.projectInto(projected);
            assertArrayEquals(x, projected, 1.0e-15);
        }
    }

    /**
     * The classical form subtracts an absolute threshold and cancels here:
     * {@code 4e200 - 1} rounds back to {@code 4e200}, so the point would come
     * back unshrunk and far outside the simplex.
     */
    @Test
    public void testSimplexResolvesATargetFarBelowTheCoordinateScale() {
        double[] x = { 3.0e200, 4.0e200 };
        Projection.simplex(2, 1.0).projectInto(x);
        assertArrayEquals(new double[] { 0.0, 1.0 }, x, 0.0);

        double[] y = { -1.0e300, -1.0e300 };
        Projection.simplex(2, 1.0).projectInto(y);
        assertArrayEquals(new double[] { 0.5, 0.5 }, y, 0.0);
    }

    // -----------------------------------------------------------------
    // balls
    // -----------------------------------------------------------------

    @Test
    public void testEuclideanBallLeavesTheInteriorAloneAndScalesTheOutside() {
        Projection p = Projection.euclideanBall(2, 5.0);
        double[] inside = { 3.0, 4.0 };
        p.projectInto(inside);
        assertArrayEquals(new double[] { 3.0, 4.0 }, inside, 0.0);

        double[] outside = { 30.0, 40.0 };
        p.projectInto(outside);
        assertArrayEquals(new double[] { 3.0, 4.0 }, outside, 1.0e-15);
    }

    /** A direct sum of squares would overflow to infinity on this input. */
    @Test
    public void testEuclideanBallDoesNotOverflow() {
        double[] x = { 3.0e200, 4.0e200 };
        Projection.euclideanBall(2, 1.0).projectInto(x);
        assertArrayEquals(new double[] { 0.6, 0.8 }, x, 1.0e-15);
    }

    @Test
    public void testL1BallShrinksZeroesAndKeepsTheSigns() {
        double[] x = { -5.0, 0.0, 0.25 };
        Projection.l1Ball(3, 1.0).projectInto(x);
        assertArrayEquals(new double[] { -1.0, 0.0, 0.0 }, x, 1.0e-15);
        assertEquals(1.0, Math.abs(x[0]) + Math.abs(x[1]) + Math.abs(x[2]), 1.0e-15);
    }

    @Test
    public void testL1BallLeavesTheInteriorAlone() {
        double[] x = { -0.5, 0.25, 0.125 };
        double[] expected = x.clone();
        Projection.l1Ball(3, 1.0).projectInto(x);
        assertArrayEquals(expected, x, 0.0);
    }

    @Test
    public void testL1BallResolvesARadiusFarBelowTheCoordinateScale() {
        double[] x = { -3.0e200, 4.0e200 };
        Projection.l1Ball(2, 1.0).projectInto(x);
        assertArrayEquals(new double[] { 0.0, 1.0 }, x, 0.0);
    }

    // -----------------------------------------------------------------
    // the defining properties, over all five sets
    // -----------------------------------------------------------------

    /**
     * The projection is the <i>nearest</i> point of the set. Checked against
     * brute force rather than against a formula: no feasible point drawn at
     * random may be closer to {@code x} than the projection is. This tests the
     * definition, which a spot check of a few values does not.
     */
    @Test
    public void testEveryProjectionReturnsTheNearestFeasiblePoint() {
        Lcg rng = new Lcg(4242L);
        int n = 4;
        assertNearest("box", Projection.box(LOWER, UPPER), n, r -> {
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = r.next(LOWER[i], UPPER[i]);
            }
            return y;
        }, rng);
        assertNearest("nonNegative", Projection.nonNegative(n), n, r -> {
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = r.next(0.0, 8.0);
            }
            return y;
        }, rng);
        assertNearest("simplex", Projection.simplex(n, 1.0), n, r -> scaleTo(draw(r, n, 0.0, 1.0), 1.0, 1), rng);
        assertNearest("euclideanBall", Projection.euclideanBall(n, 2.0), n,
                r -> scaleTo(draw(r, n, -1.0, 1.0), 2.0 * r.next(), 2), rng);
        assertNearest("l1Ball", Projection.l1Ball(n, 2.0), n,
                r -> scaleTo(draw(r, n, -1.0, 1.0), 2.0 * r.next(), 1), rng);
    }

    /**
     * Projecting twice is projecting once. Exact for the box and the orthant,
     * which only ever copy a bound in; up to the rounding of one shift or one
     * scaling for the other three.
     */
    @Test
    public void testEveryProjectionIsIdempotent() {
        Lcg rng = new Lcg(31337L);
        int n = 6;
        assertIdempotent(Projection.box(LOWER, UPPER), 4, 0.0, rng);
        assertIdempotent(Projection.nonNegative(n), n, 0.0, rng);
        assertIdempotent(Projection.simplex(n, 1.0), n, 1.0e-15, rng);
        assertIdempotent(Projection.euclideanBall(n, 2.0), n, 1.0e-14, rng);
        assertIdempotent(Projection.l1Ball(n, 2.0), n, 1.0e-14, rng);
    }

    // -----------------------------------------------------------------
    // validation
    // -----------------------------------------------------------------

    @Test
    public void testBoxValidation() {
        double[] two = { 0.0, 0.0 };
        expectIae("null lower", () -> Projection.box(null, two));
        expectIae("null upper", () -> Projection.box(two, null));
        expectIae("length mismatch", () -> Projection.box(two, new double[] { 0.0, 0.0, 0.0 }));
        expectIae("empty", () -> Projection.box(new double[0], new double[0]));
        expectIae("NaN lower", () -> Projection.box(new double[] { Double.NaN, 0.0 }, two));
        expectIae("NaN upper", () -> Projection.box(two, new double[] { Double.NaN, 0.0 }));
        expectIae("lower above upper", () -> Projection.box(new double[] { 1.0, 0.0 }, two));
    }

    @Test
    public void testDimensionAndParameterValidation() {
        expectIae("dimension 0", () -> Projection.nonNegative(0));
        expectIae("negative dimension", () -> Projection.simplex(-1, 1.0));
        expectIae("sum 0", () -> Projection.simplex(3, 0.0));
        expectIae("negative sum", () -> Projection.simplex(3, -1.0));
        expectIae("infinite sum", () -> Projection.simplex(3, Double.POSITIVE_INFINITY));
        expectIae("NaN sum", () -> Projection.simplex(3, Double.NaN));
        expectIae("radius 0", () -> Projection.euclideanBall(3, 0.0));
        expectIae("negative radius", () -> Projection.l1Ball(3, -2.0));
    }

    @Test
    public void testWrongLengthIsRejectedAtProjectionTime() {
        double[] wrong = new double[7];
        expectIae("box", () -> Projection.box(LOWER, UPPER).projectInto(wrong));
        expectIae("nonNegative", () -> Projection.nonNegative(3).projectInto(wrong));
        expectIae("simplex", () -> Projection.simplex(3, 1.0).projectInto(wrong));
        expectIae("euclideanBall", () -> Projection.euclideanBall(3, 1.0).projectInto(wrong));
        expectIae("l1Ball", () -> Projection.l1Ball(3, 1.0).projectInto(wrong));
        expectIae("null point", () -> Projection.nonNegative(3).projectInto(null));
    }

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    private static void assertNearest(String name, Projection p, int n, Sampler feasible, Lcg rng) {
        for (int trial = 0; trial < 40; trial++) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = rng.next(-6.0, 6.0);
            }
            double[] projected = x.clone();
            p.projectInto(projected);
            double best = distance(projected, x);
            for (int k = 0; k < 200; k++) {
                double[] y = feasible.sample(rng);
                double d = distance(y, x);
                assertTrue(name + ": a feasible point at distance " + d + " beats the projection at " + best,
                        d >= best - 1.0e-12 * (1.0 + best));
            }
        }
    }

    private static void assertIdempotent(Projection p, int n, double tolerance, Lcg rng) {
        for (int trial = 0; trial < 100; trial++) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = rng.next(-6.0, 6.0);
            }
            p.projectInto(x);
            double[] once = x.clone();
            p.projectInto(x);
            assertArrayEquals(once, x, tolerance);
        }
    }

    private static double[] draw(Lcg rng, int n, double lo, double hi) {
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = rng.next(lo, hi);
        }
        return y;
    }

    /** Rescales {@code y} so that its {@code p}-norm is exactly {@code norm}. */
    private static double[] scaleTo(double[] y, double norm, int p) {
        double current = 0.0;
        for (int i = 0; i < y.length; i++) {
            current += (p == 1) ? Math.abs(y[i]) : y[i] * y[i];
        }
        if (p == 2) {
            current = Math.sqrt(current);
        }
        if (current == 0.0) {
            return y;
        }
        double factor = norm / current;
        for (int i = 0; i < y.length; i++) {
            y[i] *= factor;
        }
        return y;
    }

    private static double distance(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) {
            double d = a[i] - b[i];
            s += d * d;
        }
        return Math.sqrt(s);
    }

    private static void expectIae(String what, Runnable r) {
        try {
            r.run();
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
