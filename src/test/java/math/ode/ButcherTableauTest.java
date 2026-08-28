package math.ode;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
 * <b>Those seventeen are written out by hand, and beyond them they are
 * generated.</b> Order eight needs 200 conditions, which nobody should type,
 * so the trees are enumerated recursively and the elementary weight and the
 * density of each are computed from its shape. The hand written table is what
 * makes the generator believable: it must produce exactly 1, 1, 2, 4, 9 trees
 * at orders one through five and the same seventeen residuals, which it does to
 * {@code 4.2e-16} over both tableaux and both weight vectors. Deleting the hand
 * table would leave the generator checking itself.
 * <p>
 * The measurements the tolerances below rest on: over all three tableaux the
 * row sums deviate from the stage times by at most {@code 2.2e-16}, every
 * condition a method is supposed to satisfy comes out below {@code 5.0e-15},
 * and every condition it is supposed to <em>fail</em> is missed by at least
 * {@code 1.4e-05}. Ten orders of magnitude separate the two, so nothing here is
 * close to its threshold.
 * <p>
 * <b>Failing the next order is not the same as failing all of it.</b>
 * Dormand-Prince meets nine of the twenty conditions of order six by accident,
 * which is why "of order five" is asserted as all seventeen met and <em>some</em>
 * of the twenty missed. For the two DOP853 error rows it happens to be all of
 * them, and that is recorded where it is measured rather than assumed here.
 */
public final class ButcherTableauTest {

    /** Measured worst case 5.0e-15 over all three tableaux and all conditions. */
    private static final double SATISFIED = 1.0e-14;
    /** Measured smallest violation 1.4e-05, on the DOP853 fifth order row at order six. */
    private static final double VIOLATED = 1.0e-06;
    /**
     * The interpolant reaches the ends of the step to {@code 6.0e-14} for
     * DOP853 and its slope there to {@code 4.6e-13}, both over sixteen stages
     * of coefficients that run to {@code 5.3e+02}. Dormand-Prince, over seven
     * stages, does the same at {@link #SATISFIED}.
     */
    private static final double INTERPOLANT_ENDS = 5.0e-12;

    @Test
    public void testTheRowsOfTheStageMatrixSumToTheStageTimes() {
        assertEquals(0.0, worstRowSum(ButcherTableau.CLASSIC_RK4), 0.0);
        assertTrue(worstRowSum(ButcherTableau.DORMAND_PRINCE_45) <= SATISFIED);
    }

    @Test
    public void testTheTreeCountsAreTheKnownOnes() {
        List<Tree> trees = trees(8);
        int[] count = new int[9];
        for (int i = 0; i < trees.size(); ++i) {
            ++count[trees.get(i).order];
        }
        assertArrayEquals(new int[] { 0, 1, 1, 2, 4, 9, 20, 48, 115 }, count);
        assertEquals(200, trees.size());
        // and they come out grouped by order, which is what lets a slice of the
        // residual array be read as "the conditions of order p and below"
        for (int i = 1; i < trees.size(); ++i) {
            assertTrue(trees.get(i - 1).order <= trees.get(i).order);
        }
        assertArrayEquals(new int[] { 1, 2, 4, 8, 17, 37, 85, 200 }, THROUGH_ORDER);
    }

    @Test
    public void testTheGeneratorReproducesTheSeventeenWrittenOutByHand() {
        agreeWithHand(ButcherTableau.CLASSIC_RK4, ButcherTableau.CLASSIC_RK4.b());
        agreeWithHand(ButcherTableau.DORMAND_PRINCE_45, ButcherTableau.DORMAND_PRINCE_45.b());
        agreeWithHand(ButcherTableau.DORMAND_PRINCE_45, ButcherTableau.DORMAND_PRINCE_45.bStar());
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

    /**
     * The one that matters, and the reason the generator exists: 151
     * coefficients were copied out of a Fortran file, and all 200 conditions
     * coming out at zero is a statement about every one of them at once. A
     * single mistyped digit anywhere in {@code a}, {@code b} or {@code c} moves
     * at least one of the 200.
     */
    @Test
    public void testDop853SatisfiesAllTwoHundredConditionsThroughOrderEight() {
        double[] r = generated(ButcherTableau.DOP853, ButcherTableau.DOP853.b(), 8);
        assertEquals(200, r.length);
        assertSatisfied(r, 0, r.length);
    }

    @Test
    public void testTheTwoDop853ErrorRowsAreOfOrderFiveAndOfOrderThree() {
        ButcherTableau t = ButcherTableau.DOP853;
        double[] five = generated(t, t.bStar(), 6);
        assertSatisfied(five, 0, THROUGH_ORDER[4]);
        // all twenty of order six are missed here, by between 1.4e-05 and
        // 4.5e-04, where Dormand-Prince meets nine of its twenty by accident
        assertViolated(five, THROUGH_ORDER[4], THROUGH_ORDER[5]);

        double[] three = generated(t, t.bStarSecondary(), 4);
        assertSatisfied(three, 0, THROUGH_ORDER[2]);
        // and all four of order four, by between 4.2e-03 and 2.5e-02
        assertViolated(three, THROUGH_ORDER[2], THROUGH_ORDER[3]);
    }

    /**
     * That the two thresholds are the property they claim to be. Each is
     * compared against the point where the method's own stability function
     * leaves the unit disc going along the negative real axis, computed here
     * and copied from nowhere: {@code 3.31} for Dormand-Prince against a
     * threshold of {@code 3.25}, and {@code 6.39} for DOP853 against
     * {@code 6.10}. Both sit just inside, which is what a threshold for "this
     * step is being held down by stability" has to do.
     */
    @Test
    public void testTheStiffnessThresholdsSitJustInsideTheStabilityBoundaries() {
        checkThreshold(ButcherTableau.DORMAND_PRINCE_45, 3.3065);
        checkThreshold(ButcherTableau.DOP853, 6.3936);
    }

    private static void checkThreshold(ButcherTableau t, double expectedBoundary) {
        double boundary = realStabilityBoundary(t);
        assertEquals(t.name() + " stability boundary", expectedBoundary, boundary, 1.0e-3);
        assertTrue(t.name() + " threshold " + t.stiffnessThreshold() + " is not inside " + boundary,
                t.stiffnessThreshold() < boundary);
        // and not so far inside that it would fire on a step that is fine
        assertTrue(t.name() + " threshold is needlessly cautious",
                t.stiffnessThreshold() > 0.9 * boundary);
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
        reachesTheEnd(ButcherTableau.DORMAND_PRINCE_45, SATISFIED);
        reachesTheEnd(ButcherTableau.DOP853, INTERPOLANT_ENDS);
    }

    private static void reachesTheEnd(ButcherTableau t, double tolerance) {
        double[][] d = t.dense();
        double[] b = t.b();
        for (int i = 0; i < d.length; ++i) {
            double atEnd = 0.0;
            for (int j = 0; j < d[i].length; ++j) {
                atEnd += d[i][j];
            }
            // a stage the interpolant needs but the solution does not carries
            // no weight at the end of the step, which is what makes the two
            // stage counts consistent
            double weight = (i < b.length) ? b[i] : 0.0;
            assertEquals(t.name() + " stage " + i + " at the end of the step", weight, atEnd, tolerance);
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
        leavesAndArrives(ButcherTableau.DORMAND_PRINCE_45, SATISFIED);
        leavesAndArrives(ButcherTableau.DOP853, INTERPOLANT_ENDS);
    }

    private static void leavesAndArrives(ButcherTableau t, double tolerance) {
        double[][] d = t.dense();
        // the derivative at the end of the step is the last stage that
        // propagates, not the last stage there is
        int last = t.stages() - 1;
        for (int i = 0; i < d.length; ++i) {
            double slopeAtStart = d[i][0];
            double slopeAtEnd = 0.0;
            for (int j = 0; j < d[i].length; ++j) {
                slopeAtEnd += (j + 1) * d[i][j];
            }
            // exact, because the linear coefficient is written as one or zero
            assertEquals(t.name() + " stage " + i + " at the start", (i == 0) ? 1.0 : 0.0, slopeAtStart,
                    0.0);
            assertEquals(t.name() + " stage " + i + " at the end", (i == last) ? 1.0 : 0.0, slopeAtEnd,
                    tolerance);
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
        assertFalse(dp.hasSecondaryEmbedded());
        assertEquals(0, dp.secondaryEmbeddedOrder());
        assertNull(dp.bStarSecondary());
        assertEquals(3.25, dp.stiffnessThreshold(), 0.0);

        ButcherTableau d8 = ButcherTableau.DOP853;
        assertEquals(8, d8.order());
        // twelve stages plus the one at the new point, which the next step
        // starts from, and three more that only the interpolant needs
        assertEquals(13, d8.stages());
        assertEquals(16, d8.denseStages());
        assertTrue(d8.isFsal());
        assertTrue(d8.hasEmbedded());
        assertEquals(5, d8.embeddedOrder());
        assertTrue(d8.hasSecondaryEmbedded());
        assertEquals(3, d8.secondaryEmbeddedOrder());
        assertEquals(13, d8.bStarSecondary().length);
        assertTrue(d8.hasDenseOutput());
        assertEquals(7, d8.denseDegree());
        assertEquals(16, d8.dense().length);
        assertEquals(6.1, d8.stiffnessThreshold(), 0.0);
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

    /**
     * How far the stability region of a method reaches along the negative real
     * axis. On <code>y' = lambda y</code> one step multiplies the state by
     * <code>R(z) = 1 + z b.v</code> with <code>v[i] = 1 + z sum(a[i][j] v[j])</code>,
     * so the whole thing is a forward substitution through the stages, and the
     * boundary is where that leaves the unit disc.
     */
    private static double realStabilityBoundary(ButcherTableau t) {
        double step = 1.0e-4;
        double z = step;
        while (z < 100.0 && Math.abs(stabilityFunction(t, -z)) <= 1.0) {
            z += step;
        }
        return z - step;
    }

    private static double stabilityFunction(ButcherTableau t, double z) {
        double[][] a = t.a();
        double[] b = t.b();
        double[] v = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            double sum = 1.0;
            for (int j = 0; j < i; ++j) {
                sum += z * a[i][j] * v[j];
            }
            v[i] = sum;
        }
        double dot = 0.0;
        for (int i = 0; i < b.length; ++i) {
            dot += b[i] * v[i];
        }
        return 1.0 + z * dot;
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

    // ---- the same conditions, generated rather than written out ----

    /**
     * How many conditions there are through order {@code p}, at index
     * {@code p - 1}: the running sum of 1, 1, 2, 4, 9, 20, 48, 115.
     */
    private static final int[] THROUGH_ORDER = { 1, 2, 4, 8, 17, 37, 85, 200 };

    /** Measured worst disagreement with the hand written table, 4.2e-16. */
    private static final double SAME_AS_HAND = 1.0e-15;

    /**
     * A rooted tree, as its subtrees. The leaf is the tree with none, and the
     * order is the number of nodes.
     */
    private static final class Tree {

        final Tree[] children;
        final int order;

        Tree(Tree[] children) {
            this.children = children;
            int o = 1;
            for (int i = 0; i < children.length; ++i) {
                o += children[i].order;
            }
            this.order = o;
        }
    }

    /**
     * Every rooted tree with at most {@code maxOrder} nodes, grouped by order
     * and smallest first. A tree of order {@code n} is a root over a multiset
     * of trees of total order {@code n - 1}, so the whole set is built upwards
     * from the leaf and each tree is reached once.
     */
    private static List<Tree> trees(int maxOrder) {
        List<Tree> pool = new ArrayList<Tree>();
        pool.add(new Tree(new Tree[0]));
        for (int n = 2; n <= maxOrder; ++n) {
            List<Tree> made = new ArrayList<Tree>();
            multisets(pool, n - 1, 0, new ArrayList<Tree>(), made);
            pool.addAll(made);
        }
        return pool;
    }

    /**
     * The multisets of pooled trees whose orders sum to {@code m}, as the roots
     * they become. Only non decreasing index sequences are walked, which is
     * what makes each multiset appear exactly once rather than once per
     * permutation.
     */
    private static void multisets(List<Tree> pool, int m, int from, List<Tree> chosen, List<Tree> out) {
        if (m == 0) {
            out.add(new Tree(chosen.toArray(new Tree[chosen.size()])));
            return;
        }
        for (int i = from; i < pool.size(); ++i) {
            Tree t = pool.get(i);
            if (t.order <= m) {
                chosen.add(t);
                multisets(pool, m - t.order, i, chosen, out);
                chosen.remove(chosen.size() - 1);
            }
        }
    }

    /**
     * The elementary weight of a tree, per stage:
     * <code>psi(leaf) = 1</code> and <code>psi(t) = prod over the subtrees u of
     * (A psi(u))</code>, the product being componentwise.
     */
    private static double[] elementaryWeight(Tree t, double[][] a) {
        double[] v = new double[a.length];
        Arrays.fill(v, 1.0);
        for (int k = 0; k < t.children.length; ++k) {
            double[] au = apply(a, elementaryWeight(t.children[k], a));
            for (int i = 0; i < v.length; ++i) {
                v[i] *= au[i];
            }
        }
        return v;
    }

    /**
     * The density of a tree: its order times the densities of its subtrees. The
     * condition a method of that order must meet is that the weights dotted
     * with the elementary weight come to its reciprocal.
     */
    private static double density(Tree t) {
        double g = t.order;
        for (int k = 0; k < t.children.length; ++k) {
            g *= density(t.children[k]);
        }
        return g;
    }

    /**
     * The order conditions through {@code maxOrder}, as the amount by which
     * each is missed, in the order {@link #trees(int)} produces them -- so the
     * first {@link #THROUGH_ORDER}{@code [p - 1]} of them are the conditions of
     * order {@code p} and below.
     */
    private static double[] generated(ButcherTableau t, double[] w, int maxOrder) {
        double[][] a = t.a();
        List<Tree> trees = trees(maxOrder);
        double[] out = new double[trees.size()];
        for (int i = 0; i < trees.size(); ++i) {
            Tree tree = trees.get(i);
            out[i] = dot(w, elementaryWeight(tree, a)) - 1.0 / density(tree);
        }
        return out;
    }

    /**
     * That the generated conditions through order five are the seventeen the
     * hand written table holds. They come out in a different order, so what is
     * compared is the two sets of residuals sorted.
     */
    private static void agreeWithHand(ButcherTableau t, double[] w) {
        double[] generated = generated(t, w, 5);
        double[] hand = conditions(t, w);
        assertEquals(hand.length, generated.length);
        Arrays.sort(generated);
        Arrays.sort(hand);
        for (int i = 0; i < hand.length; ++i) {
            assertEquals(t.name() + " condition " + i, hand[i], generated[i], SAME_AS_HAND);
        }
    }
}
