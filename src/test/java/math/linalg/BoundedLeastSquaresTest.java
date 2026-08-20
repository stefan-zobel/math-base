package math.linalg;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import math.fun.DiffDMultiFunction;
import math.optim.Projection;
import math.optim.SpectralProjectedGradient;

/**
 * Tests for {@link BoundedLeastSquares}.
 * <p>
 * The correctness argument does not rest on golden values. It rests on three
 * independent statements: the exhaustive {@code 3^n} enumeration of active sets,
 * which is exact for small problems; the Karush-Kuhn-Tucker conditions
 * recomputed from the returned coefficients, which works at any size; and
 * {@link SpectralProjectedGradient}, which answers the same question by a method
 * that shares no line of code with this one.
 */
public class BoundedLeastSquaresTest {

    // -----------------------------------------------------------------
    // feasibility and the exactness of the bounds
    // -----------------------------------------------------------------

    @Test
    public void testTheSolutionIsAlwaysInsideTheBox() {
        Lcg rng = new Lcg(20260820L);
        for (int trial = 0; trial < 40; trial++) {
            int m = 30;
            int n = 8;
            double[] lower = fill(n, -0.5);
            double[] upper = fill(n, 0.75);
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(design(m, n, rng), rhs(m, 3.0, rng), lower,
                    upper);
            assertTrue("trial " + trial, r.converged);
            for (int j = 0; j < n; j++) {
                assertTrue("trial " + trial + " coefficient " + j + " = " + r.solution[j],
                        r.solution[j] >= lower[j] && r.solution[j] <= upper[j]);
            }
        }
    }

    /**
     * The property that separates this from a projected gradient method: a
     * coefficient that is held at a bound carries that bound bit for bit,
     * because it is assigned rather than approached.
     */
    @Test
    public void testACoefficientAtABoundEqualsThatBoundExactly() {
        Lcg rng = new Lcg(4711L);
        int atLower = 0;
        int atUpper = 0;
        for (int trial = 0; trial < 40; trial++) {
            int m = 25;
            int n = 10;
            double[] lower = new double[n];
            double[] upper = new double[n];
            for (int j = 0; j < n; j++) {
                lower[j] = -0.3 - 0.1 * j;
                upper[j] = 0.2 + 0.05 * j;
            }
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(design(m, n, rng), rhs(m, 4.0, rng), lower,
                    upper);
            for (int j = 0; j < n; j++) {
                if (r.activeSet[j] == BoundedLeastSquares.Bound.AT_LOWER) {
                    assertEquals("trial " + trial + " coefficient " + j, lower[j], r.solution[j], 0.0);
                    atLower++;
                } else if (r.activeSet[j] == BoundedLeastSquares.Bound.AT_UPPER) {
                    assertEquals("trial " + trial + " coefficient " + j, upper[j], r.solution[j], 0.0);
                    atUpper++;
                }
            }
        }
        assertTrue("no lower bound ever bound, the test proves nothing", atLower > 0);
        assertTrue("no upper bound ever bound, the test proves nothing", atUpper > 0);
    }

    @Test
    public void testTheActiveSetAndTheSolutionAgreeInBothDirections() {
        Lcg rng = new Lcg(987654321L);
        for (int trial = 0; trial < 30; trial++) {
            int m = 20;
            int n = 6;
            double[] lower = fill(n, 0.0);
            double[] upper = fill(n, Double.POSITIVE_INFINITY);
            BoundedLeastSquares.Result r = BoundedLeastSquares.nonNegative(design(m, n, rng), rhs(m, 2.0, rng));
            int bound = 0;
            for (int j = 0; j < n; j++) {
                boolean onBound = (r.solution[j] == lower[j]) || (r.solution[j] == upper[j]);
                boolean claimed = r.activeSet[j] != BoundedLeastSquares.Bound.FREE;
                if (claimed) {
                    bound++;
                    assertTrue("trial " + trial + " coefficient " + j + " claims a bound but sits at "
                            + r.solution[j], onBound);
                }
            }
            assertEquals("trial " + trial, bound, r.atBound);
        }
    }

    @Test
    public void testNonNegativeIsExactlyBoundedWithZeroAndInfinity() {
        Lcg rng = new Lcg(13571357L);
        for (int trial = 0; trial < 20; trial++) {
            int m = 30;
            int n = 7;
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 3.0, rng);
            BoundedLeastSquares.Result viaNonNegative = BoundedLeastSquares
                    .nonNegative(new DMatrix(m, n, a.clone()), new DMatrix(m, 1, b.clone()));
            BoundedLeastSquares.Result viaBounded = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), fill(n, 0.0), fill(n, Double.POSITIVE_INFINITY));
            // same code path, so this has to hold to the last bit
            assertArrayEquals("trial " + trial, viaNonNegative.solution, viaBounded.solution, 0.0);
            assertArrayEquals("trial " + trial, viaNonNegative.residuals, viaBounded.residuals, 0.0);
        }
    }

    @Test
    public void testAnInactiveRestrictionReproducesTheUnconstrainedFit() {
        Lcg rng = new Lcg(112233L);
        for (int trial = 0; trial < 20; trial++) {
            int m = 40;
            int n = 5;
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 1.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), fill(n, -1.0e6), fill(n, 1.0e6));
            for (int j = 0; j < n; j++) {
                assertEquals("trial " + trial + " coefficient " + j, BoundedLeastSquares.Bound.FREE,
                        r.activeSet[j]);
            }
            double[] ols = ols(new DMatrix(m, n, a.clone()), new DMatrix(m, 1, b.clone()));
            // OLS shrinks by d/(d*d) where this divides by d, so the two are
            // the same formula in different association, not the same bits
            assertArrayEquals("trial " + trial, ols, r.solution, 1.0e-12);
        }
    }

    // -----------------------------------------------------------------
    // independent oracles
    // -----------------------------------------------------------------

    /**
     * Every coefficient is at its lower bound, at its upper bound or free, so
     * enumerating all {@code 3^n} patterns, solving the reduced normal
     * equations for each and keeping the feasible one of least objective gives
     * the exact minimizer by a route that has nothing in common with an active
     * set search. The designs are well conditioned on purpose: the oracle forms
     * {@code A'A} and would otherwise be the less accurate of the two.
     */
    @Test
    public void testAgreesWithTheExhaustiveActiveSetSolution() {
        int n = 4;
        int m = 12;
        for (int variant = 0; variant < 3; variant++) {
            double[] lower;
            double[] upper;
            String name;
            if (variant == 0) {
                lower = fill(n, 0.0);
                upper = fill(n, Double.POSITIVE_INFINITY);
                name = "non-negative";
            } else if (variant == 1) {
                lower = fill(n, -1.0);
                upper = fill(n, 1.0);
                name = "box [-1, 1]";
            } else {
                lower = new double[] { 0.0, -2.0, 0.5, Double.NEGATIVE_INFINITY };
                upper = new double[] { 1.0, Double.POSITIVE_INFINITY, 0.5, 3.0 };
                name = "mixed, one pinned, two half open";
            }
            Lcg rng = new Lcg(4711L + variant);
            int compared = 0;
            for (int trial = 0; trial < 150; trial++) {
                double[] a = randomArray(m * n, 1.0, rng);
                double[] b = randomArray(m, 3.0, rng);
                double[] expected = exhaustive(a, m, n, b, lower, upper);
                if (expected == null) {
                    continue;
                }
                compared++;
                BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                        new DMatrix(m, 1, b.clone()), lower, upper);
                assertTrue(name + " trial " + trial, r.converged);
                assertArrayEquals(name + " trial " + trial, expected, r.solution, 1.0e-9);
            }
            assertTrue(name + ": the oracle never produced an answer", compared > 100);
        }
    }

    /**
     * The Karush-Kuhn-Tucker conditions of the bounded problem, recomputed from
     * the coefficients that came back: the gradient {@code w = A'(b - Ax)} has
     * to vanish where the coefficient is free, point outwards at a lower bound
     * and inwards at an upper one. Unlike the enumeration above this scales, so
     * it is the check that covers the sizes a caller would actually use.
     */
    @Test
    public void testTheKktConditionsHoldAtTheSolution() {
        int[][] shapes = { { 60, 8 }, { 200, 20 }, { 400, 40 }, { 30, 45 } };
        Lcg rng = new Lcg(20260821L);
        for (int s = 0; s < shapes.length; s++) {
            int m = shapes[s][0];
            int n = shapes[s][1];
            for (int trial = 0; trial < 10; trial++) {
                double[] a = randomArray(m * n, 1.0, rng);
                double[] b = randomArray(m, 2.0, rng);
                double[] lower = fill(n, 0.0);
                double[] upper = fill(n, Double.POSITIVE_INFINITY);
                BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                        new DMatrix(m, 1, b.clone()), lower, upper);
                assertTrue(m + "x" + n + " trial " + trial, r.converged);
                double violation = kktViolation(a, m, n, r);
                double scale = frobenius(a) * frobenius(b);
                assertTrue(m + "x" + n + " trial " + trial + " : KKT violation " + violation + " against scale "
                        + scale, violation <= 1.0e-10 * scale);
            }
        }
    }

    /**
     * {@link SpectralProjectedGradient} solves the same problem by projecting a
     * gradient step onto the same set, which shares no code with an active set
     * search. It is iterative, so the agreement is to its tolerance rather than
     * exact -- and it is the only one of the three oracles that scales.
     */
    @Test
    public void testAgreesWithSpectralProjectedGradient() {
        Lcg rng = new Lcg(31337L);
        for (int trial = 0; trial < 15; trial++) {
            int m = 50;
            int n = 8;
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 3.0, rng);
            double[] lower = fill(n, 0.0);
            double[] upper = fill(n, Double.POSITIVE_INFINITY);
            BoundedLeastSquares.Result exact = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), lower, upper);
            SpectralProjectedGradient.Result iterative = new SpectralProjectedGradient()
                    .minimize(sumOfSquares(a, m, n, b), new double[n], Projection.nonNegative(n));
            assertTrue("trial " + trial + " : " + iterative.status, iterative.converged);
            assertArrayEquals("trial " + trial, iterative.point, exact.solution, 1.0e-6);
        }
    }

    @Test
    public void testDetectionIsScaleInvariant() {
        Lcg rng = new Lcg(24680L);
        int m = 30;
        int n = 6;
        double[] a = randomArray(m * n, 1.0, rng);
        double[] b = randomArray(m, 3.0, rng);
        BoundedLeastSquares.Result base = BoundedLeastSquares.nonNegative(new DMatrix(m, n, a.clone()),
                new DMatrix(m, 1, b.clone()));
        double[] factors = { 1.0e-200, 1.0e-8, 1.0e8, 1.0e200 };
        for (int f = 0; f < factors.length; f++) {
            double factor = factors[f];
            double[] scaledA = a.clone();
            double[] scaledB = b.clone();
            for (int i = 0; i < scaledA.length; i++) {
                scaledA[i] *= factor;
            }
            for (int i = 0; i < scaledB.length; i++) {
                scaledB[i] *= factor;
            }
            BoundedLeastSquares.Result scaled = BoundedLeastSquares.nonNegative(new DMatrix(m, n, scaledA),
                    new DMatrix(m, 1, scaledB));
            // A and b scaled together leave the coefficients untouched
            for (int j = 0; j < n; j++) {
                assertEquals("factor " + factor + " coefficient " + j, base.activeSet[j], scaled.activeSet[j]);
                assertEquals("factor " + factor + " coefficient " + j, base.solution[j], scaled.solution[j],
                        1.0e-9 * Math.max(1.0, Math.abs(base.solution[j])));
            }
        }
    }

    /**
     * The non-negative orthant has corners, so the restriction produces exact
     * zeros on its own -- no penalty, no tuning parameter. This is the
     * structural difference from {@link Lasso}, and it is the reason the class
     * is useful as an estimator rather than only as a repair.
     */
    @Test
    public void testNonNegativityIsSparsifyingWithoutAnyPenalty() {
        Lcg rng = new Lcg(555777L);
        int m = 40;
        int n = 30;
        int totalZeros = 0;
        for (int trial = 0; trial < 20; trial++) {
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 1.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.nonNegative(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()));
            assertTrue("trial " + trial, r.converged);
            int zeros = 0;
            for (int j = 0; j < n; j++) {
                if (r.solution[j] == 0.0) {
                    zeros++;
                    assertEquals(BoundedLeastSquares.Bound.AT_LOWER, r.activeSet[j]);
                }
            }
            assertTrue("trial " + trial + " produced no exact zero at all", zeros > 0);
            totalZeros += zeros;
        }
        // an unconstrained fit of 30 columns produces none of these
        assertTrue("only " + totalZeros + " zeros over 20 fits of 30 coefficients", totalZeros > 20 * 5);
    }

    /**
     * The estimator doing the job it exists for: a mixture of non-negative
     * components observed through a non-negative design, with noise. The
     * coefficients have to come back near the truth, and the components that
     * are not in the mixture have to come back at exactly zero without any
     * penalty having been chosen.
     */
    @Test
    public void testItRecoversAKnownNonNegativeMixture() {
        for (int trial = 0; trial < 5; trial++) {
            Lcg rng = new Lcg(987L + trial);
            int m = 120;
            int n = 40;
            double[] a = new double[m * n];
            for (int i = 0; i < a.length; i++) {
                a[i] = Math.abs(rng.centred());
            }
            double[] truth = new double[n];
            for (int j = 0; j < n; j += 7) {
                truth[j] = 1.0 + 3.0 * Math.abs(rng.centred());
            }
            double[] b = new double[m];
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    b[i] += a[j * m + i] * truth[j];
                }
            }
            for (int i = 0; i < m; i++) {
                b[i] += 0.01 * rng.centred();
            }
            BoundedLeastSquares.Result r = BoundedLeastSquares.nonNegative(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()));
            assertTrue("trial " + trial, r.converged);
            // the noise is 0.01 per observation, so this is the level the
            // coefficients can be pinned down to at all
            assertArrayEquals("trial " + trial, truth, r.solution, 0.05);
            int zeros = 0;
            for (int j = 0; j < n; j++) {
                if (r.solution[j] == 0.0) {
                    zeros++;
                }
            }
            // no penalty, no lambda: the corners of the orthant do this
            assertTrue("trial " + trial + " produced only " + zeros + " exact zeros of " + n, zeros >= n / 2);
            // and nothing that really is in the mixture may be zeroed out
            for (int j = 0; j < n; j++) {
                if (truth[j] > 0.0) {
                    assertTrue("trial " + trial + " dropped component " + j, r.solution[j] > 0.0);
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // degenerate inputs
    // -----------------------------------------------------------------

    @Test
    public void testAPinnedCoefficientKeepsItsValueExactly() {
        Lcg rng = new Lcg(864213L);
        int m = 30;
        int n = 6;
        double[] lower = { -1.0, 0.25, -1.0, -1.0, 0.0, -1.0 };
        double[] upper = { 1.0, 0.25, 1.0, 1.0, 0.0, 1.0 };
        for (int trial = 0; trial < 20; trial++) {
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(design(m, n, rng), rhs(m, 2.0, rng), lower,
                    upper);
            assertTrue("trial " + trial, r.converged);
            assertEquals("trial " + trial, 0.25, r.solution[1], 0.0);
            assertEquals("trial " + trial, 0.0, r.solution[4], 0.0);
            assertTrue("trial " + trial, r.activeSet[1] != BoundedLeastSquares.Bound.FREE);
            assertTrue("trial " + trial, r.activeSet[4] != BoundedLeastSquares.Bound.FREE);
        }
    }

    @Test
    public void testACoefficientUnboundedOnBothSidesIsAlwaysFree() {
        Lcg rng = new Lcg(1928374L);
        int m = 30;
        int n = 5;
        double[] lower = { 0.0, Double.NEGATIVE_INFINITY, 0.0, 0.0, Double.NEGATIVE_INFINITY };
        double[] upper = { 1.0, Double.POSITIVE_INFINITY, 1.0, 1.0, Double.POSITIVE_INFINITY };
        for (int trial = 0; trial < 20; trial++) {
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(design(m, n, rng), rhs(m, 5.0, rng), lower,
                    upper);
            assertTrue("trial " + trial, r.converged);
            assertEquals("trial " + trial, BoundedLeastSquares.Bound.FREE, r.activeSet[1]);
            assertEquals("trial " + trial, BoundedLeastSquares.Bound.FREE, r.activeSet[4]);
        }
    }

    /** A right hand side pointing straight out of the box binds everything. */
    @Test
    public void testEveryCoefficientCanEndUpAtABound() {
        int m = 12;
        int n = 4;
        double[] a = new double[m * n];
        for (int j = 0; j < n; j++) {
            a[j * m + j] = 1.0;
        }
        double[] b = new double[m];
        for (int i = 0; i < n; i++) {
            b[i] = 100.0;
        }
        BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a), new DMatrix(m, 1, b),
                fill(n, 0.0), fill(n, 2.0));
        assertTrue(r.converged);
        assertEquals(n, r.atBound);
        for (int j = 0; j < n; j++) {
            assertEquals(BoundedLeastSquares.Bound.AT_UPPER, r.activeSet[j]);
            assertEquals(2.0, r.solution[j], 0.0);
        }
    }

    @Test
    public void testMoreColumnsThanRows() {
        Lcg rng = new Lcg(99991L);
        int m = 25;
        int n = 40;
        for (int trial = 0; trial < 10; trial++) {
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 1.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.nonNegative(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()));
            assertTrue("trial " + trial, r.converged);
            int free = 0;
            for (int j = 0; j < n; j++) {
                if (r.activeSet[j] == BoundedLeastSquares.Bound.FREE) {
                    free++;
                }
            }
            // the free block is factorized as m x k, so k can never exceed m
            assertTrue("trial " + trial + " : " + free + " free of " + n, free <= m);
            double violation = kktViolation(a, m, n, r);
            assertTrue("trial " + trial + " : KKT violation " + violation,
                    violation <= 1.0e-10 * frobenius(a) * frobenius(b));
        }
    }

    /**
     * A design with duplicated and linearly dependent columns has no unique
     * minimizer. Refusing it would be the wrong answer -- unmixing and spline
     * designs are collinear by nature -- so it is solved rather than rejected,
     * and the answer still satisfies the optimality conditions.
     * <p>
     * Note what does <em>not</em> happen here: the free block never becomes
     * rank deficient. An exactly duplicated column carries exactly the gradient
     * of its twin, so once one of the pair is free the other has
     * {@code w_j == 0} and can never be selected. Measured over 60 designs with
     * both a duplicated and a summed column: the free set never exceeded four
     * of the six.
     */
    @Test
    public void testARankDeficientDesignIsSolvedRatherThanRejected() {
        Lcg rng = new Lcg(741852L);
        int m = 40;
        int n = 6;
        for (int trial = 0; trial < 30; trial++) {
            double[] a = randomArray(m * n, 1.0, rng);
            System.arraycopy(a, 0, a, 5 * m, m);
            for (int i = 0; i < m; i++) {
                a[4 * m + i] = a[i] + a[2 * m + i];
            }
            double[] b = randomArray(m, 1.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), fill(n, -1.0e6), fill(n, 1.0e6));
            assertTrue("trial " + trial, r.converged);
            for (int j = 0; j < n; j++) {
                assertFalse("trial " + trial, Double.isNaN(r.solution[j]));
            }
            double violation = kktViolation(a, m, n, r);
            assertTrue("trial " + trial + " : KKT violation " + violation,
                    violation <= 1.0e-8 * frobenius(a) * frobenius(b));
        }
    }

    /**
     * The one route to a rank deficient free block, and therefore the one case
     * where {@code rankDeficient} can fire: coefficients that are unbounded on
     * both sides have no bound to be held at, so they start free, and a
     * dependent pair among them is free before the first gradient is looked at.
     */
    @Test
    public void testDependentColumnsThatStartFreeAreReportedAsRankDeficient() {
        Lcg rng = new Lcg(20260825L);
        int m = 40;
        int n = 6;
        double[] lower = fill(n, -1.0e6);
        double[] upper = fill(n, 1.0e6);
        lower[0] = Double.NEGATIVE_INFINITY;
        upper[0] = Double.POSITIVE_INFINITY;
        lower[5] = Double.NEGATIVE_INFINITY;
        upper[5] = Double.POSITIVE_INFINITY;
        for (int trial = 0; trial < 20; trial++) {
            double[] a = randomArray(m * n, 1.0, rng);
            System.arraycopy(a, 0, a, 5 * m, m);
            double[] b = randomArray(m, 1.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), lower, upper);
            assertTrue("trial " + trial, r.converged);
            assertTrue("trial " + trial + " : the deficiency went unnoticed", r.rankDeficient);
            for (int j = 0; j < n; j++) {
                assertFalse("trial " + trial, Double.isNaN(r.solution[j]));
            }
            double violation = kktViolation(a, m, n, r);
            assertTrue("trial " + trial + " : KKT violation " + violation,
                    violation <= 1.0e-8 * frobenius(a) * frobenius(b));
        }
    }

    /**
     * With every bound infinite there is nothing to hold and every coefficient
     * starts free, so the loop finds no constraint to release on its first
     * pass. It used to report convergence there and hand back the starting
     * point -- all zeros -- because the free set had never been solved at all.
     */
    @Test
    public void testEveryBoundInfiniteReproducesTheUnconstrainedFit() {
        Lcg rng = new Lcg(20260826L);
        int m = 40;
        int n = 5;
        for (int trial = 0; trial < 15; trial++) {
            double[] a = randomArray(m * n, 1.0, rng);
            double[] b = randomArray(m, 2.0, rng);
            BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(new DMatrix(m, n, a.clone()),
                    new DMatrix(m, 1, b.clone()), fill(n, Double.NEGATIVE_INFINITY),
                    fill(n, Double.POSITIVE_INFINITY));
            assertTrue("trial " + trial, r.converged);
            assertEquals("trial " + trial, 0, r.atBound);
            double[] expected = ols(new DMatrix(m, n, a.clone()), new DMatrix(m, 1, b.clone()));
            assertArrayEquals("trial " + trial, expected, r.solution, 1.0e-12);
        }
    }

    @Test
    public void testASingleColumn() {
        int m = 5;
        double[] a = { 1.0, 1.0, 1.0, 1.0, 1.0 };
        double[] b = { 2.0, 2.0, 2.0, 2.0, 2.0 };
        BoundedLeastSquares.Result free = BoundedLeastSquares.nonNegative(new DMatrix(m, 1, a.clone()),
                new DMatrix(m, 1, b.clone()));
        assertEquals(2.0, free.solution[0], 1.0e-14);
        assertEquals(BoundedLeastSquares.Bound.FREE, free.activeSet[0]);

        BoundedLeastSquares.Result capped = BoundedLeastSquares.bounded(new DMatrix(m, 1, a.clone()),
                new DMatrix(m, 1, b.clone()), new double[] { 0.0 }, new double[] { 0.5 });
        assertEquals(0.5, capped.solution[0], 0.0);
        assertEquals(BoundedLeastSquares.Bound.AT_UPPER, capped.activeSet[0]);

        double[] negated = { -2.0, -2.0, -2.0, -2.0, -2.0 };
        BoundedLeastSquares.Result clamped = BoundedLeastSquares.nonNegative(new DMatrix(m, 1, a.clone()),
                new DMatrix(m, 1, negated));
        assertEquals(0.0, clamped.solution[0], 0.0);
        assertEquals(BoundedLeastSquares.Bound.AT_LOWER, clamped.activeSet[0]);
    }

    /**
     * Budget exhaustion must never be reported as convergence. This is the one
     * failure mode the class has, so it is constructed rather than hoped for.
     */
    @Test
    public void testAnExhaustedIterationBudgetIsNotReportedAsConvergence() {
        Lcg rng = new Lcg(20260822L);
        int m = 40;
        int n = 12;
        double[] a = randomArray(m * n, 1.0, rng);
        double[] b = randomArray(m, 5.0, rng);
        Bvls.Fit generous = Bvls.solve(a, m, n, b, fill(n, 0.0), fill(n, Double.POSITIVE_INFINITY),
                Bvls.defaultMaxIterations(m, n));
        assertTrue("the problem is too easy to starve", generous.iterations > 2);
        assertTrue(generous.converged);

        Bvls.Fit starved = Bvls.solve(a, m, n, b, fill(n, 0.0), fill(n, Double.POSITIVE_INFINITY), 1);
        assertFalse(starved.converged);
        assertEquals(1, starved.iterations);
    }

    @Test
    public void testTheDefaultIterationBudgetIsNeverReached() {
        Lcg rng = new Lcg(20260823L);
        int[][] shapes = { { 60, 8 }, { 200, 20 }, { 40, 45 } };
        for (int s = 0; s < shapes.length; s++) {
            int m = shapes[s][0];
            int n = shapes[s][1];
            int cap = Bvls.defaultMaxIterations(m, n);
            int worst = 0;
            for (int trial = 0; trial < 20; trial++) {
                double[] a = randomArray(m * n, 1.0, rng);
                double[] b = randomArray(m, 3.0, rng);
                Bvls.Fit fit = Bvls.solve(a, m, n, b, fill(n, -0.2), fill(n, 0.2), cap);
                assertTrue(m + "x" + n + " trial " + trial, fit.converged);
                worst = Math.max(worst, fit.iterations);
            }
            // the cap is meant to be a backstop, not a working limit
            assertTrue(m + "x" + n + " used " + worst + " of " + cap, worst < cap / 2);
        }
    }

    // -----------------------------------------------------------------
    // argument validation
    // -----------------------------------------------------------------

    @Test
    public void testMismatchedShapesAreRejected() {
        expectIae("rows disagree", new Runnable() {
            public void run() {
                BoundedLeastSquares.nonNegative(new DMatrix(5, 2, new double[10]), new DMatrix(4, 1, new double[4]));
            }
        });
        expectIae("b has two columns", new Runnable() {
            public void run() {
                BoundedLeastSquares.nonNegative(new DMatrix(5, 2, new double[10]), new DMatrix(5, 2, new double[10]));
            }
        });
        expectIae("bounds too short", new Runnable() {
            public void run() {
                BoundedLeastSquares.bounded(new DMatrix(5, 2, new double[10]), new DMatrix(5, 1, new double[5]),
                        new double[1], new double[2]);
            }
        });
    }

    @Test
    public void testNonFiniteAndImpossibleInputIsRejected() {
        final double[] a = new double[10];
        expectIae("NaN in A", new Runnable() {
            public void run() {
                double[] bad = a.clone();
                bad[7] = Double.NaN;
                BoundedLeastSquares.nonNegative(new DMatrix(5, 2, bad), new DMatrix(5, 1, new double[5]));
            }
        });
        expectIae("infinity in b", new Runnable() {
            public void run() {
                double[] bad = new double[5];
                bad[2] = Double.POSITIVE_INFINITY;
                BoundedLeastSquares.nonNegative(new DMatrix(5, 2, a.clone()), new DMatrix(5, 1, bad));
            }
        });
        expectIae("lower > upper", new Runnable() {
            public void run() {
                BoundedLeastSquares.bounded(new DMatrix(5, 2, a.clone()), new DMatrix(5, 1, new double[5]),
                        new double[] { 1.0, 0.0 }, new double[] { 0.0, 1.0 });
            }
        });
        expectIae("NaN bound", new Runnable() {
            public void run() {
                BoundedLeastSquares.bounded(new DMatrix(5, 2, a.clone()), new DMatrix(5, 1, new double[5]),
                        new double[] { Double.NaN, 0.0 }, new double[] { 1.0, 1.0 });
            }
        });
        expectIae("more doubly unbounded coefficients than rows", new Runnable() {
            public void run() {
                int m = 2;
                int n = 3;
                BoundedLeastSquares.bounded(new DMatrix(m, n, new double[m * n]), new DMatrix(m, 1, new double[m]),
                        fill(n, Double.NEGATIVE_INFINITY), fill(n, Double.POSITIVE_INFINITY));
            }
        });
    }

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /**
     * The exact minimizer over the box, by trying every assignment of every
     * coefficient to lower bound, free or upper bound, solving the reduced
     * normal equations and keeping the feasible pattern of least objective.
     */
    private static double[] exhaustive(double[] a, int m, int n, double[] b, double[] lower, double[] upper) {
        double[][] q = new double[n][n];
        double[] c = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int t = 0; t < m; t++) {
                    s += a[i * m + t] * a[j * m + t];
                }
                q[i][j] = s;
            }
            double s = 0.0;
            for (int t = 0; t < m; t++) {
                s += a[i * m + t] * b[t];
            }
            c[i] = -s;
        }
        int patterns = 1;
        for (int i = 0; i < n; i++) {
            patterns *= 3;
        }
        double bestObjective = Double.MAX_VALUE;
        double[] best = null;
        for (int pattern = 0; pattern < patterns; pattern++) {
            int[] state = new int[n];
            int rest = pattern;
            boolean usable = true;
            for (int i = 0; i < n; i++) {
                state[i] = rest % 3;
                rest /= 3;
                if (state[i] == 0 && Double.isInfinite(lower[i])) {
                    usable = false;
                }
                if (state[i] == 2 && Double.isInfinite(upper[i])) {
                    usable = false;
                }
            }
            if (!usable) {
                continue;
            }
            double[] x = new double[n];
            int freeCount = 0;
            for (int i = 0; i < n; i++) {
                if (state[i] == 0) {
                    x[i] = lower[i];
                } else if (state[i] == 2) {
                    x[i] = upper[i];
                } else {
                    freeCount++;
                }
            }
            if (freeCount > 0) {
                int[] index = new int[freeCount];
                int k = 0;
                for (int i = 0; i < n; i++) {
                    if (state[i] == 1) {
                        index[k++] = i;
                    }
                }
                double[][] reduced = new double[freeCount][freeCount];
                double[] rightHandSide = new double[freeCount];
                for (int r = 0; r < freeCount; r++) {
                    for (int cc = 0; cc < freeCount; cc++) {
                        reduced[r][cc] = q[index[r]][index[cc]];
                    }
                    double s = c[index[r]];
                    for (int i = 0; i < n; i++) {
                        if (state[i] != 1) {
                            s += q[index[r]][i] * x[i];
                        }
                    }
                    rightHandSide[r] = -s;
                }
                double[] solution = solveLinear(reduced, rightHandSide);
                if (solution == null) {
                    continue;
                }
                for (int r = 0; r < freeCount; r++) {
                    x[index[r]] = solution[r];
                }
            }
            boolean feasible = true;
            for (int i = 0; i < n && feasible; i++) {
                if (x[i] < lower[i] - 1.0e-9 || x[i] > upper[i] + 1.0e-9) {
                    feasible = false;
                }
            }
            if (!feasible) {
                continue;
            }
            double objective = 0.0;
            for (int t = 0; t < m; t++) {
                double s = -b[t];
                for (int j = 0; j < n; j++) {
                    s += a[j * m + t] * x[j];
                }
                objective += s * s;
            }
            if (objective < bestObjective) {
                bestObjective = objective;
                best = x;
            }
        }
        return best;
    }

    /** Gauss-Jordan with partial pivoting, for the oracle only. */
    private static double[] solveLinear(double[][] matrix, double[] rightHandSide) {
        int n = rightHandSide.length;
        double[][] a = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, a[i], 0, n);
            a[i][n] = rightHandSide[i];
        }
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(a[r][col]) > Math.abs(a[pivot][col])) {
                    pivot = r;
                }
            }
            if (Math.abs(a[pivot][col]) < 1.0e-13) {
                return null;
            }
            double[] swap = a[col];
            a[col] = a[pivot];
            a[pivot] = swap;
            for (int r = 0; r < n; r++) {
                if (r == col) {
                    continue;
                }
                double factor = a[r][col] / a[col][col];
                for (int cc = col; cc <= n; cc++) {
                    a[r][cc] -= factor * a[col][cc];
                }
            }
        }
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = a[i][n] / a[i][i];
        }
        return out;
    }

    /**
     * The largest amount by which any Karush-Kuhn-Tucker condition is broken:
     * {@code |w_j|} where the coefficient is free, {@code max(0, w_j)} at a
     * lower bound and {@code max(0, -w_j)} at an upper one.
     */
    private static double kktViolation(double[] a, int m, int n, BoundedLeastSquares.Result r) {
        double worst = 0.0;
        for (int j = 0; j < n; j++) {
            double w = 0.0;
            int column = j * m;
            for (int i = 0; i < m; i++) {
                w += a[column + i] * r.residuals[i];
            }
            double violation;
            if (r.activeSet[j] == BoundedLeastSquares.Bound.FREE) {
                violation = Math.abs(w);
            } else if (r.activeSet[j] == BoundedLeastSquares.Bound.AT_LOWER) {
                violation = Math.max(0.0, w);
            } else {
                violation = Math.max(0.0, -w);
            }
            worst = Math.max(worst, violation);
        }
        return worst;
    }

    private static double[] ols(DMatrix x, DMatrix y) {
        LSSummary summary = OLS.estimate(0.05, x, y);
        double[] beta = new double[x.numColumns()];
        for (int j = 0; j < beta.length; j++) {
            beta[j] = summary.getBeta().get(j);
        }
        return beta;
    }

    private static DiffDMultiFunction sumOfSquares(final double[] a, final int m, final int n, final double[] b) {
        return new DiffDMultiFunction() {

            @Override
            public double apply(double[] x) {
                double sum = 0.0;
                for (int i = 0; i < m; i++) {
                    double r = -b[i];
                    for (int j = 0; j < n; j++) {
                        r += a[j * m + i] * x[j];
                    }
                    sum += r * r;
                }
                return sum;
            }

            @Override
            public void derivativeAt(double[] x, double[] gradient) {
                Arrays.fill(gradient, 0.0);
                for (int i = 0; i < m; i++) {
                    double r = -b[i];
                    for (int j = 0; j < n; j++) {
                        r += a[j * m + i] * x[j];
                    }
                    for (int j = 0; j < n; j++) {
                        gradient[j] += 2.0 * r * a[j * m + i];
                    }
                }
            }
        };
    }

    private static double frobenius(double[] v) {
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            sum += v[i] * v[i];
        }
        return Math.sqrt(sum);
    }

    private static double[] fill(int n, double value) {
        double[] a = new double[n];
        Arrays.fill(a, value);
        return a;
    }

    private static double[] randomArray(int length, double amplitude, Lcg rng) {
        double[] a = new double[length];
        for (int i = 0; i < length; i++) {
            a[i] = amplitude * rng.centred();
        }
        return a;
    }

    private static DMatrix design(int m, int n, Lcg rng) {
        return new DMatrix(m, n, randomArray(m * n, 1.0, rng));
    }

    private static DMatrix rhs(int m, double amplitude, Lcg rng) {
        return new DMatrix(m, 1, randomArray(m, amplitude, rng));
    }

    private static void expectIae(String what, Runnable r) {
        try {
            r.run();
            fail("expected an IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }

        double centred() {
            return 2.0 * next() - 1.0;
        }
    }
}
