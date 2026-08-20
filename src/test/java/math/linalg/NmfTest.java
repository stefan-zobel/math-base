package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Tests for {@link Nmf}.
 * <p>
 * The problem is not convex, so there is no golden answer to compare against and
 * the correctness argument is built out of invariants and independent solvers
 * instead. Three of those carry it: every column of {@code H} at a converged point
 * <em>is</em> the non-negative least squares solution given {@code W}, which
 * {@link BoundedLeastSquares} answers exactly by an active set method sharing no
 * line of code with the engine here; the rank-one factorization of a non-negative
 * matrix is its leading singular triplet by Perron-Frobenius, which
 * {@link FlatParallelJacobiSVD} supplies; and the objective of an exact coordinate
 * minimization can never increase.
 */
public class NmfTest {

    // -----------------------------------------------------------------
    // structure of the returned factors
    // -----------------------------------------------------------------

    @Test
    public void testBothFactorsAreNonNegativeAndTheirZerosAreExact() {
        double[] x = parts(120, 70, 6, 0.05, 11L);
        Nmf.Result r = new Nmf().factor(x, 120, 70, 6);
        for (int i = 0; i < r.W.length; i++) {
            assertTrue("W[" + i + "] = " + r.W[i], r.W[i] >= 0.0);
        }
        for (int i = 0; i < r.H.length; i++) {
            assertTrue("H[" + i + "] = " + r.H[i], r.H[i] >= 0.0);
        }
        // a penalized fit has genuine zeros, and they are exactly zero
        Nmf.Result sparse = new Nmf(1.0e-6, 2000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.lasso(0.5)).factor(x, 120, 70, 6);
        int zeros = 0;
        for (int i = 0; i < sparse.H.length; i++) {
            assertTrue(sparse.H[i] >= 0.0);
            if (sparse.H[i] == 0.0) {
                zeros++;
            }
        }
        assertTrue("expected exact zeros under an L1 penalty, saw " + zeros, zeros > 0);
    }

    @Test
    public void testEveryColumnOfWHasUnitNorm() {
        int m = 150;
        int n = 80;
        int k = 7;
        double[] x = parts(m, n, k, 0.05, 12L);
        for (Nmf.Init init : Nmf.Init.values()) {
            Nmf.Result r = new Nmf(1.0e-5, 500, init, 3L, Nmf.Penalty.none(), Nmf.Penalty.none())
                    .factor(x, m, n, k);
            for (int c = 0; c < k; c++) {
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    s += r.W[c * m + i] * r.W[c * m + i];
                }
                assertEquals(init + ", column " + c, 1.0, Math.sqrt(s), 1.0e-14);
            }
        }
    }

    @Test
    public void testTheComponentEnergiesAreTheRowNormsOfHAndDescend() {
        int m = 100;
        int n = 90;
        int k = 6;
        double[] x = parts(m, n, k, 0.05, 13L);
        Nmf.Result r = new Nmf().factor(x, m, n, k);
        for (int c = 0; c < k; c++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                s += r.H[j * k + c] * r.H[j * k + c];
            }
            assertEquals("component " + c, Math.sqrt(s), r.componentEnergy[c], 1.0e-12);
            if (c > 0) {
                assertTrue("energies must descend", r.componentEnergy[c] <= r.componentEnergy[c - 1]);
            }
        }
    }

    @Test
    public void testTheReportedErrorIsTheOneReconstructProduces() {
        int m = 90;
        int n = 110;
        int k = 5;
        double[] x = parts(m, n, k, 0.05, 14L);
        Nmf.Result r = new Nmf().factor(x, m, n, k);
        assertEquals(r.relativeError, Nmf.reconstructionError(x, r), 1.0e-15);

        double[] product = Nmf.reconstruct(r);
        double num = 0.0;
        double den = 0.0;
        for (int i = 0; i < x.length; i++) {
            num += (x[i] - product[i]) * (x[i] - product[i]);
            den += x[i] * x[i];
        }
        assertEquals(r.relativeError, Math.sqrt(num / den), 1.0e-12);
    }

    // -----------------------------------------------------------------
    // the independent oracles
    // -----------------------------------------------------------------

    @Test
    public void testEveryColumnOfHIsTheNonNegativeLeastSquaresSolutionGivenW() {
        int m = 120;
        int n = 50;
        int k = 6;
        double[] x = parts(m, n, k, 0.05, 8200L);
        Nmf.Result r = new Nmf(1.0e-12, 20000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);

        DMatrix w = new DMatrix(m, k, r.W.clone());
        double[] column = new double[m];
        for (int j = 0; j < n; j++) {
            System.arraycopy(x, j * m, column, 0, m);
            BoundedLeastSquares.Result nnls = BoundedLeastSquares.nonNegative(w,
                    new DMatrix(m, 1, column.clone()));
            double scale = 0.0;
            for (int c = 0; c < k; c++) {
                scale = Math.max(scale, Math.abs(nnls.solution[c]));
            }
            for (int c = 0; c < k; c++) {
                // the active sets must agree exactly: which parts a sample uses is the
                // statement, and it is discrete
                boolean freeHere = r.H[j * k + c] > 0.0;
                boolean freeThere = nnls.activeSet[c] == BoundedLeastSquares.Bound.FREE;
                assertEquals("column " + j + ", component " + c, freeThere, freeHere);
                assertEquals("column " + j + ", component " + c, nnls.solution[c], r.H[j * k + c],
                        1.0e-6 * Math.max(scale, 1.0));
            }
        }
    }

    @Test
    public void testEveryRowOfWIsTheNonNegativeLeastSquaresSolutionGivenH() {
        int m = 80;
        int n = 60;
        int k = 5;
        double[] x = parts(m, n, k, 0.05, 8300L);
        Nmf.Result r = new Nmf(1.0e-12, 20000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);

        double[] ht = new double[n * k];
        for (int c = 0; c < k; c++) {
            for (int j = 0; j < n; j++) {
                ht[c * n + j] = r.H[j * k + c];
            }
        }
        DMatrix h = new DMatrix(n, k, ht);
        double[] row = new double[n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                row[j] = x[j * m + i];
            }
            BoundedLeastSquares.Result nnls = BoundedLeastSquares.nonNegative(h,
                    new DMatrix(n, 1, row.clone()));
            double scale = 0.0;
            for (int c = 0; c < k; c++) {
                scale = Math.max(scale, Math.abs(nnls.solution[c]));
            }
            for (int c = 0; c < k; c++) {
                assertEquals("row " + i + ", component " + c, nnls.solution[c], r.W[c * m + i],
                        1.0e-6 * Math.max(scale, 1.0));
            }
        }
    }

    @Test
    public void testTheRankOneFactorizationIsTheLeadingSingularTriplet() {
        int[][] shapes = { { 120, 60, 5 }, { 60, 140, 5 }, { 40, 40, 4 } };
        for (int[] shape : shapes) {
            int m = shape[0];
            int n = shape[1];
            double[] x = parts(m, n, shape[2], 0.05, 313L + m);
            Nmf.Result r = new Nmf(1.0e-12, 20000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                    Nmf.Penalty.none()).factor(x, m, n, 1);

            // the leading singular triplet of a non-negative matrix is non-negative
            // (Perron-Frobenius), so the rank-one NMF must be the rank-one SVD
            boolean wide = m < n;
            double[] a;
            int rows = wide ? n : m;
            int cols = wide ? m : n;
            if (wide) {
                a = new double[n * m];
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < m; i++) {
                        a[i * n + j] = x[j * m + i];
                    }
                }
            } else {
                a = x.clone();
            }
            FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a, rows, cols);
            double[] left = new double[m];
            double[] right = new double[n];
            for (int i = 0; i < rows; i++) {
                if (wide) {
                    right[i] = Math.abs(svd.U[i]);
                } else {
                    left[i] = Math.abs(svd.U[i]);
                }
            }
            for (int i = 0; i < cols; i++) {
                if (wide) {
                    left[i] = Math.abs(svd.V[i]);
                } else {
                    right[i] = Math.abs(svd.V[i]);
                }
            }
            double diff = 0.0;
            double scale = 0.0;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    double mine = r.W[i] * r.H[j];
                    double theirs = svd.sigma[0] * left[i] * right[j];
                    diff += (mine - theirs) * (mine - theirs);
                    scale += theirs * theirs;
                }
            }
            assertEquals(m + "x" + n, 0.0, Math.sqrt(diff / scale), 1.0e-12);
        }
    }

    @Test
    public void testTheObjectiveNeverIncreases() {
        int[][] shapes = { { 80, 50, 5 }, { 60, 90, 4 } };
        for (int[] shape : shapes) {
            int m = shape[0];
            int n = shape[1];
            int k = shape[2];
            for (int noisy = 0; noisy < 2; noisy++) {
                double[] x = parts(m, n, k, (noisy == 1) ? 0.05 : 0.0, 606L + m);
                double normSquared = 0.0;
                for (int i = 0; i < x.length; i++) {
                    normSquared += x[i] * x[i];
                }
                // a tolerance of zero is never met, so a cap of t runs exactly t sweeps
                double previous = Double.MAX_VALUE;
                for (int t = 1; t <= 60; t++) {
                    Nmf.Result r = new Nmf(0.0, t, Nmf.Init.RANDOM, 4242L, Nmf.Penalty.none(),
                            Nmf.Penalty.none()).factor(x, m, n, k);
                    assertTrue("sweep " + t + " of " + m + "x" + n + ": " + r.objective + " > " + previous,
                            r.objective <= previous + 1.0e-28 * normSquared);
                    previous = r.objective;
                }
            }
        }
    }

    @Test
    public void testTheProjectedGradientVanishesAtTheReturnedFactors() {
        int m = 100;
        int n = 70;
        int k = 6;
        double[] x = parts(m, n, k, 0.05, 4400L);
        Nmf.Result r = new Nmf(1.0e-10, 20000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);

        // grad_H = W'W H - W'X, recomputed here from the returned factors alone
        double scale = 0.0;
        double worst = 0.0;
        for (int j = 0; j < n; j++) {
            for (int c = 0; c < k; c++) {
                double g = 0.0;
                for (int i = 0; i < m; i++) {
                    double fit = 0.0;
                    for (int q = 0; q < k; q++) {
                        fit += r.W[q * m + i] * r.H[j * k + q];
                    }
                    g += r.W[c * m + i] * (fit - x[j * m + i]);
                    scale = Math.max(scale, Math.abs(x[j * m + i]));
                }
                double projected = (r.H[j * k + c] > 0.0) ? g : Math.min(0.0, g);
                worst = Math.max(worst, Math.abs(projected));
            }
        }
        assertEquals("projected gradient of H", 0.0, worst / scale, 1.0e-6);
    }

    // -----------------------------------------------------------------
    // what the factorization recovers
    // -----------------------------------------------------------------

    @Test
    public void testExactlyLowRankDataIsRecovered() {
        int m = 200;
        int n = 120;
        int k = 8;
        double[] x = parts(m, n, k, 0.0, 55L);
        Nmf.Result r = new Nmf(1.0e-9, 20000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);
        assertTrue("relative error " + r.relativeError, r.relativeError < 1.0e-8);
        assertTrue(r.converged);
        assertEquals(0, r.deadComponents);
    }

    @Test
    public void testTheErrorFallsSharplyAtTheTrueRankAndIsFlatAfterwards() {
        int m = 180;
        int n = 100;
        int trueRank = 10;
        double[] x = parts(m, n, trueRank, 0.02, 1212L);
        double[] error = new double[16];
        for (int k = 1; k <= 15; k++) {
            error[k] = new Nmf(1.0e-5, 2000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(), Nmf.Penalty.none())
                    .factor(x, m, n, k).relativeError;
            if (k > 1) {
                assertTrue("error must not grow with the rank, k = " + k, error[k] <= error[k - 1]);
            }
        }
        assertTrue("drop at the true rank was " + (error[trueRank] / error[trueRank - 1]),
                error[trueRank] / error[trueRank - 1] < 0.1);
        for (int k = trueRank + 1; k <= 15; k++) {
            assertTrue("k = " + k + " should add almost nothing", error[k] / error[k - 1] > 0.9);
        }
    }

    @Test
    public void testDifferentSeedsReachTheSameErrorAtTheDefaultTolerance() {
        int m = 200;
        int n = 100;
        int k = 8;
        double[] x = parts(m, n, k, 0.02, 20260822L + m);
        double best = Double.MAX_VALUE;
        double worst = 0.0;
        for (int t = 0; t < 8; t++) {
            Nmf.Result r = new Nmf(1.0e-5, 5000, Nmf.Init.RANDOM, 1000L + 37 * t, Nmf.Penalty.none(),
                    Nmf.Penalty.none()).factor(x, m, n, k);
            assertTrue(r.converged);
            assertEquals(0, r.deadComponents);
            best = Math.min(best, r.relativeError);
            worst = Math.max(worst, r.relativeError);
        }
        assertTrue("spread over seeds was " + ((worst - best) / best), (worst - best) / best < 1.0e-3);
    }

    @Test
    public void testTheNndsvdStartIgnoresTheSeed() {
        int m = 90;
        int n = 60;
        int k = 5;
        double[] x = parts(m, n, k, 0.05, 66L);
        Nmf.Result a = new Nmf(1.0e-6, 2000, Nmf.Init.NNDSVD, 1L, Nmf.Penalty.none(), Nmf.Penalty.none())
                .factor(x, m, n, k);
        Nmf.Result b = new Nmf(1.0e-6, 2000, Nmf.Init.NNDSVD, 987654321L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);
        assertEquals(a.iterations, b.iterations);
        for (int i = 0; i < a.W.length; i++) {
            assertEquals(bits(a.W[i]), bits(b.W[i]));
        }
        for (int i = 0; i < a.H.length; i++) {
            assertEquals(bits(a.H[i]), bits(b.H[i]));
        }
    }

    @Test
    public void testAWideMatrixIsFactoredAsWellAsItsTranspose() {
        int m = 60;
        int n = 150;
        int k = 5;
        double[] x = parts(m, n, k, 0.03, 77L);
        double[] xt = new double[n * m];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                xt[i * n + j] = x[j * m + i];
            }
        }
        Nmf nmf = new Nmf(1.0e-7, 5000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(), Nmf.Penalty.none());
        Nmf.Result wide = nmf.factor(x, m, n, k);
        Nmf.Result tall = nmf.factor(xt, n, m, k);
        assertEquals("the same data either way round", wide.relativeError, tall.relativeError, 1.0e-9);
    }

    // -----------------------------------------------------------------
    // scale
    // -----------------------------------------------------------------

    @Test
    public void testTheFactorizationIsEquivariantUnderARescalingOfTheData() {
        int m = 150;
        int n = 90;
        int k = 6;
        double[] x = parts(m, n, k, 0.03, 246L);
        for (Nmf.Init init : Nmf.Init.values()) {
            Nmf nmf = new Nmf(1.0e-6, 3000, init, 99L, Nmf.Penalty.none(), Nmf.Penalty.none());
            Nmf.Result base = nmf.factor(x, m, n, k);
            for (int exponent : new int[] { 10, -12, 20 }) {
                double c = Math.scalb(1.0, exponent);
                double[] scaled = new double[x.length];
                for (int i = 0; i < x.length; i++) {
                    scaled[i] = c * x[i];
                }
                Nmf.Result r = nmf.factor(scaled, m, n, k);
                // a power of two rescales exactly, so W is untouched and H scales with it
                assertEquals(init + ", c = " + c, bits(base.relativeError), bits(r.relativeError));
                assertEquals(init + ", c = " + c, base.iterations, r.iterations);
                for (int i = 0; i < base.W.length; i++) {
                    assertEquals(bits(base.W[i]), bits(r.W[i]));
                }
                for (int i = 0; i < base.H.length; i++) {
                    assertEquals(bits(c * base.H[i]), bits(r.H[i]));
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // the penalties
    // -----------------------------------------------------------------

    @Test
    public void testAPenaltyOfStrengthZeroReproducesTheUnpenalizedRunBitForBit() {
        int m = 120;
        int n = 80;
        int k = 6;
        double[] x = parts(m, n, k, 0.03, 369L);
        Nmf.Result plain = new Nmf(1.0e-6, 200, Nmf.Init.RANDOM, 7L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, k);
        Nmf.Result zero = new Nmf(1.0e-6, 200, Nmf.Init.RANDOM, 7L, Nmf.Penalty.elasticNet(0.0, 0.5),
                Nmf.Penalty.elasticNet(0.0, 0.5)).factor(x, m, n, k);
        assertEquals(plain.iterations, zero.iterations);
        assertEquals(bits(plain.objective), bits(zero.objective));
        for (int i = 0; i < plain.W.length; i++) {
            assertEquals(bits(plain.W[i]), bits(zero.W[i]));
        }
        for (int i = 0; i < plain.H.length; i++) {
            assertEquals(bits(plain.H[i]), bits(zero.H[i]));
        }
    }

    @Test
    public void testAnL1PenaltyMakesHSparserAndTheFitWorseMonotonically() {
        int m = 200;
        int n = 120;
        int k = 8;
        double[] x = parts(m, n, k, 0.03, 3690L);
        double[] lambdas = { 0.0, 0.001, 0.01, 0.1, 1.0, 10.0 };
        int previousZeros = -1;
        double previousError = -1.0;
        double previousSum = Double.MAX_VALUE;
        for (int t = 0; t < lambdas.length; t++) {
            Nmf.Result r = new Nmf(1.0e-6, 5000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                    Nmf.Penalty.lasso(lambdas[t])).factor(x, m, n, k);
            int zeros = 0;
            double sum = 0.0;
            for (int i = 0; i < r.H.length; i++) {
                sum += r.H[i];
                if (r.H[i] == 0.0) {
                    zeros++;
                }
            }
            assertTrue("zeros must not fall, lambda = " + lambdas[t], zeros >= previousZeros);
            assertTrue("the fit must not improve, lambda = " + lambdas[t],
                    r.relativeError >= previousError);
            assertTrue("sum H must not grow, lambda = " + lambdas[t], sum <= previousSum);
            previousZeros = zeros;
            previousError = r.relativeError;
            previousSum = sum;
        }
        assertTrue("an L1 penalty must produce zeros", previousZeros > 0);
    }

    @Test
    public void testAnL2PenaltyShrinksHWithoutTheFitImproving() {
        int m = 150;
        int n = 100;
        int k = 6;
        double[] x = parts(m, n, k, 0.03, 3691L);
        double previousSquares = Double.MAX_VALUE;
        double previousError = -1.0;
        for (double lambda : new double[] { 0.0, 0.01, 0.1, 1.0, 10.0 }) {
            Nmf.Result r = new Nmf(1.0e-6, 5000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                    Nmf.Penalty.ridge(lambda)).factor(x, m, n, k);
            double squares = 0.0;
            for (int i = 0; i < r.H.length; i++) {
                squares += r.H[i] * r.H[i];
            }
            assertTrue("||H||^2 must not grow, lambda = " + lambda, squares <= previousSquares);
            assertTrue("the fit must not improve, lambda = " + lambda, r.relativeError >= previousError);
            previousSquares = squares;
            previousError = r.relativeError;
        }
    }

    @Test
    public void testAPenaltyCannotBeEvadedByInflatingTheOtherFactor() {
        // W H = (W D)(D^-1 H) leaves the fit alone, so an unnormalized run could
        // answer an L1 penalty on H by growing W instead. The columns of W must stay
        // at unit norm however strong the penalty is.
        int m = 200;
        int n = 120;
        int k = 8;
        double[] x = parts(m, n, k, 0.03, 3692L);
        for (double lambda : new double[] { 0.01, 1.0, 100.0 }) {
            Nmf.Result r = new Nmf(1.0e-6, 5000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                    Nmf.Penalty.lasso(lambda)).factor(x, m, n, k);
            for (int c = 0; c < k; c++) {
                double squares = 0.0;
                for (int i = 0; i < m; i++) {
                    squares += r.W[c * m + i] * r.W[c * m + i];
                }
                double norm = Math.sqrt(squares);
                assertTrue("lambda = " + lambda + ", column " + c + " grew to " + norm,
                        norm == 0.0 || Math.abs(norm - 1.0) < 1.0e-12);
            }
        }
    }

    // -----------------------------------------------------------------
    // degenerate input
    // -----------------------------------------------------------------

    @Test
    public void testAnAllZeroMatrix() {
        int m = 60;
        int n = 40;
        Nmf.Result r = new Nmf().factor(new double[m * n], m, n, 3);
        assertEquals(0.0, r.relativeError, 0.0);
        assertEquals(0.0, r.objective, 0.0);
        assertTrue(r.converged);
        assertEquals(3, r.deadComponents);
        double[] product = Nmf.reconstruct(r);
        for (int i = 0; i < product.length; i++) {
            assertEquals(0.0, product[i], 0.0);
        }
    }

    @Test
    public void testAZeroRowAndAZeroColumnAreFittedWithoutTrouble() {
        int m = 60;
        int n = 40;
        int k = 5;
        double[] x = parts(m, n, k, 0.05, 1L);
        for (int j = 0; j < n; j++) {
            x[j * m + 7] = 0.0;
        }
        for (int i = 0; i < m; i++) {
            x[3 * m + i] = 0.0;
        }
        Nmf.Result r = new Nmf().factor(x, m, n, k);
        assertTrue(r.converged);
        double[] product = Nmf.reconstruct(r);
        for (int i = 0; i < product.length; i++) {
            assertTrue(product[i] >= 0.0);
        }
    }

    @Test
    public void testTheExtremeRanks() {
        int m = 60;
        int n = 40;
        double[] x = parts(m, n, 5, 0.05, 2L);
        Nmf.Result one = new Nmf().factor(x, m, n, 1);
        assertEquals(1, one.k);
        assertTrue(one.relativeError > 0.0 && one.relativeError < 1.0);
        Nmf.Result full = new Nmf(1.0e-5, 2000, Nmf.Init.NNDSVD, 0L, Nmf.Penalty.none(),
                Nmf.Penalty.none()).factor(x, m, n, Math.min(m, n));
        assertTrue("more components cannot fit worse", full.relativeError <= one.relativeError);
    }

    @Test
    public void testAnExhaustedBudgetIsNotReportedAsConvergence() {
        int m = 60;
        int n = 40;
        double[] x = parts(m, n, 5, 0.05, 3L);
        Nmf.Result r = new Nmf(1.0e-14, 2, Nmf.Init.RANDOM, 1L, Nmf.Penalty.none(), Nmf.Penalty.none())
                .factor(x, m, n, 5);
        assertEquals(2, r.iterations);
        assertEquals(Nmf.Status.TOO_MANY_ITERATIONS, r.status);
        assertFalse(r.converged);
        assertFalse(r.status.isSuccess());
    }

    @Test
    public void testTheMatrixEntryPointAgreesWithTheFlatOne() {
        int m = 80;
        int n = 50;
        int k = 5;
        double[] x = parts(m, n, k, 0.05, 4L);
        Nmf nmf = new Nmf();
        Nmf.Result flat = nmf.factor(x, m, n, k);
        Nmf.Result wrapped = nmf.factor(new DMatrix(m, n, x.clone()), k);
        assertEquals(bits(flat.relativeError), bits(wrapped.relativeError));
        for (int i = 0; i < flat.H.length; i++) {
            assertEquals(bits(flat.H[i]), bits(wrapped.H[i]));
        }
    }

    // -----------------------------------------------------------------
    // rejections
    // -----------------------------------------------------------------

    @Test
    public void testMalformedDataIsRejected() {
        int m = 6;
        int n = 4;
        double[] x = parts(m, n, 2, 0.0, 5L);
        rejects(x, m, n, 0);
        rejects(x, m, n, 5);
        rejects(x, m, n, -1);
        rejects(x, m + 1, n, 2);
        rejects(null, m, n, 2);
        rejects(new double[0], 0, 0, 1);

        double[] negative = x.clone();
        negative[7] = -1.0e-18;
        rejects(negative, m, n, 2);
        double[] notANumber = x.clone();
        notANumber[3] = Double.NaN;
        rejects(notANumber, m, n, 2);
        double[] infinite = x.clone();
        infinite[2] = Double.POSITIVE_INFINITY;
        rejects(infinite, m, n, 2);
    }

    @Test
    public void testMalformedConfigurationIsRejected() {
        try {
            new Nmf(-1.0, 10, Nmf.Init.RANDOM, 0L, Nmf.Penalty.none(), Nmf.Penalty.none());
            fail("a negative tolerance must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            new Nmf(1.0e-5, 0, Nmf.Init.RANDOM, 0L, Nmf.Penalty.none(), Nmf.Penalty.none());
            fail("a cap of zero must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            new Nmf(1.0e-5, 10, null, 0L, Nmf.Penalty.none(), Nmf.Penalty.none());
            fail("a missing start must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            new Nmf(1.0e-5, 10, Nmf.Init.RANDOM, 0L, null, Nmf.Penalty.none());
            fail("a missing penalty must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            Nmf.Penalty.lasso(-1.0);
            fail("a negative lambda must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            Nmf.Penalty.elasticNet(1.0, 1.5);
            fail("an alpha outside [0, 1] must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            Nmf.Penalty.lasso(Double.POSITIVE_INFINITY);
            fail("an infinite lambda must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    private static void rejects(double[] x, int m, int n, int k) {
        try {
            new Nmf().factor(x, m, n, k);
            fail("expected a rejection for " + m + "x" + n + ", k = " + k);
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    private static long bits(double value) {
        return Double.doubleToRawLongBits(value);
    }

    /**
     * Parts-based non-negative data of a known rank, column-major. Each factor is
     * supported on a slice of the rows, which is the situation NMF is meant for.
     */
    private static double[] parts(int m, int n, int k, double noise, long seed) {
        Lcg rng = new Lcg(seed);
        double[] w = new double[m * k];
        double[] h = new double[k * n];
        for (int r = 0; r < k; r++) {
            for (int i = 0; i < m; i++) {
                w[r * m + i] = (i % k == r) ? 1.0 + rng.next() : ((rng.next() < 0.1) ? rng.next() * 0.2 : 0.0);
            }
        }
        for (int j = 0; j < n; j++) {
            for (int r = 0; r < k; r++) {
                h[j * k + r] = (rng.next() < 0.4) ? rng.next() * 2.0 : 0.0;
            }
        }
        double[] x = new double[m * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                double v = 0.0;
                for (int r = 0; r < k; r++) {
                    v += w[r * m + i] * h[j * k + r];
                }
                v += noise * rng.next();
                x[j * m + i] = (v > 0.0) ? v : 0.0;
            }
        }
        return x;
    }

    /** The in-test generator, so that no run depends on a shared random source. */
    private static final class Lcg {

        private long state;

        Lcg(long seed) {
            this.state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return (state >>> 11) * 0x1.0p-53;
        }
    }
}
