package math.optim;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.DMatrix;
import math.linalg.Lasso;

/**
 * Tests for {@link OrthantWiseLimitedMemoryBFGS}, cross-checked against
 * {@link Lasso}. Both minimize the same objective once the data is prepared so
 * that the lasso's internal standardization is the identity, and the lasso
 * solution is verified directly against the KKT conditions, so it serves as an
 * independent oracle for the L1 path.
 */
public class OwlQnTest {

    private static final int N = 200;
    private static final int P = 12;

    /** Deterministic uniform noise in {@code [-0.5, 0.5]}. */
    static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) * 0x1.0p-53) - 0.5;
        }
    }

    /**
     * The lasso objective in the shape MALLET expects. The optimizer maximizes
     * {@code getValue()} and minimizes {@code -getValue() + l1Weight * sum|w|},
     * so {@code getValue()} is the negated mean squared error and
     * {@code l1Weight} is glmnet's {@code lambda} at {@code alpha == 1}.
     */
    static final class LeastSquares implements Optimizable.ByGradientValue {

        private final double[] x;
        private final double[] y;
        private final int n;
        private final int p;
        private final double[] w;

        LeastSquares(DMatrix X, DMatrix Y) {
            this.n = X.numRows();
            this.p = X.numColumns();
            this.x = X.getArrayUnsafe().clone();
            this.y = Y.getArrayUnsafe().clone();
            this.w = new double[p];
        }

        /** Residuals {@code y - X w}. */
        private double[] residuals() {
            double[] r = new double[n];
            System.arraycopy(y, 0, r, 0, n);
            for (int j = 0; j < p; j++) {
                double wj = w[j];
                if (wj != 0.0) {
                    int off = j * n;
                    for (int i = 0; i < n; i++) {
                        r[i] -= wj * x[off + i];
                    }
                }
            }
            return r;
        }

        @Override
        public int getNumParameters() {
            return p;
        }

        @Override
        public void getParameters(double[] buffer) {
            System.arraycopy(w, 0, buffer, 0, p);
        }

        @Override
        public double getParameter(int index) {
            return w[index];
        }

        @Override
        public void setParameters(double[] params) {
            System.arraycopy(params, 0, w, 0, p);
        }

        @Override
        public void setParameter(int index, double value) {
            w[index] = value;
        }

        @Override
        public double getValue() {
            double[] r = residuals();
            double rss = 0.0;
            for (int i = 0; i < n; i++) {
                rss += r[i] * r[i];
            }
            return -rss / (2.0 * n);
        }

        @Override
        public void getValueGradient(double[] buffer) {
            double[] r = residuals();
            for (int j = 0; j < p; j++) {
                int off = j * n;
                double d = 0.0;
                for (int i = 0; i < n; i++) {
                    d += x[off + i] * r[i];
                }
                buffer[j] = d / n;
            }
        }
    }

    /**
     * The Rosenbrock function, smooth but not convex, so the curvature
     * condition {@code s . y > 0} the BFGS update relies on can be violated.
     * MALLET maximizes {@code getValue()}, hence the negation.
     */
    static final class Rosenbrock implements Optimizable.ByGradientValue {

        private final double[] w;

        Rosenbrock(double[] start) {
            this.w = start.clone();
        }

        @Override
        public int getNumParameters() {
            return w.length;
        }

        @Override
        public void getParameters(double[] buffer) {
            System.arraycopy(w, 0, buffer, 0, w.length);
        }

        @Override
        public double getParameter(int index) {
            return w[index];
        }

        @Override
        public void setParameters(double[] params) {
            System.arraycopy(params, 0, w, 0, w.length);
        }

        @Override
        public void setParameter(int index, double value) {
            w[index] = value;
        }

        @Override
        public double getValue() {
            double s = 0.0;
            for (int i = 0; i + 1 < w.length; i++) {
                double a = w[i + 1] - w[i] * w[i];
                double b = 1.0 - w[i];
                s += 100.0 * a * a + b * b;
            }
            return -s;
        }

        @Override
        public void getValueGradient(double[] buffer) {
            for (int i = 0; i < buffer.length; i++) {
                buffer[i] = 0.0;
            }
            for (int i = 0; i + 1 < w.length; i++) {
                double a = w[i + 1] - w[i] * w[i];
                buffer[i] += 400.0 * w[i] * a + 2.0 * (1.0 - w[i]);
                buffer[i + 1] += -200.0 * a;
            }
        }
    }

    /**
     * A design whose columns are centred and scaled to
     * {@code sum_i x_ij^2 == n}, and a centred response, so that {@link Lasso}
     * re-standardizes to a no-op and fits a zero intercept. Returns
     * {@code {X, y}}.
     */
    static DMatrix[] standardizedProblem(long seed) {
        Lcg rng = new Lcg(seed);
        DMatrix X = new DMatrix(N, P);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < P; j++) {
                X.set(i, j, rng.next());
            }
        }
        double[] beta = new double[P];
        beta[1] = 8.0;
        beta[5] = -6.0;
        beta[9] = 3.0;
        DMatrix y = new DMatrix(N, 1);
        for (int i = 0; i < N; i++) {
            double v = 0.0;
            for (int j = 0; j < P; j++) {
                v += beta[j] * X.get(i, j);
            }
            y.set(i, 0, v + 0.35 * rng.next());
        }

        double[] xs = X.getArrayUnsafe();
        for (int j = 0; j < P; j++) {
            int off = j * N;
            double m = 0.0;
            for (int i = 0; i < N; i++) {
                m += xs[off + i];
            }
            m /= N;
            double ss = 0.0;
            for (int i = 0; i < N; i++) {
                xs[off + i] -= m;
                ss += xs[off + i] * xs[off + i];
            }
            double scale = Math.sqrt(ss / N);
            for (int i = 0; i < N; i++) {
                xs[off + i] /= scale;
            }
        }
        double[] ys = y.getArrayUnsafe();
        double m = 0.0;
        for (int i = 0; i < N; i++) {
            m += ys[i];
        }
        m /= N;
        for (int i = 0; i < N; i++) {
            ys[i] -= m;
        }
        return new DMatrix[] { X, y };
    }

    /** The largest absolute deviation between two coefficient vectors. */
    static double maxDeviation(double[] a, double[] b) {
        double max = 0.0;
        for (int j = 0; j < a.length; j++) {
            max = Math.max(max, Math.abs(a[j] - b[j]));
        }
        return max;
    }

    /** Runs OWL-QN at the given L1 weight and returns the parameters. */
    static double[] fit(LeastSquares f, double l1Weight) {
        OrthantWiseLimitedMemoryBFGS owlqn = new OrthantWiseLimitedMemoryBFGS(f, l1Weight);
        owlqn.optimize();
        double[] w = new double[f.getNumParameters()];
        f.getParameters(w);
        return w;
    }

    /**
     * OWL-QN and the lasso solve the same problem, so their coefficients must
     * agree. The tolerance is what the default stopping rule -- a relative
     * change of 1e-4 in the objective value -- actually delivers.
     */
    @Test
    public void testAgreesWithTheLassoOnTheSameObjective() {
        DMatrix[] problem = standardizedProblem(20260819L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        double[] lambdas = { 0.5, 0.1, 0.01 };
        for (int k = 0; k < lambdas.length; k++) {
            double lambda = lambdas[k];
            double[] expected = Lasso.estimate(X, y, lambda, 1.0).beta;
            double[] actual = fit(new LeastSquares(X, y), lambda);
            double deviation = maxDeviation(expected, actual);
            assertTrue("lambda = " + lambda + ", deviation = " + deviation, deviation < 2.0e-3);
        }
    }

    /**
     * The same cross-check with the stopping rules turned up. The defaults are
     * MALLET's and stop at a relative change of 1e-4 in the objective value;
     * reaching the lasso solution more closely than that was impossible while
     * the tolerances had no way in.
     */
    @Test
    public void testTightenedTolerancesReachTheLassoSolutionFarMoreClosely() {
        DMatrix[] problem = standardizedProblem(20260819L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        double[] lambdas = { 0.5, 0.1, 0.01 };
        for (int k = 0; k < lambdas.length; k++) {
            double lambda = lambdas[k];
            double[] expected = Lasso.estimate(X, y, lambda, 1.0).beta;
            LeastSquares f = new LeastSquares(X, y);
            OrthantWiseLimitedMemoryBFGS owlqn = new OrthantWiseLimitedMemoryBFGS(f, lambda, 1000, 1.0e-10,
                    1.0e-10, 4);
            assertTrue("lambda = " + lambda, owlqn.optimize());
            double[] actual = new double[f.getNumParameters()];
            f.getParameters(actual);
            double deviation = maxDeviation(expected, actual);
            assertTrue("lambda = " + lambda + ", deviation = " + deviation, deviation < 1.0e-5);
        }
    }

    /** The stopping rules are checked, not silently accepted. */
    @Test
    public void testArgumentValidation() {
        DMatrix[] problem = standardizedProblem(3L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        expectIae(null, 0.0, 1000, 1.0e-4, 1.0e-3, 4, X, y);
        expectIae("l1wt", -1.0e-9, 1000, 1.0e-4, 1.0e-3, 4, X, y);
        expectIae("l1wt", Double.NaN, 1000, 1.0e-4, 1.0e-3, 4, X, y);
        expectIae("maxIterations", 0.0, 0, 1.0e-4, 1.0e-3, 4, X, y);
        expectIae("tolerance", 0.0, 1000, 0.0, 1.0e-3, 4, X, y);
        expectIae("gradientTolerance", 0.0, 1000, 1.0e-4, -1.0e-9, 4, X, y);
        expectIae("m", 0.0, 1000, 1.0e-4, 1.0e-3, 0, X, y);
        expectIae("m", 0.0, 1000, 1.0e-4, 1.0e-3, 101, X, y);
    }

    private static void expectIae(String what, double l1wt, int maxIterations, double tolerance,
            double gradientTolerance, int m, DMatrix X, DMatrix y) {
        LeastSquares f = (what == null) ? null : new LeastSquares(X, y);
        try {
            new OrthantWiseLimitedMemoryBFGS(f, l1wt, maxIterations, tolerance, gradientTolerance, m);
            fail("expected IllegalArgumentException for " + what);
        } catch (IllegalArgumentException expected) {
            // as it should be
        }
    }

    /**
     * Without a penalty the problem is an ordinary least squares fit, which the
     * lasso reproduces at {@code lambda == 0}.
     */
    @Test
    public void testWithoutRegularizationReproducesOrdinaryLeastSquares() {
        DMatrix[] problem = standardizedProblem(4711L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        double[] expected = Lasso.estimate(X, y, 0.0, 1.0).beta;
        double[] actual = fit(new LeastSquares(X, y), 0.0);
        double deviation = maxDeviation(expected, actual);
        assertTrue("deviation = " + deviation, deviation < 2.0e-3);
    }

    /**
     * Why the gradient stopping criterion cannot work as written: at an L1
     * solution the gradient of the smooth part equals {@code -l1Weight} times
     * the sign on every active coordinate, so its norm is bounded away from
     * zero. The pseudo-gradient is the quantity that vanishes.
     */
    @Test
    public void testSmoothGradientStaysLargeAtAnL1SolutionWhilePseudoGradientVanishes() {
        DMatrix[] problem = standardizedProblem(99L);
        double lambda = 0.5;
        LeastSquares f = new LeastSquares(problem[0], problem[1]);
        double[] w = fit(f, lambda);

        double[] buffer = new double[w.length];
        f.getValueGradient(buffer);

        double smooth = 0.0;
        double pseudo = 0.0;
        int active = 0;
        for (int j = 0; j < w.length; j++) {
            // the minimized smooth part has the negated getValueGradient
            double g = -buffer[j];
            smooth += g * g;
            double pg;
            if (Math.abs(w[j]) > 1.0e-6) {
                pg = g + lambda * Math.signum(w[j]);
                active++;
            } else {
                pg = Math.max(0.0, Math.abs(g) - lambda);
            }
            pseudo += pg * pg;
        }
        smooth = Math.sqrt(smooth);
        pseudo = Math.sqrt(pseudo);

        assertTrue("expected active coordinates, got " + active, active >= 3);
        assertTrue("smooth gradient norm = " + smooth, smooth > 0.9 * lambda * Math.sqrt(active));
        assertTrue("pseudo gradient norm = " + pseudo, pseudo < 1.0e-2);
    }

    /** Increasing the penalty must not increase the L1 norm of the solution. */
    @Test
    public void testL1NormNeverGrowsWithThePenalty() {
        DMatrix[] problem = standardizedProblem(13L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        double[] lambdas = { 0.01, 0.05, 0.1, 0.25, 0.5, 1.0 };
        double previous = Double.MAX_VALUE;
        for (int k = 0; k < lambdas.length; k++) {
            double[] w = fit(new LeastSquares(X, y), lambdas[k]);
            double l1 = 0.0;
            for (int j = 0; j < w.length; j++) {
                l1 += Math.abs(w[j]);
            }
            assertTrue("lambda = " + lambdas[k] + ", |w|_1 = " + l1 + " > " + previous, l1 <= previous + 1.0e-6);
            previous = l1;
        }
    }

    /**
     * A penalty above {@code lambdaMax} makes the zero vector the solution, and
     * the zero vector is where the search starts. The optimizer must recognise
     * that the optimality conditions already hold instead of treating the
     * vanishing direction as a caller error.
     */
    @Test
    public void testPenaltyAboveLambdaMaxConvergesAtTheStartingPoint() {
        DMatrix[] problem = standardizedProblem(13L);
        DMatrix X = problem[0];
        DMatrix y = problem[1];
        double lambda = 1.5 * Lasso.lambdaMax(X, y, 1.0);
        LeastSquares f = new LeastSquares(X, y);
        OrthantWiseLimitedMemoryBFGS owlqn = new OrthantWiseLimitedMemoryBFGS(f, lambda);

        assertTrue("optimize() must report convergence", owlqn.optimize());
        assertTrue("isConverged() must agree", owlqn.isConverged());

        double[] w = new double[f.getNumParameters()];
        f.getParameters(w);
        for (int j = 0; j < w.length; j++) {
            assertEquals("coefficient " + j, 0.0, w[j], 0.0);
        }
    }

    /**
     * A non-convex objective can violate the curvature condition the BFGS
     * update needs. Skipping the update keeps the search going on the pairs
     * already stored; treating it as a caller error stops it dead. The classic
     * Rosenbrock start is such a case in its very first iterations.
     */
    @Test
    public void testNonConvexObjectiveSurvivesAViolatedCurvatureCondition() {
        Rosenbrock f = new Rosenbrock(new double[] { -1.2, 1.0 });
        OrthantWiseLimitedMemoryBFGS owlqn = new OrthantWiseLimitedMemoryBFGS(f, 0.0);

        assertTrue("optimize() must report convergence", owlqn.optimize());

        double[] w = new double[2];
        f.getParameters(w);
        assertEquals("w[0]", 1.0, w[0], 1.0e-4);
        assertEquals("w[1]", 1.0, w[1], 1.0e-4);
        assertEquals("minimum value", 0.0, -f.getValue(), 1.0e-8);
    }

    /** The optimizable is left on the parameters the optimizer settled on. */
    @Test
    public void testOptimizableEndsOnTheReportedPoint() {
        DMatrix[] problem = standardizedProblem(7L);
        LeastSquares f = new LeastSquares(problem[0], problem[1]);
        double[] w = fit(f, 0.1);
        double[] inside = new double[w.length];
        f.getParameters(inside);
        assertEquals(0.0, maxDeviation(w, inside), 0.0);
    }
}
