package math.minpack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.linalg.VectorOps;

/**
 * Tests for the MINPACK port, which had none at all: 3590 lines, no test and no
 * caller. The oracle is external -- the problem set of More, Garbow and
 * Hillstrom, <i>Testing Unconstrained Optimization Software</i>, ACM TOMS 7(1),
 * 1981 -- with the starting points and minimal sums of squares as published
 * there.
 * <p>
 * The published numbers are not the load-bearing assertion, though. Two checks
 * come first and depend on nothing but the problem definitions in this file:
 * that the analytic Jacobians here agree with finite differences, so that a
 * wrong derivation cannot be mistaken for a defect in the port, and that the
 * gradient of the sum of squares vanishes at every reported solution, which is
 * the defining equation of a least squares minimum. The published minima then
 * catch what those two cannot -- a mistyped data constant, which is
 * self-consistent and would pass both.
 */
public class MinpackTest {

    /** Recommended lower limit for the MINPACK tolerances. */
    private static final double SQRT_EPS = Math.sqrt(Math.ulp(1.0));

    private static final double[] BARD_Y = { 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96,
            1.34, 2.10, 4.39 };

    private static final double[] MEYER_Y = { 34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0,
            8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0 };

    private static final double[] KOWALIK_Y = { 0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323,
            0.0235, 0.0246 };

    private static final double[] KOWALIK_U = { 4.0, 2.0, 1.0, 0.5, 0.25, 0.167, 0.125, 0.1, 0.0833, 0.0714, 0.0625 };

    private static final double[] OSBORNE1_Y = { 0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818, 0.784,
            0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558, 0.538, 0.522, 0.506, 0.490, 0.478, 0.467, 0.457,
            0.448, 0.438, 0.431, 0.424, 0.420, 0.414, 0.411, 0.406 };

    /**
     * One problem of the collection, stated zero-based and independent of the
     * FORTRAN calling convention the port uses.
     */
    private abstract static class Problem {

        final String name;
        final int m;
        final int n;
        final double[] start;
        /** Published sum of squares at the minimum. */
        final double minimum;
        /** Tolerance on it, set by how many digits the source publishes. */
        final double minimumTolerance;
        /** Published minimizer, or {@code null} when only the value is given. */
        final double[] solution;
        /** Tolerance on the minimizer, relative for components above one. */
        final double solutionTolerance;
        final int maxEvaluations;

        Problem(String name, int m, int n, double[] start, double minimum, double minimumTolerance, double[] solution,
                double solutionTolerance, int maxEvaluations) {
            this.name = name;
            this.m = m;
            this.n = n;
            this.start = start;
            this.minimum = minimum;
            this.minimumTolerance = minimumTolerance;
            this.solution = solution;
            this.solutionTolerance = solutionTolerance;
            this.maxEvaluations = maxEvaluations;
        }

        abstract void residuals(double[] x, double[] r);

        abstract void jacobian(double[] x, double[][] jac);

        boolean hasZeroResidual() {
            return minimum == 0.0;
        }
    }

    private static Problem[] problems() {
        return new Problem[] {

                // (1) Rosenbrock
                new Problem("Rosenbrock", 2, 2, new double[] { -1.2, 1.0 }, 0.0, 1.0e-14, new double[] { 1.0, 1.0 },
                        1.0e-8, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        r[0] = 10.0 * (x[1] - x[0] * x[0]);
                        r[1] = 1.0 - x[0];
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        jac[0][0] = -20.0 * x[0];
                        jac[0][1] = 10.0;
                        jac[1][0] = -1.0;
                        jac[1][1] = 0.0;
                    }
                },

                // (3) Powell badly scaled
                new Problem("Powell badly scaled", 2, 2, new double[] { 0.0, 1.0 }, 0.0, 1.0e-14,
                        new double[] { 1.098159e-5, 9.106146 }, 1.0e-5, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        r[0] = 1.0e4 * x[0] * x[1] - 1.0;
                        r[1] = Math.exp(-x[0]) + Math.exp(-x[1]) - 1.0001;
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        jac[0][0] = 1.0e4 * x[1];
                        jac[0][1] = 1.0e4 * x[0];
                        jac[1][0] = -Math.exp(-x[0]);
                        jac[1][1] = -Math.exp(-x[1]);
                    }
                },

                // (4) Brown badly scaled
                new Problem("Brown badly scaled", 3, 2, new double[] { 1.0, 1.0 }, 0.0, 1.0e-8,
                        new double[] { 1.0e6, 2.0e-6 }, 1.0e-6, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        r[0] = x[0] - 1.0e6;
                        r[1] = x[1] - 2.0e-6;
                        r[2] = x[0] * x[1] - 2.0;
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        jac[0][0] = 1.0;
                        jac[0][1] = 0.0;
                        jac[1][0] = 0.0;
                        jac[1][1] = 1.0;
                        jac[2][0] = x[1];
                        jac[2][1] = x[0];
                    }
                },

                // (5) Beale
                new Problem("Beale", 3, 2, new double[] { 1.0, 1.0 }, 0.0, 1.0e-14, new double[] { 3.0, 0.5 }, 1.0e-7,
                        1000) {
                    private final double[] y = { 1.5, 2.25, 2.625 };

                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 3; i++) {
                            r[i] = y[i] - x[0] * (1.0 - Math.pow(x[1], i + 1));
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 3; i++) {
                            jac[i][0] = -(1.0 - Math.pow(x[1], i + 1));
                            jac[i][1] = (i + 1) * x[0] * Math.pow(x[1], i);
                        }
                    }
                },

                // (7) Helical valley
                new Problem("Helical valley", 3, 3, new double[] { -1.0, 0.0, 0.0 }, 0.0, 1.0e-14,
                        new double[] { 1.0, 0.0, 0.0 }, 1.0e-8, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        // theta as the source states it, not atan2: for x0 < 0
                        // the two differ by a full turn on one side of the
                        // negative axis, which is where this problem starts.
                        // atan2 jumps there, the definition below does not.
                        double theta = Math.atan(x[1] / x[0]) / (2.0 * Math.PI);
                        if (x[0] < 0.0) {
                            theta += 0.5;
                        }
                        r[0] = 10.0 * (x[2] - 10.0 * theta);
                        r[1] = 10.0 * (Math.sqrt(x[0] * x[0] + x[1] * x[1]) - 1.0);
                        r[2] = x[2];
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        double tmp = x[0] * x[0] + x[1] * x[1];
                        double root = Math.sqrt(tmp);
                        double tpi = 2.0 * Math.PI;
                        jac[0][0] = 100.0 * x[1] / (tpi * tmp);
                        jac[0][1] = -100.0 * x[0] / (tpi * tmp);
                        jac[0][2] = 10.0;
                        jac[1][0] = 10.0 * x[0] / root;
                        jac[1][1] = 10.0 * x[1] / root;
                        jac[1][2] = 0.0;
                        jac[2][0] = 0.0;
                        jac[2][1] = 0.0;
                        jac[2][2] = 1.0;
                    }
                },

                // (8) Bard
                new Problem("Bard", 15, 3, new double[] { 1.0, 1.0, 1.0 }, 8.21487e-3, 1.0e-8, null, 0.0, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 15; i++) {
                            double u = i + 1;
                            double v = 15 - i;
                            double w = Math.min(u, v);
                            r[i] = BARD_Y[i] - (x[0] + u / (v * x[1] + w * x[2]));
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 15; i++) {
                            double u = i + 1;
                            double v = 15 - i;
                            double w = Math.min(u, v);
                            double d = v * x[1] + w * x[2];
                            jac[i][0] = -1.0;
                            jac[i][1] = u * v / (d * d);
                            jac[i][2] = u * w / (d * d);
                        }
                    }
                },

                // (10) Meyer
                new Problem("Meyer", 16, 3, new double[] { 0.02, 4000.0, 250.0 }, 87.9458, 1.0e-3, null, 0.0, 4000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 16; i++) {
                            double t = 45.0 + 5.0 * (i + 1);
                            r[i] = x[0] * Math.exp(x[1] / (t + x[2])) - MEYER_Y[i];
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 16; i++) {
                            double t = 45.0 + 5.0 * (i + 1);
                            double d = t + x[2];
                            double e = Math.exp(x[1] / d);
                            jac[i][0] = e;
                            jac[i][1] = x[0] * e / d;
                            jac[i][2] = -x[0] * e * x[1] / (d * d);
                        }
                    }
                },

                // (13) Powell singular -- the Jacobian is rank deficient at the
                // minimum, which is what makes it worth having here
                new Problem("Powell singular", 4, 4, new double[] { 3.0, -1.0, 0.0, 1.0 }, 0.0, 1.0e-12,
                        new double[] { 0.0, 0.0, 0.0, 0.0 }, 1.0e-3, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        r[0] = x[0] + 10.0 * x[1];
                        r[1] = Math.sqrt(5.0) * (x[2] - x[3]);
                        r[2] = (x[1] - 2.0 * x[2]) * (x[1] - 2.0 * x[2]);
                        r[3] = Math.sqrt(10.0) * (x[0] - x[3]) * (x[0] - x[3]);
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 4; i++) {
                            for (int j = 0; j < 4; j++) {
                                jac[i][j] = 0.0;
                            }
                        }
                        jac[0][0] = 1.0;
                        jac[0][1] = 10.0;
                        jac[1][2] = Math.sqrt(5.0);
                        jac[1][3] = -Math.sqrt(5.0);
                        jac[2][1] = 2.0 * (x[1] - 2.0 * x[2]);
                        jac[2][2] = -4.0 * (x[1] - 2.0 * x[2]);
                        jac[3][0] = 2.0 * Math.sqrt(10.0) * (x[0] - x[3]);
                        jac[3][3] = -2.0 * Math.sqrt(10.0) * (x[0] - x[3]);
                    }
                },

                // (15) Kowalik and Osborne
                new Problem("Kowalik and Osborne", 11, 4, new double[] { 0.25, 0.39, 0.415, 0.39 }, 3.07505e-4, 1.0e-9,
                        null, 0.0, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 11; i++) {
                            double u = KOWALIK_U[i];
                            double num = u * u + u * x[1];
                            double den = u * u + u * x[2] + x[3];
                            r[i] = KOWALIK_Y[i] - x[0] * num / den;
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 11; i++) {
                            double u = KOWALIK_U[i];
                            double num = u * u + u * x[1];
                            double den = u * u + u * x[2] + x[3];
                            jac[i][0] = -num / den;
                            jac[i][1] = -x[0] * u / den;
                            jac[i][2] = x[0] * num * u / (den * den);
                            jac[i][3] = x[0] * num / (den * den);
                        }
                    }
                },

                // (16) Brown and Dennis
                new Problem("Brown and Dennis", 20, 4, new double[] { 25.0, 5.0, -5.0, -1.0 }, 85822.2, 0.1, null, 0.0,
                        4000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 20; i++) {
                            double t = (i + 1) / 5.0;
                            double a = x[0] + t * x[1] - Math.exp(t);
                            double b = x[2] + x[3] * Math.sin(t) - Math.cos(t);
                            r[i] = a * a + b * b;
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 20; i++) {
                            double t = (i + 1) / 5.0;
                            double a = x[0] + t * x[1] - Math.exp(t);
                            double b = x[2] + x[3] * Math.sin(t) - Math.cos(t);
                            jac[i][0] = 2.0 * a;
                            jac[i][1] = 2.0 * a * t;
                            jac[i][2] = 2.0 * b;
                            jac[i][3] = 2.0 * b * Math.sin(t);
                        }
                    }
                },

                // (17) Osborne 1
                new Problem("Osborne 1", 33, 5, new double[] { 0.5, 1.5, -1.0, 0.01, 0.02 }, 5.46489e-5, 1.0e-9, null,
                        0.0, 1000) {
                    @Override
                    void residuals(double[] x, double[] r) {
                        for (int i = 0; i < 33; i++) {
                            double t = 10.0 * i;
                            r[i] = OSBORNE1_Y[i]
                                    - (x[0] + x[1] * Math.exp(-t * x[3]) + x[2] * Math.exp(-t * x[4]));
                        }
                    }

                    @Override
                    void jacobian(double[] x, double[][] jac) {
                        for (int i = 0; i < 33; i++) {
                            double t = 10.0 * i;
                            jac[i][0] = -1.0;
                            jac[i][1] = -Math.exp(-t * x[3]);
                            jac[i][2] = -Math.exp(-t * x[4]);
                            jac[i][3] = x[1] * t * Math.exp(-t * x[3]);
                            jac[i][4] = x[2] * t * Math.exp(-t * x[4]);
                        }
                    }
                } };
    }

    /** What a run of the port returns, unpacked back to zero-based arrays. */
    private static final class Outcome {

        double[] x;
        double[] residuals;
        double sumOfSquares;
        int info;
        int functionEvaluations;
        int jacobianEvaluations;
    }

    /** Bridges a {@link Problem} to the one-based FORTRAN calling convention. */
    private static final class Bridge implements Lmder_fcn, Lmdif_fcn {

        private final Problem p;
        private final double[] x;
        private final double[] r;
        private final double[][] jac;

        Bridge(Problem p) {
            this.p = p;
            this.x = new double[p.n];
            this.r = new double[p.m];
            this.jac = new double[p.m][p.n];
        }

        private void unpack(double[] xOneBased) {
            for (int j = 0; j < p.n; j++) {
                x[j] = xOneBased[j + 1];
            }
        }

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, double[][] fjac, int[] iflag) {
            unpack(xOneBased);
            if (iflag[1] == 1) {
                p.residuals(x, r);
                for (int i = 0; i < m; i++) {
                    fvec[i + 1] = r[i];
                }
            } else {
                p.jacobian(x, jac);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        fjac[i + 1][j + 1] = jac[i][j];
                    }
                }
            }
        }

        @Override
        public void fcn(int m, int n, double[] xOneBased, double[] fvec, int[] iflag) {
            unpack(xOneBased);
            p.residuals(x, r);
            for (int i = 0; i < m; i++) {
                fvec[i + 1] = r[i];
            }
        }
    }

    /** Runs {@code lmder_f77}, the full driver -- {@code lmder1_f77} keeps its
     * evaluation counts to itself. */
    private static Outcome runLmder(Problem p, double tol) {
        int m = p.m;
        int n = p.n;
        double[] x = new double[n + 1];
        for (int j = 0; j < n; j++) {
            x[j + 1] = p.start[j];
        }
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        Minpack_f77.lmder_f77(new Bridge(p), m, n, x, fvec, fjac, tol, tol, 0.0, p.maxEvaluations, diag, 1, 100.0, 0,
                info, nfev, njev, ipvt, qtf);

        return unpack(p, x, fvec, info, nfev, njev);
    }

    /** Runs {@code lmdif_f77}, which builds the Jacobian by forward differences. */
    private static Outcome runLmdif(Problem p, double tol) {
        int m = p.m;
        int n = p.n;
        double[] x = new double[n + 1];
        for (int j = 0; j < n; j++) {
            x[j + 1] = p.start[j];
        }
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];

        Minpack_f77.lmdif_f77(new Bridge(p), m, n, x, fvec, tol, tol, 0.0, p.maxEvaluations, 0.0, diag, 1, 100.0, 0,
                info, nfev, fjac, ipvt, qtf);

        return unpack(p, x, fvec, info, nfev, new int[2]);
    }

    private static Outcome unpack(Problem p, double[] x, double[] fvec, int[] info, int[] nfev, int[] njev) {
        Outcome o = new Outcome();
        o.x = new double[p.n];
        for (int j = 0; j < p.n; j++) {
            o.x[j] = x[j + 1];
        }
        o.residuals = new double[p.m];
        double ssq = 0.0;
        for (int i = 0; i < p.m; i++) {
            o.residuals[i] = fvec[i + 1];
            ssq += fvec[i + 1] * fvec[i + 1];
        }
        o.sumOfSquares = ssq;
        o.info = info[1];
        o.functionEvaluations = nfev[1];
        o.jacobianEvaluations = njev[1];
        return o;
    }

    /**
     * The success codes. One and two are the two tolerances met, three is both,
     * and four is a residual orthogonal to the Jacobian -- a stationary point.
     * <p>
     * Eight belongs here too, and only looks like a failure: it says the same
     * thing four does and adds that the requested {@code gtol} was below what
     * can be achieved, which is trivially the case because these runs ask for
     * {@code gtol == 0}. {@code lmder1_f77} makes the same judgement in the
     * port itself, remapping eight to four before it returns. The distinction
     * matters for the facade -- read as a plain number, eight is larger than
     * every success code and smaller than nothing.
     */
    private static boolean isSuccess(int info) {
        return (info >= 1 && info <= 4) || info == 8;
    }

    private static double sumOfSquares(Problem p, double[] x) {
        double[] r = new double[p.m];
        p.residuals(x, r);
        double s = 0.0;
        for (int i = 0; i < p.m; i++) {
            s += r[i] * r[i];
        }
        return s;
    }

    /**
     * Every analytic Jacobian in this file, against central differences. This
     * runs first on purpose: a wrong derivation here would send the solver
     * somewhere odd and look exactly like a defect in the port.
     */
    @Test
    public void testAnalyticJacobiansAgreeWithFiniteDifferences() {
        double h = Math.cbrt(Math.ulp(1.0));
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            double[][] analytic = new double[p.m][p.n];
            p.jacobian(p.start, analytic);

            double[] plus = new double[p.m];
            double[] minus = new double[p.m];
            double[] x = p.start.clone();
            for (int j = 0; j < p.n; j++) {
                double step = h * Math.max(Math.abs(p.start[j]), 1.0);
                x[j] = p.start[j] + step;
                p.residuals(x, plus);
                x[j] = p.start[j] - step;
                p.residuals(x, minus);
                x[j] = p.start[j];

                for (int i = 0; i < p.m; i++) {
                    double numeric = (plus[i] - minus[i]) / (2.0 * step);
                    double scale = Math.max(Math.abs(numeric), 1.0);
                    assertEquals(p.name + ": d f[" + i + "] / d x[" + j + "]", numeric, analytic[i][j],
                            1.0e-5 * scale);
                }
            }
        }
    }

    /**
     * The problems whose residuals vanish at the solution are the ones where
     * both the minimizer and the minimum are known exactly.
     */
    @Test
    public void testZeroResidualProblemsReachTheirKnownSolution() {
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            if (!p.hasZeroResidual()) {
                continue;
            }
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": info = " + o.info, isSuccess(o.info));
            assertEquals(p.name + ": sum of squares", 0.0, o.sumOfSquares, p.minimumTolerance);
            for (int j = 0; j < p.n; j++) {
                double scale = Math.max(Math.abs(p.solution[j]), 1.0);
                assertEquals(p.name + ": x[" + j + "]", p.solution[j], o.x[j], p.solutionTolerance * scale);
            }
        }
    }

    /**
     * The rest are data fitting problems with a residual that does not vanish.
     * Only the minimal sum of squares is published for them, to six digits, and
     * the tolerance here is that of the published number rather than of the
     * algorithm.
     */
    @Test
    public void testNonZeroResidualProblemsReachThePublishedMinimum() {
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            if (p.hasZeroResidual()) {
                continue;
            }
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": info = " + o.info, isSuccess(o.info));
            assertEquals(p.name + ": sum of squares", p.minimum, o.sumOfSquares, p.minimumTolerance);
        }
    }

    /**
     * The check that needs no published number: at a least squares minimum the
     * residual is orthogonal to every column of the Jacobian, which is the
     * defining equation. Stated as a cosine so that it is scale free and one
     * bound fits problems whose residuals differ by ten orders of magnitude.
     * <p>
     * Only the problems with a residual that does not vanish are asked. Where
     * it does, the cosine is the angle of a vector of length 1e-16 and carries
     * no information; for those the sum of squares itself is the stronger
     * statement, and the test above makes it.
     * <p>
     * The bound has to be stated relative to the tolerance that was asked for,
     * because MINPACK stops on the value and the step, never on the gradient
     * unless {@code gtol} says so, and these runs pass {@code gtol == 0}. The
     * gradient is therefore only as small as the other two force it to be, and
     * the rate differs by regime:
     * <ul>
     * <li>Kowalik and Osborne has a minimal sum of squares of 3e-4. The Gauss
     * Newton model is good there, convergence is quadratic, and the cosine
     * falls <em>linearly</em> with the tolerance -- 5.8e-6 at 1.5e-8, 5.5e-8 at
     * 1e-12, a factor of about 400 throughout.</li>
     * <li>Brown and Dennis has a minimal sum of squares of 85822. Only the
     * value test binds there, and the value is quadratic at the minimum, so a
     * relative tolerance of {@code tol} in it pins the point to
     * {@code sqrt(tol)} and the gradient with it -- 1.6e-4 at 1.5e-8, 9.3e-7 at
     * 1e-12, both within a factor of two of {@code sqrt(tol)}.</li>
     * </ul>
     * The second is the weaker of the two and sets the bound. Running at two
     * tolerances three orders of magnitude apart asserts the coupling rather
     * than a single number, which is the part that would be hard to reproduce
     * by accident.
     */
    @Test
    public void testTheGradientVanishesAtEverySolution() {
        double[] tolerances = { SQRT_EPS, 1.0e-12 };
        Problem[] all = problems();
        for (int t = 0; t < tolerances.length; t++) {
            double tol = tolerances[t];
            for (int k = 0; k < all.length; k++) {
                Problem p = all[k];
                if (p.hasZeroResidual()) {
                    continue;
                }
                Outcome o = runLmder(p, tol);

                double[][] jac = new double[p.m][p.n];
                p.jacobian(o.x, jac);
                double residualNorm = VectorOps.twoNorm(o.residuals);
                for (int j = 0; j < p.n; j++) {
                    double dot = 0.0;
                    double columnNorm = 0.0;
                    for (int i = 0; i < p.m; i++) {
                        dot += jac[i][j] * o.residuals[i];
                        columnNorm += jac[i][j] * jac[i][j];
                    }
                    columnNorm = Math.sqrt(columnNorm);
                    if (columnNorm == 0.0) {
                        continue;
                    }
                    double cosine = Math.abs(dot) / (columnNorm * residualNorm);
                    assertTrue(p.name + " at tol " + tol + ": column " + j
                            + " is not orthogonal to the residual, cosine = " + cosine, cosine < 10.0 * Math.sqrt(tol));
                }
            }
        }
    }

    /**
     * {@code fvec} comes back as an output argument, so it can disagree with
     * the returned point without anything complaining. It must not.
     */
    @Test
    public void testReturnedResidualsBelongToTheReturnedPoint() {
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome o = runLmder(p, SQRT_EPS);

            double[] expected = new double[p.m];
            p.residuals(o.x, expected);
            for (int i = 0; i < p.m; i++) {
                assertEquals(p.name + ": residual " + i, expected[i], o.residuals[i], 0.0);
            }
        }
    }

    /**
     * The solver never returns a point worse than the one it was given. Trivial
     * to state, and the kind of thing a mistranslated index breaks first.
     */
    @Test
    public void testNoProblemEndsWorseThanItStarted() {
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome o = runLmder(p, SQRT_EPS);

            assertTrue(p.name + ": " + o.sumOfSquares + " > " + sumOfSquares(p, p.start),
                    o.sumOfSquares <= sumOfSquares(p, p.start));
        }
    }

    /**
     * {@code lmdif} builds the Jacobian by forward differences instead of
     * asking for it. That is a different algorithm, so what is asserted here is
     * the minimum reached and the price paid for it, never the path or the
     * status code -- see the test below for the problem where those two part
     * company. The budgets each problem carries are sized for this path rather
     * than for the analytic one, which is the cheaper of the two throughout:
     * Brown and Dennis, for instance, converges in 254 evaluations with a
     * Jacobian and 1229 without.
     */
    @Test
    public void testLmdifReachesTheSameMinimaWithoutAnAnalyticJacobian() {
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            Problem p = all[k];
            Outcome withJacobian = runLmder(p, SQRT_EPS);
            Outcome withoutJacobian = runLmdif(p, SQRT_EPS);

            double scale = Math.max(Math.abs(withJacobian.sumOfSquares), 1.0e-8);
            assertEquals(p.name + ": sum of squares", withJacobian.sumOfSquares, withoutJacobian.sumOfSquares,
                    1.0e-4 * scale);
            assertTrue(p.name + ": not asking for the Jacobian cannot be cheaper, "
                    + withoutJacobian.functionEvaluations + " against " + withJacobian.functionEvaluations,
                    withoutJacobian.functionEvaluations >= withJacobian.functionEvaluations);
            assertTrue(withJacobian.jacobianEvaluations >= 0);
        }
    }

    /**
     * Powell singular has a rank deficient Jacobian at its minimum, and it is
     * the one problem here that separates <em>having</em> the answer from
     * <em>being able to say so</em>. With an analytic Jacobian the orthogonality
     * test fires and the run ends at once. Without one, the noise of the
     * forward difference keeps that test from firing, and MINPACK falls back on
     * the relative step test -- which converges linearly on a singular problem
     * and needs about 98000 evaluations to be satisfied, having reached a sum
     * of squares of 2e-41 within the first thousand.
     * <p>
     * So {@code info = 5} here reports an exhausted budget on a point that is
     * exact to eleven digits. Any facade over this port that maps five to
     * "failed" without also handing back the point would be throwing away a
     * correct answer.
     */
    @Test
    public void testTheDerivativeFreePathPaysForARankDeficientJacobian() {
        Problem powellSingular = null;
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            if ("Powell singular".equals(all[k].name)) {
                powellSingular = all[k];
            }
        }
        assertTrue("Powell singular is in the set", powellSingular != null);

        Outcome withJacobian = runLmder(powellSingular, SQRT_EPS);
        Outcome withoutJacobian = runLmdif(powellSingular, SQRT_EPS);

        assertEquals("with an analytic Jacobian the orthogonality test fires", 8, withJacobian.info);
        assertTrue("and it is cheap: " + withJacobian.functionEvaluations, withJacobian.functionEvaluations < 100);

        assertEquals("without one the budget runs out", 5, withoutJacobian.info);
        assertTrue("yet the point is already exact, sum of squares = " + withoutJacobian.sumOfSquares,
                withoutJacobian.sumOfSquares < 1.0e-30);
        for (int j = 0; j < powellSingular.n; j++) {
            assertEquals("x[" + j + "]", 0.0, withoutJacobian.x[j], 1.0e-8);
        }
    }

    /**
     * {@code enorm_f77} carries its own three-accumulator scaling, and the
     * expected values here are stated in closed form rather than taken from
     * another implementation, because the obvious candidate for one is wrong:
     * {@link VectorOps#twoNorm} accumulates the squares directly, so it
     * overflows to infinity above about {@code 1e154} and underflows to zero
     * below about {@code 1e-162}. On the vectors that stay in range the two
     * agree, and on the ones that do not, {@code enorm_f77} is the one telling
     * the truth. Asserted in both directions so the day {@code twoNorm} is
     * rewritten in the BLAS {@code dnrm2} style this test notices.
     */
    @Test
    public void testEnormIsCorrectWhereTheNaiveNormOverflows() {
        double root3 = Math.sqrt(3.0);
        double[][] vectors = { { 3.0, 4.0 }, { 1.0e300, 1.0e300, 1.0e300 }, { 1.0e-300, 1.0e-300, 1.0e-300 },
                { 0.0, 0.0, 0.0 }, { 1.0e300, 1.0e-300, 1.0 }, { -3.0, 0.0, 4.0 } };
        double[] expected = { 5.0, root3 * 1.0e300, root3 * 1.0e-300, 0.0, 1.0e300, 5.0 };
        boolean[] withinNaiveRange = { true, false, false, true, false, true };

        for (int k = 0; k < vectors.length; k++) {
            double[] v = vectors[k];
            double[] oneBased = new double[v.length + 1];
            System.arraycopy(v, 0, oneBased, 1, v.length);

            double enorm = Minpack_f77.enorm_f77(v.length, oneBased);
            assertEquals("vector " + k, expected[k], enorm, 1.0e-14 * Math.max(expected[k], Double.MIN_NORMAL));

            double naive = VectorOps.twoNorm(v);
            if (withinNaiveRange[k]) {
                assertEquals("vector " + k + ": both agree in range", enorm, naive, 1.0e-14 * Math.max(enorm, 1.0));
            } else {
                assertTrue("vector " + k + ": twoNorm is expected to fail here, and returned " + naive,
                        naive == 0.0 || Double.isInfinite(naive));
            }
        }
    }

    /**
     * The trap this port sets for a caller who does not read the FORTRAN
     * comments: {@code info = 5} means the evaluation budget ran out. It is a
     * larger number than the four success codes and looks like a better one.
     */
    @Test
    public void testBudgetExhaustionIsReportedAsFiveRatherThanSuccess() {
        Problem meyer = null;
        Problem[] all = problems();
        for (int k = 0; k < all.length; k++) {
            if ("Meyer".equals(all[k].name)) {
                meyer = all[k];
            }
        }
        assertTrue("Meyer is in the set", meyer != null);

        int m = meyer.m;
        int n = meyer.n;
        double[] x = new double[n + 1];
        for (int j = 0; j < n; j++) {
            x[j + 1] = meyer.start[j];
        }
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] info = new int[2];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        Minpack_f77.lmder_f77(new Bridge(meyer), m, n, x, fvec, fjac, SQRT_EPS, SQRT_EPS, 0.0, 3, diag, 1, 100.0, 0,
                info, nfev, njev, ipvt, qtf);

        assertEquals("info", 5, info[1]);
        assertFalse("5 is not one of the success codes", info[1] >= 1 && info[1] <= 4);
    }

    /** Improper input is reported as {@code info = 0}, not thrown. */
    @Test
    public void testImproperInputIsReportedAsZero() {
        Problem p = problems()[0];
        int m = p.m;
        int n = p.n;
        double[] fvec = new double[m + 1];
        double[][] fjac = new double[m + 1][n + 1];
        double[] diag = new double[n + 1];
        double[] qtf = new double[n + 1];
        int[] ipvt = new int[n + 1];
        int[] nfev = new int[2];
        int[] njev = new int[2];

        // a negative tolerance
        int[] info = new int[2];
        Minpack_f77.lmder_f77(new Bridge(p), m, n, new double[n + 1], fvec, fjac, -1.0, SQRT_EPS, 0.0, 100, diag, 1,
                100.0, 0, info, nfev, njev, ipvt, qtf);
        assertEquals("negative ftol", 0, info[1]);

        // more variables than functions
        info = new int[2];
        Minpack_f77.lmder_f77(new Bridge(p), 1, n, new double[n + 1], fvec, fjac, SQRT_EPS, SQRT_EPS, 0.0, 100, diag,
                1, 100.0, 0, info, nfev, njev, ipvt, qtf);
        assertEquals("m < n", 0, info[1]);
    }
}
