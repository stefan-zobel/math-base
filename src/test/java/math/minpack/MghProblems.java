package math.minpack;

import math.fun.DiffDVectorFunction;

/**
 * The nonlinear least squares problems of More, Garbow and Hillstrom,
 * <i>Testing Unconstrained Optimization Software</i>, ACM TOMS 7(1), 1981, with
 * the starting points and minimal sums of squares as published there. The
 * numbers in the comments are the numbering of that paper.
 * <p>
 * Eleven of the thirty five are here, chosen so that each way of being
 * difficult occurs once: a residual that vanishes at the minimum and one that
 * does not, conditioning spread over twelve orders of magnitude, a Jacobian
 * that is rank deficient where the minimum is, and a problem large enough for
 * the scaling to matter.
 * <p>
 * They are stated as {@link DiffDVectorFunction}, which is the shape the
 * library uses, so the same definitions serve both the FORTRAN calling
 * convention of {@link Minpack_f77} and the facade over it. Two doors, one set
 * of problems, and the results have to agree exactly.
 */
public abstract class MghProblems {

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

    private MghProblems() {
    }

    /** One problem of the collection, together with what is known about it. */
    public abstract static class Problem implements DiffDVectorFunction {

        /** Its name in the source. */
        public final String name;
        /** Number of residuals. */
        public final int m;
        /** Number of parameters. */
        public final int n;
        /** The published starting point. */
        public final double[] start;
        /** Published sum of squares at the minimum. */
        public final double minimum;
        /** Tolerance on it, set by how many digits the source publishes. */
        public final double minimumTolerance;
        /** Published minimizer, or {@code null} when only the value is given. */
        public final double[] solution;
        /** Tolerance on the minimizer, relative for components above one. */
        public final double solutionTolerance;
        /** A budget that suffices for the derivative-free path, the costlier one. */
        public final int maxEvaluations;

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

        /** Whether the residual vanishes at the minimum. */
        public final boolean hasZeroResidual() {
            return minimum == 0.0;
        }

        /**
         * The sum of the squared residuals at {@code x}.
         *
         * @param x
         *            the parameters
         * @return the sum of squares there
         */
        public final double sumOfSquaresAt(double[] x) {
            double[] r = new double[m];
            valueAt(x, r);
            double s = 0.0;
            for (int i = 0; i < m; i++) {
                s += r[i] * r[i];
            }
            return s;
        }

        @Override
        public final String toString() {
            return name;
        }
    }

    /**
     * A fresh copy of the collection. Fresh because a caller may want to count
     * evaluations, and a shared instance would carry the previous count.
     *
     * @return the eleven problems
     */
    public static Problem[] all() {
        return new Problem[] {

                // (1) Rosenbrock
                new Problem("Rosenbrock", 2, 2, new double[] { -1.2, 1.0 }, 0.0, 1.0e-14, new double[] { 1.0, 1.0 },
                        1.0e-8, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        r[0] = 10.0 * (x[1] - x[0] * x[0]);
                        r[1] = 1.0 - x[0];
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        jac[0 * 2 + 0] = -20.0 * x[0];
                        jac[1 * 2 + 0] = 10.0;
                        jac[0 * 2 + 1] = -1.0;
                        jac[1 * 2 + 1] = 0.0;
                    }
                },

                // (3) Powell badly scaled
                new Problem("Powell badly scaled", 2, 2, new double[] { 0.0, 1.0 }, 0.0, 1.0e-14,
                        new double[] { 1.098159e-5, 9.106146 }, 1.0e-5, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        r[0] = 1.0e4 * x[0] * x[1] - 1.0;
                        r[1] = Math.exp(-x[0]) + Math.exp(-x[1]) - 1.0001;
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        jac[0 * 2 + 0] = 1.0e4 * x[1];
                        jac[1 * 2 + 0] = 1.0e4 * x[0];
                        jac[0 * 2 + 1] = -Math.exp(-x[0]);
                        jac[1 * 2 + 1] = -Math.exp(-x[1]);
                    }
                },

                // (4) Brown badly scaled
                new Problem("Brown badly scaled", 3, 2, new double[] { 1.0, 1.0 }, 0.0, 1.0e-8,
                        new double[] { 1.0e6, 2.0e-6 }, 1.0e-6, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        r[0] = x[0] - 1.0e6;
                        r[1] = x[1] - 2.0e-6;
                        r[2] = x[0] * x[1] - 2.0;
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        jac[0 * 3 + 0] = 1.0;
                        jac[1 * 3 + 0] = 0.0;
                        jac[0 * 3 + 1] = 0.0;
                        jac[1 * 3 + 1] = 1.0;
                        jac[0 * 3 + 2] = x[1];
                        jac[1 * 3 + 2] = x[0];
                    }
                },

                // (5) Beale
                new Problem("Beale", 3, 2, new double[] { 1.0, 1.0 }, 0.0, 1.0e-14, new double[] { 3.0, 0.5 }, 1.0e-7,
                        1000) {
                    private final double[] y = { 1.5, 2.25, 2.625 };

                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 3; i++) {
                            r[i] = y[i] - x[0] * (1.0 - Math.pow(x[1], i + 1));
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 3; i++) {
                            jac[0 * 3 + i] = -(1.0 - Math.pow(x[1], i + 1));
                            jac[1 * 3 + i] = (i + 1) * x[0] * Math.pow(x[1], i);
                        }
                    }
                },

                // (7) Helical valley
                new Problem("Helical valley", 3, 3, new double[] { -1.0, 0.0, 0.0 }, 0.0, 1.0e-14,
                        new double[] { 1.0, 0.0, 0.0 }, 1.0e-8, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
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
                    public void jacobianAt(double[] x, double[] jac) {
                        double tmp = x[0] * x[0] + x[1] * x[1];
                        double root = Math.sqrt(tmp);
                        double tpi = 2.0 * Math.PI;
                        jac[0 * 3 + 0] = 100.0 * x[1] / (tpi * tmp);
                        jac[1 * 3 + 0] = -100.0 * x[0] / (tpi * tmp);
                        jac[2 * 3 + 0] = 10.0;
                        jac[0 * 3 + 1] = 10.0 * x[0] / root;
                        jac[1 * 3 + 1] = 10.0 * x[1] / root;
                        jac[2 * 3 + 1] = 0.0;
                        jac[0 * 3 + 2] = 0.0;
                        jac[1 * 3 + 2] = 0.0;
                        jac[2 * 3 + 2] = 1.0;
                    }
                },

                // (8) Bard
                new Problem("Bard", 15, 3, new double[] { 1.0, 1.0, 1.0 }, 8.21487e-3, 1.0e-8, null, 0.0, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 15; i++) {
                            double u = i + 1;
                            double v = 15 - i;
                            double w = Math.min(u, v);
                            r[i] = BARD_Y[i] - (x[0] + u / (v * x[1] + w * x[2]));
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 15; i++) {
                            double u = i + 1;
                            double v = 15 - i;
                            double w = Math.min(u, v);
                            double d = v * x[1] + w * x[2];
                            jac[0 * 15 + i] = -1.0;
                            jac[1 * 15 + i] = u * v / (d * d);
                            jac[2 * 15 + i] = u * w / (d * d);
                        }
                    }
                },

                // (10) Meyer
                new Problem("Meyer", 16, 3, new double[] { 0.02, 4000.0, 250.0 }, 87.9458, 1.0e-3, null, 0.0, 4000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 16; i++) {
                            double t = 45.0 + 5.0 * (i + 1);
                            r[i] = x[0] * Math.exp(x[1] / (t + x[2])) - MEYER_Y[i];
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 16; i++) {
                            double t = 45.0 + 5.0 * (i + 1);
                            double d = t + x[2];
                            double e = Math.exp(x[1] / d);
                            jac[0 * 16 + i] = e;
                            jac[1 * 16 + i] = x[0] * e / d;
                            jac[2 * 16 + i] = -x[0] * e * x[1] / (d * d);
                        }
                    }
                },

                // (13) Powell singular -- the Jacobian is rank deficient at the
                // minimum, which is what makes it worth having here
                new Problem("Powell singular", 4, 4, new double[] { 3.0, -1.0, 0.0, 1.0 }, 0.0, 1.0e-12,
                        new double[] { 0.0, 0.0, 0.0, 0.0 }, 1.0e-3, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        r[0] = x[0] + 10.0 * x[1];
                        r[1] = Math.sqrt(5.0) * (x[2] - x[3]);
                        r[2] = (x[1] - 2.0 * x[2]) * (x[1] - 2.0 * x[2]);
                        r[3] = Math.sqrt(10.0) * (x[0] - x[3]) * (x[0] - x[3]);
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int k = 0; k < 16; k++) {
                            jac[k] = 0.0;
                        }
                        jac[0 * 4 + 0] = 1.0;
                        jac[1 * 4 + 0] = 10.0;
                        jac[2 * 4 + 1] = Math.sqrt(5.0);
                        jac[3 * 4 + 1] = -Math.sqrt(5.0);
                        jac[1 * 4 + 2] = 2.0 * (x[1] - 2.0 * x[2]);
                        jac[2 * 4 + 2] = -4.0 * (x[1] - 2.0 * x[2]);
                        jac[0 * 4 + 3] = 2.0 * Math.sqrt(10.0) * (x[0] - x[3]);
                        jac[3 * 4 + 3] = -2.0 * Math.sqrt(10.0) * (x[0] - x[3]);
                    }
                },

                // (15) Kowalik and Osborne
                new Problem("Kowalik and Osborne", 11, 4, new double[] { 0.25, 0.39, 0.415, 0.39 }, 3.07505e-4, 1.0e-9,
                        null, 0.0, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 11; i++) {
                            double u = KOWALIK_U[i];
                            double num = u * u + u * x[1];
                            double den = u * u + u * x[2] + x[3];
                            r[i] = KOWALIK_Y[i] - x[0] * num / den;
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 11; i++) {
                            double u = KOWALIK_U[i];
                            double num = u * u + u * x[1];
                            double den = u * u + u * x[2] + x[3];
                            jac[0 * 11 + i] = -num / den;
                            jac[1 * 11 + i] = -x[0] * u / den;
                            jac[2 * 11 + i] = x[0] * num * u / (den * den);
                            jac[3 * 11 + i] = x[0] * num / (den * den);
                        }
                    }
                },

                // (16) Brown and Dennis
                new Problem("Brown and Dennis", 20, 4, new double[] { 25.0, 5.0, -5.0, -1.0 }, 85822.2, 0.1, null, 0.0,
                        4000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 20; i++) {
                            double t = (i + 1) / 5.0;
                            double a = x[0] + t * x[1] - Math.exp(t);
                            double b = x[2] + x[3] * Math.sin(t) - Math.cos(t);
                            r[i] = a * a + b * b;
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 20; i++) {
                            double t = (i + 1) / 5.0;
                            double a = x[0] + t * x[1] - Math.exp(t);
                            double b = x[2] + x[3] * Math.sin(t) - Math.cos(t);
                            jac[0 * 20 + i] = 2.0 * a;
                            jac[1 * 20 + i] = 2.0 * a * t;
                            jac[2 * 20 + i] = 2.0 * b;
                            jac[3 * 20 + i] = 2.0 * b * Math.sin(t);
                        }
                    }
                },

                // (17) Osborne 1
                new Problem("Osborne 1", 33, 5, new double[] { 0.5, 1.5, -1.0, 0.01, 0.02 }, 5.46489e-5, 1.0e-9, null,
                        0.0, 1000) {
                    @Override
                    public void valueAt(double[] x, double[] r) {
                        for (int i = 0; i < 33; i++) {
                            double t = 10.0 * i;
                            r[i] = OSBORNE1_Y[i] - (x[0] + x[1] * Math.exp(-t * x[3]) + x[2] * Math.exp(-t * x[4]));
                        }
                    }

                    @Override
                    public void jacobianAt(double[] x, double[] jac) {
                        for (int i = 0; i < 33; i++) {
                            double t = 10.0 * i;
                            jac[0 * 33 + i] = -1.0;
                            jac[1 * 33 + i] = -Math.exp(-t * x[3]);
                            jac[2 * 33 + i] = -Math.exp(-t * x[4]);
                            jac[3 * 33 + i] = x[1] * t * Math.exp(-t * x[3]);
                            jac[4 * 33 + i] = x[2] * t * Math.exp(-t * x[4]);
                        }
                    }
                } };
    }
}
