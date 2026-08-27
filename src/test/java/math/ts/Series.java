package math.ts;

import math.linalg.DMatrix;

/**
 * Deterministic test data for {@link KalmanFilter} and {@link RtsSmoother}:
 * models to run and series drawn from them.
 * <p>
 * The generator is the inline LCG the rest of the test tree uses, so a seed
 * fixes the series exactly. Nothing here knows about either class under test --
 * a series drawn from the wrong model would still have to match
 * {@link StackedJoint}, since the likelihood of any data is what it is.
 */
final class Series {

    private long lcg;

    private Series(long seed) {
        this.lcg = seed;
    }

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) + 0.5) * 0x1.0p-53;
    }

    private double gauss() {
        double u1 = next();
        double u2 = next();
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    }

    static DMatrix matrix(int rows, int cols, double... values) {
        DMatrix a = new DMatrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                a.set(i, j, values[i * cols + j]);
            }
        }
        return a;
    }

    /** A random walk seen through noise, the model behind exponential smoothing. */
    static LinearGaussianModel localLevel(double q, double r, double p0) {
        return new LinearGaussianModel(matrix(1, 1, 1.0), matrix(1, 1, q), matrix(1, 1, 1.0),
                matrix(1, 1, r), new double[] { 0.0 }, matrix(1, 1, p0));
    }

    /** A level that moves at a slope that drifts. */
    static LinearGaussianModel trend() {
        return new LinearGaussianModel(matrix(2, 2, 1.0, 1.0, 0.0, 1.0),
                matrix(2, 2, 0.04, 0.0, 0.0, 0.01), matrix(1, 2, 1.0, 0.0), matrix(1, 1, 0.25),
                new double[] { 3.0, 0.5 }, matrix(2, 2, 10.0, 0.0, 0.0, 4.0));
    }

    /** Two states, two correlated sensors, everything full rank. */
    static LinearGaussianModel coupled() {
        return new LinearGaussianModel(matrix(2, 2, 0.9, 0.2, -0.1, 0.8),
                matrix(2, 2, 0.3, 0.1, 0.1, 0.2), matrix(2, 2, 1.0, 0.5, 0.0, 1.0),
                matrix(2, 2, 0.4, 0.15, 0.15, 0.6), new double[] { 1.0, -2.0 },
                matrix(2, 2, 2.0, 0.5, 0.5, 3.0));
    }

    /**
     * A level plus a three-period season through one sensor. Its {@code Q} is
     * singular -- the last seasonal coordinate is a pure carry with no shock of
     * its own, which is the ordinary shape of a structural model.
     */
    static LinearGaussianModel seasonal() {
        return new LinearGaussianModel(matrix(3, 3, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 1.0, 0.0),
                matrix(3, 3, 0.05, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0),
                matrix(1, 3, 1.0, 1.0, 0.0), matrix(1, 1, 0.3), new double[] { 5.0, 1.0, -1.0 },
                matrix(3, 3, 4.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0));
    }

    /** A state that never moves: {@code F = I} and no process noise at all. */
    static LinearGaussianModel staticState() {
        return new LinearGaussianModel(matrix(2, 2, 1.0, 0.0, 0.0, 1.0), new DMatrix(2, 2),
                matrix(2, 2, 1.0, 0.5, 0.0, 1.0), matrix(2, 2, 0.4, 0.15, 0.15, 0.6),
                new double[] { 1.0, -2.0 }, matrix(2, 2, 2.0, 0.5, 0.5, 3.0));
    }

    /** Draws {@code n} observations from the model, with its matrices of the step. */
    static DMatrix draw(StateSpaceModel model, int n, long seed) {
        Series rng = new Series(seed);
        int d = model.stateDimension();
        int m = model.observationDimension();
        DMatrix lp = factor(model.initialCovariance());

        double[] x = new double[d];
        model.initialMean(x);
        for (int i = d - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int k = 0; k <= i; ++k) {
                sum += lp.get(i, k) * rng.gauss();
            }
            x[i] += sum;
        }

        DMatrix y = new DMatrix(n, m);
        double[] moved = new double[d];
        for (int t = 0; t < n; ++t) {
            if (t > 0) {
                DMatrix f = model.transition(t);
                DMatrix lq = factor(model.processNoise(t));
                for (int i = 0; i < d; ++i) {
                    double sum = 0.0;
                    for (int j = 0; j < d; ++j) {
                        sum += f.get(i, j) * x[j];
                    }
                    for (int k = 0; k <= i; ++k) {
                        sum += lq.get(i, k) * rng.gauss();
                    }
                    moved[i] = sum;
                }
                System.arraycopy(moved, 0, x, 0, d);
            }
            DMatrix h = model.observation(t);
            DMatrix lr = factor(model.observationNoise(t));
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int j = 0; j < d; ++j) {
                    sum += h.get(i, j) * x[j];
                }
                for (int k = 0; k <= i; ++k) {
                    sum += lr.get(i, k) * rng.gauss();
                }
                y.set(t, i, sum);
            }
        }
        return y;
    }

    /**
     * A local level model sampled at irregular times: the state drifts in
     * continuous time, so the process noise of a step is proportional to the gap
     * it spans. The motivating case for time-varying matrices, and one where
     * only {@code Q} varies.
     */
    static TimeVaryingModel irregular(double qPerUnit, double r, double p0, double[] gaps) {
        int n = gaps.length;
        DMatrix[] q = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            q[t] = matrix(1, 1, qPerUnit * gaps[t]);
        }
        return new TimeVaryingModel(TimeVaryingModel.repeat(matrix(1, 1, 1.0), n), q,
                TimeVaryingModel.repeat(matrix(1, 1, 1.0), n),
                TimeVaryingModel.repeat(matrix(1, 1, r), n), new double[] { 0.0 },
                matrix(1, 1, p0));
    }

    /** The gaps of an irregular schedule, deterministic and uneven. */
    static double[] gaps(int n, long seed) {
        Series rng = new Series(seed);
        double[] out = new double[n];
        for (int t = 0; t < n; ++t) {
            // between a quarter and four time units, so Q spans four orders
            out[t] = Math.pow(4.0, 2.0 * rng.next() - 1.0);
        }
        return out;
    }

    /**
     * Two states seen through a sensor that changes: on even steps it reads the
     * first coordinate accurately, on odd steps a combination of both, badly.
     * Only {@code H} and {@code R} vary.
     */
    static TimeVaryingModel switchingSensor(int n) {
        DMatrix[] h = new DMatrix[n];
        DMatrix[] r = new DMatrix[n];
        DMatrix sharp = matrix(1, 2, 1.0, 0.0);
        DMatrix blunt = matrix(1, 2, 0.5, 0.5);
        DMatrix small = matrix(1, 1, 0.05);
        DMatrix large = matrix(1, 1, 1.5);
        for (int t = 0; t < n; ++t) {
            h[t] = (t % 2 == 0) ? sharp : blunt;
            r[t] = (t % 2 == 0) ? small : large;
        }
        return new TimeVaryingModel(TimeVaryingModel.repeat(matrix(2, 2, 1.0, 1.0, 0.0, 1.0), n),
                TimeVaryingModel.repeat(matrix(2, 2, 0.04, 0.0, 0.0, 0.01), n), h, r,
                new double[] { 3.0, 0.5 }, matrix(2, 2, 10.0, 0.0, 0.0, 4.0));
    }

    /**
     * A state rotated by a different angle at every step, with a process noise
     * that grows. Everything varies at once, and {@code F} is not symmetric, so
     * a transposed transition would show up.
     */
    static TimeVaryingModel rotating(int n) {
        DMatrix[] f = new DMatrix[n];
        DMatrix[] q = new DMatrix[n];
        DMatrix[] h = new DMatrix[n];
        DMatrix[] r = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            double angle = 0.3 + 0.11 * t;
            double c = Math.cos(angle);
            double s = Math.sin(angle);
            f[t] = matrix(2, 2, 0.95 * c, -0.95 * s, 0.95 * s, 0.95 * c);
            q[t] = matrix(2, 2, 0.02 * (1 + t), 0.0, 0.0, 0.01 * (1 + t));
            h[t] = matrix(2, 2, 1.0, 0.2 * t, 0.0, 1.0);
            r[t] = matrix(2, 2, 0.3 + 0.05 * t, 0.05, 0.05, 0.4);
        }
        return new TimeVaryingModel(f, q, h, r, new double[] { 1.0, -2.0 },
                matrix(2, 2, 2.0, 0.5, 0.5, 3.0));
    }

    /** The same model as {@code constant}, expressed as a schedule that repeats. */
    static TimeVaryingModel repeated(LinearGaussianModel constant, int n) {
        return new TimeVaryingModel(TimeVaryingModel.repeat(constant.transition(), n),
                TimeVaryingModel.repeat(constant.processNoise(), n),
                TimeVaryingModel.repeat(constant.observation(), n),
                TimeVaryingModel.repeat(constant.observationNoise(), n), meanOf(constant),
                constant.initialCovariance());
    }

    static double[] meanOf(StateSpaceModel model) {
        double[] out = new double[model.stateDimension()];
        model.initialMean(out);
        return out;
    }

    /**
     * A Cholesky factor that tolerates a singular argument, which
     * {@link math.linalg.CholeskyDecomp} does not and should not -- but a
     * singular {@code Q} is exactly what has to be simulated from.
     */
    private static DMatrix factor(DMatrix a) {
        int n = a.numRows();
        DMatrix l = new DMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = a.get(i, j);
                for (int k = 0; k < j; ++k) {
                    sum -= l.get(i, k) * l.get(j, k);
                }
                if (i == j) {
                    l.set(i, i, sum > 0.0 ? Math.sqrt(sum) : 0.0);
                } else {
                    double djj = l.get(j, j);
                    l.set(i, j, djj > 0.0 ? sum / djj : 0.0);
                }
            }
        }
        return l;
    }
}
