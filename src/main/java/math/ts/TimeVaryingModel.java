package math.ts;

import math.linalg.DMatrix;

/**
 * A {@link StateSpaceModel} whose four matrices are read from arrays, one entry
 * per time step.
 * <p>
 * What this buys over {@link LinearGaussianModel} is the two cases constant
 * matrices cannot express:
 * <ul>
 * <li><b>Unequal time steps.</b> A state that drifts in continuous time has a
 * transition and a process noise that depend on the gap -- a random walk
 * observed after two days has twice the variance of one observed after one --
 * so a series sampled irregularly needs {@code F(t)} and {@code Q(t)} even
 * though the underlying model never changes.</li>
 * <li><b>A measurement setup that changes.</b> A sensor that is recalibrated, a
 * survey whose question changes, an instrument whose precision improves: that
 * is {@code H(t)} and {@code R(t)}.</li>
 * </ul>
 * <p>
 * <b>It is not for a sensor that drops out.</b> That is a {@code NaN} in the
 * data, which {@link KalmanFilter} handles at no cost and without any change of
 * model -- and which is why the observation dimension is fixed here. Varying
 * matrices are for a model that changes; a {@code NaN} is for data that is
 * missing, and reaching for the first when the second is meant produces a model
 * that is harder to write and no more general.
 * <p>
 * <b>{@code F[0]} and {@code Q[0]} are never read.</b> {@code F(t)} is the map
 * that leads <em>into</em> time {@code t} and nothing leads into the first time
 * point; the prior does that job. They still have to be present and valid, so
 * that all four arrays have one length and no caller has to hold an offset of
 * one in their head. {@link #repeat(DMatrix, int)} makes filling them free.
 * <p>
 * Entries are copied when the model is built, and consecutive entries that are
 * the <em>same instance</em> are stored once -- so {@link #repeat(DMatrix, int)}
 * costs one matrix and not {@code n} of them, and a model where only {@code H}
 * varies pays for the varying part alone. What comes back out is the model's
 * own and must not be written to, as {@link StateSpaceModel} states.
 *
 * @since 1.5.3
 */
public final class TimeVaryingModel implements StateSpaceModel {

    private final DMatrix[] f;
    private final DMatrix[] q;
    private final DMatrix[] h;
    private final DMatrix[] r;
    private final double[] m0;
    private final DMatrix p0;
    private final int d;
    private final int m;
    private final int n;

    /**
     * The model with the given schedules of transition, process noise,
     * observation and observation noise, and the given initial state.
     * <p>
     * All four arrays have the same length, which is how many time points the
     * model is defined for. The dimensions do not vary: every {@code F[t]} is
     * square of the order {@code F[0]} has, and every {@code H[t]} has the rows
     * and columns {@code H[0]} has.
     *
     * @param F
     *            the transition matrices, one per time step, each square of
     *            order {@code d} and finite. {@code F[0]} is validated and
     *            never used. Neither the array nor its entries are modified
     * @param Q
     *            the process noise covariances, one per time step, each
     *            {@code d} by {@code d}, symmetric, with a diagonal that is not
     *            negative. Each may be singular. {@code Q[0]} is validated and
     *            never used. Not modified
     * @param H
     *            the observation matrices, one per time step, each {@code m} by
     *            {@code d} and finite. Not modified
     * @param R
     *            the observation noise covariances, one per time step, each
     *            {@code m} by {@code m}, symmetric, with a diagonal that is not
     *            negative. Each may be singular. Not modified
     * @param m0
     *            the mean of the initial state, of length {@code d}, finite.
     *            Not modified
     * @param P0
     *            the covariance of the initial state, {@code d} by {@code d},
     *            symmetric and positive definite. Not modified
     * @throws IllegalArgumentException
     *             if any argument or any entry is {@code null}, if the four
     *             arrays are not all of the same non-zero length, if any
     *             dimension disagrees with another or varies with the time
     *             step, if any entry is not finite, if a covariance is not
     *             symmetric within tolerance or has a negative variance on its
     *             diagonal, or if {@code P0} is not positive definite
     */
    public TimeVaryingModel(DMatrix[] F, DMatrix[] Q, DMatrix[] H, DMatrix[] R, double[] m0,
            DMatrix P0) {
        if (F == null || Q == null || H == null || R == null || P0 == null) {
            throw new IllegalArgumentException("no matrix argument may be null");
        }
        int length = F.length;
        if (length == 0) {
            throw new IllegalArgumentException("the series must have at least one time point");
        }
        if (Q.length != length || H.length != length || R.length != length) {
            throw new IllegalArgumentException("F, Q, H and R are of lengths " + length + ", "
                    + Q.length + ", " + H.length + " and " + R.length + ", which are not all equal");
        }
        checkEntries(F, "F");
        checkEntries(Q, "Q");
        checkEntries(H, "H");
        checkEntries(R, "R");

        if (!F[0].isSquareMatrix()) {
            throw new IllegalArgumentException("F[0] is not square");
        }
        int dim = F[0].numRows();
        if (dim == 0) {
            throw new IllegalArgumentException("the state must have at least one coordinate");
        }
        int obs = H[0].numRows();
        if (obs == 0) {
            throw new IllegalArgumentException("the observation must have at least one coordinate");
        }

        this.f = new DMatrix[length];
        this.q = new DMatrix[length];
        this.h = new DMatrix[length];
        this.r = new DMatrix[length];
        for (int t = 0; t < length; ++t) {
            // a schedule that repeats an instance is stored once, which is what
            // makes repeat(matrix, n) cost one matrix instead of n
            f[t] = shared(F, f, t);
            if (f[t] == null) {
                f[t] = Models.transition(F[t], dim, "F[" + t + "]");
            }
            q[t] = shared(Q, q, t);
            if (q[t] == null) {
                q[t] = Models.symmetrized(Q[t], dim, "Q[" + t + "]");
            }
            h[t] = shared(H, h, t);
            if (h[t] == null) {
                h[t] = Models.observation(H[t], obs, dim, "H[" + t + "]");
            }
            r[t] = shared(R, r, t);
            if (r[t] == null) {
                r[t] = Models.symmetrized(R[t], obs, "R[" + t + "]");
            }
        }
        this.m0 = Models.initialMean(m0, dim);
        this.p0 = Models.initialCovariance(P0, dim);
        this.d = dim;
        this.m = obs;
        this.n = length;
    }

    /**
     * An array holding the same matrix at every one of {@code n} time steps.
     * <p>
     * What keeps "only {@code H} varies" to one call per constant matrix rather
     * than three loops, and what fills the {@code F[0]} and {@code Q[0]} slots
     * that are never read. The entries are the same instance, which the
     * constructor recognizes and stores once.
     *
     * @param a
     *            the matrix to repeat. Not modified
     * @param n
     *            how many time steps, one or more
     * @return an array of length {@code n}, every entry {@code a}
     * @throws IllegalArgumentException
     *             if {@code a} is {@code null} or {@code n} is not positive
     */
    public static DMatrix[] repeat(DMatrix a, int n) {
        if (a == null) {
            throw new IllegalArgumentException("a must not be null");
        }
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive, not " + n);
        }
        DMatrix[] out = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            out[t] = a;
        }
        return out;
    }

    @Override
    public int stateDimension() {
        return d;
    }

    @Override
    public int observationDimension() {
        return m;
    }

    @Override
    public int length() {
        return n;
    }

    @Override
    public DMatrix transition(int t) {
        return f[t];
    }

    @Override
    public DMatrix processNoise(int t) {
        return q[t];
    }

    @Override
    public DMatrix observation(int t) {
        return h[t];
    }

    @Override
    public DMatrix observationNoise(int t) {
        return r[t];
    }

    @Override
    public void initialMean(double[] out) {
        Models.checkOut(out, d);
        System.arraycopy(m0, 0, out, 0, d);
    }

    @Override
    public DMatrix initialCovariance() {
        return p0;
    }

    /**
     * The already validated copy to reuse if this entry is the same instance as
     * the one before it, {@code null} if it has to be checked on its own.
     */
    private static DMatrix shared(DMatrix[] given, DMatrix[] stored, int t) {
        return t > 0 && given[t] == given[t - 1] ? stored[t - 1] : null;
    }

    private static void checkEntries(DMatrix[] a, String name) {
        for (int t = 0; t < a.length; ++t) {
            if (a[t] == null) {
                throw new IllegalArgumentException(name + "[" + t + "] must not be null");
            }
        }
    }
}
