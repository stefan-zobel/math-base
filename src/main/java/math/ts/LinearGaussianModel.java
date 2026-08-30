package math.ts;

import math.linalg.DMatrix;

/**
 * The linear Gaussian state space model with matrices that do not vary over
 * time:
 * <p>
 * <code>x[t] = F x[t-1] + w[t]</code>, <code>w[t] ~ N(0, Q)</code><br>
 * <code>y[t] = H x[t] + v[t]</code>, <code>v[t] ~ N(0, R)</code><br>
 * <code>x[0] ~ N(m0, P0)</code>
 * <p>
 * A hidden state that moves on its own and an observation that sees a linear
 * function of it through noise. The state is what a caller cares about and
 * never gets to see; the observation is what arrives. Everything that follows
 * from this model is Gaussian, which is why {@link KalmanFilter} can be exact
 * where a general posterior needs quadrature or a sampler.
 * <p>
 * Constant matrices cover the structural models -- a random walk plus noise, a
 * local linear trend, a seasonal component, an ARMA model in state space form
 * -- and they cover a sensor that drops out and comes back, which is a
 * {@code NaN} in the data rather than a change to the model. What they do not
 * cover is unequal time steps or a measurement setup that genuinely changes;
 * for those, {@link TimeVaryingModel} implements the same
 * {@link StateSpaceModel} interface from arrays.
 * <p>
 * <b>{@code Q} may be singular and usually is.</b> The smooth trend model puts
 * its only shock on the slope, and the state space form of an ARMA model has
 * {@code Q = s r r'} of rank one whatever its order. So {@code Q} and {@code R}
 * are checked for symmetry and for a diagonal that is not negative, and are
 * never factored; {@code P0} is factored, which is what proves it positive
 * definite. A {@code Q} or {@code R} that is symmetric but still not positive
 * semidefinite is not caught here -- it surfaces at the first step where the
 * innovation covariance stops being positive definite, and that message names
 * the time index. The check is left out on purpose: a maximum likelihood fit
 * builds one of these per function evaluation.
 * <p>
 * Instances are immutable in the sense that matters: every matrix is copied
 * when the model is built, so nothing a caller keeps a reference to can change
 * it afterwards. <b>What comes back out is the model's own and must not be
 * written to</b> -- the filter reads all four matrices at every step, and
 * copying them there would allocate for nothing. This is the convention
 * {@link KalmanFilter.Result} already follows with its public arrays.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/State-space_representation">
 * Wikipedia State-space representation</a>.
 *
 * @since 1.5.3
 */
public final class LinearGaussianModel implements StateSpaceModel {

    private final DMatrix f;
    private final DMatrix q;
    private final DMatrix h;
    private final DMatrix r;
    private final double[] m0;
    private final DMatrix p0;
    private final int d;
    private final int m;

    /**
     * The model with the given transition, process noise, observation,
     * observation noise, initial mean and initial covariance.
     *
     * @param F
     *            the transition matrix, square of order {@code d}, finite. Not
     *            modified
     * @param Q
     *            the process noise covariance, {@code d} by {@code d},
     *            symmetric, with a diagonal that is not negative. May be
     *            singular. Not modified
     * @param H
     *            the observation matrix, {@code m} by {@code d}, finite. Not
     *            modified
     * @param R
     *            the observation noise covariance, {@code m} by {@code m},
     *            symmetric, with a diagonal that is not negative. May be
     *            singular. Not modified
     * @param m0
     *            the mean of the initial state, of length {@code d}, finite.
     *            Not modified
     * @param P0
     *            the covariance of the initial state, {@code d} by {@code d},
     *            symmetric and positive definite. Not modified
     * @throws IllegalArgumentException
     *             if any argument is {@code null}, if any dimension disagrees
     *             with another, if {@code d} or {@code m} is zero, if any entry
     *             is not finite, if {@code Q}, {@code R} or {@code P0} is not
     *             symmetric within tolerance or has a negative variance on its
     *             diagonal, or if {@code P0} is not positive definite
     */
    public LinearGaussianModel(DMatrix F, DMatrix Q, DMatrix H, DMatrix R, double[] m0, DMatrix P0) {
        if (F == null || Q == null || H == null || R == null || P0 == null) {
            throw new IllegalArgumentException("no matrix argument may be null");
        }
        if (!F.isSquareMatrix()) {
            throw new IllegalArgumentException("F is not square");
        }
        int dim = F.numRows();
        if (dim == 0) {
            throw new IllegalArgumentException("the state must have at least one coordinate");
        }
        int obs = H.numRows();
        if (obs == 0) {
            throw new IllegalArgumentException("the observation must have at least one coordinate");
        }

        this.f = Models.transition(F, dim, "F");
        this.q = Models.symmetrized(Q, dim, "Q");
        this.h = Models.observation(H, obs, dim, "H");
        this.r = Models.symmetrized(R, obs, "R");
        this.m0 = Models.initialMean(m0, dim);
        this.p0 = Models.initialCovariance(P0, dim);
        this.d = dim;
        this.m = obs;
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
    public DMatrix transition(int t) {
        return f;
    }

    @Override
    public DMatrix processNoise(int t) {
        return q;
    }

    @Override
    public DMatrix observation(int t) {
        return h;
    }

    @Override
    public DMatrix observationNoise(int t) {
        return r;
    }

    /**
     * The transition matrix {@code F}, the same one {@link #transition(int)}
     * returns at every time step.
     *
     * @return {@code F}, square of order {@link #stateDimension()}, not to be
     *         modified
     */
    public DMatrix transition() {
        return f;
    }

    /**
     * The process noise covariance {@code Q}, symmetrized as it is used and the
     * same one {@link #processNoise(int)} returns at every time step.
     *
     * @return {@code Q}, square of order {@link #stateDimension()}, not to be
     *         modified
     */
    public DMatrix processNoise() {
        return q;
    }

    /**
     * The observation matrix {@code H}, the same one
     * {@link #observation(int)} returns at every time step.
     *
     * @return {@code H}, {@link #observationDimension()} by
     *         {@link #stateDimension()}, not to be modified
     */
    public DMatrix observation() {
        return h;
    }

    /**
     * The observation noise covariance {@code R}, symmetrized as it is used and
     * the same one {@link #observationNoise(int)} returns at every time step.
     *
     * @return {@code R}, square of order {@link #observationDimension()}, not
     *         to be modified
     */
    public DMatrix observationNoise() {
        return r;
    }

    @Override
    public void initialMean(double[] out) {
        Models.checkOut(out, d);
        System.arraycopy(m0, 0, out, 0, d);
    }

    /**
     * The covariance of the initial state, symmetrized as it is used.
     *
     * @return {@code P0}, square of order {@link #stateDimension()}, not to be
     *         modified
     */
    @Override
    public DMatrix initialCovariance() {
        return p0;
    }
}
