package math.ts;

import math.distribution.MultivariateNormal;
import math.linalg.DMatrix;

/**
 * A linear Gaussian state space model, whose matrices may depend on the time
 * step:
 * <p>
 * <code>x[t] = F(t) x[t-1] + w[t]</code>, <code>w[t] ~ N(0, Q(t))</code><br>
 * <code>y[t] = H(t) x[t] + v[t]</code>, <code>v[t] ~ N(0, R(t))</code><br>
 * <code>x[0] ~ N(m0, P0)</code>
 * <p>
 * Two implementations come with the package -- {@link LinearGaussianModel},
 * whose matrices do not vary, and {@link TimeVaryingModel}, which reads them
 * from arrays. A third case is what this interface exists for: a caller who
 * derives <code>F(dt)</code> and <code>Q(dt)</code> from a continuous time
 * model can implement it directly, and does not have to hand the library a time
 * grid it would have no use for.
 * <p>
 * <b>{@code F(t)} is the map that leads <em>into</em> time {@code t}</b>, so it
 * is read for {@code t >= 1} and never at zero: nothing transitions into the
 * first time point, which is what the prior is for. The same holds for
 * {@code Q(t)}. {@code H(t)} and {@code R(t)} are read at every index including
 * zero.
 * <p>
 * <b>The dimensions do not vary with {@code t}.</b> Every {@code F(t)} is
 * square of order {@link #stateDimension()} and every {@code H(t)} is
 * {@link #observationDimension()} by that. A sensor that only starts reporting
 * halfway through the series -- or stops halfway through -- is <em>not</em> a
 * change of shape: it is a {@code NaN} in the data, which
 * {@link KalmanFilter} already handles for free. The first of those two is the
 * one worth naming, because it is the one that looks like a taller {@code H}:
 * count the sensor from the beginning and mark the readings it did not yet
 * produce as missing. Varying
 * matrices are for a model that changes; a {@code NaN} is for data that is
 * missing, and the two are worth keeping apart.
 * <p>
 * <b>The matrices come back live.</b> An implementation may return the same
 * instance on every call and is expected to, since the filter reads all four at
 * every step; a caller must not write into what it is handed. This is the
 * convention {@link KalmanFilter.Result} already follows with its public
 * arrays.
 *
 * @since 1.5.3
 */
public interface StateSpaceModel {

    /**
     * The number of coordinates of the hidden state.
     *
     * @return {@code d}, one or more, the same at every time step
     */
    int stateDimension();

    /**
     * The number of coordinates of an observation.
     *
     * @return {@code m}, one or more, the same at every time step
     */
    int observationDimension();

    /**
     * The transition matrix that leads into time {@code t}.
     *
     * @param t
     *            the time step, from one to {@link #length()} minus one
     * @return {@code F(t)}, square of order {@link #stateDimension()}, not to
     *         be modified
     * @throws IndexOutOfBoundsException
     *             if the model is not defined at {@code t}
     */
    DMatrix transition(int t);

    /**
     * The covariance of the process noise added in the step into time
     * {@code t}. May be singular, and for the usual structural models is.
     *
     * @param t
     *            the time step, from one to {@link #length()} minus one
     * @return {@code Q(t)}, square of order {@link #stateDimension()}, not to
     *         be modified
     * @throws IndexOutOfBoundsException
     *             if the model is not defined at {@code t}
     */
    DMatrix processNoise(int t);

    /**
     * The matrix through which the state is observed at time {@code t}.
     *
     * @param t
     *            the time step, from zero to {@link #length()} minus one
     * @return {@code H(t)}, {@link #observationDimension()} by
     *         {@link #stateDimension()}, not to be modified
     * @throws IndexOutOfBoundsException
     *             if the model is not defined at {@code t}
     */
    DMatrix observation(int t);

    /**
     * The covariance of the observation noise at time {@code t}.
     *
     * @param t
     *            the time step, from zero to {@link #length()} minus one
     * @return {@code R(t)}, square of order {@link #observationDimension()},
     *         not to be modified
     * @throws IndexOutOfBoundsException
     *             if the model is not defined at {@code t}
     */
    DMatrix observationNoise(int t);

    /**
     * The mean of the initial state.
     *
     * @param out
     *            where the mean is written, of length
     *            {@link #stateDimension()}. Its previous contents are
     *            overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or of the wrong length
     */
    void initialMean(double[] out);

    /**
     * The covariance of the initial state, which has to be positive definite.
     *
     * @return {@code P0}, square of order {@link #stateDimension()}, not to be
     *         modified
     */
    DMatrix initialCovariance();

    /**
     * How many time points this model is defined for.
     * <p>
     * The default is {@link Integer#MAX_VALUE}, which is what a model with
     * constant matrices should report: it is defined for as long as the data
     * lasts, and a bound check against it is one comparison that never fails.
     * An implementation backed by a finite schedule <b>has to override this</b>
     * -- otherwise the filter cannot tell a caller that the series outran the
     * model and the schedule will simply run off its own end.
     *
     * @return the number of time points, or {@link Integer#MAX_VALUE} if the
     *         model does not run out
     */
    default int length() {
        return Integer.MAX_VALUE;
    }

    /**
     * The initial state as the distribution it is, {@code N(m0, P0)}.
     *
     * @return the prior over the state before any observation
     * @throws IllegalArgumentException
     *             if the initial covariance is singular
     */
    default MultivariateNormal initialState() {
        double[] mean = new double[stateDimension()];
        initialMean(mean);
        return new MultivariateNormal(mean, initialCovariance());
    }
}
