package math.ts;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;

/**
 * The Rauch-Tung-Striebel smoother: the exact posterior of the hidden state
 * given the <b>whole</b> series, <code>p(x[t] | y[1..n])</code>.
 * <p>
 * {@link KalmanFilter} answers with everything observed so far, which is all a
 * live feed has. Once the series is finished that restriction is artificial:
 * what happened at time twenty says something about where the level was at time
 * ten, and a filter is not allowed to use it. The smoother is that second
 * reading of the past, and it is not a heuristic revision -- it is the exact
 * conditional distribution, obtained by running backwards through the
 * quantities the forward pass already produced.
 * <p>
 * The recursion is
 * <code>C[t] = P[t] F(t+1)' P[t+1|t]^-1</code>,
 * <code>m[t] &lt;- m[t] + C[t] (m[t+1] - m[t+1|t])</code>,
 * <code>P[t] &lt;- P[t] + C[t] (P[t+1] - P[t+1|t]) C[t]'</code>,
 * started from the filtered state at the last index, which is already smoothed
 * -- there is nothing after it. So the two passes agree exactly there and
 * nowhere else, and everywhere else the smoothed covariance is the smaller of
 * the two: more evidence cannot make a Gaussian posterior wider.
 * <p>
 * <code>C[t]</code> is kept rather than discarded. It costs nothing, it is what
 * the lag-one covariance
 * <code>Cov(x[t], x[t+1] | y) = P[t+1] C[t]'</code> is built from, and that in
 * turn is what an expectation-maximization fit of the model needs.
 * <p>
 * Missing observations need no mention here: they were handled in the forward
 * pass, and a time point where nothing was seen is simply one where the
 * filtered and predicted states coincide.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Kalman_filter#Rauch%E2%80%93Tung%E2%80%93Striebel">
 * Wikipedia Rauch-Tung-Striebel</a>.
 *
 * @since 1.5.3
 */
public final class RtsSmoother {

    /**
     * Filters the series and smooths it in one call.
     *
     * @param model
     *            the model. Not modified
     * @param y
     *            the observations, as
     *            {@link KalmanFilter#filter(StateSpaceModel, DMatrix)}
     *            takes them. Not modified
     * @return the smoothed history
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link KalmanFilter#filter(StateSpaceModel, DMatrix)}
     *             states
     * @throws ArithmeticException
     *             if an innovation covariance or a predicted covariance is not
     *             positive definite
     */
    public static Result smooth(StateSpaceModel model, DMatrix y) {
        return smooth(model, KalmanFilter.filter(model, y));
    }

    /**
     * Runs the backward pass over a forward pass that has already been done.
     * <p>
     * The form to use when the filtered history is wanted as well: filtering
     * twice would be the only alternative.
     *
     * @param model
     *            the model the forward pass was run with. Not modified
     * @param filtered
     *            what {@link KalmanFilter#filter(StateSpaceModel, DMatrix)}
     *            returned. Not modified
     * @return the smoothed history
     * @throws IllegalArgumentException
     *             if either argument is {@code null}, or if the two do not
     *             agree on the dimension of the state
     * @throws ArithmeticException
     *             if a predicted covariance is not positive definite, which
     *             leaves the backward gain undefined
     */
    public static Result smooth(StateSpaceModel model, KalmanFilter.Result filtered) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        if (filtered == null) {
            throw new IllegalArgumentException("filtered must not be null");
        }
        int d = model.stateDimension();
        int n = filtered.length;
        if (filtered.filteredMean[0].length != d) {
            throw new IllegalArgumentException("the filtered state has "
                    + filtered.filteredMean[0].length + " coordinates, not " + d);
        }
        if (n > model.length()) {
            throw new IllegalArgumentException("the filtered history has " + n
                    + " time points and the model is defined for " + model.length());
        }

        double[][] mean = new double[n][];
        DMatrix[] covariance = new DMatrix[n];
        DMatrix[] gain = new DMatrix[n];

        // nothing comes after the last index, so it is already smoothed
        mean[n - 1] = filtered.filteredMean[n - 1].clone();
        covariance[n - 1] = filtered.filteredCovariance[n - 1].copy();

        double[] rhs = new double[d];
        double[] sol = new double[d];
        for (int t = n - 2; t >= 0; --t) {
            DMatrix predicted = filtered.predictedCovariance[t + 1];
            DMatrix factor;
            try {
                factor = CholeskyDecomp.cholesky(predicted);
            } catch (RuntimeException ex) {
                throw new ArithmeticException("the predicted covariance at step " + (t + 1)
                        + " is not positive definite : " + ex.getMessage());
            }
            for (int i = 0; i < d; ++i) {
                if (!(factor.getUnsafe(i, i) > 0.0)) {
                    throw new ArithmeticException("the predicted covariance at step " + (t + 1)
                            + " is singular: its factor has a zero at (" + i + ", " + i + ")");
                }
            }

            // C = P F' Pp^-1, so C' = Pp^-1 (P F')' : one solve per row of P F'.
            // F(t+1) and not F(t) -- the same transition that produced the
            // predicted covariance this inverts
            DMatrix pft = filtered.filteredCovariance[t].mulBTrans(model.transition(t + 1));
            DMatrix c = new DMatrix(d, d);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    rhs[j] = pft.getUnsafe(i, j);
                }
                CholeskyDecomp.solve(factor, rhs, sol);
                for (int j = 0; j < d; ++j) {
                    c.setUnsafe(i, j, sol[j]);
                }
            }
            gain[t] = c;

            mean[t] = filtered.filteredMean[t].clone();
            for (int i = 0; i < d; ++i) {
                double sum = 0.0;
                for (int j = 0; j < d; ++j) {
                    sum += c.getUnsafe(i, j) * (mean[t + 1][j] - filtered.predictedMean[t + 1][j]);
                }
                mean[t][i] += sum;
            }

            DMatrix difference = covariance[t + 1].minus(predicted);
            covariance[t] = symmetrize(
                    filtered.filteredCovariance[t].add(c.mul(difference).mulBTrans(c)));
        }
        return new Result(mean, covariance, gain);
    }

    /**
     * The whole smoothed history of one series.
     * <p>
     * Index {@code t} runs over the same time points the forward pass used.
     *
     * @since 1.5.3
     */
    public static final class Result {

        /** How many time points the series has. */
        public final int length;
        /** The mean of <code>x[t] | y[1..n]</code>, one row per time point. */
        public final double[][] mean;
        /** The covariance of <code>x[t] | y[1..n]</code>. */
        public final DMatrix[] covariance;
        /**
         * The backward gain <code>C[t]</code>, from which the lag-one
         * covariance <code>Cov(x[t], x[t+1] | y) = C[t] P[t+1]</code> follows.
         * {@code null} at the last index, where there is no next state.
         */
        public final DMatrix[] gain;

        Result(double[][] mean, DMatrix[] covariance, DMatrix[] gain) {
            this.length = mean.length;
            this.mean = mean;
            this.covariance = covariance;
            this.gain = gain;
        }

        /**
         * The smoothed state at one time point, as the distribution it is.
         *
         * @param t
         *            the time point, from zero to {@link #length} minus one
         * @return the normal distribution of <code>x[t] | y[1..n]</code>
         * @throws IndexOutOfBoundsException
         *             if {@code t} is outside the series
         * @throws IllegalArgumentException
         *             if that covariance is singular
         */
        public MultivariateNormal smoothedState(int t) {
            return new MultivariateNormal(mean[t].clone(), covariance[t]);
        }

        /**
         * The covariance between the state at {@code t} and the state at
         * {@code t + 1}, given the whole series.
         * <p>
         * <code>C[t] P[t+1]</code>, in that order: this is not symmetric, and
         * the transpose is the covariance of the pair the other way round.
         *
         * @param t
         *            the time point, from zero to {@link #length} minus two
         * @return <code>Cov(x[t], x[t+1] | y[1..n])</code>, a fresh matrix
         * @throws IndexOutOfBoundsException
         *             if {@code t} is outside the series
         * @throws IllegalArgumentException
         *             if {@code t} is the last index, where there is no next
         *             state
         */
        public DMatrix lagOneCovariance(int t) {
            if (t == length - 1) {
                throw new IllegalArgumentException(
                        "there is no state after the last index " + t);
            }
            return gain[t].mul(covariance[t + 1]);
        }
    }

    private RtsSmoother() {
        throw new AssertionError();
    }

    /** Averages a matrix with its own transpose, in place. */
    private static DMatrix symmetrize(DMatrix A) {
        int n = A.numRows();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double v = 0.5 * (A.getUnsafe(i, j) + A.getUnsafe(j, i));
                A.setUnsafe(i, j, v);
                A.setUnsafe(j, i, v);
            }
        }
        return A;
    }
}
