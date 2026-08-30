package math.ts;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;

/**
 * The Kalman filter: the exact posterior of the hidden state of a
 * {@link StateSpaceModel} given the observations up to now.
 * <p>
 * Everything in the model is Gaussian and every operation on it is linear, so
 * <code>p(x[t] | y[1..t])</code> is a normal distribution and the filter
 * computes its mean and covariance in closed form. There is nothing to
 * approximate and nothing to sample. That is what separates this from the rest
 * of the Bayesian machinery in the library: {@code math.stats.bayes} normalizes
 * a posterior by quadrature or fits a Gaussian at its mode because it has to,
 * and here the answer <em>is</em> a Gaussian.
 * <p>
 * Each time step is a prediction, <code>m &lt;- F(t) m</code> and
 * <code>P &lt;- F(t) P F(t)' + Q(t)</code>, followed by a correction against
 * the observation. The matrices are read at every step, so a
 * {@link TimeVaryingModel} costs nothing extra here and a
 * {@link LinearGaussianModel} simply returns the same ones every time. The
 * correction runs through the innovation
 * <code>e = y - H m</code> and its covariance
 * <code>S = H P H' + R</code>, and those two are worth more than the state they
 * produce: {@code e' S^-1 e} is a squared Mahalanobis distance with a known
 * chi-squared law, which makes it an anomaly score rather than a residual, and
 * the sum of {@code log N(e; 0, S)} over the series is the log likelihood of
 * the model -- the prediction error decomposition, and what a maximum
 * likelihood fit of the model maximizes.
 * <p>
 * <b>A missing observation costs nothing.</b> A {@code NaN} in a component of
 * {@code y} means that component was not observed. If the whole vector is
 * {@code NaN} the correction is skipped and only the prediction happens, so the
 * uncertainty simply grows until the next observation arrives; if some
 * components are present the correction runs through the rows of {@code H} and
 * {@code R} that were observed. Nothing is imputed and nothing is interpolated,
 * which is why irregular sampling and sensors of different rates are ordinary
 * here rather than a preprocessing problem.
 * <p>
 * There are two ways in. {@link #filter(StateSpaceModel, DMatrix)} runs a
 * whole series and hands back the entire history, which is what
 * {@link RtsSmoother} needs and what a caller plotting an interval wants; an
 * instance of this class steps one observation at a time and keeps only the
 * current state, which is what a live feed wants. The batch form is a recording
 * loop over the stepping form, so the two agree bit for bit.
 * <p>
 * <b>The covariance update is in Joseph form</b>,
 * <code>P &lt;- (I - K H) P (I - K H)' + K R K'</code>, and the result is
 * symmetrized at every step. The shorter <code>(I - K H) P</code> is
 * arithmetically the same and numerically is not: it is a difference of two
 * nearly equal matrices, and over a long series with accurate observations it
 * loses symmetry and then positive definiteness. Joseph is a sum of two
 * quadratic forms and cannot.
 * <p>
 * Instances are stateful and cannot be shared between threads. The static forms
 * hold no state and can.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Kalman_filter">Wikipedia Kalman
 * filter</a>.
 *
 * @since 1.5.3
 */
public final class KalmanFilter {

    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private final StateSpaceModel model;
    private final int d;
    private final int m;
    /**
     * How many time points the model is defined for.
     * {@link Integer#MAX_VALUE} if it does not run out, which makes every bound
     * check below one comparison that a constant model never fails.
     */
    private final int length;

    private final double[] mean;
    private DMatrix cov;
    private double logLikelihood;
    private int steps;
    private boolean corrected;

    private double[] lastInnovation;
    private DMatrix lastInnovationCovariance;
    private double lastSquaredMahalanobis;
    private int lastObserved;

    /**
     * A filter over the given model, positioned at the prior.
     *
     * @param model
     *            the model to filter with. Not modified
     * @throws IllegalArgumentException
     *             if {@code model} is {@code null}
     */
    public KalmanFilter(StateSpaceModel model) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        this.model = model;
        this.d = model.stateDimension();
        this.m = model.observationDimension();
        this.length = model.length();
        this.mean = new double[d];
        reset();
    }

    /**
     * The model this filter runs.
     *
     * @return the model handed to the constructor
     */
    public StateSpaceModel model() {
        return model;
    }

    /**
     * Returns the filter to the prior: the initial state of the model, no
     * observations, a log likelihood of zero.
     */
    public void reset() {
        model.initialMean(mean);
        cov = model.initialCovariance();
        logLikelihood = 0.0;
        steps = 0;
        corrected = false;
        lastInnovation = new double[0];
        lastInnovationCovariance = null;
        lastSquaredMahalanobis = Double.NaN;
        lastObserved = 0;
    }

    /**
     * Advances the state by one transition without observing anything:
     * <code>m &lt;- F m</code> and <code>P &lt;- F P F' + Q</code>.
     * <p>
     * Called twice in a row it is a two step ahead prediction, which is how a
     * forecast is produced: the mean stops moving once the transition has
     * settled and the covariance keeps growing, which is the widening interval.
     * <p>
     * The matrices used are {@code F(steps)} and {@code Q(steps)} -- the ones
     * leading into the time point about to be observed -- so a model with a
     * finite schedule can only be predicted forward while it lasts.
     *
     * @throws IllegalStateException
     *             if the model is defined for fewer time points than this would
     *             need
     */
    public void predict() {
        checkDefined();
        DMatrix f = model.transition(steps);
        double[] moved = new double[d];
        matVec(f, mean, moved);
        System.arraycopy(moved, 0, mean, 0, d);
        cov = symmetrize(f.mul(cov).mulBTrans(f).add(model.processNoise(steps)));
        corrected = false;
    }

    /**
     * Corrects the state against one observation and accumulates its
     * contribution to the log likelihood.
     * <p>
     * A {@code NaN} component of {@code y} was not observed. If every component
     * is {@code NaN} the state is left alone, the likelihood is unchanged, and
     * the step still counts.
     *
     * @param y
     *            the observation, of length
     *            {@link StateSpaceModel#observationDimension()}. A
     *            component may be {@code NaN} to mean it was not observed;
     *            anything else has to be finite. Not modified
     * @throws IllegalArgumentException
     *             if {@code y} is {@code null}, is of the wrong length, or
     *             holds an infinity
     * @throws IllegalStateException
     *             if the previous call was also an update -- two corrections
     *             against the same predicted state would count the same
     *             evidence twice -- or if the model is defined for fewer time
     *             points than have now been offered
     * @throws ArithmeticException
     *             if the innovation covariance at this step is not positive
     *             definite, which is a statement about the model and not about
     *             the data
     */
    public void update(double[] y) {
        if (y == null) {
            throw new IllegalArgumentException("y must not be null");
        }
        if (y.length != m) {
            throw new IllegalArgumentException("y is of length " + y.length + ", not " + m);
        }
        if (corrected) {
            throw new IllegalStateException(
                    "update() twice without predict() in between, at step " + steps);
        }
        checkDefined();
        int[] observed = observedIndices(y, steps);
        int k = observed.length;
        lastObserved = k;
        corrected = true;
        if (k == 0) {
            lastInnovation = new double[0];
            lastInnovationCovariance = null;
            lastSquaredMahalanobis = Double.NaN;
            ++steps;
            return;
        }

        DMatrix ht = selectRows(model.observation(steps), observed, k == m);
        DMatrix rt = selectBoth(model.observationNoise(steps), observed, k == m);

        // e = y_obs - H m
        double[] e = new double[k];
        for (int i = 0; i < k; ++i) {
            double sum = 0.0;
            for (int j = 0; j < d; ++j) {
                sum += ht.getUnsafe(i, j) * mean[j];
            }
            e[i] = y[observed[i]] - sum;
        }

        DMatrix pht = cov.mulBTrans(ht); // d by k
        DMatrix s = symmetrize(ht.mul(pht).add(rt)); // k by k
        DMatrix ls;
        try {
            ls = CholeskyDecomp.cholesky(s);
        } catch (RuntimeException ex) {
            throw new ArithmeticException("the innovation covariance at step " + steps
                    + " is not positive definite : " + ex.getMessage());
        }
        // cholesky() accepts a positive *semi*definite matrix, so a direction
        // with no uncertainty left gets through it and would fail four frames
        // down in the first substitution, naming neither the step nor the cause
        for (int i = 0; i < k; ++i) {
            if (!(ls.getUnsafe(i, i) > 0.0)) {
                throw new ArithmeticException("the innovation covariance at step " + steps
                        + " is singular: its factor has a zero at (" + i + ", " + i + ")");
            }
        }

        // the quadratic form and the determinant, both from the factor already
        // computed -- building a MultivariateNormal here would decompose again
        double[] z = new double[k];
        CholeskyDecomp.solveLower(ls, e, z);
        double quadratic = 0.0;
        for (int i = 0; i < k; ++i) {
            quadratic += z[i] * z[i];
        }
        logLikelihood += -0.5 * (k * LOG_TWO_PI + CholeskyDecomp.logDeterminant(ls) + quadratic);

        // K = (P H') S^-1, so K' = S^-1 (P H')' : one solve per row of P H'
        DMatrix gain = new DMatrix(d, k);
        double[] rhs = new double[k];
        double[] sol = new double[k];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < k; ++j) {
                rhs[j] = pht.getUnsafe(i, j);
            }
            CholeskyDecomp.solve(ls, rhs, sol);
            for (int j = 0; j < k; ++j) {
                gain.setUnsafe(i, j, sol[j]);
            }
        }

        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            for (int j = 0; j < k; ++j) {
                sum += gain.getUnsafe(i, j) * e[j];
            }
            mean[i] += sum;
        }

        // Joseph: P <- (I - K H) P (I - K H)' + K R K'
        DMatrix a = gain.mul(ht);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                a.setUnsafe(i, j, (i == j ? 1.0 : 0.0) - a.getUnsafe(i, j));
            }
        }
        cov = symmetrize(a.mul(cov).mulBTrans(a).add(gain.mul(rt).mulBTrans(gain)));

        lastInnovation = e;
        lastInnovationCovariance = s;
        lastSquaredMahalanobis = quadratic;
        ++steps;
    }

    /**
     * The mean of the current state.
     *
     * @param out
     *            where the mean is written, of length
     *            {@link StateSpaceModel#stateDimension()}. Its previous
     *            contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or of the wrong length
     */
    public void mean(double[] out) {
        if (out == null) {
            throw new IllegalArgumentException("out must not be null");
        }
        if (out.length != d) {
            throw new IllegalArgumentException("out is of length " + out.length + ", not " + d);
        }
        System.arraycopy(mean, 0, out, 0, d);
    }

    /**
     * The covariance of the current state.
     *
     * @return a copy, square of order
     *         {@link StateSpaceModel#stateDimension()}
     */
    public DMatrix covariance() {
        return cov.copy();
    }

    /**
     * The current state as the distribution it is.
     *
     * @return the normal distribution of the state given everything observed so
     *         far
     * @throws IllegalArgumentException
     *             if the covariance has become singular, which happens only for
     *             a model that leaves a direction of the state with no
     *             uncertainty at all
     */
    public MultivariateNormal state() {
        double[] copy = new double[d];
        System.arraycopy(mean, 0, copy, 0, d);
        return new MultivariateNormal(copy, cov);
    }

    /**
     * The innovation of the most recent update: what the observation was minus
     * what the model expected, over the components that were observed, in
     * increasing order of component.
     *
     * @return a copy, of length {@link #observedComponents()}. Empty before the
     *         first update and after one that observed nothing
     */
    public double[] innovation() {
        return lastInnovation.clone();
    }

    /**
     * The covariance of the most recent innovation,
     * <code>S = H P H' + R</code> over the observed components.
     *
     * @return a copy, square of order {@link #observedComponents()}, or
     *         {@code null} before the first update and after one that observed
     *         nothing
     */
    public DMatrix innovationCovariance() {
        return lastInnovationCovariance == null ? null : lastInnovationCovariance.copy();
    }

    /**
     * The squared Mahalanobis distance {@code e' S^-1 e} of the most recent
     * innovation.
     * <p>
     * This is the anomaly score. Under the model it is chi-squared on
     * {@link #observedComponents()} degrees of freedom, so a value far in that
     * tail says the observation does not fit what the model expected -- and
     * unlike a raw residual it is already scaled by how uncertain the
     * prediction was.
     *
     * @return the squared distance, or {@code NaN} before the first update and
     *         after one that observed nothing. {@code NaN} compares false
     *         against any threshold, so a missing observation raises no alarm
     */
    public double squaredMahalanobisDistance() {
        return lastSquaredMahalanobis;
    }

    /**
     * How many components the most recent update actually observed.
     *
     * @return a count between zero and
     *         {@link StateSpaceModel#observationDimension()}
     */
    public int observedComponents() {
        return lastObserved;
    }

    /**
     * The accumulated log likelihood of everything observed so far.
     * <p>
     * The prediction error decomposition: the sum of
     * <code>log N(e[t]; 0, S[t])</code> over the steps, which is the exact log
     * density of the observations under the model with the state integrated
     * out. Maximizing it over the entries of {@code F}, {@code Q} and {@code R}
     * is how the model is fitted.
     *
     * @return the log likelihood, zero at the prior
     */
    public double logLikelihood() {
        return logLikelihood;
    }

    /**
     * How many observations have been offered, including the ones that turned
     * out to be entirely missing.
     *
     * @return the number of calls to {@link #update(double[])} since the last
     *         {@link #reset()}
     */
    public int steps() {
        return steps;
    }

    /**
     * Runs the whole series and keeps everything.
     * <p>
     * The prior of the model is the predicted state at the first index: no
     * transition is applied before the first observation, so filtering an empty
     * series would return the prior unchanged.
     *
     * @param model
     *            the model. Not modified
     * @param y
     *            the observations, one per row, with as many columns as the
     *            model has observation coordinates. A {@code NaN} entry means
     *            that component was not observed at that time. Not modified
     * @return the whole filtered history
     * @throws IllegalArgumentException
     *             if either argument is {@code null}, if {@code y} does not
     *             have as many columns as the model expects, or if {@code y}
     *             holds an infinity
     * @throws ArithmeticException
     *             if an innovation covariance is not positive definite
     */
    public static Result filter(StateSpaceModel model, DMatrix y) {
        checkSeries(model, y);
        int n = y.numRows();
        int obs = model.observationDimension();
        KalmanFilter kf = new KalmanFilter(model);

        double[][] predictedMean = new double[n][];
        DMatrix[] predictedCovariance = new DMatrix[n];
        double[][] filteredMean = new double[n][];
        DMatrix[] filteredCovariance = new DMatrix[n];
        double[][] innovation = new double[n][];
        DMatrix[] innovationCovariance = new DMatrix[n];
        double[] squaredMahalanobis = new double[n];
        int[] observedComponents = new int[n];

        double[] row = new double[obs];
        for (int t = 0; t < n; ++t) {
            if (t > 0) {
                kf.predict();
            }
            predictedMean[t] = new double[model.stateDimension()];
            kf.mean(predictedMean[t]);
            predictedCovariance[t] = kf.covariance();
            for (int i = 0; i < obs; ++i) {
                row[i] = y.get(t, i);
            }
            kf.update(row);
            filteredMean[t] = new double[model.stateDimension()];
            kf.mean(filteredMean[t]);
            filteredCovariance[t] = kf.covariance();
            innovation[t] = kf.innovation();
            innovationCovariance[t] = kf.innovationCovariance();
            squaredMahalanobis[t] = kf.squaredMahalanobisDistance();
            observedComponents[t] = kf.observedComponents();
        }
        return new Result(predictedMean, predictedCovariance, filteredMean, filteredCovariance,
                innovation, innovationCovariance, squaredMahalanobis, observedComponents,
                kf.logLikelihood());
    }

    /**
     * The log likelihood of the series alone, without keeping the history.
     * <p>
     * The same number {@link #filter(StateSpaceModel, DMatrix)} reports,
     * bit for bit, and the form to call from inside an optimizer: it allocates
     * per step rather than per series.
     *
     * @param model
     *            the model. Not modified
     * @param y
     *            the observations, as {@link #filter(StateSpaceModel,
     *            DMatrix)} takes them. Not modified
     * @return the log likelihood of {@code y} under {@code model}
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #filter(StateSpaceModel, DMatrix)} states
     * @throws ArithmeticException
     *             if an innovation covariance is not positive definite
     */
    public static double logLikelihood(StateSpaceModel model, DMatrix y) {
        checkSeries(model, y);
        int n = y.numRows();
        int obs = model.observationDimension();
        KalmanFilter kf = new KalmanFilter(model);
        double[] row = new double[obs];
        for (int t = 0; t < n; ++t) {
            if (t > 0) {
                kf.predict();
            }
            for (int i = 0; i < obs; ++i) {
                row[i] = y.get(t, i);
            }
            kf.update(row);
        }
        return kf.logLikelihood();
    }

    /**
     * The whole filtered history of one series.
     * <p>
     * Index {@code t} runs over the rows of the observation matrix. The
     * predicted quantities are the state before the observation at {@code t}
     * was seen and the filtered ones are the state after; at {@code t == 0} the
     * predicted state is the prior of the model itself.
     * <p>
     * <b>The innovation is ragged on purpose.</b> Where components were
     * missing, {@link #innovation} holds only the ones that were observed, in
     * increasing order of component, and {@link #innovationCovariance} is of
     * that same reduced order -- and is {@code null} where nothing at all was
     * observed, since a matrix of order zero does not exist. Which components
     * those were is the pattern of {@code NaN} in the caller's own data.
     *
     * @since 1.5.3
     */
    public static final class Result {

        /** How many time points the series has. */
        public final int length;
        /** The mean of <code>x[t] | y[1..t-1]</code>, one row per time point. */
        public final double[][] predictedMean;
        /** The covariance of <code>x[t] | y[1..t-1]</code>. */
        public final DMatrix[] predictedCovariance;
        /** The mean of <code>x[t] | y[1..t]</code>, one row per time point. */
        public final double[][] filteredMean;
        /** The covariance of <code>x[t] | y[1..t]</code>. */
        public final DMatrix[] filteredCovariance;
        /** The innovation at each time point, over the observed components. */
        public final double[][] innovation;
        /** Its covariance, of the reduced order, {@code null} where nothing was observed. */
        public final DMatrix[] innovationCovariance;
        /**
         * The squared Mahalanobis distance of each innovation, the anomaly
         * score, {@code NaN} where nothing was observed.
         */
        public final double[] squaredMahalanobis;
        /** How many components were observed at each time point. */
        public final int[] observedComponents;
        /** The log likelihood of the whole series under the model. */
        public final double logLikelihood;

        Result(double[][] predictedMean, DMatrix[] predictedCovariance, double[][] filteredMean,
                DMatrix[] filteredCovariance, double[][] innovation, DMatrix[] innovationCovariance,
                double[] squaredMahalanobis, int[] observedComponents, double logLikelihood) {
            this.length = predictedMean.length;
            this.predictedMean = predictedMean;
            this.predictedCovariance = predictedCovariance;
            this.filteredMean = filteredMean;
            this.filteredCovariance = filteredCovariance;
            this.innovation = innovation;
            this.innovationCovariance = innovationCovariance;
            this.squaredMahalanobis = squaredMahalanobis;
            this.observedComponents = observedComponents;
            this.logLikelihood = logLikelihood;
        }

        /**
         * The filtered state at one time point, as the distribution it is.
         *
         * @param t
         *            the time point, from zero to {@link #length} minus one
         * @return the normal distribution of <code>x[t] | y[1..t]</code>
         * @throws IndexOutOfBoundsException
         *             if {@code t} is outside the series
         * @throws IllegalArgumentException
         *             if that covariance is singular
         */
        public MultivariateNormal filteredState(int t) {
            return new MultivariateNormal(filteredMean[t].clone(), filteredCovariance[t]);
        }

        /**
         * The predicted state at one time point, as the distribution it is.
         *
         * @param t
         *            the time point, from zero to {@link #length} minus one
         * @return the normal distribution of <code>x[t] | y[1..t-1]</code>
         * @throws IndexOutOfBoundsException
         *             if {@code t} is outside the series
         * @throws IllegalArgumentException
         *             if that covariance is singular
         */
        public MultivariateNormal predictedState(int t) {
            return new MultivariateNormal(predictedMean[t].clone(), predictedCovariance[t]);
        }
    }

    private static void checkSeries(StateSpaceModel model, DMatrix y) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y must not be null");
        }
        if (y.numColumns() != model.observationDimension()) {
            throw new IllegalArgumentException("y has " + y.numColumns() + " columns, not "
                    + model.observationDimension());
        }
        if (y.numRows() > model.length()) {
            throw new IllegalArgumentException("y has " + y.numRows()
                    + " time points and the model is defined for " + model.length());
        }
    }

    /**
     * Refuses a step the model has no matrices for. A constant model reports
     * {@link Integer#MAX_VALUE} and never fails this.
     */
    private void checkDefined() {
        if (steps >= length) {
            throw new IllegalStateException("the model is defined for " + length
                    + " time points and this would be step " + steps);
        }
    }

    /**
     * Which components of {@code y} were observed. A {@code NaN} was not; an
     * infinity is neither an observation nor a way of saying there was none.
     */
    private static int[] observedIndices(double[] y, int step) {
        int k = 0;
        for (int i = 0; i < y.length; ++i) {
            if (Double.isInfinite(y[i])) {
                throw new IllegalArgumentException(
                        "y[" + i + "] at step " + step + " is infinite : " + y[i]);
            }
            if (!Double.isNaN(y[i])) {
                ++k;
            }
        }
        int[] out = new int[k];
        int at = 0;
        for (int i = 0; i < y.length; ++i) {
            if (!Double.isNaN(y[i])) {
                out[at++] = i;
            }
        }
        return out;
    }

    private static DMatrix selectRows(DMatrix A, int[] rows, boolean all) {
        if (all) {
            return A;
        }
        DMatrix out = new DMatrix(rows.length, A.numColumns());
        for (int i = 0; i < rows.length; ++i) {
            for (int j = 0; j < A.numColumns(); ++j) {
                out.setUnsafe(i, j, A.getUnsafe(rows[i], j));
            }
        }
        return out;
    }

    private static DMatrix selectBoth(DMatrix A, int[] idx, boolean all) {
        if (all) {
            return A;
        }
        DMatrix out = new DMatrix(idx.length, idx.length);
        for (int i = 0; i < idx.length; ++i) {
            for (int j = 0; j < idx.length; ++j) {
                out.setUnsafe(i, j, A.getUnsafe(idx[i], idx[j]));
            }
        }
        return out;
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

    private static void matVec(DMatrix A, double[] x, double[] out) {
        int rows = A.numRows();
        int cols = A.numColumns();
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += A.getUnsafe(i, j) * x[j];
            }
            out[i] = sum;
        }
    }
}
