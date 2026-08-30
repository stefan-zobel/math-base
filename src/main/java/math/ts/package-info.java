/**
 * Time series models.
 * <p>
 * Not {@link math.probe}, whose comment is "Sample metrics" and whose
 * {@code ACF} and {@code PartialACF} describe the data in front of them and
 * assume nothing about where it came from. What is here assumes a model, and in
 * exchange answers questions the data alone cannot: where the level really was
 * at a moment nobody measured it accurately, how wide the interval around that
 * is, and how surprising the observation that just arrived was.
 * <p>
 * The model is the linear Gaussian state space model,
 * {@link math.ts.StateSpaceModel}: a hidden state that moves by a linear map
 * plus noise, and an observation that is a linear function of the state plus
 * more noise. It is a small model that a large number of familiar things turn
 * out to be -- a random walk plus noise, which is exponential smoothing; a
 * local linear trend; a seasonal decomposition; an ARMA model; a regression
 * whose coefficients drift; a vehicle whose position is measured badly and
 * often. {@link math.ts.LinearGaussianModel} is the case where the four
 * matrices do not change, {@link math.ts.TimeVaryingModel} the case where they
 * do.
 * <p>
 * <b>The smallest useful model, in full.</b> A quantity that drifts on its own,
 * read by an instrument that is off by about one unit, with the readings one
 * per row in {@code y}:
 *
 * <pre>
 * DMatrix keep  = new DMatrix(1, 1).set(0, 0, 1.0);   // the level stays where it was
 * DMatrix drift = new DMatrix(1, 1).set(0, 0, 0.04);  // give or take this much variance
 * DMatrix error = new DMatrix(1, 1).set(0, 0, 1.0);   // and the reading is off by this
 * DMatrix vague = new DMatrix(1, 1).set(0, 0, 100.0); // knowing nothing to begin with
 *
 * LinearGaussianModel model = new LinearGaussianModel(keep, drift, keep, error,
 *         new double[] { y.get(0, 0) }, vague);
 *
 * RtsSmoother.Result cleaned = RtsSmoother.smooth(model, y);
 * // cleaned.mean[t][0]              where the quantity really was at time t
 * // cleaned.covariance[t].get(0, 0) how sure that is
 * </pre>
 *
 * That is a random walk plus noise, and it is exponential smoothing with the
 * constant worked out rather than guessed. Only the <em>ratio</em> of the two
 * variances matters: {@code 0.04} against {@code 1.0} settles at a gain of
 * about {@code 0.18}, so eighteen percent of each new reading is taken as a
 * real move and the rest as error. Turn the ratio down and the answer
 * approaches a straight line through the data; turn it up and it follows every
 * reading. Leave a {@code Double#NaN} in {@code y} and that day is simply
 * skipped. Everything else in this package is this idea with more components in
 * the state, and {@link math.ts.KalmanFilter#logLikelihood} is how the ratio
 * stops being a guess at all.
 * <p>
 * <b>Which of those two is a question about the model, not about the data.</b>
 * A sensor that drops out, a month with no reading, a series sampled at
 * whatever times someone happened to record -- none of that needs varying
 * matrices, because a missing observation is a {@code NaN} and costs nothing.
 * What needs them is a step whose <em>length</em> changes, so that the process
 * noise it carries changes with it, or a measurement setup that is genuinely
 * different later than it was earlier. Reaching for a schedule when a
 * {@code NaN} was meant produces a model that is harder to write and no more
 * general.
 * <p>
 * What it buys is that <b>everything stays Gaussian</b>. The posterior
 * {@code p(x[t] | y[1..t])} that {@link math.ts.KalmanFilter} computes is a
 * {@link math.distribution.MultivariateNormal}, exactly, with no approximation
 * and no sampler anywhere -- which is what distinguishes this from
 * {@code math.stats.bayes}, where quadrature and a Gaussian fit at the mode are
 * necessary because the posterior of a general model is not any known
 * distribution. Here it is, and the recursion is what evaluates it.
 * <p>
 * Three things follow from that and are worth naming, because they are what
 * makes the model worth reaching for over a moving average:
 * <ul>
 * <li><b>A missing observation costs nothing.</b> A {@code NaN} means a
 * component was not observed; the correction step for it is skipped and the
 * uncertainty simply keeps growing until something arrives. Nothing is imputed.
 * That is what makes irregular sampling, mixed frequencies and a sensor that
 * drops out ordinary rather than a preprocessing problem.</li>
 * <li><b>The innovation is an anomaly score</b>, not a residual. It comes with
 * its own covariance, so {@code e' S^-1 e} is already scaled by how uncertain
 * the prediction was and is chi-squared on the number of components observed --
 * a number that can be compared against a threshold rather than eyeballed.</li>
 * <li><b>The sum of the innovation densities is the exact log likelihood</b> of
 * the series with the state integrated out. That is what a fit of the model
 * maximizes, and it is why an ARMA estimator is this package plus an
 * optimizer.</li>
 * </ul>
 * <p>
 * {@link math.ts.RtsSmoother} is the second reading of the same series. A
 * filter may only use the past, which is right for a live feed and artificial
 * once the series is over: what happened at the end says something about where
 * the level was at the beginning. The backward pass is that revision, and like
 * the forward one it is exact rather than a smoothing heuristic.
 */
package math.ts;
