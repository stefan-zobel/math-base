/**
 * Maximum likelihood estimation of the parameters of a distribution from a
 * sample.
 * <p>
 * {@link math.stats.mle.MLE} fits nine families and returns each one as its own
 * small value type -- {@link math.stats.mle.ParNormal},
 * {@link math.stats.mle.ParGamma} and the rest -- carrying public final fields
 * named after the parameters rather than getters. Every one of them implements
 * {@link math.stats.Validity}, so a fit that produced something unusable says so
 * through {@code isValid()} instead of being found out later.
 * <p>
 * The nine divide into three kinds of work, which is worth knowing before
 * calling one in a loop. The normal, log-normal, exponential and chi-squared
 * fits are closed forms and cost one pass over the sample. The gamma and
 * Weibull shapes are the root of a one-dimensional equation in digamma
 * respectively in logarithms, bracketed and then solved with
 * {@code math.solve.RootFinder.brentDekker}. The beta, Student t and Cauchy
 * fits iterate -- Newton on a two by two system for the beta, a reweighting for
 * the other two -- and those three, and only those three, carry a
 * {@code converged} flag beside their parameters. <b>An unconverged fit is
 * still returned</b>, because its numbers are usually close and throwing away a
 * near miss helps nobody; the flag is there so that a caller who cares can
 * ask.
 * <p>
 * <b>The normal fit divides by {@code n}, not by {@code n - 1}.</b> That is
 * what maximum likelihood means and it is the right thing for the use these
 * estimates are put to here, but it is not what the textbooks tabulate: the
 * published critical values for a goodness of fit test against a fitted normal
 * assume the {@code n - 1} standard deviation, and reading a maximum likelihood
 * fit off those tables rejects too often. {@code math.stats.gof.Lilliefors}
 * carries the measurement. Nothing goes wrong as long as the estimator that
 * produced a statistic is the estimator its null distribution was drawn with,
 * which is exactly why the fitted goodness of fit tests in
 * {@code math.stats.gof} estimate through this package rather than around it.
 */
package math.stats.mle;
