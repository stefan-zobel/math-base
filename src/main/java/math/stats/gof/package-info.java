/**
 * Null distributions of the goodness of fit statistics of an ordered sample
 * of {@code n} independent uniforms {@code U_i} over {@code (0,1)}, which is
 * what a continuous sample becomes once its own distribution function is
 * applied to it.
 * <p>
 * {@link math.stats.gof.KolmogorovSmirnov} carries the two-sided statistic
 * {@code D_n}, {@link math.stats.gof.KolmogorovSmirnovPlus} the one-sided
 * {@code D_n^+} -- which by symmetry serves {@code D_n^-} as well -- and
 * {@link math.stats.gof.AndersonDarling} the statistic {@code A_n^2}, which
 * weights the tails, where the Kolmogorov-Smirnov statistic is least
 * sensitive.
 * <p>
 * Three of them are not distributions of a uniform sample at all, because
 * the situations they describe have no such reduction:
 * {@link math.stats.gof.KolmogorovSmirnovTwoSample} conditions on the pooled
 * order of two samples, {@link math.stats.gof.Lilliefors} draws the null of a
 * statistic measured against a distribution fitted to the same sample, and
 * {@link math.stats.gof.DurbinWatson} inverts the characteristic function of
 * a ratio of quadratic forms whose eigenvalues come from a design matrix.
 * <p>
 * The tests themselves are in {@link math.stats.HypothesisTests}. This
 * package holds the distributions and nothing else.
 */
package math.stats.gof;
