/**
 * Classical hypothesis tests, and the small vocabulary they answer in.
 * <p>
 * {@link math.stats.HypothesisTests} fronts the package the way
 * {@code math.solve.Quadrature} fronts the integrators: the t family, the F
 * variance ratio, chi-squared goodness of fit and independence, Fisher's exact
 * test, one- and two-sample Kolmogorov-Smirnov, Cramer-von Mises and
 * Anderson-Darling each in a fully specified and a fitted form, Lilliefors,
 * Durbin-Watson, and the rank tests -- Mann-Whitney U, Wilcoxon signed rank
 * and Kruskal-Wallis -- for the samples the t family cannot be asked about.
 * <p>
 * {@link math.stats.MultipleTesting} answers the question none of those can
 * answer on its own: twenty tests at five percent are wrong once on average,
 * which is not a defect of any one of them. Given the whole family of
 * p-values it returns adjusted ones that control the false discovery rate.
 * <p>
 * Three rules shape everything here.
 * <p>
 * <b>Every test returns a {@link math.stats.TestResult}</b>, or a type holding
 * one where the test yields more than a statistic and a p-value:
 * {@link math.stats.TTestResult} adds a confidence interval,
 * {@link math.stats.FTestResult} a second degrees of freedom, and
 * {@link math.stats.SimulatedTestResult} says how much of a drawn p-value is
 * simulation noise. A wide type whose fields were {@link java.lang.Double#NaN}
 * for most tests would say less, not more.
 * <p>
 * <b>A test that cannot be computed throws</b>, naming what was wrong with the
 * sample. None of them returns {@code NaN} to mean "no answer", because a
 * caller who does not check gets a number that spreads rather than a failure
 * that stops.
 * <p>
 * <b>The significance level belongs to the reader</b>, not to the computation,
 * so it is an argument of {@link math.stats.TestResult#rejectsAt(double)} and
 * can be asked as often as one likes without running the test again. Choosing
 * it after seeing the p-value is not a test, and no signature here makes that
 * convenient.
 * <p>
 * The three neighbouring packages divide the work: {@code math.stats.gof}
 * holds the null distributions of the goodness of fit statistics,
 * {@code math.stats.rank} those of the rank statistics together with the
 * ranking they are computed from, and {@code math.stats.mle} the maximum
 * likelihood fitting that the fitted tests estimate with.
 * {@link math.stats.Validity} is the small interface those fits report through.
 */
package math.stats;
