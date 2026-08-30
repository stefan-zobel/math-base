/**
 * Classical hypothesis tests, and the small vocabulary they answer in.
 * <p>
 * {@link math.stats.HypothesisTests} fronts the package the way
 * {@code math.solve.Quadrature} fronts the integrators: the t family, the F
 * variance ratio, chi-squared goodness of fit and independence, Fisher's exact
 * test, one- and two-sample Kolmogorov-Smirnov, Cramer-von Mises and
 * Anderson-Darling each in a fully specified and a fitted form, Lilliefors,
 * Durbin-Watson, the rank tests -- Mann-Whitney U, Wilcoxon signed rank and
 * Kruskal-Wallis -- for the samples the t family cannot be asked about, the
 * three tests of zero correlation, Pearson beside Spearman and Kendall, and
 * the parametric k-sample family: one-way ANOVA, Welch's ANOVA for groups
 * that are not equally variable, and Bartlett, Levene and Brown-Forsythe for
 * asking whether they are.
 * <p>
 * <b>Three groups are not three pairs of groups.</b> Running the t-test on
 * every pair is the mistake this package is built not to make convenient:
 * {@link math.stats.HypothesisTests#oneWayAnova(double[][])} and
 * {@link math.stats.HypothesisTests#kruskalWallis(double[][])} ask the
 * question once, and a rejection is followed by pairwise comparisons
 * corrected through {@link math.stats.MultipleTesting}.
 * <p>
 * {@link math.stats.MultipleTesting} answers the question none of those can
 * answer on its own: twenty tests at five percent are wrong once on average,
 * which is not a defect of any one of them. Given the whole family of
 * p-values it returns adjusted ones -- controlling the false discovery rate
 * where there are many tests, or the family-wise error rate where there are
 * few and none of them may be wrong.
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
 * ranking and the two rank correlations they are computed from, and
 * {@code math.stats.mle} the maximum likelihood fitting that the fitted
 * tests estimate with.
 * {@link math.stats.Validity} is the small interface those fits report through.
 */
package math.stats;
