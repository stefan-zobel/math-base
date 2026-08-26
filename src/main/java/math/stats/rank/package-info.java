/**
 * The rank statistics: the ranking itself, the coefficients that take more
 * than a rank sum to compute, and the null distributions of all of them.
 * <p>
 * A rank test throws the values away and keeps only their order, which is what
 * makes it available where the t-test's normality is not: its null
 * distribution is a count of equally likely orderings and holds for every
 * continuous distribution the sample might have come from.
 * <p>
 * {@link math.stats.rank.Ranks} is the ranking, with tied values sharing the
 * average of the ranks they span. Its {@code tieSum} is what all the variance
 * corrections here take, and its {@code hasTies} is what decides between an
 * exact null and its normal approximation -- exactness is a count of
 * orderings, and a tie means there is no ordering to count.
 * <p>
 * {@link math.stats.rank.MannWhitneyU} carries the null of the two-sample
 * statistic and {@link math.stats.rank.WilcoxonSignedRank} that of the
 * one-sample one, each exactly below a measured size limit and by the normal
 * approximation above it. Kruskal-Wallis has no class here because it needs
 * none: its statistic is chi-squared on {@code k - 1} degrees of freedom, and
 * {@code math.distribution.ChiSquare} is that distribution.
 * <p>
 * {@link math.stats.rank.SpearmanRho} and {@link math.stats.rank.KendallTau}
 * are the two rank correlations, and they carry their coefficients as well as
 * their nulls -- unlike the rank sums, which are three lines once the ranking
 * is in hand, {@code tau} takes a merge sort that counts its own exchanges.
 * Their exact nulls are also the two extremes of what that costs: Kendall
 * counts orderings by a recursion and reaches {@code n = 200}, while Spearman
 * has no recursion to be had and enumerates, which stops at {@code n = 9}.
 * <p>
 * The tests themselves are in {@link math.stats.HypothesisTests}, the way the
 * goodness of fit tests are. This package holds the distributions and what
 * they are computed from.
 */
package math.stats.rank;
