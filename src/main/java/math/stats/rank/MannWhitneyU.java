package math.stats.rank;

import math.cern.ProbabilityFuncs;

/**
 * The null distribution of the Mann-Whitney statistic
 * {@code U = R_x - m (m + 1) / 2}, the number of pairs {@code (x_i, y_j)} in
 * which the first observation is the larger, when both samples come from the
 * same distribution.
 * <p>
 * https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
 *
 * @since 1.5.3
 */
public final class MannWhitneyU {

    /**
     * The largest product {@code m * n} the exact method will accept.
     * <p>
     * The recursion walks the {@code m + n} ranks once and carries a table of
     * {@code min(m, n) + 1} rows, so the cost turns on the shape and not only
     * on the product: at a product of 2500 it is measured at 11 ms for
     * {@code m = n = 50}, the worst shape there, and at 1.7 ms for
     * {@code m = 1}. Doubling the product to 5000 would take the worst shape
     * past 40 ms.
     */
    public static final long EXACT_LIMIT = 2500L;

    /**
     * The exact upper tail {@code P[U >= u]} under the null hypothesis.
     * <p>
     * Every one of the {@code C(m + n, m)} ways of dealing the ranks
     * {@code 1 .. m + n} out to the two samples is equally likely, so the tail
     * is the fraction of them whose rank sum reaches {@code u}. Those are
     * counted by walking the ranks once and carrying, for every count of
     * ranks taken so far, how many ways each attainable sum can be reached.
     * The recursion adds and never subtracts, so nothing cancels, and the
     * counts are turned into probabilities by their own total -- which is
     * {@code C(m + n, m)} and stays an ordinary number for every shape
     * {@link #EXACT_LIMIT} admits.
     * <p>
     * The statistic is a whole number under the null, so a {@code u} between
     * two attainable values is answered for the next one up.
     * <p>
     * The lower tail needs no method of its own: the null is symmetric about
     * {@code m n / 2}, so {@code P[U <= u]} is this function at
     * {@code m n - u}.
     *
     * @param m
     *            the size of the first sample, at least one
     * @param n
     *            the size of the second sample, at least one
     * @param u
     *            the value the statistic must reach
     * @return the probability that {@code U} is at least {@code u}, in
     *         {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if either size is smaller than one, if {@code u} is
     *             {@code NaN}, or if {@code m * n} is above
     *             {@link #EXACT_LIMIT}
     */
    public static double barFExact(int m, int n, double u) {
        requireSizes(m, n, u);
        if ((long) m * n > EXACT_LIMIT) {
            throw new IllegalArgumentException(
                    "m * n is " + (long) m * n + ", above the exact limit of " + EXACT_LIMIT);
        }
        double mn = (double) m * n;
        if (u <= 0.0) {
            return 1.0;
        }
        if (u > mn) {
            return 0.0;
        }
        // the distribution of U depends on the two sizes only as an unordered
        // pair, and the smaller one is both the cheaper table and the shorter
        // walk
        double[] pmf = pmf(Math.min(m, n), Math.max(m, n));
        int from = (int) Math.ceil(u);
        double tail = 0.0;
        for (int k = pmf.length - 1; k >= from; k--) {
            tail += pmf[k];
        }
        return Math.min(1.0, tail);
    }

    /**
     * The normal approximation to the upper tail {@code P[U >= u]}, which is
     * what serves when the samples are large or the ranking has ties.
     * <p>
     * {@code U} has mean {@code m n / 2} and variance
     * {@code m n ((N + 1) - tieSum / (N (N - 1))) / 12} with {@code N = m + n},
     * and the tail is read off a normal with a continuity correction of one
     * half. The correction is not decoration: measured against the exact tail
     * over the p-values between 0.001 and 0.2 it halves the error at every
     * shape tried, from 4.7e-2 to 1.5e-2 at {@code m = n = 4} and from 1.8e-3
     * to 8.2e-4 at {@code m = n = 50}. It is the deep tail where it is the
     * worse of the two in relative terms, and the deep tail of a small sample
     * is what {@link #barFExact} is for.
     * <p>
     * A tie sum of {@code N^3 - N} says every observation in both samples is
     * the same value. The variance is then exactly zero, the ranking carries
     * no information at all, and the tail is {@code 1}.
     *
     * @param m
     *            the size of the first sample, at least one
     * @param n
     *            the size of the second sample, at least one
     * @param u
     *            the value the statistic must reach
     * @param tieSum
     *            the sum of {@code t^3 - t} over the tie runs of the pooled
     *            ranking, which is {@link Ranks.Result#tieSum}, or {@code 0}
     *            when there are no ties
     * @return the approximate probability that {@code U} is at least
     *         {@code u}, in {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if either size is smaller than one, if {@code u} is
     *             {@code NaN}, or if {@code tieSum} is negative, not finite,
     *             or larger than {@code N^3 - N}
     */
    public static double barFAsymptotic(int m, int n, double u, double tieSum) {
        requireSizes(m, n, u);
        double total = (double) m + n;
        double most = total * total * total - total;
        if (!(tieSum >= 0.0) || tieSum > most) {
            throw new IllegalArgumentException("tieSum must lie in [0, " + most + "] : " + tieSum);
        }
        double mn = (double) m * n;
        double variance = mn / 12.0 * ((total + 1.0) - tieSum / (total * (total - 1.0)));
        if (!(variance > 0.0)) {
            return 1.0;
        }
        double z = (u - mn / 2.0 - 0.5) / Math.sqrt(variance);
        return Math.min(1.0, Math.max(0.0, ProbabilityFuncs.normal(-z)));
    }

    /**
     * {@code P[U = u]} for {@code u = 0 .. m n}, with {@code m} the smaller of
     * the two sizes. {@code counts[j][s]} is the number of ways to have taken
     * {@code j} of the ranks seen so far with rank sum {@code s}.
     */
    private static double[] pmf(int m, int n) {
        int ranks = m + n;
        int offset = m * (m + 1) / 2;
        int maxSum = offset + m * n;
        double[][] counts = new double[m + 1][maxSum + 1];
        counts[0][0] = 1.0;
        for (int rank = 1; rank <= ranks; rank++) {
            // downwards in j, so that every row is read before it is written
            for (int j = Math.min(rank, m); j >= 1; j--) {
                double[] taken = counts[j];
                double[] fewer = counts[j - 1];
                for (int s = maxSum; s >= rank; s--) {
                    double ways = fewer[s - rank];
                    if (ways != 0.0) {
                        taken[s] += ways;
                    }
                }
            }
        }
        double[] pmf = new double[m * n + 1];
        double total = 0.0;
        for (int u = 0; u < pmf.length; u++) {
            pmf[u] = counts[m][offset + u];
            total += pmf[u];
        }
        for (int u = 0; u < pmf.length; u++) {
            pmf[u] /= total;
        }
        return pmf;
    }

    private static void requireSizes(int m, int n, double u) {
        if (m < 1) {
            throw new IllegalArgumentException("m must be at least one : " + m);
        }
        if (n < 1) {
            throw new IllegalArgumentException("n must be at least one : " + n);
        }
        if (Double.isNaN(u)) {
            throw new IllegalArgumentException("u must not be NaN");
        }
    }

    private MannWhitneyU() {
        throw new AssertionError();
    }
}
