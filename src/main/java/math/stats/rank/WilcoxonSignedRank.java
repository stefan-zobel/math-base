package math.stats.rank;

import math.cern.ProbabilityFuncs;

/**
 * The null distribution of the Wilcoxon statistic {@code W+}, the sum of the
 * ranks of the positive differences when the differences are ranked by
 * magnitude and each of them is as likely to be positive as negative.
 * <p>
 * https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
 *
 * @since 1.5.3
 */
public final class WilcoxonSignedRank {

    /**
     * The largest sample size the exact method will accept.
     * <p>
     * The recursion is one pass per rank over a table of {@code n (n + 1) / 2}
     * entries, so it grows as the cube: measured at 0.06 ms for {@code n = 50},
     * 2.1 ms for {@code n = 200} and 7.5 ms for {@code n = 300}. At
     * {@code n = 800} it is already 138 ms, and by then the approximation is
     * accurate to better than {@code 4e-4}.
     */
    public static final int EXACT_LIMIT = 300;

    /**
     * The exact upper tail {@code P[W+ >= w]} under the null hypothesis.
     * <p>
     * Under the null each of the {@code n} ranks joins the sum or does not,
     * independently and with probability one half, so the distribution is
     * built one rank at a time as
     * {@code p(s) <- (p(s) + p(s - i)) / 2}. Carrying probabilities rather
     * than the counts of the {@code 2^n} sign patterns keeps every number in
     * {@code [0, 1]}: the counts leave the exact integers of a {@code double}
     * at {@code n = 54} and its range entirely at {@code n = 1024}.
     * <p>
     * The statistic is a whole number under the null, so a {@code w} between
     * two attainable values is answered for the next one up.
     * <p>
     * The lower tail needs no method of its own: the null is symmetric about
     * {@code n (n + 1) / 4}, so {@code P[W+ <= w]} is this function at
     * {@code n (n + 1) / 2 - w}.
     *
     * @param n
     *            the number of differences, at least one and at most
     *            {@link #EXACT_LIMIT}
     * @param w
     *            the value the statistic must reach
     * @return the probability that {@code W+} is at least {@code w}, in
     *         {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code n} is smaller than one or above
     *             {@link #EXACT_LIMIT}, or if {@code w} is {@code NaN}
     */
    public static double barFExact(int n, double w) {
        requireSize(n, w);
        if (n > EXACT_LIMIT) {
            throw new IllegalArgumentException("n is " + n + ", above the exact limit of " + EXACT_LIMIT);
        }
        double most = n * (n + 1.0) / 2.0;
        if (w <= 0.0) {
            return 1.0;
        }
        if (w > most) {
            return 0.0;
        }
        double[] pmf = pmf(n);
        int from = (int) Math.ceil(w);
        double tail = 0.0;
        for (int k = pmf.length - 1; k >= from; k--) {
            tail += pmf[k];
        }
        return Math.min(1.0, tail);
    }

    /**
     * The normal approximation to the upper tail {@code P[W+ >= w]}, which is
     * what serves when the sample is large or the magnitudes have ties.
     * <p>
     * {@code W+} has mean {@code n (n + 1) / 4} and variance
     * {@code n (n + 1) (2 n + 1) / 24 - tieSum / 48}, and the tail is read off
     * a normal with a continuity correction of one half. Measured against the
     * exact tail over the p-values between 0.001 and 0.2 the correction halves
     * the error at every size tried, from 4.4e-2 to 1.6e-2 at {@code n = 5}
     * and from 3.0e-3 to 1.6e-3 at {@code n = 50}.
     * <p>
     * Unlike the two-sample statistic this one keeps a variance whatever the
     * ties do: even if every magnitude is the same the signs are still free,
     * and the correction can take the variance down to
     * {@code n (n + 1)^2 / 16} but no further.
     *
     * @param n
     *            the number of differences, at least one
     * @param w
     *            the value the statistic must reach
     * @param tieSum
     *            the sum of {@code t^3 - t} over the tie runs of the ranked
     *            magnitudes, which is {@link Ranks.Result#tieSum}, or
     *            {@code 0} when there are no ties
     * @return the approximate probability that {@code W+} is at least
     *         {@code w}, in {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code n} is smaller than one, if {@code w} is
     *             {@code NaN}, or if {@code tieSum} is negative, not finite,
     *             or larger than {@code n^3 - n}
     */
    public static double barFAsymptotic(int n, double w, double tieSum) {
        requireSize(n, w);
        double nd = n;
        double most = nd * nd * nd - nd;
        if (!(tieSum >= 0.0) || tieSum > most) {
            throw new IllegalArgumentException("tieSum must lie in [0, " + most + "] : " + tieSum);
        }
        double mean = nd * (nd + 1.0) / 4.0;
        double variance = nd * (nd + 1.0) * (2.0 * nd + 1.0) / 24.0 - tieSum / 48.0;
        double z = (w - mean - 0.5) / Math.sqrt(variance);
        return Math.min(1.0, Math.max(0.0, ProbabilityFuncs.normal(-z)));
    }

    /**
     * {@code P[W+ = w]} for {@code w = 0 .. n (n + 1) / 2}, one rank at a
     * time, in probability form throughout.
     */
    private static double[] pmf(int n) {
        int most = n * (n + 1) / 2;
        double[] pmf = new double[most + 1];
        pmf[0] = 1.0;
        for (int i = 1; i <= n; i++) {
            // downwards, so that every entry is read before it is written
            for (int s = most; s >= i; s--) {
                pmf[s] = 0.5 * (pmf[s] + pmf[s - i]);
            }
            // and the sums this rank cannot reach only lose the half that
            // takes it
            for (int s = Math.min(i - 1, most); s >= 0; s--) {
                pmf[s] = 0.5 * pmf[s];
            }
        }
        return pmf;
    }

    private static void requireSize(int n, double w) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be at least one : " + n);
        }
        if (Double.isNaN(w)) {
            throw new IllegalArgumentException("w must not be NaN");
        }
    }

    private WilcoxonSignedRank() {
        throw new AssertionError();
    }
}
