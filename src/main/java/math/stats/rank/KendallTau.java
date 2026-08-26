package math.stats.rank;

import math.cern.ProbabilityFuncs;

/**
 * Kendall's {@code tau-b}, the excess of concordant over discordant pairs
 * scaled by what the ties leave attainable, and its null distribution when the
 * two orderings are paired at random.
 * <p>
 * https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient
 *
 * @since 1.5.3
 */
public final class KendallTau {

    /**
     * The largest sample the exact method will accept.
     * <p>
     * The recursion is one pass per observation over a table of
     * {@code n (n - 1) / 2} entries, so it grows as the cube: measured at 1.5 ms
     * for {@code n = 100}, 11 ms for {@code n = 200} and 39 ms for
     * {@code n = 300}.
     */
    public static final int EXACT_LIMIT = 200;

    /**
     * The coefficient, the concordance excess it is scaled from, and what the
     * ties leave of the null variance.
     * <p>
     * The variance travels in the result rather than as an argument of
     * {@link #barFAsymptotic(Result, double)} because correcting it for ties
     * takes six sums over the runs of equal values in the two samples, and
     * carrying six numbers to the place they are used would be worse than
     * carrying the one they produce.
     */
    public static final class Result {

        /**
         * {@code tau-b}, in {@code [-1, 1]}: {@link #s} over the geometric
         * mean of the pairs each sample leaves untied. It reaches {@code 1}
         * only if the two orderings agree and neither has a tie.
         */
        public final double tau;

        /**
         * The number of concordant pairs less the number of discordant ones. A
         * whole number, and the statistic the null distribution is over --
         * {@link #tau} divides it by a quantity that depends on the ties, which
         * would put the atoms of that distribution in different places for
         * every sample.
         */
        public final double s;

        /**
         * The variance of {@link #s} under the null hypothesis, corrected for
         * the ties in both samples. Without ties it is
         * {@code n (n - 1) (2 n + 5) / 18}.
         */
        public final double nullVariance;

        /**
         * Whether either sample repeats a value. Without ties {@link #s} moves
         * in steps of two and the exact null is available; with them it moves
         * in steps of one and it is not.
         */
        public final boolean hasTies;

        Result(double tau, double s, double nullVariance, boolean hasTies) {
            this.tau = tau;
            this.s = s;
            this.nullVariance = nullVariance;
            this.hasTies = hasTies;
        }
    }

    /**
     * Kendall's {@code tau-b} for two paired samples, in {@code O(n log n)}.
     * <p>
     * The pairs are sorted by the first sample and the second is then sorted by
     * a merge sort that counts the exchanges it makes; those exchanges are the
     * discordant pairs, which is Knight's observation and the whole reason this
     * is not the obvious quadratic loop over pairs.
     *
     * @param x
     *            the first sample, at least two observations, none of them
     *            {@code NaN}, not all equal
     * @param y
     *            the second sample, as many observations as {@code x}, under
     *            the same conditions
     * @return the coefficient, the concordance excess and the null variance
     * @throws IllegalArgumentException
     *             if either sample is {@code null}, holds fewer than two
     *             observations or a {@code NaN}, if the two are of different
     *             lengths, or if either takes only one value
     */
    public static Result of(double[] x, double[] y) {
        if (x == null || y == null) {
            throw new IllegalArgumentException("neither sample may be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException(
                    "x and y must be paired, but their lengths are " + x.length + " and " + y.length);
        }
        int n = x.length;
        if (n < 2) {
            throw new IllegalArgumentException("a correlation needs at least two observations, got " + n);
        }
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(x[i]) || Double.isNaN(y[i])) {
                throw new IllegalArgumentException("observation " + i + " is NaN, which has no place in an ordering");
            }
        }

        double[] xs = x.clone();
        double[] ys = y.clone();
        double[] bufferA = new double[n];
        double[] bufferB = new double[n];
        sortPairs(xs, ys, 0, n, bufferA, bufferB);

        double total = (double) n * (n - 1) / 2.0;
        // the runs of equal values in the first sample, and the runs in which
        // both samples repeat, which the pairs are sorted for
        double[] first = tieAggregates(xs, null);
        double[] both = tieAggregates(xs, ys);

        long discordant = sortAndCount(ys, 0, n, bufferA);
        double[] second = tieAggregates(ys, null);

        double s = total - first[0] - second[0] + both[0] - 2.0 * discordant;
        double untiedFirst = total - first[0];
        double untiedSecond = total - second[0];
        if (!(untiedFirst > 0.0)) {
            throw new IllegalArgumentException("the first sample takes one value only, so it leaves no pair "
                    + "to be concordant or discordant");
        }
        if (!(untiedSecond > 0.0)) {
            throw new IllegalArgumentException("the second sample takes one value only, so it leaves no pair "
                    + "to be concordant or discordant");
        }
        double tau = Math.max(-1.0, Math.min(1.0, s / Math.sqrt(untiedFirst * untiedSecond)));

        double nd = n;
        double v0 = nd * (nd - 1.0) * (2.0 * nd + 5.0);
        // the aggregates are held as sums of t (t - 1) / 2, and the two
        // correction terms want sums of t (t - 1)
        double tiedPairs = (2.0 * first[0]) * (2.0 * second[0]) / (2.0 * nd * (nd - 1.0));
        // a run of three is needed to contribute here, so at n = 2 the numerator
        // is zero and only the divisor would object
        double tiedTriples = (n > 2) ? first[2] * second[2] / (9.0 * nd * (nd - 1.0) * (nd - 2.0)) : 0.0;
        double variance = (v0 - first[1] - second[1]) / 18.0 + tiedPairs + tiedTriples;
        boolean hasTies = first[0] > 0.0 || second[0] > 0.0;
        return new Result(tau, s, variance, hasTies);
    }

    /**
     * The exact upper tail {@code P[S >= s]} under the null hypothesis, for
     * samples with no ties.
     * <p>
     * Without ties every pair is concordant or discordant, so {@code S} is
     * {@code n (n - 1) / 2 - 2 K} where {@code K} counts the pairs the second
     * ordering has out of step -- the inversions of one ranking against the
     * other. All {@code n!} pairings being equally likely, the null of
     * {@code K} is the number of permutations with each inversion count, built
     * one observation at a time: a new observation joins ahead of any number of
     * those already placed, from none to all of them, each equally likely, so
     * each row of the table is a running average of the last over a window that
     * grows by one.
     * <p>
     * That window is carried rather than re-summed, which costs a subtraction
     * and destroys the far end of the table -- at {@code n = 100} the
     * probabilities above about five sixths of the way up come back as noise
     * around {@code -4e-17} instead of the {@code 1e-158} they should be. It is
     * safe here because this method only ever sums the table from its
     * <b>near</b> end, where nothing has been subtracted yet and
     * {@code 1 / n!} is reproduced to every digit, and because the null is
     * symmetric: the smaller of the two tails is always the one this direction
     * gives. A {@code s} far enough below zero to need the other end is
     * answered as a number within {@code 5e-14} of one, which is the right
     * answer to four figures more than anyone reads.
     * <p>
     * {@code S} is a whole number of the same parity as {@code n (n - 1) / 2},
     * and a value between two attainable ones is answered for the next one up.
     * <p>
     * The lower tail needs no method of its own: the null is symmetric about
     * zero, so {@code P[S <= s]} is this function at {@code -s}.
     *
     * @param n
     *            the number of pairs, at least two and at most
     *            {@link #EXACT_LIMIT}
     * @param s
     *            the concordance excess the statistic must reach
     * @return the probability that {@code S} is at least {@code s}, in
     *         {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code n} is smaller than two or above
     *             {@link #EXACT_LIMIT}, or if {@code s} is {@code NaN}
     */
    public static double barFExact(int n, double s) {
        if (n < 2) {
            throw new IllegalArgumentException("n must be at least two : " + n);
        }
        if (n > EXACT_LIMIT) {
            throw new IllegalArgumentException("n is " + n + ", above the exact limit of " + EXACT_LIMIT);
        }
        if (Double.isNaN(s)) {
            throw new IllegalArgumentException("s must not be NaN");
        }
        double total = (double) n * (n - 1) / 2.0;
        if (s > total) {
            return 0.0;
        }
        double[] pmf = pmf(n);
        // S >= s is K <= (total - s) / 2
        double reach = Math.floor((total - s) / 2.0);
        if (reach < 0.0) {
            return 0.0;
        }
        int upTo = (reach >= pmf.length) ? pmf.length - 1 : (int) reach;
        double tail = 0.0;
        for (int k = 0; k <= upTo; k++) {
            tail += pmf[k];
        }
        return Math.min(1.0, Math.max(0.0, tail));
    }

    /**
     * The normal approximation to the upper tail {@code P[S >= s]}, which is
     * what serves above {@link #EXACT_LIMIT} and whenever the ranking has ties.
     * <p>
     * {@code S} has mean zero and the variance the result carries, and the tail
     * is read off a normal with half a step of continuity correction. The
     * correction is measured, not assumed: against the exact tail over the
     * p-values between 0.001 and 0.2 it takes the error from 4.6e-2 to 6.3e-3
     * at {@code n = 5} and from 3.3e-3 to 9.7e-4 at {@code n = 50}.
     * <p>
     * A step is two without ties, because every pair is then concordant or
     * discordant and moving one of them across changes {@code S} by two; with
     * ties it is one. {@link Result#hasTies} is what says which, and is why
     * this method takes the whole result rather than a variance.
     *
     * @param result
     *            the coefficient this is the null distribution of
     * @param s
     *            the concordance excess the statistic must reach
     * @return the approximate probability that {@code S} is at least
     *         {@code s}, in {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code result} is {@code null}, if {@code s} is
     *             {@code NaN}, or if the result carries no positive variance
     */
    public static double barFAsymptotic(Result result, double s) {
        if (result == null) {
            throw new IllegalArgumentException("result must not be null");
        }
        if (Double.isNaN(s)) {
            throw new IllegalArgumentException("s must not be NaN");
        }
        if (!(result.nullVariance > 0.0)) {
            throw new IllegalArgumentException("the null variance is not positive : " + result.nullVariance);
        }
        double halfAStep = result.hasTies ? 0.5 : 1.0;
        double z = (s - halfAStep) / Math.sqrt(result.nullVariance);
        return Math.min(1.0, Math.max(0.0, ProbabilityFuncs.normal(-z)));
    }

    /**
     * {@code P[K = k]} for {@code k = 0 .. n (n - 1) / 2}, the number of
     * permutations with each inversion count over {@code n!}.
     */
    private static double[] pmf(int n) {
        int maxK = n * (n - 1) / 2;
        double[] previous = new double[maxK + 1];
        double[] current = new double[maxK + 1];
        previous[0] = 1.0;
        for (int i = 2; i <= n; i++) {
            double window = 0.0;
            for (int k = 0; k <= maxK; k++) {
                window += previous[k];
                if (k >= i) {
                    window -= previous[k - i];
                }
                current[k] = window / i;
            }
            double[] swap = previous;
            previous = current;
            current = swap;
        }
        return previous;
    }

    /**
     * Over the runs of equal values in a sorted sample -- or of equal pairs,
     * when {@code second} is given -- the three sums the tie corrections need:
     * {@code t (t - 1) / 2}, {@code t (t - 1) (2 t + 5)} and
     * {@code t (t - 1) (t - 2)}.
     */
    private static double[] tieAggregates(double[] sorted, double[] second) {
        double pairs = 0.0;
        double weighted = 0.0;
        double triples = 0.0;
        int n = sorted.length;
        int i = 0;
        while (i < n) {
            int j = i;
            while (j + 1 < n && sorted[j + 1] == sorted[i] && (second == null || second[j + 1] == second[i])) {
                j++;
            }
            double t = j - i + 1;
            if (t > 1.0) {
                pairs += t * (t - 1.0) / 2.0;
                weighted += t * (t - 1.0) * (2.0 * t + 5.0);
                triples += t * (t - 1.0) * (t - 2.0);
            }
            i = j + 1;
        }
        return new double[] { pairs, weighted, triples };
    }

    /** Sorts the pairs by the first value, ties broken by the second. */
    private static void sortPairs(double[] a, double[] b, int from, int to, double[] ta, double[] tb) {
        if (to - from < 2) {
            return;
        }
        int mid = (from + to) >>> 1;
        sortPairs(a, b, from, mid, ta, tb);
        sortPairs(a, b, mid, to, ta, tb);
        int i = from;
        int j = mid;
        int k = from;
        while (i < mid && j < to) {
            if (a[i] < a[j] || (a[i] == a[j] && b[i] <= b[j])) {
                ta[k] = a[i];
                tb[k] = b[i];
                i++;
            } else {
                ta[k] = a[j];
                tb[k] = b[j];
                j++;
            }
            k++;
        }
        while (i < mid) {
            ta[k] = a[i];
            tb[k] = b[i];
            i++;
            k++;
        }
        while (j < to) {
            ta[k] = a[j];
            tb[k] = b[j];
            j++;
            k++;
        }
        System.arraycopy(ta, from, a, from, to - from);
        System.arraycopy(tb, from, b, from, to - from);
    }

    /**
     * Sorts and returns how many exchanges it took, which is the number of
     * pairs that were out of order -- the discordant ones.
     */
    private static long sortAndCount(double[] a, int from, int to, double[] t) {
        if (to - from < 2) {
            return 0L;
        }
        int mid = (from + to) >>> 1;
        long exchanges = sortAndCount(a, from, mid, t) + sortAndCount(a, mid, to, t);
        int i = from;
        int j = mid;
        int k = from;
        while (i < mid && j < to) {
            if (a[i] <= a[j]) {
                t[k++] = a[i++];
            } else {
                t[k++] = a[j++];
                // everything still waiting in the left half is larger than the
                // value just taken from the right one
                exchanges += mid - i;
            }
        }
        while (i < mid) {
            t[k++] = a[i++];
        }
        while (j < to) {
            t[k++] = a[j++];
        }
        System.arraycopy(t, from, a, from, to - from);
        return exchanges;
    }

    private KendallTau() {
        throw new AssertionError();
    }
}
