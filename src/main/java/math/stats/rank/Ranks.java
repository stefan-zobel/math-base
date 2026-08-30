package math.stats.rank;

import java.util.Arrays;

/**
 * The ranks of a sample, where values that occur more than once share the
 * average of the ranks they span -- the midranks {@code 1, 2.5, 2.5, 4}
 * rather than an arbitrary order among the equals.
 * <p>
 * https://en.wikipedia.org/wiki/Ranking#Ranking_in_statistics
 *
 * @since 1.5.3
 */
public final class Ranks {

    /**
     * The ranking of one sample: the midranks themselves, and what the tests
     * built on them need to know about the ties.
     */
    public static final class Result {

        /**
         * The midranks, in the order the values were given. They sum to
         * {@code n (n + 1) / 2} whether or not there are ties.
         */
        public final double[] ranks;

        /**
         * The sum of {@code t^3 - t} over the runs of equal values, where
         * {@code t} is the length of a run. It is {@code 0} without ties, and
         * it is exactly the correction the variance of every rank statistic
         * here takes.
         */
        public final double tieSum;

        /**
         * Whether any value occurs more than once. This is what decides
         * between an exact null distribution and its normal approximation:
         * the exact ones count equally likely orderings, and a tie means
         * there is no ordering to count.
         */
        public final boolean hasTies;

        Result(double[] ranks, double tieSum, boolean hasTies) {
            this.ranks = ranks;
            this.tieSum = tieSum;
            this.hasTies = hasTies;
        }
    }

    /**
     * Ranks a sample from smallest to largest, giving every run of equal
     * values the average of the ranks it spans.
     * <p>
     * Infinities are ranked like any other value -- an infinite observation
     * has a perfectly well defined position in an ordering -- but
     * {@link Double#NaN} is not, and is rejected.
     *
     * @param values
     *            the sample, at least one observation, none of them
     *            {@code NaN}
     * @return the midranks in the caller's order, and the tie information
     * @throws IllegalArgumentException
     *             if {@code values} is {@code null}, is empty, or holds a
     *             {@code NaN}
     */
    public static Result of(double[] values) {
        if (values == null) {
            throw new IllegalArgumentException("values must not be null");
        }
        if (values.length == 0) {
            throw new IllegalArgumentException("values must not be empty");
        }
        for (int i = 0; i < values.length; i++) {
            if (Double.isNaN(values[i])) {
                throw new IllegalArgumentException("values[" + i + "] is NaN, which has no place in an ordering");
            }
        }
        int n = values.length;
        double[] sorted = values.clone();
        Arrays.sort(sorted);

        double[] byPosition = new double[n];
        double tieSum = 0.0;
        boolean hasTies = false;
        int i = 0;
        while (i < n) {
            int j = i;
            while (j + 1 < n && sorted[j + 1] == sorted[i]) {
                j++;
            }
            // the ranks i+1 .. j+1 averaged, which is their midpoint
            double midrank = 0.5 * (i + j) + 1.0;
            for (int k = i; k <= j; k++) {
                byPosition[k] = midrank;
            }
            long t = j - i + 1L;
            if (t > 1L) {
                hasTies = true;
                tieSum += (double) t * t * t - t;
            }
            i = j + 1;
        }

        double[] ranks = new double[n];
        for (int k = 0; k < n; k++) {
            // every position inside a run of equal values carries the same
            // midrank, so it does not matter which of them the search returns
            // -- which is what makes the walk back to the caller's order a
            // binary search rather than a sort of boxed indices
            ranks[k] = byPosition[Arrays.binarySearch(sorted, values[k])];
        }
        return new Result(ranks, tieSum, hasTies);
    }

    private Ranks() {
        throw new AssertionError();
    }
}
