package math.stats.rank;

import math.distribution.StudentT;

/**
 * Spearman's rank correlation coefficient {@code rho} -- Pearson's correlation
 * of the midranks -- and its null distribution when the two rankings are
 * paired at random.
 * <p>
 * https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
 *
 * @since 1.5.3
 */
public final class SpearmanRho {

    /**
     * The largest sample the exact method will accept.
     * <p>
     * The exact null is an enumeration of the {@code n!} pairings, and there
     * is no cheaper way to it: measured at 0.66 ms for {@code n = 8}, 7.2 ms
     * for {@code n = 9} and 67 ms for {@code n = 10}, each step about ten
     * times the last. At {@code n = 10} the approximation is already within
     * 3.4e-3 of the exact tail over the p-values a reader acts on, which is
     * what makes nine the place to stop rather than ten.
     */
    public static final int EXACT_LIMIT = 9;

    /**
     * Spearman's {@code rho} for two paired samples.
     * <p>
     * The definition is Pearson's correlation of the two midrankings. The
     * familiar {@code 1 - 6 sum d^2 / (n^3 - n)} is the same number only when
     * neither sample has a tie; with ties it is not, and this method computes
     * the definition.
     *
     * @param x
     *            the first sample, at least two observations, none of them
     *            {@code NaN}, not all equal
     * @param y
     *            the second sample, as many observations as {@code x}, under
     *            the same conditions
     * @return the coefficient, in {@code [-1, 1]}
     * @throws IllegalArgumentException
     *             if either sample is {@code null}, holds fewer than two
     *             observations or a {@code NaN}, if the two are of different
     *             lengths, or if either takes only one value
     */
    public static double coefficient(double[] x, double[] y) {
        if (x == null || y == null) {
            throw new IllegalArgumentException("neither sample may be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException(
                    "x and y must be paired, but their lengths are " + x.length + " and " + y.length);
        }
        return coefficient(Ranks.of(x), Ranks.of(y));
    }

    /**
     * Spearman's {@code rho} from two rankings that have already been formed,
     * for a caller who needs the ranks themselves as well.
     *
     * @param x
     *            the ranking of the first sample
     * @param y
     *            the ranking of the second sample, of the same length
     * @return the coefficient, in {@code [-1, 1]}
     * @throws IllegalArgumentException
     *             if either ranking is {@code null}, they are of different
     *             lengths, there are fewer than two observations, or either
     *             ranking has no spread, which is one value repeated
     */
    public static double coefficient(Ranks.Result x, Ranks.Result y) {
        if (x == null || y == null) {
            throw new IllegalArgumentException("neither ranking may be null");
        }
        int n = x.ranks.length;
        if (n != y.ranks.length) {
            throw new IllegalArgumentException(
                    "the two rankings must be paired, but their lengths are " + n + " and " + y.ranks.length);
        }
        if (n < 2) {
            throw new IllegalArgumentException("a correlation needs at least two observations, got " + n);
        }
        // a midranking sums to n (n + 1) / 2 whatever its ties do, so its mean
        // is (n + 1) / 2 exactly and there is nothing to estimate here
        double mean = (n + 1.0) / 2.0;
        double sxy = 0.0;
        double sxx = 0.0;
        double syy = 0.0;
        for (int i = 0; i < n; i++) {
            double a = x.ranks[i] - mean;
            double b = y.ranks[i] - mean;
            sxy += a * b;
            sxx += a * a;
            syy += b * b;
        }
        if (!(sxx > 0.0)) {
            throw new IllegalArgumentException("the first sample takes one value only, so it has no ranking to "
                    + "correlate");
        }
        if (!(syy > 0.0)) {
            throw new IllegalArgumentException("the second sample takes one value only, so it has no ranking to "
                    + "correlate");
        }
        double rho = sxy / Math.sqrt(sxx * syy);
        // the arithmetic can leave the interval by an ulp where the ranks agree
        // exactly, and a coefficient of 1.0000000000000002 is not one
        return Math.max(-1.0, Math.min(1.0, rho));
    }

    /**
     * The exact upper tail {@code P[rho >= r]} under the null hypothesis that
     * the two rankings are paired at random.
     * <p>
     * All {@code n!} pairings are equally likely, and the tail is the fraction
     * of them that reaches {@code r}. There is no recursion to be had here as
     * there is for the rank sums: {@code rho} is a sum of squared rank
     * differences over a permutation, which has no structure to build on, so
     * the pairings are enumerated. That is what makes {@link #EXACT_LIMIT}
     * nine rather than a few hundred.
     * <p>
     * The attainable values of {@code rho} are {@code 12 / (n^3 - n)} apart,
     * and {@code r} is answered for the nearest of them -- unlike the rank sum
     * statistics next door, which are whole numbers and are answered for the
     * next value up.
     * <p>
     * The lower tail needs no method of its own: the null is symmetric about
     * zero, so {@code P[rho <= r]} is this function at {@code -r}.
     *
     * @param n
     *            the number of pairs, at least two and at most
     *            {@link #EXACT_LIMIT}
     * @param rho
     *            the value the coefficient must reach
     * @return the probability that {@code rho} is at least that, in
     *         {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code n} is smaller than two or above
     *             {@link #EXACT_LIMIT}, or if {@code rho} is {@code NaN}
     */
    public static double barFExact(int n, double rho) {
        if (n < 2) {
            throw new IllegalArgumentException("n must be at least two : " + n);
        }
        if (n > EXACT_LIMIT) {
            throw new IllegalArgumentException("n is " + n + ", above the exact limit of " + EXACT_LIMIT);
        }
        if (Double.isNaN(rho)) {
            throw new IllegalArgumentException("rho must not be NaN");
        }
        double span = n * (double) n * n - n;
        double halfAnAtom = 6.0 / span;
        double[] pmf = pmf(n);
        double tail = 0.0;
        for (int k = 0; k < pmf.length; k++) {
            // the atoms come in descending order of rho, so this walks the
            // upper tail from its far end inwards
            if (1.0 - 6.0 * (2 * k) / span < rho - halfAnAtom) {
                break;
            }
            tail += pmf[k];
        }
        return Math.min(1.0, tail);
    }

    /**
     * The approximation to the upper tail {@code P[rho >= r]}, which is what
     * serves above {@link #EXACT_LIMIT} and whenever the ranking has ties.
     * <p>
     * {@code t = rho sqrt((n - 2) / (1 - rho^2))} against a {@link StudentT}
     * on {@code n - 2} degrees of freedom, the same approximation that carries
     * the correlation coefficient of a simple regression.
     * <p>
     * <b>There is no continuity correction here</b>, and that is measured
     * rather than assumed. Half an atom of {@code rho} does improve the tail
     * at the very small sizes -- 3.3e-2 to 7.6e-3 at {@code n = 5} -- but from
     * {@code n = 9} upwards it makes the tail worse, 2.7e-3 against 3.5e-3 at
     * {@code n = 11}, and above {@link #EXACT_LIMIT} is the only place this
     * method has to be right.
     * <p>
     * With ties at a small {@code n} this is the only path there is, and it is
     * the weakest thing in the class: at {@code n = 10} it is within 3.4e-3 of
     * the exact tail, at {@code n = 5} only within 3.3e-2. Measured over 20000
     * replications of a true null with the observations rounded to half units,
     * a test run at five percent rejects 5.8 percent of the time at
     * {@code n = 5} and 5.7 at {@code n = 8}, and is back to 5.0 by
     * {@code n = 12}. It errs towards rejecting, which is the direction that
     * costs the reader something.
     *
     * @param n
     *            the number of pairs, at least three
     * @param rho
     *            the value the coefficient must reach
     * @return the approximate probability that {@code rho} is at least that,
     *         in {@code [0, 1]}
     * @throws IllegalArgumentException
     *             if {@code n} is smaller than three, or if {@code rho} is
     *             {@code NaN}
     */
    public static double barFAsymptotic(int n, double rho) {
        if (n < 3) {
            throw new IllegalArgumentException("n must be at least three for n - 2 degrees of freedom : " + n);
        }
        if (Double.isNaN(rho)) {
            throw new IllegalArgumentException("rho must not be NaN");
        }
        if (!(rho * rho < 1.0)) {
            // a perfect agreement or a perfect reversal: the t is infinite and
            // the tail is one end or the other
            return rho > 0.0 ? 0.0 : 1.0;
        }
        double t = rho * Math.sqrt((n - 2.0) / (1.0 - rho * rho));
        return Math.min(1.0, Math.max(0.0, new StudentT(n - 2.0).cdf(-t)));
    }

    /**
     * {@code P[S = 2k]} for {@code k = 0 .. (n^3 - n) / 6}, where {@code S} is
     * the sum of the squared rank differences -- which is always even, since
     * the differences sum to zero. Index {@code k} is the atom at
     * {@code rho = 1 - 12k / (n^3 - n)}, the largest first.
     */
    private static double[] pmf(int n) {
        int maxS = (n * n * n - n) / 3;
        double[] counts = new double[maxS / 2 + 1];
        walk(1, n, 0, new boolean[n + 1], counts);
        double total = 0.0;
        for (int k = 0; k < counts.length; k++) {
            total += counts[k];
        }
        for (int k = 0; k < counts.length; k++) {
            counts[k] /= total;
        }
        return counts;
    }

    /** Pairs rank {@code i} with each rank still free, carrying the sum along. */
    private static void walk(int i, int n, int s, boolean[] taken, double[] counts) {
        if (i > n) {
            counts[s / 2] += 1.0;
            return;
        }
        for (int v = 1; v <= n; v++) {
            if (!taken[v]) {
                taken[v] = true;
                int d = i - v;
                walk(i + 1, n, s + d * d, taken, counts);
                taken[v] = false;
            }
        }
    }

    private SpearmanRho() {
        throw new AssertionError();
    }
}
