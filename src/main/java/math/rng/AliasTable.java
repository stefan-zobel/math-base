package math.rng;

/**
 * A prepared categorical distribution over {@code 0 .. n-1}, drawn from in
 * constant time by Walker's alias method in Vose's formulation.
 * <p>
 * Inversion over a cumulative table costs a binary search, so it grows with the
 * number of outcomes. This costs one {@code int} and one {@code double} per
 * draw whatever {@code n} is, in exchange for an {@code O(n)} setup and two
 * arrays: the table is worth building when it will be drawn from more than a
 * handful of times, which is the case a sampler is usually in.
 * <p>
 * The construction pairs each outcome that is lighter than average with one
 * that is heavier, until every outcome carries at most a full share of its own
 * and the remainder of one other. A draw then picks a column uniformly and
 * flips a biased coin to decide which of the two it lands on.
 * <p>
 * Instances are immutable and can be shared between threads. The generator
 * passed to {@link #sample(PseudoRandom)} cannot -- that restriction belongs to
 * the generator, not to this table.
 * <p>
 * https://en.wikipedia.org/wiki/Alias_method
 *
 * @since 1.5.3
 */
public final class AliasTable {

    /**
     * The probability of staying in column {@code i} rather than following its
     * alias, each in {@code [0, 1]}.
     */
    private final double[] stay;

    /** Where column {@code i} sends a draw that does not stay. */
    private final int[] alias;

    private AliasTable(double[] stay, int[] alias) {
        this.stay = stay;
        this.alias = alias;
    }

    /**
     * Builds the table for the given weights.
     * <p>
     * The weights need not sum to one; they are normalized, so counts serve as
     * well as probabilities. A weight of zero is allowed and its outcome is
     * then never drawn.
     *
     * @param weights
     *            the weight of each outcome, at least one of them, each finite
     *            and not negative, and not all zero. Not modified
     * @return a table over {@code 0 .. weights.length - 1}
     * @throws IllegalArgumentException
     *             if {@code weights} is {@code null}, is empty, holds a value
     *             that is negative or not finite, or sums to zero
     */
    public static AliasTable of(double[] weights) {
        if (weights == null) {
            throw new IllegalArgumentException("weights must not be null");
        }
        int n = weights.length;
        if (n == 0) {
            throw new IllegalArgumentException("weights must not be empty");
        }
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            double w = weights[i];
            if (!(w >= 0.0) || w == Double.POSITIVE_INFINITY) {
                throw new IllegalArgumentException(
                        "weights[" + i + "] is not a finite non-negative number : " + w);
            }
            sum += w;
        }
        if (!(sum > 0.0) || sum == Double.POSITIVE_INFINITY) {
            throw new IllegalArgumentException("the weights must sum to a finite positive number, not " + sum);
        }

        // each outcome's share of n, so that the average is exactly one and
        // "light" and "heavy" mean below and above it
        double[] stay = new double[n];
        double scale = n / sum;
        for (int i = 0; i < n; i++) {
            stay[i] = weights[i] * scale;
        }

        // the two worklists share one array: the light ones grow from the
        // bottom and the heavy ones from the top, and they cannot meet because
        // together they hold at most n entries
        int[] work = new int[n];
        int light = 0;
        int heavy = n;
        for (int i = n - 1; i >= 0; i--) {
            if (stay[i] < 1.0) {
                work[light++] = i;
            } else {
                work[--heavy] = i;
            }
        }

        int[] alias = new int[n];
        while (light > 0 && heavy < n) {
            int small = work[--light];
            int large = work[heavy++];
            alias[small] = large;
            // written as (large + small) - 1 rather than large - (1 - small),
            // which is Vose's own note: the first form does not lose the
            // digits of a large residue against a tiny one
            stay[large] = (stay[large] + stay[small]) - 1.0;
            if (stay[large] < 1.0) {
                work[light++] = large;
            } else {
                work[--heavy] = large;
            }
        }
        // whatever is left over carries a full share, up to the rounding that
        // put it on the wrong list in the first place
        while (heavy < n) {
            stay[work[heavy++]] = 1.0;
        }
        while (light > 0) {
            stay[work[--light]] = 1.0;
        }
        return new AliasTable(stay, alias);
    }

    /**
     * Draws one outcome.
     *
     * @param prng
     *            the generator to draw from
     * @return an outcome in {@code 0 .. outcomes() - 1}
     * @throws IllegalArgumentException
     *             if {@code prng} is {@code null}
     */
    public int sample(PseudoRandom prng) {
        if (prng == null) {
            throw new IllegalArgumentException("prng must not be null");
        }
        return draw(prng);
    }

    /** The same draw without the null check, for a caller that has made it already. */
    int draw(PseudoRandom prng) {
        int column = prng.nextInt(stay.length);
        return prng.nextDouble() < stay[column] ? column : alias[column];
    }

    /**
     * The number of outcomes this table covers.
     *
     * @return the number of outcomes, one or more
     */
    public int outcomes() {
        return stay.length;
    }

    /**
     * The probability that column {@code i} keeps a draw that lands on it,
     * exposed for the tests that check the construction rather than its
     * consequences.
     *
     * @param i
     *            the column
     * @return the probability of staying in column {@code i}
     */
    double stayProbability(int i) {
        return stay[i];
    }

    /**
     * Where column {@code i} sends a draw that does not stay.
     *
     * @param i
     *            the column
     * @return the alias of column {@code i}
     */
    int aliasOf(int i) {
        return alias[i];
    }
}
