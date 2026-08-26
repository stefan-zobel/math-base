package math.stats.gof;

/**
 * The null distribution of the Cramer-von Mises statistic
 * {@code W_n^2 = 1/(12n) + sum (u_(i) - (2i - 1)/(2n))^2} for an ordered sample
 * of {@code n} independent uniforms {@code U_i} over {@code (0,1)}.
 * <p>
 * The limit of {@code W_n^2} is a weighted sum of chi-squares whose weights are
 * known in closed form, {@code lambda_j = 1 / (pi^2 j^2)}, so the distribution
 * is {@link WeightedChiSquare} and no table is needed. The series is truncated
 * and what it drops is <b>subtracted from the argument</b> rather than ignored:
 * the discarded weights sum to a known number and contribute that much to the
 * statistic on average. With the correction, thirty terms land within
 * {@code 5.2e-7} of twenty thousand; without it the same thirty terms would be
 * short by {@code 3.3e-3}, four orders of magnitude worse, for a correction
 * that costs one subtraction.
 * <p>
 * Finite {@code n} is covered by Stephens' modification, which reads
 * {@code (W_n^2 - 0.4/n + 0.6/n^2)(1 + 1/n)} against the limit.
 * <b>It is a tail device, and it is worth knowing which half of the
 * distribution it was built for.</b> Measured against a drawn null of 100000
 * samples at each size:
 * <table border="1">
 * <caption>Error of the computed tail, by sample size and level</caption>
 * <tr><th>{@code n}</th><th>{@code 0.75}</th><th>{@code 0.50}</th>
 *     <th>{@code 0.10}</th><th>{@code 0.05}</th><th>{@code 0.01}</th></tr>
 * <tr><td>3</td><td>+0.250</td><td>+0.200</td><td>-0.0089</td><td>-0.0094</td>
 *     <td>-0.0010</td></tr>
 * <tr><td>10</td><td>+0.172</td><td>+0.106</td><td>+0.0017</td><td>-0.0011</td>
 *     <td>-0.0004</td></tr>
 * <tr><td>30</td><td>+0.062</td><td>+0.034</td><td>+0.0013</td><td>-0.0001</td>
 *     <td>-0.0002</td></tr>
 * <tr><td>200</td><td>+0.010</td><td>+0.005</td><td>-0.0003</td><td>+0.0003</td>
 *     <td>+0.0000</td></tr>
 * </table>
 * <p>
 * From the ten percent point downwards it is right to within {@code 0.01} at
 * every size tried and to within {@code 0.0025} from {@code n = 5} up, which is
 * the whole range in which a p-value decides anything. Above it the modification
 * carries a small sample too far and the error grows to a quarter at
 * {@code n = 3}, shrinking like {@code 1/n}. A p-value of {@code 0.6} that is
 * really {@code 0.4} costs nothing; one of {@code 0.04} that is really
 * {@code 0.06} would cost everything, and that is the one this is accurate for.
 * <p>
 * Leaving the modification out is not the safer choice. The unmodified limit is
 * <em>liberal</em> in exactly the place that matters and rejects {@code 0.0137}
 * of the time at the one percent level for {@code n = 5}.
 * <p>
 * A single observation is the one case the modification cannot serve, and it
 * does not have to: {@code W_1^2} is {@code 1/12 + (u - 1/2)^2}, so it lives on
 * {@code [1/12, 1/3]} and its distribution function is {@code 2 sqrt(x - 1/12)}
 * exactly. Read through the limit instead, a statistic of {@code 0.5} at
 * {@code n = 1} would come back with a p-value of {@code 2.9e-4} for a value
 * the statistic cannot reach at all.
 * <p>
 * https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93von_Mises_criterion
 *
 * @since 1.5.3
 */
public final class CramerVonMises {

    /**
     * How many terms of {@code sum 1 / (pi^2 j^2)} are carried.
     * <p>
     * The rest is folded into the argument, which is what makes the number
     * small. Measured against a twenty thousand term reference at the four
     * published points, the worst absolute error is {@code 2.4e-5} at eight
     * terms, {@code 1.7e-6} at twenty and {@code 5.2e-7} here. The cost rises in
     * <em>both</em> directions -- with too few weights Imhof's remainder bound
     * pushes the upper limit of the integral much further out -- and is flattest
     * between twelve and fifty.
     * <p>
     * A p-value costs about ten milliseconds, which is four orders of magnitude
     * more than a table lookup would. That is the price of not having a table:
     * the integrand is evaluated a few thousand times and each evaluation is an
     * arc tangent and a logarithm per weight. It is invisible for a test and it
     * is not invisible in a loop, so a caller measuring a null distribution by
     * simulation should compare statistics against one critical value rather
     * than compute a p-value per replication.
     */
    private static final int SERIES_TERMS = 30;

    /** The weights carried, {@code 1 / (pi^2 j^2)} for {@code j = 1 .. SERIES_TERMS}. */
    private static final double[] LAMBDA = weights();

    /**
     * What the discarded weights contribute on average, {@code 1/6} being the
     * sum of the whole series.
     */
    private static final double DROPPED = dropped();

    /** The smallest value {@code W_1^2} can take, and the largest. */
    private static final double ONE_TWELFTH = 1.0 / 12.0;
    private static final double ONE_THIRD = 1.0 / 3.0;

    /**
     * The accuracy asked of the inversion.
     * <p>
     * Four orders of magnitude looser than the default, which costs
     * {@code 4.9e-11} and saves three quarters of the running time. There is
     * nothing to spend that accuracy on: the finite-sample correction below is
     * itself good to about {@code 1e-3}, so a tail computed to {@code 1e-13}
     * would be quoting digits the correction has already lost.
     */
    private static final double TOLERANCE = 1.0e-9;

    /**
     * Computes the <b>complementary</b> distribution function
     * {@code P[W_n^2 >= x]}.
     *
     * @param n
     *            the sample size, {@code 1} or greater
     * @param x
     *            the value of the statistic
     * @return {@code P[W_n^2 >= x]}
     * @throws IllegalArgumentException
     *             if {@code n} is not strictly positive or {@code x} is not a
     *             number
     */
    public static double barF(int n, double x) {
        requireArguments(n, x);
        if (n == 1) {
            // written so that the tail does not cancel: with
            // w = x - 1/12, one minus 2 sqrt(w) is (1 - 4w) / (1 + 2 sqrt(w))
            double w = x - ONE_TWELFTH;
            if (w <= 0.0) {
                return 1.0;
            }
            double root = Math.sqrt(w);
            return root >= 0.5 ? 0.0 : 4.0 * (ONE_THIRD - x) / (1.0 + 2.0 * root);
        }
        if (x <= 0.0) {
            return 1.0;
        }
        return WeightedChiSquare.barF(LAMBDA, modified(n, x) - DROPPED, TOLERANCE);
    }

    /**
     * Computes the distribution function {@code P[W_n^2 <= x]}.
     * <p>
     * It is read from {@link WeightedChiSquare#cdf} directly rather than as one
     * minus {@link #barF}, so that a small lower tail keeps its digits.
     *
     * @param n
     *            the sample size, {@code 1} or greater
     * @param x
     *            the value of the statistic
     * @return {@code P[W_n^2 <= x]}
     * @throws IllegalArgumentException
     *             if {@code n} is not strictly positive or {@code x} is not a
     *             number
     */
    public static double cdf(int n, double x) {
        requireArguments(n, x);
        if (n == 1) {
            double w = x - ONE_TWELFTH;
            if (w <= 0.0) {
                return 0.0;
            }
            double p = 2.0 * Math.sqrt(w);
            return p >= 1.0 ? 1.0 : p;
        }
        if (x <= 0.0) {
            return 0.0;
        }
        return WeightedChiSquare.cdf(LAMBDA, modified(n, x) - DROPPED, TOLERANCE);
    }

    /**
     * Stephens' modification, which carries the statistic of a sample of
     * {@code n} onto the scale the limit distribution is stated on.
     * <p>
     * It can come out negative for a statistic near zero, which is not a
     * problem: every weight is positive, so the sum is positive with
     * probability one and both tails answer such an argument without being
     * asked twice.
     */
    private static double modified(int n, double x) {
        return (x - 0.4 / n + 0.6 / (n * (double) n)) * (1.0 + 1.0 / n);
    }

    private static double[] weights() {
        double[] lambda = new double[SERIES_TERMS];
        for (int j = 1; j <= SERIES_TERMS; j++) {
            lambda[j - 1] = 1.0 / (Math.PI * Math.PI * j * (double) j);
        }
        return lambda;
    }

    private static double dropped() {
        double kept = 0.0;
        for (int j = 0; j < LAMBDA.length; j++) {
            kept += LAMBDA[j];
        }
        // sum_{j>=1} 1 / (pi^2 j^2) = (1/pi^2) (pi^2 / 6)
        return 1.0 / 6.0 - kept;
    }

    private static void requireArguments(int n, double x) {
        if (n <= 0) {
            throw new IllegalArgumentException("n <= 0 : " + n);
        }
        if (Double.isNaN(x)) {
            throw new IllegalArgumentException("x must not be NaN");
        }
    }

    private CramerVonMises() {
        throw new AssertionError();
    }
}
