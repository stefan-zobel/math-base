package math.stats;

import java.util.Arrays;

/**
 * Corrections for having run more than one test.
 * <p>
 * A single test at {@code alpha = 0.05} is wrong one time in twenty by
 * construction. Twenty of them are wrong once on average, which is not a defect
 * of any one test and cannot be fixed inside one: it has to be answered by
 * looking at the whole family of p-values at once, which is what these methods
 * do.
 * <p>
 * There are two promises to choose between, and the choice is about how many
 * tests there are.
 * <p>
 * {@link #benjaminiHochberg} and {@link #benjaminiYekutieli} control the
 * <b>false discovery rate</b>, the expected share of the rejections that are
 * wrong. That is the useful promise when there are many tests and some of them
 * are expected to be real -- which is the situation that produces many tests in
 * the first place. Controlling the false discovery rate at ten percent means
 * accepting that about one in ten of what is found is noise.
 * <p>
 * {@link #holmBonferroni} controls the <b>family-wise error rate</b>, the
 * probability of <em>any</em> wrong rejection at all. That is the promise worth
 * having for the handful of comparisons that follow a rejected omnibus test --
 * three groups make three pairs -- where one wrong answer out of three is not a
 * rate one is willing to run at. Over a thousand tests the same promise means
 * finding almost nothing, which is why it is not the only one here.
 * <p>
 * Each method returns <b>adjusted p-values in the order they were given in</b>.
 * Rejecting every hypothesis whose adjusted value is at most {@code alpha} is
 * exactly the original sequential decision at {@code alpha}, so nothing is lost
 * by handing back numbers rather than a verdict, and the level no longer has to
 * be chosen before the correction is applied.
 * <p>
 * https://en.wikipedia.org/wiki/False_discovery_rate
 * <p>
 * https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method
 *
 * @since 1.5.3
 */
public final class MultipleTesting {

    /**
     * Benjamini-Hochberg adjusted p-values, which control the false discovery
     * rate when the p-values are independent.
     * <p>
     * The guarantee also survives positive regression dependence, which covers
     * the common case of tests sharing a control group or a variance estimate,
     * but not dependence in general. When that cannot be argued,
     * {@link #benjaminiYekutieli} pays a known price for a guarantee that holds
     * whatever the dependence is.
     * <p>
     * The adjusted values are monotone in the p-values and equal p-values get
     * equal adjusted values. None is ever smaller than the p-value it came
     * from, and none is ever larger than the Bonferroni value
     * {@code min(1, m p)}.
     *
     * @param pValues
     *            the p-values of the family, each in {@code [0, 1]}, at least
     *            one of them. Not modified
     * @return the adjusted p-values, in the order {@code pValues} were given in
     * @throws IllegalArgumentException
     *             if {@code pValues} is {@code null}, is empty, or holds a
     *             value that is not a probability
     */
    public static double[] benjaminiHochberg(double[] pValues) {
        int m = requireProbabilities(pValues, "pValues");
        return stepUp(pValues, m, 1.0);
    }

    /**
     * Benjamini-Yekutieli adjusted p-values, which control the false discovery
     * rate whatever the dependence between the p-values is.
     * <p>
     * This is {@link #benjaminiHochberg} multiplied by
     * {@code 1 + 1/2 + ... + 1/m} and clamped at one, and that factor is the
     * whole price: about {@code 5.2} for ten tests, {@code 7.5} for a hundred
     * and {@code 14.4} for a million, since it grows only as the logarithm.
     * Paying it buys a guarantee that needs no argument about how the tests are
     * related, which is worth having when they clearly are related and it is
     * not clear how.
     *
     * @param pValues
     *            the p-values of the family, each in {@code [0, 1]}, at least
     *            one of them. Not modified
     * @return the adjusted p-values, in the order {@code pValues} were given in
     * @throws IllegalArgumentException
     *             if {@code pValues} is {@code null}, is empty, or holds a
     *             value that is not a probability
     */
    public static double[] benjaminiYekutieli(double[] pValues) {
        int m = requireProbabilities(pValues, "pValues");
        return stepUp(pValues, m, harmonic(m));
    }

    /**
     * Holm-Bonferroni adjusted p-values, which control the family-wise error
     * rate whatever the dependence between the p-values is.
     * <p>
     * This is the step-<em>down</em> counterpart of the two above: the smallest
     * p-value is judged against {@code alpha / m}, the next against
     * {@code alpha / (m - 1)}, and so on until one of them survives, after which
     * nothing is rejected. In adjusted form that is
     * {@code q_(k) = max over j <= k of (m - j + 1) p_(j)}, clamped at one.
     * <p>
     * <b>It is uniformly better than Bonferroni and there is no reason to prefer
     * that one.</b> Every adjusted value here is at most
     * {@code min(1, m p)} -- equal for the smallest p-value and smaller for all
     * the others -- and the guarantee is the same, with no assumption about how
     * the tests are related. The price against the false discovery rate
     * procedures is the one the class comment describes: a stricter promise,
     * fewer rejections.
     * <p>
     * The adjusted values are monotone in the p-values, equal p-values get equal
     * adjusted values, and none is ever smaller than the p-value it came from.
     *
     * @param pValues
     *            the p-values of the family, each in {@code [0, 1]}, at least
     *            one of them. Not modified
     * @return the adjusted p-values, in the order {@code pValues} were given in
     * @throws IllegalArgumentException
     *             if {@code pValues} is {@code null}, is empty, or holds a
     *             value that is not a probability
     */
    public static double[] holmBonferroni(double[] pValues) {
        int m = requireProbabilities(pValues, "pValues");
        double[] ascending = pValues.clone();
        Arrays.sort(ascending);

        double[] byRank = new double[m];
        double running = Double.NEGATIVE_INFINITY;
        for (int k = 1; k <= m; k++) {
            // walking the ranks upwards turns the inner maximum into a running
            // one, and that running maximum is also what makes the result
            // monotone and what makes tied p-values come out equal without a
            // special case: where p_(k) equals p_(k+1) the later term is the
            // smaller of the two, so the maximum does not move
            double candidate = ascending[k - 1] * (m - k + 1);
            if (candidate > running) {
                running = candidate;
            }
            byRank[k - 1] = running;
        }

        double[] adjusted = new double[m];
        for (int i = 0; i < m; i++) {
            // as in stepUp: the p-values were copied, so the search always finds
            // one of the positions holding this value, and which one does not
            // matter. The multiplier is a whole number and never below one, so
            // no candidate can round below the p-value it came from and there is
            // nothing to repair here
            double q = byRank[Arrays.binarySearch(ascending, pValues[i])];
            adjusted[i] = q > 1.0 ? 1.0 : q;
        }
        return adjusted;
    }

    /**
     * The step-up sweep the two false discovery rate procedures are:
     * {@code q_(k) = min(1, min over j >= k of p_(j) factor m / j)}.
     * <p>
     * Walking the ranks downwards turns the inner minimum into a running one,
     * so the whole thing is a single backwards pass. That running minimum is
     * also what makes the result monotone, and what makes tied p-values come
     * out equal without a special case: where {@code p_(k)} equals
     * {@code p_(k+1)}, the term at {@code j = k+1} is the smaller of the two,
     * so the minimum from {@code k} down is the one from {@code k+1} down.
     * <p>
     * The clamp at one is applied last rather than inside the sweep, which is
     * what makes {@link #benjaminiYekutieli} exactly
     * {@link #benjaminiHochberg} scaled and then clamped.
     */
    private static double[] stepUp(double[] pValues, int m, double factor) {
        double[] ascending = pValues.clone();
        Arrays.sort(ascending);

        double[] byRank = new double[m];
        double running = Double.POSITIVE_INFINITY;
        for (int k = m; k >= 1; k--) {
            double candidate = ascending[k - 1] * factor * m / k;
            if (candidate < running) {
                running = candidate;
            }
            // every term of the minimum is at least p_(k), since j never
            // exceeds m and p_(j) never falls below p_(k) -- but the largest
            // rank computes p * factor * m / m, and two roundings of that can
            // land one unit in the last place below p. Repairing it here is
            // what makes "never smaller than the p-value it came from" true
            // rather than nearly true
            byRank[k - 1] = running < ascending[k - 1] ? ascending[k - 1] : running;
        }

        double[] adjusted = new double[m];
        for (int i = 0; i < m; i++) {
            // the p-values were copied, so the search always finds one of the
            // positions holding this value. Which one does not matter: equal
            // p-values earn equal adjusted values, as above
            int rank = Arrays.binarySearch(ascending, pValues[i]);
            double q = byRank[rank];
            adjusted[i] = q > 1.0 ? 1.0 : q;
        }
        return adjusted;
    }

    /**
     * {@code 1 + 1/2 + ... + 1/m}, summed from the smallest term upwards so
     * that the small ones are not lost against a partial sum that has already
     * grown past them.
     */
    private static double harmonic(int m) {
        double sum = 0.0;
        for (int j = m; j >= 1; j--) {
            sum += 1.0 / j;
        }
        return sum;
    }

    private static int requireProbabilities(double[] pValues, String name) {
        if (pValues == null) {
            throw new IllegalArgumentException(name + " must not be null");
        }
        if (pValues.length == 0) {
            throw new IllegalArgumentException(name + " must not be empty");
        }
        for (int i = 0; i < pValues.length; i++) {
            // written this way round so that NaN is refused with the rest
            if (!(pValues[i] >= 0.0 && pValues[i] <= 1.0)) {
                throw new IllegalArgumentException(
                        name + "[" + i + "] is not a probability : " + pValues[i]);
            }
        }
        return pValues.length;
    }

    private MultipleTesting() {
        throw new AssertionError();
    }
}
