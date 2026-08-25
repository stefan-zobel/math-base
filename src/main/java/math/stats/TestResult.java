package math.stats;

import java.util.Locale;

/**
 * The outcome of one hypothesis test: a statistic, the p-value it earns under
 * the null hypothesis, and which tail that p-value covers.
 * <p>
 * The type carries only what every test has. Where a test yields more -- a
 * confidence interval, a second degrees of freedom -- it returns a type of its
 * own holding one of these, rather than adding a field that would be
 * {@link Double#NaN} for everything else.
 * <p>
 * The significance level is deliberately not a field. It belongs to the reader,
 * not to the computation, so it is an argument of
 * {@link #rejectsAt(double)} and can be asked as often as one likes without
 * running the test again.
 *
 * @since 1.5.3
 */
public final class TestResult {

    /** What was tested, for example {@code "one-sample t"}. */
    public final String test;
    /** The value of the test statistic. */
    public final double statistic;
    /**
     * The probability, under the null hypothesis, of a statistic at least as
     * extreme as {@link #statistic} in the direction {@link #alternative}
     * names. It is a number in {@code [0, 1]} and never {@link Double#NaN}: a
     * test that cannot be computed throws instead.
     */
    public final double pValue;
    /** Which departure from the null hypothesis {@link #pValue} covers. */
    public final Alternative alternative;
    /**
     * Degrees of freedom of the null distribution, or {@link Double#NaN} for a
     * test whose null distribution has none. It is a {@code double} because
     * Welch's approximation is not an integer.
     */
    public final double degreesOfFreedom;

    TestResult(String test, double statistic, double pValue, Alternative alternative, double degreesOfFreedom) {
        this.test = test;
        this.statistic = statistic;
        this.pValue = pValue;
        this.alternative = alternative;
        this.degreesOfFreedom = degreesOfFreedom;
    }

    /**
     * Whether the null hypothesis is rejected at the significance level
     * {@code alpha}, which is {@link #pValue} {@code < alpha}.
     * <p>
     * Asking after the fact is cheap; choosing {@code alpha} after seeing
     * {@link #pValue} is not a test.
     *
     * @param alpha
     *            the significance level, strictly between {@code 0} and
     *            {@code 1}
     * @return {@code true} if the null hypothesis is rejected at that level
     * @throws IllegalArgumentException
     *             if {@code alpha} does not lie in {@code (0, 1)}
     */
    public boolean rejectsAt(double alpha) {
        if (!(alpha > 0.0 && alpha < 1.0)) {
            throw new IllegalArgumentException("alpha must lie in (0, 1) : " + alpha);
        }
        return pValue < alpha;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder b = new StringBuilder(96);
        b.append(test).append(": statistic = ").append(String.format(Locale.ROOT, "%.6g", statistic));
        if (!Double.isNaN(degreesOfFreedom)) {
            b.append(", df = ");
            if (degreesOfFreedom == Math.rint(degreesOfFreedom) && Math.abs(degreesOfFreedom) < 1.0e15) {
                // Long.toString carries no locale, and "11" reads better than
                // "11.0000" for the degrees of freedom that are whole numbers
                b.append((long) degreesOfFreedom);
            } else {
                b.append(String.format(Locale.ROOT, "%.6g", degreesOfFreedom));
            }
        }
        b.append(", p = ").append(String.format(Locale.ROOT, "%.6g", pValue));
        b.append(" (").append(describe(alternative)).append(')');
        return b.toString();
    }

    /** The alternative in words, so the printed line says what the p-value means. */
    private static String describe(Alternative alternative) {
        switch (alternative) {
        case LESS:
            return "one-sided, less";
        case GREATER:
            return "one-sided, greater";
        default:
            return "two-sided";
        }
    }
}
