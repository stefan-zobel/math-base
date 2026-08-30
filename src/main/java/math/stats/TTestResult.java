package math.stats;

import java.util.Locale;

/**
 * The outcome of a t-test: the {@link TestResult} every test returns, plus the
 * quantity that was estimated and a confidence interval for it.
 * <p>
 * The interval and the decision are two views of the same computation. It is
 * built to match {@link TestResult#alternative}, so it has one infinite end for
 * a one-sided test, and it excludes the null value exactly when the test
 * rejects at {@code 1 - }{@link #confidenceLevel}.
 *
 * @since 1.5.3
 */
public final class TTestResult {

    /** The statistic, its p-value, the alternative and the degrees of freedom. */
    public final TestResult test;
    /**
     * The quantity the null hypothesis is about: the sample mean for the
     * one-sample test, the difference of the means for the two-sample and
     * paired tests.
     */
    public final double estimate;
    /** The standard error of {@link #estimate}. */
    public final double standardError;
    /** The level {@link #lower} and {@link #upper} were computed at. */
    public final double confidenceLevel;
    /**
     * Lower end of the confidence interval for {@link #estimate}, or
     * {@link Double#NEGATIVE_INFINITY} where the alternative is one-sided in
     * the other direction.
     */
    public final double lower;
    /**
     * Upper end of the confidence interval for {@link #estimate}, or
     * {@link Double#POSITIVE_INFINITY} where the alternative is one-sided in
     * the other direction.
     */
    public final double upper;

    TTestResult(TestResult test, double estimate, double standardError, double confidenceLevel, double lower,
            double upper) {
        this.test = test;
        this.estimate = estimate;
        this.standardError = standardError;
        this.confidenceLevel = confidenceLevel;
        this.lower = lower;
        this.upper = upper;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return test.toString() + String.format(Locale.ROOT, "%n  estimate = %.6g, standard error = %.6g%n"
                + "  %.4g confidence interval = [%.6g, %.6g]", estimate, standardError, confidenceLevel, lower, upper);
    }
}
