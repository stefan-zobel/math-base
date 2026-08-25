package math.stats;

import java.util.Locale;

/**
 * The outcome of an F-test: the {@link TestResult} every test returns, plus the
 * two degrees of freedom its null distribution needs. The statistic there is
 * the ratio of the two sample variances, in the order the samples were given.
 * <p>
 * {@link TestResult#degreesOfFreedom} is {@link Double#NaN} for these tests --
 * one number does not describe an F distribution -- so the pair lives here.
 *
 * @since 1.5.3
 */
public final class FTestResult {

    /** The statistic, its p-value and the alternative. */
    public final TestResult test;
    /** Degrees of freedom of the numerator. */
    public final int numeratorDf;
    /** Degrees of freedom of the denominator. */
    public final int denominatorDf;

    FTestResult(TestResult test, int numeratorDf, int denominatorDf) {
        this.test = test;
        this.numeratorDf = numeratorDf;
        this.denominatorDf = denominatorDf;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        // the statistic is the variance ratio, so it is already on the line
        // above; only the pair of degrees of freedom is missing from it
        return test.toString() + String.format(Locale.ROOT, "%n  df = (%d, %d)", numeratorDf, denominatorDf);
    }
}
