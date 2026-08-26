package math.stats;

import java.util.Locale;

/**
 * The outcome of an F-test: the {@link TestResult} every test returns, plus the
 * two degrees of freedom its null distribution needs.
 * <p>
 * {@link TestResult#degreesOfFreedom} is {@link Double#NaN} for these tests --
 * one number does not describe an F distribution -- so the pair lives here.
 * <p>
 * They are {@code double} for the same reason
 * {@link TestResult#degreesOfFreedom} is: Welch's approximation to a k-sample
 * test of equal means puts a fractional value in the denominator, and every
 * other test here happens to put a whole number there.
 *
 * @since 1.5.3
 */
public final class FTestResult {

    /** The statistic, its p-value and the alternative. */
    public final TestResult test;
    /** Degrees of freedom of the numerator. */
    public final double numeratorDf;
    /** Degrees of freedom of the denominator, fractional for Welch's test. */
    public final double denominatorDf;

    FTestResult(TestResult test, double numeratorDf, double denominatorDf) {
        this.test = test;
        this.numeratorDf = numeratorDf;
        this.denominatorDf = denominatorDf;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        // the statistic is already on the line above; only the pair of
        // degrees of freedom is missing from it
        return test.toString() + String.format(Locale.ROOT, "%n  df = (%s, %s)", format(numeratorDf),
                format(denominatorDf));
    }

    /**
     * One degrees of freedom in the shape {@link TestResult#toString()} prints
     * one: as a whole number where it is one, since {@code (2, 27)} reads better
     * than {@code (2.00000, 27.0000)} and Welch's value is the only one here
     * that is not whole.
     */
    private static String format(double df) {
        if (df == Math.rint(df) && Math.abs(df) < 1.0e15) {
            return Long.toString((long) df);
        }
        return String.format(Locale.ROOT, "%.6g", df);
    }
}
