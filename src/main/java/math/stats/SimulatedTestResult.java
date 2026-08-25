package math.stats;

import java.util.Locale;

/**
 * The result of a test whose null distribution was drawn rather than computed:
 * the test itself, and how much of its p-value is simulation noise.
 *
 * @since 1.5.3
 */
public final class SimulatedTestResult {

    /** The statistic, its p-value and what was fitted. */
    public final TestResult test;

    /** How many samples the null distribution was drawn from. */
    public final int replications;

    /**
     * The standard error of {@link TestResult#pValue} as an estimate,
     * {@code sqrt(p (1 - p) / replications)}. It is the distance the p-value
     * would be expected to move if the test were run again from another seed,
     * and it shrinks only as the square root of the work.
     */
    public final double monteCarloStandardError;

    SimulatedTestResult(TestResult test, int replications, double monteCarloStandardError) {
        this.test = test;
        this.replications = replications;
        this.monteCarloStandardError = monteCarloStandardError;
    }

    @Override
    public String toString() {
        return String.format(Locale.ROOT, "%s, %d replications (+- %.2g)", test, replications,
                monteCarloStandardError);
    }
}
