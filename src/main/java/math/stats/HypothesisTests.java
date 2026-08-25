package math.stats;

import java.util.Arrays;

import math.MathConsts;
import math.cern.ProbabilityFuncs;
import math.distribution.ContinuousDistribution;
import math.distribution.FisherF;
import math.distribution.Hypergeometric;
import math.distribution.StudentT;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.list.DoubleList;
import math.stats.gof.AndersonDarling;
import math.stats.gof.DurbinWatson;
import math.stats.gof.KolmogorovSmirnov;
import math.stats.gof.KolmogorovSmirnovPlus;
import math.stats.gof.KolmogorovSmirnovTwoSample;
import math.stats.gof.Lilliefors;

/**
 * Classical hypothesis tests, fronting the distributions in
 * {@code math.distribution} the way {@code math.solve.Quadrature} fronts the
 * integrators.
 * <p>
 * Every method returns a {@link TestResult}, or a type holding one where the
 * test yields more than a statistic and a p-value. A test that cannot be
 * computed throws {@link IllegalArgumentException} naming what was wrong with
 * the sample; none of them returns {@link Double#NaN} to mean "no answer".
 *
 * @since 1.5.3
 */
public final class HypothesisTests {

    /**
     * Tests whether the mean of {@code x} equals {@code mu}, assuming the
     * observations are independent and roughly normal, and returns a
     * confidence interval for that mean.
     * <p>
     * The statistic is {@code (mean(x) - mu) / (s / sqrt(n))} against a
     * {@link StudentT} on {@code n - 1} degrees of freedom, where {@code s} is
     * the sample standard deviation. The interval is built to match
     * {@code alternative}, so it has one infinite end for a one-sided test and
     * excludes {@code mu} exactly when the test rejects at
     * {@code 1 - confidenceLevel}.
     *
     * @param x
     *            the sample, at least two finite observations
     * @param mu
     *            the mean under the null hypothesis
     * @param alternative
     *            which departure from {@code mu} to look for
     * @param confidenceLevel
     *            the level of the returned interval, strictly between
     *            {@code 0} and {@code 1}
     * @return the statistic, its p-value, the estimate and its interval
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, holds fewer than two
     *             observations, holds a value that is not finite, or has no
     *             spread at all; if {@code mu} is not finite; if
     *             {@code alternative} is {@code null}; or if
     *             {@code confidenceLevel} does not lie in {@code (0, 1)}
     */
    public static TTestResult tOneSample(double[] x, double mu, Alternative alternative, double confidenceLevel) {
        int n = requireFiniteSample(x, "x");
        requireAtLeastTwo(n, "x");
        if (!isFinite(mu)) {
            throw new IllegalArgumentException("mu must be finite : " + mu);
        }
        requireAlternative(alternative);
        requireLevel(confidenceLevel);

        double[] moments = meanAndDeviation(x);
        double mean = moments[0];
        // s / sqrt(n) rather than sqrt(s^2 / n): a variance leaves the double
        // range where a deviation, and the scale invariant statistic built from
        // it, are still perfectly representable
        double standardError = moments[1] / Math.sqrt(n);
        if (!(standardError > 0.0)) {
            throw new IllegalArgumentException("x has no spread: every observation is " + x[0]);
        }

        return tTest("one-sample t", mean, mu, standardError, n - 1.0, alternative, confidenceLevel);
    }

    /**
     * Tests whether the means of two independent samples are equal, without
     * assuming that their variances are, and returns a confidence interval for
     * the difference.
     * <p>
     * This is Welch's test, and it is the default on purpose: the pooled
     * alternative buys a little power for an assumption of equal variances that
     * callers rarely check and that Welch does not need. Use
     * {@link #tTwoSamplePooled(double[], double[], Alternative, double)} where
     * the assumption is actually justified, so that the choice is visible at the
     * call site.
     * <p>
     * The estimate is {@code mean(x) - mean(y)}, and the alternative is read in
     * that direction: {@link Alternative#GREATER} asks whether the mean of
     * {@code x} exceeds the mean of {@code y}. The degrees of freedom come from
     * the Welch-Satterthwaite approximation and are not a whole number.
     *
     * @param x
     *            the first sample, at least two finite observations
     * @param y
     *            the second sample, at least two finite observations
     * @param alternative
     *            which departure from equality to look for
     * @param confidenceLevel
     *            the level of the returned interval, strictly between
     *            {@code 0} and {@code 1}
     * @return the statistic, its p-value, the difference of the means and an
     *         interval for it
     * @throws IllegalArgumentException
     *             if either sample is {@code null}, holds fewer than two
     *             observations or holds a value that is not finite; if neither
     *             sample has any spread; if {@code alternative} is {@code null};
     *             or if {@code confidenceLevel} does not lie in {@code (0, 1)}
     */
    public static TTestResult tTwoSample(double[] x, double[] y, Alternative alternative,
            double confidenceLevel) {
        int nx = requireFiniteSample(x, "x");
        requireAtLeastTwo(nx, "x");
        int ny = requireFiniteSample(y, "y");
        requireAtLeastTwo(ny, "y");
        requireAlternative(alternative);
        requireLevel(confidenceLevel);

        double[] mx = meanAndDeviation(x);
        double[] my = meanAndDeviation(y);
        // the standard error of each mean, never the variance: hypot is the
        // overflow-safe way to add them in quadrature
        double ex = mx[1] / Math.sqrt(nx);
        double ey = my[1] / Math.sqrt(ny);
        double standardError = Math.hypot(ex, ey);
        if (!(standardError > 0.0)) {
            throw new IllegalArgumentException("neither sample has any spread");
        }
        return tTest("Welch two-sample t", mx[0] - my[0], 0.0, standardError,
                welchDegreesOfFreedom(ex, ey, nx, ny), alternative, confidenceLevel);
    }

    /**
     * Tests whether the means of two independent samples are equal, assuming
     * their variances are, and returns a confidence interval for the
     * difference.
     * <p>
     * The assumption is in the name because it is an assumption. Where it does
     * not hold this test is not conservative: it borrows degrees of freedom from
     * the larger sample and reports a p-value that is too small.
     * {@link #tTwoSample(double[], double[], Alternative, double)} is the one to
     * reach for by default.
     *
     * @param x
     *            the first sample, at least two finite observations
     * @param y
     *            the second sample, at least two finite observations
     * @param alternative
     *            which departure from equality to look for
     * @param confidenceLevel
     *            the level of the returned interval, strictly between
     *            {@code 0} and {@code 1}
     * @return the statistic, its p-value, the difference of the means and an
     *         interval for it
     * @throws IllegalArgumentException
     *             under the same conditions as
     *             {@link #tTwoSample(double[], double[], Alternative, double)}
     */
    public static TTestResult tTwoSamplePooled(double[] x, double[] y, Alternative alternative,
            double confidenceLevel) {
        int nx = requireFiniteSample(x, "x");
        requireAtLeastTwo(nx, "x");
        int ny = requireFiniteSample(y, "y");
        requireAtLeastTwo(ny, "y");
        requireAlternative(alternative);
        requireLevel(confidenceLevel);

        double[] mx = meanAndDeviation(x);
        double[] my = meanAndDeviation(y);
        double largest = Math.max(mx[1], my[1]);
        if (!(largest > 0.0)) {
            throw new IllegalArgumentException("neither sample has any spread");
        }
        // square deviations that have been scaled to at most one and put the
        // scale back afterwards, so a spread beyond 1e154 still has a pooled
        // standard deviation
        double p = mx[1] / largest;
        double q = my[1] / largest;
        double pooled = largest * Math.sqrt(((nx - 1.0) * p * p + (ny - 1.0) * q * q) / (nx + ny - 2.0));
        double standardError = pooled * Math.sqrt(1.0 / nx + 1.0 / ny);
        return tTest("pooled two-sample t", mx[0] - my[0], 0.0, standardError, nx + ny - 2.0, alternative,
                confidenceLevel);
    }

    /**
     * Tests whether the mean of the differences between two paired samples is
     * zero, and returns a confidence interval for it.
     * <p>
     * This is {@link #tOneSample(double[], double, Alternative, double)} on
     * {@code x[i] - y[i]} against zero, and it is that literally: the numbers it
     * returns agree with that call bit for bit. Pairing is what makes it a
     * different test from
     * {@link #tTwoSample(double[], double[], Alternative, double)} -- it removes
     * whatever the two observations of a pair have in common, which is usually
     * most of the variance.
     *
     * @param x
     *            the first of the two paired samples, at least two finite
     *            observations
     * @param y
     *            the second, of the same length as {@code x}
     * @param alternative
     *            which departure from zero to look for
     * @param confidenceLevel
     *            the level of the returned interval, strictly between
     *            {@code 0} and {@code 1}
     * @return the statistic, its p-value, the mean difference and an interval
     *         for it
     * @throws IllegalArgumentException
     *             if either sample is {@code null}, they are of different
     *             lengths, they hold fewer than two observations, a value or a
     *             difference is not finite, the differences have no spread,
     *             {@code alternative} is {@code null}, or
     *             {@code confidenceLevel} does not lie in {@code (0, 1)}
     */
    public static TTestResult tPaired(double[] x, double[] y, Alternative alternative, double confidenceLevel) {
        int n = requireFiniteSample(x, "x");
        requireFiniteSample(y, "y");
        if (y.length != n) {
            throw new IllegalArgumentException(
                    "x and y must be paired, but their lengths are " + n + " and " + y.length);
        }
        requireAtLeastTwo(n, "x");
        requireAlternative(alternative);
        requireLevel(confidenceLevel);

        double[] differences = new double[n];
        for (int i = 0; i < n; i++) {
            differences[i] = x[i] - y[i];
            if (!isFinite(differences[i])) {
                throw new IllegalArgumentException(
                        "x[" + i + "] - y[" + i + "] is not a finite number : " + differences[i]);
            }
        }
        double[] moments = meanAndDeviation(differences);
        double standardError = moments[1] / Math.sqrt(n);
        if (!(standardError > 0.0)) {
            throw new IllegalArgumentException(
                    "the differences have no spread: every one of them is " + differences[0]);
        }
        return tTest("paired t", moments[0], 0.0, standardError, n - 1.0, alternative, confidenceLevel);
    }

    /**
     * Tests whether two independent samples have the same variance, and
     * returns the ratio of their sample variances.
     * <p>
     * The statistic is {@code var(x) / var(y)} against an {@link FisherF} on
     * {@code (nx - 1, ny - 1)} degrees of freedom, and the alternative is read
     * in that direction: {@link Alternative#GREATER} asks whether {@code x} is
     * the more variable of the two.
     * <p>
     * The test is worth using only where normality is credible. It is not
     * robust to a heavy tail -- a departure from normality moves its level far
     * more than a departure from equal variances does -- and where the question
     * is really "may I pool", the answer is to use
     * {@link #tTwoSample(double[], double[], Alternative, double)} and not ask.
     *
     * @param x
     *            the first sample, at least two finite observations
     * @param y
     *            the second sample, at least two finite observations
     * @param alternative
     *            which departure from equal variances to look for
     * @return the statistic, its p-value and the two degrees of freedom
     * @throws IllegalArgumentException
     *             if either sample is {@code null}, holds fewer than two
     *             observations or holds a value that is not finite; if
     *             {@code y} has no spread, which would put a zero in the
     *             denominator; or if {@code alternative} is {@code null}
     */
    public static FTestResult fVarianceRatio(double[] x, double[] y, Alternative alternative) {
        int nx = requireFiniteSample(x, "x");
        requireAtLeastTwo(nx, "x");
        int ny = requireFiniteSample(y, "y");
        requireAtLeastTwo(ny, "y");
        requireAlternative(alternative);

        double sx = meanAndDeviation(x)[1];
        double sy = meanAndDeviation(y)[1];
        if (!(sy > 0.0)) {
            throw new IllegalArgumentException("y has no spread, so the variance ratio is not defined");
        }
        // the ratio of the deviations first and the square afterwards: the two
        // variances on their own leave the double range from a spread of 1e154
        // upwards, where their quotient is still an ordinary number
        double ratio = sx / sy;
        double statistic = ratio * ratio;
        int numeratorDf = nx - 1;
        int denominatorDf = ny - 1;

        double lower = new FisherF(numeratorDf, denominatorDf).cdf(statistic);
        // the upper tail through the reciprocal rather than as 1 - cdf: if
        // F ~ F(a, b) then 1/F ~ F(b, a), and the subtraction would throw away
        // every significant digit of a small p-value
        double upper = (statistic > 0.0) ? new FisherF(denominatorDf, numeratorDf).cdf(1.0 / statistic) : 1.0;
        double pValue;
        switch (alternative) {
        case LESS:
            pValue = lower;
            break;
        case GREATER:
            pValue = upper;
            break;
        default:
            pValue = Math.min(1.0, 2.0 * Math.min(lower, upper));
            break;
        }
        TestResult test = new TestResult("F variance ratio", statistic, pValue, alternative, Double.NaN);
        return new FTestResult(test, numeratorDf, denominatorDf);
    }

    /**
     * Tests whether counts follow a distribution given as expected counts.
     * <p>
     * The statistic is {@code sum (observed - expected)^2 / expected} against a
     * {@link math.distribution.ChiSquare} on
     * {@code categories - 1 - estimatedParameters} degrees of freedom, and its
     * p-value is the upper tail: a poor fit makes the statistic large. The
     * result therefore reports {@link Alternative#GREATER}, which is a statement
     * about the statistic and not about a parameter.
     * <p>
     * {@code estimatedParameters} is not decoration. Fitting the expected
     * counts to the same sample they are tested against uses up degrees of
     * freedom, and a test that does not subtract them reports a p-value that is
     * too large -- it fails to reject a model it should. Pass {@code 0} only
     * when the expected counts came from somewhere other than this sample.
     * <p>
     * The approximation is a large-sample one. The usual rule of thumb asks for
     * an expected count of at least five in every category; this method does
     * not enforce it, because the rule is a rule of thumb and the right answer
     * to a sparse table is often to pool categories rather than to give up.
     *
     * @param observed
     *            the counts that were seen, at least two categories, none
     *            negative
     * @param expected
     *            the counts the null hypothesis predicts, as many as there are
     *            observed, every one strictly positive, summing to the same
     *            total as {@code observed}
     * @param estimatedParameters
     *            how many parameters of the null distribution were estimated
     *            from this same sample
     * @return the statistic, its p-value and the degrees of freedom
     * @throws IllegalArgumentException
     *             if either array is {@code null}, they are of different
     *             lengths, there are fewer than two categories, an observed
     *             count is negative or not finite, an expected count is not
     *             strictly positive and finite, the two totals disagree,
     *             {@code estimatedParameters} is negative, or nothing is left
     *             of the degrees of freedom
     */
    public static TestResult chiSquaredGoodnessOfFit(long[] observed, double[] expected,
            int estimatedParameters) {
        if (observed == null) {
            throw new IllegalArgumentException("observed must not be null");
        }
        if (expected == null) {
            throw new IllegalArgumentException("expected must not be null");
        }
        if (observed.length != expected.length) {
            throw new IllegalArgumentException("observed and expected must match, but their lengths are "
                    + observed.length + " and " + expected.length);
        }
        int categories = observed.length;
        if (categories < 2) {
            throw new IllegalArgumentException("at least two categories are needed, got " + categories);
        }
        if (estimatedParameters < 0) {
            throw new IllegalArgumentException("estimatedParameters must not be negative : "
                    + estimatedParameters);
        }
        double degreesOfFreedom = categories - 1.0 - estimatedParameters;
        if (degreesOfFreedom < 1.0) {
            throw new IllegalArgumentException("estimating " + estimatedParameters + " parameters from "
                    + categories + " categories leaves " + degreesOfFreedom + " degrees of freedom");
        }

        double observedTotal = 0.0;
        double expectedTotal = 0.0;
        double statistic = 0.0;
        for (int i = 0; i < categories; i++) {
            if (observed[i] < 0L) {
                throw new IllegalArgumentException("observed[" + i + "] must not be negative : " + observed[i]);
            }
            if (!(expected[i] > 0.0) || !isFinite(expected[i])) {
                throw new IllegalArgumentException(
                        "expected[" + i + "] must be strictly positive and finite : " + expected[i]);
            }
            observedTotal += observed[i];
            expectedTotal += expected[i];
            double residual = observed[i] - expected[i];
            statistic += residual * residual / expected[i];
        }
        if (Math.abs(observedTotal - expectedTotal) > TOTALS_TOLERANCE * Math.max(observedTotal, 1.0)) {
            throw new IllegalArgumentException("expected holds counts, not probabilities: observed sums to "
                    + observedTotal + " and expected to " + expectedTotal);
        }
        return chiSquaredResult("chi-squared goodness of fit", statistic, degreesOfFreedom);
    }

    /**
     * Tests whether the row and the column of a contingency table are
     * independent.
     * <p>
     * The expected count of a cell is the product of its margins over the grand
     * total, the statistic is {@code sum (observed - expected)^2 / expected},
     * and the degrees of freedom are {@code (rows - 1) * (columns - 1)}. As for
     * the goodness of fit test the p-value is the upper tail, and the caveat
     * about small expected counts applies here too -- more so, because a table
     * has more cells to spread the same total over.
     *
     * @param table
     *            the counts, at least two rows of at least two columns, all of
     *            the same length, none of them negative, no margin empty
     * @return the statistic, its p-value and the degrees of freedom
     * @throws IllegalArgumentException
     *             if the table is {@code null}, holds a {@code null} row, is
     *             ragged, is smaller than two by two, holds a negative count,
     *             or has an empty row or column
     */
    public static TestResult chiSquaredIndependence(long[][] table) {
        if (table == null) {
            throw new IllegalArgumentException("table must not be null");
        }
        int rows = table.length;
        if (rows < 2) {
            throw new IllegalArgumentException("at least two rows are needed, got " + rows);
        }
        if (table[0] == null) {
            throw new IllegalArgumentException("table[0] must not be null");
        }
        int columns = table[0].length;
        if (columns < 2) {
            throw new IllegalArgumentException("at least two columns are needed, got " + columns);
        }

        double[] rowTotals = new double[rows];
        double[] columnTotals = new double[columns];
        double total = 0.0;
        for (int i = 0; i < rows; i++) {
            if (table[i] == null) {
                throw new IllegalArgumentException("table[" + i + "] must not be null");
            }
            if (table[i].length != columns) {
                throw new IllegalArgumentException("table[" + i + "] has " + table[i].length
                        + " columns where table[0] has " + columns);
            }
            for (int k = 0; k < columns; k++) {
                long count = table[i][k];
                if (count < 0L) {
                    throw new IllegalArgumentException(
                            "table[" + i + "][" + k + "] must not be negative : " + count);
                }
                rowTotals[i] += count;
                columnTotals[k] += count;
                total += count;
            }
        }
        for (int i = 0; i < rows; i++) {
            if (!(rowTotals[i] > 0.0)) {
                throw new IllegalArgumentException("row " + i + " is empty, so it has no expected counts");
            }
        }
        for (int k = 0; k < columns; k++) {
            if (!(columnTotals[k] > 0.0)) {
                throw new IllegalArgumentException("column " + k + " is empty, so it has no expected counts");
            }
        }

        double statistic = 0.0;
        for (int i = 0; i < rows; i++) {
            for (int k = 0; k < columns; k++) {
                // the margins over the grand total, grouped so that a large
                // table cannot overflow the product on its way
                double expected = (rowTotals[i] / total) * columnTotals[k];
                double residual = table[i][k] - expected;
                statistic += residual * residual / expected;
            }
        }
        return chiSquaredResult("chi-squared independence", statistic, (rows - 1.0) * (columns - 1.0));
    }

    /**
     * Tests a two by two table for independence exactly, conditioning on both
     * margins.
     * <p>
     * With the margins held fixed a single cell determines the whole table, and
     * under independence that cell follows a {@link Hypergeometric}
     * distribution, so the p-value is a sum over the null rather than an
     * approximation to it. That is what makes this the test to reach for on the
     * small and sparse tables where
     * {@link #chiSquaredIndependence(long[][])} is not to be trusted -- it needs
     * no expected count to be large, and it needs no continuity correction.
     * <p>
     * Two properties of an exact test of a discrete statistic are worth knowing
     * before the result is read. It is <b>conservative</b>: only finitely many
     * p-values are attainable, so the rejection rate at a given level sits below
     * that level rather than on it, measurably so -- at {@code alpha = 0.05} and
     * a total of 20 it is near 0.019. And it costs a sum whose length is the
     * smaller of the two margins, so a table of a million counts takes
     * milliseconds where the chi-squared test takes none.
     * <p>
     * {@link Alternative#GREATER} sums the upper tail from the observed value of
     * {@code table[0][0]} and so asks whether that cell is <em>larger</em> than
     * independence predicts, {@link Alternative#LESS} sums the lower tail, and
     * {@link Alternative#TWO_SIDED} sums the probability of every table no more
     * likely than the one observed. The last is the convention that leaves the
     * answer unchanged when the table is transposed or its columns are
     * exchanged, so it does not matter which cell is taken to carry the test.
     *
     * @param table
     *            the counts, two rows of two columns, none of them negative, no
     *            margin empty and no count and no total above
     *            {@link Integer#MAX_VALUE}
     * @param alternative
     *            which departure from independence to look for
     * @return the count in {@code table[0][0]} as the statistic -- given the
     *         margins it is the whole table -- its p-value, and no degrees of
     *         freedom
     * @throws IllegalArgumentException
     *             if {@code table} is {@code null}, is not two by two, holds a
     *             {@code null} row or a negative count, has an empty row or
     *             column, or counts past {@link Integer#MAX_VALUE}; or if
     *             {@code alternative} is {@code null}
     */
    public static TestResult fisherExact(long[][] table, Alternative alternative) {
        requireAlternative(alternative);
        if (table == null) {
            throw new IllegalArgumentException("table must not be null");
        }
        if (table.length != 2 || table[0] == null || table[0].length != 2 || table[1] == null
                || table[1].length != 2) {
            throw new IllegalArgumentException(
                    "the exact test needs a two by two table; larger ones go to chiSquaredIndependence");
        }

        long total = 0L;
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                long count = table[i][k];
                if (count < 0L) {
                    throw new IllegalArgumentException(
                            "table[" + i + "][" + k + "] must not be negative : " + count);
                }
                if (count > Integer.MAX_VALUE) {
                    throw new IllegalArgumentException(
                            "table[" + i + "][" + k + "] is past Integer.MAX_VALUE : " + count);
                }
                total += count;
            }
        }
        if (total > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("the grand total is past Integer.MAX_VALUE : " + total);
        }
        long firstRow = table[0][0] + table[0][1];
        long firstColumn = table[0][0] + table[1][0];
        if (firstRow == 0L || firstRow == total) {
            throw new IllegalArgumentException("a row is empty, so the margins already fix the table");
        }
        if (firstColumn == 0L || firstColumn == total) {
            throw new IllegalArgumentException("a column is empty, so the margins already fix the table");
        }

        int a = (int) table[0][0];
        Hypergeometric nullDistribution = new Hypergeometric((int) total, (int) firstColumn, (int) firstRow);
        int lo = nullDistribution.supportLowerBound();
        int hi = nullDistribution.supportUpperBound();

        // every tail is summed from its far end inwards, never as one minus the
        // other: the terms then grow as the sum grows, so no small one is lost
        // under a large one, and a tail of 1e-300 comes back as itself rather
        // than as the round-off of a subtraction from one
        double pValue = 0.0;
        switch (alternative) {
        case LESS:
            for (int k = lo; k <= a; k++) {
                pValue += nullDistribution.pmf(k);
            }
            break;
        case GREATER:
            for (int k = hi; k >= a; k--) {
                pValue += nullDistribution.pmf(k);
            }
            break;
        default:
            // the mirror image of the observed table has the same probability in
            // exact arithmetic but not in floating point, and dropping it would
            // halve the p-value, which is what the tolerance is for
            double observed = nullDistribution.pmf(a) * (1.0 + LIKELIHOOD_TOLERANCE);
            for (int k = lo; k <= hi; k++) {
                double mass = nullDistribution.pmf(k);
                if (mass <= observed) {
                    pValue += mass;
                }
            }
            break;
        }
        return new TestResult("Fisher exact", a, Math.min(1.0, pValue), alternative, Double.NaN);
    }

    /**
     * Tests whether a sample was drawn from a given continuous distribution.
     * <p>
     * The sample is carried into {@code (0, 1)} by the distribution function --
     * if the hypothesis holds, the transformed values are uniform there -- and
     * the statistic is the largest gap between their empirical distribution
     * function and the diagonal. The p-value is exact for
     * {@link Alternative#TWO_SIDED} up to at least ten decimal digits for
     * {@code n <= 500} and six up to {@code n = 200000}.
     * <p>
     * The alternatives select the one-sided statistics and are stated about the
     * <em>sample</em>, not about its distribution function, which is the
     * opposite of one convention in the field and worth reading twice:
     * {@link Alternative#GREATER} asks whether the sample is stochastically
     * <em>larger</em> than {@code nullDistribution}, which shows up as its
     * distribution function lying <em>below</em> the hypothesized one.
     * <p>
     * <b>The distribution has to be fully specified before the sample is
     * seen.</b> Estimating its parameters from the same data and testing
     * against the result is a different test with a different null
     * distribution -- Lilliefors -- and the p-value this method returns is then
     * much too large. {@code math.stats.mle} sitting in the next package makes
     * that mistake easy to make; a fit of the sample is not a hypothesis about
     * it.
     *
     * @param x
     *            the sample, at least one finite observation
     * @param nullDistribution
     *            the distribution the sample is tested against, with every
     *            parameter fixed independently of {@code x}
     * @param alternative
     *            which departure to look for
     * @return the statistic, its p-value and no degrees of freedom, the null
     *         distribution of {@code D_n} having none
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is empty or holds a value that
     *             is not finite; or if {@code nullDistribution} or
     *             {@code alternative} is {@code null}
     */
    public static TestResult kolmogorovSmirnov(double[] x, ContinuousDistribution nullDistribution,
            Alternative alternative) {
        int n = requireFiniteSample(x, "x");
        if (nullDistribution == null) {
            throw new IllegalArgumentException("nullDistribution must not be null");
        }
        requireAlternative(alternative);

        double[] uniform = new double[n];
        for (int i = 0; i < n; i++) {
            uniform[i] = nullDistribution.cdf(x[i]);
        }
        Arrays.sort(uniform);

        // D+ is how far the empirical distribution function rises above the
        // diagonal and D- how far it falls below it. A sample shifted upwards
        // has large transformed values, which is a large D-
        double dPlus = 0.0;
        double dMinus = 0.0;
        for (int i = 0; i < n; i++) {
            dPlus = Math.max(dPlus, (i + 1.0) / n - uniform[i]);
            dMinus = Math.max(dMinus, uniform[i] - i / (double) n);
        }

        double statistic;
        double pValue;
        switch (alternative) {
        case LESS:
            statistic = dPlus;
            pValue = KolmogorovSmirnovPlus.barF(n, statistic);
            break;
        case GREATER:
            statistic = dMinus;
            pValue = KolmogorovSmirnovPlus.barF(n, statistic);
            break;
        default:
            statistic = Math.max(dPlus, dMinus);
            pValue = KolmogorovSmirnov.barF(n, statistic);
            break;
        }
        return new TestResult("Kolmogorov-Smirnov", statistic, Math.min(1.0, Math.max(0.0, pValue)),
                alternative, Double.NaN);
    }

    /**
     * Tests whether a sample was drawn from a given continuous distribution,
     * weighting the tails.
     * <p>
     * The sample is carried into {@code (0, 1)} by the distribution function
     * the way {@link #kolmogorovSmirnov} carries it, but the statistic is the
     * integrated squared distance divided by {@code u (1 - u)} rather than the
     * largest gap:
     * <p>
     * {@code A_n^2 = -n - (1/n) sum (2i - 1) (ln u_i + ln(1 - u_(n+1-i)))}
     * <p>
     * That divisor is the whole point. {@code D_n} is a supremum of an
     * unweighted difference, and near the ends of the range there is almost no
     * room for a difference to be large, so a sample that departs only in its
     * tails moves {@code D_n} very little. Measured on a unit exponential
     * contaminated with five percent of a five times longer one, this test
     * rejects 24 percent of samples of 200 at the five percent level where
     * {@link #kolmogorovSmirnov} rejects 11 percent. The advantage is not
     * universal -- for a standardized {@code t(5)} at {@code n <= 100} it was
     * measured slightly the other way -- but where the two disagree it is
     * usually this one that sees the departure.
     * <p>
     * There is no {@code alternative} parameter: the statistic squares the
     * difference, so it has no one-sided form, and the result carries
     * {@link Alternative#TWO_SIDED}.
     * <p>
     * <b>The distribution has to be fully specified before the sample is
     * seen</b>, for the reason {@link #kolmogorovSmirnov} gives at more length;
     * {@link #andersonDarling(double[], Lilliefors.Family)} is the test for a
     * distribution fitted to the same sample. Against a fitted normal this
     * method rejected none of 6000 samples that really were normal.
     * <p>
     * The p-value comes from {@link AndersonDarling#barF}, which is an
     * interpolation of a numerically integrated asymptotic distribution with an
     * empirical correction in {@code 1/n}, good to about three decimal digits.
     * It is not exact the way the Kolmogorov-Smirnov p-value is.
     *
     * @param x
     *            the sample, at least one finite observation
     * @param nullDistribution
     *            the distribution the sample is tested against, with every
     *            parameter fixed independently of {@code x}
     * @return the statistic, its p-value and no degrees of freedom, the null
     *         distribution of {@code A_n^2} having none
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is empty or holds a value that
     *             is not finite; or if {@code nullDistribution} is {@code null}
     */
    public static TestResult andersonDarling(double[] x, ContinuousDistribution nullDistribution) {
        int n = requireFiniteSample(x, "x");
        if (nullDistribution == null) {
            throw new IllegalArgumentException("nullDistribution must not be null");
        }

        double[] uniform = new double[n];
        for (int i = 0; i < n; i++) {
            uniform[i] = nullDistribution.cdf(x[i]);
        }
        Arrays.sort(uniform);

        // the i-th term pairs the i-th smallest transformed value with the i-th
        // largest, which is what puts the weight on both tails at once
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            double lower = uniform[i];
            double upper = 1.0 - uniform[n - 1 - i];
            if (lower < UNIFORM_FLOOR) {
                lower = UNIFORM_FLOOR;
            }
            if (upper < UNIFORM_FLOOR) {
                upper = UNIFORM_FLOOR;
            }
            sum += (2 * i + 1) * (Math.log(lower) + Math.log(upper));
        }

        double statistic = -n - sum / n;
        double pValue = AndersonDarling.barF(n, statistic);
        return new TestResult("Anderson-Darling", statistic, Math.min(1.0, Math.max(0.0, pValue)),
                Alternative.TWO_SIDED, Double.NaN);
    }

    /**
     * Tests whether two samples were drawn from the same continuous
     * distribution, naming neither of them.
     * <p>
     * The statistic is the largest gap between the two empirical distribution
     * functions. Under the null hypothesis every interleaving of the two
     * samples is equally likely, which is what makes the p-value computable
     * without knowing the distribution at all -- the price is that a shift and
     * a change of shape are the same finding here, and the test says only that
     * the two differ.
     * <p>
     * The p-value is <b>exact</b> when the pooled sample has no repeated value
     * and {@code x.length * y.length} is at most
     * {@link KolmogorovSmirnovTwoSample#EXACT_LIMIT}, and the Kolmogorov limit
     * otherwise; {@link TestResult#test} says which of the two produced the
     * number, since nothing else in the result would show it. Ties force the
     * approximation because the exact distribution counts interleavings of
     * distinct values, and the limit is conservative in both cases: measured
     * against the exact value it always came out too large.
     * <p>
     * The alternatives are stated about {@code x}, in the same direction the
     * one-sample test uses: {@link Alternative#GREATER} asks whether {@code x}
     * is stochastically <em>larger</em> than {@code y}, which shows up as its
     * distribution function lying <em>below</em> the other one.
     *
     * @param x
     *            the first sample, at least one finite observation
     * @param y
     *            the second sample, at least one finite observation
     * @param alternative
     *            which departure to look for
     * @return the statistic, its p-value and no degrees of freedom, the null
     *         distribution of {@code D_(m,n)} having none
     * @throws IllegalArgumentException
     *             if {@code x} or {@code y} is {@code null}, is empty or holds a
     *             value that is not finite, or if {@code alternative} is
     *             {@code null}
     */
    public static TestResult kolmogorovSmirnovTwoSample(double[] x, double[] y, Alternative alternative) {
        int m = requireFiniteSample(x, "x");
        int n = requireFiniteSample(y, "y");
        requireAlternative(alternative);

        double[] first = x.clone();
        double[] second = y.clone();
        Arrays.sort(first);
        Arrays.sort(second);

        // walk the pooled sample once, stepping past every repetition of a
        // value before reading the two empirical distribution functions off:
        // between the repetitions neither function is defined yet
        double dPlus = 0.0;
        double dMinus = 0.0;
        boolean tied = false;
        int i = 0;
        int k = 0;
        while (i < m && k < n) {
            double value = Math.min(first[i], second[k]);
            int consumed = 0;
            while (i < m && first[i] == value) {
                i++;
                consumed++;
            }
            while (k < n && second[k] == value) {
                k++;
                consumed++;
            }
            if (consumed > 1) {
                tied = true;
            }
            double gap = i / (double) m - k / (double) n;
            dPlus = Math.max(dPlus, gap);
            dMinus = Math.max(dMinus, -gap);
        }
        // whichever sample is left over only closes the gap again, so the
        // maximum has already been seen -- but its repetitions have not
        tied = tied || hasARepetition(first, i) || hasARepetition(second, k);

        boolean exact = !tied && (long) m * n <= KolmogorovSmirnovTwoSample.EXACT_LIMIT;
        String name = exact ? "two-sample Kolmogorov-Smirnov, exact"
                : "two-sample Kolmogorov-Smirnov, asymptotic";
        double statistic;
        double pValue;
        switch (alternative) {
        case LESS:
            statistic = dPlus;
            pValue = exact ? KolmogorovSmirnovTwoSample.barFExactOneSided(m, n, statistic)
                    : KolmogorovSmirnovTwoSample.barFAsymptoticOneSided(m, n, statistic);
            break;
        case GREATER:
            statistic = dMinus;
            pValue = exact ? KolmogorovSmirnovTwoSample.barFExactOneSided(m, n, statistic)
                    : KolmogorovSmirnovTwoSample.barFAsymptoticOneSided(m, n, statistic);
            break;
        default:
            statistic = Math.max(dPlus, dMinus);
            pValue = exact ? KolmogorovSmirnovTwoSample.barFExact(m, n, statistic)
                    : KolmogorovSmirnovTwoSample.barFAsymptotic(m, n, statistic);
            break;
        }
        return new TestResult(name, statistic, Math.min(1.0, Math.max(0.0, pValue)), alternative, Double.NaN);
    }

    /**
     * Tests whether a sample was drawn from a family of distributions, with the
     * parameters of that family estimated from the sample itself.
     * <p>
     * This is the test {@link #kolmogorovSmirnov} is not. Fitting a
     * distribution to a sample and then measuring how far the sample sits from
     * the fit understates the distance -- the fit has already moved towards the
     * sample -- so the ordinary null distribution returns a p-value that is much
     * too large. Here the null is drawn instead: samples of the same size from
     * the same family, each fitted and measured the same way, which is exactly
     * the situation the real one is in.
     * <p>
     * That makes the p-value a <b>simulation</b>, with an uncertainty of its own
     * that {@link SimulatedTestResult#monteCarloStandardError} reports. The short
     * form uses {@link Lilliefors#DEFAULT_REPLICATIONS} replications and a fixed
     * seed, so the same sample always gives the same answer; the long form is
     * there to buy more accuracy, or to see how much the answer moves when the
     * seed does.
     *
     * @param x
     *            the sample, at least {@link Lilliefors#MINIMUM_SAMPLE} finite
     *            observations, strictly positive for the two families that need
     *            it
     * @param family
     *            the family to fit and test against
     * @return the statistic, its simulated p-value, and how many replications
     *         went into it
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is too short or holds a value
     *             that is not finite or not in the support of the family; or if
     *             {@code family} is {@code null}
     */
    public static SimulatedTestResult lilliefors(double[] x, Lilliefors.Family family) {
        return lilliefors(x, family, Lilliefors.DEFAULT_REPLICATIONS, Lilliefors.DEFAULT_SEED);
    }

    /**
     * Tests whether a sample was drawn from a fitted family, drawing the null
     * distribution from {@code replications} samples seeded by {@code seed}.
     * <p>
     * See {@link #lilliefors(double[], Lilliefors.Family)} for what the test
     * does. The p-value is reproducible from {@code seed} alone, whether or not
     * the replications end up spread over several threads.
     *
     * @param x
     *            the sample, at least {@link Lilliefors#MINIMUM_SAMPLE} finite
     *            observations, strictly positive for the two families that need
     *            it
     * @param family
     *            the family to fit and test against
     * @param replications
     *            how many samples the null distribution is drawn from,
     *            {@code 1} or more
     * @param seed
     *            the seed the replications are derived from
     * @return the statistic, its simulated p-value, and how many replications
     *         went into it
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is too short or holds a value
     *             that is not finite or not in the support of the family; if
     *             {@code family} is {@code null}; or if {@code replications} is
     *             not strictly positive
     */
    public static SimulatedTestResult lilliefors(double[] x, Lilliefors.Family family, int replications,
            long seed) {
        return fitted(Lilliefors.Statistic.KOLMOGOROV_SMIRNOV, "Lilliefors", x, family, replications, seed);
    }

    /**
     * Tests whether a sample was drawn from a family of distributions,
     * weighting the tails, with the parameters of that family estimated from
     * the sample itself.
     * <p>
     * This is to {@link #andersonDarling(double[], ContinuousDistribution)} what
     * {@link #lilliefors(double[], Lilliefors.Family)} is to
     * {@link #kolmogorovSmirnov}, and it is the test usually meant by "the
     * Anderson-Darling normality test": the fit has already moved towards the
     * sample, so the tabulated null distribution is far too generous and the
     * null has to be drawn instead. Against a fitted normal the fully specified
     * test rejected none of 6000 samples that really were normal.
     * <p>
     * The p-value is a <b>simulation</b>, with an uncertainty of its own that
     * {@link SimulatedTestResult#monteCarloStandardError} reports. The short form
     * uses {@link Lilliefors#DEFAULT_REPLICATIONS} replications and a fixed seed,
     * so the same sample always gives the same answer.
     *
     * @param x
     *            the sample, at least {@link Lilliefors#MINIMUM_SAMPLE} finite
     *            observations, strictly positive for the two families that need
     *            it
     * @param family
     *            the family to fit and test against
     * @return the statistic, its simulated p-value, and how many replications
     *         went into it
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is too short or holds a value
     *             that is not finite or not in the support of the family; or if
     *             {@code family} is {@code null}
     */
    public static SimulatedTestResult andersonDarling(double[] x, Lilliefors.Family family) {
        return andersonDarling(x, family, Lilliefors.DEFAULT_REPLICATIONS, Lilliefors.DEFAULT_SEED);
    }

    /**
     * Tests whether a sample was drawn from a fitted family, weighting the tails,
     * drawing the null distribution from {@code replications} samples seeded by
     * {@code seed}.
     * <p>
     * See {@link #andersonDarling(double[], Lilliefors.Family)} for what the test
     * does. The p-value is reproducible from {@code seed} alone, whether or not
     * the replications end up spread over several threads.
     *
     * @param x
     *            the sample, at least {@link Lilliefors#MINIMUM_SAMPLE} finite
     *            observations, strictly positive for the two families that need
     *            it
     * @param family
     *            the family to fit and test against
     * @param replications
     *            how many samples the null distribution is drawn from,
     *            {@code 1} or more
     * @param seed
     *            the seed the replications are derived from
     * @return the statistic, its simulated p-value, and how many replications
     *         went into it
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null}, is too short or holds a value
     *             that is not finite or not in the support of the family; if
     *             {@code family} is {@code null}; or if {@code replications} is
     *             not strictly positive
     */
    public static SimulatedTestResult andersonDarling(double[] x, Lilliefors.Family family, int replications,
            long seed) {
        return fitted(Lilliefors.Statistic.ANDERSON_DARLING, "Anderson-Darling", x, family, replications,
                seed);
    }

    /**
     * The body both fitted tests share: measure the sample against its own fit,
     * then draw the null distribution of that measurement.
     */
    private static SimulatedTestResult fitted(Lilliefors.Statistic which, String name, double[] x,
            Lilliefors.Family family, int replications, long seed) {
        int n = requireFiniteSample(x, "x");
        double statistic = Lilliefors.statistic(which, family, x);
        double pValue = Lilliefors.barF(which, family, n, statistic, replications, seed);
        TestResult test = new TestResult(name + ", " + fittedName(family) + " fitted", statistic, pValue,
                Alternative.TWO_SIDED, Double.NaN);
        return new SimulatedTestResult(test, replications,
                Math.sqrt(pValue * (1.0 - pValue) / replications));
    }

    /**
     * Tests whether the errors behind a linear least squares fit are
     * autocorrelated at lag one.
     * <p>
     * The statistic is
     * {@code d = sum (e_i - e_(i-1))^2 / sum e_i^2}, which is near {@code 2}
     * when successive residuals are unrelated, near {@code 0} when each one
     * follows the last and near {@code 4} when each one reverses it.
     * <p>
     * <b>The alternatives are about the autocorrelation, not about the
     * statistic, and the two run opposite ways.</b>
     * {@link Alternative#GREATER} asks whether the errors are <em>positively</em>
     * autocorrelated, which is a <em>small</em> {@code d} and therefore the
     * lower tail of the null distribution; {@link Alternative#LESS} asks the
     * other question and reads the upper tail. {@link Alternative#TWO_SIDED} is
     * twice the smaller of the two.
     * <p>
     * The p-value is exact rather than a pair of bounds, and that is what the
     * design matrix is for: {@code d} is a ratio of quadratic forms whose null
     * distribution depends on the column space of {@code X}, so there is no
     * table to look it up in. The work is dominated by an eigendecomposition of
     * an {@code n x n} matrix, which is cubic in the number of observations and
     * is repeated on every call. A caller testing many fits against one design
     * should compute {@link DurbinWatson#nullEigenvalues} once and read the
     * tails off {@link DurbinWatson#cdf} and {@link DurbinWatson#barF} instead.
     * <p>
     * <b>The test assumes normal errors and a design that does not contain a
     * lagged regressand.</b> The first is what makes the null distribution the
     * one computed here; the second is the classical way this test misleads,
     * and no computation on {@code X} can detect it.
     *
     * @param residuals
     *            the residuals of the fit, in the order the observations were
     *            taken, at least two of them and all finite
     * @param designMatrix
     *            the design matrix the fit went through, with as many rows as
     *            there are residuals and full column rank
     * @param alternative
     *            which departure to look for, stated about the autocorrelation
     * @return the statistic, its p-value and no degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code residuals} is {@code null}, is shorter than two or
     *             holds a value that is not finite; if every residual is zero;
     *             if {@code designMatrix} is {@code null}, has the wrong number
     *             of rows, has at least as many columns as rows, or does not
     *             have full column rank; or if {@code alternative} is
     *             {@code null}
     */
    public static TestResult durbinWatson(double[] residuals, DMatrix designMatrix,
            Alternative alternative) {
        int n = requireFiniteSample(residuals, "residuals");
        requireAtLeastTwo(n, "residuals");
        if (designMatrix == null) {
            throw new IllegalArgumentException("designMatrix must not be null");
        }
        requireAlternative(alternative);
        if (designMatrix.numRows() != n) {
            throw new IllegalArgumentException("the design matrix has " + designMatrix.numRows()
                    + " rows against " + n + " residuals");
        }

        double statistic = durbinWatsonStatistic(residuals);
        double[] nullEigenvalues = DurbinWatson.nullEigenvalues(designMatrix.getArrayUnsafe(), n,
                designMatrix.numColumns());
        double pValue;
        switch (alternative) {
        case GREATER:
            // positive autocorrelation drives the statistic down
            pValue = DurbinWatson.cdf(nullEigenvalues, statistic);
            break;
        case LESS:
            pValue = DurbinWatson.barF(nullEigenvalues, statistic);
            break;
        default:
            pValue = Math.min(1.0, 2.0 * Math.min(DurbinWatson.cdf(nullEigenvalues, statistic),
                    DurbinWatson.barF(nullEigenvalues, statistic)));
            break;
        }
        return new TestResult("Durbin-Watson", statistic, pValue, alternative, Double.NaN);
    }

    /**
     * Tests the residuals of a fit that {@code math.linalg.OLS} produced.
     * <p>
     * See {@link #durbinWatson(double[], DMatrix, Alternative)} for what the
     * test does and which way the alternatives point.
     * <p>
     * A <b>weighted</b> fit is refused. {@code LSSummary.getResiduals()} returns
     * the raw {@code y - X beta} for a weighted fit as well as for an unweighted
     * one, and those are not the residuals whose null distribution this is --
     * the weighted ones are, against the row-scaled design. Passing the raw
     * pair would give a number that looks like a p-value and is not one.
     *
     * @param fit
     *            an unweighted least squares summary that still holds its
     *            residuals and its design matrix
     * @param alternative
     *            which departure to look for, stated about the autocorrelation
     * @return the statistic, its p-value and no degrees of freedom
     * @throws IllegalArgumentException
     *             if {@code fit} is {@code null}, came from a weighted fit, or
     *             has had {@code clearTemporaries()} called on it; or for any
     *             of the reasons the other overload throws
     */
    public static TestResult durbinWatson(LSSummary fit, Alternative alternative) {
        if (fit == null) {
            throw new IllegalArgumentException("fit must not be null");
        }
        if (fit.getWeights() != null) {
            throw new IllegalArgumentException("this is a weighted fit, whose raw residuals are not the "
                    + "ones this null distribution belongs to");
        }
        DoubleList residuals = fit.getResiduals();
        DMatrix designMatrix = fit.getXMatrix();
        if (residuals == null || designMatrix == null) {
            throw new IllegalArgumentException("the fit no longer holds its residuals or its design "
                    + "matrix; clearTemporaries() releases both");
        }
        double[] e = new double[residuals.size()];
        for (int i = 0; i < e.length; i++) {
            e[i] = residuals.get(i);
        }
        return durbinWatson(e, designMatrix, alternative);
    }

    // ------------------------------------------------------- the assembly --

    /**
     * Wraps a chi-squared statistic in its result. The upper tail comes from
     * {@link ProbabilityFuncs#chiSquareComplemented(double, double)}, which
     * integrates it directly rather than subtracting the lower tail from one
     * and losing the answer to cancellation.
     */
    private static TestResult chiSquaredResult(String name, double statistic, double degreesOfFreedom) {
        double pValue = Math.min(1.0, Math.max(0.0,
                ProbabilityFuncs.chiSquareComplemented(degreesOfFreedom, statistic)));
        return new TestResult(name, statistic, pValue, Alternative.GREATER, degreesOfFreedom);
    }

    /**
     * Builds the result from the four pieces every t-test has: the estimate, the
     * value the null hypothesis puts it at, the standard error of the estimate,
     * and the degrees of freedom.
     */
    private static TTestResult tTest(String name, double estimate, double nullValue, double standardError,
            double degreesOfFreedom, Alternative alternative, double confidenceLevel) {
        double statistic = (estimate - nullValue) / standardError;
        StudentT nullDistribution = new StudentT(degreesOfFreedom);
        double pValue = tailProbability(nullDistribution, statistic, alternative);
        TestResult test = new TestResult(name, statistic, pValue, alternative, degreesOfFreedom);
        return interval(test, estimate, standardError, confidenceLevel, nullDistribution);
    }

    /**
     * The Welch-Satterthwaite degrees of freedom from the two standard errors
     * {@code ex} and {@code ey}.
     * <p>
     * The expression squares them twice, which would overflow from about
     * {@code 1e77} upwards. The result is a ratio in which the common scale
     * cancels exactly, so it is divided out before the powers are formed.
     */
    private static double welchDegreesOfFreedom(double ex, double ey, int nx, int ny) {
        double largest = Math.max(Math.abs(ex), Math.abs(ey));
        double p = ex / largest;
        double q = ey / largest;
        double p2 = p * p;
        double q2 = q * q;
        double sum = p2 + q2;
        return (sum * sum) / (p2 * p2 / (nx - 1.0) + q2 * q2 / (ny - 1.0));
    }

    // ------------------------------------------------------------ the tails --

    /**
     * The p-value of {@code statistic} against a symmetric null distribution.
     * <p>
     * The upper tail is taken as {@code cdf(-statistic)} rather than as
     * {@code 1 - cdf(statistic)}: the two are equal in exact arithmetic, but
     * the subtraction loses every significant digit of a small p-value.
     */
    private static double tailProbability(StudentT nullDistribution, double statistic, Alternative alternative) {
        switch (alternative) {
        case LESS:
            return nullDistribution.cdf(statistic);
        case GREATER:
            return nullDistribution.cdf(-statistic);
        default:
            return Math.min(1.0, 2.0 * nullDistribution.cdf(-Math.abs(statistic)));
        }
    }

    /**
     * The confidence interval that matches the alternative the test was run
     * with, so that the interval excludes the null value exactly when the test
     * rejects.
     */
    private static TTestResult interval(TestResult test, double estimate, double standardError,
            double confidenceLevel, StudentT nullDistribution) {
        double lower;
        double upper;
        switch (test.alternative) {
        case LESS:
            lower = Double.NEGATIVE_INFINITY;
            upper = estimate + nullDistribution.inverseCdf(confidenceLevel) * standardError;
            break;
        case GREATER:
            lower = estimate - nullDistribution.inverseCdf(confidenceLevel) * standardError;
            upper = Double.POSITIVE_INFINITY;
            break;
        default:
            double half = nullDistribution.inverseCdf(1.0 - 0.5 * (1.0 - confidenceLevel)) * standardError;
            lower = estimate - half;
            upper = estimate + half;
            break;
        }
        return new TTestResult(test, estimate, standardError, confidenceLevel, lower, upper);
    }

    // ----------------------------------------------------------- the moments --

    /**
     * The mean and the sample standard deviation -- the {@code n - 1}
     * denominator -- of {@code x}, in that order.
     * <p>
     * A deviation rather than a variance because the square of one leaves the
     * double range far sooner: data at {@code 1e300} has a variance of about
     * {@code 1e596} and data at {@code 1e-320} one of about {@code 1e-640},
     * neither representable, while both deviations and the scale invariant
     * statistic built from them are. This is why
     * {@code math.list.DoubleArrayList} returns a deviation too.
     * <p>
     * The variance underneath is the corrected two-pass: a second pass over the
     * deviations, minus the square of their residual sum, which is what carries
     * the rounding error left behind in the mean. The one-pass
     * {@code E[x^2] - mean^2} is unusable once the data sits far from zero, and
     * Welford is the streaming algorithm and pays for being one; both were
     * measured elsewhere in this library.
     * <p>
     * The rescaling is optimistic, in the shape
     * {@code math.list.DoubleArrayList.norm2} uses: a square leaves the double
     * range long before a deviation does, so the direct route is taken first and
     * the data is scaled to about one -- by a power of two, which is exact --
     * only once the sum of squares has overflowed or dropped into the
     * subnormals. The common case pays nothing for it.
     */
    static double[] meanAndDeviation(double[] x) {
        final int n = x.length;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i];
        }
        if (isFinite(sum)) {
            double mean = sum / n;
            double sumDev = 0.0;
            double sumSqrDev = 0.0;
            for (int i = 0; i < n; i++) {
                double d = x[i] - mean;
                sumDev += d;
                sumSqrDev += d * d;
            }
            if (sumSqrDev >= n * Double.MIN_NORMAL && isFinite(sumSqrDev)) {
                return new double[] { mean, Math.sqrt(corrected(sumSqrDev, sumDev, n)) };
            }
        }

        double maxAbs = 0.0;
        for (int i = 0; i < n; i++) {
            maxAbs = Math.max(maxAbs, Math.abs(x[i]));
        }
        if (maxAbs == 0.0) {
            return new double[] { 0.0, 0.0 };
        }
        int k = -Math.getExponent(maxAbs);
        if (k < -1022) {
            // the factor itself would turn subnormal and stop being exact
            k = -1022;
        }
        final double scale = Math.scalb(1.0, k);
        double scaledSum = 0.0;
        for (int i = 0; i < n; i++) {
            scaledSum += x[i] * scale;
        }
        final double scaledMean = scaledSum / n;
        double sumDev = 0.0;
        double sumSqrDev = 0.0;
        for (int i = 0; i < n; i++) {
            double d = x[i] * scale - scaledMean;
            sumDev += d;
            sumSqrDev += d * d;
        }
        // the mean and the deviation both scale by the factor
        return new double[] { Math.scalb(scaledMean, -k),
                Math.scalb(Math.sqrt(corrected(sumSqrDev, sumDev, n)), -k) };
    }

    /** {@code sum (e_i - e_(i-1))^2 / sum e_i^2}, in one pass. */
    private static double durbinWatsonStatistic(double[] e) {
        double numerator = 0.0;
        double denominator = e[0] * e[0];
        for (int i = 1; i < e.length; i++) {
            double step = e[i] - e[i - 1];
            numerator += step * step;
            denominator += e[i] * e[i];
        }
        if (denominator == 0.0) {
            throw new IllegalArgumentException("every residual is zero, so the statistic is 0/0");
        }
        return numerator / denominator;
    }

    /** The two-pass variance, corrected by the residual sum of the deviations. */
    private static double corrected(double sumSqrDev, double sumDev, int n) {
        double variance = (sumSqrDev - (sumDev * sumDev) / n) / (n - 1);
        // the subtraction carries no sign guarantee, and a negative variance
        // would turn the standard error into NaN
        return variance > 0.0 ? variance : 0.0;
    }

    // -------------------------------------------------------- the guard rail --

    /** Rejects a sample that is absent, empty or not made of numbers. */
    private static int requireFiniteSample(double[] x, String name) {
        if (x == null) {
            throw new IllegalArgumentException(name + " must not be null");
        }
        if (x.length == 0) {
            throw new IllegalArgumentException(name + " must not be empty");
        }
        for (int i = 0; i < x.length; i++) {
            if (!isFinite(x[i])) {
                throw new IllegalArgumentException(name + "[" + i + "] is not a finite number : " + x[i]);
            }
        }
        return x.length;
    }

    private static void requireAtLeastTwo(int n, String name) {
        if (n < 2) {
            throw new IllegalArgumentException(name + " needs at least two observations, got " + n);
        }
    }

    private static void requireAlternative(Alternative alternative) {
        if (alternative == null) {
            throw new IllegalArgumentException("alternative must not be null");
        }
    }

    private static void requireLevel(double confidenceLevel) {
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException("confidenceLevel must lie in (0, 1) : " + confidenceLevel);
        }
    }

    /**
     * How far the expected counts of a goodness of fit test may sum away from
     * the observed total before they are taken to be something else, such as
     * probabilities.
     */
    private static final double TOTALS_TOLERANCE = 1.0e-9;

    /**
     * How far above the probability of the observed table another table may
     * still be counted as no more likely than it. Two tables of equal
     * probability reach that probability by different sums of logarithms, and
     * the largest disagreement measured over the tables in the test is 5e-13
     * relative.
     */
    private static final double LIKELIHOOD_TOLERANCE = 1.0e-9;

    /**
     * How close to the ends of {@code (0, 1)} a transformed observation is
     * allowed to come before {@link #andersonDarling} pulls it back.
     * <p>
     * A distribution function written in {@code double} arithmetic returns
     * exactly {@code 0} or {@code 1} far more often than the mathematics says
     * it should -- 1310 of 80000 draws from a {@code t(1.5)} read through a
     * standard normal did -- and the logarithm of either end is infinite. The
     * floor turns an infinite statistic into a merely enormous one, which is
     * the same decision with a number attached.
     */
    private static final double UNIFORM_FLOOR = MathConsts.BIG_INV / 2.0;

    /**
     * Whether the sorted {@code sample} repeats a value from {@code from}
     * onwards. The merge walk of the two-sample test sees every repetition up
     * to the point where one sample runs out, and this covers the tail of the
     * other one.
     */
    private static boolean hasARepetition(double[] sample, int from) {
        for (int i = Math.max(1, from); i < sample.length; i++) {
            if (sample[i] == sample[i - 1]) {
                return true;
            }
        }
        return false;
    }

    /** What the fitted family is called in the name of the test. */
    private static String fittedName(Lilliefors.Family family) {
        if (family == null) {
            throw new IllegalArgumentException("family must not be null");
        }
        switch (family) {
        case EXPONENTIAL:
            return "exponential";
        case LOGNORMAL:
            return "log-normal";
        default:
            return "normal";
        }
    }

    /** {@code Double.isFinite} of Java 8, spelled out for the release 8 target. */
    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }

    private HypothesisTests() {
        throw new AssertionError();
    }
}
