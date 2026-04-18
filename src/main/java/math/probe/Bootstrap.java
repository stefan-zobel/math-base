package math.probe;

import java.util.Arrays;
import java.util.stream.IntStream;
import math.cern.Arithmetic;
import math.cern.ProbabilityFuncs;
import math.rng.SplittablePseudoRandom;
import math.rng.XorShiftRot256StarStar;

/**
 * A decent-performance, parallelized implementation of the Bootstrap resampling
 * method.
 *
 * This class provides mechanisms to estimate the sampling distribution of
 * nearly any statistic by performing random resampling with replacement (Monte
 * Carlo bootstrap). It supports both the standard percentile method and the
 * more accurate <b>Bias-Corrected and Accelerated (BCa)</b> interval
 * calculation.
 *
 * <p>
 * Features:
 * <ul>
 * <li><b>Parallel Execution:</b> Resampling and Jackknife calculations are
 * parallelized using Java Streams to maximize CPU utilization.</li>
 * <li><b>Thread-Safe RNG:</b> Utilizes {@link math.rng.SplittablePseudoRandom}
 * to provide independent, non-blocking random streams for each thread.</li>
 * <li><b>BCa Support:</b> Adjusts for bias and skewness in the bootstrap
 * distribution, providing higher accuracy for non-symmetrical statistics.</li>
 * </ul>
 *
 * <p>
 * The complexity is approximately O(B * N) for resampling and O(N^2) for the
 * Jackknife acceleration phase, where B is the number of iterations and N is
 * the sample size.
 *
 * @since 1.4.1
 */
public final class Bootstrap {

    private final double[] originalSample;
    private final SampleStatistic statistic;
    private final double observedStatistic;
    private final double[] bootResults;
    private final SplittablePseudoRandom rng;

    public Bootstrap(double[] sample, SampleStatistic statistic, int iterations) {
        this.originalSample = sample.clone();
        this.statistic = statistic;
        this.observedStatistic = statistic.apply(originalSample);
        this.bootResults = new double[iterations];
        this.rng = XorShiftRot256StarStar.getDefault();

        runResampling(iterations);
    }

    private void runResampling(int iterations) {
        int n = originalSample.length;

        // We use Parallel Streams to utilize the hardware cores
        IntStream.range(0, iterations).parallel().forEach(b -> {
            // Each iteration needs its own RNG state, to avoid
            // thread interference (SplittablePseudoRandom is ideal)
            SplittablePseudoRandom localRng = rng.split();

            double[] resample = new double[n];
            for (int i = 0; i < n; i++) {
                resample[i] = originalSample[localRng.nextInt(n)];
            }

            // Since we write into an array at a fixed index, this
            // is thread-safe (no structural modification of the array)
            bootResults[b] = statistic.apply(resample);
        });

        // After the parallel phase: sort once for the quantiles
        Arrays.sort(bootResults);
    }

    public double getMean() {
        double sum = 0.0;
        for (double val : bootResults) {
            sum += val;
        }
        return sum / bootResults.length;
    }

    public double getStdDev() {
        double avg = getMean();
        double sumSq = 0.0;
        for (double val : bootResults) {
            double diff = val - avg;
            sumSq += diff * diff;
        }
        return Math.sqrt(sumSq / (bootResults.length - 1));
    }

    /**
     * @param confidenceLevel
     *            e.g., 0.95 for 95% CI
     * @return [lowerBound, upperBound]
     */
    public double[] getConfidenceInterval(double confidenceLevel) {
        if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0) {
            throw new IllegalArgumentException("confidenceLevel must be in (0,1)");
        }
        double alphaHalf = (1.0 - confidenceLevel) / 2.0;
        return new double[] { getQuantile(alphaHalf), getQuantile(1.0 - alphaHalf) };
    }

    /**
     * Computes the Bias-Corrected and Accelerated (BCa) confidence interval.
     * 
     * The BCa interval is a robust bootstrap method that adjusts for both:
     * <ul>
     * <li><b>Bias:</b> Discrepancy between the median of the bootstrap
     * distribution and the observed sample statistic.</li>
     * <li><b>Acceleration:</b> Changes in the standard error of the statistic
     * with respect to its true value (skewness).</li>
     * </ul>
     * 
     * This method is generally more accurate than the simple percentile method,
     * especially when the bootstrap distribution is skewed or biased. It uses a
     * Jackknife-based acceleration factor and a standard normal distribution
     * transformation to adjust the quantile points.
     * 
     * @param confidenceLevel
     *            the desired confidence level (e.g., 0.95 for 95%)
     * @return an array containing [lower_bound, upper_bound]
     * @throws IllegalArgumentException
     *             if confidenceLevel is not between 0 and 1
     */
    public double[] getConfidenceIntervalBCa(double confidenceLevel) {
        if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0) {
            throw new IllegalArgumentException("confidenceLevel must be in (0,1)");
        }
        int B = bootResults.length;
        double alpha = 1.0 - confidenceLevel;
        double alpha1 = alpha / 2.0;
        double alpha2 = 1.0 - alpha1;

        // Bias correction (z0) with epsilon protection
        long countLower = 0;
        for (double val : bootResults) {
            if (val < observedStatistic) {
                countLower++;
            }
        }

        // Stabilization for z0: proportion must not be exactly 0 or 1
        double eps = 0.5 / B;
        double proportionLower = (double) countLower / B;
        proportionLower = Math.max(eps, Math.min(1.0 - eps, proportionLower));
        double z0 = ProbabilityFuncs.normalInverse(proportionLower);

        // Acceleration (a)
        double a = calculateAcceleration();

        // Corrected quantile points
        double zAlpha1 = ProbabilityFuncs.normalInverse(alpha1);
        double zAlpha2 = ProbabilityFuncs.normalInverse(alpha2);

        // Helper variables for the BCa formula (avoiding division by zero)
        double den1 = 1.0 - a * (z0 + zAlpha1);
        double den2 = 1.0 - a * (z0 + zAlpha2);

        // Protection against division by 0 (very rare, but possible)
        if (Math.abs(den1) < 1e-10) {
            den1 = 1e-10;
        }
        if (Math.abs(den2) < 1e-10) {
            den2 = 1e-10;
        }

        double pLower = ProbabilityFuncs.normal(z0 + (z0 + zAlpha1) / den1);
        double pUpper = ProbabilityFuncs.normal(z0 + (z0 + zAlpha2) / den2);

        // Numerical stabilization (Clamping) for getQuantile
        pLower = Math.max(eps, Math.min(1.0 - eps, pLower));
        pUpper = Math.max(eps, Math.min(1.0 - eps, pUpper));

        return new double[] { getQuantile(pLower), getQuantile(pUpper) };
    }

    public double getQuantile(double p) {
        if (p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0,1]");
        }

        double index = p * (bootResults.length - 1);
        int i0 = (int) Math.floor(index);
        int i1 = (int) Math.ceil(index);

        if (i0 == i1) {
            return bootResults[i0];
        }

        return Arithmetic.interpolateY((double) i0, (double) i1, bootResults[i0], bootResults[i1], index);
    }

    private double calculateAcceleration() { // for BCa
        int n = originalSample.length;
        double[] jackknifeValues = new double[n];

        // Parallelization over the sample indices
        IntStream.range(0, n).parallel().forEach(i -> {
            double[] subset = new double[n - 1];
            // leave i out
            System.arraycopy(originalSample, 0, subset, 0, i);
            System.arraycopy(originalSample, i + 1, subset, i, n - i - 1);

            // Calculation of the "Leave-one-out" statistic
            jackknifeValues[i] = statistic.apply(subset);
        });

        // Calculate the mean of the jackknife values
        double sumJack = 0.0;
        for (double val : jackknifeValues) {
            sumJack += val;
        }
        double meanJack = sumJack / n;

        // Calculate skewness
        double num = 0;
        double den = 0;
        for (double val : jackknifeValues) {
            double diff = meanJack - val;
            num += (diff * diff * diff);
            den += (diff * diff);
        }

        // If the denominator is near zero, acceleration is zero
        if (Math.abs(den) < 1e-18) {
            return 0.0;
        }

        // Equation (11.40) in Efron / Hastie (2016):
        // "Computer Age Statistical Inference", page 194
        return num / (6.0 * Math.pow(den, 1.5));
    }

    /**
     * Returns a formatted summary of the bootstrap results, including mean, 
     * standard deviation, and a comparison of confidence intervals.
     * 
     * @param confidenceLevel the desired confidence level (e.g., 0.95)
     * @return a string containing the statistical summary
     */
    public String summary(double confidenceLevel) {
        double[] ciPercentile = getConfidenceInterval(confidenceLevel);
        double[] ciBCa = getConfidenceIntervalBCa(confidenceLevel);
        double stdErr = getStdDev();
        double mean = getMean();

        StringBuilder sb = new StringBuilder();
        sb.append("Bootstrap Summary (").append(bootResults.length).append(" iterations)\n");
        sb.append("--------------------------------------------------\n");
        sb.append(String.format("Observed Statistic : %.8f\n", observedStatistic));
        sb.append(String.format("Bootstrap Mean     : %.8f\n", mean));
        sb.append(String.format("Bootstrap StdError : %.8f\n", stdErr));
        sb.append(String.format("Bias               : %.8f\n", mean - observedStatistic));
        sb.append("--------------------------------------------------\n");
        sb.append(String.format("%d%% Confidence Intervals:\n", (int) (confidenceLevel * 100)));
        sb.append(String.format("  Percentile       : [%.8f, %.8f]\n", ciPercentile[0], ciPercentile[1]));
        sb.append(String.format("  BCa (corrected)  : [%.8f, %.8f]\n", ciBCa[0], ciBCa[1]));

        // --- Warning for the case of instability ---
        if (Math.abs(ciBCa[0] - ciBCa[1]) < 1e-12) {
            sb.append("\nWARNING: BCa interval has collapsed (lower == upper).\n");
            sb.append("The sample size (N=").append(originalSample.length).append(") might be too small or\n");
            sb.append("the data is non-smooth or (less likely) the statistic is extremely biased for this data.\n");
        }
        sb.append("--------------------------------------------------\n");

        return sb.toString();
    }

    @Override
    public String toString() {
        return String.format("Bootstrap(n=%d, iterations=%d, observed=%.4f)", 
                             originalSample.length, bootResults.length, observedStatistic);
    }
}
