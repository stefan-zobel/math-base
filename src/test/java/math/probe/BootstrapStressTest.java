package math.probe;

import org.junit.Assert;
import org.junit.Test;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;

public class BootstrapStressTest {

    @Test
    public void stressTestParallelBCa() {
        // Setup: a dataset where we know the theoretical values
        // We use a normal distribution (MU=0, SIGMA=1)
        int n = 2_000;
        int iterations = 100_000; // High load for the parallel streams

        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        double[] data = rng.normal(n, 0.0, 1.0).toArray();

        // We test variance as the statistic (slightly more expensive than the mean)
        SampleStatistic varianceSource = sample -> {
            double avg = 0;
            for (double d : sample) avg += d;
            avg /= sample.length;
            double sumSq = 0;
            for (double d : sample) {
                double diff = d - avg;
                sumSq += diff * diff;
            }
            return sumSq / (sample.length - 1);
        };

        // Start time measurement
        long start = System.currentTimeMillis();

        // Run bootstrap
        Bootstrap bs = new Bootstrap(data, varianceSource, iterations);

        // Calculate BCa confidence interval (also triggers the parallel Jackknife)
        double[] ci = bs.getConfidenceIntervalBCa(0.95);

        long end = System.currentTimeMillis();
        System.out.println("Stress Test finished in " + (end - start) + " ms");
        System.out.println("95% BCa CI for Variance: [" + ci[0] + ", " + ci[1] + "]");

        // Validation
        // A) Basic checks
        Assert.assertEquals("CI must have 2 values", 2, ci.length);
        Assert.assertTrue("Lower bound must be < upper bound", ci[0] < ci[1]);

        // B) Numerical plausibility
        // The variance of a standard normal distribution is 1.0 
        // The CI should contain 1.0
        Assert.assertTrue("CI should contain the true variance 1.0", ci[0] < 1.0 && ci[1] > 1.0);

        // C) Check for NaN (sign of threading issues or division by zero)
        Assert.assertFalse("Lower bound is NaN", Double.isNaN(ci[0]));
        Assert.assertFalse("Upper bound is NaN", Double.isNaN(ci[1]));

        // D) Consistency check: The mean of the bootstrap distribution
        // should be close to the variance of the original sample
        double originalVar = varianceSource.apply(data);
        Assert.assertEquals(originalVar, bs.getMean(), 0.01);
        System.out.println("stressTestParallelBCa:");
        System.out.println(bs.summary(0.95));
    }

    @Test
    public void testExponentialSkewnessBCa() {
        int n = 500;
        int iterations = 50_000;

        // Setup: Exponential distribution (theoretical skewness is always = 2.0, independent of lambda)
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        // Exponential(lambda=1.0)
        double[] data = rng.exponential(n, 1.0).toArray();

        // Statistic: Skewness
        SampleStatistic skewness = sample -> {
            double m = 0, m2 = 0, m3 = 0;
            int len = sample.length;
            for (double d : sample) m += d;
            m /= len;
            for (double d : sample) {
                double diff = d - m;
                m2 += diff * diff;
                m3 += diff * diff * diff;
            }
            // Moment-based skewness
            return (m3 / len) / Math.pow(m2 / len, 1.5);
        };

        // Execution
        Bootstrap bs = new Bootstrap(data, skewness, iterations);
        System.out.println("Exponential Skewness Test:");
        System.out.println(bs.summary(0.95));

        // Validation
        double[] ciBCa = bs.getConfidenceIntervalBCa(0.95);
        double[] ciPercentile = bs.getConfidenceInterval(0.95);

        // BCa should be significantly different from Percentile for skewed distributions
        double diffLower = Math.abs(ciBCa[0] - ciPercentile[0]);
        double diffUpper = Math.abs(ciBCa[1] - ciPercentile[1]);

        System.out.format("Difference BCa vs Percentile: Lower=%.4f, Upper=%.4f\n", diffLower, diffUpper);

        Assert.assertTrue("BCa should correct for skewness", (diffLower + diffUpper) > 0.0);
        Assert.assertTrue("CI should be around 2.0", ciBCa[0] < 2.0 && ciBCa[1] > 1.0);
    }

    @Test
    public void testHighContentionSmallSample() {
        // Tests whether the system remains stable even with extremely small
        // samples (many Jackknife threads) and a high number of iterations
        double[] tinyData = {1.0, 2.0, 3.0, 4.0, 5.0};
        Bootstrap bs = new Bootstrap(tinyData, sample -> sample[0], 50_000);

        double[] ci = bs.getConfidenceIntervalBCa(0.99);
        Assert.assertNotNull(ci);
        System.out.println("testHighContentionSmallSample:");
        System.out.println(bs.summary(0.95));
    }
}
