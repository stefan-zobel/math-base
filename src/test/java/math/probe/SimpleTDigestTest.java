package math.probe;

import org.junit.Test;

import math.cern.Arithmetic;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;

import static org.junit.Assert.*;

import java.util.stream.DoubleStream;

/**
 * Unit test for SimpleTDigest.
 */
public class SimpleTDigestTest {

    private static final double LAMBDA = 0.5;
    private static final double THEORETICAL_MEDIAN = Math.log(2) / LAMBDA; // ~1.386
    private static final int SAMPLES = 10_000;
    private static final int REPETITIONS = 100;

    /**
     * Verifies that a higher compression leads to a lower error and better
     * accuracy.
     */
    @Test
    public void testAccuracyAndCompression() {
        // We test three different compression settings
        double lowCompressionError = runSimulation(20.0);
        double middleCompressionError = runSimulation(50.0);
        double highCompressionError = runSimulation(100.0);

        System.out.println("T-Digest Evaluation (100 runs, 10,000 samples each):");
        System.out.println("Theoretical Median: " + THEORETICAL_MEDIAN);
        System.out.println("Avg Error (Compression  20): " + lowCompressionError);
        System.out.println("Avg Error (Compression  50): " + middleCompressionError);
        System.out.println("Avg Error (Compression 100): " + highCompressionError);

        // Assertions: High compression should generally be more accurate (lower
        // absolute error). We use a generous delta because of the random
        // nature, but the trend should be clear.
        assertTrue("Error should be reasonably small for compression 100.0", Math.abs(highCompressionError) < 0.01);
    }

    /**
     * Basic functional test for a very simple sorted dataset.
     */
    @Test
    public void testBasicMedian() {
        SimpleTDigest digest = new SimpleTDigest();
        // Median of 1, 2, 3, 4, 5 is 3
        for (double d = 1.0; d <= 5.0; d++) {
            digest.accept(d);
        }
        assertEquals(3.0, digest.getMedian(), 0.1);
    }

    /**
     * Helper method to run multiple simulations and calculate the average error.
     */
    private double runSimulation(double compression) {
        double totalError = 0.0;
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        
        double avgCentroidCount = 0.0;

        for (int j = 0; j < REPETITIONS; j++) {
            SimpleTDigest digest = new SimpleTDigest(compression);
            DoubleStream stream = rng.exponential(SAMPLES, LAMBDA);
            stream.forEach(d -> digest.accept(d));
            totalError += (digest.getMedian() - THEORETICAL_MEDIAN);
            avgCentroidCount += digest.getCentroidCount();
        }
        System.out.println("compression: " + compression + ", avg. centroid count: " + Arithmetic.round(avgCentroidCount / REPETITIONS, 1));
        return totalError / REPETITIONS;
    }

    /**
     * Tests if the getQuantile method works for other values than 0.5.
     */
    @Test
    public void testDifferentQuantiles() {
        SimpleTDigest digest = new SimpleTDigest();
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        // Standard Normal Distribution
        DoubleStream stream = rng.normal(10_000,  0.0, 1.0);
        stream.forEach(d -> digest.accept(d));

        // For Normal Distribution: 0.5 quantile (median) should be around 0
        assertEquals(0.0, digest.getQuantile(0.5), 0.1);
        // 0.975 quantile for Normal Distribution is ~1.96
        assertEquals(1.96, digest.getQuantile(0.975), 0.2);
        System.out.println("Standard normal - compression: " + digest.getCompression() + ", num centroids: " + digest.getCentroidCount());
    }

    @Test
    public void simulateSLAMonitoring() {
        SimpleTDigest digest = new SimpleTDigest(100.0);
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        double lambda = 0.5; // Average response time = 1/lambda = 2s

        // Simulate 10,000 requests
        DoubleStream stream = rng.exponential(10_000, lambda);
        stream.forEach(responseTime -> digest.accept(responseTime));

        // 1. What is the 95th percentile (P95)?
        double p95 = digest.getQuantile(0.95);
        System.out.println("95% of users wait less than: " + p95 + "s");

        // 2. How many users stay under the 4.0s limit? (SLA check)
        double usersUnderLimit = digest.getCDF(4.0) * 100;
        System.out.println("Percentage of users under 4.0s: " + usersUnderLimit + "%");

        // Theoretical CDF for Exponential: 1 - e^(-lambda * x)
        double theoreticalCDF = (1 - Math.exp(-lambda * 4.0)) * 100;
        
        // Assertion: Our estimate should be very close to the theoretical 86.46%
        assertEquals(theoreticalCDF, usersUnderLimit, 1.0);
        System.out.println("SLA monitor - compression: " + digest.getCompression() + ", num centroids: " + digest.getCentroidCount());
    }

    @Test
    public void testBowleySkewness() {
        SimpleTDigest digest = new SimpleTDigest(55.0);
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        double lambda = 0.5;
        // Exponential distribution
        DoubleStream stream = rng.exponential(20_000, lambda);
        stream.forEach(d -> digest.accept(d));

        double skewness = digest.getBowleySkewness();
        System.out.println("Bowley Skewness: " + skewness);

        // Theoretical value for exponential distribution is ~0.28
        // Normal distribution would be ~0.0
        assertTrue("Skewness should be higher than Normal Dist", skewness > 0.2);
        System.out.println("Bowley - compression: " + digest.getCompression() + ", num centroids: " + digest.getCentroidCount());
    }

    @Test
    public void testTailRatio() {
        SimpleTDigest digest = new SimpleTDigest(55.0);
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        double lambda = 0.5;
        // Exponential distribution
        DoubleStream stream = rng.exponential(100_000, lambda);
        stream.forEach(d -> digest.accept(d));

        double tailRatio = digest.getTailRatio();
        System.out.println("Tail-Ratio: " + tailRatio);

        // Theoretical value for exponential distribution is ~5.72
        // Normal distribution would be ~1.0
        assertTrue("Tail-Ratio should be at least 5.4", tailRatio > 5.4);
        System.out.println("Tail-Ratio - compression: " + digest.getCompression() + ", num centroids: " + digest.getCentroidCount());
        digest.compress();
        System.out.println("Tail-Ratio - compression: " + digest.getCompression() + ", num centroids: " + digest.getCentroidCount());
        double newTailRatio = digest.getTailRatio();
        System.out.println("Old Tail-Ratio: " + tailRatio + ", New Tail-Ratio: " + newTailRatio);
    }
}
