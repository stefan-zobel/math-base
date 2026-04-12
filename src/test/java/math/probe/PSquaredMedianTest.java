package math.probe;

import org.junit.Before;
import org.junit.Test;

import math.rng.DefaultRng;
import math.rng.PseudoRandom;

import static org.junit.Assert.*;

import java.util.stream.DoubleStream;

/**
 * Unit tests for PSquaredMedian.
 */
public class PSquaredMedianTest {

    private PSquaredMedian p2;
    private static final double EPSILON = 0.05; // 5% error margin for approximation

    @Before
    public void setUp() {
        p2 = new PSquaredMedian();
    }

    /**
     * Tests the median with a small, fixed set of numbers.
     */
    @Test
    public void testSmallFixedDataset() {
        double[] values = {10.0, 2.0, 38.0, 23.0, 1.0, 7.0, 15.0};
        for (double v : values) {
            p2.accept(v);
        }

        // Exact median of {1, 2, 7, 10, 15, 23, 38} is 10.0
        assertEquals(10.0, p2.getMedian(), 1.0);
    }

    /**
     * Tests the algorithm against an Exponential Distribution.
     * The theoretical median of an exponential distribution is ln(2)/lambda.
     */
    @Test
    public void testExponentialDistribution() {
        double lambda = 0.5;
        double theoreticalMedian = Math.log(2) / lambda; // approx 1.386

        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        int samples = 10_000;
        DoubleStream stream = rng.exponential(samples, lambda);
        stream.forEach(d -> p2.accept(d));

        double estimatedMedian = p2.getMedian();

        System.out.println("Exponential Distribution (lambda=" + lambda + ")");
        System.out.println("Theoretical Median: " + theoreticalMedian);
        System.out.println("P-Square Estimate: " + estimatedMedian);

        // Check if estimate is within 5% of theoretical value
        assertEquals(theoreticalMedian, estimatedMedian, theoreticalMedian * EPSILON);
    }

    @Test
    public void testAverageError() {
        final double lambda = 0.5;
        final double theoreticalMedian = Math.log(2) / lambda;
        final int iterations = 100;
        double totalError = 0.0;

        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        int samples = 10_000;

        for (int j = 0; j < iterations; j++) {
            PSquaredMedian p2Test = new PSquaredMedian();
            DoubleStream stream = rng.exponential(samples, lambda);
            stream.forEach(d -> p2Test.accept(d));

            totalError += (p2Test.getMedian() - theoreticalMedian);
        }
        System.out.println("Average Error: " + (totalError / iterations));
        assertTrue("Average Error <= 0.015", (totalError / iterations) <= 0.015);
    }

    /**
     * Ensures the count is tracked correctly.
     */
    @Test
    public void testCountTracking() {
        int iterations = 123;
        for (int i = 0; i < iterations; i++) {
            p2.accept(i);
        }
        assertEquals(iterations, p2.getCount());
    }
}
