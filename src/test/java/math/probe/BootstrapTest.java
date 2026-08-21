package math.probe;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.concurrent.ForkJoinPool;

import org.junit.Test;

/**
 * Unit tests for {@link Bootstrap}. The bounds are not guesses: they were
 * measured over the seeds and sample sizes used here and then given at least
 * threefold headroom.
 */
public final class BootstrapTest {

    private static final SampleStatistic MEAN = new SampleStatistic() {
        @Override
        public double apply(double[] sample) {
            double total = 0.0;
            for (int i = 0; i < sample.length; ++i) {
                total += sample[i];
            }
            return total / sample.length;
        }
    };

    private static final SampleStatistic VARIANCE = new SampleStatistic() {
        @Override
        public double apply(double[] sample) {
            double m = MEAN.apply(sample);
            double total = 0.0;
            for (int i = 0; i < sample.length; ++i) {
                double d = sample[i] - m;
                total += d * d;
            }
            return total / (sample.length - 1);
        }
    };

    /** Deterministic standard normal, Box-Muller on an inline LCG. */
    private static double[] normal(int n, long seed) {
        double[] x = new double[n];
        long s = seed;
        for (int i = 0; i < n; i += 2) {
            s = s * 6364136223846793005L + 1442695040888963407L;
            double u1 = (s >>> 11) * 0x1.0p-53;
            s = s * 6364136223846793005L + 1442695040888963407L;
            double u2 = (s >>> 11) * 0x1.0p-53;
            if (u1 < 1.0e-300) {
                u1 = 1.0e-300;
            }
            double r = Math.sqrt(-2.0 * Math.log(u1));
            double a = 2.0 * Math.PI * u2;
            x[i] = r * Math.cos(a);
            if (i + 1 < n) {
                x[i + 1] = r * Math.sin(a);
            }
        }
        return x;
    }

    private static double stdDev(double[] x) {
        return Math.sqrt(VARIANCE.apply(x));
    }

    @Test
    public void testTheSameSeedReproducesEverythingExactly() {
        // the reason the seeded constructor exists. Tolerance 0.0: the
        // replications are derived from the replication index, so the two runs
        // do the same arithmetic in the same order
        double[] x = normal(300, 42L);
        int[] iterations = { 1000, 5000, 20000, 100000 };
        for (int k = 0; k < iterations.length; ++k) {
            int b = iterations[k];
            Bootstrap first = new Bootstrap(x, MEAN, b, 12345L);
            Bootstrap second = new Bootstrap(x, MEAN, b, 12345L);
            assertEquals("mean, B=" + b, first.getMean(), second.getMean(), 0.0);
            assertEquals("stddev, B=" + b, first.getStdDev(), second.getStdDev(), 0.0);
            double[] p1 = first.getConfidenceInterval(0.95);
            double[] p2 = second.getConfidenceInterval(0.95);
            assertEquals("percentile lower, B=" + b, p1[0], p2[0], 0.0);
            assertEquals("percentile upper, B=" + b, p1[1], p2[1], 0.0);
            double[] b1 = first.getConfidenceIntervalBCa(0.95);
            double[] b2 = second.getConfidenceIntervalBCa(0.95);
            assertEquals("BCa lower, B=" + b, b1[0], b2[0], 0.0);
            assertEquals("BCa upper, B=" + b, b1[1], b2[1], 0.0);
        }
    }

    @Test
    public void testTheResultDoesNotDependOnHowTheStreamIsScheduled() throws Exception {
        // resampling runs on a parallel stream. Executed from inside a
        // ForkJoinPool task the stream uses that pool, so this really is a
        // different number of workers and a different work-stealing pattern.
        // The seed is derived per replication index rather than from state
        // shared between the threads, so the outcome must not move
        double[] x = normal(400, 7L);
        double[] common = new Bootstrap(x, MEAN, 20000, 999L).getConfidenceIntervalBCa(0.95);
        int[] workers = { 1, 2, 3, 7 };
        for (int k = 0; k < workers.length; ++k) {
            final int w = workers[k];
            ForkJoinPool pool = new ForkJoinPool(w);
            try {
                double[] got = pool.submit(new java.util.concurrent.Callable<double[]>() {
                    @Override
                    public double[] call() {
                        return new Bootstrap(x, MEAN, 20000, 999L).getConfidenceIntervalBCa(0.95);
                    }
                }).get();
                assertEquals("BCa lower on " + w + " workers", common[0], got[0], 0.0);
                assertEquals("BCa upper on " + w + " workers", common[1], got[1], 0.0);
            } finally {
                pool.shutdown();
            }
        }
    }

    @Test
    public void testDifferentSeedsDisagreeInTheNoiseAndAgreeOnTheAnswer() {
        double[] x = normal(400, 555L);
        double[] previous = null;
        int differing = 0;
        for (long seed = 1L; seed <= 6L; ++seed) {
            double[] ci = new Bootstrap(x, MEAN, 20000, seed).getConfidenceInterval(0.95);
            assertTrue("interval must bracket the sample mean, seed " + seed,
                    ci[0] < MEAN.apply(x) && MEAN.apply(x) < ci[1]);
            if (previous != null) {
                if (ci[0] != previous[0]) {
                    ++differing;
                }
                // Monte Carlo noise only: measured worst shift was 0.014 of
                // the interval width
                double width = ci[1] - ci[0];
                assertTrue("seed " + seed + " moved the interval too far",
                        Math.abs(ci[0] - previous[0]) / width < 0.05);
            }
            previous = ci;
        }
        assertTrue("different seeds must give different replications", differing >= 4);
    }

    @Test
    public void testThePercentileIntervalMatchesTheClassicalOne() {
        // for the mean of a normal sample the bootstrap has a closed-form
        // counterpart, which is the only external check available here.
        // Measured worst deviation of the half width was 1.2 percent
        for (int n = 100; n <= 800; n *= 2) {
            for (long seed = 1L; seed <= 4L; ++seed) {
                double[] x = normal(n, seed * 7919L);
                double[] ci = new Bootstrap(x, MEAN, 20000, seed).getConfidenceInterval(0.95);
                double classical = 1.96 * stdDev(x) / Math.sqrt(n);
                double bootstrap = 0.5 * (ci[1] - ci[0]);
                assertEquals("half width, n=" + n + ", seed " + seed, 1.0, bootstrap / classical, 0.05);
            }
        }
    }

    @Test
    public void testTheBootstrapMeanTracksTheSampleStatistic() {
        // the bootstrap distribution of the mean is centred on the sample
        // mean; any bias must be far below one standard error. Measured worst
        // was 0.017 standard errors
        for (long seed = 1L; seed <= 8L; ++seed) {
            double[] x = normal(400, seed * 104729L);
            Bootstrap bs = new Bootstrap(x, MEAN, 20000, seed);
            double standardError = stdDev(x) / Math.sqrt(400);
            double bias = Math.abs(bs.getMean() - MEAN.apply(x));
            assertTrue("bias " + bias + " too large at seed " + seed, bias / standardError < 0.05);
        }
    }

    @Test
    public void testBCaLeavesASymmetricStatisticAloneAndCorrectsASkewedOne() {
        // the mean of a normal sample is symmetric and unbiased, so BCa has
        // almost nothing to correct: measured worst shift was 0.7 percent of
        // the interval width
        for (long seed = 1L; seed <= 8L; ++seed) {
            double[] x = normal(400, seed * 15485863L);
            Bootstrap bs = new Bootstrap(x, MEAN, 20000, seed);
            double[] percentile = bs.getConfidenceInterval(0.95);
            double[] bca = bs.getConfidenceIntervalBCa(0.95);
            double width = percentile[1] - percentile[0];
            assertTrue("BCa moved a symmetric statistic, seed " + seed,
                    Math.abs(percentile[0] - bca[0]) / width < 0.05
                            && Math.abs(percentile[1] - bca[1]) / width < 0.05);
        }
        // the sample variance is right-skewed, and there BCa has to move both
        // ends upwards; that is the whole point of the correction
        for (long seed = 1L; seed <= 4L; ++seed) {
            double[] x = normal(200, seed * 32452843L);
            Bootstrap bs = new Bootstrap(x, VARIANCE, 20000, seed);
            double[] percentile = bs.getConfidenceInterval(0.95);
            double[] bca = bs.getConfidenceIntervalBCa(0.95);
            assertTrue("BCa lower must exceed the percentile lower, seed " + seed, bca[0] > percentile[0]);
            assertTrue("BCa upper must exceed the percentile upper, seed " + seed, bca[1] > percentile[1]);
        }
    }

    @Test
    public void testTheSeedlessConstructorStillProducesAUsableInterval() {
        // it draws its own seed, so nothing about its exact values can be
        // asserted; what must hold is that it answers at all
        double[] x = normal(300, 2024L);
        Bootstrap bs = new Bootstrap(x, MEAN, 5000);
        double[] percentile = bs.getConfidenceInterval(0.95);
        double[] bca = bs.getConfidenceIntervalBCa(0.95);
        double m = MEAN.apply(x);
        assertTrue("percentile must bracket the sample mean", percentile[0] < m && m < percentile[1]);
        assertTrue("BCa must bracket the sample mean", bca[0] < m && m < bca[1]);
        assertFalse("mean must be finite", Double.isNaN(bs.getMean()) || Double.isInfinite(bs.getMean()));
        assertFalse("stddev must be finite", Double.isNaN(bs.getStdDev()) || Double.isInfinite(bs.getStdDev()));
    }

    @Test
    public void testTheConfidenceLevelIsChecked() {
        double[] x = normal(64, 3L);
        Bootstrap bs = new Bootstrap(x, MEAN, 1000, 17L);
        double[] levels = { 0.0, 1.0, -0.5, 1.5 };
        for (int i = 0; i < levels.length; ++i) {
            try {
                bs.getConfidenceInterval(levels[i]);
                org.junit.Assert.fail("percentile accepted " + levels[i]);
            } catch (IllegalArgumentException expected) {
                // that is the contract
            }
            try {
                bs.getConfidenceIntervalBCa(levels[i]);
                org.junit.Assert.fail("BCa accepted " + levels[i]);
            } catch (IllegalArgumentException expected) {
                // that is the contract
            }
        }
    }
}
