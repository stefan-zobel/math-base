package math.gemm;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.AbstractExecutorService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.Assume;
import org.junit.Test;

/**
 * Correctness and threshold tests for the parallel variants of
 * {@link SgemmMRxNR} and {@link DgemmMRxNR}.
 * <p>
 * Both BLIS-style kernels use {@code NC = 3840} as the column-panel width and
 * {@code PARALLEL_WORK_THRESHOLD = 45_000_000}. Parallelization therefore
 * requires:
 * <ul>
 *   <li>{@code taskCountUncapped(colsB, 3840) >= 3}, i.e. {@code colsB >= 7681}
 *       (three full NC-blocks), <em>and</em></li>
 *   <li>{@code rowsA * colsB * colsA >= 45_000_000}.</li>
 * </ul>
 * Unlike the baseline Gemm kernels, both BLIS kernels handle only the
 * no-transpose, column-major case ({@code incRow = 1}, {@code incCol = lda}).
 * There are therefore no transpose-combination loops in these tests.
 */
public class GemmMRxNRParallelTest {

    // NC used internally by both SgemmMRxNR and DgemmMRxNR
    private static final int NC = 3840;

    // -----------------------------------------------------------------------
    // Basic correctness: parallel result must equal sequential result exactly
    // -----------------------------------------------------------------------

    @Test
    public void sgemmParallelMatchesSequential() throws Exception {
        // n=7681 -> 3 NC-blocks; work = 128 * 7681 * 128 = 125,861,888 >> 45M
        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            verifySgemmCase(128, 7681, 128, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void dgemmParallelMatchesSequential() throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            verifyDgemmCase(128, 7681, 128, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
    }

    // -----------------------------------------------------------------------
    // Verify that tasks are actually submitted when all conditions are met
    // -----------------------------------------------------------------------

    @Test
    public void sgemmActualParallelSubmitsTasks() throws Exception {
        // taskCountUncapped(7681, 3840) = min(availProcs, 3); need >= 3 CPUs
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCase(128, 7681, 128, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmActualParallelSubmitsTasks() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyDgemmCase(128, 7681, 128, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    // -----------------------------------------------------------------------
    // taskCount < 3: n = 7680 gives exactly 2 NC-blocks -> always sequential
    // -----------------------------------------------------------------------

    @Test
    public void sgemmExactNcBoundaryUsesSingleTaskFallback() throws Exception {
        // ceil(7680 / 3840) = 2 -> taskCount <= 2 < 3 -> sequential regardless of work
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCase(128, 7680, 128, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmExactNcBoundaryUsesSingleTaskFallback() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyDgemmCase(128, 7680, 128, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void sgemmNcPlusOneSubmitsParallelTasks() throws Exception {
        // ceil(7681 / 3840) = 3 blocks; work = 128 * 7681 * 128 >> 45M
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCase(128, 7681, 128, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmNcPlusOneSubmitsParallelTasks() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyDgemmCase(128, 7681, 128, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    // -----------------------------------------------------------------------
    // Work-threshold boundary (n = 7681 gives 3 NC-blocks).
    // The k values that straddle the threshold are computed at runtime from
    // getParallelWorkThreshold() so that the same test method works correctly
    // under both the Java-8 scalar kernel (threshold = 45M) and the Java-25
    // Vector API kernel (threshold = 25M) in the multi-release jar.
    //
    // With m = 64, n = 7681:
    //   kLo = floor((threshold - 1) / (m * n))  ->  m * n * kLo  < threshold
    //   kHi = kLo + 1                            ->  m * n * kHi  > threshold
    // -----------------------------------------------------------------------

    @Test
    public void sgemmWorkThresholdLowerBoundaryRemainsSequential() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        final int m = 64, n = 7681;
        final int kLo = (int) ((SgemmMRxNR.getParallelWorkThreshold() - 1) / ((long) m * n));
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCase(m, n, kLo, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmWorkThresholdLowerBoundaryRemainsSequential() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        final int m = 64, n = 7681;
        final int kLo = (int) ((DgemmMRxNR.getParallelWorkThreshold() - 1) / ((long) m * n));
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyDgemmCase(m, n, kLo, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void sgemmWorkThresholdUpperBoundaryGoesParallel() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        final int m = 64, n = 7681;
        final int kHi = (int) ((SgemmMRxNR.getParallelWorkThreshold() - 1) / ((long) m * n)) + 1;
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCase(m, n, kHi, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmWorkThresholdUpperBoundaryGoesParallel() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        final int m = 64, n = 7681;
        final int kHi = (int) ((DgemmMRxNR.getParallelWorkThreshold() - 1) / ((long) m * n)) + 1;
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyDgemmCase(m, n, kHi, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    // -----------------------------------------------------------------------
    // Null executor falls back to sequential
    // -----------------------------------------------------------------------

    @Test
    public void sgemmNullExecutorMatchesSequential() throws Exception {
        verifySgemmCase(128, 7681, 64, 0.75f, -0.25f, null);
    }

    @Test
    public void dgemmNullExecutorMatchesSequential() throws Exception {
        verifyDgemmCase(128, 7681, 64, -1.25, 0.5, null);
    }

    // -----------------------------------------------------------------------
    // Shutdown executor falls back to sequential
    // -----------------------------------------------------------------------

    @Test
    public void sgemmShutdownExecutorFallsBackToSequential() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        executor.shutdown();
        try {
            verifySgemmCase(128, 7681, 128, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmShutdownExecutorFallsBackToSequential() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        executor.shutdown();
        try {
            verifyDgemmCase(128, 7681, 128, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    // -----------------------------------------------------------------------
    // alpha == 0 early-exit: no tasks must be submitted
    // -----------------------------------------------------------------------

    @Test
    public void sgemmAlphaZeroBetaZeroClearsC() throws Exception {
        verifySgemmAlphaZeroCase(0.0f);
    }

    @Test
    public void dgemmAlphaZeroBetaZeroClearsC() throws Exception {
        verifyDgemmAlphaZeroCase(0.0);
    }

    @Test
    public void sgemmAlphaZeroScalesC() throws Exception {
        verifySgemmAlphaZeroCase(-0.5f);
    }

    @Test
    public void dgemmAlphaZeroScalesC() throws Exception {
        verifyDgemmAlphaZeroCase(-0.5);
    }

    // -----------------------------------------------------------------------
    // Special alpha / beta combinations
    // -----------------------------------------------------------------------

    @Test
    public void sgemmAlphaOneBetaZeroMatchesSequential() throws Exception {
        verifySgemmCase(33, 7681, 57, 1.0f, 0.0f, null);
    }

    @Test
    public void sgemmAlphaOneBetaOneMatchesSequential() throws Exception {
        verifySgemmCase(33, 7681, 57, 1.0f, 1.0f, null);
    }

    @Test
    public void dgemmAlphaOneBetaZeroMatchesSequential() throws Exception {
        verifyDgemmCase(33, 7681, 57, 1.0, 0.0, null);
    }

    @Test
    public void dgemmAlphaOneBetaOneMatchesSequential() throws Exception {
        verifyDgemmCase(33, 7681, 57, 1.0, 1.0, null);
    }

    // -----------------------------------------------------------------------
    // Input arrays and padding outside the C matrix must not be modified
    // -----------------------------------------------------------------------

    @Test
    public void sgemmDoesNotModifyInputsOrPadding() throws Exception {
        verifySgemmPreservesInputAndPadding(13, 7681, 11, 0.75f, -0.25f);
    }

    @Test
    public void dgemmDoesNotModifyInputsOrPadding() throws Exception {
        verifyDgemmPreservesInputAndPadding(13, 7681, 11, -1.25, 0.5);
    }

    // -----------------------------------------------------------------------
    // Repeated calls with a shared executor must produce identical results
    // -----------------------------------------------------------------------

    @Test
    public void sgemmRepeatedCallsWithSharedExecutorRemainDeterministic() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifyRepeatedSgemmCase(43, 7681, 256, 0.75f, -0.25f, executor);
            verifyRepeatedSgemmCase(43, 7681, 256, 1.0f, 1.0f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmRepeatedCallsWithSharedExecutorRemainDeterministic() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            // 128 * 7681 * 64 = 62,861,312 >> 45M threshold
            verifyRepeatedDgemmCase(128, 7681, 64, -1.25, 0.5, executor);
            verifyRepeatedDgemmCase(128, 7681, 64, 1.0, 1.0, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    // -----------------------------------------------------------------------
    // Randomized small cases against naive reference
    // -----------------------------------------------------------------------

    @Test
    public void sgemmRandomizedSmallCasesMatchNaiveReference() throws Exception {
        Random random = new Random(0x5EED5EEDL);
        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            for (int caseIndex = 0; caseIndex < 20; caseIndex++) {
                int m = random.nextInt(9);
                int n = random.nextInt(9);
                int k = random.nextInt(9);
                float alpha = pickFloatScalar(random);
                float beta  = pickFloatScalar(random);
                verifySgemmCaseAgainstReference(m, n, k, alpha, beta,
                        (caseIndex & 1) == 0 ? null : executor, 1.0e-4f);
            }
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void dgemmRandomizedSmallCasesMatchNaiveReference() throws Exception {
        Random random = new Random(0xD6EED6EEL);
        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            for (int caseIndex = 0; caseIndex < 20; caseIndex++) {
                int m = random.nextInt(9);
                int n = random.nextInt(9);
                int k = random.nextInt(9);
                double alpha = pickDoubleScalar(random);
                double beta  = pickDoubleScalar(random);
                verifyDgemmCaseAgainstReference(m, n, k, alpha, beta,
                        (caseIndex & 1) == 0 ? null : executor, 1.0e-10);
            }
        } finally {
            shutdown(executor);
        }
    }

    // -----------------------------------------------------------------------
    // Degenerate dimensions (m=0, n=0, or k=0) against naive reference
    // -----------------------------------------------------------------------

    @Test
    public void sgemmDegenerateDimensionsMatchNaiveReference() {
        verifySgemmCaseAgainstReference(0, 7, 5, 0.75f, -0.25f, null, 0.0f);
        verifySgemmCaseAgainstReference(5, 0, 7, 0.75f, -0.25f, null, 0.0f);
        verifySgemmCaseAgainstReference(5, 7, 0, 0.75f, -0.25f, null, 0.0f);
    }

    @Test
    public void dgemmDegenerateDimensionsMatchNaiveReference() {
        verifyDgemmCaseAgainstReference(0, 7, 5, -1.25, 0.5, null, 0.0);
        verifyDgemmCaseAgainstReference(5, 0, 7, -1.25, 0.5, null, 0.0);
        verifyDgemmCaseAgainstReference(5, 7, 0, -1.25, 0.5, null, 0.0);
    }

    // -----------------------------------------------------------------------
    // High-work cases against naive reference (parallel path exercised)
    // -----------------------------------------------------------------------

    @Test
    public void sgemmHighWorkCasesMatchNaiveReference() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            verifySgemmCaseAgainstReference(43, 7681, 256, 0.75f, -0.25f, executor, 1.0e-3f);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmHighWorkCasesMatchNaiveReference() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() >= 3);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            // 64 * 7681 * 92 = 45,225,728 > 45M threshold; naive reference stays manageable
            verifyDgemmCaseAgainstReference(64, 7681, 92, -1.25, 0.5, executor, 1.0e-10);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    // -----------------------------------------------------------------------
    // taskCountUncapped structural tests
    // -----------------------------------------------------------------------

    @Test
    public void taskCountUncappedIsOneForNonPositiveColumns() {
        assertEquals(1, GemmParallelSupport.taskCountUncapped(0, NC));
        assertEquals(1, GemmParallelSupport.taskCountUncapped(-1, NC));
    }

    @Test
    public void taskCountUncappedIsOneForExactBlockSize() {
        assertEquals(1, GemmParallelSupport.taskCountUncapped(NC, NC));
    }

    @Test
    public void taskCountUncappedIsNotArtificiallyCapped() {
        // With a very large column count the result must exceed MAX_PARALLEL_TASKS (4)
        // on any machine with more than 4 logical CPUs.
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > GemmParallelSupport.MAX_PARALLEL_TASKS);
        int largeN = NC * (Runtime.getRuntime().availableProcessors() + 1);
        assertTrue(GemmParallelSupport.taskCountUncapped(largeN, NC)
                > GemmParallelSupport.MAX_PARALLEL_TASKS);
    }

    @Test
    public void blockRangesCoverAllColumnsWithoutOverlapForBLISPartitioning() {
        // 2-task split: n = NC + 1
        assertCoverageWithoutOverlap(NC + 1, NC, 2);
        // 3-task split: n = 2*NC + 1
        int n3 = 2 * NC + 1;
        assertCoverageWithoutOverlap(n3, NC, Math.min(Runtime.getRuntime().availableProcessors(),
                GemmParallelSupport.taskCountUncapped(n3, NC)));
        // 4-task split: n = 3*NC + 1
        int n4 = 3 * NC + 1;
        assertCoverageWithoutOverlap(n4, NC, Math.min(Runtime.getRuntime().availableProcessors(),
                GemmParallelSupport.taskCountUncapped(n4, NC)));
    }

    // =======================================================================
    // Private helpers
    // =======================================================================

    /**
     * Compares the result of the sequential {@code sgemm} call against the
     * parallel (executor-gated) call. The two results must be bit-identical.
     */
    private static void verifySgemmCase(int m, int n, int k,
            float alpha, float beta, ExecutorService executor) {
        int lda = m + 2;
        int ldb = k + 3;
        int ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;

        float[] a        = new float[offA + lda * k];
        float[] b        = new float[offB + ldb * n];
        float[] expected = new float[offC + ldc * n + 9];
        float[] actual   = new float[expected.length];

        fill(a);
        fill(b);
        fill(expected);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, expected, 1, ldc);
        SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual,   1, ldc,
                executor);

        assertArrayEquals(expected, actual, 0.0f);
    }

    private static void verifyDgemmCase(int m, int n, int k,
            double alpha, double beta, ExecutorService executor) {
        int lda = m + 2;
        int ldb = k + 3;
        int ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;

        double[] a        = new double[offA + lda * k];
        double[] b        = new double[offB + ldb * n];
        double[] expected = new double[offC + ldc * n + 9];
        double[] actual   = new double[expected.length];

        fill(a);
        fill(b);
        fill(expected);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, expected, 1, ldc);
        DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual,   1, ldc,
                executor);

        assertArrayEquals(expected, actual, 0.0);
    }

    private static void verifySgemmCaseAgainstReference(int m, int n, int k,
            float alpha, float beta, ExecutorService executor, float delta) {
        int lda = m + 2;
        int ldb = k + 3;
        int ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(seedForCase(m, n, k, Float.floatToIntBits(alpha), Float.floatToIntBits(beta)));

        float[] a        = new float[offA + lda * k + 5];
        float[] b        = new float[offB + ldb * n + 5];
        float[] expected = new float[offC + ldc * n + 11];
        float[] actual   = new float[expected.length];

        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(expected, random);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        referenceSgemm(m, n, k, alpha, a, offA, lda, b, offB, ldb, beta, expected, offC, ldc);
        SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual, 1, ldc,
                executor);

        assertArrayEquals(expected, actual, delta);
    }

    private static void verifyDgemmCaseAgainstReference(int m, int n, int k,
            double alpha, double beta, ExecutorService executor, double delta) {
        int lda = m + 2;
        int ldb = k + 3;
        int ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(seedForCase(m, n, k, Double.doubleToLongBits(alpha), Double.doubleToLongBits(beta)));

        double[] a        = new double[offA + lda * k + 5];
        double[] b        = new double[offB + ldb * n + 5];
        double[] expected = new double[offC + ldc * n + 11];
        double[] actual   = new double[expected.length];

        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(expected, random);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        referenceDgemm(m, n, k, alpha, a, offA, lda, b, offB, ldb, beta, expected, offC, ldc);
        DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual, 1, ldc,
                executor);

        assertArrayEquals(expected, actual, delta);
    }

    private static void verifySgemmAlphaZeroCase(float beta) throws Exception {
        int m = 17, n = 33, k = 19;
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;

        float[] a      = new float[offA + lda * k];
        float[] b      = new float[offB + ldb * n];
        float[] actual = new float[offC + ldc * n + 11];
        float[] expect = new float[actual.length];

        fill(a);
        fill(b);
        fill(actual);
        System.arraycopy(actual, 0, expect, 0, actual.length);
        applyExpectedAlphaZero(expect, m, n, beta, offC, ldc);

        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            SgemmMRxNR.sgemm(m, n, k, 0.0f, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual, 1, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }
        // alpha == 0 is an early exit; the executor must never be touched
        assertEquals(0, executor.getSubmittedTaskCount());
        assertArrayEquals(expect, actual, 0.0f);
    }

    private static void verifyDgemmAlphaZeroCase(double beta) throws Exception {
        int m = 17, n = 33, k = 19;
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;

        double[] a      = new double[offA + lda * k];
        double[] b      = new double[offB + ldb * n];
        double[] actual = new double[offC + ldc * n + 11];
        double[] expect = new double[actual.length];

        fill(a);
        fill(b);
        fill(actual);
        System.arraycopy(actual, 0, expect, 0, actual.length);
        applyExpectedAlphaZero(expect, m, n, beta, offC, ldc);

        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(8));
        try {
            DgemmMRxNR.dgemm(m, n, k, 0.0, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, actual, 1, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
        assertArrayEquals(expect, actual, 0.0);
    }

    private static void verifySgemmPreservesInputAndPadding(int m, int n, int k,
            float alpha, float beta) throws Exception {
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(17L);

        float[] a = new float[offA + lda * k + 9];
        float[] b = new float[offB + ldb * n + 9];
        float[] c = new float[offC + ldc * n + 13];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(c, random);

        float[] origA = a.clone();
        float[] origB = b.clone();
        float[] origC = c.clone();

        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, c, 1, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }

        assertArrayEquals(origA, a, 0.0f);
        assertArrayEquals(origB, b, 0.0f);
        assertPaddingUntouched(origC, c, m, n, offC, ldc);
    }

    private static void verifyDgemmPreservesInputAndPadding(int m, int n, int k,
            double alpha, double beta) throws Exception {
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(23L);

        double[] a = new double[offA + lda * k + 9];
        double[] b = new double[offB + ldb * n + 9];
        double[] c = new double[offC + ldc * n + 13];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(c, random);

        double[] origA = a.clone();
        double[] origB = b.clone();
        double[] origC = c.clone();

        ExecutorService executor = Executors.newFixedThreadPool(8);
        try {
            DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, c, 1, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }

        assertArrayEquals(origA, a, 0.0);
        assertArrayEquals(origB, b, 0.0);
        assertPaddingUntouched(origC, c, m, n, offC, ldc);
    }

    private static void verifyRepeatedSgemmCase(int m, int n, int k,
            float alpha, float beta, ExecutorService executor) {
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(seedForCase(m, n, k, Float.floatToIntBits(alpha), Float.floatToIntBits(beta)));

        float[] a     = new float[offA + lda * k + 5];
        float[] b     = new float[offB + ldb * n + 5];
        float[] baseC = new float[offC + ldc * n + 11];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(baseC, random);

        float[] r1 = baseC.clone();
        float[] r2 = baseC.clone();

        SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, r1, 1, ldc, executor);
        SgemmMRxNR.sgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, r2, 1, ldc, executor);

        assertArrayEquals(r1, r2, 0.0f);
    }

    private static void verifyRepeatedDgemmCase(int m, int n, int k,
            double alpha, double beta, ExecutorService executor) {
        int lda = m + 2, ldb = k + 3, ldc = m + 4;
        int offA = 3, offB = 5, offC = 7;
        Random random = new Random(seedForCase(m, n, k, Double.doubleToLongBits(alpha), Double.doubleToLongBits(beta)));

        double[] a     = new double[offA + lda * k + 5];
        double[] b     = new double[offB + ldb * n + 5];
        double[] baseC = new double[offC + ldc * n + 11];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(baseC, random);

        double[] r1 = baseC.clone();
        double[] r2 = baseC.clone();

        DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, r1, 1, ldc, executor);
        DgemmMRxNR.dgemm(m, n, k, alpha, offA, a, 1, lda, offB, b, 1, ldb, beta, offC, r2, 1, ldc, executor);

        assertArrayEquals(r1, r2, 0.0);
    }

    // -----------------------------------------------------------------------
    // Naive reference implementations (column-major, no transpose)
    // -----------------------------------------------------------------------

    private static void referenceSgemm(int m, int n, int k,
            float alpha, float[] a, int offA, int lda,
            float[] b, int offB, int ldb,
            float beta, float[] c, int offC, int ldc) {
        for (int col = 0; col < n; col++) {
            int cBase = offC + col * ldc;
            for (int row = 0; row < m; row++) {
                float sum = 0.0f;
                for (int p = 0; p < k; p++) {
                    sum += a[offA + row + p * lda] * b[offB + p + col * ldb];
                }
                c[cBase + row] = alpha * sum + beta * c[cBase + row];
            }
        }
    }

    private static void referenceDgemm(int m, int n, int k,
            double alpha, double[] a, int offA, int lda,
            double[] b, int offB, int ldb,
            double beta, double[] c, int offC, int ldc) {
        for (int col = 0; col < n; col++) {
            int cBase = offC + col * ldc;
            for (int row = 0; row < m; row++) {
                double sum = 0.0;
                for (int p = 0; p < k; p++) {
                    sum += a[offA + row + p * lda] * b[offB + p + col * ldb];
                }
                c[cBase + row] = alpha * sum + beta * c[cBase + row];
            }
        }
    }

    // -----------------------------------------------------------------------
    // Assertion helpers
    // -----------------------------------------------------------------------

    private static void applyExpectedAlphaZero(float[] c, int m, int n,
            float beta, int offC, int ldc) {
        for (int col = 0; col < n; col++) {
            int base = offC + col * ldc;
            for (int row = 0; row < m; row++) {
                c[base + row] = (beta == 0.0f) ? 0.0f : beta * c[base + row];
            }
        }
    }

    private static void applyExpectedAlphaZero(double[] c, int m, int n,
            double beta, int offC, int ldc) {
        for (int col = 0; col < n; col++) {
            int base = offC + col * ldc;
            for (int row = 0; row < m; row++) {
                c[base + row] = (beta == 0.0) ? 0.0 : beta * c[base + row];
            }
        }
    }

    private static void assertPaddingUntouched(float[] expected, float[] actual,
            int m, int n, int offC, int ldc) {
        for (int index = 0; index < actual.length; index++) {
            if (!isInsideMatrix(index, m, n, offC, ldc)) {
                assertEquals(expected[index], actual[index], 0.0f);
            }
        }
    }

    private static void assertPaddingUntouched(double[] expected, double[] actual,
            int m, int n, int offC, int ldc) {
        for (int index = 0; index < actual.length; index++) {
            if (!isInsideMatrix(index, m, n, offC, ldc)) {
                assertEquals(expected[index], actual[index], 0.0);
            }
        }
    }

    private static boolean isInsideMatrix(int index, int m, int n, int offC, int ldc) {
        int relative = index - offC;
        if (relative < 0) {
            return false;
        }
        int col = relative / ldc;
        int row = relative % ldc;
        return col >= 0 && col < n && row >= 0 && row < m;
    }

    private static void assertCoverageWithoutOverlap(int columns, int blockSize, int taskCount) {
        int previousEnd = 0;
        for (int taskIndex = 0; taskIndex < taskCount; taskIndex++) {
            int start = GemmParallelSupport.blockStart(taskIndex, taskCount, columns, blockSize);
            int end   = GemmParallelSupport.blockEnd(taskIndex, taskCount, columns, blockSize);
            assertEquals(previousEnd, start);
            assertTrue(end >= start);
            previousEnd = end;
        }
        assertEquals(columns, previousEnd);
    }

    // -----------------------------------------------------------------------
    // Utility
    // -----------------------------------------------------------------------

    private static long seedForCase(int m, int n, int k, long scalarBits0, long scalarBits1) {
        long seed = 1469598103934665603L;
        seed = mix(seed, m);
        seed = mix(seed, n);
        seed = mix(seed, k);
        seed = mix(seed, scalarBits0);
        seed = mix(seed, scalarBits1);
        return seed;
    }

    private static long mix(long seed, long value) {
        return (seed ^ value) * 1099511628211L;
    }

    private static float pickFloatScalar(Random random) {
        switch (random.nextInt(6)) {
        case 0:  return 0.0f;
        case 1:  return 1.0f;
        case 2:  return -1.0f;
        case 3:  return 0.75f;
        case 4:  return -0.25f;
        default: return 1.5f;
        }
    }

    private static double pickDoubleScalar(Random random) {
        switch (random.nextInt(6)) {
        case 0:  return 0.0;
        case 1:  return 1.0;
        case 2:  return -1.0;
        case 3:  return -1.25;
        case 4:  return 0.5;
        default: return 2.0;
        }
    }

    private static void fill(float[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (((i * 37) % 101) - 50) / 17.0f;
        }
    }

    private static void fill(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (((i * 53) % 211) - 105) / 31.0;
        }
    }

    private static void fillRandom(float[] array, Random random) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (random.nextInt(2001) - 1000) / 137.0f;
        }
    }

    private static void fillRandom(double[] array, Random random) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (random.nextInt(4001) - 2000) / 257.0;
        }
    }

    private static void shutdown(ExecutorService executor) throws Exception {
        if (executor == null) {
            return;
        }
        executor.shutdown();
        if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
            executor.shutdownNow();
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                throw new AssertionError("executor did not terminate");
            }
        }
    }

    // -----------------------------------------------------------------------
    // CountingExecutor
    // -----------------------------------------------------------------------

    private static final class CountingExecutor extends AbstractExecutorService {

        private final ExecutorService delegate;
        private final AtomicInteger submittedTaskCount = new AtomicInteger();

        private CountingExecutor(ExecutorService delegate) {
            this.delegate = delegate;
        }

        int getSubmittedTaskCount() {
            return submittedTaskCount.get();
        }

        @Override
        public void shutdown() {
            delegate.shutdown();
        }

        @Override
        public List<Runnable> shutdownNow() {
            List<Runnable> tasks = delegate.shutdownNow();
            return (tasks != null) ? tasks : Collections.<Runnable>emptyList();
        }

        @Override
        public boolean isShutdown() {
            return delegate.isShutdown();
        }

        @Override
        public boolean isTerminated() {
            return delegate.isTerminated();
        }

        @Override
        public boolean awaitTermination(long timeout, TimeUnit unit) throws InterruptedException {
            return delegate.awaitTermination(timeout, unit);
        }

        @Override
        public void execute(Runnable command) {
            submittedTaskCount.incrementAndGet();
            delegate.execute(command);
        }
    }
}
