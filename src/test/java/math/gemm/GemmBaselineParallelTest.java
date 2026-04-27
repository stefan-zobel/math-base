package math.gemm;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

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

public class GemmBaselineParallelTest {

    private static final int M = 9;
    private static final int N = 2051;
    private static final int K = 57;

    @Test
    public void sgemmParallelMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifySgemm(true, true);
        verifySgemm(true, false);
        verifySgemm(false, true);
        verifySgemm(false, false);
    }

    @Test
    public void dgemmParallelMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifyDgemm(true, true);
        verifyDgemm(true, false);
        verifyDgemm(false, true);
        verifyDgemm(false, false);
    }

    @Test
    public void sgemmActualParallelMatchesSequentialForAllTransposeCombinations() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllSgemmCases(128, 4097, 86, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmActualParallelMatchesSequentialForAllTransposeCombinations() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllDgemmCases(M, N, K, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void sgemmExecutorFallsBackToSequentialBelowThreshold() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllSgemmCases(M, N, K, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmExecutorFallsBackToSequentialBelowThreshold() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllDgemmCases(9, 511, 57, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void sgemmNullExecutorMatchesSequential() throws Exception {
        verifyAllSgemmCases(31, 2051, 57, 0.75f, -0.25f, null);
    }

    @Test
    public void dgemmNullExecutorMatchesSequential() throws Exception {
        verifyAllDgemmCases(31, 2051, 57, -1.25, 0.5, null);
    }

    @Test
    public void sgemmShutdownExecutorFallsBackToSequential() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        executor.shutdown();
        try {
            verifyAllSgemmCases(128, 2051, 86, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmShutdownExecutorFallsBackToSequential() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        executor.shutdown();
        try {
            verifyAllDgemmCases(128, 2051, 86, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

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

    @Test
    public void sgemmAlphaOneBetaZeroMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifyAllSgemmCases(33, 2051, 57, 1.0f, 0.0f, null);
    }

    @Test
    public void sgemmAlphaOneBetaOneMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifyAllSgemmCases(33, 2051, 57, 1.0f, 1.0f, null);
    }

    @Test
    public void dgemmAlphaOneBetaZeroMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifyAllDgemmCases(33, 2051, 57, 1.0, 0.0, null);
    }

    @Test
    public void dgemmAlphaOneBetaOneMatchesSequentialForAllTransposeCombinations() throws Exception {
        verifyAllDgemmCases(33, 2051, 57, 1.0, 1.0, null);
    }

    @Test
    public void sgemmBlockedPathMatchesSequentialForNnAndNt() throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            verifySgemmCase(129, 2051, 257, true, true, 0.75f, -0.25f, executor);
            verifySgemmCase(129, 2051, 257, true, false, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void dgemmBlockedPathMatchesSequentialForNnAndNt() throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            verifyDgemmCase(129, 2051, 257, true, true, -1.25, 0.5, executor);
            verifyDgemmCase(129, 2051, 257, true, false, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void taskCountIsCappedAtFour() {
        assertTrue(GemmParallelSupport.taskCount(Integer.MAX_VALUE, 1) <= GemmParallelSupport.MAX_PARALLEL_TASKS);
    }

    @Test
    public void taskCountIsOneForNonPositiveColumns() {
        assertEquals(1, GemmParallelSupport.taskCount(0, 2048));
        assertEquals(1, GemmParallelSupport.taskCount(-1, 2048));
    }

    @Test
    public void taskCountIsOneForExactBlockSize() {
        assertEquals(1, GemmParallelSupport.taskCount(2048, 2048));
        assertEquals(1, GemmParallelSupport.taskCount(1024, 1024));
    }

    @Test
    public void sgemmExactNcBoundaryUsesSingleTaskFallback() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifySgemmCase(128, 2048, 86, true, true, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void dgemmExactNcBoundaryUsesSingleTaskFallback() throws Exception {
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyDgemmCase(128, 1024, 8, true, true, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
    }

    @Test
    public void sgemmNcPlusOneSubmitsParallelTasks() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifySgemmCase(128, 4097, 86, true, true, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmNcPlusOneSubmitsParallelTasks() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyDgemmCase(128, 2049, 8, true, true, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void blockRangesCoverAllColumnsWithoutOverlapForSgemmLikePartitioning() {
        assertExactPartition(2049, 2048, 2, new int[] { 0, 2048, 2048, 2049 });
        assertCoverageWithoutOverlap(2049, 2048, 2);
        assertCoverageWithoutOverlap(6145, 2048, Math.min(4, GemmParallelSupport.taskCount(6145, 2048)));
    }

    @Test
    public void blockRangesCoverAllColumnsWithoutOverlapForDgemmLikePartitioning() {
        assertExactPartition(2051, 1024, 3, new int[] { 0, 1024, 1024, 2048, 2048, 2051 });
        assertCoverageWithoutOverlap(2051, 1024, 3);
        assertCoverageWithoutOverlap(6145, 1024, Math.min(4, GemmParallelSupport.taskCount(6145, 1024)));
    }

    @Test
    public void sgemmDegenerateDimensionsMatchNaiveReference() {
        verifyAllSgemmCasesAgainstReference(0, 7, 5, 0.75f, -0.25f, null, 0.0f);
        verifyAllSgemmCasesAgainstReference(5, 0, 7, 0.75f, -0.25f, null, 0.0f);
        verifyAllSgemmCasesAgainstReference(5, 7, 0, 0.75f, -0.25f, null, 0.0f);
    }

    @Test
    public void dgemmDegenerateDimensionsMatchNaiveReference() {
        verifyAllDgemmCasesAgainstReference(0, 7, 5, -1.25, 0.5, null, 0.0);
        verifyAllDgemmCasesAgainstReference(5, 0, 7, -1.25, 0.5, null, 0.0);
        verifyAllDgemmCasesAgainstReference(5, 7, 0, -1.25, 0.5, null, 0.0);
    }

    @Test
    public void sgemmDoesNotModifyInputsOrPadding() throws Exception {
        verifySgemmPreservesInputAndPadding(13, 17, 11, false, false, 0.75f, -0.25f);
    }

    @Test
    public void dgemmDoesNotModifyInputsOrPadding() throws Exception {
        verifyDgemmPreservesInputAndPadding(13, 17, 11, false, false, -1.25, 0.5);
    }

    @Test
    public void sgemmRepeatedCallsWithSharedExecutorRemainDeterministic() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyRepeatedSgemmCase(43, 4097, 256, true, false, 0.75f, -0.25f, executor);
            verifyRepeatedSgemmCase(43, 4097, 256, false, true, 1.0f, 1.0f, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmRepeatedCallsWithSharedExecutorRemainDeterministic() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyRepeatedDgemmCase(31, 2049, 32, true, false, -1.25, 0.5, executor);
            verifyRepeatedDgemmCase(31, 2049, 32, false, true, 1.0, 1.0, executor);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void sgemmRandomizedSmallCasesMatchNaiveReference() throws Exception {
        Random random = new Random(0x5EED5EEDL);
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            for (int caseIndex = 0; caseIndex < 16; caseIndex++) {
                int m = random.nextInt(8);
                int n = random.nextInt(8);
                int k = random.nextInt(8);
                float alpha = pickFloatScalar(random);
                float beta = pickFloatScalar(random);
                verifyAllSgemmCasesAgainstReference(m, n, k, alpha, beta, (caseIndex & 1) == 0 ? null : executor, 1.0e-4f);
            }
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void dgemmRandomizedSmallCasesMatchNaiveReference() throws Exception {
        Random random = new Random(0xD6EED6EEL);
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            for (int caseIndex = 0; caseIndex < 16; caseIndex++) {
                int m = random.nextInt(8);
                int n = random.nextInt(8);
                int k = random.nextInt(8);
                double alpha = pickDoubleScalar(random);
                double beta = pickDoubleScalar(random);
                verifyAllDgemmCasesAgainstReference(m, n, k, alpha, beta, (caseIndex & 1) == 0 ? null : executor, 1.0e-10);
            }
        } finally {
            shutdown(executor);
        }
    }

    @Test
    public void sgemmHighWorkCasesMatchNaiveReferenceForAllTransposeCombinations() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllSgemmCasesAgainstReference(43, 4097, 256, 0.75f, -0.25f, executor, 1.0e-3f);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    @Test
    public void dgemmHighWorkCasesMatchNaiveReferenceForAllTransposeCombinations() throws Exception {
        Assume.assumeTrue(Runtime.getRuntime().availableProcessors() > 1);
        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            verifyAllDgemmCasesAgainstReference(31, 2049, 32, -1.25, 0.5, executor, 1.0e-10);
        } finally {
            shutdown(executor);
        }
        assertTrue(executor.getSubmittedTaskCount() > 0);
    }

    private static void verifySgemm(boolean notATransposed, boolean notBTransposed) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            verifySgemmCase(M, N, K, notATransposed, notBTransposed, 0.75f, -0.25f, executor);
        } finally {
            shutdown(executor);
        }
    }

    private static void verifyDgemm(boolean notATransposed, boolean notBTransposed) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            verifyDgemmCase(M, N, K, notATransposed, notBTransposed, -1.25, 0.5, executor);
        } finally {
            shutdown(executor);
        }
    }

    private static void verifyAllSgemmCases(int m, int n, int k, float alpha, float beta, ExecutorService executor) throws Exception {
        verifySgemmCase(m, n, k, true, true, alpha, beta, executor);
        verifySgemmCase(m, n, k, true, false, alpha, beta, executor);
        verifySgemmCase(m, n, k, false, true, alpha, beta, executor);
        verifySgemmCase(m, n, k, false, false, alpha, beta, executor);
    }

    private static void verifyAllDgemmCases(int m, int n, int k, double alpha, double beta, ExecutorService executor) throws Exception {
        verifyDgemmCase(m, n, k, true, true, alpha, beta, executor);
        verifyDgemmCase(m, n, k, true, false, alpha, beta, executor);
        verifyDgemmCase(m, n, k, false, true, alpha, beta, executor);
        verifyDgemmCase(m, n, k, false, false, alpha, beta, executor);
    }

    private static void verifyAllSgemmCasesAgainstReference(int m, int n, int k, float alpha, float beta, ExecutorService executor,
            float delta) {
        verifySgemmCaseAgainstReference(m, n, k, true, true, alpha, beta, executor, delta);
        verifySgemmCaseAgainstReference(m, n, k, true, false, alpha, beta, executor, delta);
        verifySgemmCaseAgainstReference(m, n, k, false, true, alpha, beta, executor, delta);
        verifySgemmCaseAgainstReference(m, n, k, false, false, alpha, beta, executor, delta);
    }

    private static void verifyAllDgemmCasesAgainstReference(int m, int n, int k, double alpha, double beta, ExecutorService executor,
            double delta) {
        verifyDgemmCaseAgainstReference(m, n, k, true, true, alpha, beta, executor, delta);
        verifyDgemmCaseAgainstReference(m, n, k, true, false, alpha, beta, executor, delta);
        verifyDgemmCaseAgainstReference(m, n, k, false, true, alpha, beta, executor, delta);
        verifyDgemmCaseAgainstReference(m, n, k, false, false, alpha, beta, executor, delta);
    }

    private static void verifySgemmCase(int m, int n, int k, boolean notATransposed, boolean notBTransposed, float alpha, float beta,
            ExecutorService executor) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;

        float[] a = new float[aOffset + lda * aColumns];
        float[] b = new float[bOffset + ldb * bColumns];
        float[] expected = new float[cOffset + ldc * n + 9];
        float[] actual = new float[expected.length];

        fill(a);
        fill(b);
        fill(expected);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, expected, cOffset, ldc);
        SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                executor);

        assertArrayEquals(expected, actual, 0.0f);
    }

    private static void verifyDgemmCase(int m, int n, int k, boolean notATransposed, boolean notBTransposed, double alpha, double beta,
            ExecutorService executor) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;

        double[] a = new double[aOffset + lda * aColumns];
        double[] b = new double[bOffset + ldb * bColumns];
        double[] expected = new double[cOffset + ldc * n + 9];
        double[] actual = new double[expected.length];

        fill(a);
        fill(b);
        fill(expected);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, expected, cOffset, ldc);
        DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                executor);

        assertArrayEquals(expected, actual, 0.0);
    }

    private static void verifySgemmCaseAgainstReference(int m, int n, int k, boolean notATransposed, boolean notBTransposed, float alpha,
            float beta, ExecutorService executor, float delta) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;

        float[] a = new float[aOffset + lda * aColumns + 5];
        float[] b = new float[bOffset + ldb * bColumns + 5];
        float[] expected = new float[cOffset + ldc * n + 11];
        float[] actual = new float[expected.length];
        Random random = new Random(seedForCase(m, n, k, notATransposed, notBTransposed, Float.floatToIntBits(alpha), Float.floatToIntBits(beta)));

        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(expected, random);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        referenceSgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, expected, cOffset, ldc);
        SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                executor);

        assertArrayEquals(expected, actual, delta);
    }

    private static void verifyDgemmCaseAgainstReference(int m, int n, int k, boolean notATransposed, boolean notBTransposed, double alpha,
            double beta, ExecutorService executor, double delta) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;

        double[] a = new double[aOffset + lda * aColumns + 5];
        double[] b = new double[bOffset + ldb * bColumns + 5];
        double[] expected = new double[cOffset + ldc * n + 11];
        double[] actual = new double[expected.length];
        Random random = new Random(seedForCase(m, n, k, notATransposed, notBTransposed, Double.doubleToLongBits(alpha), Double.doubleToLongBits(beta)));

        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(expected, random);
        System.arraycopy(expected, 0, actual, 0, expected.length);

        referenceDgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, expected, cOffset, ldc);
        DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                executor);

        assertArrayEquals(expected, actual, delta);
    }

    private static void verifySgemmPreservesInputAndPadding(int m, int n, int k, boolean notATransposed, boolean notBTransposed, float alpha,
            float beta) throws Exception {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        Random random = new Random(17L);

        float[] a = new float[aOffset + lda * aColumns + 9];
        float[] b = new float[bOffset + ldb * bColumns + 9];
        float[] c = new float[cOffset + ldc * n + 13];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(c, random);

        float[] originalA = a.clone();
        float[] originalB = b.clone();
        float[] originalC = c.clone();

        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }

        assertArrayEquals(originalA, a, 0.0f);
        assertArrayEquals(originalB, b, 0.0f);
        assertPaddingUntouched(originalC, c, m, n, cOffset, ldc);
    }

    private static void verifyDgemmPreservesInputAndPadding(int m, int n, int k, boolean notATransposed, boolean notBTransposed, double alpha,
            double beta) throws Exception {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        Random random = new Random(23L);

        double[] a = new double[aOffset + lda * aColumns + 9];
        double[] b = new double[bOffset + ldb * bColumns + 9];
        double[] c = new double[cOffset + ldc * n + 13];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(c, random);

        double[] originalA = a.clone();
        double[] originalB = b.clone();
        double[] originalC = c.clone();

        ExecutorService executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
        try {
            DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }

        assertArrayEquals(originalA, a, 0.0);
        assertArrayEquals(originalB, b, 0.0);
        assertPaddingUntouched(originalC, c, m, n, cOffset, ldc);
    }

    private static void verifyRepeatedSgemmCase(int m, int n, int k, boolean notATransposed, boolean notBTransposed, float alpha, float beta,
            ExecutorService executor) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        Random random = new Random(seedForCase(m, n, k, notATransposed, notBTransposed, Float.floatToIntBits(alpha), Float.floatToIntBits(beta)));

        float[] a = new float[aOffset + lda * aColumns + 5];
        float[] b = new float[bOffset + ldb * bColumns + 5];
        float[] baseC = new float[cOffset + ldc * n + 11];
        float[] actual1 = new float[baseC.length];
        float[] actual2 = new float[baseC.length];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(baseC, random);
        System.arraycopy(baseC, 0, actual1, 0, baseC.length);
        System.arraycopy(baseC, 0, actual2, 0, baseC.length);

        SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual1, cOffset, ldc,
                executor);
        SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual2, cOffset, ldc,
                executor);

        assertArrayEquals(actual1, actual2, 0.0f);
    }

    private static void verifyRepeatedDgemmCase(int m, int n, int k, boolean notATransposed, boolean notBTransposed, double alpha, double beta,
            ExecutorService executor) {
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        Random random = new Random(seedForCase(m, n, k, notATransposed, notBTransposed, Double.doubleToLongBits(alpha), Double.doubleToLongBits(beta)));

        double[] a = new double[aOffset + lda * aColumns + 5];
        double[] b = new double[bOffset + ldb * bColumns + 5];
        double[] baseC = new double[cOffset + ldc * n + 11];
        double[] actual1 = new double[baseC.length];
        double[] actual2 = new double[baseC.length];
        fillRandom(a, random);
        fillRandom(b, random);
        fillRandom(baseC, random);
        System.arraycopy(baseC, 0, actual1, 0, baseC.length);
        System.arraycopy(baseC, 0, actual2, 0, baseC.length);

        DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual1, cOffset, ldc,
                executor);
        DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, actual2, cOffset, ldc,
                executor);

        assertArrayEquals(actual1, actual2, 0.0);
    }

    private static void verifySgemmAlphaZeroCase(float beta) throws Exception {
        boolean notATransposed = false;
        boolean notBTransposed = false;
        int m = 17;
        int n = 33;
        int k = 19;
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;

        float[] a = new float[aOffset + lda * aColumns];
        float[] b = new float[bOffset + ldb * bColumns];
        float[] actual = new float[cOffset + ldc * n + 11];
        float[] expected = new float[actual.length];

        fill(a);
        fill(b);
        fill(actual);
        System.arraycopy(actual, 0, expected, 0, actual.length);
        applyExpectedAlphaZero(expected, m, n, beta, cOffset, ldc);

        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            SgemmBaseline.sgemm(notATransposed, notBTransposed, m, n, k, 0.0f, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
        assertArrayEquals(expected, actual, 0.0f);
    }

    private static void verifyDgemmAlphaZeroCase(double beta) throws Exception {
        boolean notATransposed = false;
        boolean notBTransposed = false;
        int m = 17;
        int n = 33;
        int k = 19;
        int lda = (notATransposed ? m : k) + 2;
        int ldb = (notBTransposed ? k : n) + 3;
        int ldc = m + 4;
        int aOffset = 3;
        int bOffset = 5;
        int cOffset = 7;
        int aColumns = notATransposed ? k : m;
        int bColumns = notBTransposed ? n : k;

        double[] a = new double[aOffset + lda * aColumns];
        double[] b = new double[bOffset + ldb * bColumns];
        double[] actual = new double[cOffset + ldc * n + 11];
        double[] expected = new double[actual.length];

        fill(a);
        fill(b);
        fill(actual);
        System.arraycopy(actual, 0, expected, 0, actual.length);
        applyExpectedAlphaZero(expected, m, n, beta, cOffset, ldc);

        CountingExecutor executor = new CountingExecutor(Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS));
        try {
            DgemmBaseline.dgemm(notATransposed, notBTransposed, m, n, k, 0.0, a, aOffset, lda, b, bOffset, ldb, beta, actual, cOffset, ldc,
                    executor);
        } finally {
            shutdown(executor);
        }
        assertEquals(0, executor.getSubmittedTaskCount());
        assertArrayEquals(expected, actual, 0.0);
    }

    private static void applyExpectedAlphaZero(float[] c, int m, int n, float beta, int cOffset, int ldc) {
        for (int col = 0; col < n; col++) {
            int base = cOffset + col * ldc;
            for (int row = 0; row < m; row++) {
                c[base + row] = (beta == 0.0f) ? 0.0f : beta * c[base + row];
            }
        }
    }

    private static void applyExpectedAlphaZero(double[] c, int m, int n, double beta, int cOffset, int ldc) {
        for (int col = 0; col < n; col++) {
            int base = cOffset + col * ldc;
            for (int row = 0; row < m; row++) {
                c[base + row] = (beta == 0.0) ? 0.0 : beta * c[base + row];
            }
        }
    }

    private static void assertExactPartition(int columns, int blockSize, int taskCount, int[] expectedRanges) {
        if (expectedRanges.length != taskCount * 2) {
            fail("expectedRanges must contain start/end pairs for every task");
        }
        for (int taskIndex = 0; taskIndex < taskCount; taskIndex++) {
            assertEquals(expectedRanges[taskIndex * 2], GemmParallelSupport.blockStart(taskIndex, taskCount, columns, blockSize));
            assertEquals(expectedRanges[taskIndex * 2 + 1], GemmParallelSupport.blockEnd(taskIndex, taskCount, columns, blockSize));
        }
    }

    private static void assertCoverageWithoutOverlap(int columns, int blockSize, int taskCount) {
        int previousEnd = 0;
        for (int taskIndex = 0; taskIndex < taskCount; taskIndex++) {
            int start = GemmParallelSupport.blockStart(taskIndex, taskCount, columns, blockSize);
            int end = GemmParallelSupport.blockEnd(taskIndex, taskCount, columns, blockSize);
            assertEquals(previousEnd, start);
            assertTrue(end >= start);
            previousEnd = end;
        }
        assertEquals(columns, previousEnd);
    }

    private static void assertPaddingUntouched(float[] expected, float[] actual, int m, int n, int cOffset, int ldc) {
        for (int index = 0; index < actual.length; index++) {
            if (!isInsideMatrix(index, m, n, cOffset, ldc)) {
                assertEquals(expected[index], actual[index], 0.0f);
            }
        }
    }

    private static void assertPaddingUntouched(double[] expected, double[] actual, int m, int n, int cOffset, int ldc) {
        for (int index = 0; index < actual.length; index++) {
            if (!isInsideMatrix(index, m, n, cOffset, ldc)) {
                assertEquals(expected[index], actual[index], 0.0);
            }
        }
    }

    private static boolean isInsideMatrix(int index, int m, int n, int cOffset, int ldc) {
        int relative = index - cOffset;
        if (relative < 0) {
            return false;
        }
        int col = relative / ldc;
        int row = relative % ldc;
        return col >= 0 && col < n && row >= 0 && row < m;
    }

    private static void referenceSgemm(boolean notATransposed, boolean notBTransposed, int m, int n, int k, float alpha, float[] a,
            int aOffset, int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc) {
        for (int col = 0; col < n; col++) {
            int cColOffset = cOffset + col * ldc;
            for (int row = 0; row < m; row++) {
                float sum = 0.0f;
                for (int p = 0; p < k; p++) {
                    float av = notATransposed ? a[aOffset + row + p * lda] : a[aOffset + p + row * lda];
                    float bv = notBTransposed ? b[bOffset + p + col * ldb] : b[bOffset + col + p * ldb];
                    sum += av * bv;
                }
                c[cColOffset + row] = alpha * sum + beta * c[cColOffset + row];
            }
        }
    }

    private static void referenceDgemm(boolean notATransposed, boolean notBTransposed, int m, int n, int k, double alpha, double[] a,
            int aOffset, int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        for (int col = 0; col < n; col++) {
            int cColOffset = cOffset + col * ldc;
            for (int row = 0; row < m; row++) {
                double sum = 0.0;
                for (int p = 0; p < k; p++) {
                    double av = notATransposed ? a[aOffset + row + p * lda] : a[aOffset + p + row * lda];
                    double bv = notBTransposed ? b[bOffset + p + col * ldb] : b[bOffset + col + p * ldb];
                    sum += av * bv;
                }
                c[cColOffset + row] = alpha * sum + beta * c[cColOffset + row];
            }
        }
    }

    private static long seedForCase(int m, int n, int k, boolean notATransposed, boolean notBTransposed, long scalarBits0, long scalarBits1) {
        long seed = 1469598103934665603L;
        seed = mix(seed, m);
        seed = mix(seed, n);
        seed = mix(seed, k);
        seed = mix(seed, notATransposed ? 1 : 0);
        seed = mix(seed, notBTransposed ? 1 : 0);
        seed = mix(seed, scalarBits0);
        seed = mix(seed, scalarBits1);
        return seed;
    }

    private static long mix(long seed, long value) {
        return (seed ^ value) * 1099511628211L;
    }

    private static float pickFloatScalar(Random random) {
        switch (random.nextInt(6)) {
        case 0:
            return 0.0f;
        case 1:
            return 1.0f;
        case 2:
            return -1.0f;
        case 3:
            return 0.75f;
        case 4:
            return -0.25f;
        default:
            return 1.5f;
        }
    }

    private static double pickDoubleScalar(Random random) {
        switch (random.nextInt(6)) {
        case 0:
            return 0.0;
        case 1:
            return 1.0;
        case 2:
            return -1.0;
        case 3:
            return -1.25;
        case 4:
            return 0.5;
        default:
            return 2.0;
        }
    }

    private static void fill(float[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (((i * 37) % 101) - 50) / 17.0f;
        }
    }

    private static void fillRandom(float[] array, Random random) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (random.nextInt(2001) - 1000) / 137.0f;
        }
    }

    private static void fill(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = (((i * 53) % 211) - 105) / 31.0;
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
