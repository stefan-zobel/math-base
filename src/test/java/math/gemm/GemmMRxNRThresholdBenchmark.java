package math.gemm;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Threads;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Benchmarks the sequential BLIS-style GEMM kernels against the forced-parallel
 * N-splitting path to validate the parallelization heuristic.
 * <p>
 * Both {@link SgemmMRxNR} and {@link DgemmMRxNR} share {@code NC = 3840} as the
 * column-panel width. The {@code PARALLEL_WORK_THRESHOLD} differs between the
 * two jar versions:
 * <ul>
 *   <li><strong>Java-8 scalar kernel</strong>: threshold = 45&thinsp;000&thinsp;000</li>
 *   <li><strong>Java-25 Vector-API kernel</strong>: threshold = 25&thinsp;000&thinsp;000
 *       (vectorized micro-kernel is faster, so parallelism pays off earlier)</li>
 * </ul>
 * Parallelization requires both {@code taskCountUncapped(colsB, 3840) >= 3}
 * (i.e. {@code colsB >= 7681}) <em>and</em>
 * {@code rowsA * colsB * colsA >= PARALLEL_WORK_THRESHOLD}.
 * Unlike the baseline kernels there are no transpose combinations; all cases
 * use column-major, no-transpose storage ({@code incRow=1}, {@code incCol=lda}).
 *
 * <p><strong>Main benchmark set ({@link BenchmarkState}, {@code @Param}):</strong>
 * Fixed dimension strings that are meaningful on both Java versions:
 * <ol>
 *   <li>{@code 128x3841x128} (~63.5M ops) &mdash; 2 NC-blocks: {@code taskCount < 3}
 *       &rarr; always sequential regardless of threshold &mdash; negative control.</li>
 *   <li>{@code 128x7681x64}  (~63.1M ops) &mdash; 3 NC-blocks, moderate work
 *       (well above both thresholds).</li>
 *   <li>{@code 128x7681x128} (~126M ops)  &mdash; 3 NC-blocks, large work.</li>
 *   <li>{@code 128x11521x128} (~189M ops) &mdash; 4 NC-blocks (requires &ge;4 CPUs).</li>
 * </ol>
 *
 * <p><strong>Threshold-boundary benchmark set ({@link ThresholdBoundaryState}):</strong>
 * Dimensions are computed at runtime via {@code getParallelWorkThreshold()} so
 * that {@code kLo} and {@code kHi} always straddle the active threshold
 * regardless of which jar version is running.  With {@code m=64}, {@code n=7681}:
 * <pre>
 *   kLo = floor((threshold - 1) / (m * n))  ->  m * n * kLo  &lt; threshold (sequential)
 *   kHi = kLo + 1                            ->  m * n * kHi  &gt; threshold (parallel)
 * </pre>
 *
 * <p><strong>JMH settings rationale:</strong> measurement time is set to 2 s
 * to obtain enough samples for the slower large-matrix cases ({@literal >}10
 * ms/op) and thereby reduce error variance.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Threads(1)
@Fork(1)
@Warmup(iterations = 4, time = 2)
@Measurement(iterations = 5, time = 2)
public class GemmMRxNRThresholdBenchmark {

    private static final Method SGEMM_PARALLELIZE_OVER_N =
            lookupParallelizeOverN(SgemmMRxNR.class, float.class, float[].class);
    private static final Method DGEMM_PARALLELIZE_OVER_N =
            lookupParallelizeOverN(DgemmMRxNR.class, double.class, double[].class);

    // -----------------------------------------------------------------------
    // sgemm benchmarks — main set
    // -----------------------------------------------------------------------

    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void sgemmSequential(BenchmarkState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.k, state.sAlpha,
                state.offA, state.sa, 1, state.lda,
                state.offB, state.sb, 1, state.ldb,
                state.sBeta, state.offC, state.sc, 1, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void sgemmForcedParallel(BenchmarkState state, Blackhole bh) {
        invokeParallelizeOverN(SGEMM_PARALLELIZE_OVER_N,
                state.executor, state.m, state.n, state.k, state.sAlpha,
                state.offA, state.sa, 1, state.lda,
                state.offB, state.sb, 1, state.ldb,
                state.sBeta, state.offC, state.sc, 1, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void sgemmAuto(BenchmarkState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.k, state.sAlpha,
                state.offA, state.sa, 1, state.lda,
                state.offB, state.sb, 1, state.ldb,
                state.sBeta, state.offC, state.sc, 1, state.ldc,
                state.executor);
        consume(state.sc, bh);
    }

    // -----------------------------------------------------------------------
    // dgemm benchmarks — main set
    // -----------------------------------------------------------------------

    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void dgemmSequential(BenchmarkState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.k, state.dAlpha,
                state.offA, state.da, 1, state.lda,
                state.offB, state.db, 1, state.ldb,
                state.dBeta, state.offC, state.dc, 1, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmForcedParallel(BenchmarkState state, Blackhole bh) {
        invokeParallelizeOverN(DGEMM_PARALLELIZE_OVER_N,
                state.executor, state.m, state.n, state.k, state.dAlpha,
                state.offA, state.da, 1, state.lda,
                state.offB, state.db, 1, state.ldb,
                state.dBeta, state.offC, state.dc, 1, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmAuto(BenchmarkState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.k, state.dAlpha,
                state.offA, state.da, 1, state.lda,
                state.offB, state.db, 1, state.ldb,
                state.dBeta, state.offC, state.dc, 1, state.ldc,
                state.executor);
        consume(state.dc, bh);
    }

    // -----------------------------------------------------------------------
    // sgemm benchmarks — threshold boundary (version-adaptive dimensions)
    // -----------------------------------------------------------------------

    /** Sequential path with work just below the active threshold. */
    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void sgemmThresholdLower(ThresholdBoundaryState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.sKlo, state.sAlpha,
                state.offA, state.saLo, 1, state.sLdaLo,
                state.offB, state.sbLo, 1, state.sLdbLo,
                state.sBeta, state.offC, state.scLo, 1, state.sLdcLo,
                state.executor);
        consume(state.scLo, bh);
    }

    /** Auto path with work just above the active threshold. */
    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void sgemmThresholdUpper(ThresholdBoundaryState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.sKhi, state.sAlpha,
                state.offA, state.saHi, 1, state.sLdaHi,
                state.offB, state.sbHi, 1, state.sLdbHi,
                state.sBeta, state.offC, state.scHi, 1, state.sLdcHi,
                state.executor);
        consume(state.scHi, bh);
    }

    // -----------------------------------------------------------------------
    // dgemm benchmarks — threshold boundary (version-adaptive dimensions)
    // -----------------------------------------------------------------------

    /** Sequential path with work just below the active threshold. */
    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void dgemmThresholdLower(ThresholdBoundaryState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.dKlo, state.dAlpha,
                state.offA, state.daLo, 1, state.dLdaLo,
                state.offB, state.dbLo, 1, state.dLdbLo,
                state.dBeta, state.offC, state.dcLo, 1, state.dLdcLo,
                state.executor);
        consume(state.dcLo, bh);
    }

    /** Auto path with work just above the active threshold. */
    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void dgemmThresholdUpper(ThresholdBoundaryState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.dKhi, state.dAlpha,
                state.offA, state.daHi, 1, state.dLdaHi,
                state.offB, state.dbHi, 1, state.dLdbHi,
                state.dBeta, state.offC, state.dcHi, 1, state.dLdcHi,
                state.executor);
        consume(state.dcHi, bh);
    }

    // -----------------------------------------------------------------------
    // Main benchmark state (fixed @Param dimensions, valid for both JVM versions)
    // -----------------------------------------------------------------------

    @State(Scope.Thread)
    public static class BenchmarkState {

        /**
         * Fixed dimension strings in the form {@code MxNxK}.
         * <ul>
         *   <li>{@code 128x3841x128}: 2 NC-blocks &rarr; taskCount &lt; 3
         *       &rarr; always sequential (negative control for both versions).</li>
         *   <li>{@code 128x7681x64}:  3 NC-blocks, ~63M ops &mdash; well above
         *       both thresholds (25M and 45M).</li>
         *   <li>{@code 128x7681x128}: 3 NC-blocks, ~126M ops.</li>
         *   <li>{@code 128x11521x128}: 4 NC-blocks (requires &ge;4 CPUs).</li>
         * </ul>
         */
        @Param({
            "128x3841x128",   // 2 NC-blocks -> sequential (negative control)
            "128x7681x64",    // ~63.1M ops, 3 NC-blocks, above both thresholds
            "128x7681x128",   // ~125.9M ops, 3 NC-blocks
            "128x11521x128"   // ~188.7M ops, 4 NC-blocks
        })
        public String dimensions;

        int m;
        int n;
        int k;

        int lda;
        int ldb;
        int ldc;

        int offA;
        int offB;
        int offC;

        final float  sAlpha = 0.75f;
        final float  sBeta  = -0.25f;
        final double dAlpha = -1.25;
        final double dBeta  = 0.5;

        float[]  sa;
        float[]  sb;
        float[]  sBaseC;
        float[]  sc;

        double[] da;
        double[] db;
        double[] dBaseC;
        double[] dc;

        ExecutorService executor;

        @Setup(Level.Trial)
        public void setupTrial() {
            parseDimensions();

            offA = 3;
            offB = 5;
            offC = 7;

            lda = m + 2;
            ldb = k + 3;
            ldc = m + 4;

            sa     = new float [offA + lda * k];
            sb     = new float [offB + ldb * n];
            sBaseC = new float [offC + ldc * n];
            sc     = new float [sBaseC.length];

            da     = new double[offA + lda * k];
            db     = new double[offB + ldb * n];
            dBaseC = new double[offC + ldc * n];
            dc     = new double[dBaseC.length];

            fill(sa);
            fill(sb);
            fill(sBaseC);
            fill(da);
            fill(db);
            fill(dBaseC);

            int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
            executor = Executors.newFixedThreadPool(threads);
        }

        @Setup(Level.Invocation)
        public void setupInvocation() {
            System.arraycopy(sBaseC, 0, sc, 0, sBaseC.length);
            System.arraycopy(dBaseC, 0, dc, 0, dBaseC.length);
        }

        @TearDown(Level.Trial)
        public void tearDown() throws InterruptedException {
            executor.shutdown();
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                executor.shutdownNow();
                if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                    throw new AssertionError("executor did not terminate");
                }
            }
        }

        private void parseDimensions() {
            String[] tokens = dimensions.split("x");
            if (tokens.length != 3) {
                throw new IllegalArgumentException(
                        "Expected dimensions in the form MxNxK but got: " + dimensions);
            }
            m = Integer.parseInt(tokens[0]);
            n = Integer.parseInt(tokens[1]);
            k = Integer.parseInt(tokens[2]);
        }
    }

    // -----------------------------------------------------------------------
    // Threshold-boundary state (no @Param; dimensions computed at runtime)
    // -----------------------------------------------------------------------

    /**
     * Benchmark state for the two work-threshold boundary cases.
     * <p>
     * Dimensions are derived at setup time from the threshold reported by
     * {@link SgemmMRxNR#getParallelWorkThreshold()} and
     * {@link DgemmMRxNR#getParallelWorkThreshold()}, so the same benchmark
     * method correctly probes the boundary under both the Java-8 scalar kernel
     * (threshold = 45M) and the Java-25 Vector-API kernel (threshold = 25M).
     * <p>
     * Fixed: {@code m = 64}, {@code n = 7681} (3 NC-blocks).
     * <pre>
     *   kLo = floor((threshold - 1) / (m * n))  ->  work &lt; threshold (sequential)
     *   kHi = kLo + 1                            ->  work &gt; threshold (parallel)
     * </pre>
     */
    @State(Scope.Thread)
    public static class ThresholdBoundaryState {

        // fixed row / column counts
        final int m = 64;
        final int n = 7681;

        final float  sAlpha = 0.75f;
        final float  sBeta  = -0.25f;
        final double dAlpha = -1.25;
        final double dBeta  = 0.5;

        final int offA = 3;
        final int offB = 5;
        final int offC = 7;

        // sgemm boundary dimensions
        int sKlo;
        int sKhi;
        int sLdaLo, sLdaHi;
        int sLdbLo, sLdbHi;
        int sLdcLo, sLdcHi;
        float[] saLo, saHi;
        float[] sbLo, sbHi;
        float[] sBaseCLo, sBaseCHi;
        float[] scLo, scHi;

        // dgemm boundary dimensions
        int dKlo;
        int dKhi;
        int dLdaLo, dLdaHi;
        int dLdbLo, dLdbHi;
        int dLdcLo, dLdcHi;
        double[] daLo, daHi;
        double[] dbLo, dbHi;
        double[] dBaseCLo, dBaseCHi;
        double[] dcLo, dcHi;

        ExecutorService executor;

        @Setup(Level.Trial)
        public void setupTrial() {
            long sThreshold = SgemmMRxNR.getParallelWorkThreshold();
            sKlo = (int) ((sThreshold - 1) / ((long) m * n));
            sKhi = sKlo + 1;

            long dThreshold = DgemmMRxNR.getParallelWorkThreshold();
            dKlo = (int) ((dThreshold - 1) / ((long) m * n));
            dKhi = dKlo + 1;

            sLdaLo = m + 2;  sLdbLo = sKlo + 3;  sLdcLo = m + 4;
            sLdaHi = m + 2;  sLdbHi = sKhi + 3;  sLdcHi = m + 4;
            dLdaLo = m + 2;  dLdbLo = dKlo + 3;  dLdcLo = m + 4;
            dLdaHi = m + 2;  dLdbHi = dKhi + 3;  dLdcHi = m + 4;

            saLo = new float [offA + sLdaLo * sKlo];
            sbLo = new float [offB + sLdbLo * n];
            sBaseCLo = new float [offC + sLdcLo * n];
            scLo = new float [sBaseCLo.length];

            saHi = new float [offA + sLdaHi * sKhi];
            sbHi = new float [offB + sLdbHi * n];
            sBaseCHi = new float [offC + sLdcHi * n];
            scHi = new float [sBaseCHi.length];

            daLo = new double[offA + dLdaLo * dKlo];
            dbLo = new double[offB + dLdbLo * n];
            dBaseCLo = new double[offC + dLdcLo * n];
            dcLo = new double[dBaseCLo.length];

            daHi = new double[offA + dLdaHi * dKhi];
            dbHi = new double[offB + dLdbHi * n];
            dBaseCHi = new double[offC + dLdcHi * n];
            dcHi = new double[dBaseCHi.length];

            fill(saLo);  fill(saHi);
            fill(sbLo);  fill(sbHi);
            fill(sBaseCLo);  fill(sBaseCHi);
            fill(daLo);  fill(daHi);
            fill(dbLo);  fill(dbHi);
            fill(dBaseCLo);  fill(dBaseCHi);

            int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
            executor = Executors.newFixedThreadPool(threads);
        }

        @Setup(Level.Invocation)
        public void setupInvocation() {
            System.arraycopy(sBaseCLo, 0, scLo, 0, sBaseCLo.length);
            System.arraycopy(sBaseCHi, 0, scHi, 0, sBaseCHi.length);
            System.arraycopy(dBaseCLo, 0, dcLo, 0, dBaseCLo.length);
            System.arraycopy(dBaseCHi, 0, dcHi, 0, dBaseCHi.length);
        }

        @TearDown(Level.Trial)
        public void tearDown() throws InterruptedException {
            executor.shutdown();
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                executor.shutdownNow();
                if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                    throw new AssertionError("executor did not terminate");
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Reflection helpers for accessing the private parallelizeOverN method
    // -----------------------------------------------------------------------

    /**
     * Looks up the private static {@code parallelizeOverN} method of a
     * BLIS-style GEMM class. The method signature is:
     * <pre>
     * void parallelizeOverN(ExecutorService, int, int, int,
     *                       T, int, T[], int, int,
     *                       int, T[], int, int,
     *                       T, int, T[], int, int)
     * </pre>
     * where {@code T} is either {@code float} or {@code double}.
     */
    private static Method lookupParallelizeOverN(Class<?> owner,
            Class<?> scalarType, Class<?> arrayType) {
        try {
            Method method = owner.getDeclaredMethod("parallelizeOverN",
                    ExecutorService.class,
                    int.class, int.class, int.class,
                    scalarType, int.class, arrayType, int.class, int.class,
                    int.class, arrayType, int.class, int.class,
                    scalarType, int.class, arrayType, int.class, int.class);
            method.setAccessible(true);
            return method;
        } catch (ReflectiveOperationException e) {
            throw new ExceptionInInitializerError(e);
        }
    }

    private static void invokeParallelizeOverN(Method method, Object... args) {
        try {
            method.invoke(null, args);
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        } catch (InvocationTargetException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new RuntimeException(cause);
        }
    }

    // -----------------------------------------------------------------------
    // Blackhole consumption helpers
    // -----------------------------------------------------------------------

    private static void consume(float[] values, Blackhole bh) {
        bh.consume(values[0]);
        bh.consume(values[values.length >>> 1]);
        bh.consume(values[values.length - 1]);
    }

    private static void consume(double[] values, Blackhole bh) {
        bh.consume(values[0]);
        bh.consume(values[values.length >>> 1]);
        bh.consume(values[values.length - 1]);
    }

    // -----------------------------------------------------------------------
    // Array initializers
    // -----------------------------------------------------------------------

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
}
