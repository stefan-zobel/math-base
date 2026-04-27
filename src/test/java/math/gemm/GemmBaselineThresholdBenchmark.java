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
 * Benchmarks the sequential baseline GEMM kernels against the forced-parallel
 * N-splitting path to validate the parallelization heuristic.
 * <p>
 * <strong>Key finding:</strong> parallelism only pays off when at least
 * <em>3 tasks</em> can be created. With exactly 2 tasks the thread-dispatch
 * overhead consistently exceeds the speedup regardless of the total work
 * volume. The {@code shouldParallelize} guard therefore requires
 * {@code taskCount >= 3}, which implies:
 * <ul>
 * <li>{@code dgemm} (NC=1024): parallelizes only for N &gt; 2&times;NC = 2048</li>
 * <li>{@code sgemm} (NC=2048): parallelizes only for N &gt; 2&times;NC = 4096</li>
 * </ul>
 *
 * <p><strong>Transpose cases:</strong> only {@code NN} (best case, sequential
 * memory access) and {@code TT} (worst case, strided access on both operands)
 * are included. {@code TN} and {@code NT} consistently fall between these two
 * extremes and add no additional insight for heuristic validation.
 *
 * <p><strong>Dimension strings</strong> exercise the task-count boundaries for
 * both kernels and probe the {@code sgemm} work threshold at N=4097:
 * <ol>
 * <li>{@code 128x1025x128} (~16.8M ops) &mdash; dgemm 2-task <em>negative</em>
 *     control: {@code Auto} must equal {@code Sequential} after the fix.</li>
 * <li>{@code 64x2049x16} (~2.1M ops) &mdash; dgemm 3-task low-work boundary
 *     (N just above 2&times;NC for dgemm); sgemm 2-task negative control.</li>
 * <li>{@code 128x2049x64} (~16.8M ops) &mdash; dgemm 3-task moderate work;
 *     sgemm 2-task negative control.</li>
 * <li>{@code 64x4097x85} (~22.3M ops) &mdash; sgemm work-threshold lower
 *     boundary (N=4097 gives 3 sgemm tasks, but work &lt; 22.5M threshold):
 *     {@code Auto} must equal {@code Sequential} for sgemm.</li>
 * <li>{@code 64x4097x86} (~22.5M ops) &mdash; sgemm work-threshold upper
 *     boundary (3 sgemm tasks, work just above 22.5M threshold):
 *     {@code Auto} must beat {@code Sequential} for sgemm.</li>
 * <li>{@code 128x4097x128} (~67.1M ops) &mdash; sgemm 3-task large case
 *     (N just above 2&times;NC for sgemm); also 4 dgemm tasks.</li>
 * <li>{@code 128x6145x64} (~50.4M ops) &mdash; large-scale validation
 *     (4 dgemm tasks, 3 sgemm tasks).</li>
 * </ol>
 *
 * <p><strong>JMH settings rationale:</strong> measurement time is set to 2 s
 * (instead of 1 s) to obtain enough samples for the slower large-matrix cases
 * (&gt;10 ms/op) and thereby reduce error variance. The number of transpose
 * cases is halved (NN + TT only) to keep the total wall-clock time comparable
 * to earlier runs (~35 min).
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Threads(1)
@Fork(1)
@Warmup(iterations = 4, time = 2)
@Measurement(iterations = 5, time = 2)
public class GemmBaselineThresholdBenchmark {

    private static final Method SGEMM_PARALLELIZE_OVER_N = lookupParallelizeOverN(SgemmBaseline.class, float.class, float[].class);
    private static final Method DGEMM_PARALLELIZE_OVER_N = lookupParallelizeOverN(DgemmBaseline.class, double.class, double[].class);

    @Benchmark
    public void sgemmSequential(BenchmarkState state, Blackhole bh) {
        SgemmBaseline.sgemm(state.notATransposed, state.notBTransposed, state.m, state.n, state.k, state.sAlpha, state.sa,
                state.aOffset, state.lda, state.sb, state.bOffset, state.ldb, state.sBeta, state.sc, state.cOffset, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    public void sgemmForcedParallel(BenchmarkState state, Blackhole bh) {
        invokeParallelizeOverN(SGEMM_PARALLELIZE_OVER_N, state.executor, state.notATransposed, state.notBTransposed, state.m, state.n,
                state.k, state.sAlpha, state.sa, state.aOffset, state.lda, state.sb, state.bOffset, state.ldb, state.sBeta, state.sc,
                state.cOffset, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    public void sgemmAuto(BenchmarkState state, Blackhole bh) {
        SgemmBaseline.sgemm(state.notATransposed, state.notBTransposed, state.m, state.n, state.k, state.sAlpha, state.sa,
                state.aOffset, state.lda, state.sb, state.bOffset, state.ldb, state.sBeta, state.sc, state.cOffset, state.ldc,
                state.executor);
        consume(state.sc, bh);
    }

    @Benchmark
    @Measurement(iterations = 15, time = 2)
    public void dgemmSequential(BenchmarkState state, Blackhole bh) {
        DgemmBaseline.dgemm(state.notATransposed, state.notBTransposed, state.m, state.n, state.k, state.dAlpha, state.da,
                state.aOffset, state.lda, state.db, state.bOffset, state.ldb, state.dBeta, state.dc, state.cOffset, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmForcedParallel(BenchmarkState state, Blackhole bh) {
        invokeParallelizeOverN(DGEMM_PARALLELIZE_OVER_N, state.executor, state.notATransposed, state.notBTransposed, state.m, state.n,
                state.k, state.dAlpha, state.da, state.aOffset, state.lda, state.db, state.bOffset, state.ldb, state.dBeta, state.dc,
                state.cOffset, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmAuto(BenchmarkState state, Blackhole bh) {
        DgemmBaseline.dgemm(state.notATransposed, state.notBTransposed, state.m, state.n, state.k, state.dAlpha, state.da,
                state.aOffset, state.lda, state.db, state.bOffset, state.ldb, state.dBeta, state.dc, state.cOffset, state.ldc,
                state.executor);
        consume(state.dc, bh);
    }

    @State(Scope.Thread)
    public static class BenchmarkState {

        @Param({"NN", "TT"})
        public String transposeCase;

        @Param({
                // dgemm 2-task negative control (N=1025 -> 2 dgemm tasks):
                // taskCount < 3, so Auto must equal Sequential after the fix
                "128x1025x128",  // ~16.8M ops, 2 dgemm tasks -> sequential
                // dgemm 3-task boundary (N=2049 -> 3 dgemm tasks, just above 2*NC_dgemm=2048):
                // sgemm has only 2 tasks here -> sgemm negative control as well
                "64x2049x16",    //  ~2.1M ops, 3 dgemm tasks (low-work boundary)
                "128x2049x64",   // ~16.8M ops, 3 dgemm tasks / 2 sgemm tasks
                // sgemm work-threshold boundary (N=4097 -> 3 sgemm tasks, just above 2*NC_sgemm=4096):
                // probe the 22.5M work threshold; also 4 dgemm tasks for both cases
                "64x4097x85",    // ~22.3M ops - just below sgemm threshold -> Auto == Sequential for sgemm
                "64x4097x86",    // ~22.5M ops - just above sgemm threshold -> Auto beats Sequential for sgemm
                // large-scale validation
                "128x4097x128",  // ~67.1M ops, 3 sgemm tasks / 4 dgemm tasks
                "128x6145x64"    // ~50.4M ops, 3 sgemm tasks / 4 dgemm tasks
        })
        public String dimensions;

        boolean notATransposed;
        boolean notBTransposed;

        int m;
        int n;
        int k;

        int lda;
        int ldb;
        int ldc;

        int aOffset;
        int bOffset;
        int cOffset;

        final float sAlpha = 0.75f;
        final float sBeta = -0.25f;
        final double dAlpha = -1.25;
        final double dBeta = 0.5;

        float[] sa;
        float[] sb;
        float[] sBaseC;
        float[] sc;

        double[] da;
        double[] db;
        double[] dBaseC;
        double[] dc;

        ExecutorService executor;

        @Setup(Level.Trial)
        public void setupTrial() {
            parseTransposeCase();
            parseDimensions();

            aOffset = 3;
            bOffset = 5;
            cOffset = 7;

            lda = (notATransposed ? m : k) + 2;
            ldb = (notBTransposed ? k : n) + 3;
            ldc = m + 4;

            int aColumns = notATransposed ? k : m;
            int bColumns = notBTransposed ? n : k;

            sa = new float[aOffset + lda * aColumns];
            sb = new float[bOffset + ldb * bColumns];
            sBaseC = new float[cOffset + ldc * n];
            sc = new float[sBaseC.length];

            da = new double[aOffset + lda * aColumns];
            db = new double[bOffset + ldb * bColumns];
            dBaseC = new double[cOffset + ldc * n];
            dc = new double[dBaseC.length];

            fill(sa);
            fill(sb);
            fill(sBaseC);
            fill(da);
            fill(db);
            fill(dBaseC);

            executor = Executors.newFixedThreadPool(GemmParallelSupport.MAX_PARALLEL_TASKS);
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

        private void parseTransposeCase() {
            if ("NN".equals(transposeCase)) {
                notATransposed = true;
                notBTransposed = true;
            } else if ("TN".equals(transposeCase)) {
                notATransposed = false;
                notBTransposed = true;
            } else if ("NT".equals(transposeCase)) {
                notATransposed = true;
                notBTransposed = false;
            } else if ("TT".equals(transposeCase)) {
                notATransposed = false;
                notBTransposed = false;
            } else {
                throw new IllegalArgumentException("Unsupported transposeCase: " + transposeCase);
            }
        }

        private void parseDimensions() {
            String[] tokens = dimensions.split("x");
            if (tokens.length != 3) {
                throw new IllegalArgumentException("Expected dimensions in the form MxNxK but got: " + dimensions);
            }
            m = Integer.parseInt(tokens[0]);
            n = Integer.parseInt(tokens[1]);
            k = Integer.parseInt(tokens[2]);
        }
    }

    private static Method lookupParallelizeOverN(Class<?> owner, Class<?> scalarType, Class<?> arrayType) {
        try {
            Method method = owner.getDeclaredMethod("parallelizeOverN", ExecutorService.class, boolean.class, boolean.class, int.class,
                    int.class, int.class, scalarType, arrayType, int.class, int.class, arrayType, int.class, int.class, scalarType,
                    arrayType, int.class, int.class);
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
