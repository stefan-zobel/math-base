package math.gemm;

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
 * Compares the two production APIs that may be selected by a future dispatcher:
 * <ul>
 *   <li>{@link SgemmBaseline}/{@link DgemmBaseline} with their internal
 *       threshold-based parallelization (capped at 4 worker tasks)</li>
 *   <li>{@link SgemmMRxNR}/{@link DgemmMRxNR} with their internal
 *       threshold-based parallelization (uncapped task count, limited by cores)
 *   </li>
 * </ul>
 *
 * <p>
 * This benchmark currently covers the hardcoded transpose cases NN and TT.
 *
 * <p>
 * The dimension set is grouped by practical use cases:
 * <ul>
 *   <li><strong>Tiny / small</strong>: expected Baseline-favored region where
 *       packing overhead can dominate.</li>
 *   <li><strong>Medium</strong>: transition region where either implementation
 *       may win depending on JVM and hardware.</li>
 *   <li><strong>Large / NC-sensitive</strong>: expected BLIS-favored region,
 *       including NC boundary cases ({@code n=3841} and {@code n=7681}).</li>
 *   <li><strong>Shape variants</strong>: skinny-N and skinny-K cases to avoid
 *       overfitting thresholds to square matrices only.</li>
 * </ul>
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Threads(1)
@Fork(1)
@Warmup(iterations = 4, time = 2)
@Measurement(iterations = 10, time = 2)
public class GemmApiSwitchBenchmark {

    @Benchmark
    public void sgemmBaselineSequential(BenchmarkState state, Blackhole bh) {
        SgemmBaseline.sgemm(state.notATransposed, state.notBTransposed,
                state.m, state.n, state.k, state.sAlpha,
                state.sa, state.offA, state.lda,
                state.sb, state.offB, state.ldb,
                state.sBeta, state.sc, state.offC, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    public void sgemmBlisSequential(BenchmarkState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.k, state.sAlpha,
                state.offA, state.sa, state.incRowA, state.incColA,
                state.offB, state.sb, state.incRowB, state.incColB,
                state.sBeta, state.offC, state.sc, 1, state.ldc);
        consume(state.sc, bh);
    }

    @Benchmark
    public void sgemmBaselineAuto(BenchmarkState state, Blackhole bh) {
        SgemmBaseline.sgemm(state.notATransposed, state.notBTransposed,
                state.m, state.n, state.k, state.sAlpha,
                state.sa, state.offA, state.lda,
                state.sb, state.offB, state.ldb,
                state.sBeta, state.sc, state.offC, state.ldc,
                state.executor);
        consume(state.sc, bh);
    }

    @Benchmark
    public void sgemmBlisAuto(BenchmarkState state, Blackhole bh) {
        SgemmMRxNR.sgemm(state.m, state.n, state.k, state.sAlpha,
                state.offA, state.sa, state.incRowA, state.incColA,
                state.offB, state.sb, state.incRowB, state.incColB,
                state.sBeta, state.offC, state.sc, 1, state.ldc,
                state.executor);
        consume(state.sc, bh);
    }

    @Benchmark
    public void dgemmBaselineSequential(BenchmarkState state, Blackhole bh) {
        DgemmBaseline.dgemm(state.notATransposed, state.notBTransposed,
                state.m, state.n, state.k, state.dAlpha,
                state.da, state.offA, state.lda,
                state.db, state.offB, state.ldb,
                state.dBeta, state.dc, state.offC, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmBlisSequential(BenchmarkState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.k, state.dAlpha,
                state.offA, state.da, state.incRowA, state.incColA,
                state.offB, state.db, state.incRowB, state.incColB,
                state.dBeta, state.offC, state.dc, 1, state.ldc);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmBaselineAuto(BenchmarkState state, Blackhole bh) {
        DgemmBaseline.dgemm(state.notATransposed, state.notBTransposed,
                state.m, state.n, state.k, state.dAlpha,
                state.da, state.offA, state.lda,
                state.db, state.offB, state.ldb,
                state.dBeta, state.dc, state.offC, state.ldc,
                state.executor);
        consume(state.dc, bh);
    }

    @Benchmark
    public void dgemmBlisAuto(BenchmarkState state, Blackhole bh) {
        DgemmMRxNR.dgemm(state.m, state.n, state.k, state.dAlpha,
                state.offA, state.da, state.incRowA, state.incColA,
                state.offB, state.db, state.incRowB, state.incColB,
                state.dBeta, state.offC, state.dc, 1, state.ldc,
                state.executor);
        consume(state.dc, bh);
    }

    @State(Scope.Thread)
    public static class BenchmarkState {

        @Param({"NN", "TT"})
        public String transposeCase;

        /**
         * Reduced high-signal set in the form MxNxK.
         * Focuses on medium/large and shape-sensitive cases to validate BLIS behavior
         * without spending time on many tiny matrices.
         */
        @Param({
                "128x256x128",
                "256x256x64",
                "64x1024x64",
                "128x3841x64",
                "128x7681x64",
                "128x7681x128",
                "32x32x32",
                "64x64x64",
                "96x96x96",
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

        boolean notATransposed;
        boolean notBTransposed;

        int incRowA;
        int incColA;
        int incRowB;
        int incColB;

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

        private static void checkVersion() {
            String jvmVersion = System.getProperty("java.vm.version");
            String location = DgemmMRxNR.class.getClassLoader().getResource(DgemmMRxNR.class.getName().replace('.', '/') + ".class").toString();
            System.err.println("\nJVM version: " + jvmVersion);
            System.err.println("DgemmMRxNR location: " + location);
            System.err.flush();
        }

        @Setup(Level.Trial)
        public void setupTrial() {
            checkVersion();
            parseDimensions();
            parseTransposeCase();

            offA = 3;
            offB = 5;
            offC = 7;

            configureStridesAndLeadingDimensions();

            sa = new float[offA + lda * physicalColsA()];
            sb = new float[offB + ldb * physicalColsB()];
            sBaseC = new float[offC + ldc * n];
            sc = new float[sBaseC.length];

            da = new double[offA + lda * physicalColsA()];
            db = new double[offB + ldb * physicalColsB()];
            dBaseC = new double[offC + ldc * n];
            dc = new double[dBaseC.length];

            fill(sa);
            fill(sb);
            fill(sBaseC);
            fill(da);
            fill(db);
            fill(dBaseC);

            // Use machine parallelism; each implementation applies its own policy.
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
                throw new IllegalArgumentException("Expected dimensions in the form MxNxK but got: " + dimensions);
            }
            m = Integer.parseInt(tokens[0]);
            n = Integer.parseInt(tokens[1]);
            k = Integer.parseInt(tokens[2]);
        }

        private void parseTransposeCase() {
            switch (transposeCase) {
            case "NN":
                notATransposed = true;
                notBTransposed = true;
                break;
            case "TT":
                notATransposed = false;
                notBTransposed = false;
                break;
            default:
                throw new IllegalArgumentException("Unsupported transposeCase: " + transposeCase);
            }
        }

        private void configureStridesAndLeadingDimensions() {
            int physicalRowsA = notATransposed ? m : k;
            int physicalRowsB = notBTransposed ? k : n;

            lda = physicalRowsA + 2;
            ldb = physicalRowsB + 3;
            ldc = m + 4;

            if (notATransposed) {
                incRowA = 1;
                incColA = lda;
            } else {
                incRowA = lda;
                incColA = 1;
            }

            if (notBTransposed) {
                incRowB = 1;
                incColB = ldb;
            } else {
                incRowB = ldb;
                incColB = 1;
            }
        }

        private int physicalColsA() {
            return notATransposed ? k : m;
        }

        private int physicalColsB() {
            return notBTransposed ? n : k;
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
