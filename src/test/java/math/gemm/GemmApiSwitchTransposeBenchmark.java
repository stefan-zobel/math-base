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
 * Focused dispatcher benchmark for transpose-sensitive cases NT and TN.
 *
 * <p>This benchmark is intentionally smaller and faster than
 * {@link GemmApiSwitchBenchmark}: it measures AUTO paths only, so the runtime
 * impact of adding NT/TN coverage stays manageable.</p>
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Threads(1)
@Fork(1)
@Warmup(iterations = 2, time = 1)
@Measurement(iterations = 5, time = 1)
public class GemmApiSwitchTransposeBenchmark {

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

        @Param({"NT", "TN"})
        public String transposeCase;

        /**
         * Small, medium, and NC-sensitive shapes to keep runtime low while
         * preserving dispatcher signal.
         */
        @Param({
                "32x32x32",
                "96x96x96",
                "128x256x128",
                "64x1024x64",
                "128x3841x64",
                "128x7681x64"
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
            String location = DgemmMRxNR.class.getClassLoader()
                    .getResource(DgemmMRxNR.class.getName().replace('.', '/') + ".class").toString();
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
            case "NT":
                notATransposed = true;
                notBTransposed = false;
                break;
            case "TN":
                notATransposed = false;
                notBTransposed = true;
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
