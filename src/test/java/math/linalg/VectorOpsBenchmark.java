package math.linalg;

import org.openjdk.jmh.annotations.*;

import java.util.Random;
import java.util.concurrent.TimeUnit;

@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@State(Scope.Thread)
@Fork(1)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 7, time = 1)
public class VectorOpsBenchmark {

    @Param({"1024", "65536", "1048576"}) // 1K ... 1M
    private int size;

    private double[] baseA;
    private double[] baseB;
    private double[] baseC;
    private double[] baseD;

    private double[] arrayA;
    private double[] arrayB;
    private double[] arrayC;
    private double[] arrayD;

    private static void checkVersion() {
        String jvmVersion = System.getProperty("java.vm.version");
        String location = VectorOps.class.getClassLoader().getResource(VectorOps.class.getName().replace('.', '/') + ".class").toString();
        System.err.println("\nJVM version: " + jvmVersion);
        System.err.println("VectorOps location: " + location);
        System.err.flush();
    }

    @Setup(Level.Trial)
    public void setupTrial() {
        checkVersion();

        baseA = new double[size];
        baseB = new double[size];
        baseC = new double[size];
        baseD = new double[size];

        arrayA = new double[size];
        arrayB = new double[size];
        arrayC = new double[size];
        arrayD = new double[size];

        Random rand = new Random(123456789L);
        for (int i = 0; i < size; i++) {
            baseA[i] = rand.nextDouble();
            baseB[i] = rand.nextDouble();
            baseC[i] = rand.nextDouble();
            baseD[i] = rand.nextDouble();

            if (i % 20 == 0) {
                baseC[i] = Double.POSITIVE_INFINITY;
                baseD[i] = Double.NEGATIVE_INFINITY;
            }
        }
    }

    @Setup(Level.Invocation)
    public void setupInvocation() {
        System.arraycopy(baseA, 0, arrayA, 0, size);
        System.arraycopy(baseB, 0, arrayB, 0, size);
        System.arraycopy(baseC, 0, arrayC, 0, size);
        System.arraycopy(baseD, 0, arrayD, 0, size);
    }

    @Benchmark
    public double benchmarkDotProduct() {
        return VectorOps.dotProduct(arrayA, arrayB);
    }

    @Benchmark
    public void benchmarkTimesEquals() {
        VectorOps.timesEquals(arrayA, 1.0001);
    }

    @Benchmark
    public void benchmarkPlusEquals() {
        VectorOps.plusEquals(arrayC, arrayD);
    }

    @Benchmark
    public void benchmarkPlusEqualsWithFactor() {
        VectorOps.plusEquals(arrayC, arrayD, 2.0);
    }
}
