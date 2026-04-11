package math.linalg;

import org.openjdk.jmh.annotations.*;
import java.util.concurrent.TimeUnit;
import java.util.Random;

@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@State(Scope.Thread)
@Fork(1)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 7, time = 1)
public class VectorOpsBenchmark {

    @Param({"1024", "65536", "1048576"}) // 1K ... 1M
    private int size;

    private double[] arrayA;
    private double[] arrayB;

    @Setup
    public void setup() {
        arrayA = new double[size];
        arrayB = new double[size];
        Random rand = new Random();
        for (int i = 0; i < size; i++) {
            arrayA[i] = rand.nextDouble();
            arrayB[i] = rand.nextDouble();
        }
    }

    @Benchmark
    public double benchmarkDotProduct() {
        return VectorOps.dotProduct(arrayA, arrayB);
    }

    @Benchmark
    public void benchmarkTimesEquals() {
        VectorOps.timesEquals(arrayA, 1.0001);
    }
}
