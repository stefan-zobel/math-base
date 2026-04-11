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
    private double[] arrayC;
    private double[] arrayD;

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
    
    @Setup
    public void setupPlusEquals() {
        arrayC = new double[size];
        arrayD = new double[size];
        java.util.Random rand = new java.util.Random();
        
        for (int i = 0; i < size; i++) {
            arrayC[i] = rand.nextDouble();
            arrayD[i] = rand.nextDouble();
            
            if (i % 20 == 0) {
                arrayC[i] = Double.POSITIVE_INFINITY;
                arrayD[i] = Double.NEGATIVE_INFINITY;
            }
        }
    }

    /**
     * Benchmark für plusEquals(double[] m1, double[] m2)
     */
    @Benchmark
    public void benchmarkPlusEquals() {
        VectorOps.plusEquals(arrayC, arrayD);
    }

    /**
     * Benchmark für plusEquals(double[] m1, double[] m2, double factor)
     */
    @Benchmark
    public void benchmarkPlusEqualsWithFactor() {
        VectorOps.plusEquals(arrayC, arrayD, 2.0);
    }
}
