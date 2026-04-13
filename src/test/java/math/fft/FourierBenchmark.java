package math.fft;

import org.openjdk.jmh.annotations.*;
import java.util.concurrent.TimeUnit;
import java.util.Random;

@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@State(Scope.Thread)
@Fork(1)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 7, time = 1)
public class FourierBenchmark {

    @Param({"256", "1000", "1024", "4000", "4096"})
    private int size;

    private double[] arrayReal;
    private double[] arrayImag;

    @Setup
    public void setup() {
        arrayReal = new double[size];
        arrayImag = new double[size];
        Random rand = new Random(42);
        for (int i = 0; i < size; i++) {
            arrayReal[i] = rand.nextDouble();
            arrayImag[i] = rand.nextDouble();
        }
    }

    @Benchmark
    public Object benchmarkForwardReal() {
        return Fourier.forwardDFT(arrayReal);
    }

    @Benchmark
    public Object benchmarkForwardComplex() {
        return Fourier.forwardDFT(arrayReal, arrayImag);
    }
}
