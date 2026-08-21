package math.fft;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Assert;
import org.junit.Test;

public final class FourierTest {

    private static final double TOL = 1.0e-9;

    @Test
    public void testIsVectorized() {
        String versionStr = System.getProperty("java.specification.version");
        double version = Double.parseDouble(versionStr.startsWith("1.") ? versionStr.substring(2) : versionStr);
        if (version >= 25.0) {
            assertTrue("Java 25 should use the Vector API", Fourier.isVectorized());
        } else {
            assertFalse("Java 8 must not use the Vector API", Fourier.isVectorized());
        }
    }

    @Test
    public void testForwardRealMatchesNaiveForPowerOfTwoAndBluesteinLengths() {
        assertForwardRealMatchesNaive(sampleReal(4), 1.0e-12);
        assertForwardRealMatchesNaive(sampleReal(8), 1.0e-12);
        assertForwardRealMatchesNaive(sampleReal(7), TOL);
        assertForwardRealMatchesNaive(sampleReal(13), TOL);
    }

    @Test
    public void testForwardComplexMatchesNaiveForMixedLengths() {
        assertForwardComplexMatchesNaive(sampleReal(4), sampleImag(4), 1.0e-12);
        assertForwardComplexMatchesNaive(sampleReal(8), sampleImag(8), 1.0e-12);
        assertForwardComplexMatchesNaive(sampleReal(7), sampleImag(7), TOL);
        assertForwardComplexMatchesNaive(sampleReal(11), sampleImag(11), TOL);
    }

    @Test
    public void testVectorMaskBoundaries() {
        // Test array lengths that trigger VectorMask tail loops (e.g. not a multiple of 4/8)
        assertForwardRealMatchesNaive(sampleReal(14), TOL);
        assertForwardRealMatchesNaive(sampleReal(31), TOL);
        assertForwardRealMatchesNaive(sampleReal(63), TOL);
        
        assertForwardComplexMatchesNaive(sampleReal(14), sampleImag(14), TOL);
        assertForwardComplexMatchesNaive(sampleReal(31), sampleImag(31), TOL);
        assertForwardComplexMatchesNaive(sampleReal(63), sampleImag(63), TOL);
        
        assertRoundTrip(sampleReal(14), sampleImag(14), TOL);
        assertRoundTrip(sampleReal(31), sampleImag(31), TOL);
        assertRoundTrip(sampleReal(63), sampleImag(63), TOL);
    }

    @Test
    public void testLargePowerOfTwo() {
        // Ensures that the unrolled / pure Radix-4 paths are exhaustively tested
        assertForwardRealMatchesNaive(sampleReal(256), 1.0e-10);
        assertForwardComplexMatchesNaive(sampleReal(256), sampleImag(256), 1.0e-10);
        assertRoundTrip(sampleReal(256), sampleImag(256), 1.0e-10);
    }

    @Test
    public void testLargeNonPowerOfTwo() {
        // Exhaustively tests the Bluestein paths with length >= 128
        assertForwardRealMatchesNaive(sampleReal(137), TOL);
        assertForwardComplexMatchesNaive(sampleReal(137), sampleImag(137), TOL);
        assertRoundTrip(sampleReal(137), sampleImag(137), TOL);
    }

    @Test
    public void testVectorizedDotAndElementwiseProduct() {
        ComplexArray a = new ComplexArray(sampleReal(42), sampleImag(42));
        ComplexArray b = new ComplexArray(sampleImag(42), sampleReal(42)); // cross inputs 
        
        // Elementwise
        ComplexArray actProd = ComplexArray.elementwiseProduct(a, b);
        ComplexArray expProd = naiveElementwiseProduct(a, b);
        Assert.assertArrayEquals(expProd.re(), actProd.re(), TOL);
        Assert.assertArrayEquals(expProd.im(), actProd.im(), TOL);

        // Dot
        double[] actDot = ComplexArray.dot(a, b);
        double[] expDot = naiveDot(a, b);
        Assert.assertArrayEquals(expDot, actDot, TOL);
    }

    @Test
    public void testAbsSquared() {
        ComplexArray ca = new ComplexArray(sampleReal(31), sampleImag(31)); // Using 31 to trigger tail loop
        
        // Unscaled
        double[] actualUnscaled = ca.absSquared();
        double[] expectedUnscaled = naiveAbsSquared(ca, false);
        Assert.assertArrayEquals(expectedUnscaled, actualUnscaled, TOL);

        // Scaled
        double[] actualScaled = ca.absSquaredScaled();
        double[] expectedScaled = naiveAbsSquared(ca, true);
        Assert.assertArrayEquals(expectedScaled, actualScaled, TOL);
    }

    private static double[] naiveAbsSquared(ComplexArray ca, boolean scaled) {
        int n = ca.length();
        double[] res = new double[n];
        double invScale = scaled ? (1.0 / n) : 1.0;
        double max = 0.0;
        for (int i = 0; i < n; ++i) {
            double rX = ca.re()[i];
            double iY = ca.im()[i];
            double square = (rX * rX + iY * iY) * invScale;
            res[i] = square;
            max = Math.max(max, square);
        }
        // the same relative rule the class uses, squared because these are
        // squares
        double factor = 4.0 * n * math.MathConsts.MACH_EPS_DBL;
        double cutoff = max * factor * factor;
        for (int i = 0; i < n; ++i) {
            if (res[i] <= cutoff) {
                res[i] = 0.0;
            }
        }
        return res;
    }

    private static double[] naiveDot(ComplexArray a, ComplexArray b) {
        double res_re = 0.0;
        double res_im = 0.0;
        for (int i = 0; i < a.length(); ++i) {
            double a_re = a.re()[i];
            double b_re = b.re()[i];
            double a_im = a.im()[i];
            double b_im = b.im()[i];
            res_re += a_re * b_re - a_im * b_im;
            res_im += a_re * b_im + a_im * b_re;
        }
        return new double[] { res_re, res_im };
    }

    private static ComplexArray naiveElementwiseProduct(ComplexArray a, ComplexArray b) {
        int n = a.length();
        double[] real = new double[n];
        double[] imag = new double[n];
        for (int i = 0; i < n; ++i) {
            double a_re = a.re()[i];
            double b_re = b.re()[i];
            double a_im = a.im()[i];
            double b_im = b.im()[i];
            real[i] = a_re * b_re - a_im * b_im;
            imag[i] = a_re * b_im + a_im * b_re;
        }
        return new ComplexArray(real, imag, false);
    }

    @Test
    public void testInverseRoundTripForMixedLengths() {
        assertRoundTrip(sampleReal(4), sampleImag(4), 1.0e-12);
        assertRoundTrip(sampleReal(16), sampleImag(16), 1.0e-12);
        assertRoundTrip(sampleReal(7), sampleImag(7), TOL);
        assertRoundTrip(sampleReal(15), sampleImag(15), TOL);
    }

    private static void assertForwardRealMatchesNaive(double[] data, double tol) {
        ComplexArray actual = Fourier.forwardDFT(data);
        ComplexArray expected = ComplexArray.naiveForwarDFT(data);
        Assert.assertArrayEquals(expected.re(), actual.re(), tol);
        Assert.assertArrayEquals(expected.im(), actual.im(), tol);
    }

    private static void assertForwardComplexMatchesNaive(double[] real, double[] imag, double tol) {
        ComplexArray actual = Fourier.forwardDFT(real, imag);
        ComplexArray expected = new ComplexArray(real, imag, true).naiveForwardDFT();
        Assert.assertArrayEquals(expected.re(), actual.re(), tol);
        Assert.assertArrayEquals(expected.im(), actual.im(), tol);
    }

    private static void assertRoundTrip(double[] real, double[] imag, double tol) {
        ComplexArray transformed = Fourier.forwardDFT(real, imag);
        ComplexArray restored = Fourier.inverseDFT(transformed);
        Assert.assertArrayEquals(real, restored.re(), tol);
        Assert.assertArrayEquals(imag, restored.im(), tol);
    }

    // ---------------------------------------------------------------------
    // The tests above compare the shipped code against reimplementations that
    // live in this file, which can only catch a wrong transcription of a rule,
    // never a wrong rule. The ones below are properties the transform has to
    // satisfy whatever its implementation, and they are what the cleanup
    // threshold is actually held to.
    // ---------------------------------------------------------------------

    /** A constant series: only the DC bin is nonzero. */
    private static double[] constant(int n) {
        double[] x = new double[n];
        for (int i = 0; i < n; ++i) {
            x[i] = 2.5;
        }
        return x;
    }

    /** Real and even, so the whole imaginary part of the transform is zero. */
    private static double[] realEven(int n) {
        double[] x = new double[n];
        for (int i = 0; i < n; ++i) {
            int j = Math.min(i, n - i);
            x[i] = 1.0 + Math.cos(0.9 * j) + 0.3 * j;
        }
        return x;
    }

    /** A cosine exactly on a bin: two nonzero coefficients, the rest zero. */
    private static double[] binSinusoid(int n) {
        double[] x = new double[n];
        int k0 = Math.max(1, n / 8);
        for (int i = 0; i < n; ++i) {
            x[i] = Math.cos(2.0 * Math.PI * k0 * i / n);
        }
        return x;
    }

    private static double[] scaled(double[] x, double c) {
        double[] y = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            y[i] = c * x[i];
        }
        return y;
    }

    private static int exactZeros(ComplexArray f) {
        int count = 0;
        for (int i = 0; i < f.length(); ++i) {
            if (f.re()[i] == 0.0) {
                ++count;
            }
            if (f.im()[i] == 0.0) {
                ++count;
            }
        }
        return count;
    }

    @Test
    public void testCoefficientsThatAreExactlyZeroComeBackAsZero() {
        // an absolute threshold delivered 3 of 1024 of these at unit scale and
        // all 1024 of them plus the signal at 1e-8; a relative one has to
        // deliver exactly the right set at every scale
        int[] lengths = { 64, 256, 1024, 500, 1000 };
        double[] scales = { 1.0, 1.0e-8, 1.0e-14, 1.0e8 };
        for (int k = 0; k < lengths.length; ++k) {
            int n = lengths[k];
            for (int s = 0; s < scales.length; ++s) {
                double c = scales[s];
                String at = " at n=" + n + ", scale " + c;
                Assert.assertEquals("constant" + at, 2 * n - 1,
                        exactZeros(Fourier.forwardDFT(scaled(constant(n), c))));
                Assert.assertEquals("real even" + at, n,
                        exactZeros(Fourier.forwardDFT(scaled(realEven(n), c))));
                Assert.assertEquals("bin sinusoid" + at, 2 * n - 2,
                        exactZeros(Fourier.forwardDFT(scaled(binSinusoid(n), c))));
            }
        }
    }

    @Test
    public void testTheResultDoesNotDependOnTheUnitTheDataIsMeasuredIn() {
        int n = 512;
        double[] base = sampleReal(n);
        ComplexArray reference = Fourier.forwardDFT(base);
        double[] refSquares = reference.absSquared();
        // round-off in a transform is proportional to the largest coefficient,
        // not to each coefficient separately -- which is the whole point of
        // this change, so the tolerance here is stated the same way
        double magnitude = 0.0;
        double largestSquare = 0.0;
        for (int i = 0; i < n; ++i) {
            magnitude = Math.max(magnitude, Math.max(Math.abs(reference.re()[i]), Math.abs(reference.im()[i])));
            largestSquare = Math.max(largestSquare, refSquares[i]);
        }
        for (int e = -12; e <= 12; e += 4) {
            if (e == 0) {
                continue;
            }
            double c = Math.pow(10.0, e);
            ComplexArray f = Fourier.forwardDFT(scaled(base, c));
            double[] squares = f.absSquared();
            for (int i = 0; i < n; ++i) {
                Assert.assertEquals("re[" + i + "] at 1e" + e, reference.re()[i], f.re()[i] / c,
                        1.0e-12 * magnitude);
                Assert.assertEquals("im[" + i + "] at 1e" + e, reference.im()[i], f.im()[i] / c,
                        1.0e-12 * magnitude);
                Assert.assertEquals("|X|^2[" + i + "] at 1e" + e, refSquares[i], squares[i] / (c * c),
                        1.0e-12 * largestSquare);
            }
        }
    }

    @Test
    public void testTheRoundTripSurvivesAtEveryScale() {
        // this failed from 1e-12 downwards while the threshold was absolute,
        // and returned all zeros from 1e-16
        for (int e = 0; e <= 20; e += 4) {
            double c = Math.pow(10.0, -e);
            assertScaledRoundTrip(1024, c);
            assertScaledRoundTrip(1000, c);
        }
    }

    private static void assertScaledRoundTrip(int n, double c) {
        double[] real = scaled(sampleReal(n), c);
        double[] imag = scaled(sampleImag(n), c);
        ComplexArray restored = Fourier.inverseDFT(Fourier.forwardDFT(real, imag));
        double magnitude = 0.0;
        for (int i = 0; i < n; ++i) {
            magnitude = Math.max(magnitude, Math.max(Math.abs(real[i]), Math.abs(imag[i])));
        }
        for (int i = 0; i < n; ++i) {
            Assert.assertEquals("n=" + n + ", scale " + c + ", re[" + i + "]", real[i], restored.re()[i],
                    1.0e-10 * magnitude);
            Assert.assertEquals("n=" + n + ", scale " + c + ", im[" + i + "]", imag[i], restored.im()[i],
                    1.0e-10 * magnitude);
        }
    }

    @Test
    public void testAGenuineWeakLineIsNotMistakenForRoundOff() {
        // the price of clearing the round-off floor is that a coefficient
        // below the floor cannot be told from it; this pins where that is
        int n = 1024;
        for (int e = 6; e <= 12; e += 2) {
            double ratio = Math.pow(10.0, -e);
            double[] x = new double[n];
            for (int i = 0; i < n; ++i) {
                x[i] = Math.cos(2.0 * Math.PI * 100 * i / n) + ratio * Math.cos(2.0 * Math.PI * 300 * i / n);
            }
            ComplexArray f = Fourier.forwardDFT(x);
            Assert.assertTrue("a line at 1e-" + e + " of the strongest was discarded", f.re()[300] != 0.0);
        }
    }

    @Test
    public void testDegenerateInputsAreNotTurnedIntoZeros() {
        // the reference for a relative threshold is the largest magnitude; if
        // that is zero or not finite there is nothing to be relative to, and
        // an unguarded cutoff would zero the whole array
        ComplexArray allZero = Fourier.forwardDFT(new double[8]);
        for (int i = 0; i < allZero.length(); ++i) {
            Assert.assertEquals(0.0, allZero.re()[i], 0.0);
            Assert.assertEquals(0.0, allZero.im()[i], 0.0);
        }
        ComplexArray single = Fourier.forwardDFT(new double[] { 3.0 });
        Assert.assertEquals(3.0, single.re()[0], 0.0);

        ComplexArray withInfinity = Fourier.forwardDFT(new double[] { 1.0, 2.0, Double.POSITIVE_INFINITY, 4.0 });
        boolean anyInfinite = false;
        for (int i = 0; i < withInfinity.length(); ++i) {
            anyInfinite |= Double.isInfinite(withInfinity.re()[i]) || Double.isInfinite(withInfinity.im()[i]);
        }
        Assert.assertTrue("an infinite input was silently turned into zeros", anyInfinite);

        ComplexArray withNaN = Fourier.forwardDFT(new double[] { 1.0, 2.0, Double.NaN, 4.0 });
        boolean anyNaN = false;
        for (int i = 0; i < withNaN.length(); ++i) {
            anyNaN |= Double.isNaN(withNaN.re()[i]) || Double.isNaN(withNaN.im()[i]);
        }
        Assert.assertTrue("a NaN input was silently turned into zeros", anyNaN);
    }

    @Test
    public void testDotAndElementwiseProductAreScaleInvariantToo() {
        int n = 31;
        ComplexArray a = new ComplexArray(sampleReal(n), sampleImag(n));
        ComplexArray b = new ComplexArray(sampleImag(n), sampleReal(n));
        double[] referenceDot = ComplexArray.dot(a, b);
        ComplexArray referenceProduct = ComplexArray.elementwiseProduct(a, b);
        for (int e = -10; e <= 10; e += 5) {
            if (e == 0) {
                continue;
            }
            double c = Math.pow(10.0, e);
            ComplexArray sa = new ComplexArray(scaled(sampleReal(n), c), scaled(sampleImag(n), c));
            double[] dot = ComplexArray.dot(sa, b);
            Assert.assertEquals("dot re at 1e" + e, referenceDot[0], dot[0] / c,
                    1.0e-12 * Math.abs(referenceDot[0]));
            Assert.assertEquals("dot im at 1e" + e, referenceDot[1], dot[1] / c,
                    1.0e-12 * Math.abs(referenceDot[1]));
            ComplexArray product = ComplexArray.elementwiseProduct(sa, b);
            for (int i = 0; i < n; ++i) {
                Assert.assertEquals("product re at 1e" + e, referenceProduct.re()[i], product.re()[i] / c,
                        1.0e-12 * Math.abs(referenceProduct.re()[i]) + 1.0e-300);
            }
        }
    }

    private static double[] sampleReal(int n) {
        double[] data = new double[n];
        for (int i = 0; i < n; ++i) {
            data[i] = Math.sin(0.37 * i) + Math.cos(0.19 * (i + 1)) + ((i % 3) - 1) * 0.125;
        }
        return data;
    }

    private static double[] sampleImag(int n) {
        double[] data = new double[n];
        for (int i = 0; i < n; ++i) {
            data[i] = Math.cos(0.23 * i) - Math.sin(0.41 * (i + 2)) + ((i % 4) - 1.5) * 0.0625;
        }
        return data;
    }
}
