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
        // ComplexArray.TOL
        double tol = 5.0 * 2.220446049250313E-16; // 5.0 * MathConsts.MACH_EPS_DBL

        for (int i = 0; i < n; ++i) {
            double rX = ca.re()[i];
            double iY = ca.im()[i];
            double square = (rX * rX + iY * iY) * invScale;
            if (square <= tol) {
                square = 0.0;
            }
            res[i] = square;
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
