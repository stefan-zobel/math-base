package math.linalg;

import static org.junit.Assert.*;
import org.junit.Test;

import math.minpack.Minpack_f77;

public class VectorOpsTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testIsVectorized() {
        String versionStr = System.getProperty("java.specification.version");
        double version = Double.parseDouble(versionStr.startsWith("1.") ? versionStr.substring(2) : versionStr);
        if (version >= 25.0) {
            assertTrue("Java 25 should use the Vector API", VectorOps.isVectorized());
        } else {
            assertFalse("Java 8 must not use the Vector API", VectorOps.isVectorized());
        }
    }

    @Test
    public void testTimesEquals() {
        double[] m = {1.0, 2.0, 3.0, 4.0, 5.0};
        VectorOps.timesEquals(m, 2.5);

        double[] expected = {2.5, 5.0, 7.5, 10.0, 12.5};
        assertArrayEquals(expected, m, EPSILON);
    }

    @Test
    public void testDotProduct() {
        double[] m1 = {1.0, 2.0, 3.0};
        double[] m2 = {4.0, 5.0, 6.0};
        // (1*4) + (2*5) + (3*6) = 4 + 10 + 18 = 32
        double result = VectorOps.dotProduct(m1, m2);
        assertEquals(32.0, result, EPSILON);
    }

    @Test
    public void testPlusEquals() {
        double[] m1 = {
            10.0,
            Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY,
            Double.POSITIVE_INFINITY
        };
        double[] m2 = {
            5.5,
            Double.NEGATIVE_INFINITY, // Inf + -Inf -> 0.0
            Double.POSITIVE_INFINITY, // -Inf + Inf -> 0.0
            Double.POSITIVE_INFINITY  // Inf + Inf -> Inf
        };

        VectorOps.plusEquals(m1, m2);

        assertEquals(15.5, m1[0], EPSILON);
        assertEquals("Inf + -Inf must be 0.0", 0.0, m1[1], EPSILON);
        assertEquals("-Inf + Inf must be 0.0", 0.0, m1[2], EPSILON);
        assertEquals("Inf + Inf must stay Inf", Double.POSITIVE_INFINITY, m1[3], 0.0);
    }

    @Test
    public void testPlusEqualsWithFactor() {
        double[] m1 = {1.0, Double.POSITIVE_INFINITY, 2.0};
        double[] m2 = {3.0, Double.NEGATIVE_INFINITY, 4.0};
        double factor = 2.0;

        // Operation: m1[i] += m2[i] * factor
        VectorOps.plusEquals(m1, m2, factor);

        assertEquals(1.0 + (3.0 * 2.0), m1[0], EPSILON); // 7.0
        assertEquals("Special case Infinity with a factor", 0.0, m1[1], EPSILON);
        assertEquals(2.0 + (4.0 * 2.0), m1[2], EPSILON); // 10.0
    }

    @Test
    public void testTailCleanup() {
        double[] m1 = {1, 1, 1, 1, 1, 1, 1};
        double[] m2 = {2, 2, 2, 2, 2, 2, 2};

        VectorOps.plusEquals(m1, m2);

        for (double val : m1) {
            assertEquals(3.0, val, EPSILON);
        }
    }

    // ------------------------------------------------------------------
    // the norms
    // ------------------------------------------------------------------

    /** Lengths that straddle every plausible vector lane count. */
    private static final int[] LENGTHS = { 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 33, 64, 129, 1000 };

    private long lcg = 987654321987654321L;

    /** Deterministic pseudo random number in [-1, 1). */
    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) / (double) (1L << 53)) * 2.0 - 1.0;
    }

    /**
     * The Euclidean norm from MINPACK, an independent implementation of the
     * same quantity by a different method: three accumulators and a split by
     * magnitude rather than one scaling pass. Note that the f2j port indexes
     * from one, so the argument has to be shifted.
     */
    private static double enorm(double[] v) {
        double[] shifted = new double[v.length + 1];
        System.arraycopy(v, 0, shifted, 1, v.length);
        return Minpack_f77.enorm_f77(v.length, shifted);
    }

    /** The norm as the straightforward loop computes it. */
    private static double naiveTwoNorm(double[] v) {
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            sum += v[i] * v[i];
        }
        return Math.sqrt(sum);
    }

    private static void assertSameBits(String message, double expected, double actual) {
        assertEquals(message + " : expected " + expected + " but was " + actual,
                Double.doubleToRawLongBits(expected), Double.doubleToRawLongBits(actual));
    }

    /**
     * Over the whole exponent range, against MINPACK. The vector holds
     * {@code +/-2^k} in every position, so the answer is {@code 2^k * sqrt(n)}
     * and no rounding enters through the data.
     */
    @Test
    public void testTheTwoNormAgreesWithMinpackOverTheWholeExponentRange() {
        int[] lengths = { 1, 2, 3, 7, 8, 17, 100 };
        for (int k = -1074; k <= 1023; ++k) {
            double s = Math.scalb(1.0, k);
            for (int li = 0; li < lengths.length; ++li) {
                int n = lengths[li];
                double[] v = new double[n];
                for (int i = 0; i < n; ++i) {
                    v[i] = ((i & 1) == 0) ? s : -s;
                }
                double expected = enorm(v);
                if (Double.isNaN(expected) || Double.isInfinite(expected) || expected == 0.0) {
                    // outside the range where MINPACK itself is the authority
                    continue;
                }
                double actual = VectorOps.twoNorm(v);
                assertEquals("n = " + n + " at 2^" + k, expected, actual, 1.0e-15 * expected);
            }
        }
    }

    /**
     * The three magnitudes the naive accumulation of squares gets wrong. It
     * returns {@code Infinity} from about {@code 8.4e152} upwards, {@code 0.0}
     * below about {@code 1.1e-162}, and is already an order of magnitude short
     * of full precision before that.
     */
    @Test
    public void testTheTwoNormNeitherOverflowsNorUnderflows() {
        double[] scales = { 1.0e154, 1.0e160, 1.0e300, Double.MAX_VALUE / 4.0, 1.0e-160, 1.0e-162, 1.0e-200,
                1.0e-300, Double.MIN_NORMAL, Double.MIN_VALUE };
        for (int k = 0; k < scales.length; ++k) {
            double s = scales[k];
            double[] v = { s, s, s, s };
            double actual = VectorOps.twoNorm(v);
            assertTrue("the norm of a finite vector must be finite, at " + s, !Double.isInfinite(actual));
            assertTrue("the norm of a non-zero vector must not be zero, at " + s, actual != 0.0);
            // the exact answer is 2 * s
            assertEquals("at " + s, 2.0 * s, actual, 1.0e-14 * 2.0 * s);
        }
    }

    /**
     * For every vector whose magnitude leaves the squares in the normal range
     * -- which is every vector this library is likely to see -- the result is
     * exactly the one the straightforward loop gives, so nothing that was
     * already working can move.
     * <p>
     * On the Java 25 tree that equality holds up to the reduction order, which
     * {@code DoubleVector.reduceLanes(ADD)} deliberately leaves unspecified:
     * the same call site there produces two different results one unit in the
     * last place apart depending on whether it is running interpreted or
     * compiled.
     */
    @Test
    public void testTheOrdinaryCaseIsExactlyTheStraightforwardLoop() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            for (int trial = 0; trial < 60; ++trial) {
                double s = Math.scalb(1.0, -120 + trial * 4);
                double[] v = new double[n];
                for (int i = 0; i < n; ++i) {
                    v[i] = next() * s;
                }
                double naive = naiveTwoNorm(v);
                if (naive == 0.0 || Double.isInfinite(naive)) {
                    continue;
                }
                double actual = VectorOps.twoNorm(v);
                if (VectorOps.isVectorized()) {
                    assertEquals("n = " + n + " trial " + trial, naive, actual, 1.0e-14 * naive);
                } else {
                    assertSameBits("n = " + n + " trial " + trial, naive, actual);
                }
            }
        }
    }

    /**
     * {@code ||c * v|| == |c| * ||v||} for {@code c} a power of two, which
     * scales every squared term by the same exact factor and so cannot change
     * the rounding -- exactly on the scalar tree, and up to the unspecified
     * reduction order on the vectorized one, where the two calls being
     * compared need not have been compiled alike.
     */
    @Test
    public void testTheTwoNormIsAbsolutelyHomogeneous() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            double[] v = new double[n];
            for (int i = 0; i < n; ++i) {
                v[i] = next();
            }
            for (int k = -300; k <= 300; k += 20) {
                double c = Math.scalb(1.0, k);
                double[] scaled = new double[n];
                for (int i = 0; i < n; ++i) {
                    scaled[i] = v[i] * c;
                }
                double expected = VectorOps.twoNorm(v) * c;
                if (expected == 0.0 || Double.isInfinite(expected)) {
                    continue;
                }
                double actual = VectorOps.twoNorm(scaled);
                String where = "n = " + n + " scaled by 2^" + k;
                if (VectorOps.isVectorized()) {
                    assertEquals(where, expected, actual, 1.0e-15 * expected);
                } else {
                    assertSameBits(where, expected, actual);
                }
            }
        }
    }

    /** {@code max |v_i| <= ||v||_2 <= sum |v_i|}, the ordering of the three. */
    @Test
    public void testTheThreeNormsAreOrdered() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            for (int trial = 0; trial < 40; ++trial) {
                double[] v = new double[n];
                for (int i = 0; i < n; ++i) {
                    v[i] = next() * Math.scalb(1.0, -30 + trial);
                }
                double inf = VectorOps.infinityNorm(v);
                double two = VectorOps.twoNorm(v);
                double one = VectorOps.absNorm(v);
                assertTrue("infinityNorm > twoNorm at n = " + n, inf <= two * (1.0 + 1.0e-15));
                assertTrue("twoNorm > absNorm at n = " + n, two <= one * (1.0 + 1.0e-15));
            }
        }
    }

    @Test
    public void testTheTriangleInequalityHolds() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            double[] a = new double[n];
            double[] b = new double[n];
            double[] s = new double[n];
            for (int i = 0; i < n; ++i) {
                a[i] = next() * 1.0e12;
                b[i] = next() * 1.0e-12;
                s[i] = a[i] + b[i];
            }
            double lhs = VectorOps.twoNorm(s);
            double rhs = VectorOps.twoNorm(a) + VectorOps.twoNorm(b);
            assertTrue("n = " + n + " : " + lhs + " > " + rhs, lhs <= rhs * (1.0 + 1.0e-14));
        }
    }

    /** A unit basis vector has norm one, wherever it sits on the exponent axis. */
    @Test
    public void testAScaledBasisVectorHasTheNormOfItsOneEntry() {
        for (int k = -1000; k <= 1000; k += 10) {
            double value = Math.scalb(1.0, k);
            for (int li = 0; li < LENGTHS.length; ++li) {
                int n = LENGTHS[li];
                double[] v = new double[n];
                v[n / 2] = -value;
                assertSameBits("n = " + n + " at 2^" + k, value, VectorOps.twoNorm(v));
            }
        }
    }

    /** All three norms answer {@code 0.0} for an empty vector. */
    @Test
    public void testTheNormsOfAnEmptyVectorAreZero() {
        double[] empty = new double[0];
        assertSameBits("twoNorm", 0.0, VectorOps.twoNorm(empty));
        assertSameBits("absNorm", 0.0, VectorOps.absNorm(empty));
        assertSameBits("infinityNorm", 0.0, VectorOps.infinityNorm(empty));
        assertSameBits("sum", 0.0, VectorOps.sum(empty));
    }

    /** And for a vector of zeros, whichever sign the zeros carry. */
    @Test
    public void testTheNormsOfAZeroVectorAreZero() {
        double[] zeros = { 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0 };
        assertSameBits("twoNorm", 0.0, VectorOps.twoNorm(zeros));
        assertSameBits("absNorm", 0.0, VectorOps.absNorm(zeros));
        assertSameBits("infinityNorm", 0.0, VectorOps.infinityNorm(zeros));
    }

    /** A {@code NaN} in any position reaches every norm. */
    @Test
    public void testANaNAnywhereReachesEveryNorm() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            for (int pos = 0; pos < n; ++pos) {
                double[] v = new double[n];
                for (int i = 0; i < n; ++i) {
                    v[i] = i + 1.0;
                }
                v[pos] = Double.NaN;
                String where = "n = " + n + ", NaN at " + pos;
                assertTrue(where + " : twoNorm", Double.isNaN(VectorOps.twoNorm(v)));
                assertTrue(where + " : absNorm", Double.isNaN(VectorOps.absNorm(v)));
                assertTrue(where + " : infinityNorm", Double.isNaN(VectorOps.infinityNorm(v)));
                assertTrue(where + " : isNaN", VectorOps.isNaN(v));
            }
        }
    }

    /**
     * An infinity reaches every norm as an infinity, and a {@code NaN} beside
     * one still wins -- {@code NaN + Infinity} is {@code NaN}, and that is what
     * the straightforward loop has always returned here.
     */
    @Test
    public void testAnInfinityReachesEveryNorm() {
        double[] v = { 1.0, Double.POSITIVE_INFINITY, 2.0 };
        assertEquals(Double.POSITIVE_INFINITY, VectorOps.twoNorm(v), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, VectorOps.absNorm(v), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, VectorOps.infinityNorm(v), 0.0);

        double[] w = { 1.0, Double.NEGATIVE_INFINITY, 2.0 };
        assertEquals(Double.POSITIVE_INFINITY, VectorOps.twoNorm(w), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, VectorOps.infinityNorm(w), 0.0);

        double[] both = { Double.NaN, Double.POSITIVE_INFINITY };
        assertTrue("NaN beside an infinity must reach the result",
                Double.isNaN(VectorOps.twoNorm(both)));
    }

    /** Subnormal input, where the squares are far below the smallest normal. */
    @Test
    public void testSubnormalInputStillHasANorm() {
        double[] v = { Double.MIN_VALUE, Double.MIN_VALUE, Double.MIN_VALUE, Double.MIN_VALUE };
        double norm = VectorOps.twoNorm(v);
        assertTrue("the norm of a non-zero vector must not be zero", norm > 0.0);
        assertEquals(enorm(v), norm, 0.0);

        double[] w = { Double.MIN_NORMAL, Double.MIN_NORMAL };
        assertEquals(enorm(w), VectorOps.twoNorm(w), 0.0);
    }

    // ------------------------------------------------------------------
    // sum, and the one-norm it is not
    // ------------------------------------------------------------------

    /**
     * {@code sum} adds the elements as they are, {@code absNorm} adds their
     * absolute values. The deprecated {@code oneNorm} is the former under a
     * name that promises the latter, which is what the deprecation says.
     */
    @Test
    @SuppressWarnings("deprecation")
    public void testSumIsSignedAndAbsNormIsTheOneNorm() {
        double[] v = { 3.0, -4.0, 5.0, -5.0 };
        assertEquals("the signed sum", -1.0, VectorOps.sum(v), 0.0);
        assertEquals("the one-norm", 17.0, VectorOps.absNorm(v), 0.0);
        assertEquals("oneNorm still returns what it always returned", -1.0, VectorOps.oneNorm(v), 0.0);

        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            double[] w = new double[n];
            for (int i = 0; i < n; ++i) {
                w[i] = next() * 1.0e3;
            }
            assertSameBits("n = " + n, VectorOps.sum(w), VectorOps.oneNorm(w));
        }

        double[] positive = { 1.0, 2.0, 3.0, 4.0 };
        assertSameBits("with no negative entry the two coincide", VectorOps.absNorm(positive),
                VectorOps.sum(positive));
    }

    // ------------------------------------------------------------------
    // the remaining members, none of which had a test
    // ------------------------------------------------------------------

    @Test
    public void testSetAll() {
        for (int li = 0; li < LENGTHS.length; ++li) {
            double[] v = new double[LENGTHS[li]];
            VectorOps.setAll(v, -2.75);
            for (int i = 0; i < v.length; ++i) {
                assertSameBits("element " + i, -2.75, v[i]);
            }
        }
    }

    @Test
    public void testSet() {
        double[] source = { 1.0, -2.0, 3.0 };
        double[] dest = new double[3];
        VectorOps.set(dest, source);
        assertArrayEquals(source, dest, 0.0);
        // a copy, not an alias
        source[0] = 99.0;
        assertEquals(1.0, dest[0], 0.0);

        try {
            VectorOps.set(new double[2], new double[3]);
            fail("a length mismatch was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
    }

    @Test
    public void testAbsNormalize() {
        double[] v = { 3.0, -4.0, 5.0, -5.0 };
        double returned = VectorOps.absNormalize(v);
        assertEquals("the return value is the norm it divided by", 17.0, returned, 0.0);
        assertEquals("and afterwards the one-norm is one", 1.0, VectorOps.absNorm(v), 1.0e-15);
        assertEquals(3.0 / 17.0, v[0], 1.0e-15);
        assertEquals(-4.0 / 17.0, v[1], 1.0e-15);

        double[] zeros = { 0.0, 0.0, 0.0 };
        assertEquals("a zero vector is left alone", 0.0, VectorOps.absNormalize(zeros), 0.0);
        assertArrayEquals(new double[] { 0.0, 0.0, 0.0 }, zeros, 0.0);
    }

    @Test
    public void testIsNaN() {
        assertFalse(VectorOps.isNaN(new double[0]));
        assertFalse(VectorOps.isNaN(new double[] { 1.0, 2.0, Double.POSITIVE_INFINITY }));
        assertTrue(VectorOps.isNaN(new double[] { 1.0, Double.NaN }));
        assertTrue(VectorOps.isNaN(new double[] { Double.NaN }));
    }

    @Test
    public void testMeanMaxMin() {
        double[] v = { 4.0, -2.0, 10.0, 0.0 };
        assertEquals(3.0, VectorOps.mean(v), 1.0e-15);
        assertEquals(10.0, VectorOps.max(v), 0.0);
        assertEquals(-2.0, VectorOps.min(v), 0.0);

        for (int li = 0; li < LENGTHS.length; ++li) {
            int n = LENGTHS[li];
            double[] w = new double[n];
            double sum = 0.0;
            double lo = Double.POSITIVE_INFINITY;
            double hi = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < n; ++i) {
                w[i] = next() * 1.0e6;
                sum += w[i];
                lo = Math.min(lo, w[i]);
                hi = Math.max(hi, w[i]);
            }
            assertEquals("max at n = " + n, hi, VectorOps.max(w), 0.0);
            assertEquals("min at n = " + n, lo, VectorOps.min(w), 0.0);
            assertEquals("mean at n = " + n, sum / n, VectorOps.mean(w), 1.0e-9);
            assertSameBits("infinityNorm at n = " + n, Math.max(Math.abs(hi), Math.abs(lo)),
                    VectorOps.infinityNorm(w));
        }
    }
}
