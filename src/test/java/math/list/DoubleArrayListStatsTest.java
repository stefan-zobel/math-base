package math.list;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.NoSuchElementException;

import org.junit.Test;

/**
 * Tests for the reductions of {@link DoubleArrayList}: the median, which has to
 * agree with the same values held in a root list, the standard deviation, which
 * has to be invariant under a shift of the data, and the euclidean norm, which
 * has to survive data whose squares do not fit into a {@code double}.
 */
public class DoubleArrayListStatsTest {

    /** 60 decimal digits, far beyond anything a double can hold. */
    private static final MathContext EXACT = new MathContext(60);

    /**
     * Relative tolerance for a standard deviation. The worst error measured
     * against the reference below, over 32000 sets of up to 41 points spread
     * across magnitudes from 1e-150 to 1e150 and centers up to 1e15, is
     * 6.2e-16.
     */
    private static final double SD_TOL = 1.0e-14;

    private long lcg = 20260823L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return (lcg >>> 11) * 0x1.0p-53;
    }

    // ----- median -------------------------------------------------------

    @Test
    public void medianOfKnownSubLists() {
        DoubleArrayList list = DoubleArrayList.of(1.0, 2.0, 3.0, 4.0, 5.0);
        assertEquals(3.0, list.median(), 0.0);
        assertEquals(3.0, list.subList(1, 4).median(), 0.0);
        assertEquals(4.0, list.subList(2, 5).median(), 0.0);
        assertEquals(2.5, list.subList(1, 3).median(), 0.0);
        assertEquals(1.0, list.subList(0, 1).median(), 0.0);
        assertEquals(5.0, list.subList(4, 5).median(), 0.0);
        assertEquals(4.5, list.subList(3, 5).median(), 0.0);
    }

    @Test
    public void medianOfASubListMatchesTheSameValuesInARootList() {
        for (int trial = 0; trial < 300; ++trial) {
            int n = 1 + (int) (next() * 24);
            double[] values = new double[n];
            for (int i = 0; i < n; ++i) {
                values[i] = Math.floor(next() * 100.0);
            }
            DoubleArrayList root = DoubleArrayList.of(values);
            for (int from = 0; from < n; ++from) {
                for (int to = from + 1; to <= n; ++to) {
                    double[] range = Arrays.copyOfRange(values, from, to);
                    // the same values, but at offset zero -- the offset is the
                    // only thing that differs, so the two must agree exactly
                    assertEquals(Arrays.toString(range), DoubleArrayList.of(range).median(),
                            root.subList(from, to).median(), 0.0);
                }
            }
        }
    }

    @Test
    public void medianMatchesTheDefinitionOnASortedCopy() {
        for (int trial = 0; trial < 500; ++trial) {
            int n = 1 + (int) (next() * 20);
            double[] values = new double[n];
            for (int i = 0; i < n; ++i) {
                values[i] = next() * 1000.0 - 500.0;
            }
            DoubleArrayList root = DoubleArrayList.of(values);
            assertEquals(medianOfSortedCopy(values), root.median(), 0.0);
            int from = (int) (next() * n);
            int to = from + 1 + (int) (next() * (n - from));
            double[] range = Arrays.copyOfRange(values, from, to);
            assertEquals(medianOfSortedCopy(range), root.subList(from, to).median(), 0.0);
        }
    }

    @Test
    public void medianOfNestedSubLists() {
        DoubleArrayList root = DoubleArrayList.of(9.0, 8.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 6.0);
        DoubleList outer = root.subList(1, 8);
        DoubleList inner = outer.subList(1, 5);
        assertEquals(medianOfSortedCopy(new double[] { 8.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0 }), outer.median(), 0.0);
        assertEquals(medianOfSortedCopy(new double[] { 1.0, 2.0, 3.0, 4.0 }), inner.median(), 0.0);
    }

    @Test
    public void medianOfAnEmptyListThrows() {
        try {
            DoubleArrayList.of().median();
            fail("no exception for an empty root list");
        } catch (NoSuchElementException expected) {
            // expected
        }
        try {
            DoubleArrayList.of(1.0, 2.0, 3.0).subList(1, 1).median();
            fail("no exception for an empty sub-list");
        } catch (NoSuchElementException expected) {
            // expected
        }
    }

    // ----- standard deviation -------------------------------------------

    @Test
    public void stddevIsBitIdenticalUnderAShiftOfRepresentableData() {
        double[] base = { 1.0, 2.0, 3.0, 4.0, 5.0, 3.0, 2.0, 4.0 };
        double[] shifts = { 0.0, 1.0e3, 1.0e6, 1.0e8, 1.0e9, 1.0e10, 1.0e12, 1.0e14, 1.0e15 };
        double unshifted = DoubleArrayList.of(base).stddev();
        for (int k = 0; k < shifts.length; ++k) {
            double[] shifted = new double[base.length];
            for (int i = 0; i < base.length; ++i) {
                shifted[i] = shifts[k] + base[i];
                // every shift is exact for this data, so the standard
                // deviation of the stored values does not move at all
                assertEquals(base[i], shifted[i] - shifts[k], 0.0);
            }
            assertEquals("shift " + shifts[k], unshifted, DoubleArrayList.of(shifted).stddev(), 0.0);
        }
    }

    @Test
    public void stddevIsInvariantUnderAShiftOfRandomData() {
        double[] shifts = { 1.0e3, 1.0e6, 1.0e8, 1.0e9, 1.0e10, 1.0e12, 1.0e14, 1.0e15 };
        for (int trial = 0; trial < 400; ++trial) {
            int n = 2 + (int) (next() * 30);
            double[] base = new double[n];
            for (int i = 0; i < n; ++i) {
                // integers, so that adding the shift loses nothing
                base[i] = Math.floor(next() * 1000.0);
            }
            double unshifted = DoubleArrayList.of(base).stddev();
            if (unshifted == 0.0) {
                continue;
            }
            for (int k = 0; k < shifts.length; ++k) {
                double[] shifted = new double[n];
                for (int i = 0; i < n; ++i) {
                    shifted[i] = shifts[k] + base[i];
                }
                assertClose("shift " + shifts[k], unshifted, DoubleArrayList.of(shifted).stddev(), SD_TOL);
            }
        }
    }

    @Test
    public void stddevScalesExactlyWithAPowerOfTwo() {
        int[] exponents = { -300, -100, -20, -1, 1, 20, 100, 300 };
        for (int trial = 0; trial < 300; ++trial) {
            int n = 2 + (int) (next() * 20);
            double[] base = new double[n];
            for (int i = 0; i < n; ++i) {
                base[i] = next() - 0.5;
            }
            double unscaled = DoubleArrayList.of(base).stddev();
            for (int k = 0; k < exponents.length; ++k) {
                double[] scaled = new double[n];
                for (int i = 0; i < n; ++i) {
                    scaled[i] = Math.scalb(base[i], exponents[k]);
                }
                // scaling by a power of two is exact, and so is the way the
                // computation removes it again
                assertEquals("2^" + exponents[k], Math.scalb(unscaled, exponents[k]),
                        DoubleArrayList.of(scaled).stddev(), 0.0);
            }
        }
    }

    @Test
    public void stddevMatchesAHighPrecisionReference() {
        double[] magnitudes = { 1.0e-150, 1.0e-6, 1.0, 1.0e6, 1.0e150 };
        double[] centers = { 0.0, 1.0, 1.0e8, 1.0e15 };
        for (int m = 0; m < magnitudes.length; ++m) {
            for (int c = 0; c < centers.length; ++c) {
                for (int trial = 0; trial < 60; ++trial) {
                    int n = 2 + (int) (next() * 30);
                    double[] values = new double[n];
                    for (int i = 0; i < n; ++i) {
                        values[i] = (centers[c] + next()) * magnitudes[m];
                    }
                    double reference = referenceStddev(values);
                    if (reference > 0.0 && !Double.isInfinite(reference)) {
                        assertClose("magnitude " + magnitudes[m] + " center " + centers[c], reference,
                                DoubleArrayList.of(values).stddev(), SD_TOL);
                    }
                }
            }
        }
    }

    @Test
    public void stddevOfConstantDataIsExactlyZero() {
        double[] constants = { 0.0, 0.1, 1.0 / 3.0, 123.456, 1.0e8, 1.0e300, -1.0e300, 1.0e-300, Double.MAX_VALUE,
                Double.MIN_VALUE };
        for (int c = 0; c < constants.length; ++c) {
            for (int n = 2; n <= 64; ++n) {
                double[] values = new double[n];
                Arrays.fill(values, constants[c]);
                // never NaN, and never a residue of the constant either
                assertEquals(constants[c] + " x " + n, 0.0, DoubleArrayList.of(values).stddev(), 0.0);
            }
        }
    }

    @Test
    public void stddevSurvivesTheEndsOfTheDoubleRange() {
        // the population standard deviation of {c, 2c, 3c} is c * sqrt(2/3),
        // whatever c is -- an oracle that needs no reference implementation
        double[] cs = { 1.0, 1.0e160, 1.0e300, 1.0e-160, 1.0e-300, 1.0e-320 };
        for (int i = 0; i < cs.length; ++i) {
            double[] values = { cs[i], 2.0 * cs[i], 3.0 * cs[i] };
            assertClose("c = " + cs[i], cs[i] * Math.sqrt(2.0 / 3.0), DoubleArrayList.of(values).stddev(), SD_TOL);
        }
        // the widest spread a double can hold, and one that cannot be squared
        assertEquals(Double.MAX_VALUE, DoubleArrayList.of(-Double.MAX_VALUE, Double.MAX_VALUE).stddev(), 0.0);
        assertEquals(Double.MAX_VALUE / 2.0, DoubleArrayList.of(0.0, Double.MAX_VALUE).stddev(), 0.0);
        // subnormal data, where squaring a deviation gives zero
        assertEquals(Double.MIN_VALUE, DoubleArrayList.of(0.0, 2.0 * Double.MIN_VALUE).stddev(), 0.0);
    }

    @Test
    public void stddevPropagatesNonFiniteValues() {
        assertTrue(Double.isNaN(DoubleArrayList.of(1.0, Double.NaN, 3.0).stddev()));
        assertTrue(Double.isNaN(DoubleArrayList.of(1.0, Double.POSITIVE_INFINITY, 3.0).stddev()));
        assertTrue(Double.isNaN(
                DoubleArrayList.of(Double.NEGATIVE_INFINITY, 0.0, Double.POSITIVE_INFINITY).stddev()));
    }

    @Test
    public void stddevDividesByNAndNeedsTwoElements() {
        // the documented convention is the population standard deviation, so
        // {0, 2} is 1.0 and not sqrt(2)
        assertEquals(1.0, DoubleArrayList.of(0.0, 2.0).stddev(), 0.0);
        assertEquals(0.5, DoubleArrayList.of(1.0, 2.0).stddev(), 0.0);
        // two is the smallest number of observations in which a dispersion
        // is observable at all, so the refusal is deliberate: the 0.0 that
        // the formula would produce would claim there is no spread, where in
        // truth there is nothing to go on
        try {
            DoubleArrayList.of(1.0).stddev();
            fail("no exception for a single element");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("2"));
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("1"));
        }
        // and the same through a sub-list, which shares the check
        try {
            DoubleArrayList.of(1.0, 2.0, 3.0).subList(1, 2).stddev();
            fail("no exception for a single element sub-list");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("1"));
        }
        // while the reductions that are defined on one observation give it
        assertEquals(3.5, DoubleArrayList.of(3.5).avg(), 0.0);
        assertEquals(3.5, DoubleArrayList.of(3.5).median(), 0.0);
    }

    @Test
    public void theEmptyListSaysWhichWayItIsWrong() {
        DoubleArrayList empty = DoubleArrayList.of();
        // "there is no element to return" and "there are not enough elements
        // to compute with" are two different situations, and the interface
        // documents a different exception for each
        try {
            empty.min();
            fail("no exception for min on an empty list");
        } catch (NoSuchElementException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            empty.avg();
            fail("no exception for avg on an empty list");
        } catch (NoSuchElementException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            empty.stddev();
            fail("no exception for stddev on an empty list");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("0"));
        }
        try {
            empty.iqr();
            fail("no exception for iqr on an empty list");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().endsWith("0"));
        }
    }

    @Test
    public void stddevOfASubListMatchesTheSameValuesInARootList() {
        for (int trial = 0; trial < 60; ++trial) {
            int n = 2 + (int) (next() * 20);
            double[] values = new double[n];
            for (int i = 0; i < n; ++i) {
                values[i] = 1.0e9 + next();
            }
            DoubleArrayList root = DoubleArrayList.of(values);
            for (int from = 0; from < n - 1; ++from) {
                for (int to = from + 2; to <= n; ++to) {
                    double[] range = Arrays.copyOfRange(values, from, to);
                    assertEquals(from + ".." + to, DoubleArrayList.of(range).stddev(),
                            root.subList(from, to).stddev(), 0.0);
                }
            }
        }
    }

    // ----- euclidean norm -----------------------------------------------

    @Test
    public void norm2OfKnownVectors() {
        assertEquals(0.0, DoubleArrayList.of().norm2(), 0.0);
        assertEquals(0.0, DoubleArrayList.of(0.0, 0.0, 0.0).norm2(), 0.0);
        assertEquals(5.0, DoubleArrayList.of(3.0, 4.0).norm2(), 0.0);
        assertEquals(13.0, DoubleArrayList.of(3.0, 4.0, 12.0).norm2(), 0.0);
        assertEquals(2.0, DoubleArrayList.of(1.0, 1.0, 1.0, 1.0).norm2(), 0.0);
    }

    @Test
    public void norm2IsUnchangedWhereTheDirectFormWorks() {
        for (int trial = 0; trial < 3000; ++trial) {
            int n = 1 + (int) (next() * 40);
            double[] values = new double[n];
            double magnitude = Math.pow(10.0, (next() - 0.5) * 200.0);
            for (int i = 0; i < n; ++i) {
                values[i] = (next() - 0.5) * magnitude;
            }
            double squares = 0.0;
            for (int i = 0; i < n; ++i) {
                squares += values[i] * values[i];
            }
            // wherever the sum of squares stays inside the double range, the
            // answer must be the one the direct route gives, to the last bit
            assertEquals("magnitude " + magnitude, Math.sqrt(squares), DoubleArrayList.of(values).norm2(), 0.0);
        }
    }

    @Test
    public void norm2SurvivesTheEndsOfTheDoubleRange() {
        // the norm of {c, 2c, 3c} is c * sqrt(14), whatever c is
        double[] cs = { 1.0, 1.0e160, 1.0e300, 1.0e-160, 1.0e-300, 1.0e-310 };
        for (int i = 0; i < cs.length; ++i) {
            double[] values = { cs[i], 2.0 * cs[i], 3.0 * cs[i] };
            assertClose("c = " + cs[i], cs[i] * Math.sqrt(14.0), DoubleArrayList.of(values).norm2(), SD_TOL);
        }
        // a single element decides the answer, and its square does not fit
        assertEquals(1.0e200, DoubleArrayList.of(1.0e200, 1.0).norm2(), 0.0);
        assertEquals(Double.MAX_VALUE, DoubleArrayList.of(Double.MAX_VALUE, 0.0).norm2(), 0.0);
        assertEquals(5.0e-200, DoubleArrayList.of(3.0e-200, 4.0e-200).norm2(), 0.0);
        assertEquals(Double.MIN_VALUE, DoubleArrayList.of(Double.MIN_VALUE).norm2(), 0.0);
    }

    @Test
    public void norm2MatchesAHighPrecisionReference() {
        double[] magnitudes = { 1.0e-300, 1.0e-160, 1.0e-6, 1.0, 1.0e6, 1.0e160, 1.0e300 };
        for (int m = 0; m < magnitudes.length; ++m) {
            for (int trial = 0; trial < 100; ++trial) {
                int n = 1 + (int) (next() * 40);
                double[] values = new double[n];
                for (int i = 0; i < n; ++i) {
                    values[i] = (next() - 0.5) * magnitudes[m];
                }
                double reference = referenceNorm2(values);
                if (reference > 0.0) {
                    assertClose("magnitude " + magnitudes[m], reference, DoubleArrayList.of(values).norm2(), SD_TOL);
                }
            }
        }
    }

    @Test
    public void norm2ScalesExactlyWithAPowerOfTwo() {
        int[] exponents = { -300, -100, -20, -1, 1, 20, 100, 300 };
        for (int trial = 0; trial < 300; ++trial) {
            int n = 1 + (int) (next() * 20);
            double[] base = new double[n];
            for (int i = 0; i < n; ++i) {
                base[i] = next() - 0.5;
            }
            double unscaled = DoubleArrayList.of(base).norm2();
            for (int k = 0; k < exponents.length; ++k) {
                double[] scaled = new double[n];
                for (int i = 0; i < n; ++i) {
                    scaled[i] = Math.scalb(base[i], exponents[k]);
                }
                assertEquals("2^" + exponents[k], Math.scalb(unscaled, exponents[k]),
                        DoubleArrayList.of(scaled).norm2(), 0.0);
            }
        }
    }

    @Test
    public void norm2OfASubListMatchesTheSameValuesInARootList() {
        for (int trial = 0; trial < 100; ++trial) {
            int n = 1 + (int) (next() * 20);
            double[] values = new double[n];
            for (int i = 0; i < n; ++i) {
                // small enough that every square is subnormal, which is where
                // the sub-list and the root list could most easily part ways
                values[i] = (next() - 0.5) * 1.0e-170;
            }
            DoubleArrayList root = DoubleArrayList.of(values);
            for (int from = 0; from < n; ++from) {
                for (int to = from + 1; to <= n; ++to) {
                    double[] range = Arrays.copyOfRange(values, from, to);
                    assertEquals(from + ".." + to, DoubleArrayList.of(range).norm2(),
                            root.subList(from, to).norm2(), 0.0);
                }
            }
        }
    }

    @Test
    public void norm2PropagatesNonFiniteValues() {
        assertTrue(Double.isNaN(DoubleArrayList.of(1.0, Double.NaN).norm2()));
        assertEquals(Double.POSITIVE_INFINITY, DoubleArrayList.of(1.0, Double.POSITIVE_INFINITY).norm2(), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, DoubleArrayList.of(Double.NEGATIVE_INFINITY, 1.0).norm2(), 0.0);
    }

    // ----- helpers ------------------------------------------------------

    private static double medianOfSortedCopy(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        int n = sorted.length;
        if (n % 2 == 0) {
            return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
        }
        return sorted[n / 2];
    }

    /** The population standard deviation, carried at 60 decimal digits. */
    private static double referenceStddev(double[] values) {
        BigDecimal n = new BigDecimal(values.length);
        BigDecimal sum = BigDecimal.ZERO;
        for (int i = 0; i < values.length; ++i) {
            sum = sum.add(new BigDecimal(values[i]));
        }
        BigDecimal mean = sum.divide(n, EXACT);
        BigDecimal squares = BigDecimal.ZERO;
        for (int i = 0; i < values.length; ++i) {
            BigDecimal deviation = new BigDecimal(values[i]).subtract(mean);
            squares = squares.add(deviation.multiply(deviation));
        }
        return Math.sqrt(squares.divide(n, EXACT).doubleValue());
    }

    /**
     * The euclidean norm at 60 decimal digits, with a power of two factored
     * out so that the reference itself does not leave the double range.
     */
    private static double referenceNorm2(double[] values) {
        double max = 0.0;
        for (int i = 0; i < values.length; ++i) {
            max = Math.max(max, Math.abs(values[i]));
        }
        if (max == 0.0) {
            return 0.0;
        }
        int exponent = Math.getExponent(max);
        BigDecimal squares = BigDecimal.ZERO;
        for (int i = 0; i < values.length; ++i) {
            BigDecimal scaled = new BigDecimal(Math.scalb(values[i], -exponent));
            squares = squares.add(scaled.multiply(scaled));
        }
        return Math.scalb(Math.sqrt(squares.doubleValue()), exponent);
    }

    private static void assertClose(String what, double expected, double actual, double tolerance) {
        double error = Math.abs(actual - expected) / Math.abs(expected);
        if (!(error <= tolerance)) {
            fail(what + ": expected " + expected + ", got " + actual + ", relative error " + error);
        }
    }
}
