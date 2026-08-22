package math.linalg;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

public class StandardizationTest {

    private long lcg = 13572468135724681L;

    /** Deterministic pseudo random number in [-1, 1). */
    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) / (double) (1L << 53)) * 2.0 - 1.0;
    }

    /** A table whose columns sit on wildly different scales and offsets. */
    private double[][] table(int n, int p) {
        double[][] x = new double[n][p];
        for (int j = 0; j < p; ++j) {
            double offset = Math.pow(10.0, (j % 7) - 3);
            double spread = Math.pow(10.0, (j % 5) - 2);
            for (int i = 0; i < n; ++i) {
                x[i][j] = offset + spread * next();
            }
        }
        return x;
    }

    private static double columnMean(double[][] x, int j) {
        double sum = 0.0;
        for (int i = 0; i < x.length; ++i) {
            sum += x[i][j];
        }
        return sum / x.length;
    }

    private static double columnVariance(double[][] x, int j, double divisor) {
        double mean = columnMean(x, j);
        double ss = 0.0;
        for (int i = 0; i < x.length; ++i) {
            ss += (x[i][j] - mean) * (x[i][j] - mean);
        }
        return ss / divisor;
    }

    // ------------------------------------------------------------------

    /**
     * The defining property: after the transform every column has mean zero and
     * variance one, measured at the divisor the transform was fitted with.
     * <p>
     * The tolerance carries the conditioning of the centering rather than being
     * a flat number. Subtracting a mean of {@code 100} from values that spread
     * over {@code 0.01} throws away four digits before anything else happens,
     * and columns exactly like that are in the table below on purpose, so the
     * bound is {@code |mean| / scale} units of rounding error. A flat
     * {@code 1e-12} passes on the well-scaled columns and fails on those, which
     * says nothing about the code.
     */
    @Test
    public void testEveryColumnComesOutWithZeroMeanAndUnitVariance() {
        int[] rows = { 2, 3, 10, 178, 1000 };
        int[] cols = { 1, 2, 13, 30 };
        for (int ri = 0; ri < rows.length; ++ri) {
            for (int ci = 0; ci < cols.length; ++ci) {
                int n = rows[ri];
                int p = cols[ci];
                double[][] x = table(n, p);
                Standardization s = Standardization.of(x);
                double[] mean = s.mean();
                double[] scale = s.scale();
                double[][] z = s.transform(x);
                for (int j = 0; j < p; ++j) {
                    String where = "n = " + n + ", column " + j;
                    double conditioning = Math.max(1.0, Math.abs(mean[j]) / scale[j]);
                    assertEquals(where + " mean", 0.0, columnMean(z, j), 1.0e-12 * conditioning);
                    if (n > 1) {
                        assertEquals(where + " variance", 1.0, columnVariance(z, j, n - 1),
                                1.0e-12 * conditioning);
                    }
                }
                Standardization sp = Standardization.ofPopulation(x);
                double[] meanP = sp.mean();
                double[] scaleP = sp.scale();
                double[][] zp = sp.transform(x);
                for (int j = 0; j < p; ++j) {
                    double conditioning = Math.max(1.0, Math.abs(meanP[j]) / scaleP[j]);
                    assertEquals("population variance, n = " + n + ", column " + j, 1.0,
                            columnVariance(zp, j, n), 1.0e-12 * conditioning);
                }
            }
        }
    }

    /**
     * The two divisors differ by exactly {@code sqrt((n - 1) / n)} and by
     * nothing else, which is the whole of the distinction between them.
     */
    @Test
    public void testThePopulationDivisorDiffersOnlyByTheKnownFactor() {
        int[] rows = { 2, 5, 100, 1000 };
        for (int ri = 0; ri < rows.length; ++ri) {
            int n = rows[ri];
            double[][] x = table(n, 6);
            double[] sample = Standardization.of(x).scale();
            double[] population = Standardization.ofPopulation(x).scale();
            double factor = Math.sqrt((n - 1.0) / n);
            for (int j = 0; j < sample.length; ++j) {
                assertEquals("n = " + n + ", column " + j, sample[j] * factor, population[j],
                        1.0e-14 * sample[j]);
            }
            // and the means do not depend on the divisor at all
            assertArrayEquals("n = " + n, Standardization.of(x).mean(),
                    Standardization.ofPopulation(x).mean(), 0.0);
        }
    }

    /** {@code inverse(transform(x))} is {@code x} again. */
    @Test
    public void testTheInverseRecoversTheOriginal() {
        double[][] x = table(200, 9);
        Standardization s = Standardization.of(x);
        double[][] back = s.inverse(s.transform(x));
        for (int i = 0; i < x.length; ++i) {
            for (int j = 0; j < x[i].length; ++j) {
                assertEquals("(" + i + ", " + j + ")", x[i][j], back[i][j],
                        1.0e-12 * Math.max(Math.abs(x[i][j]), 1.0e-12));
            }
        }
    }

    /** The row forms are the matrix forms, one row at a time. */
    @Test
    public void testTheRowFormsAgreeWithTheMatrixForms() {
        double[][] x = table(50, 7);
        Standardization s = Standardization.of(x);
        double[][] z = s.transform(x);
        double[][] back = s.inverse(z);
        for (int i = 0; i < x.length; ++i) {
            assertArrayEquals("row " + i, z[i], s.transformRow(x[i]), 0.0);
            assertArrayEquals("row " + i, back[i], s.inverseRow(z[i]), 0.0);
        }
    }

    /**
     * The point of it being a fitted transform: applied to data it was not
     * fitted on it uses the mean and scale it was fitted with, not the ones
     * that batch happens to have. Standardizing each batch on its own is what
     * this class exists to make hard, so a test says which of the two it does.
     */
    @Test
    public void testTheTransformIsTheFittedOneNotARecomputedOne() {
        double[][] training = table(300, 5);
        double[][] holdout = table(40, 5);
        // move the hold-out batch decisively off the training scale
        for (int i = 0; i < holdout.length; ++i) {
            for (int j = 0; j < holdout[i].length; ++j) {
                holdout[i][j] = holdout[i][j] * 3.0 + 17.0;
            }
        }
        Standardization s = Standardization.of(training);
        double[] mean = s.mean();
        double[] scale = s.scale();
        double[][] z = s.transform(holdout);
        for (int i = 0; i < holdout.length; ++i) {
            for (int j = 0; j < holdout[i].length; ++j) {
                assertEquals("(" + i + ", " + j + ")", (holdout[i][j] - mean[j]) / scale[j], z[i][j], 0.0);
            }
        }
        // and it is emphatically not the hold-out batch's own standardization
        for (int j = 0; j < holdout[0].length; ++j) {
            assertTrue("column " + j + " came out with mean 0, so it was refitted",
                    Math.abs(columnMean(z, j)) > 1.0e-6);
        }
    }

    /** Shifting a column by a constant does not move its z-scores. */
    @Test
    public void testShiftingAColumnLeavesTheResultUnchanged() {
        double[][] x = table(120, 4);
        double[][] z = Standardization.of(x).transform(x);
        double[] shifts = { 1.0, -25.0, 1.0e4 };
        for (int k = 0; k < shifts.length; ++k) {
            double[][] y = new double[x.length][];
            for (int i = 0; i < x.length; ++i) {
                y[i] = x[i].clone();
                y[i][1] += shifts[k];
            }
            double[][] zy = Standardization.of(y).transform(y);
            for (int i = 0; i < x.length; ++i) {
                assertEquals("shift " + shifts[k] + ", row " + i, z[i][1], zy[i][1], 1.0e-9);
            }
        }
    }

    /** And neither does scaling it by a positive constant. */
    @Test
    public void testScalingAColumnLeavesTheResultUnchanged() {
        double[][] x = table(120, 4);
        double[][] z = Standardization.of(x).transform(x);
        double[] factors = { 2.0, 0.001, 1.0e6 };
        for (int k = 0; k < factors.length; ++k) {
            double[][] y = new double[x.length][];
            for (int i = 0; i < x.length; ++i) {
                y[i] = x[i].clone();
                y[i][2] *= factors[k];
            }
            double[][] zy = Standardization.of(y).transform(y);
            for (int i = 0; i < x.length; ++i) {
                assertEquals("factor " + factors[k] + ", row " + i, z[i][2], zy[i][2], 1.0e-9);
            }
        }
    }

    /** Nothing handed in is written to. */
    @Test
    public void testTheArgumentIsNotModified() {
        double[][] x = table(30, 5);
        double[][] copy = new double[x.length][];
        for (int i = 0; i < x.length; ++i) {
            copy[i] = x[i].clone();
        }
        Standardization s = Standardization.of(x);
        double[][] z = s.transform(x);
        s.inverse(z);
        s.transformRow(x[0]);
        s.inverseRow(z[0]);
        Standardization.standardize(x[0]);
        for (int i = 0; i < x.length; ++i) {
            assertArrayEquals("row " + i, copy[i], x[i], 0.0);
        }
        assertNotSame("transform must not hand back its argument", x, z);
    }

    /** {@code mean()} and {@code scale()} hand out copies, not the state. */
    @Test
    public void testMeanAndScaleAreCopies() {
        double[][] x = table(20, 3);
        Standardization s = Standardization.of(x);
        double[] mean = s.mean();
        double[] scale = s.scale();
        assertNotSame(mean, s.mean());
        assertNotSame(scale, s.scale());
        mean[0] = 999.0;
        scale[0] = 999.0;
        assertTrue("the state was reachable through mean()", s.mean()[0] != 999.0);
        assertTrue("the state was reachable through scale()", s.scale()[0] != 999.0);
        assertEquals(3, s.dimension());
    }

    /** A constant column has no scale, and saying so beats dividing by zero. */
    @Test
    public void testAConstantColumnIsRefused() {
        double[][] constant = { { 1.0, 7.0 }, { 2.0, 7.0 }, { 3.0, 7.0 } };
        try {
            Standardization.of(constant);
            fail("a constant column was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("the message must name the column: " + expected.getMessage(),
                    expected.getMessage().contains("column 1"));
        }
        double[][] zeros = { { 0.0, 1.0 }, { 0.0, 2.0 } };
        try {
            Standardization.ofPopulation(zeros);
            fail("an all-zero column was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("column 0"));
        }
        try {
            Standardization.standardize(new double[] { 4.0, 4.0, 4.0 });
            fail("a constant vector was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        // a single row has no sample spread either, whatever it holds
        try {
            Standardization.of(new double[][] { { 1.0, 2.0 } });
            fail("a single row was accepted for a sample standard deviation");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        // but the population form is defined there, and refuses for the other
        // reason: one row is constant in every column
        try {
            Standardization.ofPopulation(new double[][] { { 1.0, 2.0 } });
            fail("a single row was accepted for a population standard deviation");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("column 0"));
        }
    }

    /**
     * A non-finite entry is refused, and the message says so rather than
     * calling the column constant. It arrives as a {@code NaN} scale, which is
     * the same branch a constant column would take if the two were not
     * separated, so the wrong message here is easy to ship and hard to read
     * later.
     */
    @Test
    public void testANonFiniteEntryIsNamedAsOneRatherThanCalledConstant() {
        double[] poison = { Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
        for (int k = 0; k < poison.length; ++k) {
            double[][] x = { { 1.0, 4.0 }, { 2.0, 5.0 }, { 3.0, 6.0 } };
            x[1][1] = poison[k];
            try {
                Standardization.of(x);
                fail(poison[k] + " was accepted");
            } catch (IllegalArgumentException expected) {
                String message = expected.getMessage();
                assertTrue(poison[k] + " : " + message, message.contains("column 1"));
                assertTrue(poison[k] + " : " + message, message.contains("not finite"));
                assertTrue(poison[k] + " must not be called constant : " + message,
                        !message.contains("constant"));
            }
            try {
                Standardization.standardize(new double[] { 1.0, poison[k], 3.0 });
                fail(poison[k] + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage(), expected.getMessage().contains("not finite"));
            }
        }
    }

    /** The one-variable convenience is the one-column table. */
    @Test
    public void testStandardizeAgreesWithTheOneColumnTable() {
        int[] lengths = { 2, 3, 17, 500 };
        for (int li = 0; li < lengths.length; ++li) {
            int n = lengths[li];
            double[] v = new double[n];
            for (int i = 0; i < n; ++i) {
                v[i] = 1.0e3 + next() * 5.0;
            }
            double[][] asTable = new double[n][1];
            for (int i = 0; i < n; ++i) {
                asTable[i][0] = v[i];
            }
            double[][] z = Standardization.of(asTable).transform(asTable);
            double[] direct = Standardization.standardize(v);
            for (int i = 0; i < n; ++i) {
                assertEquals("n = " + n + ", element " + i, z[i][0], direct[i], 0.0);
            }
        }
    }

    @Test
    public void testArgumentValidation() {
        try {
            Standardization.of(null);
            fail("null was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            Standardization.of(new double[0][]);
            fail("an empty table was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            Standardization.of(new double[][] { {}, {} });
            fail("a table with no columns was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            Standardization.of(new double[][] { { 1.0, 2.0 }, { 3.0 } });
            fail("a ragged table was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("row 1"));
        }
        try {
            Standardization.of(new double[][] { { 1.0, 2.0 }, null });
            fail("a null row was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("row 1"));
        }
        try {
            Standardization.standardize(null);
            fail("null was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            Standardization.standardize(new double[] { 1.0 });
            fail("a single value was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }

        Standardization s = Standardization.of(table(10, 4));
        try {
            s.transform(null);
            fail("null was accepted");
        } catch (IllegalArgumentException expected) {
            assertNotNull(expected.getMessage());
        }
        try {
            s.transformRow(new double[3]);
            fail("a row of the wrong length was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("fitted on 4"));
        }
        try {
            s.inverseRow(new double[5]);
            fail("a row of the wrong length was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("fitted on 4"));
        }
        try {
            s.transform(new double[][] { new double[4], new double[2] });
            fail("a block with a wrong row was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("row 1"));
        }
    }
}
