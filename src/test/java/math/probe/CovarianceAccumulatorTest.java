package math.probe;

import static org.junit.Assert.*;

import org.junit.Test;

public class CovarianceAccumulatorTest {

    private static final double EPSILON = 1e-12;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /** Deterministic m x n row-major data with entries in [offset - 0.5, offset + 0.5). */
    private static double[][] random(int m, int n, long seed, double offset) {
        double[][] d = new double[m][n];
        long lcg = seed;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                d[i][j] = ((lcg >>> 11) / (double) (1L << 53)) - 0.5 + offset;
            }
        }
        return d;
    }

    /** Straightforward two-pass covariance, column-major, for cross-checking. */
    private static double[] twoPassCovariance(double[][] data) {
        int m = data.length;
        int n = data[0].length;
        double[] mean = new double[n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) mean[j] += data[i][j];
        for (int j = 0; j < n; j++) mean[j] /= m;

        double[] c = new double[n * n];
        for (int p = 0; p < n; p++) {
            for (int q = 0; q < n; q++) {
                double s = 0.0;
                for (int i = 0; i < m; i++) s += (data[i][p] - mean[p]) * (data[i][q] - mean[q]);
                c[q * n + p] = s / (m - 1);
            }
        }
        return c;
    }

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testMatchesTwoPassCovariance() {
        double[][] data = random(200, 5, 4711L, 0.0);
        CovarianceAccumulator acc = new CovarianceAccumulator(5);
        acc.addAll(data);
        assertEquals(200L, acc.count());
        assertEquals(5, acc.dimension());
        assertArrayEquals(twoPassCovariance(data), acc.covariance(), 1e-13);
    }

    @Test
    public void testEmittedMatrixIsExactlySymmetric() {
        // the symmetric eigensolver downstream relies on both triangles being
        // bit-identical, not merely equal up to rounding
        double[][] data = random(97, 9, 13579L, 3.5);
        CovarianceAccumulator acc = new CovarianceAccumulator(9);
        acc.addAll(data);
        double[] c = acc.covariance();
        for (int j = 0; j < 9; j++)
            for (int i = 0; i < 9; i++)
                assertEquals("entry (" + i + "," + j + ")", c[j * 9 + i], c[i * 9 + j], 0.0);
    }

    @Test
    public void testKnownCovarianceOfTwoFeatures() {
        // y = 2x exactly -> var(y) = 4 var(x), cov(x,y) = 2 var(x), correlation 1
        double[][] data = { { 1.0, 2.0 }, { 2.0, 4.0 }, { 3.0, 6.0 }, { 4.0, 8.0 } };
        CovarianceAccumulator acc = new CovarianceAccumulator(2);
        acc.addAll(data);
        double[] c = acc.covariance();
        assertEquals(5.0 / 3.0, c[0], EPSILON);         // var(x) for 1..4
        assertEquals(2.0 * 5.0 / 3.0, c[2], EPSILON);   // cov(x,y)
        assertEquals(4.0 * 5.0 / 3.0, c[3], EPSILON);   // var(y)
        assertArrayEquals(new double[] { 2.5, 5.0 }, acc.mean(), EPSILON);
        double[] r = acc.correlation();
        assertArrayEquals(new double[] { 1.0, 1.0, 1.0, 1.0 }, r, 1e-14);
    }

    @Test
    public void testAccurateForDataWithHugeOffset() {
        // the point of Welford: a mean of 1e9 against unit spread destroys the
        // naive sum(x*y) - n*mx*my formula, but must not affect this one
        double[][] centered = random(500, 4, 24680L, 0.0);
        double[][] shifted = new double[500][4];
        for (int i = 0; i < 500; i++)
            for (int j = 0; j < 4; j++) shifted[i][j] = centered[i][j] + 1.0e9;

        CovarianceAccumulator a = new CovarianceAccumulator(4);
        a.addAll(centered);
        CovarianceAccumulator b = new CovarianceAccumulator(4);
        b.addAll(shifted);

        // Covariance is shift invariant, so both must agree. The residual (~1e-8
        // against entries of ~8e-2) is the representation error of the shifted
        // data itself - ulp(1e9) is 1.2e-7 - not an error of the algorithm. The
        // naive formula would be off by ~1e4 here, since sum(x*y) is 5e20 while
        // the wanted difference is only 41.
        assertArrayEquals(a.covariance(), b.covariance(), 1e-7);
        for (int j = 0; j < 4; j++)
            assertEquals(1.0e9, b.mean()[j] - a.mean()[j], 1e-6);
    }

    @Test
    public void testMergeMatchesSequentialAccumulation() {
        double[][] data = random(300, 6, 99991L, 7.0);
        CovarianceAccumulator whole = new CovarianceAccumulator(6);
        whole.addAll(data);

        CovarianceAccumulator left = new CovarianceAccumulator(6);
        CovarianceAccumulator right = new CovarianceAccumulator(6);
        for (int i = 0; i < 120; i++) left.add(data[i]);
        for (int i = 120; i < 300; i++) right.add(data[i]);
        left.merge(right);

        assertEquals(whole.count(), left.count());
        assertArrayEquals(whole.mean(), left.mean(), 1e-13);
        assertArrayEquals(whole.covariance(), left.covariance(), 1e-12);
        // the absorbed accumulator must be untouched
        assertEquals(180L, right.count());
    }

    @Test
    public void testMergeWithEmptyOnEitherSide() {
        double[][] data = random(50, 3, 2024L, 1.0);
        CovarianceAccumulator filled = new CovarianceAccumulator(3);
        filled.addAll(data);

        CovarianceAccumulator emptyRight = new CovarianceAccumulator(3);
        emptyRight.addAll(data);
        emptyRight.merge(new CovarianceAccumulator(3));
        assertArrayEquals(filled.covariance(), emptyRight.covariance(), 0.0);

        CovarianceAccumulator emptyLeft = new CovarianceAccumulator(3);
        emptyLeft.merge(filled);
        assertEquals(50L, emptyLeft.count());
        assertArrayEquals(filled.covariance(), emptyLeft.covariance(), 0.0);
        assertArrayEquals(filled.mean(), emptyLeft.mean(), 0.0);
    }

    @Test
    public void testVarianceIsTheDiagonalOfTheCovariance() {
        double[][] data = random(150, 7, 555L, -2.0);
        CovarianceAccumulator acc = new CovarianceAccumulator(7);
        acc.addAll(data);
        double[] c = acc.covariance();
        double[] v = acc.variance();
        for (int j = 0; j < 7; j++) assertEquals(c[j * 7 + j], v[j], 0.0);
    }

    @Test
    public void testPopulationVersusSampleCovariance() {
        double[][] data = random(64, 3, 31415L, 0.0);
        CovarianceAccumulator acc = new CovarianceAccumulator(3);
        acc.addAll(data);
        double[] sample = acc.covariance();
        double[] pop = acc.populationCovariance();
        for (int idx = 0; idx < sample.length; idx++)
            assertEquals(sample[idx] * 63.0 / 64.0, pop[idx], 1e-15);
    }

    @Test
    public void testCorrelationHasUnitDiagonalAndIsBounded() {
        double[][] data = random(120, 5, 8080L, 12.0);
        CovarianceAccumulator acc = new CovarianceAccumulator(5);
        acc.addAll(data);
        double[] r = acc.correlation();
        for (int j = 0; j < 5; j++) {
            assertEquals(1.0, r[j * 5 + j], 1e-14);
            for (int i = 0; i < 5; i++) assertTrue("|r| <= 1", Math.abs(r[j * 5 + i]) <= 1.0 + 1e-14);
        }
    }

    @Test
    public void testCorrelationOfAConstantFeature() {
        // feature 1 never varies -> undefined correlation, reported as 0 off-diagonal
        double[][] data = new double[40][2];
        long lcg = 777L;
        for (int i = 0; i < 40; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            data[i][0] = ((lcg >>> 11) / (double) (1L << 53));
            data[i][1] = 5.0;
        }
        CovarianceAccumulator acc = new CovarianceAccumulator(2);
        acc.addAll(data);
        double[] r = acc.correlation();
        assertEquals(1.0, r[0], 1e-14);
        assertEquals(1.0, r[3], 0.0);
        assertEquals(0.0, r[1], 0.0);
        assertEquals(0.0, r[2], 0.0);
    }

    @Test
    public void testSingleFeature() {
        double[][] data = { { 1.0 }, { 2.0 }, { 3.0 }, { 4.0 }, { 5.0 } };
        CovarianceAccumulator acc = new CovarianceAccumulator(1);
        acc.addAll(data);
        assertEquals(2.5, acc.covariance()[0], EPSILON);   // sample variance of 1..5
        assertEquals(3.0, acc.mean()[0], EPSILON);
    }

    @Test(expected = IllegalStateException.class)
    public void testCovarianceNeedsTwoSamples() {
        CovarianceAccumulator acc = new CovarianceAccumulator(2);
        acc.add(new double[] { 1.0, 2.0 });
        acc.covariance();
    }

    @Test(expected = IllegalStateException.class)
    public void testPopulationCovarianceNeedsOneSample() {
        new CovarianceAccumulator(2).populationCovariance();
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWrongSampleLength() {
        new CovarianceAccumulator(3).add(new double[] { 1.0, 2.0 });
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsMergeOfDifferentDimension() {
        new CovarianceAccumulator(3).merge(new CovarianceAccumulator(4));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsZeroDimension() {
        new CovarianceAccumulator(0);
    }
}
