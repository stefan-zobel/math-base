package math.linalg;

import static org.junit.Assert.*;

import org.junit.Test;

import math.probe.CovarianceAccumulator;

public class CovariancePCATest {

    private static final double EPSILON = 1e-12;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /**
     * Deterministic m x n row-major data with three dominant, clearly separated
     * directions on top of a large offset, so that the leading components are
     * well conditioned and the mean actually matters.
     */
    private static double[][] makeData(int m, int n, long seed) {
        double[][] d = new double[m][n];
        long lcg = seed;
        for (int i = 0; i < m; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f1 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f2 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double f3 = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            for (int j = 0; j < n; j++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double noise = (((lcg >>> 11) / (double) (1L << 53)) - 0.5) * 0.02;
                d[i][j] = 6.0 * f1 * Math.cos(j) + 3.0 * f2 * Math.sin(j)
                        + 1.5 * f3 * Math.cos(2 * j) + 10.0 + noise;
            }
        }
        return d;
    }

    private static double[] covarianceOf(double[][] data) {
        CovarianceAccumulator acc = new CovarianceAccumulator(data[0].length);
        acc.addAll(data);
        return acc.covariance();
    }

    /** Column-major diagonal matrix. */
    private static double[] diag(double[] d) {
        int n = d.length;
        double[] a = new double[n * n];
        for (int j = 0; j < n; j++) a[j * n + j] = d[j];
        return a;
    }

    private static double orthonormalityError(double[][] c) {
        double max = 0.0;
        for (int p = 0; p < c.length; p++)
            for (int q = 0; q < c.length; q++) {
                double dot = 0.0;
                for (int j = 0; j < c[p].length; j++) dot += c[p][j] * c[q][j];
                max = Math.max(max, Math.abs(dot - (p == q ? 1.0 : 0.0)));
            }
        return max;
    }

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testAgreesWithJacobiPCAOnTheSameData() {
        // the whole point: the covariance route must reproduce the SVD route for
        // the leading, well separated components
        int m = 400, n = 8, k = 3;
        double[][] data = makeData(m, n, 111L);

        CovarianceAccumulator acc = new CovarianceAccumulator(n);
        acc.addAll(data);
        CovariancePCA cpca = new CovariancePCA().fit(acc.covariance(), n, acc.mean(), k);

        JacobiPCA rpca = new JacobiPCA();
        rpca.pca(data, k);

        // both canonicalize signs the same way, so the directions must match outright
        for (int c = 0; c < k; c++)
            assertArrayEquals("component " + c,
                    rpca.getComponents()[c], cpca.getComponents()[c], 1e-8);

        double[] expected = rpca.getExplainedVariance();
        double[] actual = cpca.getExplainedVariance();
        for (int c = 0; c < k; c++)
            assertEquals(expected[c], actual[c], 1e-10 * expected[0]);
    }

    @Test
    public void testProjectionAgreesWithJacobiPCA() {
        int m = 200, n = 6, k = 2;
        double[][] data = makeData(m, n, 2222L);
        CovarianceAccumulator acc = new CovarianceAccumulator(n);
        acc.addAll(data);

        double[][] byCov = new CovariancePCA()
                .fit(acc.covariance(), n, acc.mean(), k).transform(data);
        double[][] byRaw = new JacobiPCA().pca(data, k);

        for (int i = 0; i < m; i++)
            assertArrayEquals("sample " + i, byRaw[i], byCov[i], 1e-7);
    }

    @Test
    public void testDiagonalCovarianceGivesSortedEigenvaluesAndUnitVectors() {
        double[] c = diag(new double[] { 2.0, 9.0, 5.0 });
        CovariancePCA pca = new CovariancePCA().fit(c, 3, 3);
        assertArrayEquals(new double[] { 9.0, 5.0, 2.0 }, pca.getEigenvalues(), EPSILON);
        // components are the unit vectors of the corresponding axes, sign canonicalized
        assertArrayEquals(new double[] { 0.0, 1.0, 0.0 }, pca.getComponents()[0], EPSILON);
        assertArrayEquals(new double[] { 0.0, 0.0, 1.0 }, pca.getComponents()[1], EPSILON);
        assertArrayEquals(new double[] { 1.0, 0.0, 0.0 }, pca.getComponents()[2], EPSILON);
        assertEquals(16.0, pca.getTotalVariance(), EPSILON);
        assertArrayEquals(new double[] { 9.0 / 16.0, 5.0 / 16.0, 2.0 / 16.0 },
                pca.getExplainedVarianceRatio(), EPSILON);
    }

    @Test
    public void testIdentityCovariance() {
        int n = 5;
        double[] one = new double[n];
        for (int j = 0; j < n; j++) one[j] = 1.0;
        CovariancePCA pca = new CovariancePCA().fit(diag(one), n, n);
        for (int j = 0; j < n; j++) assertEquals(1.0, pca.getEigenvalues()[j], EPSILON);
        assertEquals(0.0, orthonormalityError(pca.getComponents()), EPSILON);
        assertEquals((double) n, pca.getTotalVariance(), EPSILON);
    }

    @Test
    public void testComponentsAreOrthonormalAndRatiosSumToOne() {
        int n = 7;
        double[][] data = makeData(250, n, 333L);
        CovariancePCA pca = new CovariancePCA().fit(covarianceOf(data), n, n);
        assertEquals(0.0, orthonormalityError(pca.getComponents()), 1e-13);

        double sum = 0.0;
        for (int c = 0; c < n; c++) sum += pca.getExplainedVarianceRatio()[c];
        assertEquals(1.0, sum, 1e-12);

        double traceFromEigenvalues = 0.0;
        for (int c = 0; c < n; c++) traceFromEigenvalues += pca.getEigenvalues()[c];
        assertEquals(pca.getTotalVariance(), traceFromEigenvalues, 1e-12 * pca.getTotalVariance());
    }

    @Test
    public void testEigenvaluesAreDescending() {
        int n = 6;
        double[][] data = makeData(180, n, 4444L);
        CovariancePCA pca = new CovariancePCA().fit(covarianceOf(data), n, n);
        for (int c = 1; c < n; c++)
            assertTrue("eigenvalues must be descending",
                    pca.getEigenvalues()[c - 1] >= pca.getEigenvalues()[c]);
    }

    @Test
    public void testTransformCentersWithTheSuppliedMean() {
        int n = 4;
        double[][] data = makeData(100, n, 5555L);
        CovarianceAccumulator acc = new CovarianceAccumulator(n);
        acc.addAll(data);
        double[] mean = acc.mean();

        CovariancePCA withMean = new CovariancePCA().fit(acc.covariance(), n, mean, 2);
        CovariancePCA withoutMean = new CovariancePCA().fit(acc.covariance(), n, 2);

        double[][] centered = new double[100][n];
        for (int i = 0; i < 100; i++)
            for (int j = 0; j < n; j++) centered[i][j] = data[i][j] - mean[j];

        // centering inside transform must equal centering the samples beforehand
        double[][] a = withMean.transform(data);
        double[][] b = withoutMean.transform(centered);
        for (int i = 0; i < 100; i++) assertArrayEquals(b[i], a[i], 1e-11);

        // the projections of the centered data must average out to zero
        for (int c = 0; c < 2; c++) {
            double s = 0.0;
            for (int i = 0; i < 100; i++) s += a[i][c];
            assertEquals(0.0, s / 100.0, 1e-12);
        }
    }

    @Test
    public void testWhitenedDataHasUnitVariance() {
        int m = 300, n = 5, k = 3;
        double[][] data = makeData(m, n, 6666L);
        CovarianceAccumulator acc = new CovarianceAccumulator(n);
        acc.addAll(data);
        double[][] w = new CovariancePCA()
                .fit(acc.covariance(), n, acc.mean(), k).whiten(data);

        for (int c = 0; c < k; c++) {
            double s2 = 0.0;
            for (int i = 0; i < m; i++) s2 += w[i][c] * w[i][c];
            assertEquals("component " + c, 1.0, s2 / (m - 1), 1e-10);
        }
        // and the whitened components stay uncorrelated
        for (int c = 1; c < k; c++) {
            double dot = 0.0;
            for (int i = 0; i < m; i++) dot += w[i][0] * w[i][c];
            assertEquals(0.0, dot / (m - 1), 1e-10);
        }
    }

    @Test
    public void testWhiteningIsRefusedForANonPositiveEigenvalue() {
        double[] c = diag(new double[] { 4.0, 0.0 });
        CovariancePCA pca = new CovariancePCA().fit(c, 2, 2);
        try {
            pca.whiten(new double[][] { { 1.0, 1.0 } });
            fail("whitening a singular covariance must be refused");
        } catch (IllegalStateException expected) {
            // the message must name the offending component
            assertTrue(expected.getMessage().contains("eigenvalue"));
        }
        // keeping only the healthy component is still fine
        double[][] w = new CovariancePCA().fit(c, 2, 1).whiten(new double[][] { { 2.0, 1.0 } });
        assertEquals(1.0, w[0][0], EPSILON);
    }

    @Test
    public void testDetectsAnIndefiniteInput() {
        // a matrix that is symmetric but not positive semidefinite: the SVD would
        // report |lambda| and hide this, the eigensolver does not
        double[] c = { 1.0, 2.0, 2.0, 1.0 };   // eigenvalues 3 and -1
        CovariancePCA pca = new CovariancePCA().fit(c, 2, 2);
        assertArrayEquals(new double[] { 3.0, -1.0 }, pca.getEigenvalues(), EPSILON);
        assertEquals(-1.0, pca.getMinEigenvalue(), EPSILON);
        assertEquals(2.0, pca.getTotalVariance(), EPSILON);
    }

    @Test
    public void testRankDeficientCovarianceKeepsACompleteBasis() {
        // feature 3 is the sum of features 0 and 1 -> exactly singular
        int m = 150, n = 4;
        double[][] base = makeData(m, n - 1, 7777L);
        double[][] data = new double[m][n];
        for (int i = 0; i < m; i++) {
            System.arraycopy(base[i], 0, data[i], 0, n - 1);
            data[i][n - 1] = base[i][0] + base[i][1];
        }
        CovariancePCA pca = new CovariancePCA().fit(covarianceOf(data), n, n);
        double lambdaMax = pca.getEigenvalues()[0];
        assertEquals("the null direction must show up as a vanishing eigenvalue",
                0.0, pca.getMinEigenvalue(), 1e-12 * lambdaMax);
        // unlike the SVD of a singular matrix, the basis stays complete
        assertEquals(0.0, orthonormalityError(pca.getComponents()), 1e-13);
    }

    @Test
    public void testInputIsNotModified() {
        double[][] data = makeData(60, 5, 8888L);
        double[] cov = covarianceOf(data);
        double[] copy = cov.clone();
        new CovariancePCA().fit(cov, 5, 3);
        assertArrayEquals("fit() must not touch its input", copy, cov, 0.0);
    }

    @Test
    public void testSingleFeature() {
        CovariancePCA pca = new CovariancePCA().fit(new double[] { 7.0 }, 1, 1);
        assertArrayEquals(new double[] { 7.0 }, pca.getEigenvalues(), EPSILON);
        assertArrayEquals(new double[] { 1.0 }, pca.getComponents()[0], EPSILON);
        assertArrayEquals(new double[] { 1.0 }, pca.getExplainedVarianceRatio(), EPSILON);
    }

    @Test
    public void testMeanIsCopiedDefensively() {
        double[] mean = { 1.0, 2.0 };
        CovariancePCA pca = new CovariancePCA().fit(diag(new double[] { 3.0, 1.0 }), 2, mean, 1);
        mean[0] = 99.0;
        assertArrayEquals(new double[] { 1.0, 2.0 }, pca.getMean(), 0.0);
        assertNull(new CovariancePCA().fit(diag(new double[] { 3.0, 1.0 }), 2, 1).getMean());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsAsymmetricInput() {
        new CovariancePCA().fit(new double[] { 1.0, 0.5, -0.5, 1.0 }, 2, 2);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWrongLength() {
        new CovariancePCA().fit(new double[5], 2, 2);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsTooManyComponents() {
        new CovariancePCA().fit(diag(new double[] { 1.0, 1.0 }), 2, 3);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsWrongMeanLength() {
        new CovariancePCA().fit(diag(new double[] { 1.0, 1.0 }), 2, new double[3], 1);
    }

    @Test(expected = IllegalStateException.class)
    public void testTransformBeforeFit() {
        new CovariancePCA().transform(new double[][] { { 1.0 } });
    }
}
