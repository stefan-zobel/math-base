package math.linalg;

import static org.junit.Assert.*;

import java.io.OutputStream;
import java.io.PrintStream;

import org.junit.Test;

public class JacobiPCATest {

    private static final double EPSILON = 1e-10;

    // -----------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------

    /**
     * Deterministic row-major data with {@code latent} dominant directions of clearly
     * separated strength plus a little noise, so that the leading components are
     * well conditioned and their order is unambiguous.
     */
    private static double[][] data(int m, int n, int latent, long seed) {
        double[][] d = new double[m][n];
        long lcg = seed;
        double[] f = new double[latent];
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < latent; k++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                f[k] = ((lcg >>> 11) / (double) (1L << 53)) - 0.5;
            }
            for (int j = 0; j < n; j++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                double v = (((lcg >>> 11) / (double) (1L << 53)) - 0.5) * 0.01;
                for (int k = 0; k < latent; k++) {
                    v += (8.0 - 2.5 * k) * f[k] * Math.cos((k + 1) * (j + 1) * 0.7);
                }
                d[i][j] = v;
            }
        }
        return d;
    }

    private static double[] columnMeans(double[][] d) {
        int m = d.length, n = d[0].length;
        double[] mean = new double[n];
        for (double[] row : d) {
            for (int j = 0; j < n; j++) {
                mean[j] += row[j];
            }
        }
        for (int j = 0; j < n; j++) {
            mean[j] /= m;
        }
        return mean;
    }

    /** Sample covariance matrix (divisor m-1) as an n x n column-major array. */
    private static double[] covariance(double[][] d) {
        int m = d.length, n = d[0].length;
        double[] mean = columnMeans(d);
        double[] c = new double[n * n];
        for (int p = 0; p < n; p++) {
            for (int q = p; q < n; q++) {
                double sum = 0.0;
                for (int i = 0; i < m; i++) {
                    sum += (d[i][p] - mean[p]) * (d[i][q] - mean[q]);
                }
                double cov = sum / (m - 1);
                c[q * n + p] = cov;
                c[p * n + q] = cov;
            }
        }
        return c;
    }

    /** Flips the sign so that the entry of largest magnitude is positive. */
    private static double[] canonicalize(double[] v) {
        double[] out = v.clone();
        int argmax = 0;
        double best = -1.0;
        for (int j = 0; j < out.length; j++) {
            double abs = Math.abs(out[j]);
            if (abs > best) {
                best = abs;
                argmax = j;
            }
        }
        if (out[argmax] < 0.0) {
            for (int j = 0; j < out.length; j++) {
                out[j] = -out[j];
            }
        }
        return out;
    }

    /** max |C C^T - I| for components C ([k][n]). */
    private static double componentOrthoError(double[][] c) {
        double max = 0.0;
        for (int p = 0; p < c.length; p++) {
            for (int q = 0; q < c.length; q++) {
                double dot = 0.0;
                for (int j = 0; j < c[p].length; j++) {
                    dot += c[p][j] * c[q][j];
                }
                max = Math.max(max, Math.abs(dot - (p == q ? 1.0 : 0.0)));
            }
        }
        return max;
    }

    private static double[] column(double[][] a, int k) {
        double[] col = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            col[i] = a[i][k];
        }
        return col;
    }

    private static double mean(double[] x) {
        double s = 0.0;
        for (double v : x) {
            s += v;
        }
        return s / x.length;
    }

    /** Sample variance with divisor n-1. */
    private static double variance(double[] x) {
        double mu = mean(x);
        double s = 0.0;
        for (double v : x) {
            s += (v - mu) * (v - mu);
        }
        return s / (x.length - 1);
    }

    /**
     * Runs the checks that must hold for every shape: orthonormal components,
     * descending singular values, canonical signs, centered and uncorrelated
     * projections, and agreement with the eigen decomposition of the sample
     * covariance matrix.
     */
    private static void assertGoodPca(double[][] d, int noDims, double tol) {
        int m = d.length, n = d[0].length;
        JacobiPCA pca = new JacobiPCA();
        double[][] y = pca.pca(d, noDims);

        assertEquals(m, y.length);
        assertEquals(noDims, y[0].length);
        assertArrayEquals("per-feature mean", columnMeans(d), pca.getMean(), tol);
        assertEquals("components orthonormal", 0.0, componentOrthoError(pca.getComponents()), tol);

        double[] sv = pca.getSingularValues();
        assertEquals(noDims, sv.length);
        for (int k = 1; k < noDims; k++) {
            assertTrue("singular values must be descending", sv[k - 1] >= sv[k]);
        }
        for (double[] comp : pca.getComponents()) {
            assertArrayEquals("sign must be canonical", canonicalize(comp), comp, 0.0);
        }

        // projections are centered, and their variance is the explained variance
        double[] ev = pca.getExplainedVariance();
        for (int k = 0; k < noDims; k++) {
            double[] yk = column(y, k);
            assertEquals("projection must be centered", 0.0, mean(yk), tol);
            assertEquals("variance of component " + k, ev[k], variance(yk), tol);
        }
        // different components are uncorrelated
        for (int p = 0; p < noDims; p++) {
            for (int q = p + 1; q < noDims; q++) {
                double[] yp = column(y, p);
                double[] yq = column(y, q);
                double dot = 0.0;
                for (int i = 0; i < m; i++) {
                    dot += yp[i] * yq[i];
                }
                assertEquals("components must be uncorrelated", 0.0, dot, tol);
            }
        }

        // cross-check against the eigen decomposition of the covariance matrix
        SymmetricJacobiEigen.Result eig = new SymmetricJacobiEigen().decompose(covariance(d), n);
        for (int k = 0; k < noDims; k++) {
            assertEquals("explained variance " + k, eig.lambda[k], ev[k], tol);
            assertArrayEquals("principal direction " + k,
                    canonicalize(eig.eigenvector(k)), pca.getComponents()[k], tol);
        }
    }

    // -----------------------------------------------------------------
    // tests
    // -----------------------------------------------------------------

    @Test
    public void testTallCase() {
        // m > n -> QR preconditioned path
        assertGoodPca(data(200, 8, 3, 111L), 3, 1e-8);
    }

    @Test
    public void testSquareCase() {
        // m == n -> direct SVD path
        assertGoodPca(data(9, 9, 3, 222L), 3, 1e-9);
    }

    @Test
    public void testWideCase() {
        // m < n -> the transposed path; the centered data has rank <= m-1 = 5
        assertGoodPca(data(6, 40, 3, 333L), 3, 1e-8);
    }

    @Test
    public void testSmallestPossibleShapes() {
        assertGoodPca(data(2, 1, 1, 444L), 1, 1e-9);
        assertGoodPca(data(3, 2, 2, 666L), 2, 1e-9);
    }

    @Test
    public void testSingleSampleIsDegenerate() {
        // A single sample has no variance at all: after centering the data is the
        // zero matrix, so there is no direction to find. The current behavior is to
        // return a zero component (and thus a zero projection) rather than to reject
        // the input - pinned here so a future change to that contract is noticed.
        double[][] d = { { 1.0, 2.0, 3.0 } };
        JacobiPCA pca = new JacobiPCA();
        double[][] y = pca.pca(d, 1);

        assertArrayEquals(new double[] { 1.0, 2.0, 3.0 }, pca.getMean(), EPSILON);
        assertEquals(0.0, pca.getSingularValues()[0], 0.0);
        assertArrayEquals(new double[] { 0.0, 0.0, 0.0 }, pca.getComponents()[0], 0.0);
        assertArrayEquals(new double[] { 0.0 }, y[0], 0.0);
    }

    @Test
    public void testExactlyCollinearData() {
        // points on the line y = 2x + 1 -> one direction carries all the variance
        int m = 5;
        double[][] d = new double[m][2];
        for (int i = 0; i < m; i++) {
            d[i][0] = i;
            d[i][1] = 2.0 * i + 1.0;
        }
        JacobiPCA pca = new JacobiPCA();
        pca.pca(d, 2);

        double[] sv = pca.getSingularValues();
        assertTrue("first singular value must dominate", sv[0] > 1.0);
        assertEquals("second direction carries no variance", 0.0, sv[1], 1e-13 * sv[0]);

        double inv = 1.0 / Math.sqrt(5.0);
        assertArrayEquals(new double[] { inv, 2.0 * inv }, pca.getComponents()[0], 1e-12);
        assertArrayEquals(new double[] { 2.0, 5.0 }, pca.getMean(), EPSILON);
    }

    @Test
    public void testAxisAlignedDataKeepsTheDominantAxis() {
        // feature 0 varies a lot, feature 1 barely -> first component is the x axis
        int m = 40;
        double[][] d = new double[m][2];
        for (int i = 0; i < m; i++) {
            d[i][0] = 100.0 * (i - m / 2.0);
            d[i][1] = ((i % 2 == 0) ? 1.0 : -1.0) * 0.001;
        }
        JacobiPCA pca = new JacobiPCA();
        double[][] y = pca.pca(d, 1);
        // the residual tilt away from the x axis is real (the two features are not
        // exactly uncorrelated), it is just tiny - hence 1e-6 and not 1e-12 here
        assertArrayEquals(new double[] { 1.0, 0.0 }, pca.getComponents()[0], 1e-6);
        assertEquals(m, y.length);
        assertEquals(1, y[0].length);
    }

    @Test
    public void testExplainedVarianceSumsToTotalVariance() {
        // keeping all components must preserve the total variance of the data
        double[][] d = data(50, 6, 4, 777L);
        JacobiPCA pca = new JacobiPCA();
        pca.pca(d, 6);

        double totalVariance = 0.0;
        for (int j = 0; j < 6; j++) {
            totalVariance += variance(column(d, j));
        }
        double explained = 0.0;
        for (double e : pca.getExplainedVariance()) {
            explained += e;
        }
        assertEquals(totalVariance, explained, 1e-9 * totalVariance);
    }

    @Test
    public void testTransformReproducesTheTrainingProjection() {
        double[][] d = data(60, 5, 3, 888L);
        JacobiPCA pca = new JacobiPCA();
        double[][] fitted = pca.pca(d, 3);
        double[][] transformed = pca.transform(d);
        assertEquals(fitted.length, transformed.length);
        for (int i = 0; i < fitted.length; i++) {
            assertArrayEquals(fitted[i], transformed[i], 1e-12);
        }
    }

    @Test
    public void testTransformOfUnseenSamples() {
        double[][] d = data(60, 5, 3, 999L);
        JacobiPCA pca = new JacobiPCA();
        pca.pca(d, 2);

        // a sample equal to the mean projects onto the origin
        double[][] atMean = { pca.getMean().clone() };
        assertArrayEquals(new double[] { 0.0, 0.0 }, pca.transform(atMean)[0], 1e-12);

        // mean + component_0 projects onto (1, 0)
        double[] shifted = pca.getMean().clone();
        for (int j = 0; j < shifted.length; j++) {
            shifted[j] += pca.getComponents()[0][j];
        }
        assertArrayEquals(new double[] { 1.0, 0.0 }, pca.transform(new double[][] { shifted })[0], 1e-12);
    }

    @Test
    public void testRefitOnDifferentData() {
        JacobiPCA pca = new JacobiPCA();
        pca.pca(data(30, 4, 2, 1234L), 2);
        double[][] second = data(20, 6, 3, 5678L);
        double[][] y = pca.pca(second, 3);
        assertEquals(20, y.length);
        assertEquals(3, y[0].length);
        assertEquals(6, pca.getMean().length);
        assertEquals(3, pca.getComponents().length);
        assertEquals(6, pca.getComponents()[0].length);
        assertArrayEquals(columnMeans(second), pca.getMean(), 1e-12);
    }

    @Test
    public void testCustomSvdInstanceIsUsed() {
        // the SVD passed in must actually drive the decomposition: a solver that is
        // not allowed to converge makes the PCA fail instead of returning garbage
        double[][] d = data(20, 6, 3, 4321L);
        JacobiPCA pca = new JacobiPCA(new FlatParallelJacobiSVD(1e-15, 1, 1L << 16));
        // the solver writes its non-convergence warning to System.err
        // - expected here, but it has no business in the build log
        PrintStream err = System.err;
        System.setErr(new PrintStream(new OutputStream() {
            @Override
            public void write(int b) {
                // discard
            }
        }));
        try {
            pca.pca(d, 2);
            fail("expected IllegalStateException for a non-converged SVD");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("did not converge"));
        } finally {
            System.setErr(err);
        }
    }

    @Test(expected = IllegalStateException.class)
    public void testTransformBeforeFitFails() {
        new JacobiPCA().transform(new double[][] { { 1.0, 2.0 } });
    }

    @Test(expected = IllegalArgumentException.class)
    public void testTransformWithWrongSampleLengthFails() {
        JacobiPCA pca = new JacobiPCA();
        pca.pca(data(10, 4, 2, 1L), 2);
        pca.transform(new double[][] { { 1.0, 2.0 } });
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsNullMatrix() {
        new JacobiPCA().pca(null, 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsEmptyMatrix() {
        new JacobiPCA().pca(new double[0][0], 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsRowsWithoutFeatures() {
        new JacobiPCA().pca(new double[3][0], 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsRaggedMatrix() {
        new JacobiPCA().pca(new double[][] { { 1.0, 2.0 }, { 3.0 } }, 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsTooFewDimensions() {
        new JacobiPCA().pca(data(10, 4, 2, 2L), 0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testRejectsTooManyDimensions() {
        // no_dims must not exceed min(m, n) = 4
        new JacobiPCA().pca(data(10, 4, 2, 3L), 5);
    }
}
