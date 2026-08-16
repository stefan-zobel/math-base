package math.linalg;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.junit.Test;

/**
 * Tests for {@link TruncatedPCA}, which computes only the leading components and therefore has to
 * agree with {@link JacobiPCA} wherever the data determines those components at all.
 */
public class TruncatedPCATest {

    private static final int ROWS = 400;
    private static final int COLS = 120;

    /**
     * Data with a decaying spectrum: 20 latent directions with weights 1/(k+1) plus noise. Every
     * leading component is well separated from the next, so both methods have to find the same one.
     */
    private static double[][] decayingSpectrum(int m, int n) {
        Random random = new Random(23);
        int factors = Math.min(20, n);
        double[][] directions = new double[factors][n];
        for (int k = 0; k < factors; k++) {
            for (int j = 0; j < n; j++) directions[k][j] = random.nextGaussian();
        }
        double[][] data = new double[m][n];
        for (int i = 0; i < m; i++) {
            double[] coefficients = new double[factors];
            for (int k = 0; k < factors; k++) coefficients[k] = random.nextGaussian() / (k + 1.0);
            for (int j = 0; j < n; j++) {
                double value = random.nextGaussian() * 0.3;
                for (int k = 0; k < factors; k++) value += coefficients[k] * directions[k][j];
                data[i][j] = value;
            }
        }
        return data;
    }

    /** One dominant direction (the all-ones vector) and nothing but noise behind it. */
    private static double[][] oneDominantDirection(int m, int n) {
        Random random = new Random(7);
        double[][] data = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) data[i][j] = random.nextGaussian() + (i % 10) * 3.0;
        }
        return data;
    }

    /** Largest deviation in component {@code c}, relative to its spread and up to a sign flip. */
    private static double deviation(double[][] expected, double[][] actual, int c) {
        double scale = 0.0, plus = 0.0, minus = 0.0;
        for (int i = 0; i < expected.length; i++) scale = Math.max(scale, Math.abs(expected[i][c]));
        for (int i = 0; i < expected.length; i++) {
            plus = Math.max(plus, Math.abs(expected[i][c] - actual[i][c]));
            minus = Math.max(minus, Math.abs(expected[i][c] + actual[i][c]));
        }
        return Math.min(plus, minus) / scale;
    }

    @Test
    public void agreesWithTheExactPcaWhereTheSpectrumIsSeparated() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        double[][] exact = new JacobiPCA().pca(data, 2);

        TruncatedPCA truncated = new TruncatedPCA();
        double[][] approximate = truncated.pca(data, 2);

        assertTrue("should have converged", truncated.converged());
        assertTrue("component 1: " + deviation(exact, approximate, 0), deviation(exact, approximate, 0) < 1e-4);
        assertTrue("component 2: " + deviation(exact, approximate, 1), deviation(exact, approximate, 1) < 1e-4);
    }

    @Test
    public void agreesWithTheExactPcaOnTheSingularValues() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        JacobiPCA exact = new JacobiPCA();
        exact.pca(data, 3);
        TruncatedPCA truncated = new TruncatedPCA();
        truncated.pca(data, 3);

        double[] expected = exact.getSingularValues();
        double[] actual = truncated.getSingularValues();
        for (int k = 0; k < expected.length; k++) {
            assertEquals("singular value " + k, expected[k], actual[k], 1e-6 * expected[0]);
        }
    }

    @Test
    public void findsTheLeadingComponentWhenTheRestIsDegenerate() {
        double[][] data = oneDominantDirection(ROWS, COLS);
        double[][] exact = new JacobiPCA().pca(data, 2);

        TruncatedPCA truncated = new TruncatedPCA();
        double[][] approximate = truncated.pca(data, 2);

        assertTrue("component 1: " + deviation(exact, approximate, 0), deviation(exact, approximate, 0) < 1e-8);
        // component 2 is deliberately not asserted: behind the dominant direction the spectrum is
        // degenerate, so the second direction is not determined by the data and the exact solver
        // picks an equally arbitrary one. Its spread is smaller than the first by two orders of
        // magnitude, which the next assertion pins down instead.
        double spreadFirst = spread(approximate, 0);
        double spreadSecond = spread(approximate, 1);
        assertTrue("the dominant direction must carry far more spread",
                spreadFirst > 50.0 * spreadSecond);
    }

    private static double spread(double[][] projected, int c) {
        double mean = 0.0;
        for (double[] row : projected) mean += row[c];
        mean /= projected.length;
        double variance = 0.0;
        for (double[] row : projected) {
            double d = row[c] - mean;
            variance += d * d;
        }
        return Math.sqrt(variance / projected.length);
    }

    @Test
    public void reportsNonConvergenceOnAStructurelessSpectrum() {
        // pure noise has no dominant direction: the Ritz values are distinct but separated by so
        // little that they keep moving for very many iterations. This is the case the fallback in
        // the caller exists for, and the result really is unusable here.
        Random random = new Random(5);
        double[][] noise = new double[ROWS][COLS];
        for (int i = 0; i < ROWS; i++) {
            for (int j = 0; j < COLS; j++) noise[i][j] = random.nextGaussian();
        }
        TruncatedPCA truncated = new TruncatedPCA();
        truncated.pca(noise, 2);
        assertFalse("structureless data must not be reported as converged", truncated.converged());
    }

    @Test
    public void theProjectionMatchesTheComponentsItReports() {
        double[][] data = decayingSpectrum(200, 60);
        TruncatedPCA truncated = new TruncatedPCA();
        double[][] projected = truncated.pca(data, 2);
        double[] mean = truncated.getMean();
        double[][] components = truncated.getComponents();

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < 2; k++) {
                double dot = 0.0;
                for (int j = 0; j < data[i].length; j++) dot += (data[i][j] - mean[j]) * components[k][j];
                assertEquals(dot, projected[i][k], 1e-9);
            }
        }
    }

    @Test
    public void componentsAreOrthonormalAndSingularValuesDescend() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        TruncatedPCA truncated = new TruncatedPCA();
        truncated.pca(data, 3);

        double[][] c = truncated.getComponents();
        for (int p = 0; p < c.length; p++) {
            for (int q = 0; q < c.length; q++) {
                double dot = 0.0;
                for (int j = 0; j < c[p].length; j++) dot += c[p][j] * c[q][j];
                assertEquals(p == q ? 1.0 : 0.0, dot, 1e-9);
            }
        }
        double[] sv = truncated.getSingularValues();
        for (int k = 0; k + 1 < sv.length; k++) {
            assertTrue("singular values must descend", sv[k] >= sv[k + 1]);
        }
    }

    @Test
    public void canonicalizesTheSignOfEveryComponent() {
        TruncatedPCA truncated = new TruncatedPCA();
        truncated.pca(decayingSpectrum(ROWS, COLS), 3);
        for (double[] component : truncated.getComponents()) {
            int argmax = 0;
            for (int j = 0; j < component.length; j++) {
                if (Math.abs(component[j]) > Math.abs(component[argmax])) argmax = j;
            }
            assertTrue("largest magnitude entry must be positive", component[argmax] > 0.0);
        }
    }

    @Test
    public void isDeterministic() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        double[][] first = new TruncatedPCA().pca(data, 2);
        double[][] second = new TruncatedPCA().pca(data, 2);
        for (int i = 0; i < first.length; i++) {
            assertArrayEquals(first[i], second[i], 0.0);
        }
    }

    @Test
    public void reportsNonConvergenceWhenTheIterationLimitIsTooLow() {
        // one iteration cannot produce a second estimate to compare against, so the stability test
        // never runs - this is the branch that makes the caller fall back to the exact PCA
        TruncatedPCA truncated = new TruncatedPCA(6, 1e-15, 1, 20260815L);
        truncated.pca(decayingSpectrum(ROWS, COLS), 2);
        assertFalse(truncated.converged());
        assertEquals(1, truncated.getIterations());
    }

    @Test
    public void convergesWithoutOversampling() {
        TruncatedPCA truncated = new TruncatedPCA(0, 1e-9, 60, 20260815L);
        double[][] data = decayingSpectrum(ROWS, COLS);
        double[][] approximate = truncated.pca(data, 2);
        assertTrue("should have converged", truncated.converged());
        assertTrue(deviation(new JacobiPCA().pca(data, 2), approximate, 0) < 1e-4);
    }

    @Test
    public void fixedIterationsRunsExactlyTheRequestedNumber() {
        // no stability test decides anything here: the count alone ends the loop, which is the point
        // of the factory - see its javadoc for why that is the right design for many components
        TruncatedPCA truncated = TruncatedPCA.fixedIterations(10, 6);
        truncated.pca(decayingSpectrum(ROWS, COLS), 2);
        assertEquals(6, truncated.getIterations());
    }

    @Test
    public void fixedIterationsIsDeterministic() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        double[][] first = TruncatedPCA.fixedIterations(10, 6).pca(data, 2);
        double[][] second = TruncatedPCA.fixedIterations(10, 6).pca(data, 2);
        for (int i = 0; i < first.length; i++) {
            assertArrayEquals("row " + i, first[i], second[i], 0.0);
        }
    }

    @Test
    public void fixedIterationsResolvesTheLeadingComponents() {
        double[][] data = decayingSpectrum(ROWS, COLS);
        double[][] approximate = TruncatedPCA.fixedIterations(10, 6).pca(data, 2);
        double[][] exact = new JacobiPCA().pca(data, 2);
        assertTrue("component 1: " + deviation(exact, approximate, 0), deviation(exact, approximate, 0) < 1e-4);
        assertTrue("component 2: " + deviation(exact, approximate, 1), deviation(exact, approximate, 1) < 1e-4);
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsAnEmptyMatrix() {
        new TruncatedPCA().pca(new double[0][0], 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsARaggedMatrix() {
        new TruncatedPCA().pca(new double[][] { { 1.0, 2.0 }, { 3.0 } }, 1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void rejectsTooManyComponents() {
        new TruncatedPCA().pca(new double[][] { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } }, 3);
    }
}
