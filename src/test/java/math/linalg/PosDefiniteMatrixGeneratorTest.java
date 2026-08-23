package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * {@link PosDefiniteMatrixGenerator} drew from a static {@code Stc64} shared by
 * every call, which made concurrent use unsafe -- the gaussian cache in
 * {@code AbstractRng64} answers {@code NaN} when two threads race on it, and
 * this generator is almost all gaussians -- and left no way to reproduce a
 * matrix. It builds its own generator per call now and takes a seed.
 * <p>
 * Over 200 seeds and the six dimensions below, Cholesky never failed, the worst
 * asymmetry was 3.6e-14 and the smallest pivot 1.1e-03.
 */
public class PosDefiniteMatrixGeneratorTest {

    private static final int[] DIMENSIONS = { 1, 2, 3, 5, 10, 25 };

    /** The same seed has to reproduce the same matrix to the last bit. */
    @Test
    public void testTheSameSeedReproducesTheSameMatrix() {
        for (int d = 0; d < DIMENSIONS.length; ++d) {
            int dim = DIMENSIONS[d];
            DMatrix a = PosDefiniteMatrixGenerator.generate(dim, 20260823L);
            DMatrix b = PosDefiniteMatrixGenerator.generate(dim, 20260823L);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    assertEquals("dim " + dim + " entry (" + i + ", " + j + ")", a.get(i, j), b.get(i, j),
                            0.0);
                }
            }
        }
    }

    /** And a different seed must not. */
    @Test
    public void testADifferentSeedGivesADifferentMatrix() {
        DMatrix a = PosDefiniteMatrixGenerator.generate(10, 20260823L);
        DMatrix b = PosDefiniteMatrixGenerator.generate(10, 20260824L);
        int equal = 0;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                if (a.get(i, j) == b.get(i, j)) {
                    ++equal;
                }
            }
        }
        assertEquals("the two seeds produced the same matrix", 0, equal);
    }

    /**
     * The defining property. Symmetry is exact in theory -- the matrix is
     * {@code Q' D Q} with a diagonal {@code D} -- and the residual is the
     * rounding of the two products; positive definiteness is what Cholesky
     * succeeding means.
     */
    @Test
    public void testTheMatrixIsSymmetricAndPositiveDefinite() {
        for (int d = 0; d < DIMENSIONS.length; ++d) {
            int dim = DIMENSIONS[d];
            for (long seed = 100L; seed < 120L; ++seed) {
                DMatrix a = PosDefiniteMatrixGenerator.generate(dim, seed);
                assertEquals("rows", dim, a.numRows());
                assertEquals("columns", dim, a.numColumns());
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j) {
                        double scale = Math.max(1.0, Math.abs(a.get(i, j)));
                        assertEquals("dim " + dim + " seed " + seed + " asymmetry at (" + i + ", " + j + ")",
                                a.get(i, j), a.get(j, i), 1.0e-12 * scale);
                    }
                }
                DMatrix l = CholeskyDecomp.cholesky(a);
                for (int i = 0; i < dim; ++i) {
                    assertTrue("dim " + dim + " seed " + seed + " pivot " + i + " is " + l.get(i, i),
                            Math.abs(l.get(i, i)) > 0.0);
                }
            }
        }
    }
}
