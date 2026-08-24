package math.trmm;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.Test;

import math.gemm.Trans;

/**
 * {@link Dtrmm} against a naive dense triangular multiply, and the two
 * implementations behind it against each other. The public entry point picks
 * the blocked Lehn path above a size threshold and the netlib translation
 * below it, so a test that only goes through {@code dtrmm} cannot compare the
 * two at one size. This class lives in the package and reaches both directly.
 */
public class DtrmmTest {

    /** Shapes crossing the LEFT threshold at 191 and the RIGHT one at 281. */
    private static final int[][] SHAPES = { { 1, 1 }, { 2, 3 }, { 4, 4 }, { 7, 5 }, { 16, 16 }, { 31, 17 },
            { 64, 64 }, { 191, 191 }, { 200, 130 }, { 281, 281 } };

    private static final double[] ALPHAS = { 1.0, -0.75 };

    private static final Trans[] TRANS = { Trans.NO_TRANS, Trans.TRANS };

    /**
     * Norm-wise per configuration, {@code max|got - want| / max|want|}. The
     * entry-wise relative error is the wrong measure here: entries whose true
     * value nearly cancels reach 5.3e-11 on a correct implementation. Measured
     * over the 320 configurations below, the worst norm-wise residual is
     * 1.633e-15.
     */
    private static final double TOL = 1e-13;

    // ----- the cross-validation ---------------------------------------------

    @Test
    public void theProductAgreesWithANaiveTriangularMultiply() {
        int configs = 0;
        double worst = 0.0;
        String worstAt = "";
        for (Side side : Side.values()) {
            for (UpLo uplo : UpLo.values()) {
                for (int t = 0; t < TRANS.length; ++t) {
                    for (Diag diag : Diag.values()) {
                        for (int ai = 0; ai < ALPHAS.length; ++ai) {
                            for (int s = 0; s < SHAPES.length; ++s) {
                                int m = SHAPES[s][0];
                                int n = SHAPES[s][1];
                                int k = (side == Side.LEFT) ? m : n;
                                double[] a = triangular(k, 100L + s);
                                double[] b = dense(m * n, 200L + s);
                                double[] want = naive(side, uplo, TRANS[t], diag, m, n, ALPHAS[ai], a, k, b);
                                double[] got = b.clone();
                                Dtrmm.dtrmm(side, uplo, TRANS[t], diag, m, n, ALPHAS[ai], a, 0, k, got, 0, m);
                                ++configs;
                                double rel = residual(got, want);
                                if (rel > worst) {
                                    worst = rel;
                                    worstAt = m + "x" + n + " " + side + "/" + uplo + "/" + TRANS[t] + "/" + diag
                                            + " alpha=" + ALPHAS[ai];
                                }
                            }
                        }
                    }
                }
            }
        }
        assertEquals("every combination was covered", 320, configs);
        assertTrue("worst relative residual " + worst + " at " + worstAt, worst < TOL);
    }

    // ----- the two paths against each other ---------------------------------

    @Test
    public void alphaZeroGivesTheSameZerosOnBothPaths() {
        // alpha == 0 means B is set to zero whatever it held, so the two
        // implementations have to agree bit for bit -- including the sign of
        // the zero, and including entries that were not finite
        for (Side side : Side.values()) {
            for (UpLo uplo : UpLo.values()) {
                for (Diag diag : Diag.values()) {
                    for (int n : new int[] { 8, 64 }) {
                        double[] a = triangular(n, 5L);
                        double[] seed = dense(n * n, 6L);
                        // every entry negative, so a multiplication would give -0.0
                        for (int i = 0; i < seed.length; ++i) {
                            seed[i] = -Math.abs(seed[i]) - 1.0;
                        }
                        seed[0] = Double.NaN;
                        seed[1] = Double.POSITIVE_INFINITY;
                        seed[2] = Double.NEGATIVE_INFINITY;
                        double[] blocked = seed.clone();
                        double[] netlib = seed.clone();
                        blocked(side, uplo, Trans.NO_TRANS, diag, n, n, 0.0, a, n, blocked, n);
                        netlib(side, uplo, Trans.NO_TRANS, diag, n, n, 0.0, a, n, netlib, n);
                        String at = " (" + side + "/" + uplo + "/" + diag + ", n = " + n + ")";
                        for (int i = 0; i < blocked.length; ++i) {
                            assertSameBits("entry " + i + at, netlib[i], blocked[i]);
                            assertSameBits("entry " + i + " is +0.0" + at, 0.0, blocked[i]);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void theDispatchThresholdDoesNotChangeTheAnswer() {
        // through the public entry point, at the sizes where it switches
        // implementation: LEFT above 190, RIGHT above 280
        int[][] cases = { { 190, 191 }, { 280, 281 } };
        for (int c = 0; c < cases.length; ++c) {
            Side side = (c == 0) ? Side.LEFT : Side.RIGHT;
            for (int k = 0; k < cases[c].length; ++k) {
                int n = cases[c][k];
                double[] a = triangular(n, 21L);
                double[] b = new double[n * n];
                Arrays.fill(b, -2.5);
                b[0] = Double.NaN;
                b[1] = Double.POSITIVE_INFINITY;
                Dtrmm.dtrmm(side, UpLo.UPPER, Trans.NO_TRANS, Diag.NON_UNIT, n, n, 0.0, a, 0, n, b, 0, n);
                String at = " (" + side + ", " + n + "x" + n + ")";
                for (int i = 0; i < b.length; ++i) {
                    assertSameBits("entry " + i + at, 0.0, b[i]);
                }
            }
        }
    }

    // ----- the pieces underneath --------------------------------------------

    @Test
    public void geaxpyHonoursItsStartOffsetInTheGeneralBranch() {
        // the general branch is the one taken when X and Y are neither both
        // column major nor both row major, which is what Ugemm reaches
        // whenever C is not column major. Dropping X_start cannot throw -- the
        // index only gets smaller -- it silently reads the wrong data
        int m = 3;
        int n = 2;
        int incRowX = 1;
        int incColX = BlockSizes.MR;
        int incRowY = 4;
        int incColY = 1;
        for (int xStart = 0; xStart <= 7; ++xStart) {
            double[] x = new double[xStart + 64];
            for (int i = 0; i < x.length; ++i) {
                x[i] = i + 1;
            }
            double[] got = new double[64];
            double[] want = new double[64];
            Geaxpy.geaxpy(m, n, 2.0, xStart, x, incRowX, incColX, 0, got, incRowY, incColY);
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    want[i * incRowY + j * incColY] += 2.0 * x[xStart + i * incRowX + j * incColX];
                }
            }
            assertTrue("X_start = " + xStart + ", got " + Arrays.toString(Arrays.copyOf(got, 12)) + ", want "
                    + Arrays.toString(Arrays.copyOf(want, 12)), Arrays.equals(got, want));
        }
    }

    @Test
    public void anEmptyDimensionReturnsWithoutTouchingB() {
        int[][] shapes = { { 0, 5 }, { 5, 0 }, { 0, 0 } };
        for (int s = 0; s < shapes.length; ++s) {
            int m = shapes[s][0];
            int n = shapes[s][1];
            double[] a = triangular(8, 3L);
            double[] b = new double[8];
            Arrays.fill(b, 9.0);
            double[] before = b.clone();
            Dtrmm.dtrmm(Side.LEFT, UpLo.UPPER, Trans.NO_TRANS, Diag.NON_UNIT, m, n, 2.0, a, 0, 8, b, 0,
                    Math.max(1, m));
            assertTrue("m = " + m + ", n = " + n + " left B alone", Arrays.equals(before, b));
        }
    }

    // ----- helpers -----------------------------------------------------------

    /** {@code Trmm.trmm} takes the inverted flags, exactly as Dtrmm passes them. */
    private static void blocked(Side side, UpLo uplo, Trans trans, Diag diag, int m, int n, double alpha,
            double[] a, int lda, double[] b, int ldb) {
        Trmm.trmm(side == Side.LEFT, uplo != UpLo.UPPER, trans != Trans.NO_TRANS, diag != Diag.NON_UNIT, m, n,
                alpha, a, 0, lda, b, 0, ldb);
    }

    private static void netlib(Side side, UpLo uplo, Trans trans, Diag diag, int m, int n, double alpha,
            double[] a, int lda, double[] b, int ldb) {
        DtrmmNetlib.dtrmm(side == Side.LEFT, uplo == UpLo.UPPER, trans == Trans.NO_TRANS, diag == Diag.NON_UNIT,
                m, n, alpha, a, 0, lda, b, 0, ldb);
    }

    /** {@code B := alpha * op(A) * B} or {@code B := alpha * B * op(A)}, column major. */
    private static double[] naive(Side side, UpLo uplo, Trans trans, Diag diag, int m, int n, double alpha,
            double[] a, int lda, double[] b) {
        int k = (side == Side.LEFT) ? m : n;
        double[] t = new double[k * k];
        for (int j = 0; j < k; ++j) {
            for (int i = 0; i < k; ++i) {
                boolean inTriangle = (uplo == UpLo.UPPER) ? (i <= j) : (i >= j);
                double v = inTriangle ? a[j * lda + i] : 0.0;
                if (i == j && diag == Diag.UNIT) {
                    v = 1.0;
                }
                if (trans == Trans.NO_TRANS) {
                    t[j * k + i] = v;
                } else {
                    t[i * k + j] = v;
                }
            }
        }
        double[] r = new double[m * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                double s = 0.0;
                if (side == Side.LEFT) {
                    for (int p = 0; p < m; ++p) {
                        s += t[p * m + i] * b[j * m + p];
                    }
                } else {
                    for (int p = 0; p < n; ++p) {
                        s += b[p * m + i] * t[j * n + p];
                    }
                }
                r[j * m + i] = alpha * s;
            }
        }
        return r;
    }

    private static double residual(double[] got, double[] want) {
        double maxWant = 0.0;
        double maxDiff = 0.0;
        for (int i = 0; i < got.length; ++i) {
            maxWant = Math.max(maxWant, Math.abs(want[i]));
            maxDiff = Math.max(maxDiff, Math.abs(got[i] - want[i]));
        }
        return maxDiff / Math.max(maxWant, 1e-300);
    }

    /**
     * A full n x n array, column major, with a diagonal kept away from zero.
     * Both halves carry data on purpose: whichever one {@code uplo} does not
     * name must not be referenced, and filling it is what would catch a
     * blocked path that packs the wrong triangle.
     */
    private static double[] triangular(int n, long seed) {
        double[] a = new double[n * n];
        long lcg = seed;
        for (int i = 0; i < a.length; ++i) {
            lcg = next(lcg);
            a[i] = unit(lcg);
        }
        for (int j = 0; j < n; ++j) {
            a[j * n + j] += 1.5;
        }
        return a;
    }

    private static double[] dense(int len, long seed) {
        double[] b = new double[len];
        long lcg = seed;
        for (int i = 0; i < len; ++i) {
            lcg = next(lcg);
            b[i] = unit(lcg);
        }
        return b;
    }

    private static long next(long lcg) {
        return lcg * 6364136223846793005L + 1442695040888963407L;
    }

    private static double unit(long lcg) {
        return ((lcg >>> 11) * 0x1.0p-53) - 0.5;
    }

    private static void assertSameBits(String what, double expected, double actual) {
        assertEquals(what + " : expected " + expected + " but was " + actual,
                Double.doubleToRawLongBits(expected), Double.doubleToRawLongBits(actual));
    }
}
