package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * {@link ScaledDesign}, which was package private until it acquired a consumer
 * in another package and had no test at all.
 * <p>
 * What is pinned here is the <em>parameterization</em> rather than the
 * arithmetic: a penalty or a prior means what it means only relative to this
 * divisor, so an estimator in another package that fits on this transform is
 * relying on the root mean square and not on the sample standard deviation.
 * Changing that divisor would silently reparameterize every {@code lambda} in
 * the library.
 */
public final class ScaledDesignTest {

    /** A design with columns of deliberately different scale and offset. */
    private static double[] design(int n, int p, long seed) {
        double[] x = new double[n * p];
        long lcg = seed;
        for (int j = 0; j < p; j++) {
            double offset = (j + 1) * 10.0;
            double spread = Math.pow(10.0, j - 1);
            for (int i = 0; i < n; i++) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                x[j * n + i] = offset + spread * (((lcg >>> 11) * 0x1.0p-53) * 2.0 - 1.0);
            }
        }
        return x;
    }

    private static double[] response(int n, long seed) {
        double[] y = new double[n];
        long lcg = seed;
        for (int i = 0; i < n; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            y[i] = 5.0 + 3.0 * (((lcg >>> 11) * 0x1.0p-53) * 2.0 - 1.0);
        }
        return y;
    }

    @Test
    public void everyColumnComesOutCenteredAndOfUnitRootMeanSquare() {
        // sum_i x_ij^2 == n exactly is what CoordinateDescent's inner loop
        // divides by, so it is the contract and not a coincidence
        for (int n : new int[] { 5, 20, 137 }) {
            for (int p : new int[] { 1, 3, 7 }) {
                ScaledDesign s = ScaledDesign.of(design(n, p, 42L), response(n, 7L), n, p, null);
                assertEquals("rows", n, s.n);
                assertEquals("columns", p, s.p);
                for (int j = 0; j < p; j++) {
                    double mean = 0.0;
                    double ss = 0.0;
                    for (int i = 0; i < n; i++) {
                        double v = s.x[j * n + i];
                        mean += v;
                        ss += v * v;
                    }
                    assertEquals("n=" + n + " p=" + p + ": column " + j + " is not centered", 0.0, mean / n,
                            1.0e-12);
                    assertEquals("n=" + n + " p=" + p + ": column " + j + " is not of unit mean square", 1.0,
                            ss / n, 1.0e-12);
                }
                double yMean = 0.0;
                for (int i = 0; i < n; i++) {
                    yMean += s.y[i];
                }
                assertEquals("the response is not centered", 0.0, yMean / n, 1.0e-12);
            }
        }
    }

    @Test
    public void theDivisorIsTheRootMeanSquareAndNotTheStandardDeviation() {
        // the whole reason this class is not Standardization, and the thing
        // that would reparameterize every lambda in the library if it changed
        int n = 40;
        int p = 3;
        double[] xs = design(n, p, 11L);
        ScaledDesign s = ScaledDesign.of(xs, response(n, 3L), n, p, null);
        for (int j = 0; j < p; j++) {
            double mean = 0.0;
            for (int i = 0; i < n; i++) {
                mean += xs[j * n + i];
            }
            mean /= n;
            double ss = 0.0;
            for (int i = 0; i < n; i++) {
                double c = xs[j * n + i] - mean;
                ss += c * c;
            }
            assertEquals("the column mean", mean, s.xBar[j], 1.0e-12 * Math.max(1.0, Math.abs(mean)));
            assertEquals("the divisor is not sqrt(ss/n)", Math.sqrt(ss / n), s.scale[j],
                    1.0e-12 * Math.sqrt(ss / n));
            // and it is measurably not the sample standard deviation
            assertTrue("the two divisors are indistinguishable at n=" + n,
                    Math.abs(s.scale[j] - Math.sqrt(ss / (n - 1))) > 1.0e-6 * s.scale[j]);
        }
    }

    @Test
    public void unscaleInvertsTheScaling() {
        int n = 25;
        int p = 4;
        ScaledDesign s = ScaledDesign.of(design(n, p, 5L), response(n, 9L), n, p, null);
        double[] scaled = { 1.5, -2.0, 0.25, 7.0 };
        double[] original = s.unscale(scaled);
        for (int j = 0; j < p; j++) {
            assertEquals("unscale at " + j, scaled[j] / s.scale[j], original[j], 0.0);
            // and the fitted contribution is the same either way, which is what
            // unscaling is for
            assertEquals("the contribution changed at " + j, scaled[j] * s.x[j * n],
                    original[j] * (s.x[j * n] * s.scale[j]), 1.0e-12);
        }
    }

    @Test
    public void theInterceptIsTheOneItPromises() {
        int n = 30;
        int p = 3;
        double[] xs = design(n, p, 17L);
        double[] ys = response(n, 21L);
        ScaledDesign s = ScaledDesign.of(xs, ys, n, p, null);
        double[] beta = { 0.5, -1.25, 2.0 };
        double want = s.yBar;
        for (int j = 0; j < p; j++) {
            want -= beta[j] * s.xBar[j];
        }
        assertEquals("the intercept", want, s.intercept(beta), 1.0e-12 * Math.max(1.0, Math.abs(want)));

        // and it is what makes the fitted values average to ybar
        double fittedMean = 0.0;
        for (int i = 0; i < n; i++) {
            double v = s.intercept(beta);
            for (int j = 0; j < p; j++) {
                v += beta[j] * xs[j * n + i];
            }
            fittedMean += v;
        }
        assertEquals("the fitted values do not average to yBar", s.yBar, fittedMean / n, 1.0e-10);
    }

    @Test
    public void theRowSubsetSelectsWhatItSays() {
        int rows = 40;
        int p = 2;
        double[] xs = design(rows, p, 31L);
        double[] ys = response(rows, 33L);
        int[] use = { 3, 7, 11, 12, 20, 21, 22, 30 };

        ScaledDesign sub = ScaledDesign.of(xs, ys, rows, p, use);
        assertEquals("the subset has the wrong number of rows", use.length, sub.n);

        // the same rows copied out first must give the same transform
        double[] xsCopy = new double[use.length * p];
        double[] ysCopy = new double[use.length];
        for (int i = 0; i < use.length; i++) {
            ysCopy[i] = ys[use[i]];
            for (int j = 0; j < p; j++) {
                xsCopy[j * use.length + i] = xs[j * rows + use[i]];
            }
        }
        ScaledDesign direct = ScaledDesign.of(xsCopy, ysCopy, use.length, p, null);
        for (int j = 0; j < p; j++) {
            assertEquals("the subset mean differs at " + j, direct.xBar[j], sub.xBar[j], 0.0);
            assertEquals("the subset scale differs at " + j, direct.scale[j], sub.scale[j], 0.0);
            for (int i = 0; i < use.length; i++) {
                assertEquals("the subset design differs at (" + i + ", " + j + ")",
                        direct.x[j * use.length + i], sub.x[j * use.length + i], 0.0);
            }
        }
        assertEquals("the subset response mean differs", direct.yBar, sub.yBar, 0.0);
    }

    @Test
    public void itHoldsItsOwnArrays() {
        int n = 12;
        int p = 2;
        double[] xs = design(n, p, 61L);
        double[] ys = response(n, 63L);
        ScaledDesign s = ScaledDesign.of(xs, ys, n, p, null);
        double before = s.x[0];
        double yBefore = s.y[0];
        xs[0] = 999.0;
        ys[0] = 999.0;
        assertEquals("the transform followed the caller's design", before, s.x[0], 0.0);
        assertEquals("the transform followed the caller's response", yBefore, s.y[0], 0.0);
    }

    @Test
    public void aConstantColumnIsRefused() {
        // pinned as it stands: a constant column carries no information once
        // the intercept is fitted separately, and it would divide by zero
        int n = 10;
        int p = 2;
        double[] xs = new double[n * p];
        for (int i = 0; i < n; i++) {
            xs[i] = 3.0;
            xs[n + i] = i;
        }
        try {
            ScaledDesign.of(xs, response(n, 1L), n, p, null);
            fail("a constant column was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue("threw without a message", expected.getMessage() != null);
            assertTrue("the message does not name the column : " + expected.getMessage(),
                    expected.getMessage().contains("column 0"));
        }
    }
}
