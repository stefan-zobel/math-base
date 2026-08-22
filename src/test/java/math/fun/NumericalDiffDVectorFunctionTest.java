package math.fun;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

/**
 * Covers {@link NumericalDiffDVectorFunction} and, since it had never had one,
 * {@link NumericalDiffDMultiFunction} beside it.
 */
public class NumericalDiffDVectorFunctionTest {

    /** What a forward difference can be held to, relative to the scale. */
    private static final double FORWARD_DIFFERENCE = 1.0e-6;

    /**
     * {@code F(x) = (x0^2 + x1^2 + x2^2, sin(x0) * x1, exp(x2) - x0)}, whose
     * Jacobian is written down below it.
     */
    private abstract static class Sphere {

        static void values(double[] x, double[] f) {
            f[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            f[1] = Math.sin(x[0]) * x[1];
            f[2] = Math.exp(x[2]) - x[0];
        }

        /** Column-major, three by three. */
        static void jacobian(double[] x, double[] j) {
            j[0] = 2.0 * x[0];
            j[1] = Math.cos(x[0]) * x[1];
            j[2] = -1.0;
            j[3] = 2.0 * x[1];
            j[4] = Math.sin(x[0]);
            j[5] = 0.0;
            j[6] = 2.0 * x[2];
            j[7] = 0.0;
            j[8] = Math.exp(x[2]);
        }
    }

    private static NumericalDiffDVectorFunction sphere() {
        return new NumericalDiffDVectorFunction(3) {
            @Override
            public void valueAt(double[] x, double[] values) {
                Sphere.values(x, values);
            }
        };
    }

    private static double scaleOf(double[] a) {
        double worst = 0.0;
        for (int i = 0; i < a.length; ++i) {
            worst = Math.max(worst, Math.abs(a[i]));
        }
        return Math.max(worst, 1.0);
    }

    private long lcg = 24680135791113151L;

    private double next() {
        lcg = lcg * 6364136223846793005L + 1442695040888963407L;
        return ((lcg >>> 11) / (double) (1L << 53)) * 4.0 - 2.0;
    }

    @Test
    public void testTheApproximatedJacobianMatchesTheAnalyticOne() {
        NumericalDiffDVectorFunction f = sphere();
        double[] approximated = new double[9];
        double[] exact = new double[9];
        for (int t = 0; t < 200; ++t) {
            double[] x = { next(), next(), next() };
            f.jacobianAt(x, approximated);
            Sphere.jacobian(x, exact);
            double scale = scaleOf(exact);
            for (int i = 0; i < 9; ++i) {
                assertEquals("entry " + i + " at trial " + t, exact[i], approximated[i],
                        FORWARD_DIFFERENCE * scale);
            }
        }
    }

    /**
     * Entry {@code (i, j)} of the result is the derivative of result {@code i}
     * with respect to argument {@code j} and sits at {@code j * m + i}. The
     * function below has one non-zero derivative per cell, so a transposed
     * layout could not pass.
     */
    @Test
    public void testTheLayoutIsTheColumnMajorOneTheInterfaceDocuments() {
        final int m = 4;
        final int n = 3;
        NumericalDiffDVectorFunction f = new NumericalDiffDVectorFunction(m) {
            @Override
            public void valueAt(double[] x, double[] values) {
                for (int i = 0; i < m; ++i) {
                    values[i] = 0.0;
                    for (int j = 0; j < n; ++j) {
                        values[i] += (10.0 * (i + 1) + (j + 1)) * x[j];
                    }
                }
            }
        };
        double[] jacobian = new double[m * n];
        f.jacobianAt(new double[] { 0.5, -1.5, 2.5 }, jacobian);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                assertEquals("(" + i + ", " + j + ")", 10.0 * (i + 1) + (j + 1), jacobian[j * m + i], 1.0e-6);
            }
        }
    }

    /**
     * The one place this class deviates from
     * {@link NumericalDiffDMultiFunction}: the step is floored at one, so a
     * component that happens to sit near zero still produces a column. A step
     * proportional to that component alone would be far below the last bit of
     * the function values and the whole column would come back as exactly
     * zero, which turns an ordinary point into a singular Jacobian.
     */
    @Test
    public void testAComponentNearZeroStillProducesAColumn() {
        NumericalDiffDVectorFunction circles = new NumericalDiffDVectorFunction(2) {
            @Override
            public void valueAt(double[] x, double[] values) {
                values[0] = x[0] * x[0] + x[1] * x[1] - 1.0;
                values[1] = x[0] + x[1] - 1.0;
            }
        };
        double[] jacobian = new double[4];
        double[] smaller = { 1.0e-4, 1.0e-9, 1.0e-14, 0.0 };
        for (int k = 0; k < smaller.length; ++k) {
            circles.jacobianAt(new double[] { smaller[k], 1.0 }, jacobian);
            String where = "x0 = " + smaller[k];
            // the exact first column is (2 * x0, 1), so the second entry is the
            // one that must survive
            assertEquals(where, 1.0, jacobian[1], 1.0e-6);
            assertEquals(where, 2.0, jacobian[2], 1.0e-6);
            assertEquals(where, 1.0, jacobian[3], 1.0e-6);
            double determinant = jacobian[0] * jacobian[3] - jacobian[2] * jacobian[1];
            assertTrue(where + " : the matrix is not singular, determinant " + determinant,
                    Math.abs(determinant) > 0.5);
        }
    }

    /** The argument comes back exactly as it went in. */
    @Test
    public void testTheArgumentIsRestored() {
        NumericalDiffDVectorFunction f = sphere();
        double[] x = { 0.25, -3.5, 1.75 };
        double[] copy = x.clone();
        f.jacobianAt(x, new double[9]);
        for (int i = 0; i < x.length; ++i) {
            assertEquals("component " + i, copy[i], x[i], 0.0);
        }
    }

    /** One Jacobian costs one evaluation per argument, plus the base point. */
    @Test
    public void testAJacobianCostsOneEvaluationPerArgumentPlusOne() {
        final int[] calls = new int[1];
        NumericalDiffDVectorFunction f = new NumericalDiffDVectorFunction(5) {
            @Override
            public void valueAt(double[] x, double[] values) {
                ++calls[0];
                for (int i = 0; i < 5; ++i) {
                    values[i] = x[i % x.length] * (i + 1);
                }
            }
        };
        f.jacobianAt(new double[] { 1.0, 2.0, 3.0 }, new double[15]);
        assertEquals(4, calls[0]);
        calls[0] = 0;
        f.jacobianAt(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 }, new double[35]);
        assertEquals(8, calls[0]);
    }

    @Test
    public void testResultCountIsReported() {
        assertEquals(3, sphere().resultCount());
        assertEquals(11, new NumericalDiffDVectorFunction(11, 1.0e-6) {
            @Override
            public void valueAt(double[] x, double[] values) {
                throw new AssertionError("not evaluated");
            }
        }.resultCount());
    }

    /**
     * A coarser step is coarser, which is the only thing the caller can
     * actually control here and the reason it is a constructor argument at all.
     */
    @Test
    public void testACoarserStepIsCoarser() {
        double[] exact = new double[9];
        double[] x = { 0.7, -1.3, 0.4 };
        Sphere.jacobian(x, exact);
        double previous = -1.0;
        double[] scales = { 1.5 * Math.sqrt(Math.ulp(1.0)), 1.0e-5, 1.0e-3, 1.0e-1 };
        for (int k = 0; k < scales.length; ++k) {
            final double scale = scales[k];
            NumericalDiffDVectorFunction f = new NumericalDiffDVectorFunction(3, scale) {
                @Override
                public void valueAt(double[] xx, double[] values) {
                    Sphere.values(xx, values);
                }
            };
            double[] approximated = new double[9];
            f.jacobianAt(x.clone(), approximated);
            double worst = 0.0;
            for (int i = 0; i < 9; ++i) {
                worst = Math.max(worst, Math.abs(approximated[i] - exact[i]));
            }
            assertTrue("step " + scale + " is not coarser than the one before it : " + worst + " against "
                    + previous, worst > previous);
            previous = worst;
        }
    }

    @Test
    public void testArgumentValidation() {
        int[] badCounts = { 0, -1, Integer.MIN_VALUE };
        for (int k = 0; k < badCounts.length; ++k) {
            try {
                new NumericalDiffDVectorFunction(badCounts[k]) {
                    @Override
                    public void valueAt(double[] x, double[] values) {
                        throw new AssertionError("not evaluated");
                    }
                };
                fail("result count " + badCounts[k] + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertNotNull(expected.getMessage());
            }
        }
        double[] badScales = { 0.0, -1.0e-8, Double.NaN, Double.POSITIVE_INFINITY };
        for (int k = 0; k < badScales.length; ++k) {
            try {
                new NumericalDiffDVectorFunction(3, badScales[k]) {
                    @Override
                    public void valueAt(double[] x, double[] values) {
                        throw new AssertionError("not evaluated");
                    }
                };
                fail("scale " + badScales[k] + " was accepted");
            } catch (IllegalArgumentException expected) {
                assertNotNull(expected.getMessage());
            }
        }
    }

    // --- the scalar sibling, which had no test of its own ---

    private static NumericalDiffDMultiFunction quadratic() {
        return new NumericalDiffDMultiFunction() {
            @Override
            public double apply(double[] x) {
                return 3.0 * x[0] * x[0] + 2.0 * x[0] * x[1] - x[1] * x[1] * x[1] + Math.cos(x[2]);
            }
        };
    }

    private static void quadraticGradient(double[] x, double[] g) {
        g[0] = 6.0 * x[0] + 2.0 * x[1];
        g[1] = 2.0 * x[0] - 3.0 * x[1] * x[1];
        g[2] = -Math.sin(x[2]);
    }

    @Test
    public void testTheApproximatedGradientMatchesTheAnalyticOne() {
        NumericalDiffDMultiFunction f = quadratic();
        double[] approximated = new double[3];
        double[] exact = new double[3];
        for (int t = 0; t < 200; ++t) {
            double[] x = { next(), next(), next() };
            f.derivativeAt(x, approximated);
            quadraticGradient(x, exact);
            double scale = scaleOf(exact);
            for (int i = 0; i < 3; ++i) {
                assertEquals("component " + i + " at trial " + t, exact[i], approximated[i],
                        FORWARD_DIFFERENCE * scale);
            }
        }
    }

    @Test
    public void testTheScalarSiblingRestoresItsArgumentToo() {
        NumericalDiffDMultiFunction f = quadratic();
        double[] x = { 0.25, -3.5, 1.75 };
        double[] copy = x.clone();
        f.derivativeAt(x, new double[3]);
        for (int i = 0; i < x.length; ++i) {
            assertEquals("component " + i, copy[i], x[i], 0.0);
        }
    }

    /**
     * A gradient costs one evaluation per argument plus the base point, the
     * same accounting as the Jacobian above.
     */
    @Test
    public void testAGradientCostsOneEvaluationPerArgumentPlusOne() {
        final int[] calls = new int[1];
        NumericalDiffDMultiFunction f = new NumericalDiffDMultiFunction() {
            @Override
            public double apply(double[] x) {
                ++calls[0];
                double sum = 0.0;
                for (int i = 0; i < x.length; ++i) {
                    sum += x[i] * x[i];
                }
                return sum;
            }
        };
        f.derivativeAt(new double[] { 1.0, 2.0, 3.0, 4.0 }, new double[4]);
        assertEquals(5, calls[0]);
    }
}
