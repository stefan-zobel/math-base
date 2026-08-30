package math.ts;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;

/**
 * The validation both models share.
 * <p>
 * It lives here rather than in either of them because it is <em>one policy</em>
 * -- what counts as symmetric, and what a covariance is allowed to be -- and
 * two copies of a tolerance are two tolerances waiting to drift apart.
 */
final class Models {

    /**
     * How far the two halves of a covariance may differ before it is called
     * asymmetric, relative to the largest entry on its diagonal. The same
     * policy {@link MultivariateNormal} applies, and for the same reason: a
     * covariance accumulated from data is symmetric only up to its own
     * round-off, which at a scale of {@code 1e6} is already larger than the
     * absolute tolerance {@link CholeskyDecomp} uses.
     */
    private static final double SYMMETRY_TOL = 1.0e-10;

    /**
     * Reads the matrix as symmetric: every entry finite, the two halves equal
     * within a tolerance relative to the diagonal, and the copy averaged rather
     * than taken from one side, so that what is stored is exactly symmetric
     * whichever half was the rounder one.
     */
    static DMatrix symmetrized(DMatrix A, int order, String name) {
        if (!A.isSquareMatrix()) {
            throw new IllegalArgumentException(name + " is not square");
        }
        if (A.numRows() != order) {
            throw new IllegalArgumentException(
                    name + " is of order " + A.numRows() + ", not " + order);
        }
        double scale = 0.0;
        for (int i = 0; i < order; ++i) {
            double v = A.get(i, i);
            if (!isFinite(v)) {
                throw new IllegalArgumentException(name + "[" + i + "][" + i + "] is not finite : " + v);
            }
            if (v < 0.0) {
                throw new IllegalArgumentException(
                        name + "[" + i + "][" + i + "] is a negative variance : " + v);
            }
            scale = Math.max(scale, v);
        }
        // a zero matrix is a legitimate Q -- a state that moves without noise
        double tol = SYMMETRY_TOL * Math.max(scale, 1.0);

        DMatrix out = new DMatrix(order, order);
        for (int i = 0; i < order; ++i) {
            out.set(i, i, A.get(i, i));
            for (int j = i + 1; j < order; ++j) {
                double a = A.get(i, j);
                double b = A.get(j, i);
                if (!isFinite(a) || !isFinite(b)) {
                    throw new IllegalArgumentException(
                            name + "[" + i + "][" + j + "] is not finite : " + a + ", " + b);
                }
                if (Math.abs(a - b) > tol) {
                    throw new IllegalArgumentException(name + " is not symmetric at (" + i + ", " + j
                            + ") : " + a + " against " + b);
                }
                double v = 0.5 * (a + b);
                out.set(i, j, v);
                out.set(j, i, v);
            }
        }
        return out;
    }

    /** The initial covariance, symmetrized and proved positive definite. */
    static DMatrix initialCovariance(DMatrix P0, int order) {
        DMatrix out = symmetrized(P0, order, "P0");
        DMatrix factor;
        try {
            factor = CholeskyDecomp.cholesky(out);
        } catch (RuntimeException e) {
            throw new IllegalArgumentException("P0 is not positive definite : " + e.getMessage());
        }
        for (int i = 0; i < order; ++i) {
            if (!(factor.get(i, i) > 0.0)) {
                throw new IllegalArgumentException(
                        "P0 is singular: its factor has a zero at (" + i + ", " + i + ")");
            }
        }
        return out;
    }

    /** The initial mean, checked and copied. */
    static double[] initialMean(double[] m0, int order) {
        if (m0 == null) {
            throw new IllegalArgumentException("m0 must not be null");
        }
        if (m0.length != order) {
            throw new IllegalArgumentException("m0 is of length " + m0.length + ", not " + order);
        }
        for (int i = 0; i < order; ++i) {
            if (!isFinite(m0[i])) {
                throw new IllegalArgumentException("m0[" + i + "] is not finite : " + m0[i]);
            }
        }
        return m0.clone();
    }

    /** The transition, checked for shape and finiteness, and copied. */
    static DMatrix transition(DMatrix F, int order, String name) {
        if (!F.isSquareMatrix()) {
            throw new IllegalArgumentException(name + " is not square");
        }
        if (F.numRows() != order) {
            throw new IllegalArgumentException(
                    name + " is of order " + F.numRows() + ", not " + order);
        }
        checkFinite(F, name);
        return F.copy();
    }

    /** The observation matrix, checked for shape and finiteness, and copied. */
    static DMatrix observation(DMatrix H, int rows, int cols, String name) {
        if (H.numRows() != rows) {
            throw new IllegalArgumentException(name + " has " + H.numRows() + " rows, not " + rows);
        }
        if (H.numColumns() != cols) {
            throw new IllegalArgumentException(
                    name + " has " + H.numColumns() + " columns, not " + cols);
        }
        checkFinite(H, name);
        return H.copy();
    }

    static void checkFinite(DMatrix A, String name) {
        for (int i = 0; i < A.numRows(); ++i) {
            for (int j = 0; j < A.numColumns(); ++j) {
                double v = A.get(i, j);
                if (!isFinite(v)) {
                    throw new IllegalArgumentException(name + "[" + i + "][" + j + "] is not finite : " + v);
                }
            }
        }
    }

    /** What every accessor writing into a caller's array has to check first. */
    static void checkOut(double[] out, int order) {
        if (out == null) {
            throw new IllegalArgumentException("out must not be null");
        }
        if (out.length != order) {
            throw new IllegalArgumentException("out is of length " + out.length + ", not " + order);
        }
    }

    static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private Models() {
        throw new AssertionError();
    }
}
