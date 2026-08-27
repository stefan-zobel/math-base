/*
 * Copyright 2024, 2026 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.linalg;

/**
 * Computes the Cholesky decomposition of a symmetric positive semidefinite
 * matrix <code>A = L * L<sup>T</sup></code>.
 */
public final class CholeskyDecomp {

    private static final double TOL = 1.0e-10;

    private static boolean isSymmetricMatrix(DMatrix A) {
        for (int col_ = 0; col_ < A.numColumns(); ++col_) {
            for (int row_ = 0; row_ < A.numRows(); ++row_) {
                if (Math.abs(A.getUnsafe(col_, row_) - A.getUnsafe(row_, col_)) > TOL) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Returns the Cholesky factor {@code L} of a positive semidefinite matrix
     * <code>A = L * L<sup>T</sup></code>.
     * 
     * @param A
     *            a positive semidefinite matrix
     * @return the Cholesky factor {@code L}
     */
    public static DMatrix cholesky(DMatrix A) {
        if (!A.isSquareMatrix()) {
            throw new IllegalArgumentException("matrix is not square");
        }
        if (!isSymmetricMatrix(A)) {
            throw new IllegalArgumentException("matrix is not symmetric");
        }
        DMatrix L = new DMatrix(A.numRows(), A.numRows());
        for (int i = 0; i < L.numRows(); ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                    sum += L.getUnsafe(i, k) * L.getUnsafe(j, k);
                }
                if (i == j) {
                    L.setUnsafe(i, i, Math.sqrt(A.getUnsafe(i, i) - sum));
                } else {
                    L.setUnsafe(i, j, (A.getUnsafe(i, j) - sum) * (1.0 / L.getUnsafe(j, j)));
                }
            }
            if (Double.isNaN(L.getUnsafe(i, i))) {
                throw new RuntimeException("(covariance) matrix is not positive semidefinite");
            }
        }
        return L;
    }

    /**
     * Solves <code>L * y = b</code> by forward substitution, where {@code L} is
     * a Cholesky factor as {@link #cholesky(DMatrix)} returns it.
     * <p>
     * Only the lower triangle of {@code L} is read, which is the convention
     * LAPACK uses and which makes the contract exact rather than checked: what
     * stands above the diagonal is ignored, not required to be zero.
     * <p>
     * {@code out} may be the same array as {@code b}.
     *
     * @param L
     *            a Cholesky factor, square, with a non-zero diagonal
     * @param b
     *            the right hand side, as long as {@code L} has rows. Not
     *            modified unless it is also {@code out}
     * @param out
     *            where the solution is written. Its previous contents are
     *            overwritten
     * @throws IllegalArgumentException
     *             if any argument is {@code null}, if {@code L} is not square,
     *             if {@code b} or {@code out} is not as long as {@code L} has
     *             rows, or if the diagonal of {@code L} holds a zero
     * @since 1.5.3
     */
    public static void solveLower(DMatrix L, double[] b, double[] out) {
        int n = checkFactor(L, b, out);
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int k = 0; k < i; ++k) {
                sum -= L.getUnsafe(i, k) * out[k];
            }
            out[i] = sum / L.getUnsafe(i, i);
        }
    }

    /**
     * Solves <code>A * x = b</code> for <code>A = L * L<sup>T</sup></code>, by
     * forward substitution through {@code L} and back substitution through
     * <code>L<sup>T</sup></code>.
     * <p>
     * This is what a Cholesky factor is for, and it is not the same as
     * multiplying by {@link DMatrix#inverse()}: the inverse is a separate LU
     * factorization and squares the condition number of {@code A}, for a
     * quantity {@code L} already holds. Column by column against the unit
     * vectors it also gives the inverse itself, more accurately than forming
     * one.
     * <p>
     * As in {@link #solveLower(DMatrix, double[], double[])}, only the lower
     * triangle of {@code L} is read and {@code out} may be {@code b}.
     *
     * @param L
     *            a Cholesky factor, square, with a non-zero diagonal
     * @param b
     *            the right hand side, as long as {@code L} has rows. Not
     *            modified unless it is also {@code out}
     * @param out
     *            where the solution is written. Its previous contents are
     *            overwritten
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #solveLower(DMatrix, double[], double[])} states
     * @since 1.5.3
     */
    public static void solve(DMatrix L, double[] b, double[] out) {
        int n = checkFactor(L, b, out);
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int k = 0; k < i; ++k) {
                sum -= L.getUnsafe(i, k) * out[k];
            }
            out[i] = sum / L.getUnsafe(i, i);
        }
        for (int i = n - 1; i >= 0; --i) {
            double sum = out[i];
            for (int k = i + 1; k < n; ++k) {
                // the transpose, so the factor is indexed the other way round
                sum -= L.getUnsafe(k, i) * out[k];
            }
            out[i] = sum / L.getUnsafe(i, i);
        }
    }

    /**
     * The natural logarithm of the determinant of
     * <code>A = L * L<sup>T</sup></code> -- <b>of {@code A} and not of
     * {@code L}</b>, which is twice as much, and getting the two the wrong way
     * round is a factor of two in every normalizing constant downstream.
     * <p>
     * Formed as {@code 2 * sum log(L_ii)}, which answers where the determinant
     * itself cannot: a 50 by 50 covariance of unit scale has a determinant well
     * below the smallest {@code double} while its logarithm is an ordinary
     * number.
     *
     * @param L
     *            a Cholesky factor, square, with a strictly positive diagonal
     * @return the natural logarithm of {@code det(A)}
     * @throws IllegalArgumentException
     *             if {@code L} is {@code null}, is not square, or holds a
     *             diagonal entry that is not strictly positive
     * @since 1.5.3
     */
    public static double logDeterminant(DMatrix L) {
        if (L == null) {
            throw new IllegalArgumentException("L must not be null");
        }
        if (!L.isSquareMatrix()) {
            throw new IllegalArgumentException("L is not square");
        }
        double sum = 0.0;
        for (int i = 0; i < L.numRows(); ++i) {
            double d = L.getUnsafe(i, i);
            if (!(d > 0.0)) {
                throw new IllegalArgumentException(
                        "L[" + i + "][" + i + "] is not strictly positive : " + d);
            }
            sum += Math.log(d);
        }
        return 2.0 * sum;
    }

    /** Shared validation of a factor and the two vectors it is applied to. */
    private static int checkFactor(DMatrix L, double[] b, double[] out) {
        if (L == null) {
            throw new IllegalArgumentException("L must not be null");
        }
        if (b == null || out == null) {
            throw new IllegalArgumentException("b and out must not be null");
        }
        if (!L.isSquareMatrix()) {
            throw new IllegalArgumentException("L is not square");
        }
        int n = L.numRows();
        if (b.length != n) {
            throw new IllegalArgumentException("b must hold " + n + " entries, not " + b.length);
        }
        if (out.length != n) {
            throw new IllegalArgumentException("out must hold " + n + " entries, not " + out.length);
        }
        for (int i = 0; i < n; ++i) {
            if (L.getUnsafe(i, i) == 0.0) {
                throw new IllegalArgumentException("L[" + i + "][" + i + "] is zero, so L is singular");
            }
        }
        return n;
    }

    private CholeskyDecomp() {
        throw new AssertionError();
    }
}
