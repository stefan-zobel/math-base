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
package math.solve;

import math.MathConsts;
import math.gemm.Trans;
import math.lapack.Dgels;
import math.lapack.Dgesv;
import math.linalg.DMatrix;

/**
 * Computes the solution {@code X} to a real system of linear equations
 * {@code A * X = B}.
 */
public final class LinearEquationsSolver {

    /**
     * Computes the solution ({@code X}) to a real system of linear equations
     * {@code A * X = B}, where {@code A} is either a {@code n x n} matrix and
     * {@code X} and {@code B} are {@code n x r} matrices, or where {@code A} is
     * a {@code n x m} and matrix {@code X} is a {@code m x r} matrix and
     * {@code B} is a {@code n x r} matrix.
     * 
     * @param A
     *            either a {@code n x n} matrix or a {@code n x m} matrix
     * @param B
     *            matrix with the same number of rows as matrix {@code A}, and
     *            the same number of columns as {@code X}
     * @param X
     *            matrix with number of rows equal to the number of columns of
     *            matrix {@code A}, and the same number of columns as {@code B}
     * @return {@code X}, the solution of dimension either {@code n x r} (in the
     *         {@code n x n} case) or {@code m x r} (in the {@code m x n} case)
     * @throws RuntimeException
     *             if {@code A} is numerically rank deficient, on either path.
     *             The decomposition is rejected when the smallest entry on the
     *             diagonal of its triangular factor drops to
     *             {@code max(m, n) * u} times the largest or below, where
     *             {@code u} is {@link math.MathConsts#MACH_EPS_DBL}. That ratio
     *             bounds the reciprocal of the condition number of {@code A}
     *             from above, so it detects the loss of rank rather than
     *             computing the rank: a rejection means {@code A} is
     *             ill conditioned beyond what double precision can carry, but a
     *             matrix whose triangular factor has a well spread diagonal
     *             while the matrix itself is nearly singular -- Kahan's matrix
     *             is the classical example -- passes undetected, as it must
     *             without column pivoting. Use {@code math.linalg.OLS}, which
     *             is based on a singular value decomposition, when a rank
     *             deficient problem is to be solved rather than rejected
     */
    public static DMatrix solve(DMatrix A, DMatrix B, DMatrix X) {
        checkSolve(A, B, X);
        if (A.isSquareMatrix()) {
            return lusolve(A, B, X);
        }
        return qrsolve(A, B, X);
    }

    private static DMatrix lusolve(DMatrix A, DMatrix B, DMatrix X) {
        // copy B into X
        System.arraycopy(B.getArrayUnsafe(), 0, X.getArrayUnsafe(), 0, B.getArrayUnsafe().length);

        int n = A.numRows();
        int lda = Math.max(1, n);
        // kept rather than passed inline: Dgesv overwrites it with L and U
        double[] lu = A.getArrayUnsafe().clone();

        boolean success = Dgesv.dgesv(n, B.numColumns(), lu, 0, lda, new int[n], 0, X.getArrayUnsafe(), 0, lda);

        if (!success) {
            throw new RuntimeException(
                    "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
        }

        // Dgetf2 only ever compares the pivot against exactly 0.0, which
        // floating point arithmetic hardly ever produces, so the check above
        // misses almost every singular matrix
        checkFullRank(lu, n, n, lda);

        return X;
    }

    private static DMatrix qrsolve(DMatrix A, DMatrix B, DMatrix X) {
        int rhsCount = B.numColumns();
        int mm = A.numRows();
        int nn = A.numColumns();

        DMatrix tmp = new DMatrix(Math.max(mm, nn), rhsCount);
        for (int j = 0; j < rhsCount; ++j) {
            for (int i = 0; i < mm; ++i) {
                tmp.setUnsafe(i, j, B.getUnsafe(i, j));
            }
        }

        double[] work = new double[1];
        boolean success = Dgels.dgels(Trans.NO_TRANS, mm, nn, rhsCount, new double[0], 0, Math.max(1, mm),
                new double[0], 0, Math.max(1, Math.max(mm, nn)), work, 0, -1);

        if (success) {
            work = new double[(int) work[0]];
            int lda = Math.max(1, mm);
            // kept rather than passed inline: Dgels overwrites it with the
            // triangular factor, R of the QR when mm >= nn, L of the LQ
            // otherwise
            double[] factor = A.getArrayUnsafe().clone();
            success = Dgels.dgels(Trans.NO_TRANS, mm, nn, rhsCount, factor, 0, lda, tmp.getArrayUnsafe(), 0,
                    Math.max(1, Math.max(mm, nn)), work, 0, work.length);

            if (success) {
                // Dtrtrs only ever compares the diagonal against exactly 0.0,
                // which a Householder factorization hardly ever produces, so
                // Dgels reports full rank for almost every deficient matrix
                checkFullRank(factor, mm, nn, lda);

                for (int j = 0; j < rhsCount; ++j) {
                    for (int i = 0; i < nn; ++i) {
                        X.setUnsafe(i, j, tmp.getUnsafe(i, j));
                    }
                }
            }
        }

        if (!success) {
            throw new RuntimeException("A does not have full rank. Least squares solution could not be computed.");
        }

        return X;
    }

    /**
     * Rejects a factorization whose triangular factor is numerically rank
     * deficient.
     * <p>
     * Both {@code Dgesv} and {@code Dgels} overwrite their first argument with
     * the triangular factor of the decomposition -- {@code U} of the LU,
     * {@code R} of the QR, {@code L} of the LQ -- so its diagonal is available
     * here at no arithmetic cost. A smallest entry that is tiny relative to the
     * largest is the standard indicator that the factor, and with it {@code A},
     * has lost rank.
     * <p>
     * On the QR and LQ paths this is more than an indicator. Since
     * {@code sigma_min <= min|d_i|} and {@code max|d_i| <= sigma_max}, the ratio
     * bounds {@code 1/cond_2(A)} from above, so a value below the threshold
     * proves the condition number lies above its reciprocal and the rejection
     * is never a false alarm. On the LU path element growth breaks that bound
     * and the test is a heuristic -- still far better than the comparison
     * against exactly zero that it backs up.
     *
     * @param factor
     *            the array the LAPACK routine overwrote with the triangular
     *            factor, in column major order
     * @param m
     *            the number of rows of {@code A}
     * @param n
     *            the number of columns of {@code A}
     * @param lda
     *            the leading dimension of {@code factor}
     */
    private static void checkFullRank(double[] factor, int m, int n, int lda) {
        int k = Math.min(m, n);
        if (k == 0) {
            // Dgels quick returns without factorizing, there is no diagonal
            return;
        }
        double min = Double.MAX_VALUE;
        double max = 0.0;
        for (int i = 0; i < k; ++i) {
            double d = Math.abs(factor[i + i * lda]);
            // both comparisons below are false for a NaN, which would drop it
            // out of the scan unnoticed and leave a finite looking ratio, so it
            // has to be caught here
            if (!Double.isFinite(d)) {
                throw new RuntimeException("The decomposition of A is not finite: entry " + i + " on the diagonal"
                        + " of the triangular factor is " + factor[i + i * lda]
                        + ". Solution could not be computed.");
            }
            if (d < min) {
                min = d;
            }
            if (d > max) {
                max = d;
            }
        }
        // the threshold scales with max(m, n) because the backward error of the
        // factorization does, that being the noise floor the diagonal has to
        // stand above. MACH_EPS_DBL is the unit roundoff u = eps/2, so this is
        // half of the max(m, n) * eps that numpy and MATLAB use for the same
        // purpose
        double tol = Math.max(m, n) * MathConsts.MACH_EPS_DBL;
        double ratio = min / max;
        // negated so that a NaN ratio is rejected as well. It arises from a NaN
        // in A, and from the all zero matrix, which Dgels short circuits before
        // factorizing and which therefore leaves min == max == 0
        if (!(ratio > tol)) {
            throw new RuntimeException("A does not have full rank: the smallest diagonal entry of the triangular"
                    + " factor is " + ratio + " times the largest, at or below the threshold " + tol
                    + ". Solution could not be computed.");
        }
    }

    /**
     * Rejects a non-finite entry before it reaches LAPACK.
     * <p>
     * This is not a matter of taste. {@code Dlascl}, which both {@code Dgels}
     * and {@code Dgesv} reach when the norm of their argument leaves the safe
     * range, scales by repeatedly multiplying {@code cfrom} by the smallest
     * normal number until the product is small enough. For an infinite
     * {@code cfrom} that product stays infinite, so the loop never terminates
     * and the call hangs forever rather than returning a wrong answer. A
     * {@code NaN} does not hang but silently produces a {@code NaN} solution.
     * Both are better reported here.
     *
     * @param M
     *            the matrix to scan
     * @param name
     *            the name to use in the exception message
     */
    private static void checkFinite(DMatrix M, String name) {
        double[] a = M.getArrayUnsafe();
        for (int i = 0; i < a.length; ++i) {
            if (!Double.isFinite(a[i])) {
                throw new IllegalArgumentException("Matrix " + name + " has a non-finite entry (" + a[i]
                        + ") at row " + (i % M.numRows()) + ", column " + (i / M.numRows()));
            }
        }
    }

    private static void checkSolve(DMatrix A, DMatrix B, DMatrix X) {
        checkSameRows(A, B);
        checkFinite(A, "A");
        checkFinite(B, "B");
        if (A.numColumns() != X.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != X.numRows() (" + A.numColumns() + " != " + X.numRows() + ")");
        }
        if (X.numColumns() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "X.numColumns() != B.numColumns() (" + X.numColumns() + " != " + B.numColumns() + ")");
        }
    }

    private static void checkSameRows(DMatrix A, DMatrix B) {
        if (A.numRows() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numRows() (" + A.numRows() + " != " + B.numRows() + ")");
        }
    }

    private LinearEquationsSolver() {
        throw new AssertionError();
    }
}
