/*
 * Copyright 2024 Stefan Zobel
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
     *             for exactly singular factors in the LU decomposition of a
     *             quadratic matrix or for a non-quadratic matrix that doesn't
     *             have full rank
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

        boolean success = Dgesv.dgesv(A.numRows(), B.numColumns(), A.getArrayUnsafe().clone(), 0,
                Math.max(1, A.numRows()), new int[A.numRows()], 0, X.getArrayUnsafe(), 0, Math.max(1, A.numRows()));

        if (!success) {
            throw new RuntimeException(
                    "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
        }

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
            success = Dgels.dgels(Trans.NO_TRANS, mm, nn, rhsCount, A.getArrayUnsafe().clone(), 0, Math.max(1, mm),
                    tmp.getArrayUnsafe(), 0, Math.max(1, Math.max(mm, nn)), work, 0, work.length);

            if (success) {
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

    private static void checkSolve(DMatrix A, DMatrix B, DMatrix X) {
        checkSameRows(A, B);
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
