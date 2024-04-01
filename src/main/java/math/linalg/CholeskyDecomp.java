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
package math.linalg;

/**
 * Computes the Cholesky decomposition of a symmetric positive semidefinite
 * matrix <code>A = L * L<sup>T</sup></code>.
 */
public final class CholeskyDecomp {

    private static boolean isSquareMatrix(DMatrix A) {
        return A.numRows() == A.numColumns();
    }

    private static boolean isSymmetricMatrix(DMatrix A) {
        for (int col_ = 0; col_ < A.numColumns(); ++col_) {
            for (int row_ = 0; row_ < A.numRows(); ++row_) {
                if (A.getUnsafe(col_, row_) != A.getUnsafe(row_, col_)) {
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
        if (!isSquareMatrix(A)) {
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
            if (L.getUnsafe(i, i) <= 0.0) {
                throw new RuntimeException("(covariance) matrix is not positive semidefinite");
            }
        }
        return L;
    }

    private CholeskyDecomp() {
        throw new AssertionError();
    }
}
