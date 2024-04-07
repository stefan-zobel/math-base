/*
 * Copyright 2018 Stefan Zobel
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
package math.trmm;

import java.util.Objects;

import math.gemm.Trans;

/**
 * {@code DTRMM} performs one of the matrix-matrix operations
 * {@code B := alpha*op( A )*B}, or {@code B := alpha*B*op( A )}, where alpha is
 * a scalar, B is an m by n matrix, A is a unit, or non-unit, upper or lower
 * triangular matrix and {@code op( A )} is one of {@code op( A ) = A} or
 * {@code op( A ) = A**T}.
 */
public final class Dtrmm {
    /**
     * <pre>
     * <code>
     *  Purpose
     *  =======
     *
     *  DTRMM  performs one of the matrix-matrix operations
     *
     *     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
     *
     *  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     *
     *     op( A ) = A   or   op( A ) = A'.
     *
     *  Arguments
     *  ==========
     *
     *  SIDE   - CHARACTER*1.
     *           On entry,  SIDE specifies whether  op( A ) multiplies B from
     *           the left or right as follows:
     *
     *              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
     *
     *              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
     *
     *           Unchanged on exit.
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix A is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANSA - CHARACTER*1.
     *           On entry, TRANSA specifies the form of op( A ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSA = 'N' or 'n'   op( A ) = A.
     *
     *              TRANSA = 'T' or 't'   op( A ) = A'.
     *
     *              TRANSA = 'C' or 'c'   op( A ) = A'.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit triangular
     *           as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry, M specifies the number of rows of B. M must be at
     *           least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of B.  N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
     *           zero then  A is not referenced and  B need not be set before
     *           entry.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
     *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
     *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
     *           upper triangular part of the array  A must contain the upper
     *           triangular matrix  and the strictly lower triangular part of
     *           A is not referenced.
     *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
     *           lower triangular part of the array  A must contain the lower
     *           triangular matrix  and the strictly upper triangular part of
     *           A is not referenced.
     *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
     *           A  are not referenced either,  but are assumed to be  unity.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
     *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
     *           then LDA must be at least max( 1, n ).
     *           Unchanged on exit.
     *
     *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
     *           Before entry,  the leading  m by n part of the array  B must
     *           contain the matrix  B,  and  on exit  is overwritten  by the
     *           transformed matrix.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in  the  calling  (sub)  program.   LDB  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param side
     *            specifies whether op( A ) multiplies B from the left or right
     * @param uplo
     *            specifies whether the matrix A is an upper or lower triangular
     * @param transa
     *            specifies the form of op( A ) to be used in the matrix
     *            multiplication
     * @param diag
     *            specifies whether or not A is unit triangular
     * @param m
     *            specifies the number of rows of B
     * @param n
     *            specifies the number of columns of B
     * @param alpha
     *            specifies the scalar alpha
     * @param a
     *            array containing matrix A
     * @param _a_offset
     *            offset into the array {@code a}
     * @param lda
     *            the first dimension of A
     * @param b
     *            array containing matrix B
     * @param _b_offset
     *            offset into the array {@code b}
     * @param ldb
     *            the first dimension of B
     */
    public static void dtrmm(Side side, UpLo uplo, Trans transa, Diag diag, int m, int n, double alpha,
            double[] a, int _a_offset, int lda, double[] b, int _b_offset, int ldb) {

        requireNonNull(side, uplo, transa, diag, a, b);

        boolean lside = (Side.LEFT == side); // L
        boolean nounit = (Diag.NON_UNIT == diag); // N
        boolean upper = (UpLo.UPPER == uplo); // U
        boolean notrans = (Trans.NO_TRANS == transa); // N

        if (m < 0) {
            throw new IllegalArgumentException("m < 0");
        } else if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        } else if (lda < Math.max(1, (lside) ? m : n)) {
            throw new IllegalArgumentException("lda < Math.max(1, (Side.LEFT == side) ? m : n");
        } else if (ldb < Math.max(1, m)) {
            throw new IllegalArgumentException("ldb < Math.max(1, m)");
        }

        // Quick return if possible
        if (n == 0) {
            return;
        }

        if ((lside && Math.min(m, n) > 190) || (!lside && Math.min(m, n) > 280)) {
            Trmm.trmm(lside, !upper, !notrans, !nounit, m, n, alpha, a, _a_offset, lda, b, _b_offset, ldb);
        } else {
            DtrmmNetlib.dtrmm(lside, upper, notrans, nounit, m, n, alpha, a, _a_offset, lda, b, _b_offset, ldb);
        }
    }

    private static void requireNonNull(Object... args) {
        for (Object arg : args) {
            Objects.requireNonNull(arg);
        }
    }

    private Dtrmm() {
        throw new AssertionError();
    }
}
