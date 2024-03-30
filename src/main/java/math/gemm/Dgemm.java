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
package math.gemm;

import java.util.Objects;

/**
 * Java implementation of BLAS generalized matrix multiplication (gemm) for
 * double-precision Fortran-style matrices (i.e., column-major storage layout
 * matrices holding doubles).
 */
public final class Dgemm {
    /**
     * <pre>
     * <code>
     * Purpose =======
     *
     * DGEMM performs one of the matrix-matrix operations
     *
     * C := alpha*op( A )*op( B ) + beta*C,
     *
     * where op( X ) is one of
     *
     * op( X ) = X or op( X ) = X',
     *
     * alpha and beta are scalars, and A, B and C are matrices, with op( A ) an
     * m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
     *
     * Arguments ==========
     *
     * TRANSA - CHARACTER*1. On entry, TRANSA specifies the form of op( A ) to
     * be used in the matrix multiplication as follows:
     *
     * TRANSA = 'N' or 'n', op( A ) = A.
     *
     * TRANSA = 'T' or 't', op( A ) = A'.
     *
     * TRANSA = 'C' or 'c', op( A ) = A'.
     *
     * Unchanged on exit.
     *
     * TRANSB - CHARACTER*1. On entry, TRANSB specifies the form of op( B ) to
     * be used in the matrix multiplication as follows:
     *
     * TRANSB = 'N' or 'n', op( B ) = B.
     *
     * TRANSB = 'T' or 't', op( B ) = B'.
     *
     * TRANSB = 'C' or 'c', op( B ) = B'.
     *
     * Unchanged on exit.
     *
     * M - INTEGER. On entry, M specifies the number of rows of the matrix op( A
     * ) and of the matrix C. M must be at least zero. Unchanged on exit.
     *
     * N - INTEGER. On entry, N specifies the number of columns of the matrix
     * op( B ) and the number of columns of the matrix C. N must be at least
     * zero. Unchanged on exit.
     *
     * K - INTEGER. On entry, K specifies the number of columns of the matrix
     * op( A ) and the number of rows of the matrix op( B ). K must be at least
     * zero. Unchanged on exit.
     *
     * ALPHA - DOUBLE PRECISION. On entry, ALPHA specifies the scalar alpha.
     * Unchanged on exit.
     *
     * A - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is k when
     * TRANSA = 'N' or 'n', and is m otherwise. Before entry with TRANSA = 'N'
     * or 'n', the leading m by k part of the array A must contain the matrix A,
     * otherwise the leading k by m part of the array A must contain the matrix
     * A. Unchanged on exit.
     *
     * LDA - INTEGER. On entry, LDA specifies the first dimension of A as
     * declared in the calling (sub) program. When TRANSA = 'N' or 'n' then LDA
     * must be at least max( 1, m ), otherwise LDA must be at least max( 1, k ).
     * Unchanged on exit.
     *
     * B - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is n when
     * TRANSB = 'N' or 'n', and is k otherwise. Before entry with TRANSB = 'N'
     * or 'n', the leading k by n part of the array B must contain the matrix B,
     * otherwise the leading n by k part of the array B must contain the matrix
     * B. Unchanged on exit.
     *
     * LDB - INTEGER. On entry, LDB specifies the first dimension of B as
     * declared in the calling (sub) program. When TRANSB = 'N' or 'n' then LDB
     * must be at least max( 1, k ), otherwise LDB must be at least max( 1, n ).
     * Unchanged on exit.
     *
     * BETA - DOUBLE PRECISION. On entry, BETA specifies the scalar beta. When
     * BETA is supplied as zero then C need not be set on input. Unchanged on
     * exit.
     *
     * C - DOUBLE PRECISION array of DIMENSION ( LDC, n ). Before entry, the
     * leading m by n part of the array C must contain the matrix C, except when
     * beta is zero, in which case C need not be set on entry. On exit, the
     * array C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) +
     * beta*C ).
     *
     * LDC - INTEGER. On entry, LDC specifies the first dimension of C as
     * declared in the calling (sub) program. LDC must be at least max( 1, m ).
     * Unchanged on exit.
     *
     *
     * Level 3 Blas routine.
     *
     * -- Written on 8-February-1989. Jack Dongarra, Argonne National
     * Laboratory. Iain Duff, AERE Harwell. Jeremy Du Croz, Numerical Algorithms
     * Group Ltd. Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param transa
     *            transpose matrix {@code A} or not
     * @param transb
     *            transpose matrix {@code B} or not
     * @param m
     *            the number of rows of the matrix op(A) and of the matrix C
     * @param n
     *            the number of columns of the matrix op(B) and the number of
     *            columns of the matrix C
     * @param k
     *            the number of columns of the matrix op(A) and the number of
     *            rows of the matrix op(B)
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
     * @param beta
     *            specifies the scalar beta
     * @param c
     *            array containing matrix C
     * @param _c_offset
     *            offset into the array {@code c}
     * @param ldc
     *            the first dimension of C
     */
    public static void dgemm(Trans transa, Trans transb, int m, int n, int k, double alpha, double[] a, int _a_offset,
            int lda, double[] b, int _b_offset, int ldb, double beta, double[] c, int _c_offset, int ldc) {

        requireNonNull(transa, transb, a, b, c);

        int nRowA = 0;
        int nRowB = 0;
        boolean notA = (transa == Trans.NO_TRANS); // N
        boolean notB = (transb == Trans.NO_TRANS); // N
        if (notA) {
            nRowA = m;
        } else {
            nRowA = k;
        }
        if (notB) {
            nRowB = k;
        } else {
            nRowB = n;
        }

        if (m < 0) {
            throw new IllegalArgumentException("m < 0");
        } else if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        } else if (k < 0) {
            throw new IllegalArgumentException("k < 0");
        } else if (lda < Math.max(1, nRowA)) {
            throw new IllegalArgumentException("lda < Math.max(1, nRowA)");
        } else if (ldb < Math.max(1, nRowB)) {
            throw new IllegalArgumentException("ldb < Math.max(1, nRowB)");
        } else if (ldc < Math.max(1, m)) {
            throw new IllegalArgumentException("ldc < Math.max(1, m)");
        }

        // Quick return if possible
        if ((m == 0 || n == 0) || ((alpha == 0.0 || k == 0) && beta == 1.0)) {
            return;
        }

        /*
         * @param m the number of rows of the matrix op(A) and of the matrix C
         * 
         * @param n the number of columns of the matrix op(B) and the number of
         * columns of the matrix C
         * 
         * @param k the number of columns of the matrix op(A) and the number of
         * rows of the matrix op(B)
         */
        if (Math.round(Math.sqrt(m * n)) <= 220L) {
            DgemmNetlib.dgemm(notA, notB, m, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset,
                    ldc);
        } else {
            DgemmLehn.dgemm(notA, notB, m, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
        }
    }

    private static void requireNonNull(Object... args) {
        for (Object arg : args) {
            Objects.requireNonNull(arg);
        }
    }
}
