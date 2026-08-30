/*
 * Copyright © ???? The University of Tennessee. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * · Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * · Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer listed in this license in
 *   the documentation and/or other materials provided with the distribution.
 * · Neither the name of the copyright holders nor the names of its contributors
 *   may be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * This software is provided by the copyright holders and contributors "as is" and
 * any express or implied warranties, including, but not limited to, the implied
 * warranties of merchantability and fitness for a particular purpose are disclaimed.
 * In no event shall the copyright owner or contributors be liable for any direct,
 * indirect, incidental, special, exemplary, or consequential damages (including,
 * but not limited to, procurement of substitute goods or services; loss of use,
 * data, or profits; or business interruption) however caused and on any theory of
 * liability, whether in contract, strict liability, or tort (including negligence
 * or otherwise) arising in any way out of the use of this software, even if advised
 * of the possibility of such damage. 
 */
package math.lapack;

import java.util.Objects;

/**
 * DGESV computes the solution to a real system of linear equations
 * {@code A * X = B}, where A is an N-by-N matrix and X and B are N-by-NRHS
 * matrices.
 * <p>
 * The LU decomposition with partial pivoting and row interchanges is used to
 * factor A as {@code A = P * L * U}, where P is a permutation matrix, L is unit
 * lower triangular, and U is upper triangular. The factored form of A is then
 * used to solve the system of equations {@code A * X = B}.
 */
public final class Dgesv {
    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGESV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     *
     *  The LU decomposition with partial pivoting and row interchanges is
     *  used to factor A as
     *     A = P * L * U,
     *  where P is a permutation matrix, L is unit lower triangular, and U is
     *  upper triangular.  The factored form of A is then used to solve the
     *  system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N &ge; 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS &ge; 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the N-by-N coefficient matrix A.
     *          On exit, the factors L and U from the factorization
     *          A = P*L*U; the unit diagonal elements of L are not stored.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA &ge; max(1,N).
     *
     *  IPIV    (output) INTEGER array, dimension (N)
     *          The pivot indices that define the permutation matrix P;
     *          row i of the matrix was interchanged with row IPIV(i).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS matrix of right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB &ge; max(1,N).
     *
     *  INFO    (output) INTEGER
     *          =    0:  successful exit
     *          &lt; 0:  if INFO = -i, the i-th argument had an illegal value
     *          &gt; 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                   has been completed, but the factor U is exactly
     *                   singular, so the solution could not be computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     *            the number of linear equations
     * @param nrhs
     *            the number of right hand sides
     * @param a
     *            On entry, the N-by-N coefficient matrix A. On exit, the
     *            factors L and U from the factorization {@code A = P*L*U}; the
     *            unit diagonal elements of L are not stored.
     * @param _a_offset
     *            offset into the array {@code a}
     * @param lda
     *            the leading dimension of the array A
     * @param ipiv
     *            the pivot indices that define the permutation matrix P
     * @param _ipiv_offset
     *            offset into the array {@code ipiv}
     * @param b
     *            On entry, the N-by-NRHS matrix of right hand side matrix B. On
     *            exit, the N-by-NRHS solution matrix X
     * @param _b_offset
     *            offset into the array {@code b}
     * @param ldb
     *            the leading dimension of the array B
     * @return {@code true} if the factorization has been completed
     *         successfully, {@code false} if the factor {@code U} is exactly
     *         singular, or if {@code A} holds an entry that is not finite, so
     *         that the solution could not be computed
     * @throws NullPointerException
     *             if {@code a}, {@code ipiv} or {@code b} is {@code null}
     * @throws IllegalArgumentException
     *             if a size, an offset or a leading dimension is out of range,
     *             or an array is too short for the problem described
     */
    public static boolean dgesv(int n, int nrhs, double[] a, int _a_offset, int lda, int[] ipiv, int _ipiv_offset,
            double[] b, int _b_offset, int ldb) {

        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");

        if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        } else if (nrhs < 0) {
            throw new IllegalArgumentException("nrhs < 0");
        } else if (lda < Math.max(1, n)) {
            throw new IllegalArgumentException("lda < Math.max(1, n)");
        } else if (ldb < Math.max(1, n)) {
            throw new IllegalArgumentException("ldb < Math.max(1, n)");
        }
        checkOffset(_a_offset, "_a_offset");
        checkOffset(_ipiv_offset, "_ipiv_offset");
        checkOffset(_b_offset, "_b_offset");
        if (ipiv.length - _ipiv_offset < n) {
            throw new IllegalArgumentException("ipiv must hold at least " + n + " indices from offset "
                    + _ipiv_offset + " (length = " + ipiv.length + ")");
        }
        // the tight bound: the last column needs n entries, not lda of them,
        // so a validly packed array is not rejected
        checkLen(a, _a_offset, n > 0 ? (long) lda * (n - 1) + n : 0L, "a");
        checkLen(b, _b_offset, (n > 0 && nrhs > 0) ? (long) ldb * (nrhs - 1) + n : 0L, "b");
        if (!allFinite(a, _a_offset, lda, n, n)) {
            return false;
        }

        intW info_ = new intW(0);
        // Compute the LU factorization of A. A negative info would mean an
        // illegal argument, which the checks above have already ruled out and
        // which Xerbla would throw on in any case, so only the singularity
        // report can come back from here
        Dgetrf.dgetrf(n, n, a, _a_offset, lda, ipiv, _ipiv_offset, info_);
        if (info_.val != 0) {
            // Factor U in the LU decomposition is exactly singular.
            // Solution could not be computed.
            return false;
        }
        // Solve the system A*X = B, overwriting B with X
        Dgetrs.dgetrs("No transpose", n, nrhs, a, _a_offset, lda, ipiv, _ipiv_offset, b, _b_offset, ldb, info_);
        return true;
    }

    private static void checkOffset(int offset, String name) {
        if (offset < 0) {
            throw new IllegalArgumentException(name + " must not be negative, but is " + offset);
        }
    }

    private static void checkLen(double[] array, int offset, long needed, String name) {
        if (array.length - offset < needed) {
            throw new IllegalArgumentException("Array '" + name + "' must hold at least " + needed
                    + " entries from offset " + offset + " (length = " + array.length + ")");
        }
    }

    /**
     * Whether every entry of the {@code m} by {@code n} matrix is finite.
     * <p>
     * {@code Math.abs(v) <= Double.MAX_VALUE} is false for a {@code NaN} and
     * for either infinity and true for every other {@code double}, which makes
     * it one comparison where the two {@code Double} predicates are two. It
     * measured at 2.9 to 10.5 percent of the factorization it guards, against
     * 4.7 to 22.6 percent for {@code isNaN() || isInfinite()}.
     */
    private static boolean allFinite(double[] a, int offset, int lda, int m, int n) {
        for (int j = 0; j < n; ++j) {
            int base = offset + j * lda;
            for (int i = 0; i < m; ++i) {
                if (!(Math.abs(a[base + i]) <= Double.MAX_VALUE)) {
                    return false;
                }
            }
        }
        return true;
    }

    private Dgesv() {
        throw new AssertionError();
    }
}
