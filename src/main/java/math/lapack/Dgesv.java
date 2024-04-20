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
     *         singular, so the solution could not be computed
     */
    public static boolean dgesv(int n, int nrhs, double[] a, int _a_offset, int lda, int[] ipiv, int _ipiv_offset,
            double[] b, int _b_offset, int ldb) {

        requireNonNull(a, ipiv, b);
        checkMinLen(ipiv, n, "ipiv");

        if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        } else if (nrhs < 0) {
            throw new IllegalArgumentException("nrhs < 0");
        } else if (lda < Math.max(1, n)) {
            throw new IllegalArgumentException("lda < Math.max(1, n)");
        } else if (ldb < Math.max(1, n)) {
            throw new IllegalArgumentException("ldb < Math.max(1, n)");
        }

        intW info_ = new intW(0);
        // Compute the LU factorization of A
        Dgetrf.dgetrf(n, n, a, _a_offset, lda, ipiv, _ipiv_offset, info_);
        if (info_.val < 0) {
            throwIllegalArg(info_.val);
        } else if (info_.val == 0) {
            // Solve the system A*X = B, overwriting B with X
            Dgetrs.dgetrs("No transpose", n, nrhs, a, _a_offset, lda, ipiv, _ipiv_offset, b, _b_offset, ldb, info_);
            if (info_.val < 0) {
                throwIllegalArg(info_.val);
            } else if (info_.val > 0) {
                // Factor U in the LU decomposition is exactly singular.
                // Solution could not be computed.
                return false;
            }
            return true;
        }
        return false;
    }

    private static void checkMinLen(int[] array, int minLen, String name) {
        if (array.length < minLen) {
            throw new IllegalArgumentException("Length of array '" + name + "' argument must be at least " + minLen
                    + " (length = " + array.length + ")");
        }
    }

    private static void throwIllegalArg(int pos) {
        throw new IllegalArgumentException("Illegal argument at position: " + pos);
    }

    private static void requireNonNull(Object... args) {
        for (Object arg : args) {
            Objects.requireNonNull(arg);
        }
    }

    private Dgesv() {
        throw new AssertionError();
    }
}
