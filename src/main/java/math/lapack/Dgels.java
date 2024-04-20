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

import math.gemm.Trans;

/**
 * DGELS solves overdetermined or underdetermined real linear systems involving
 * an M-by-N matrix A, or its transpose, using a QR or LQ factorization of A. It
 * is assumed that A has full rank.
 * <p>
 * The following options are provided:
 * <p>
 * 1. If TRANS = 'N' and m &ge; n: find the least squares solution of an
 * overdetermined system, i.e., solve the least squares problem
 * {@code minimize || B - A*X ||}.
 * <p>
 * 2. If TRANS = 'N' and m &lt; n: find the minimum norm solution of an
 * underdetermined system {@code A * X = B}.
 * <p>
 * 3. If TRANS = 'T' and m &ge; n: find the minimum norm solution of an
 * underdetermined system {@code A**T * X = B}.
 * <p>
 * 4. If TRANS = 'T' and m &lt; n: find the least squares solution of an
 * overdetermined system, i.e., solve the least squares problem
 * {@code minimize || B - A**T * X ||}.
 * <p>
 * Several right hand side vectors b and solution vectors x can be handled in a
 * single call; they are stored as the columns of the M-by-NRHS right hand side
 * matrix B and the N-by-NRHS solution matrix X.
 */
public final class Dgels {
    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGELS solves overdetermined or underdetermined real linear systems
     *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
     *  factorization of A.  It is assumed that A has full rank.
     *
     *  The following options are provided:
     *
     *  1. If TRANS = 'N' and m &ge; n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A*X ||.
     *
     *  2. If TRANS = 'N' and m &lt; n:  find the minimum norm solution of
     *     an underdetermined system A * X = B.
     *
     *  3. If TRANS = 'T' and m &ge; n:  find the minimum norm solution of
     *     an undetermined system A**T * X = B.
     *
     *  4. If TRANS = 'T' and m &lt; n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A**T * X ||.
     *
     *  Several right hand side vectors b and solution vectors x can be
     *  handled in a single call; they are stored as the columns of the
     *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     *  matrix X.
     *
     *  Arguments
     *  =========
     *
     *  TRANS   (input) CHARACTER*1
     *          = 'N': the linear system involves A;
     *          = 'T': the linear system involves A**T.
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M &ge; 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N &ge; 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of
     *          columns of the matrices B and X. NRHS &ge; 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *            if M &ge; N, A is overwritten by details of its QR
     *                       factorization as returned by DGEQRF;
     *            if M &lt; N, A is overwritten by details of its LQ
     *                       factorization as returned by DGELQF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA &ge; max(1,M).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the matrix B of right hand side vectors, stored
     *          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
     *          if TRANS = 'T'.
     *          On exit, if INFO = 0, B is overwritten by the solution
     *          vectors, stored columnwise:
     *          if TRANS = 'N' and m &ge; n, rows 1 to n of B contain the least
     *          squares solution vectors; the residual sum of squares for the
     *          solution in each column is given by the sum of squares of
     *          elements N+1 to M in that column;
     *          if TRANS = 'N' and m &lt; n, rows 1 to N of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m &ge; n, rows 1 to M of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m &lt; n, rows 1 to M of B contain the
     *          least squares solution vectors; the residual sum of squares
     *          for the solution in each column is given by the sum of
     *          squares of elements M+1 to N in that column.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B. LDB &ge; MAX(1,M,N).
     *
     *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,L
     *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     *
     *  LWORK   (input) INTEGER
     *          The dimension of the array WORK.
     *          LWORK &ge; max( 1, MN + max( MN, NRHS ) ).
     *          For optimal performance,
     *          LWORK &ge; max( 1, MN + max( MN, NRHS )*NB ).
     *          where MN = min(M,N) and NB is the optimum block size.
     *
     *          If LWORK = -1, then a workspace query is assumed; the routine
     *          only calculates the optimal size of the WORK array, returns
     *          this value as the first entry of the WORK array, and no error
     *          message related to LWORK is issued by XERBLA.
     *
     *  INFO    (output) INTEGER
     *          =    0:  successful exit
     *          &lt; 0:  if INFO = -i, the i-th argument had an illegal value
     *          &gt; 0:  if INFO =  i, the i-th diagonal element of the
     *                   triangular factor of A is zero, so that A does not have
     *                   full rank; the least squares solution could not be
     *                   computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param trans
     *            specifies the form of op( A )
     * @param m
     *            the number of rows of the matrix A
     * @param n
     *            the number of columns of the matrix A
     * @param nrhs
     *            the number of right hand sides, i.e., the number of columns of
     *            the matrices B and X
     * @param a
     *            (input/output) array a
     * @param _a_offset
     *            offset into the array {@code a}
     * @param lda
     *            the leading dimension of the array a
     * @param b
     *            (input/output) array b
     * @param _b_offset
     *            offset into the array {@code b}
     * @param ldb
     *            the leading dimension of the array B
     * @param work
     *            (output) workspace array
     * @param _work_offset
     *            offset into the array {@code work}
     * @param lwork
     *            the dimension of the array work
     * @return {@code true} if the computation did succeed,
     *         {@code false otherwise}
     */
    public static boolean dgels(Trans trans, int m, int n, int nrhs, double[] a, int _a_offset, int lda, double[] b,
            int _b_offset, int ldb, double[] work, int _work_offset, int lwork) {

        requireNonNull(trans, a, b, work);

        intW info = new intW(0);
        int mn = Math.min(m, n);
        boolean lquery = (lwork == -1);

        if (m < 0) {
            throw new IllegalArgumentException("m < 0");
        } else if (n < 0) {
            throw new IllegalArgumentException("n < 0");
        } else if (nrhs < 0) {
            throw new IllegalArgumentException("nrhs < 0");
        } else if (lda < Math.max(1, m)) {
            throw new IllegalArgumentException("lda < Math.max(1, m)");
        } else if (ldb < Util.max(1, m, n)) {
            throw new IllegalArgumentException("ldb < max(1, m, n)");
        } else if (!lquery && lwork < Math.max(1, mn + Math.max(mn, nrhs))) {
            info.val = -10;
        }

        // Figure out optimal block size
        int wsize = 0;
        boolean tpsd = false;
        if (info.val == 0 || info.val == -10) {
            tpsd = true;
            if (Trans.NO_TRANS == trans) {
                tpsd = false;
            }
            int nb;
            if (m >= n) {
                nb = Ilaenv.ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
                if (tpsd) {
                    nb = Math.max(nb, Ilaenv.ilaenv(1, "DORMQR", "LN", m, nrhs, n, -1));
                } else {
                    nb = Math.max(nb, Ilaenv.ilaenv(1, "DORMQR", "LT", m, nrhs, n, -1));
                }
            } else {
                nb = Ilaenv.ilaenv(1, "DGELQF", " ", m, n, -1, -1);
                if (tpsd) {
                    nb = Math.max(nb, Ilaenv.ilaenv(1, "DORMLQ", "LT", n, nrhs, m, -1));
                } else {
                    nb = Math.max(nb, Ilaenv.ilaenv(1, "DORMLQ", "LN", n, nrhs, m, -1));
                }
            }
            wsize = Math.max(1, mn + Math.max(mn, nrhs) * nb);
            work[_work_offset] = wsize;
        }
        if (info.val == -10) {
            throw new IllegalArgumentException(
                    "lwork != -1 && lwork < Math.max(1, Math.min(m, n) + Math.max(Math.min(m, n), nrhs))");
        }
        if (lquery) {
            return true;
        }
        // Quick return if possible
        if (Util.min(m, n, nrhs) == 0) {
            Dlaset.dlaset("Full", Math.max(m, n), nrhs, 0.0, 0.0, b, _b_offset, ldb);
            return true;
        }

        // Get machine parameters
        doubleW smlnum = new doubleW(1.0020841800044864E-292);
        doubleW bignum = new doubleW(9.979201547673599E291);
        // Scale A, B if max element outside range [SMLNUM,BIGNUM]
        double[] rwork = new double[1];
        double anrm = Dlange.dlange("M", m, n, a, _a_offset, lda, rwork, 0);
        int iascl = 0;
        if (anrm > 0.0 && anrm < smlnum.val) {
            // Scale matrix norm up to SMLNUM
            Dlascl.dlascl("G", 0, 0, anrm, smlnum.val, m, n, a, _a_offset, lda, info);
            iascl = 1;
        } else if (anrm > bignum.val) {
            // Scale matrix norm down to BIGNUM
            Dlascl.dlascl("G", 0, 0, anrm, bignum.val, m, n, a, _a_offset, lda, info);
            iascl = 2;
        } else if (anrm == 0.0) {
            // Matrix all zero. Return zero solution
            Dlaset.dlaset("F", Math.max(m, n), nrhs, 0.0, 0.0, b, _b_offset, ldb);
            work[_work_offset] = wsize;
            return true;
        }
        int brow = m;
        if (tpsd) {
            brow = n;
        }
        double bnrm = Dlange.dlange("M", brow, nrhs, b, _b_offset, ldb, rwork, 0);
        int ibscl = 0;
        if (bnrm > 0.0 && bnrm < smlnum.val) {
            // Scale matrix norm up to SMLNUM
            Dlascl.dlascl("G", 0, 0, bnrm, smlnum.val, brow, nrhs, b, _b_offset, ldb, info);
            ibscl = 1;
        } else if (bnrm > bignum.val) {
            // Scale matrix norm down to BIGNUM
            Dlascl.dlascl("G", 0, 0, bnrm, bignum.val, brow, nrhs, b, _b_offset, ldb, info);
            ibscl = 2;
        }
        int scllen = 0;
        if (m >= n) {
            // compute QR factorization of A
            Dgeqrf.dgeqrf(m, n, a, _a_offset, lda, work, _work_offset, work, mn + _work_offset, lwork - mn, info);
            // workspace at least N, optimally N*NB
            if (!tpsd) {
                // Least-Squares Problem min || A * X - B ||
                // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
                Dormqr.dormqr("Left", Trans.TRANS, m, nrhs, n, a, _a_offset, lda, work, _work_offset, b, _b_offset, ldb,
                        work, mn + _work_offset, lwork - mn, info);
                // workspace at least NRHS, optimally NRHS*NB
                // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
                Dtrtrs.dtrtrs("Upper", "No transpose", "Non-unit", n, nrhs, a, _a_offset, lda, b, _b_offset, ldb, info);
                if (info.val > 0) {
                    return false;
                }
                scllen = n;
            } else {
                // Underdetermined system of equations A**T * X = B
                // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
                Dtrtrs.dtrtrs("Upper", "Transpose", "Non-unit", n, nrhs, a, _a_offset, lda, b, _b_offset, ldb, info);
                if (info.val > 0) {
                    return false;
                }
                // B(N+1:M,1:NRHS) = 0
                int j = 1;
                for (int p = nrhs; p > 0; p--) {
                    int i = n + 1;
                    for (int q = m - n; q > 0; q--) {
                        b[(i - 1) + (j - 1) * ldb + _b_offset] = 0.0;
                        i++;
                    }

                    j++;
                }

                // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
                Dormqr.dormqr("Left", Trans.NO_TRANS, m, nrhs, n, a, _a_offset, lda, work, _work_offset, b, _b_offset,
                        ldb, work, mn + _work_offset, lwork - mn, info);
                // workspace at least NRHS, optimally NRHS*NB
                scllen = m;
            }
        } else {
            // Compute LQ factorization of A
            Dgelqf.dgelqf(m, n, a, _a_offset, lda, work, _work_offset, work, mn + _work_offset, lwork - mn, info);
            // workspace at least M, optimally M*NB
            if (!tpsd) {
                // underdetermined system of equations A * X = B
                // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
                Dtrtrs.dtrtrs("Lower", "No transpose", "Non-unit", m, nrhs, a, _a_offset, lda, b, _b_offset, ldb, info);
                if (info.val > 0) {
                    return false;
                }
                // B(M+1:N,1:NRHS) = 0
                int j = 1;
                for (int p = nrhs; p > 0; p--) {
                    int i = m + 1;
                    for (int q = n - m; q > 0; q--) {
                        b[(i - 1) + (j - 1) * ldb + _b_offset] = 0.0;
                        i++;
                    }

                    j++;
                }

                // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
                Dormlq.dormlq("Left", "Transpose", n, nrhs, m, a, _a_offset, lda, work, _work_offset, b, _b_offset, ldb,
                        work, mn + _work_offset, lwork - mn, info);
                // workspace at least NRHS, optimally NRHS*NB
                scllen = n;
            } else {
                // overdetermined system min || A**T * X - B ||
                // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
                Dormlq.dormlq("Left", "No transpose", n, nrhs, m, a, _a_offset, lda, work, _work_offset, b, _b_offset,
                        ldb, work, mn + _work_offset, lwork - mn, info);
                // workspace at least NRHS, optimally NRHS*NB
                // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
                Dtrtrs.dtrtrs("Lower", "Transpose", "Non-unit", m, nrhs, a, _a_offset, lda, b, _b_offset, ldb, info);
                if (info.val > 0) {
                    return false;
                }
                scllen = m;
            }
        }

        // Undo scaling
        if (iascl == 1) {
            Dlascl.dlascl("G", 0, 0, anrm, smlnum.val, scllen, nrhs, b, _b_offset, ldb, info);
        } else if (iascl == 2) {
            Dlascl.dlascl("G", 0, 0, anrm, bignum.val, scllen, nrhs, b, _b_offset, ldb, info);
        }
        if (ibscl == 1) {
            Dlascl.dlascl("G", 0, 0, smlnum.val, bnrm, scllen, nrhs, b, _b_offset, ldb, info);
        } else if (ibscl == 2) {
            Dlascl.dlascl("G", 0, 0, bignum.val, bnrm, scllen, nrhs, b, _b_offset, ldb, info);
        }

        work[_work_offset] = wsize;
        if (info.val > 0) {
            return false;
        } else if (info.val < 0) {
            throw new IllegalArgumentException();
        }
        return true;
    }

    private static void requireNonNull(Object... args) {
        for (Object arg : args) {
            Objects.requireNonNull(arg);
        }
    }

    private Dgels() {
        throw new AssertionError();
    }
}
