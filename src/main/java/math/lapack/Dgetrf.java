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

import math.gemm.Dgemm;
import math.gemm.Trans;

// DGETRF computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
final class Dgetrf {

    static void dgetrf(int m, int n, double[] a, int _a_offset, int lda, int[] ipiv, int _ipiv_offset,
            intW info) {

        info.val = 0;
        if (m < 0) {
            info.val = -1;
        } else if (n < 0) {
            info.val = -2;
        } else if (lda < Math.max(1, m)) {
            info.val = -4;
        }
        if (info.val != 0) {
            Xerbla.xerbla("DGETRF", -info.val);
            return;
        }

        // Quick return if possible
        if (m == 0 || n == 0) {
            return;
        }
        // Determine the block size for this environment
        int nb = Ilaenv.ilaenv(1, "DGETRF", " ", m, n, -1, -1);
        if (nb <= 1 || nb >= Math.min(m, n)) {
            // Use unblocked code
            Dgetf2.dgetf2(m, n, a, _a_offset, lda, ipiv, _ipiv_offset, info);
        } else {
            // Use blocked code
            intW iinfo = new intW(0);
            int j = 1;
            for (int p = ((Math.min(m, n) - 1) + nb) / nb; p > 0; p--) {
                int jb = Math.min((Math.min(m, n) - j) + 1, nb);
                // Factor diagonal and subdiagonal blocks and test
                // for exact singularity
                Dgetf2.dgetf2(m - j + 1, jb, a, j - 1 + (j - 1) * lda + _a_offset, lda, ipiv, j - 1 + _ipiv_offset,
                        iinfo);
                // Adjust INFO and the pivot indices
                if (info.val == 0 && iinfo.val > 0) {
                    info.val = (iinfo.val + j) - 1;
                }
                int i = j;
                for (int q = Math.min(m, j + jb - 1) - j + 1; q > 0; q--) {
                    ipiv[i - 1 + _ipiv_offset] = j - 1 + ipiv[i - 1 + _ipiv_offset];
                    i++;
                }

                // Apply interchanges to columns 1:J-1
                Dlaswp.dlaswp(j - 1, a, _a_offset, lda, j, j + jb - 1, ipiv, _ipiv_offset, 1);

                if (j + jb <= n) {
                    // Apply interchanges to columns J+JB:N
                    Dlaswp.dlaswp(n - j - jb + 1, a, (j + jb - 1) * lda + _a_offset, lda, j, j + jb - 1, ipiv,
                            _ipiv_offset, 1);
                    // Compute block row of U
                    Dtrsm.dtrsm("Left", "Lower", "No transpose", "Unit", jb, n - j - jb + 1, 1.0, a,
                            j - 1 + (j - 1) * lda + _a_offset, lda, a, j - 1 + (j + jb - 1) * lda + _a_offset, lda);

                    if (j + jb <= m) {
                        // Update trailing submatrix
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, m - j - jb + 1, n - j - jb + 1, jb, -1.0, a,
                                j + jb - 1 + (j - 1) * lda + _a_offset, lda, a, j - 1 + (j + jb - 1) * lda + _a_offset,
                                lda, 1.0, a, j + jb - 1 + (j + jb - 1) * lda + _a_offset, lda);
                    }
                }
                j += nb;
            }
        }
    }
}
