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

// DGETF2 computes an LU factorization of a general m-by-n matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
final class Dgetf2 {

    static void dgetf2(int m, int n, double[] a, int _a_offset, int lda, int[] ipiv, int _ipiv_offset,
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
            Xerbla.xerbla("DGETF2", -info.val);
            return;
        }
        // Quick return if possible
        if (m == 0 || n == 0) {
            return;
        }

        int j = 1;
        for (int p = Math.min(m, n); p > 0; p--) {
            // Find pivot and test for singularity
            int jp = j - 1 + Idamax.idamax(m - j + 1, a, j - 1 + (j - 1) * lda + _a_offset, 1);
            ipiv[j - 1 + _ipiv_offset] = jp;
            if (a[jp - 1 + (j - 1) * lda + _a_offset] != 0.0) {
                // Apply the interchange to columns 1:N
                if (jp != j) {
                    Dswap.dswap(n, a, j - 1 + _a_offset, lda, a, jp - 1 + _a_offset, lda);
                }
                // Compute elements J+1:M of J-th column
                if (j < m) {
                    if (Math.abs(a[j - 1 + (j - 1) * lda + _a_offset]) >= 2.2250738585072014E-308) {
                        Dscal.dscal(m - j, 1.0 / a[j - 1 + (j - 1) * lda + _a_offset], a, j + (j - 1) * lda + _a_offset,
                                1);
                    } else {
                        int i = 1;
                        for (int q = m - j; q > 0; q--) {
                            a[j + i - 1 + (j - 1) * lda + _a_offset] = a[j + i - 1 + (j - 1) * lda + _a_offset]
                                    / a[j - 1 + (j - 1) * lda + _a_offset];
                            i++;
                        }
                    }
                }
            } else if (info.val == 0) {
                info.val = j;
            }
            if (j < Math.min(m, n)) {
                // Update trailing submatrix
                Dger.dger(m - j, n - j, -1.0, a, j + (j - 1) * lda + _a_offset, 1, a, j - 1 + j * lda + _a_offset, lda,
                        a, j + j * lda + _a_offset, lda);
            }
            j++;
        }
    }
}
