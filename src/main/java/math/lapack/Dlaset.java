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

// DLASET initializes an m-by-n matrix A to BETA on the diagonal and
// ALPHA on the offdiagonals.
final class Dlaset {

    static void dlaset(String uplo, int m, int n, double alpha, double beta, double[] a, int _a_offset,
            int lda) {

        int k = 0;
        // Set the strictly upper triangular or trapezoidal part of the
        // array to ALPHA.
        if (Lsame.lsame(uplo, "U")) {
            int l = 2;
            for (int j = n - 1; j > 0; j--) {
                k = 1;
                for (int i = Math.min(l - 1, m); i > 0; i--) {
                    a[k - 1 + (l - 1) * lda + _a_offset] = alpha;
                    k++;
                }

                l++;
            }

        } else if (Lsame.lsame(uplo, "L")) {
            // Set the strictly lower triangular or trapezoidal part of the
            // array to ALPHA.
            int l = 1;
            for (int j = Math.min(m, n); j > 0; j--) {
                k = l + 1;
                for (int i = m - l; i > 0; i--) {
                    a[k - 1 + (l - 1) * lda + _a_offset] = alpha;
                    k++;
                }

                l++;
            }

        } else {
            // Set the leading m-by-n submatrix to ALPHA.
            int l = 1;
            for (int j = n; j > 0; j--) {
                k = 1;
                for (int i = m; i > 0; i--) {
                    a[k - 1 + (l - 1) * lda + _a_offset] = alpha;
                    k++;
                }

                l++;
            }

        }
        // Set the first min(M,N) diagonal elements to BETA.
        k = 1;
        for (int i = Math.min(m, n); i > 0; i--) {
            a[k - 1 + (k - 1) * lda + _a_offset] = beta;
            k++;
        }
    }
}
