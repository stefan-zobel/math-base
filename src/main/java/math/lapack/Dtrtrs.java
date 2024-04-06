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

// DTRTRS solves a triangular system of the form
//   A * X = B  or  A**T * X = B,
// where A is a triangular matrix of order N, and B is an N-by-NRHS
// matrix.  A check is made to verify that A is nonsingular.
final class Dtrtrs {

    static void dtrtrs(String uplo, String trans, String diag, int n, int nrhs, double[] a, int _a_offset,
            int lda, double[] b, int _b_offset, int ldb, intW info) {

        info.val = 0;
        boolean nounit = Lsame.lsame(diag, "N");
        if (!Lsame.lsame(uplo, "U") && !Lsame.lsame(uplo, "L")) {
            info.val = -1;
        } else if (!Lsame.lsame(trans, "N") && !Lsame.lsame(trans, "T") && !Lsame.lsame(trans, "C")) {
            info.val = -2;
        } else if (!nounit && !Lsame.lsame(diag, "U")) {
            info.val = -3;
        } else if (n < 0) {
            info.val = -4;
        } else if (nrhs < 0) {
            info.val = -5;
        } else if (lda < Math.max(1, n)) {
            info.val = -7;
        } else if (ldb < Math.max(1, n)) {
            info.val = -9;
        }
        if (info.val != 0) {
            Xerbla.xerbla("DTRTRS", -info.val);
            return;
        }
        if (n == 0) {
            return;
        }
        if (nounit) {
            info.val = 1;
            for (int p = n; p > 0; p--) {
                if (a[(info.val - 1) + (info.val - 1) * lda + _a_offset] == 0.0) {
                    return;
                }
                info.val = info.val + 1;
            }
        }
        info.val = 0;
        // Solve A * x = b or A**T * x = b
        Dtrsm.dtrsm("Left", uplo, trans, diag, n, nrhs, 1.0, a, _a_offset, lda, b, _b_offset, ldb);
    }
}
