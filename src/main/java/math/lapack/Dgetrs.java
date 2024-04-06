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

// DGETRS solves a system of linear equations
//    A * X = B  or  A**T * X = B
// with a general N-by-N matrix A using the LU factorization
// computed by DGETRF.
final class Dgetrs {

    static void dgetrs(String trans, int n, int nrhs, double[] a, int _a_offset, int lda, int[] ipiv,
            int _ipiv_offset, double[] b, int _b_offset, int ldb, intW info) {

        info.val = 0;
        boolean notran = Lsame.lsame(trans, "N");
        if (!notran && !Lsame.lsame(trans, "T") && !Lsame.lsame(trans, "C")) {
            info.val = -1;
        } else if (n < 0) {
            info.val = -2;
        } else if (nrhs < 0) {
            info.val = -3;
        } else if (lda < Math.max(1, n)) {
            info.val = -5;
        } else if (ldb < Math.max(1, n)) {
            info.val = -8;
        }

        if (info.val != 0) {
            Xerbla.xerbla("DGETRS", -info.val);
            return;
        }
        // Quick return if possible
        if (n == 0 || nrhs == 0) {
            return;
        }

        if (notran) {
            // Solve A * X = B
            // Apply row interchanges to the right hand sides
            Dlaswp.dlaswp(nrhs, b, _b_offset, ldb, 1, n, ipiv, _ipiv_offset, 1);
            // Solve L*X = B, overwriting B with X
            Dtrsm.dtrsm("Left", "Lower", "No transpose", "Unit", n, nrhs, 1.0, a, _a_offset, lda, b, _b_offset, ldb);
            // Solve U*X = B, overwriting B with X
            Dtrsm.dtrsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, 1.0, a, _a_offset, lda, b, _b_offset,
                    ldb);
        } else {
            // Solve A**T * X = B
            // Solve U**T *X = B, overwriting B with X
            Dtrsm.dtrsm("Left", "Upper", "Transpose", "Non-unit", n, nrhs, 1.0, a, _a_offset, lda, b, _b_offset, ldb);
            // Solve L**T *X = B, overwriting B with X
            Dtrsm.dtrsm("Left", "Lower", "Transpose", "Unit", n, nrhs, 1.0, a, _a_offset, lda, b, _b_offset, ldb);
            // Apply row interchanges to the solution vectors
            Dlaswp.dlaswp(nrhs, b, _b_offset, ldb, 1, n, ipiv, _ipiv_offset, -1);
        }
    }
}
