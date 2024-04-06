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

import math.gemm.Trans;
import math.trmm.Side;

// DORMLQ overwrites the general real M-by-N matrix C with
//
//                 SIDE = 'L'     SIDE = 'R'
// TRANS = 'N':      Q * C          C * Q
// TRANS = 'T':      Q**T * C       C * Q**T
//
// where Q is a real orthogonal matrix defined as the product
// of k elementary reflectors
//
// Q = H(k) . . . H(2) H(1)
//
// as returned by DGELQF. Q is of order M if SIDE = 'L' and
// of order N if SIDE = 'R'.
final class Dormlq {

    static void dormlq(String side, String trans, int m, int n, int k, double[] a, int _a_offset, int lda,
            double[] tau, int _tau_offset, double[] c, int _c_offset, int ldc, double[] work, int _work_offset,
            int lwork, intW info) {

        info.val = 0;
        boolean left = Lsame.lsame(side, "L");
        boolean notran = Lsame.lsame(trans, "N");
        boolean lquery = (lwork == -1);
        // NQ is the order of Q and NW is the minimum dimension of WORK
        int nq = 0;
        int nw = 0;
        if (left) {
            nq = m;
            nw = n;
        } else {
            nq = n;
            nw = m;
        }
        if (!left && !Lsame.lsame(side, "R")) {
            info.val = -1;
        } else if (!notran && !Lsame.lsame(trans, "T")) {
            info.val = -2;
        } else if (m < 0) {
            info.val = -3;
        } else if (n < 0) {
            info.val = -4;
        } else if (k < 0 || k > nq) {
            info.val = -5;
        } else if (lda < Math.max(1, k)) {
            info.val = -7;
        } else if (ldc < Math.max(1, m)) {
            info.val = -10;
        } else if (!lquery && lwork < Math.max(1, nw)) {
            info.val = -12;
        }

        int lwkopt = 0;
        int nb = 0;
        if (info.val == 0) {
            // Compute the workspace requirements
            nb = Math.min(nbmax, Ilaenv.ilaenv(1, "DORMLQ", side + trans, m, n, k, -1));
            lwkopt = Math.max(1, nw) * nb;
            work[_work_offset] = lwkopt;
        }
        if (info.val != 0) {
            Xerbla.xerbla("DORMLQ", -info.val);
            return;
        }
        if (lquery) {
            return;
        }
        // Quick return if possible
        if ((m == 0 || n == 0) || k == 0) {
            work[_work_offset] = 1;
            return;
        }
        int nbmin = 2;
        int ldwork = nw;
        if (nb > 1 && nb < k) {
            if (lwork < nw * nb) {
                nb = lwork / ldwork;
                nbmin = Math.max(2, Ilaenv.ilaenv(2, "DORMLQ", side + trans, m, n, k, -1));
            }
        }

        if (nb < nbmin || nb >= k) {
            // Use unblocked code
            Dorml2.dorml2(side, trans, m, n, k, a, _a_offset, lda, tau, _tau_offset, c, _c_offset, ldc, work,
                    _work_offset, refInfo);
        } else {
            // Use blocked code
            int ni = 0;
            int jc = 0;

            int mi = 0;
            int ic = 0;

            int i1;
            int i2;
            int i3;
            if ((left && notran) || (!left && !notran)) {
                i1 = 1;
                i2 = k;
                i3 = nb;
            } else {
                i1 = ((k - 1) / nb) * nb + 1;
                i2 = 1;
                i3 = -nb;
            }
            if (left) {
                ni = n;
                jc = 1;
            } else {
                mi = m;
                ic = 1;
            }
            Trans transt = notran ? Trans.TRANS : Trans.NO_TRANS;
            double[] buffer = new double[(nbmax + 1) * nbmax];
            int i = i1;
            for (int p = ((i2 - i1) + i3) / i3; p > 0; p--) {
                int ib = Math.min(nb, k - i + 1);
                // Form the triangular factor of the block reflector
                // H = H(i) H(i+1) . . . H(i+ib-1)
                Dlarft.dlarft("Forward", "Rowwise", nq - i + 1, ib, a, i - 1 + (i - 1) * lda + _a_offset, lda, tau,
                        i - 1 + _tau_offset, buffer, 0, (nbmax + 1));
                if (left) {
                    // H or H**T is applied to C(i:m,1:n)
                    mi = m - i + 1;
                    ic = i;
                } else {
                    // H or H**T is applied to C(1:m,i:n)
                    ni = n - i + 1;
                    jc = i;
                }
                // Apply H or H**T
                Side side_ = left ? Side.LEFT : Side.RIGHT;
                Dlarfb.dlarfb(side_, transt, "Forward", "Rowwise", mi, ni, ib, a, i - 1 + (i - 1) * lda + _a_offset, lda,
                        buffer, 0, (nbmax + 1), c, ic - 1 + (jc - 1) * ldc + _c_offset, ldc, work, _work_offset, ldwork);
                i += i3;
            }
        }
        work[_work_offset] = lwkopt;
    }

    private static final int nbmax = 64;
    private static final intW refInfo = new intW(0);
}
