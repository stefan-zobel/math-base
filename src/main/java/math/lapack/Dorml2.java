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

// DORML2 overwrites the general real m by n matrix C with
//
//       Q * C  if SIDE = 'L' and TRANS = 'N', or
//
//       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
//
//       C * Q  if SIDE = 'R' and TRANS = 'N', or
//
//       C * Q**T if SIDE = 'R' and TRANS = 'T',
//
// where Q is a real orthogonal matrix defined as the product
// of k elementary reflectors
//
// Q = H(k) . . . H(2) H(1)
//
// as returned by DGELQF. Q is of order m if SIDE = 'L' and
// of order n if SIDE = 'R'.
final class Dorml2 {

    static void dorml2(String side, String trans, int m, int n, int k, double[] a, int _a_offset, int lda,
            double[] tau, int _tau_offset, double[] c, int _c_offset, int ldc, double[] work, int _work_offset,
            intW info) {

        info.val = 0;
        boolean left = Lsame.lsame(side, "L");
        boolean notran = Lsame.lsame(trans, "N");
        // NQ is the order of Q
        int nq;
        if (left) {
            nq = m;
        } else {
            nq = n;
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
        }

        if (info.val != 0) {
            Xerbla.xerbla("DORML2", -info.val);
            return;
        }
        // Quick return if possible
        if ((m == 0 || n == 0) || k == 0) {
            return;
        }

        int i1;
        int i2;
        int i3;
        if ((left && notran) || (!left && !notran)) {
            i1 = 1;
            i2 = k;
            i3 = 1;
        } else {
            i1 = k;
            i2 = 1;
            i3 = -1;
        }

        int ni = 0;
        int jc = 0;
        int mi = 0;
        int ic = 0;
        if (left) {
            ni = n;
            jc = 1;
        } else {
            mi = m;
            ic = 1;
        }

        int i = i1;
        for (int p = (i2 - i1 + i3) / i3; p > 0; p--) {
            if (left) {
                // H(i) is applied to C(i:m,1:n)
                mi = m - i + 1;
                ic = i;
            } else {
                // H(i) is applied to C(1:m,i:n)
                ni = n - i + 1;
                jc = i;
            }
            double aii = a[i - 1 + (i - 1) * lda + _a_offset];
            a[i - 1 + (i - 1) * lda + _a_offset] = 1.0;
            Dlarf.dlarf(side, mi, ni, a, i - 1 + (i - 1) * lda + _a_offset, lda, tau[i - 1 + _tau_offset], c,
                    ic - 1 + (jc - 1) * ldc + _c_offset, ldc, work, _work_offset);
            a[i - 1 + (i - 1) * lda + _a_offset] = aii;
            i += i3;
        }
    }
}
