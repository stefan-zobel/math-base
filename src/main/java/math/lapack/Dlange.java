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

// DLANGE  returns the value of the one norm,  or the Frobenius norm, or
// the  infinity norm,  or the  element of  largest absolute value  of a
// real matrix A.
//
//    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
//    (
//    ( norm1(A),         NORM = '1', 'O' or 'o'
//    (
//    ( normI(A),         NORM = 'I' or 'i'
//    (
//    ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
//
// where norm1 denotes the  one norm of a matrix (maximum column sum),
// normI denotes the infinity norm  of a matrix (maximum row sum) and
// normF denotes the Frobenius norm of a matrix (square root of sum of
// squares). Note that max(abs(A(i,j))) is not a consistent matrix norm.
final class Dlange {

    static double dlange(String norm, int m, int n, double[] a, int _a_offset, int lda, double[] work,
            int _work_offset) {

        double value = 0.0;

        if (Math.min(m, n) == 0) {
            value = 0.0;
        } else if (Lsame.lsame(norm, "M")) {
            // Find max(abs(A(i,j)))
            value = 0.0;
            int j = 1;
            for (int p = n; p > 0; p--) {
                int i = 1;
                for (int q = m; q > 0; q--) {
                    value = Math.max(value, Math.abs(a[i - 1 + (j - 1) * lda + _a_offset]));
                    i++;
                }

                j++;
            }

        } else if (Lsame.lsame(norm, "O") || norm.regionMatches(0, "1", 0, 1)) {
            // Find norm1(A)
            value = 0.0;
            int j = 1;
            for (int p = n; p > 0; p--) {
                double sum = 0.0;
                int i = 1;
                for (int q = m; q > 0; q--) {
                    sum += Math.abs(a[i - 1 + (j - 1) * lda + _a_offset]);
                    i++;
                }

                value = Math.max(value, sum);
                j++;
            }

        } else if (Lsame.lsame(norm, "I")) {
            // Find normI(A)
            int i = 1;
            for (int p = m; p > 0; p--) {
                work[i - 1 + _work_offset] = 0.0;
                i++;
            }

            int j = 1;
            for (int p = n; p > 0; p--) {
                i = 1;
                for (int q = m; q > 0; q--) {
                    work[i - 1 + _work_offset] = work[i - 1 + _work_offset]
                            + Math.abs(a[i - 1 + (j - 1) * lda + _a_offset]);
                    i++;
                }

                j++;
            }

            value = 0.0;
            i = 1;
            for (int p = m; p > 0; p--) {
                value = Math.max(value, work[i - 1 + _work_offset]);
                i++;
            }

        } else if (Lsame.lsame(norm, "F") || Lsame.lsame(norm, "E")) {
            // Find normF(A)
            doubleW scale = new doubleW(0.0);
            doubleW sum = new doubleW(1.0);

            int j = 1;
            for (int p = n; p > 0; p--) {
                Dlassq.dlassq(m, a, (j - 1) * lda + _a_offset, 1, scale, sum);
                j++;
            }

            value = scale.val * Math.sqrt(sum.val);
        }
        return value;
    }
}
