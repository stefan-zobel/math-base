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

// DLASWP performs a series of row interchanges on the matrix A.
// One row interchange is initiated for each of rows K1 through K2 of A.
final class Dlaswp {

    static void dlaswp(int n, double[] a, int _a_offset, int lda, int k1, int k2, int[] ipiv, int _ipiv_offset,
            int incx) {

        int ix0;
        int i1;
        int i2;
        int inc;

        // Interchange row I with row IPIV(K1+(I-K1)*abs(INCX))
        // for each of rows K1 through K2

        if (incx > 0) {
            ix0 = k1;
            i1 = k1;
            i2 = k2;
            inc = 1;
        } else if (incx < 0) {
            ix0 = 1 + (1 - k2) * incx;
            i1 = k2;
            i2 = k1;
            inc = -1;
        } else {
            return;
        }

        int n32 = (n / 32) * 32;
        if (n32 != 0) {
            int j = 1;
            for (int p = (n32 + 31) / 32; p > 0; p--) {
                int ix = ix0;
                int i = i1;
                for (int q = (i2 - i1 + inc) / inc; q > 0; q--) {
                    int ip = ipiv[ix - 1 + _ipiv_offset];
                    if (ip != i) {
                        int k = j;
                        for (int r = 32; r > 0; r--) {
                            double temp = a[i - 1 + (k - 1) * lda + _a_offset];
                            a[i - 1 + (k - 1) * lda + _a_offset] = a[ip - 1 + (k - 1) * lda + _a_offset];
                            a[ip - 1 + (k - 1) * lda + _a_offset] = temp;
                            k++;
                        }
                    }
                    ix += incx;
                    i += inc;
                }

                j += 32;
            }

        }

        if (n32 != n) {
            n32++;
            int ix = ix0;
            int i = i1;
            for (int p = (i2 - i1 + inc) / inc; p > 0; p--) {
                int ip = ipiv[ix - 1 + _ipiv_offset];
                if (ip != i) {
                    int k = n32;
                    for (int q = n - n32 + 1; q > 0; q--) {
                        double temp = a[i - 1 + (k - 1) * lda + _a_offset];
                        a[i - 1 + (k - 1) * lda + _a_offset] = a[ip - 1 + (k - 1) * lda + _a_offset];
                        a[ip - 1 + (k - 1) * lda + _a_offset] = temp;
                        k++;
                    }
                }
                ix += incx;
                i += inc;
            }
        }
    }
}
