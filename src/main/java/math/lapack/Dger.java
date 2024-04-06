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

// DGER performs the rank 1 operation A := alpha*x*y**T + A,
// where alpha is a scalar, x is an m element vector, y is
// an n element vector and A is an m by n matrix.
final class Dger {

    static void dger(int m, int n, double alpha, double[] x, int _x_offset, int incx, double[] y, int _y_offset,
            int incy, double[] a, int _a_offset, int lda) {

        int info = 0;
        if (m < 0) {
            info = 1;
        } else if (n < 0) {
            info = 2;
        } else if (incx == 0) {
            info = 5;
        } else if (incy == 0) {
            info = 7;
        } else if (lda < Math.max(1, m)) {
            info = 9;
        }
        if (info != 0) {
            Xerbla.xerbla("DGER  ", info);
            return;
        }
        if ((m == 0 || n == 0) || alpha == 0.0) {
            return;
        }

        int jy = 0;
        if (incy > 0) {
            jy = 1;
        } else {
            jy = 1 - (n - 1) * incy;
        }
        if (incx == 1) {
            int j = 1;
            for (int p = n; p > 0; p--) {
                if (y[jy - 1 + _y_offset] != 0.0) {
                    double temp = alpha * y[jy - 1 + _y_offset];
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset]
                                + x[i - 1 + _x_offset] * temp;
                        i++;
                    }
                }
                jy += incy;
                j++;
            }
        } else {
            int kx;
            if (incx > 0) {
                kx = 1;
            } else {
                kx = 1 - (m - 1) * incx;
            }
            int j = 1;
            for (int p = n; p > 0; p--) {
                if (y[jy - 1 + _y_offset] != 0.0) {
                    double temp = alpha * y[jy - 1 + _y_offset];
                    int ix = kx;
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset]
                                + x[ix - 1 + _x_offset] * temp;
                        ix += incx;
                        i++;
                    }

                }
                jy += incy;
                j++;
            }
        }
    }
}
