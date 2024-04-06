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

// DGEMV performs one of the matrix-vector operations
// y := alpha*A*x + beta*y,  or  y := alpha*A**T*x + beta*y,
// where alpha and beta are scalars, x and y are vectors and
// A is an m by n matrix.
final class Dgemv {

    static void dgemv(String trans, int m, int n, double alpha, double[] a, int _a_offset, int lda, double[] x,
            int _x_offset, int incx, double beta, double[] y, int _y_offset, int incy) {

        int info = 0;
        if ((!Lsame.lsame(trans, "N") && !Lsame.lsame(trans, "T")) && !Lsame.lsame(trans, "C")) {
            info = 1;
        } else if (m < 0) {
            info = 2;
        } else if (n < 0) {
            info = 3;
        } else if (lda < Math.max(1, m)) {
            info = 6;
        } else if (incx == 0) {
            info = 8;
        } else if (incy == 0) {
            info = 11;
        }
        if (info != 0) {
            Xerbla.xerbla("DGEMV ", info);
            return;
        }

        if ((m == 0 || n == 0) || (alpha == 0.0 && beta == 1.0)) {
            return;
        }
        int lenx;
        int leny;
        if (Lsame.lsame(trans, "N")) {
            lenx = n;
            leny = m;
        } else {
            lenx = m;
            leny = n;
        }
        int kx;
        if (incx > 0) {
            kx = 1;
        } else {
            kx = 1 - (lenx - 1) * incx;
        }
        int ky;
        if (incy > 0) {
            ky = 1;
        } else {
            ky = 1 - (leny - 1) * incy;
        }
        // First form y := beta*y
        if (beta != 1.0) {
            if (incy == 1) {
                if (beta == 0.0) {
                    int i = 1;
                    for (int p = leny; p > 0; p--) {
                        y[(i - 1) + _y_offset] = 0.0;
                        i++;
                    }
                } else {
                    int i = 1;
                    for (int p = leny; p > 0; p--) {
                        y[(i - 1) + _y_offset] = beta * y[(i - 1) + _y_offset];
                        i++;
                    }
                }
            } else {
                int iy = ky;
                if (beta == 0.0) {
                    for (int p = leny; p > 0; p--) {
                        y[(iy - 1) + _y_offset] = 0.0;
                        iy += incy;
                    }
                } else {
                    for (int p = leny; p > 0; p--) {
                        y[(iy - 1) + _y_offset] = beta * y[(iy - 1) + _y_offset];
                        iy += incy;
                    }
                }
            }
        }
        if (alpha == 0.0) {
            return;
        }
        // Form y := alpha*A*x + y
        if (Lsame.lsame(trans, "N")) {
            int jx = kx;
            if (incy == 1) {
                int j = 1;
                for (int p = n; p > 0; p--) {
                    if (x[(jx - 1) + _x_offset] != 0.0) {
                        double temp = alpha * x[(jx - 1) + _x_offset];
                        int i = 1;
                        for (int q = m; q > 0; q--) {
                            y[(i - 1) + _y_offset] = y[(i - 1) + _y_offset]
                                    + temp * a[(i - 1) + (j - 1) * lda + _a_offset];
                            i++;
                        }
                    }
                    jx += incx;
                    j++;
                }
            } else {
                int j = 1;
                for (int p = n; p > 0; p--) {
                    if (x[(jx - 1) + _x_offset] != 0.0) {
                        double temp = alpha * x[(jx - 1) + _x_offset];
                        int iy = ky;
                        int i = 1;
                        for (int q = m; q > 0; q--) {
                            y[(iy - 1) + _y_offset] = y[(iy - 1) + _y_offset]
                                    + temp * a[(i - 1) + (j - 1) * lda + _a_offset];
                            iy += incy;
                            i++;
                        }
                    }
                    jx += incx;
                    j++;
                }
            }
        } else {
            // Form y := alpha*A**T*x + y
            int jy = ky;
            if (incx == 1) {
                int j = 1;
                for (int p = n; p > 0; p--) {
                    double temp = 0.0;
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        temp += a[(i - 1) + (j - 1) * lda + _a_offset] * x[(i - 1) + _x_offset];
                        i++;
                    }

                    y[(jy - 1) + _y_offset] = y[(jy - 1) + _y_offset] + alpha * temp;
                    jy += incy;
                    j++;
                }
            } else {
                int j = 1;
                for (int p = n; p > 0; p--) {
                    double temp = 0.0;
                    int ix = kx;
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        temp += a[(i - 1) + (j - 1) * lda + _a_offset] * x[(ix - 1) + _x_offset];
                        ix += incx;
                        i++;
                    }

                    y[(jy - 1) + _y_offset] = y[(jy - 1) + _y_offset] + alpha * temp;
                    jy += incy;
                    j++;
                }
            }
        }
    }
}
