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

// DSWAP interchanges two vectors.
// Uses unrolled loops for increments equal to 1.
final class Dswap {

    static void dswap(int n, double[] dx, int _dx_offset, int incx, double[] dy, int _dy_offset, int incy) {

        if (n <= 0) {
            return;
        }

        // code for unequal increments or equal increments
        // not equal to 1
        if (incx != 1 || incy != 1) {
            int ix = 1;
            int iy = 1;
            if (incx < 0) {
                ix = (1 - n) * incx + 1;
            }
            if (incy < 0) {
                iy = (1 - n) * incy + 1;
            }
            for (int i = n; i > 0; i--) {
                double dtemp = dx[ix - 1 + _dx_offset];
                dx[ix - 1 + _dx_offset] = dy[iy - 1 + _dy_offset];
                dy[iy - 1 + _dy_offset] = dtemp;
                ix += incx;
                iy += incy;
            }

            return;
        }

        // code for both increments equal to 1
        int m = n % 3;
        int i;
        if (m != 0) {
            i = 1;
            for (int k = m; k > 0; k--) {
                double dtemp = dx[i - 1 + _dx_offset];
                dx[i - 1 + _dx_offset] = dy[i - 1 + _dy_offset];
                dy[i - 1 + _dy_offset] = dtemp;
                i++;
            }

            if (n < 3) {
                return;
            }
        }

        int mp1 = m + 1;
        i = mp1;

        for (int k = (n - mp1 + 3) / 3; k > 0; k--) {

            double dtemp           = dx[i - 1 + _dx_offset];
            dx[i - 1 + _dx_offset] = dy[i - 1 + _dy_offset];
            dy[i - 1 + _dy_offset] = dtemp;

            dtemp                  = dx[i +     _dx_offset];
            dx[i +     _dx_offset] = dy[i +     _dy_offset];
            dy[i +     _dy_offset] = dtemp;

            dtemp                  = dx[i + 1 + _dx_offset];
            dx[i + 1 + _dx_offset] = dy[i + 1 + _dy_offset];
            dy[i + 1 + _dy_offset] = dtemp;

            i += 3;
        }
    }
}
