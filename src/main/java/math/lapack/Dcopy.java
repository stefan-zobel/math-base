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

// DCOPY copies a vector, x, to a vector, y.
final class Dcopy {

    static void dcopy(int n, double[] dx, int _dx_offset, int incx, double[] dy, int _dy_offset, int incy) {

        if (n <= 0) {
            return;
        }

        // code for unequal increments or equal increments
        // not equal to 1
        if (incx != 1 || incy != 1) {
            int ix = 1;
            int iy = 1;
            if (incx < 0) {
                ix = (-n + 1) * incx + 1;
            }
            if (incy < 0) {
                iy = (-n + 1) * incy + 1;
            }
            for (int i = n; i > 0; i--) {
                dy[iy - 1 + _dy_offset] = dx[ix - 1 + _dx_offset];
                ix += incx;
                iy += incy;
            }

            return;
        }

        // code for both increments equal to 1
        System.arraycopy(dx, _dx_offset, dy, _dy_offset, n);
    }
}
