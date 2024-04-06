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

// DSCAL scales a vector by a constant.
final class Dscal {

    static void dscal(int n, double da, double[] dx, int _dx_offset, int incx) {

        if (n <= 0 || incx <= 0) {
            return;
        }

        // code for increment not equal to 1
        if (incx != 1) {
            int nincx = n * incx;
            int l = 1;
            for (int i = ((nincx - 1) + incx) / incx; i > 0; i--) {
                dx[l - 1 + _dx_offset] = da * dx[l - 1 + _dx_offset];
                l += incx;
            }

            return;
        }

        // code for increment equal to 1
        for (int j = 0; j < n; ++j) {
            dx[j + _dx_offset] *= da;
        }
    }
}
