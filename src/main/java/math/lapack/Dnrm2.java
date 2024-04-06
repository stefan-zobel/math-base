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

// DNRM2 returns the euclidean norm of a vector via the function
// name, so that DNRM2 := sqrt( x'*x )
final class Dnrm2 {

    static double dnrm2(int n, double[] x, int _x_offset, int incx) {

        double norm = 0.0;

        if (n < 1 || incx < 1) {
            norm = 0.0;
        } else if (n == 1) {
            norm = Math.abs(x[_x_offset]);
        } else {
            double scale = 0.0;
            double ssq = 1.0;
            int ix = 1;
            for (int i = (((1 + (n - 1) * incx) - 1) + incx) / incx; i > 0; i--) {
                if (x[(ix - 1) + _x_offset] != 0.0) {
                    double absxi = Math.abs(x[(ix - 1) + _x_offset]);
                    if (scale < absxi) {
                        ssq = 1.0 + ssq * Math.pow(scale / absxi, 2);
                        scale = absxi;
                    } else {
                        ssq += Math.pow(absxi / scale, 2);
                    }
                }
                ix += incx;
            }

            norm = scale * Math.sqrt(ssq);
        }
        return norm;
    }
}
