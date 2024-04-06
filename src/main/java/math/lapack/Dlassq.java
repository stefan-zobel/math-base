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

// DLASSQ  returns the values  scl and smsq such that
//
//    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
//
// where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is
// assumed to be non-negative and scl returns the value
//
//    scl = max( scale, abs( x( i ) ) ).
//
// scale and sumsq must be supplied in SCALE and SUMSQ and
// scl and smsq are overwritten on SCALE and SUMSQ respectively.
//
// The routine makes only one pass through the vector x.
final class Dlassq {

    static void dlassq(int n, double[] x, int _x_offset, int incx, doubleW scale, doubleW sumsq) {

        if (n > 0) {
            int ix = 1;
            for (int p = ((n - 1) * incx + incx) / incx; p > 0; p--) {
                if (x[ix - 1 + _x_offset] != 0.0) {
                    double absxi = Math.abs(x[ix - 1 + _x_offset]);
                    if (scale.val < absxi) {
                        sumsq.val = 1.0 + sumsq.val * Math.pow(scale.val / absxi, 2.0);
                        scale.val = absxi;
                    } else {
                        sumsq.val = sumsq.val + Math.pow(absxi / scale.val, 2.0);
                    }
                }
                ix += incx;
            }
        }
    }
}
