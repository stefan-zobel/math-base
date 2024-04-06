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

// DLARFG generates a real elementary reflector H of order n, such
// that
//
//       H * ( alpha ) = ( beta ),   H**T * H = I.
//           (   x   )   (   0  )
//
// where alpha and beta are scalars, and x is an (n-1)-element real
// vector. H is represented in the form
//
//       H = I - tau * ( 1 ) * ( 1 v**T ) ,
//                     ( v )
//
// where tau is a real scalar and v is a real (n-1)-element
// vector.
//
// If the elements of x are all zero, then tau = 0 and H is taken
// to be the unit matrix.
//
// Otherwise  1 <= tau <= 2.
final class Dlarfg {

    static void dlarfg(int n, doubleW alpha, double[] x, int _x_offset, int incx, doubleW tau) {

        if (n <= 1) {
            tau.val = 0.0;
            return;
        }
        double xnorm = Dnrm2.dnrm2(n - 1, x, _x_offset, incx);
        if (xnorm == 0.0) {
            // H = I
            tau.val = 0.0;
        } else {
            // general case
            double beta = -Util.dsign(Dlapy2.dlapy2(alpha.val, xnorm), alpha.val);
            final double safmin = 2.0041683600089728E-292;
            if (Math.abs(beta) < safmin) {
                final double rsafmn = 4.9896007738367995E291;
                int knt = 0;
                do {
                    knt++;
                    Dscal.dscal(n - 1, rsafmn, x, _x_offset, incx);
                    beta *= rsafmn;
                    alpha.val = alpha.val * rsafmn;
                } while (Math.abs(beta) < safmin);

                // New BETA is at most 1, at least SAFMIN

                xnorm = Dnrm2.dnrm2(n - 1, x, _x_offset, incx);
                beta = -Util.dsign(Dlapy2.dlapy2(alpha.val, xnorm), alpha.val);
                tau.val = (beta - alpha.val) / beta;
                Dscal.dscal(n - 1, 1.0 / (alpha.val - beta), x, _x_offset, incx);
                alpha.val = beta;
                // If ALPHA is subnormal, it may lose relative accuracy
                for (int j = knt; j > 0; j--) {
                    alpha.val = alpha.val * safmin;
                }
            } else {
                tau.val = (beta - alpha.val) / beta;
                Dscal.dscal(n - 1, 1.0 / (alpha.val - beta), x, _x_offset, incx);
                alpha.val = beta;
            }
        }
    }
}
