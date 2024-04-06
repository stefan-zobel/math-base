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

// DLARF applies a real elementary reflector H to a
// real m by n matrix C, from either the left or the
// right. H is represented in the form
//
//       H = I - tau * v * v**T
//
// where tau is a real scalar and v is a real vector.
//
// If tau = 0, then H is taken to be the unit matrix.
final class Dlarf {

    static void dlarf(String side, int m, int n, double[] v, int _v_offset, int incv, double tau, double[] c,
            int _c_offset, int ldc, double[] work, int _work_offset) {

        if (Lsame.lsame(side, "L")) {
            if (tau != 0.0) {
                Dgemv.dgemv("Transpose", m, n, 1.0, c, _c_offset, ldc, v, _v_offset, incv, 0.0, work, _work_offset, 1);
                Dger.dger(m, n, -tau, v, _v_offset, incv, work, _work_offset, 1, c, _c_offset, ldc);
            }
        } else if (tau != 0.0) {
            Dgemv.dgemv("No transpose", m, n, 1.0, c, _c_offset, ldc, v, _v_offset, incv, 0.0, work, _work_offset, 1);
            Dger.dger(m, n, -tau, work, _work_offset, 1, v, _v_offset, incv, c, _c_offset, ldc);
        }
    }
}
