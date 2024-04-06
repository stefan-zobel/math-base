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

// This program sets problem and machine dependent parameters
// useful for xHSEQR and related subroutines for eigenvalue
// problems. It is called whenever
// IPARMQ is called with 12 <= ISPEC <= 16
final class Iparmq {

    static int iparmq(int ispec, String name, String opts, int n, int ilo, int ihi, int lwork) {

        int nh = 0;
        int ns = 0;
        int iparmq = 0;
        if (ispec == 15 || ispec == 13 || ispec == 16) {
            nh = (ihi - ilo) + 1;
            ns = 2;
            if (nh >= 30) {
                ns = 4;
            }
            if (nh >= 60) {
                ns = 10;
            }
            if (nh >= 150) {
                ns = Math.max(10, nh / Util.nint((float) Math.log(nh) / (float) Math.log(2.0f)));
            }
            if (nh >= 590) {
                ns = 64;
            }
            if (nh >= 3000) {
                ns = 128;
            }
            if (nh >= 6000) {
                ns = 256;
            }
            ns = Math.max(2, ns - ns % 2);
        }
        if (ispec == 12) {
            iparmq = 75;
        } else if (ispec == 14) {
            iparmq = 14;
        } else if (ispec == 15) {
            iparmq = ns;
        } else if (ispec == 13) {
            if (nh <= 500) {
                iparmq = ns;
            } else {
                iparmq = (3 * ns) / 2;
            }
        } else if (ispec == 16) {
            iparmq = 0;
            if (ns >= 14) {
                iparmq = 1;
            }
            if (ns >= 14) {
                iparmq = 2;
            }
        } else {
            iparmq = -1;
        }
        return iparmq;
    }
}
