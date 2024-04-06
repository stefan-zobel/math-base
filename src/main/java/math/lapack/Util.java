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

final class Util {

    static int nint(float f) {
        return (int) (f < 0.0F ? f - 0.5D : f + 0.5D);
    }

    static double dsign(double d, double d1) {
        return d1 < 0.0D ? -Math.abs(d) : Math.abs(d);
    }

    static int min(int i, int j, int k) {
        return Math.min(i >= j ? j : i, Math.min(j, k));
    }

    static int max(int i, int j, int k) {
        return Math.max(i <= j ? j : i, Math.max(j, k));
    }

    static String stringInsert(String s, String s1, int i, int j) {
        return new String(s.substring(0, i - 1) + s1.substring(0, (j - i) + 1) + s.substring(j, s.length()));
    }
}
