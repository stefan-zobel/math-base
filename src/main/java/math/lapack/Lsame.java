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

// LSAME returns true if its two arguments begin with the same letter,
// regardless of case.
final class Lsame {

    // The f2j translation wrote this as
    //     s1.substring(0, 1).equalsIgnoreCase(s.substring(0, 1))
    // which allocates two strings on every call. It is called at the entry of
    // every routine in this package and inside the loops of Dlarft and
    // Dlarfb, and removing the allocation measured 1.7x to 2.8x on
    // LuFactorization.solve at orders 3 to 12, which is where Rosenbrock
    // works.
    static boolean lsame(String s, String s1) {
        char a = s.charAt(0);
        char b = s1.charAt(0);
        if (a == b) {
            return true;
        }
        return Character.toUpperCase(a) == Character.toUpperCase(b);
    }
}
