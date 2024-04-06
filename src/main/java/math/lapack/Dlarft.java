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

//DLARFT forms the triangular factor T of a real block reflector H
//of order n, which is defined as a product of k elementary reflectors.
//
//If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
//
//If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
//
//If STOREV = 'C', the vector which defines the elementary reflector
//H(i) is stored in the i-th column of the array V, and
//
//   H  =  I - V * T * V**T
//
//If STOREV = 'R', the vector which defines the elementary reflector
//H(i) is stored in the i-th row of the array V, and
//
//   H  =  I - V**T * T * V
final class Dlarft {

    static void dlarft(String direct, String storev, int n, int k, double[] v, int _v_offset, int ldv,
            double[] tau, int _tau_offset, double[] t, int _t_offset, int ldt) {

        if (n == 0) {
            return;
        }
        if (Lsame.lsame(direct, "F")) {
            int I = 1;
            for (int p = k; p > 0; p--) {
                if (tau[(I - 1) + _tau_offset] == 0.0) {
                    int J = 1;
                    for (int q = I; q > 0; q--) {
                        t[(J - 1) + (I - 1) * ldt + _t_offset] = 0.0;
                        J++;
                    }

                } else {
                    double VII = v[(I - 1) + (I - 1) * ldv + _v_offset];
                    v[(I - 1) + (I - 1) * ldv + _v_offset] = 1.0;
                    if (Lsame.lsame(storev, "C")) {
                        Dgemv.dgemv("Transpose", (n - I) + 1, I - 1, -tau[(I - 1) + _tau_offset], v,
                                (I - 1) + _v_offset, ldv, v, (I - 1) + (I - 1) * ldv + _v_offset, 1, 0.0, t,
                                (I - 1) * ldt + _t_offset, 1);
                    } else {
                        Dgemv.dgemv("No transpose", I - 1, (n - I) + 1, -tau[(I - 1) + _tau_offset], v,
                                (I - 1) * ldv + _v_offset, ldv, v, (I - 1) + (I - 1) * ldv + _v_offset, ldv, 0.0, t,
                                (I - 1) * ldt + _t_offset, 1);
                    }
                    v[(I - 1) + (I - 1) * ldv + _v_offset] = VII;
                    Dtrmv.dtrmv("Upper", "No transpose", "Non-unit", I - 1, t, _t_offset, ldt, t,
                            (I - 1) * ldt + _t_offset, 1);
                    t[(I - 1) + (I - 1) * ldt + _t_offset] = tau[(I - 1) + _tau_offset];
                }
                I++;
            }

        } else {
            int I = k;
            for (int p = k; p > 0; p--) {
                if (tau[(I - 1) + _tau_offset] == 0.0) {
                    int J = I;
                    for (int q = (k - I) + 1; q > 0; q--) {
                        t[(J - 1) + (I - 1) * ldt + _t_offset] = 0.0;
                        J++;
                    }

                } else {
                    if (I < k) {
                        if (Lsame.lsame(storev, "C")) {
                            double VII = v[(((n - k) + I) - 1) + (I - 1) * ldv + _v_offset];
                            v[(((n - k) + I) - 1) + (I - 1) * ldv + _v_offset] = 1.0;
                            Dgemv.dgemv("Transpose", (n - k) + I, k - I, -tau[(I - 1) + _tau_offset], v,
                                    I * ldv + _v_offset, ldv, v, (I - 1) * ldv + _v_offset, 1, 0.0, t,
                                    I + (I - 1) * ldt + _t_offset, 1);
                            v[(((n - k) + I) - 1) + (I - 1) * ldv + _v_offset] = VII;
                        } else {
                            double VII = v[(I - 1) + (((n - k) + I) - 1) * ldv + _v_offset];
                            v[(I - 1) + (((n - k) + I) - 1) * ldv + _v_offset] = 1.0;
                            Dgemv.dgemv("No transpose", k - I, (n - k) + I, -tau[(I - 1) + _tau_offset], v,
                                    I + _v_offset, ldv, v, (I - 1) + _v_offset, ldv, 0.0, t,
                                    I + (I - 1) * ldt + _t_offset, 1);
                            v[(I - 1) + (((n - k) + I) - 1) * ldv + _v_offset] = VII;
                        }
                        Dtrmv.dtrmv("Lower", "No transpose", "Non-unit", k - I, t, I + I * ldt + _t_offset, ldt, t,
                                I + (I - 1) * ldt + _t_offset, 1);
                    }
                    t[(I - 1) + (I - 1) * ldt + _t_offset] = tau[(I - 1) + _tau_offset];
                }
                I--;
            }
        }
    }
}
