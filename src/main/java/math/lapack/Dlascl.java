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

// DLASCL multiplies the M by N real matrix A by the real scalar
// CTO/CFROM.  This is done without over/underflow as long as the final
// result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
// A may be full, upper triangular, lower triangular, upper Hessenberg,
// or banded.
final class Dlascl {

    static void dlascl(String type, int kl, int ku, double cfrom, double cto, int m, int n, double[] a,
            int _a_offset, int lda, intW info) {

        info.val = 0;
        int itype = 0;
        if (Lsame.lsame(type, "G")) {
            itype = 0;
        } else if (Lsame.lsame(type, "L")) {
            itype = 1;
        } else if (Lsame.lsame(type, "U")) {
            itype = 2;
        } else if (Lsame.lsame(type, "H")) {
            itype = 3;
        } else if (Lsame.lsame(type, "B")) {
            itype = 4;
        } else if (Lsame.lsame(type, "Q")) {
            itype = 5;
        } else if (Lsame.lsame(type, "Z")) {
            itype = 6;
        } else {
            itype = -1;
        }

        if (itype == -1) {
            info.val = -1;
        } else if (cfrom == 0.0) {
            info.val = -4;
        } else if (m < 0) {
            info.val = -6;
        } else if ((n < 0 || (itype == 4 && n != m)) || (itype == 5 && n != m)) {
            info.val = -7;
        } else if (itype <= 3 && lda < Math.max(1, m)) {
            info.val = -9;
        } else if (itype >= 4) {
            if (kl < 0 || kl > Math.max(m - 1, 0)) {
                info.val = -2;
            } else if ((ku < 0 || ku > Math.max(n - 1, 0)) || ((itype == 4 || itype == 5) && (kl != ku))) {
                info.val = -3;
            } else if (((itype == 4 && lda < kl + 1) || (itype == 5 && lda < ku + 1))
                    || (itype == 6 && lda < 2 * kl + ku + 1)) {
                info.val = -9;
            }
        }

        if (info.val != 0) {
            Xerbla.xerbla("DLASCL", -info.val);
            return;
        }
        // Quick return if possible
        if (n == 0 || m == 0) {
            return;
        }

        boolean done = false;
        // Get machine parameters
        final double smlnum = 2.2250738585072014E-308;
        final double bignum = 4.49423283715579E307;
        double cfromc = cfrom;
        double ctoc = cto;
        do {
            double cfrom1 = cfromc * smlnum;
            double cto1 = ctoc / bignum;
            double mul;
            if (Math.abs(cfrom1) > Math.abs(ctoc) && ctoc != 0.0) {
                mul = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if (Math.abs(cto1) > Math.abs(cfromc)) {
                mul = bignum;
                done = false;
                ctoc = cto1;
            } else {
                mul = ctoc / cfromc;
                done = true;
            }
            if (itype == 0) {
                // Full matrix
                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 1) {
                // Lower triangular matrix
                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = j;
                    for (int q = m - j + 1; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 2) {
                // Upper triangular matrix
                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = 1;
                    for (int q = Math.min(j, m); q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 3) {
                // Upper Hessenberg matrix
                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = 1;
                    for (int q = Math.min(j + 1, m); q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 4) {
                // Lower half of a symmetric band matrix
                int k3 = kl + 1;
                int k4 = n + 1;

                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = 1;
                    for (int q = Math.min(k3, k4 - j); q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 5) {
                // Upper half of a symmetric band matrix
                int k1 = ku + 2;
                int k3 = ku + 1;

                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = Math.max(k1 - j, 1);
                    for (int q = k3 - Math.max(k1 - j, 1) + 1; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            } else if (itype == 6) {
                // Band matrix
                int k1 = kl + ku + 2;
                int k2 = kl + 1;
                int k3 = 2 * kl + ku + 1;
                int k4 = kl + ku + 1 + m;

                int j = 1;
                for (int p = n; p > 0; p--) {
                    int i = Math.max(k1 - j, k2);
                    for (int q = (Math.min(k3, k4 - j) - Math.max(k1 - j, k2)) + 1; q > 0; q--) {
                        a[i - 1 + (j - 1) * lda + _a_offset] = a[i - 1 + (j - 1) * lda + _a_offset] * mul;
                        i++;
                    }

                    j++;
                }

            }
        } while (!done);
    }
}
