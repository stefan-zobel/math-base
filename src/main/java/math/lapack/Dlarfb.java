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

import math.gemm.Dgemm;
import math.gemm.Trans;
import math.trmm.Dtrmm;
import math.trmm.Side;
import math.trmm.UpLo;
import math.trmm.Diag;

// DLARFB applies a real block reflector H or its transpose H**T
// to a real m by n matrix C, from either the left or the right.
final class Dlarfb {

    static void dlarfb(Side side, Trans trans, String direct, String storev, int m, int n, int k, double[] v,
            int _v_offset, int ldv, double[] t, int _t_offset, int ldt, double[] c, int _c_offset, int ldc,
            double[] work, int _work_offset, int ldwork) {

        if (m <= 0 || n <= 0) {
            return;
        }

        Trans transt;
        if (Trans.NO_TRANS == trans) {
            transt = Trans.TRANS;
        } else {
            transt = Trans.NO_TRANS;
        }

        if (Lsame.lsame(storev, "C")) {
            if (Lsame.lsame(direct, "F")) {
                if (Side.LEFT == side) {
                    int j = 1;
                    for (int p = k; p > 0; p--) {
                        Dcopy.dcopy(n, c, (j - 1) + _c_offset, ldc, work, (j - 1) * ldwork + _work_offset, 1);
                        j++;
                    }

                    Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.NO_TRANS, Diag.UNIT, n, k, 1.0, v, _v_offset, ldv, work,
                            _work_offset, ldwork);
                    if (m > k) {
                        Dgemm.dgemm(Trans.TRANS, Trans.NO_TRANS, n, k, m - k, 1.0, c, k + _c_offset, ldc, v,
                                k + _v_offset, ldv, 1.0, work, _work_offset, ldwork);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, transt, Diag.NON_UNIT, n, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                            ldwork);
                    if (m > k) {
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m - k, n, k, -1.0, v, k + _v_offset, ldv, work,
                                _work_offset, ldwork, 1.0, c, k + _c_offset, ldc);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.TRANS, Diag.UNIT, n, k, 1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork);
                    j = 1;
                    for (int p = k; p > 0; p--) {
                        int i = 1;
                        for (int q = n; q > 0; q--) {
                            c[(j - 1) + (i - 1) * ldc + _c_offset] = c[(j - 1) + (i - 1) * ldc + _c_offset]
                                    - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                            i++;
                        }

                        j++;
                    }

                } else if (Side.RIGHT == side) {
                    int j = 1;
                    for (int p = k; p > 0; p--) {
                        Dcopy.dcopy(m, c, (j - 1) * ldc + _c_offset, 1, work, (j - 1) * ldwork + _work_offset, 1);
                        j++;
                    }

                    Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.NO_TRANS, Diag.UNIT, m, k, 1.0, v, _v_offset, ldv, work,
                            _work_offset, ldwork);
                    if (n > k) {
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, m, k, n - k, 1.0, c, k * ldc + _c_offset, ldc, v,
                                k * ldv + _v_offset, ldv, 1.0, work, _work_offset, ldwork);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, trans, Diag.NON_UNIT, m, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                            ldwork);
                    if (n > k) {
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m, n - k, k, -1.0, work, _work_offset, ldwork, v,
                                k + _v_offset, ldv, 1.0, c, k * ldc + _c_offset, ldc);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.TRANS, Diag.UNIT, m, k, 1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork);
                    j = 1;
                    for (int p = k; p > 0; p--) {
                        int i = 1;
                        for (int q = m; q > 0; q--) {
                            c[(i - 1) + (j - 1) * ldc + _c_offset] = c[(i - 1) + (j - 1) * ldc + _c_offset]
                                    - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                            i++;
                        }

                        j++;
                    }

                }
            } else if (Side.LEFT == side) {
                int j = 1;
                for (int p = k; p > 0; p--) {
                    Dcopy.dcopy(n, c, ((m - k) + j - 1) + _c_offset, ldc, work, (j - 1) * ldwork + _work_offset, 1);
                    j++;
                }

                Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.NO_TRANS, Diag.UNIT, n, k, 1.0, v, (m - k) + _v_offset, ldv, work,
                        _work_offset, ldwork);
                if (m > k) {
                    Dgemm.dgemm(Trans.TRANS, Trans.NO_TRANS, n, k, m - k, 1.0, c, _c_offset, ldc, v, _v_offset, ldv,
                            1.0, work, _work_offset, ldwork);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, transt, Diag.NON_UNIT, n, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                        ldwork);
                if (m > k) {
                    Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m - k, n, k, -1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork, 1.0, c, _c_offset, ldc);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.TRANS, Diag.UNIT, n, k, 1.0, v, (m - k) + _v_offset, ldv, work,
                        _work_offset, ldwork);
                j = 1;
                for (int p = k; p > 0; p--) {
                    int i = 1;
                    for (int q = n; q > 0; q--) {
                        c[(((m - k) + j) - 1) + (i - 1) * ldc
                                + _c_offset] = c[(((m - k) + j) - 1) + (i - 1) * ldc + _c_offset]
                                        - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                        i++;
                    }

                    j++;
                }

            } else if (Side.RIGHT == side) {
                int j = 1;
                for (int p = k; p > 0; p--) {
                    Dcopy.dcopy(m, c, (((n - k) + j) - 1) * ldc + _c_offset, 1, work, (j - 1) * ldwork + _work_offset,
                            1);
                    j++;
                }

                Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.NO_TRANS, Diag.UNIT, m, k, 1.0, v, (n - k) + _v_offset, ldv, work,
                        _work_offset, ldwork);
                if (n > k) {
                    Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, m, k, n - k, 1.0, c, _c_offset, ldc, v, _v_offset, ldv,
                            1.0, work, _work_offset, ldwork);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, trans, Diag.NON_UNIT, m, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                        ldwork);
                if (n > k) {
                    Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m, n - k, k, -1.0, work, _work_offset, ldwork, v,
                            _v_offset, ldv, 1.0, c, _c_offset, ldc);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.TRANS, Diag.UNIT, m, k, 1.0, v, (n - k) + _v_offset, ldv, work,
                        _work_offset, ldwork);
                j = 1;
                for (int p = k; p > 0; p--) {
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        c[(i - 1) + (((n - k) + j) - 1) * ldc
                                + _c_offset] = c[(i - 1) + (((n - k) + j) - 1) * ldc + _c_offset]
                                        - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                        i++;
                    }

                    j++;
                }

            }
        } else if (Lsame.lsame(storev, "R"))
            if (Lsame.lsame(direct, "F")) {
                if (Side.LEFT == side) {
                    int j = 1;
                    for (int p = k; p > 0; p--) {
                        Dcopy.dcopy(n, c, (j - 1) + _c_offset, ldc, work, (j - 1) * ldwork + _work_offset, 1);
                        j++;
                    }

                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.TRANS, Diag.UNIT, n, k, 1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork);
                    if (m > k) {
                        Dgemm.dgemm(Trans.TRANS, Trans.TRANS, n, k, m - k, 1.0, c, k + _c_offset, ldc, v,
                                k * ldv + _v_offset, ldv, 1.0, work, _work_offset, ldwork);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, transt, Diag.NON_UNIT, n, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                            ldwork);
                    if (m > k) {
                        Dgemm.dgemm(Trans.TRANS, Trans.TRANS, m - k, n, k, -1.0, v, k * ldv + _v_offset, ldv, work,
                                _work_offset, ldwork, 1.0, c, k + _c_offset, ldc);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.NO_TRANS, Diag.UNIT, n, k, 1.0, v, _v_offset, ldv, work,
                            _work_offset, ldwork);
                    j = 1;
                    for (int p = k; p > 0; p--) {
                        int i = 1;
                        for (int q = n; q > 0; q--) {
                            c[(j - 1) + (i - 1) * ldc + _c_offset] = c[(j - 1) + (i - 1) * ldc + _c_offset]
                                    - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                            i++;
                        }

                        j++;
                    }

                } else if (Side.RIGHT == side) {
                    int j = 1;
                    for (int p = k; p > 0; p--) {
                        Dcopy.dcopy(m, c, (j - 1) * ldc + _c_offset, 1, work, (j - 1) * ldwork + _work_offset, 1);
                        j++;
                    }

                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.TRANS, Diag.UNIT, m, k, 1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork);
                    if (n > k) {
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m, k, n - k, 1.0, c, k * ldc + _c_offset, ldc, v,
                                k * ldv + _v_offset, ldv, 1.0, work, _work_offset, ldwork);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, trans, Diag.NON_UNIT, m, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                            ldwork);
                    if (n > k) {
                        Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, m, n - k, k, -1.0, work, _work_offset, ldwork, v,
                                k * ldv + _v_offset, ldv, 1.0, c, k * ldc + _c_offset, ldc);
                    }
                    Dtrmm.dtrmm(Side.RIGHT, UpLo.UPPER, Trans.NO_TRANS, Diag.UNIT, m, k, 1.0, v, _v_offset, ldv, work,
                            _work_offset, ldwork);
                    j = 1;
                    for (int p = k; p > 0; p--) {
                        int i = 1;
                        for (int q = m; q > 0; q--) {
                            c[(i - 1) + (j - 1) * ldc + _c_offset] = c[(i - 1) + (j - 1) * ldc + _c_offset]
                                    - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                            i++;
                        }

                        j++;
                    }

                }
            } else if (Side.LEFT == side) {
                int j = 1;
                for (int p = k; p > 0; p--) {
                    Dcopy.dcopy(n, c, (m - k) + j - 1 + _c_offset, ldc, work, (j - 1) * ldwork + _work_offset, 1);
                    j++;
                }

                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.TRANS, Diag.UNIT, n, k, 1.0, v, (m - k) * ldv + _v_offset, ldv, work,
                        _work_offset, ldwork);
                if (m > k) {
                    Dgemm.dgemm(Trans.TRANS, Trans.TRANS, n, k, m - k, 1.0, c, _c_offset, ldc, v, _v_offset, ldv, 1.0,
                            work, _work_offset, ldwork);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, transt, Diag.NON_UNIT, n, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                        ldwork);
                if (m > k) {
                    Dgemm.dgemm(Trans.TRANS, Trans.TRANS, m - k, n, k, -1.0, v, _v_offset, ldv, work, _work_offset,
                            ldwork, 1.0, c, _c_offset, ldc);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.NO_TRANS, Diag.UNIT, n, k, 1.0, v, (m - k) * ldv + _v_offset, ldv,
                        work, _work_offset, ldwork);
                j = 1;
                for (int p = k; p > 0; p--) {
                    int i = 1;
                    for (int q = n; q > 0; q--) {
                        c[((m - k) + j - 1) + (i - 1) * ldc
                                + _c_offset] = c[((m - k) + j - 1) + (i - 1) * ldc + _c_offset]
                                        - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                        i++;
                    }

                    j++;
                }

            } else if (Side.RIGHT == side) {
                int j = 1;
                for (int p = k; p > 0; p--) {
                    Dcopy.dcopy(m, c, ((n - k) + j - 1) * ldc + _c_offset, 1, work, (j - 1) * ldwork + _work_offset, 1);
                    j++;
                }

                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.TRANS, Diag.UNIT, m, k, 1.0, v, (n - k) * ldv + _v_offset, ldv, work,
                        _work_offset, ldwork);
                if (n > k) {
                    Dgemm.dgemm(Trans.NO_TRANS, Trans.TRANS, m, k, n - k, 1.0, c, _c_offset, ldc, v, _v_offset, ldv,
                            1.0, work, _work_offset, ldwork);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, trans, Diag.NON_UNIT, m, k, 1.0, t, _t_offset, ldt, work, _work_offset,
                        ldwork);
                if (n > k) {
                    Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, m, n - k, k, -1.0, work, _work_offset, ldwork, v,
                            _v_offset, ldv, 1.0, c, _c_offset, ldc);
                }
                Dtrmm.dtrmm(Side.RIGHT, UpLo.LOWER, Trans.NO_TRANS, Diag.UNIT, m, k, 1.0, v, (n - k) * ldv + _v_offset, ldv,
                        work, _work_offset, ldwork);
                j = 1;
                for (int p = k; p > 0; p--) {
                    int i = 1;
                    for (int q = m; q > 0; q--) {
                        c[(i - 1) + ((n - k) + j - 1) * ldc
                                + _c_offset] = c[(i - 1) + ((n - k) + j - 1) * ldc + _c_offset]
                                        - work[(i - 1) + (j - 1) * ldwork + _work_offset];
                        i++;
                    }

                    j++;
                }

            }
    }
}
