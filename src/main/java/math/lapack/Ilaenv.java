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

final class Ilaenv {

    static int ilaenv(int ispec, String name, String opts, int n1, int n2, int n3, int n4) {

        label0: {
            boolean flag1;
            boolean flag2;
            String s1;
            String s2;
            String s3;

            label1: {
                label2: {
                    int i = ispec;
                    int j;

                    if (i != 1 && i != 2 && i != 3) {

                        if (i != 4) {
                            if (i != 5) {
                                if (i != 6) {
                                    if (i != 7) {
                                        if (i != 8) {
                                            if (i != 9) {
                                                if (i != 10) {
                                                    if (i != 11)
                                                        if (i != 12 && i != 13 && i != 14 && i != 15 && i != 16) {
                                                            return -1;
                                                        } else {
                                                            return Iparmq.iparmq(ispec, name, opts, n1, n2, n3, n4);
                                                        }
                                                } else {
                                                    return 1;
                                                }
                                            } else {
                                                return 25;
                                            }
                                        } else {
                                            return 50;
                                        }
                                    } else {
                                        return 1;
                                    }
                                } else {
                                    return (int) ((float) Math.min(n1, n2) * 1.6F);
                                }
                            } else {
                                return 2;
                            }
                        } else {
                            return 6;
                        }
                        break label0;
                    }

                    String s4 = name;
                    char c = s4.substring(0, 1).charAt(0);
                    char c1 = 'Z';
                    if ((c1 == 'Z') || (c1 == 'z')) {
                        if ((c >= 'a') && (c <= 'z')) {
                            s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), 1, 1);
                            j = 2;
                            for (i = 5; i > 0; i--) {
                                c = s4.substring(j + -1, j).charAt(0);
                                if ((c >= 'a') && (c <= 'z'))
                                    s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), j, j);
                                j++;
                            }
                        }
                    } else if ((c1 == '\351') || (c1 == '\251')) {
                        if ((((c >= '\201') && (c <= '\211')) || ((c >= '\221') && (c <= '\231')))
                                || ((c >= '\242') && (c <= '\251'))) {
                            s4 = Util.stringInsert(s4, Character.valueOf((char) (c + 64)).toString(), 1, 1);
                            j = 2;
                            for (i = 5; i > 0; i--) {
                                c = s4.substring(j + -1, j).charAt(0);
                                if ((((c >= '\201') && (c <= '\211')) || ((c >= '\221') && (c <= '\231')))
                                        || ((c >= '\242') && (c <= '\251')))
                                    s4 = Util.stringInsert(s4, Character.valueOf((char) (c + 64)).toString(), j, j);
                                j++;
                            }
                        }
                    } else if (((c1 == '\332') || (c1 == '\372')) && ((c >= '\341') && (c <= '\372'))) {
                        s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), 1, 1);
                        j = 2;
                        for (i = 5; i > 0; i--) {
                            c = s4.substring(j + -1, j).charAt(0);
                            if ((c >= '\341') && (c <= '\372'))
                                s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), j, j);
                            j++;
                        }
                    }
                    String s5 = s4.substring(0, 1);
                    flag2 = s5.regionMatches(0, "S", 0, 1) || s5.regionMatches(0, "D", 0, 1);
                    flag1 = s5.regionMatches(0, "C", 0, 1) || s5.regionMatches(0, "Z", 0, 1);
                    if (!(flag1 || flag2)) {
                        return 1;
                    }
                    s1 = s4.substring(1, 3);
                    s3 = s4.substring(3, 6);
                    s2 = s3.substring(1, 3);
                    i = ispec;
                    if (i != 1) {
                        if (i == 2) {
                            break label2;
                        }
                        if (i == 3) {
                            break label1;
                        }
                    }
                    int x0 = 1;
                    if (s1.regionMatches(0, "GE", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3)) {
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                        } else if (((s3.regionMatches(0, "QRF", 0, 3) || s3.regionMatches(0, "RQF", 0, 3))
                                || s3.regionMatches(0, "LQF", 0, 3)) || s3.regionMatches(0, "QLF", 0, 3)) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (s3.regionMatches(0, "HRD", 0, 3)) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (s3.regionMatches(0, "BRD", 0, 3)) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (s3.regionMatches(0, "TRI", 0, 3))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (s1.regionMatches(0, "PO", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (s1.regionMatches(0, "SY", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3)) {
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                        } else if (flag2 && s3.regionMatches(0, "TRD", 0, 3))
                            x0 = 32;
                        else if (flag2 && s3.regionMatches(0, "GST", 0, 3))
                            x0 = 64;
                    } else if (flag1 && s1.regionMatches(0, "HE", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3))
                            x0 = 64;
                        else if (s3.regionMatches(0, "TRD", 0, 3))
                            x0 = 32;
                        else if (s3.regionMatches(0, "GST", 0, 3))
                            x0 = 64;
                    } else if (flag2 && s1.regionMatches(0, "OR", 0, 2)) {
                        if (s3.substring(0, 1).regionMatches(0, "G", 0, 1)) {
                            if ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                    || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                    || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                    || s2.regionMatches(0, "BR", 0, 2))
                                x0 = 32;
                        } else if (s3.substring(0, 1).regionMatches(0, "M", 0, 1)
                                && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                        || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                        || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                        || s2.regionMatches(0, "BR", 0, 2)))
                            x0 = 32;
                    } else if (flag1 && s1.regionMatches(0, "UN", 0, 2)) {
                        if (s3.substring(0, 1).regionMatches(0, "G", 0, 1)) {
                            if ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                    || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                    || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                    || s2.regionMatches(0, "BR", 0, 2))
                                x0 = 32;
                        } else if (s3.substring(0, 1).regionMatches(0, "M", 0, 1)
                                && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                        || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                        || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                        || s2.regionMatches(0, "BR", 0, 2)))
                            x0 = 32;
                    } else if (s1.regionMatches(0, "GB", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3))
                            if (flag2) {
                                if (n4 <= 64)
                                    x0 = 1;
                                else
                                    x0 = 32;
                            } else if (n4 <= 64)
                                x0 = 1;
                            else
                                x0 = 32;
                    } else if (s1.regionMatches(0, "PB", 0, 2)) {
                        if (s3.regionMatches(0, "TRF", 0, 3))
                            if (flag2) {
                                if (n2 <= 64)
                                    x0 = 1;
                                else
                                    x0 = 32;
                            } else if (n2 <= 64)
                                x0 = 1;
                            else
                                x0 = 32;
                    } else if (s1.regionMatches(0, "TR", 0, 2)) {
                        if (s3.regionMatches(0, "TRI", 0, 3))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (s1.regionMatches(0, "LA", 0, 2)) {
                        if (s3.regionMatches(0, "UUM", 0, 3))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if ((flag2 && s1.regionMatches(0, "ST", 0, 2)) && s3.regionMatches(0, "EBZ", 0, 3))
                        x0 = 1;

                    return x0;
                } // label2
                int x1 = 2;
                if (s1.regionMatches(0, "GE", 0, 2)) {
                    if (((s3.regionMatches(0, "QRF", 0, 3) || s3.regionMatches(0, "RQF", 0, 3))
                            || s3.regionMatches(0, "LQF", 0, 3)) || s3.regionMatches(0, "QLF", 0, 3)) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (s3.regionMatches(0, "HRD", 0, 3)) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (s3.regionMatches(0, "BRD", 0, 3)) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (s3.regionMatches(0, "TRI", 0, 3))
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                } else if (s1.regionMatches(0, "SY", 0, 2)) {
                    if (s3.regionMatches(0, "TRF", 0, 3)) {
                        if (flag2)
                            x1 = 8;
                        else
                            x1 = 8;
                    } else if (flag2 && s3.regionMatches(0, "TRD", 0, 3))
                        x1 = 2;
                } else if (flag1 && s1.regionMatches(0, "HE", 0, 2)) {
                    if (s3.regionMatches(0, "TRD", 0, 3))
                        x1 = 2;
                } else if (flag2 && s1.regionMatches(0, "OR", 0, 2)) {
                    if (s3.substring(0, 1).regionMatches(0, "G", 0, 1)) {
                        if ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                || s2.regionMatches(0, "BR", 0, 2))
                            x1 = 2;
                    } else if (s3.substring(0, 1).regionMatches(0, "M", 0, 1)
                            && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                    || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                    || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                    || s2.regionMatches(0, "BR", 0, 2)))
                        x1 = 2;
                } else if (flag1 && s1.regionMatches(0, "UN", 0, 2))
                    if (s3.substring(0, 1).regionMatches(0, "G", 0, 1)) {
                        if ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                || s2.regionMatches(0, "BR", 0, 2))
                            x1 = 2;
                    } else if (s3.substring(0, 1).regionMatches(0, "M", 0, 1)
                            && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                    || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                    || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                    || s2.regionMatches(0, "BR", 0, 2)))
                        x1 = 2;

                return x1;
            } // label1
            char c2 = '\0';
            if (s1.regionMatches(0, "GE", 0, 2)) {
                if (((s3.regionMatches(0, "QRF", 0, 3) || s3.regionMatches(0, "RQF", 0, 3))
                        || s3.regionMatches(0, "LQF", 0, 3)) || s3.regionMatches(0, "QLF", 0, 3)) {
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
                } else if (s3.regionMatches(0, "HRD", 0, 3)) {
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
                } else if (s3.regionMatches(0, "BRD", 0, 3))
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
            } else if (s1.regionMatches(0, "SY", 0, 2)) {
                if (flag2 && s3.regionMatches(0, "TRD", 0, 3))
                    c2 = ' ';
            } else if (flag1 && s1.regionMatches(0, "HE", 0, 2)) {
                if (s3.regionMatches(0, "TRD", 0, 3))
                    c2 = ' ';
            } else if (flag2 && s1.regionMatches(0, "OR", 0, 2)) {
                if (s3.substring(0, 1).regionMatches(0, "G", 0, 1)
                        && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                                || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                                || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                                || s2.regionMatches(0, "BR", 0, 2)))
                    c2 = '\200';
            } else if ((flag1 && s1.regionMatches(0, "UN", 0, 2)) && s3.substring(0, 1).regionMatches(0, "G", 0, 1)
                    && ((((((s2.regionMatches(0, "QR", 0, 2) || s2.regionMatches(0, "RQ", 0, 2))
                            || s2.regionMatches(0, "LQ", 0, 2)) || s2.regionMatches(0, "QL", 0, 2))
                            || s2.regionMatches(0, "HR", 0, 2)) || s2.regionMatches(0, "TR", 0, 2))
                            || s2.regionMatches(0, "BR", 0, 2)))
                c2 = '\200';
            return c2;
        } // label0
        return 1;
    }
}
