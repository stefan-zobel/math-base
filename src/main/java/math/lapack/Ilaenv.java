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

        if (name == null || name.length() < 6) {
            // the name is decomposed as name[1..2] and name[3..5] below,
            // so anything shorter would throw. Reference LAPACK pads
            // SUBNAM to CHARACTER*6 and then fails to recognize it; 1 is
            // what this routine already returns for a name it does not
            // know, see the flag1/flag2 test further down
            return 1;
        }

        label0: {
            boolean flag1;
            boolean flag2;
            // the uppercased routine name. It is read at fixed offsets --
            // LAPACK's SUBNAM is a one-character type prefix, a two-character
            // matrix group and a three-character routine -- so it has to be
            // visible in the ispec 2 and ispec 3 blocks below, which sit
            // outside the labeled block it is assigned in
            String s4;

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

                    s4 = name;
                    char c = name.charAt(0);
                    // Reference ILAENV asks the runtime what number 'Z' has and
                    // uppercases the name with the arithmetic that fits the
                    // answer: 90 or 122 for ASCII, 233 or 169 for EBCDIC, whose
                    // letters are not contiguous, and 218 or 250 for Prime
                    // machines, whose PRIMOS stored characters as ASCII with the
                    // high bit set. Fortran 77 had neither a toupper nor a
                    // guaranteed character set.
                    //
                    // A Java char is UTF-16, so 'Z' is 90 on every JVM and the
                    // other two branches could not be reached under any
                    // configuration. They are gone; the ASCII arithmetic below is
                    // what the reference does for IZ == 90.
                    if ((c >= 'a') && (c <= 'z')) {
                        s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), 1, 1);
                        j = 2;
                        for (i = 5; i > 0; i--) {
                            c = s4.charAt(j - 1);
                            if ((c >= 'a') && (c <= 'z')) {
                                s4 = Util.stringInsert(s4, Character.valueOf((char) (c - 32)).toString(), j, j);
                            }
                            j++;
                        }
                    }
                    flag2 = prefix(s4, "S") || prefix(s4, "D");
                    flag1 = prefix(s4, "C") || prefix(s4, "Z");
                    if (!(flag1 || flag2)) {
                        return 1;
                    }
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
                    if (group(s4, "GE")) {
                        if (routine(s4, "TRF")) {
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                        } else if (((routine(s4, "QRF") || routine(s4, "RQF"))
                                || routine(s4, "LQF")) || routine(s4, "QLF")) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (routine(s4, "HRD")) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (routine(s4, "BRD")) {
                            if (flag2)
                                x0 = 32;
                            else
                                x0 = 32;
                        } else if (routine(s4, "TRI"))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (group(s4, "PO")) {
                        if (routine(s4, "TRF"))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (group(s4, "SY")) {
                        if (routine(s4, "TRF")) {
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                        } else if (flag2 && routine(s4, "TRD"))
                            x0 = 32;
                        else if (flag2 && routine(s4, "GST"))
                            x0 = 64;
                    } else if (flag1 && group(s4, "HE")) {
                        if (routine(s4, "TRF"))
                            x0 = 64;
                        else if (routine(s4, "TRD"))
                            x0 = 32;
                        else if (routine(s4, "GST"))
                            x0 = 64;
                    } else if (flag2 && group(s4, "OR")) {
                        if (action(s4, "G")) {
                            if ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                    || pair(s4, "LQ")) || pair(s4, "QL"))
                                    || pair(s4, "HR")) || pair(s4, "TR"))
                                    || pair(s4, "BR"))
                                x0 = 32;
                        } else if (action(s4, "M")
                                && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                        || pair(s4, "LQ")) || pair(s4, "QL"))
                                        || pair(s4, "HR")) || pair(s4, "TR"))
                                        || pair(s4, "BR")))
                            x0 = 32;
                    } else if (flag1 && group(s4, "UN")) {
                        if (action(s4, "G")) {
                            if ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                    || pair(s4, "LQ")) || pair(s4, "QL"))
                                    || pair(s4, "HR")) || pair(s4, "TR"))
                                    || pair(s4, "BR"))
                                x0 = 32;
                        } else if (action(s4, "M")
                                && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                        || pair(s4, "LQ")) || pair(s4, "QL"))
                                        || pair(s4, "HR")) || pair(s4, "TR"))
                                        || pair(s4, "BR")))
                            x0 = 32;
                    } else if (group(s4, "GB")) {
                        if (routine(s4, "TRF"))
                            if (flag2) {
                                if (n4 <= 64)
                                    x0 = 1;
                                else
                                    x0 = 32;
                            } else if (n4 <= 64)
                                x0 = 1;
                            else
                                x0 = 32;
                    } else if (group(s4, "PB")) {
                        if (routine(s4, "TRF"))
                            if (flag2) {
                                if (n2 <= 64)
                                    x0 = 1;
                                else
                                    x0 = 32;
                            } else if (n2 <= 64)
                                x0 = 1;
                            else
                                x0 = 32;
                    } else if (group(s4, "TR")) {
                        if (routine(s4, "TRI"))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if (group(s4, "LA")) {
                        if (routine(s4, "UUM"))
                            if (flag2)
                                x0 = 64;
                            else
                                x0 = 64;
                    } else if ((flag2 && group(s4, "ST")) && routine(s4, "EBZ"))
                        x0 = 1;

                    return x0;
                } // label2
                int x1 = 2;
                if (group(s4, "GE")) {
                    if (((routine(s4, "QRF") || routine(s4, "RQF"))
                            || routine(s4, "LQF")) || routine(s4, "QLF")) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (routine(s4, "HRD")) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (routine(s4, "BRD")) {
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                    } else if (routine(s4, "TRI"))
                        if (flag2)
                            x1 = 2;
                        else
                            x1 = 2;
                } else if (group(s4, "SY")) {
                    if (routine(s4, "TRF")) {
                        if (flag2)
                            x1 = 8;
                        else
                            x1 = 8;
                    } else if (flag2 && routine(s4, "TRD"))
                        x1 = 2;
                } else if (flag1 && group(s4, "HE")) {
                    if (routine(s4, "TRD"))
                        x1 = 2;
                } else if (flag2 && group(s4, "OR")) {
                    if (action(s4, "G")) {
                        if ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                || pair(s4, "LQ")) || pair(s4, "QL"))
                                || pair(s4, "HR")) || pair(s4, "TR"))
                                || pair(s4, "BR"))
                            x1 = 2;
                    } else if (action(s4, "M")
                            && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                    || pair(s4, "LQ")) || pair(s4, "QL"))
                                    || pair(s4, "HR")) || pair(s4, "TR"))
                                    || pair(s4, "BR")))
                        x1 = 2;
                } else if (flag1 && group(s4, "UN"))
                    if (action(s4, "G")) {
                        if ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                || pair(s4, "LQ")) || pair(s4, "QL"))
                                || pair(s4, "HR")) || pair(s4, "TR"))
                                || pair(s4, "BR"))
                            x1 = 2;
                    } else if (action(s4, "M")
                            && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                    || pair(s4, "LQ")) || pair(s4, "QL"))
                                    || pair(s4, "HR")) || pair(s4, "TR"))
                                    || pair(s4, "BR")))
                        x1 = 2;

                return x1;
            } // label1
            char c2 = '\0';
            if (group(s4, "GE")) {
                if (((routine(s4, "QRF") || routine(s4, "RQF"))
                        || routine(s4, "LQF")) || routine(s4, "QLF")) {
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
                } else if (routine(s4, "HRD")) {
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
                } else if (routine(s4, "BRD"))
                    if (flag2)
                        c2 = '\200';
                    else
                        c2 = '\200';
            } else if (group(s4, "SY")) {
                if (flag2 && routine(s4, "TRD"))
                    c2 = ' ';
            } else if (flag1 && group(s4, "HE")) {
                if (routine(s4, "TRD"))
                    c2 = ' ';
            } else if (flag2 && group(s4, "OR")) {
                if (action(s4, "G")
                        && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                                || pair(s4, "LQ")) || pair(s4, "QL"))
                                || pair(s4, "HR")) || pair(s4, "TR"))
                                || pair(s4, "BR")))
                    c2 = '\200';
            } else if ((flag1 && group(s4, "UN")) && action(s4, "G")
                    && ((((((pair(s4, "QR") || pair(s4, "RQ"))
                            || pair(s4, "LQ")) || pair(s4, "QL"))
                            || pair(s4, "HR")) || pair(s4, "TR"))
                            || pair(s4, "BR")))
                c2 = '\200';
            return c2;
        } // label0
        return 1;
    }

    // LAPACK addresses a routine by the layout of its name, SUBNAM, which is a
    // one-character type prefix, a two-character matrix group and a
    // three-character routine whose first letter is the action. The f2j
    // translation cut the name into four substrings to match those pieces,
    // which allocated five to seven strings on every call; regionMatches takes
    // an offset into the receiver, so the pieces never have to exist.

    /** the type prefix: S or D for real, C or Z for complex */
    private static boolean prefix(String s, String lit) {
        return s.regionMatches(0, lit, 0, 1);
    }

    /** the two-character matrix group, GE, SY, PO and the rest */
    private static boolean group(String s, String lit) {
        return s.regionMatches(1, lit, 0, 2);
    }

    /** the three-character routine, TRF, QRF, TRD and the rest */
    private static boolean routine(String s, String lit) {
        return s.regionMatches(3, lit, 0, 3);
    }

    /** the first letter of the routine alone: G generates a factor, M multiplies by one */
    private static boolean action(String s, String lit) {
        return s.regionMatches(3, lit, 0, 1);
    }

    /** the two characters after the action, naming the factorization: QR, LQ and the rest */
    private static boolean pair(String s, String lit) {
        return s.regionMatches(4, lit, 0, 2);
    }
}
