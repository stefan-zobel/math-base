/*
 * Copyright 2018 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.trmm;

/**
 * Trampoline into Michael Lehn's cache-friendly trlmm / trumm routines.
 */
final class Trmm {

    static void trmm(boolean leftSide, boolean lowerTriang, boolean transA, boolean unitDiag, int m, int n,
            double alpha, double[] A, int _a_offset, int lda, double[] B, int _b_offset, int ldb) {

        if (leftSide) {
            if (lowerTriang) {
                if (!transA) {
                    TrlmmLehn.trlmm(m, n, alpha, unitDiag, _a_offset, A, 1, lda, _b_offset, B, 1, ldb);
                } else {
                    TrummLehn.trumm(m, n, alpha, unitDiag, _a_offset, A, lda, 1, _b_offset, B, 1, ldb);
                }
            } else {
                if (!transA) {
                    TrummLehn.trumm(m, n, alpha, unitDiag, _a_offset, A, 1, lda, _b_offset, B, 1, ldb);
                } else {
                    TrlmmLehn.trlmm(m, n, alpha, unitDiag, _a_offset, A, lda, 1, _b_offset, B, 1, ldb);
                }
            }
        } else {
            if (lowerTriang) {
                if (!transA) {
                    TrummLehn.trumm(n, m, alpha, unitDiag, _a_offset, A, lda, 1, _b_offset, B, ldb, 1);
                } else {
                    TrlmmLehn.trlmm(n, m, alpha, unitDiag, _a_offset, A, 1, lda, _b_offset, B, ldb, 1);
                }
            } else {
                if (!transA) {
                    TrlmmLehn.trlmm(n, m, alpha, unitDiag, _a_offset, A, lda, 1, _b_offset, B, ldb, 1);
                } else {
                    TrummLehn.trumm(n, m, alpha, unitDiag, _a_offset, A, 1, lda, _b_offset, B, ldb, 1);
                }
            }
        }
    }

    private Trmm() {
        throw new AssertionError();
    }
}
