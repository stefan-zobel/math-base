/*
 * Copyright 2018 Stefan Zobel
 *
 * The original C code this was derived from is Copyright Dr. Michael Lehn,
 * Ulm University
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

final class Mgemm {

    static void mgemm(int mc, int nc, int kc, double alpha, double[] A_, double[] B_, double beta, int C_start,
            double[] C, int incRowC, int incColC, double[] AB, double[] workC) {

        final int mp = (mc + BlockSizes.MR - 1) / BlockSizes.MR;
        final int np = (nc + BlockSizes.NR - 1) / BlockSizes.NR;

        final int mr_ = (mc & (BlockSizes.MR - 1));
        final int nr_ = (nc & (BlockSizes.NR - 1));

        int mr, nr;

        for (int j = 0; j < np; ++j) {
            nr = (j != np - 1 || nr_ == 0) ? BlockSizes.NR : nr_;

            for (int i = 0; i < mp; ++i) {
                mr = (i != mp - 1 || mr_ == 0) ? BlockSizes.MR : mr_;

                if (mr == BlockSizes.MR && nr == BlockSizes.NR) {
                    Ugemm.ugemm(kc, alpha, (i * kc * BlockSizes.MR), // A_ start
                            A_, (j * kc * BlockSizes.NR), // B_ start
                            B_, beta, (C_start + i * BlockSizes.MR * incRowC + j * BlockSizes.NR * incColC), // C
                                                                                                                // start
                            C, incRowC, incColC, AB);
                } else {
                    // Call the buffered micro kernel
                    Ugemm.ugemm(mr, nr, kc, alpha, (i * kc * BlockSizes.MR), // A_
                                                                                // start
                            A_, (j * kc * BlockSizes.NR), // B_ start
                            B_, beta, (C_start + i * BlockSizes.MR * incRowC + j * BlockSizes.NR * incColC), // C
                                                                                                                // start
                            C, incRowC, incColC, AB, workC);
                }
            }
        }
    }

    private Mgemm() {
        throw new AssertionError();
    }
}
