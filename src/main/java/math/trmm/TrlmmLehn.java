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

final class TrlmmLehn {

    static void trlmm(int m, int n, double alpha, boolean unitDiag, int A_start, double[] A, int incRowA,
            int incColA, int B_start, double[] B, int incRowB, int incColB) {

        if (alpha == 0.0) {
            Gescal.gescal(m, n, 0.0, B_start, B, incRowB, incColB);
            return;
        }

        final int mb = (m + BlockSizes.MC - 1) / BlockSizes.MC;
        final int nb = (n + BlockSizes.NC - 1) / BlockSizes.NC;

        final int mc_ = m % BlockSizes.MC;
        final int nc_ = (n & (BlockSizes.NC - 1));

        final double[] A_ = new double[BlockSizes.MC * BlockSizes.MC + BlockSizes.MR];
        final double[] B_ = new double[BlockSizes.MC * BlockSizes.NC + BlockSizes.NR];
        final double[] C_ = new double[BlockSizes.MR * BlockSizes.NR];
        final double[] AB = new double[BlockSizes.MR * BlockSizes.NR];

        for (int j = 0; j < nb; ++j) {
            int nc = (j != nb - 1 || nc_ == 0) ? BlockSizes.NC : nc_;

            for (int l = mb - 1; l >= 0; --l) {
                int kc = (l != mb - 1 || mc_ == 0) ? BlockSizes.MC : mc_;

                Gepack.gepack_B(kc, nc, (B_start + l * BlockSizes.MC * incRowB + j * BlockSizes.NC * incColB), // B
                                                                                                                // start
                        B, incRowB, incColB, B_);

                Trlpack.trlpack(kc, unitDiag, (A_start + l * BlockSizes.MC * (incRowA + incColA)), // A
                                                                                                    // start
                        A, incRowA, incColA, A_);

                Mtrlmm.mtrlmm(kc, nc, alpha, A_, B_,
                        (B_start + l * BlockSizes.MC * incRowB + j * BlockSizes.NC * incColB), // B
                                                                                                // start
                        B, incRowB, incColB, AB, C_);

                for (int i = l + 1; i < mb; ++i) {
                    int mc = (i != mb - 1 || mc_ == 0) ? BlockSizes.MC : mc_;

                    Gepack.gepack_A(mc, kc, (A_start + i * BlockSizes.MC * incRowA + l * BlockSizes.MC * incColA), // A
                                                                                                                    // start
                            A, incRowA, incColA, A_);

                    Mgemm.mgemm(mc, nc, kc, alpha, A_, B_, 1.0,
                            (B_start + i * BlockSizes.MC * incRowB + j * BlockSizes.NC * incColB), // B
                                                                                                    // start
                            B, incRowB, incColB, AB, C_);
                }
            }
        }

    }

    private TrlmmLehn() {
        throw new AssertionError();
    }
}
