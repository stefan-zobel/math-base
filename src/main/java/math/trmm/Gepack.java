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

final class Gepack {

    static void gepack_A(int mc, int kc, int A_start, double[] A, int incRowA, int incColA, double[] buffer) {
        final int mp = mc / BlockSizes.MR;
        final int mr_ = (mc & (BlockSizes.MR - 1));

        int buffer_start = 0;

        for (int i = 0; i < mp; ++i) {
            packA_MRxk(kc, A_start, A, incRowA, incColA, buffer, buffer_start);
            buffer_start += kc * BlockSizes.MR;
            A_start += BlockSizes.MR * incRowA;
        }
        if (mr_ > 0) {
            for (int j = 0; j < kc; ++j) {
                for (int i = 0; i < mr_; ++i) {
                    buffer[buffer_start + i] = A[A_start + i * incRowA];
                }
                for (int i = mr_; i < BlockSizes.MR; ++i) {
                    buffer[buffer_start + i] = 0.0;
                }
                buffer_start += BlockSizes.MR;
                A_start += incColA;
            }
        }
    }

    public static void gepack_B(int kc, int nc, int B_start, double[] B, int incRowB, int incColB, double[] buffer) {
        final int np = nc / BlockSizes.NR;
        final int nr_ = (nc & (BlockSizes.NR - 1));

        int buffer_start = 0;

        for (int j = 0; j < np; ++j) {
            packB_kxNR(kc, B_start, B, incRowB, incColB, buffer, buffer_start);
            buffer_start += kc * BlockSizes.NR;
            B_start += BlockSizes.NR * incColB;
        }
        if (nr_ > 0) {
            for (int i = 0; i < kc; ++i) {
                for (int j = 0; j < nr_; ++j) {
                    buffer[buffer_start + j] = B[B_start + j * incColB];
                }
                for (int j = nr_; j < BlockSizes.NR; ++j) {
                    buffer[buffer_start + j] = 0.0;
                }
                buffer_start += BlockSizes.NR;
                B_start += incRowB;
            }
        }
    }

    private static void packB_kxNR(int k, int B_start, double[] B, int incRowB, int incColB, double[] work,
            int work_start) {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < BlockSizes.NR; ++j) {
                work[work_start + j] = B[B_start + j * incColB];
            }
            work_start += BlockSizes.NR;
            B_start += incRowB;
        }
    }

    private static void packA_MRxk(int k, int A_start, double[] A, int incRowA, int incColA, double[] work,
            int work_start) {
        for (int j = 0; j < k; ++j) {
            for (int i = 0; i < BlockSizes.MR; ++i) {
                work[work_start + i] = A[A_start + i * incRowA];
            }
            work_start += BlockSizes.MR;
            A_start += incColA;
        }
    }

    private Gepack() {
        throw new AssertionError();
    }
}
