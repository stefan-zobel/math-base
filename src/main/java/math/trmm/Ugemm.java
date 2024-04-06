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

import java.util.Arrays;

final class Ugemm {

    //
    // Buffered variant. Used for zero padded panels.
    //
    static void ugemm(int mr, int nr, int kc, double alpha, int A_start, double[] A, int B_start, double[] B,
            double beta, int C_start, double[] C, int incRowC, int incColC, double[] AB, double[] workC) {

        ugemm(kc, alpha, A_start, A, B_start, B, 0.0, 0, workC, 1, BlockSizes.MR, AB);
        Gescal.gescal(mr, nr, beta, C_start, C, incRowC, incColC);
        Geaxpy.geaxpy(mr, nr, 1.0, 0, workC, 1, BlockSizes.MR, C_start, C, incRowC, incColC);
    }

    //
    // Micro kernel for multiplying panels from A and B. Unbuffered variant.
    //
    static void ugemm(int kc, double alpha, int A_start, double[] A, int B_start, double[] B, double beta,
            int C_start, double[] C, int incRowC, int incColC, double[] AB) {

        // clear buffer
        Arrays.fill(AB, 0.0);

        //
        // Compute AB = A*B
        //
        for (int l = 0; l < kc; ++l) {

            double b0j = B[B_start++];
            double b1j = B[B_start++];
            double b2j = B[B_start++];
            double b3j = B[B_start++];

            double a0i = A[A_start++];
            double a1i = A[A_start++];
            double a2i = A[A_start++];
            double a3i = A[A_start++];

            int idx = 0;
            AB[idx++] += a0i * b0j;
            AB[idx++] += a1i * b0j;
            AB[idx++] += a2i * b0j;
            AB[idx++] += a3i * b0j;

            AB[idx++] += a0i * b1j;
            AB[idx++] += a1i * b1j;
            AB[idx++] += a2i * b1j;
            AB[idx++] += a3i * b1j;

            AB[idx++] += a0i * b2j;
            AB[idx++] += a1i * b2j;
            AB[idx++] += a2i * b2j;
            AB[idx++] += a3i * b2j;

            AB[idx++] += a0i * b3j;
            AB[idx++] += a1i * b3j;
            AB[idx++] += a2i * b3j;
            AB[idx++] += a3i * b3j;
        }

        //
        // Update C <- beta*C
        //
        if (beta != 1.0) {
            betaMulC(beta, C_start, C, incRowC, incColC);
        }

        //
        // Update C <- C + alpha*AB (note: the case alpha==0.0 must
        // already have been treated in one of the the above layers!)
        //
        plusAlphaAB(alpha, C_start, C, incRowC, incColC, AB);
    }

    //
    // Update C <- beta*C
    //
    private static void betaMulC(double beta, int C_start, double[] C, int incRowC, int incColC) {

        if (beta == 0.0) {
            for (int j = 0; j < BlockSizes.NR; ++j) {
                int base_C = C_start + j * incColC;
                for (int i = 0; i < BlockSizes.MR; ++i) {
                    C[base_C + i * incRowC] = 0.0;
                }
            }
        } else {
            for (int j = 0; j < BlockSizes.NR; ++j) {
                int base_C = C_start + j * incColC;
                for (int i = 0; i < BlockSizes.MR; ++i) {
                    C[base_C + i * incRowC] *= beta;
                }
            }
        }
    }

    //
    // Update C <- C + alpha*AB (note: the case alpha==0.0 must
    // already have been treated in one of the the above layers!)
    //
    private static void plusAlphaAB(double alpha, int C_start, double[] C, int incRowC, int incColC, double[] AB) {

        if (alpha == 1.0) {
            for (int j = 0; j < BlockSizes.NR; ++j) {
                int jIdx = j * BlockSizes.MR;
                int base_C = C_start + j * incColC;
                for (int i = 0; i < BlockSizes.MR; ++i) {
                    C[base_C + i * incRowC] += AB[i + jIdx];
                }
            }
        } else {
            for (int j = 0; j < BlockSizes.NR; ++j) {
                int jIdx = j * BlockSizes.MR;
                int base_C = C_start + j * incColC;
                for (int i = 0; i < BlockSizes.MR; ++i) {
                    C[base_C + i * incRowC] += alpha * AB[i + jIdx];
                }
            }
        }
    }

    private Ugemm() {
        throw new AssertionError();
    }
}
