/*
 * Copyright 2018, 2026 Stefan Zobel
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
package math.gemm;

import java.util.Arrays;


/**
 * A straightforward Java translation of Michael Lehn's pure ANSI C variant of a
 * cache-friendly BLIS sgemm routine.
 */
public final class SgemmMRxNR {

    // float
    private static final int MR_Height = 4; // 8 // that's extremely frustrating!
    private static final int NR_Width = 6;  // 16

    private static final int KC = 256;  // 512;  // 512 .. 1024 // that's extremely frustrating!
    private static final int MC = 96;   // 192;  // 128 .. 256 (256) // that's extremely frustrating!
    private static final int NC = 3840; // 3840 or more

    //
    // Packing complete panels from A (i.e. without padding)
    //
    private static void pack_A_4xk(int k, int A_start, float[] A, int incRowA, int incColA, float[] work,
            int work_start) {

        for (int j = 0; j < k; ++j) {
            // Packing manually unrolled
            work[work_start + 0] = A[A_start + 0 * incRowA];
            work[work_start + 1] = A[A_start + 1 * incRowA];
            work[work_start + 2] = A[A_start + 2 * incRowA];
            work[work_start + 3] = A[A_start + 3 * incRowA];

            work_start += MR_Height;
            A_start += incColA;
        }
    }

    private static void pack_A_4xk_fast(int k, int A_start, float[] A, int incColA, float[] work,
            int work_start) {

        for (int j = 0; j < k; ++j) {
            // Packing manually unrolled
            work[work_start + 0] = A[A_start + 0];
            work[work_start + 1] = A[A_start + 1];
            work[work_start + 2] = A[A_start + 2];
            work[work_start + 3] = A[A_start + 3];

            work_start += MR_Height;
            A_start += incColA;
        }
    }

    //
    // Packing complete panels from B (i.e. without padding)
    //
    private static void pack_B_kx6(int k, int B_start, float[] B, int incRowB, int incColB, float[] work,
            int work_start) {

        for (int i = 0; i < k; ++i) {
            // Packing manually unrolled
            work[work_start + 0] = B[B_start + 0 * incColB];
            work[work_start + 1] = B[B_start + 1 * incColB];
            work[work_start + 2] = B[B_start + 2 * incColB];
            work[work_start + 3] = B[B_start + 3 * incColB];
            work[work_start + 4] = B[B_start + 4 * incColB];
            work[work_start + 5] = B[B_start + 5 * incColB];

            work_start += NR_Width;
            B_start += incRowB;
        }
    }

    private static void pack_B_kx6_fast(int k, int B_start, float[] B, int incRowB, float[] work,
            int work_start) {

        for (int i = 0; i < k; ++i) {
            // Packing manually unrolled
            work[work_start + 0] = B[B_start + 0];
            work[work_start + 1] = B[B_start + 1];
            work[work_start + 2] = B[B_start + 2];
            work[work_start + 3] = B[B_start + 3];
            work[work_start + 4] = B[B_start + 4];
            work[work_start + 5] = B[B_start + 5];

            work_start += NR_Width;
            B_start += incRowB;
        }
    }

    //
    // Packing panels from A with padding if required
    //
    static void pack_A(int mc, int kc, int A_start, float[] A, int incRowA, int incColA, float[] work) {

        final int mp = mc / MR_Height;
        final int _mr = mc % MR_Height;

        int work_start = 0;
        int i;

        for (i = 0; i < mp; ++i) {
            if (incRowA == 1) {
                pack_A_4xk_fast(kc, A_start, A, incColA, work, work_start);
            } else {
                pack_A_4xk(kc, A_start, A, incRowA, incColA, work, work_start);
            }
            work_start += kc * MR_Height;
            A_start += MR_Height * incRowA;
        }
        if (_mr > 0) {
            for (int j = 0; j < kc; ++j) {
                for (i = 0; i < _mr; ++i) {
                    work[work_start + i] = A[A_start + i * incRowA];
                }
                for (i = _mr; i < MR_Height; ++i) {
                    work[work_start + i] = 0.0f;
                }
                work_start += MR_Height;
                A_start += incColA;
            }
        }
    }

    //
    // Packing panels from B with padding if required
    //
    static void pack_B(int kc, int nc, int B_start, float[] B, int incRowB, int incColB, float[] work) {

        final int np = nc / NR_Width;
        final int _nr = nc % NR_Width;

        int work_start = 0;
        int j;

        for (j = 0; j < np; ++j) {
            if (incColB == 1) {
                pack_B_kx6_fast(kc, B_start, B, incRowB, work, work_start);
            } else {
                pack_B_kx6(kc, B_start, B, incRowB, incColB, work, work_start);
            }
            work_start += kc * NR_Width;
            B_start += NR_Width * incColB;
        }
        if (_nr > 0) {
            for (int i = 0; i < kc; ++i) {
                for (j = 0; j < _nr; ++j) {
                    work[work_start + j] = B[B_start + j * incColB];
                }
                for (j = _nr; j < NR_Width; ++j) {
                    work[work_start + j] = 0.0f;
                }
                work_start += NR_Width;
                B_start += incRowB;
            }
        }
    }

    //
    // Micro kernel for multiplying panels from A and B
    //
    private static void sgemm_micro_kernel(int kc, float alpha, int A_panel_start, float[] A_panel, int B_panel_start,
            float[] B_panel, float beta, final int C_panel_start, float[] C_panel, int incRowC, int incColC,
            float[] AB) {

        // clear buffer
        Arrays.fill(AB, 0.0f);

        //
        // Compute AB = A*B
        //
        for (int l = 0; l < kc; ++l) {

            float b0j = B_panel[B_panel_start++];
            float b1j = B_panel[B_panel_start++];
            float b2j = B_panel[B_panel_start++];
            float b3j = B_panel[B_panel_start++];
            float b4j = B_panel[B_panel_start++];
            float b5j = B_panel[B_panel_start++];

            float a0i = A_panel[A_panel_start++];
            float a1i = A_panel[A_panel_start++];
            float a2i = A_panel[A_panel_start++];
            float a3i = A_panel[A_panel_start++];

            int idx = 0;
            // Column 0 of B
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

            AB[idx++] += a0i * b4j;
            AB[idx++] += a1i * b4j;
            AB[idx++] += a2i * b4j;
            AB[idx++] += a3i * b4j;

            // Column 5 of B
            AB[idx++] += a0i * b5j;
            AB[idx++] += a1i * b5j;
            AB[idx++] += a2i * b5j;
            AB[idx++] += a3i * b5j;
        }

        //
        // Update C <- beta*C
        //
        if (beta != 1.0f) {
            sgemm_micro_betaMulC_MR4(beta, C_panel_start, C_panel, incColC);
        }

        //
        // Update C <- C + alpha*AB (note: the case alpha==0.0 was already
        // treated in the above layer dgemm)
        //
        sgemm_micro_plusAlphaAB_MR4(alpha, C_panel_start, C_panel, incColC, AB);
    }

    //
    // Update C <- beta*C
    //
    private static void sgemm_micro_betaMulC_MR4(float beta, int C_panel_start, float[] C_panel,
            int incColC) {

        if (beta == 0.0f) {
            for (int j = 0; j < NR_Width; ++j) {
                int base_C = C_panel_start + j * incColC;
                C_panel[base_C + 0] = 0.0f;
                C_panel[base_C + 1] = 0.0f;
                C_panel[base_C + 2] = 0.0f;
                C_panel[base_C + 3] = 0.0f;
            }
        } else {
            for (int j = 0; j < NR_Width; ++j) {
                int base_C = C_panel_start + j * incColC;
                C_panel[base_C + 0] *= beta;
                C_panel[base_C + 1] *= beta;
                C_panel[base_C + 2] *= beta;
                C_panel[base_C + 3] *= beta;
            }
        }
    }

    //
    // Update C <- C + alpha*AB (note: the case alpha==0.0f was already
    // treated in the above layer sgemm)
    //
    private static void sgemm_micro_plusAlphaAB_MR4(float alpha, int C_panel_start, float[] C_panel,
            int incColC, float[] AB) {

        if (alpha == 1.0f) {
            for (int j = 0; j < NR_Width; ++j) {
                int jIdx = j * MR_Height;
                int base_C = C_panel_start + j * incColC;
                C_panel[base_C + 0] += AB[0 + jIdx];
                C_panel[base_C + 1] += AB[1 + jIdx];
                C_panel[base_C + 2] += AB[2 + jIdx];
                C_panel[base_C + 3] += AB[3 + jIdx];
            }
        } else {
            for (int j = 0; j < NR_Width; ++j) {
                int jIdx = j * MR_Height;
                int base_C = C_panel_start + j * incColC;
                C_panel[base_C + 0] += alpha * AB[0 + jIdx];
                C_panel[base_C + 1] += alpha * AB[1 + jIdx];
                C_panel[base_C + 2] += alpha * AB[2 + jIdx];
                C_panel[base_C + 3] += alpha * AB[3 + jIdx];
            }
        }
    }

    //
    // Macro Kernel for the multiplication of blocks of A and B. We assume that
    // these blocks were previously packed to buffers _A and _B.
    //
    static int sgemm_macro_kernel(int mc, int nc, int kc, float alpha, float beta, int C_start, float[] C,
            int incRowC, int incColC, float[] _A, float[] _B, float[] AB, float[] workC) {

        int micro_kernel_calls = 0;

        final int mp = (mc + MR_Height - 1) / MR_Height;
        final int np = (nc + NR_Width - 1) / NR_Width;

        final int _mr = mc % MR_Height;
        final int _nr = nc % NR_Width;

        int nr, mr;
        int i, j;

        for (j = 0; j < np; ++j) {
            nr = (j != np - 1 || _nr == 0) ? NR_Width : _nr;

            for (i = 0; i < mp; ++i) {
                mr = (i != mp - 1 || _mr == 0) ? MR_Height : _mr;

                if (mr == MR_Height && nr == NR_Width) {
                    sgemm_micro_kernel(kc, alpha, (i * kc * MR_Height), _A, (j * kc * NR_Width), _B, beta,
                            (C_start + i * MR_Height * incRowC + j * NR_Width * incColC), C, incRowC, incColC, AB);
                    ++micro_kernel_calls;
                } else {
                    sgemm_micro_kernel(kc, alpha, (i * kc * MR_Height), _A, (j * kc * NR_Width), _B, 0.0f, 0, workC, 1,
                            MR_Height, AB);

                    sgescal(mr, nr, beta, (C_start + i * MR_Height * incRowC + j * NR_Width * incColC), C, incRowC,
                            incColC);

                    sgeaxpy(mr, nr, 1.0f, workC, 1, MR_Height,
                            (C_start + i * MR_Height * incRowC + j * NR_Width * incColC), C, incRowC, incColC);
                    ++micro_kernel_calls;
                }
            }
        }
        return micro_kernel_calls;
    }

    //
    // Compute X *= alpha
    //
    private static void sgescal(int m, int n, float alpha, int X_start, float[] X, int incRowX, int incColX) {

        if (alpha != 0.0f) {
            for (int j = 0; j < n; ++j) {
                int base_X = X_start + j * incColX;
                for (int i = 0; i < m; ++i) {
                    X[base_X + i * incRowX] *= alpha;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                int base_X = X_start + j * incColX;
                for (int i = 0; i < m; ++i) {
                    X[base_X + i * incRowX] = 0.0f;
                }
            }
        }
    }

    //
    // Compute Y += alpha*X
    //
    private static void sgeaxpy(int m, int n, float alpha, float[] X, int incRowX, int incColX, int Y_start,
            float[] Y, int incRowY, int incColY) {

        if (alpha != 1.0f) {
            for (int j = 0; j < n; ++j) {
                int base_Y = Y_start + j * incColY;
                int _incColX = j * incColX;
                for (int i = 0; i < m; ++i) {
                    Y[base_Y + i * incRowY] += alpha * X[i * incRowX + _incColX];
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                int base_Y = Y_start + j * incColY;
                int _incColX = j * incColX;
                for (int i = 0; i < m; ++i) {
                    Y[base_Y + i * incRowY] += X[i * incRowX + _incColX];
                }
            }
        }
    }

    //
    // Compute C <- beta*C + alpha*A*B
    //
    public static int sgemm(int rowsA, int colsB, int colsA, float alpha, int offA, float[] A, int incRowA, int incColA,
            int offB, float[] B, int incRowB, int incColB, float beta, int offC, float[] C, int incRowC,
            int incColC) {

        int micro_kernel_calls = 0;

        if (alpha == 0.0f || colsA == 0) {
            sgescal(rowsA, colsB, beta, offC, C, incRowC, incColC);
            return micro_kernel_calls;
        }

        final int mb = (rowsA + MC - 1) / MC;
        final int nb = (colsB + NC - 1) / NC;
        final int kb = (colsA + KC - 1) / KC;

        final int _mc = rowsA % MC;
        final int _nc = colsB % NC;
        final int _kc = colsA % KC;

        //
        // Local buffers for storing panels from A, B and C
        //
        final float[] _A = new float[MC * KC];
        final float[] _B = new float[KC * NC];
        final float[] _C = new float[MR_Height * NR_Width];
        final float[] AB = new float[MR_Height * NR_Width];

        for (int j = 0; j < nb; ++j) {
            int nc = (j != nb - 1 || _nc == 0) ? NC : _nc;

            for (int l = 0; l < kb; ++l) {
                int kc = (l != kb - 1 || _kc == 0) ? KC : _kc;
                float _beta = (l == 0) ? beta : 1.0f;

                pack_B(kc, nc, (offB + l * KC * incRowB + j * NC * incColB), B, incRowB, incColB, _B);

                for (int i = 0; i < mb; ++i) {
                    int mc = (i != mb - 1 || _mc == 0) ? MC : _mc;

                    pack_A(mc, kc, (offA + i * MC * incRowA + l * KC * incColA), A, incRowA, incColA, _A);

                    micro_kernel_calls += sgemm_macro_kernel(mc, nc, kc, alpha, _beta,
                            (offC + i * MC * incRowC + j * NC * incColC), C, incRowC, incColC, _A, _B, AB, _C);
                }
            }
        }
        return micro_kernel_calls;
    }
}
