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
package math.gemm;

import java.util.Arrays;

/**
 * A straightforward Java translation of Michael Lehn's pure ANSI C variant of a
 * cache-friendly dgemm routine.
 * <p>
 * https://apfel.mathematik.uni-ulm.de/~lehn/sghpc/day08/page06.html
 * <p>
 * http://apfel.mathematik.uni-ulm.de/~lehn/sghpc/gemm/page02/index.html
 * <p>
 * https://apfel.mathematik.uni-ulm.de/~lehn/sghpc/day08/page02.html
 * <p>
 * https://apfel.mathematik.uni-ulm.de/~lehn/sghpc/day08/page03.html
 * <p>
 * https://apfel.mathematik.uni-ulm.de/~lehn/sghpc/day08/page04.html
 * <p>
 * https://stackoverflow.com/questions/1303182/how-does-blas-get-such-extreme-performance/11421344#11421344
 * <p>
 * http://apfel.mathematik.uni-ulm.de/~lehn/sghpc/gemm/
 * <p>
 * https://stackoverflow.com/questions/35620853/how-to-write-a-matrix-matrix-product-that-can-compete-with-eigen/35637007#35637007
 */
final class Dgemm4x4 {

    static final int MR_Height = 4; // 4
    static final int NR_Width = 4; // 4

    static final int MC = 384; // 384
    static final int KC = 384; // 384
    static final int NC = 4096; // 4096 .. 16384

    //
    // Packing complete panels from A (i.e. without padding)
    //
    private static void pack_A_MRxk(int k, int A_start, double[] A, int incRowA, int incColA, double[] work,
            int work_start) {

        for (int j = 0; j < k; ++j) {
            for (int i = 0; i < MR_Height; ++i) {
                work[work_start + i] = A[A_start + i * incRowA];
            }
            work_start += MR_Height;
            A_start += incColA;
        }
    }

    //
    // Packing complete panels from B (i.e. without padding)
    //
    private static void pack_B_kxNR(int k, int B_start, double[] B, int incRowB, int incColB, double[] work,
            int work_start) {

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < NR_Width; ++j) {
                work[work_start + j] = B[B_start + j * incColB];
            }
            work_start += NR_Width;
            B_start += incRowB;
        }
    }

    //
    // Packing panels from A with padding if required
    //
    static void pack_A(int mc, int kc, int A_start, double[] A, int incRowA, int incColA, double[] work) {

        final int mp = mc / MR_Height;
        final int _mr = mc % MR_Height;

        int work_start = 0;
        int i;

        for (i = 0; i < mp; ++i) {
            pack_A_MRxk(kc, A_start, A, incRowA, incColA, work, work_start);
            work_start += kc * MR_Height;
            A_start += MR_Height * incRowA;
        }
        if (_mr > 0) {
            for (int j = 0; j < kc; ++j) {
                for (i = 0; i < _mr; ++i) {
                    work[work_start + i] = A[A_start + i * incRowA];
                }
                for (i = _mr; i < MR_Height; ++i) {
                    work[work_start + i] = 0.0;
                }
                work_start += MR_Height;
                A_start += incColA;
            }
        }
    }

    //
    // Packing panels from B with padding if required
    //
    static void pack_B(int kc, int nc, int B_start, double[] B, int incRowB, int incColB, double[] work) {

        final int np = nc / NR_Width;
        final int _nr = nc % NR_Width;

        int work_start = 0;
        int j;

        for (j = 0; j < np; ++j) {
            pack_B_kxNR(kc, B_start, B, incRowB, incColB, work, work_start);
            work_start += kc * NR_Width;
            B_start += NR_Width * incColB;
        }
        if (_nr > 0) {
            for (int i = 0; i < kc; ++i) {
                for (j = 0; j < _nr; ++j) {
                    work[work_start + j] = B[B_start + j * incColB];
                }
                for (j = _nr; j < NR_Width; ++j) {
                    work[work_start + j] = 0.0;
                }
                work_start += NR_Width;
                B_start += incRowB;
            }
        }
    }

    //
    // Micro kernel for multiplying panels from A and B
    //
    private static void dgemm_micro_kernel(int kc, double alpha, int A_panel_start, double[] A_panel, int B_panel_start,
            double[] B_panel, double beta, final int C_panel_start, double[] C_panel, int incRowC, int incColC,
            double[] AB) {

        // clear buffer
        Arrays.fill(AB, 0.0);

        //
        // Compute AB = A*B
        //
        for (int l = 0; l < kc; ++l) {

            double b0j = B_panel[B_panel_start++];
            double b1j = B_panel[B_panel_start++];
            double b2j = B_panel[B_panel_start++];
            double b3j = B_panel[B_panel_start++];

            double a0i = A_panel[A_panel_start++];
            double a1i = A_panel[A_panel_start++];
            double a2i = A_panel[A_panel_start++];
            double a3i = A_panel[A_panel_start++];

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
            dgemm_micro_betaMulC(beta, C_panel_start, C_panel, incRowC, incColC);
        }

        //
        // Update C <- C + alpha*AB (note: the case alpha==0.0 was already
        // treated in the above layer dgemm)
        //
        dgemm_micro_plusAlphaAB(alpha, C_panel_start, C_panel, incRowC, incColC, AB);
    }

    //
    // Update C <- beta*C
    //
    private static void dgemm_micro_betaMulC(double beta, int C_panel_start, double[] C_panel, int incRowC,
            int incColC) {

        if (beta == 0.0) {
            for (int j = 0; j < NR_Width; ++j) {
                int base_C = C_panel_start + j * incColC;
                for (int i = 0; i < MR_Height; ++i) {
                    C_panel[base_C + i * incRowC] = 0.0;
                }
            }
        } else {
            for (int j = 0; j < NR_Width; ++j) {
                int base_C = C_panel_start + j * incColC;
                for (int i = 0; i < MR_Height; ++i) {
                    C_panel[base_C + i * incRowC] *= beta;
                }
            }
        }
    }

    //
    // Update C <- C + alpha*AB (note: the case alpha==0.0 was already
    // treated in the above layer dgemm)
    //
    private static void dgemm_micro_plusAlphaAB(double alpha, int C_panel_start, double[] C_panel, int incRowC,
            int incColC, double[] AB) {

        if (alpha == 1.0) {
            for (int j = 0; j < NR_Width; ++j) {
                int jIdx = j * MR_Height;
                int base_C = C_panel_start + j * incColC;
                for (int i = 0; i < MR_Height; ++i) {
                    C_panel[base_C + i * incRowC] += AB[i + jIdx];
                }
            }
        } else {
            for (int j = 0; j < NR_Width; ++j) {
                int jIdx = j * MR_Height;
                int base_C = C_panel_start + j * incColC;
                for (int i = 0; i < MR_Height; ++i) {
                    C_panel[base_C + i * incRowC] += alpha * AB[i + jIdx];
                }
            }
        }
    }

    //
    // Macro Kernel for the multiplication of blocks of A and B. We assume that
    // these blocks were previously packed to buffers _A and _B.
    //
    static int dgemm_macro_kernel(int mc, int nc, int kc, double alpha, double beta, int C_start, double[] C,
            int incRowC, int incColC, double[] _A, double[] _B, double[] AB, double[] workC) {

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
                    dgemm_micro_kernel(kc, alpha, (i * kc * MR_Height), _A, (j * kc * NR_Width), _B, beta,
                            (C_start + i * MR_Height * incRowC + j * NR_Width * incColC), C, incRowC, incColC, AB);
                    ++micro_kernel_calls;
                } else {
                    dgemm_micro_kernel(kc, alpha, (i * kc * MR_Height), _A, (j * kc * NR_Width), _B, 0.0, 0, workC, 1,
                            MR_Height, AB);

                    dgescal(mr, nr, beta, (C_start + i * MR_Height * incRowC + j * NR_Width * incColC), C, incRowC,
                            incColC);

                    dgeaxpy(mr, nr, 1.0, workC, 1, MR_Height,
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
    private static void dgescal(int m, int n, double alpha, int X_start, double[] X, int incRowX, int incColX) {

        if (alpha != 0.0) {
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
                    X[base_X + i * incRowX] = 0.0;
                }
            }
        }
    }

    //
    // Compute Y += alpha*X
    //
    private static void dgeaxpy(int m, int n, double alpha, double[] X, int incRowX, int incColX, int Y_start,
            double[] Y, int incRowY, int incColY) {

        if (alpha != 1.0) {
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
    static int dgemm(int rowsA, int colsB, int colsA, double alpha, int offA, double[] A, int incRowA, int incColA,
            int offB, double[] B, int incRowB, int incColB, double beta, int offC, double[] C, int incRowC,
            int incColC) {

        int micro_kernel_calls = 0;

        if (alpha == 0.0 || colsA == 0) {
            dgescal(rowsA, colsB, beta, offC, C, incRowC, incColC);
            return micro_kernel_calls;
        }

        final int mb = (rowsA + MC - 1) / MC;
        final int nb = (colsB + NC - 1) / NC;
        final int kb = (colsA + KC - 1) / KC;

        final int _mc = rowsA % MC;
        final int _nc = colsB % NC;
        final int _kc = colsA % KC;

        // check whether we can parallelize the computation
        if ((mb > 1 || nb > 1) && DgemmTasks.availableCores() > 1) {
            return Dgemm4x4Parallel.dgemm(nb, kb, mb, _nc, _kc, _mc, alpha, offA, A, incRowA, incColA, offB, B, incRowB,
                    incColB, beta, offC, C, incRowC, incColC);
        }

        //
        // Local buffers for storing panels from A, B and C
        //
        final double[] _A = new double[MC * KC];
        final double[] _B = new double[KC * NC];
        final double[] _C = new double[MR_Height * NR_Width];
        final double[] AB = new double[MR_Height * NR_Width];

        for (int j = 0; j < nb; ++j) {
            int nc = (j != nb - 1 || _nc == 0) ? NC : _nc;

            for (int l = 0; l < kb; ++l) {
                int kc = (l != kb - 1 || _kc == 0) ? KC : _kc;
                double _beta = (l == 0) ? beta : 1.0;

                pack_B(kc, nc, (offB + l * KC * incRowB + j * NC * incColB), B, incRowB, incColB, _B);

                for (int i = 0; i < mb; ++i) {
                    int mc = (i != mb - 1 || _mc == 0) ? MC : _mc;

                    pack_A(mc, kc, (offA + i * MC * incRowA + l * KC * incColA), A, incRowA, incColA, _A);

                    micro_kernel_calls += dgemm_macro_kernel(mc, nc, kc, alpha, _beta,
                            (offC + i * MC * incRowC + j * NC * incColC), C, incRowC, incColC, _A, _B, AB, _C);
                }
            }
        }
        return micro_kernel_calls;
    }
}
