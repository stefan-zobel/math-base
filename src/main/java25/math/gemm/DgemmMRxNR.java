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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;


/**
 * A straightforward Java translation of Michael Lehn's pure ANSI C variant of a
 * cache-friendly BLIS dgemm routine.
 */
public final class DgemmMRxNR {

    // double
    private static final int MR_Height = 4; // 4
    private static final int NR_Width = 8;  // 8

    private static final int KC = 256;  // 256 .. 512
    private static final int MC = 128;  // 128 .. 256 (128)
    private static final int NC = 3840; // 3840 or more

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    //
    // Packing complete panels from A (i.e. without padding)
    //
    private static void pack_A_4xk(int k, int A_start, double[] A, int incRowA, int incColA, double[] work,
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

    private static void pack_A_4xk_fast(int k, int A_start, double[] A, int incColA, double[] work,
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
    private static void pack_B_kx8(int k, int B_start, double[] B, int incRowB, int incColB, double[] work,
            int work_start) {

        for (int i = 0; i < k; ++i) {
            // Packing manually unrolled
            work[work_start + 0] = B[B_start + 0 * incColB];
            work[work_start + 1] = B[B_start + 1 * incColB];
            work[work_start + 2] = B[B_start + 2 * incColB];
            work[work_start + 3] = B[B_start + 3 * incColB];
            work[work_start + 4] = B[B_start + 4 * incColB];
            work[work_start + 5] = B[B_start + 5 * incColB];
            work[work_start + 6] = B[B_start + 6 * incColB];
            work[work_start + 7] = B[B_start + 7 * incColB];

            work_start += NR_Width;
            B_start += incRowB;
        }
    }

    private static void pack_B_kx8_fast(int k, int B_start, double[] B, int incRowB, double[] work,
            int work_start) {

        for (int i = 0; i < k; ++i) {
            // Packing manually unrolled
            work[work_start + 0] = B[B_start + 0];
            work[work_start + 1] = B[B_start + 1];
            work[work_start + 2] = B[B_start + 2];
            work[work_start + 3] = B[B_start + 3];
            work[work_start + 4] = B[B_start + 4];
            work[work_start + 5] = B[B_start + 5];
            work[work_start + 6] = B[B_start + 6];
            work[work_start + 7] = B[B_start + 7];

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
            if (incColB == 1) {
                pack_B_kx8_fast(kc, B_start, B, incRowB, work, work_start);
            } else {
                pack_B_kx8(kc, B_start, B, incRowB, incColB, work, work_start);
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

        // 1. Initialize 8 accumulator registers (vectors) with 0
        DoubleVector c0 = DoubleVector.zero(SPECIES);
        DoubleVector c1 = DoubleVector.zero(SPECIES);
        DoubleVector c2 = DoubleVector.zero(SPECIES);
        DoubleVector c3 = DoubleVector.zero(SPECIES);
        DoubleVector c4 = DoubleVector.zero(SPECIES);
        DoubleVector c5 = DoubleVector.zero(SPECIES);
        DoubleVector c6 = DoubleVector.zero(SPECIES);
        DoubleVector c7 = DoubleVector.zero(SPECIES);

        // 2. Main lopp (mostly compute operations)
        for (int l = 0; l < kc; ++l) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, A_panel, A_panel_start);

            // Broadcast the 8 B values
            c0 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 0]), c0);
            c1 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 1]), c1);
            c2 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 2]), c2);
            c3 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 3]), c3);
            c4 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 4]), c4);
            c5 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 5]), c5);
            c6 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 6]), c6);
            c7 = va.fma(DoubleVector.broadcast(SPECIES, B_panel[B_panel_start + 7]), c7);

            A_panel_start += MR_Height; // 4
            B_panel_start += NR_Width;  // 8
        }

        // 3. Write the results to the AB array
        // This time the increment is 4 since we can only store 4 doubles
        c0.intoArray(AB, 0);
        c1.intoArray(AB, 4);
        c2.intoArray(AB, 8);
        c3.intoArray(AB, 12);
        c4.intoArray(AB, 16);
        c5.intoArray(AB, 20);
        c6.intoArray(AB, 24);
        c7.intoArray(AB, 28);

        //
        // Update C <- beta*C
        //
        if (beta != 1.0) {
            dgemm_micro_betaMulC_MR4(beta, C_panel_start, C_panel, incColC);
        }

        //
        // Update C <- C + alpha*AB (note: the case alpha==0.0 was already
        // treated in the above layer dgemm)
        //
        dgemm_micro_plusAlphaAB_MR4(alpha, C_panel_start, C_panel, incColC, AB);
    }

    //
    // Update C <- beta*C
    //
    private static void dgemm_micro_betaMulC_MR4(double beta, int C_panel_start, double[] C_panel,
            int incColC) {

        if (beta == 0.0) {
            for (int j = 0; j < NR_Width; ++j) {
                int base_C = C_panel_start + j * incColC;
                C_panel[base_C + 0] = 0.0;
                C_panel[base_C + 1] = 0.0;
                C_panel[base_C + 2] = 0.0;
                C_panel[base_C + 3] = 0.0;
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
    // Update C <- C + alpha*AB (note: the case alpha==0.0 was already
    // treated in the above layer dgemm)
    //
    private static void dgemm_micro_plusAlphaAB_MR4(double alpha, int C_panel_start, double[] C_panel,
            int incColC, double[] AB) {

        if (alpha == 1.0) {
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
    // Compute C <- beta*C + alpha*A*B (sequential, full column range)
    //
    public static int dgemm(int rowsA, int colsB, int colsA, double alpha, int offA, double[] A, int incRowA, int incColA,
            int offB, double[] B, int incRowB, int incColB, double beta, int offC, double[] C, int incRowC,
            int incColC) {

        if (alpha == 0.0 || colsA == 0) {
            dgescal(rowsA, colsB, beta, offC, C, incRowC, incColC);
            return 0;
        }
        return dgemmRange(rowsA, 0, colsB, colsA, alpha, offA, A, incRowA, incColA,
                offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC);
    }

    //
    // Compute C <- beta*C + alpha*A*B, optionally parallelized over the N dimension.
    // Falls back to the sequential path when executor is null, shut down, or the
    // workload heuristic decides against parallel execution.
    //
    public static int dgemm(int rowsA, int colsB, int colsA, double alpha, int offA, double[] A, int incRowA, int incColA,
            int offB, double[] B, int incRowB, int incColB, double beta, int offC, double[] C, int incRowC,
            int incColC, ExecutorService executor) {

        if (alpha == 0.0 || colsA == 0) {
            dgescal(rowsA, colsB, beta, offC, C, incRowC, incColC);
            return 0;
        }
        long work = (long) rowsA * (long) colsB * (long) colsA;
        if (shouldParallelize(executor, colsB, work)) {
            parallelizeOverN(executor, rowsA, colsB, colsA, alpha, offA, A, incRowA, incColA,
                    offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC);
            return 0;
        }
        return dgemmRange(rowsA, 0, colsB, colsA, alpha, offA, A, incRowA, incColA,
                offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC);
    }

    public static boolean isVectorized() {
        return true;
    }

    public static boolean isVectorApiPresent() {
        return ModuleLayer.boot().findModule("jdk.incubator.vector").isPresent();
    }

    // The Vector API micro-kernel is faster than the Java-8 scalar kernel, so
    // parallelism pays off at a lower work volume.
    private static final long PARALLEL_WORK_THRESHOLD = 25_000_000L;

    static long getParallelWorkThreshold() {
        return PARALLEL_WORK_THRESHOLD;
    }

    private static boolean shouldParallelize(ExecutorService executor, int colsB, long work) {
        return executor != null
                && !executor.isShutdown()
                && work >= PARALLEL_WORK_THRESHOLD
                && GemmParallelSupport.taskCountUncapped(colsB, NC) >= 3;
    }

    private static void parallelizeOverN(ExecutorService executor, int rowsA, int colsB, int colsA,
            double alpha, int offA, double[] A, int incRowA, int incColA,
            int offB, double[] B, int incRowB, int incColB,
            double beta, int offC, double[] C, int incRowC, int incColC) {

        int taskCount = GemmParallelSupport.taskCountUncapped(colsB, NC);
        if (taskCount <= 1) {
            dgemmRange(rowsA, 0, colsB, colsA, alpha, offA, A, incRowA, incColA,
                    offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC);
            return;
        }

        ArrayList<Future<?>> futures = new ArrayList<>(taskCount);
        for (int taskIndex = 0; taskIndex < taskCount; taskIndex++) {
            int jFrom = GemmParallelSupport.blockStart(taskIndex, taskCount, colsB, NC);
            int jTo   = GemmParallelSupport.blockEnd(taskIndex, taskCount, colsB, NC);
            futures.add(executor.submit(new ColumnRangeTask(
                    rowsA, jFrom, jTo, colsA, alpha, offA, A, incRowA, incColA,
                    offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC)));
        }
        GemmParallelSupport.awaitAll(futures);
    }

    //
    // Processes only columns [colsB_from, colsB_to) of C and B.
    // All panel buffers are allocated locally so this method is safe to call
    // from multiple threads with disjoint column ranges.
    //
    private static int dgemmRange(int rowsA, int colsB_from, int colsB_to, int colsA,
            double alpha, int offA, double[] A, int incRowA, int incColA,
            int offB, double[] B, int incRowB, int incColB,
            double beta, int offC, double[] C, int incRowC, int incColC) {

        int micro_kernel_calls = 0;

        final int rangeN = colsB_to - colsB_from;
        final int mb = (rowsA + MC - 1) / MC;
        final int nb = (rangeN + NC - 1) / NC;
        final int kb = (colsA + KC - 1) / KC;

        final int _mc = rowsA % MC;
        final int _nc = rangeN % NC;
        final int _kc = colsA % KC;

        //
        // Local buffers for storing panels from A, B and C - one set per thread
        //
        final double[] _A = new double[MC * KC];
        final double[] _B = new double[KC * NC];
        final double[] _C = new double[MR_Height * NR_Width];
        final double[] AB = new double[MR_Height * NR_Width];

        for (int j = 0; j < nb; ++j) {
            int nc = (j != nb - 1 || _nc == 0) ? NC : _nc;
            int globalColOffset = colsB_from + j * NC;

            for (int l = 0; l < kb; ++l) {
                int kc = (l != kb - 1 || _kc == 0) ? KC : _kc;
                double _beta = (l == 0) ? beta : 1.0;

                pack_B(kc, nc, (offB + l * KC * incRowB + globalColOffset * incColB), B, incRowB, incColB, _B);

                for (int i = 0; i < mb; ++i) {
                    int mc = (i != mb - 1 || _mc == 0) ? MC : _mc;

                    pack_A(mc, kc, (offA + i * MC * incRowA + l * KC * incColA), A, incRowA, incColA, _A);

                    micro_kernel_calls += dgemm_macro_kernel(mc, nc, kc, alpha, _beta,
                            (offC + i * MC * incRowC + globalColOffset * incColC), C, incRowC, incColC, _A, _B, AB, _C);
                }
            }
        }
        return micro_kernel_calls;
    }

    private static final class ColumnRangeTask implements Runnable {

        private final int rowsA;
        private final int jFrom;
        private final int jTo;
        private final int colsA;
        private final double alpha;
        private final int offA;
        private final double[] A;
        private final int incRowA;
        private final int incColA;
        private final int offB;
        private final double[] B;
        private final int incRowB;
        private final int incColB;
        private final double beta;
        private final int offC;
        private final double[] C;
        private final int incRowC;
        private final int incColC;

        private ColumnRangeTask(int rowsA, int jFrom, int jTo, int colsA,
                double alpha, int offA, double[] A, int incRowA, int incColA,
                int offB, double[] B, int incRowB, int incColB,
                double beta, int offC, double[] C, int incRowC, int incColC) {
            this.rowsA   = rowsA;
            this.jFrom   = jFrom;
            this.jTo     = jTo;
            this.colsA   = colsA;
            this.alpha   = alpha;
            this.offA    = offA;
            this.A       = A;
            this.incRowA = incRowA;
            this.incColA = incColA;
            this.offB    = offB;
            this.B       = B;
            this.incRowB = incRowB;
            this.incColB = incColB;
            this.beta    = beta;
            this.offC    = offC;
            this.C       = C;
            this.incRowC = incRowC;
            this.incColC = incColC;
        }

        @Override
        public void run() {
            dgemmRange(rowsA, jFrom, jTo, colsA, alpha, offA, A, incRowA, incColA,
                    offB, B, incRowB, incColB, beta, offC, C, incRowC, incColC);
        }
    }
}
