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
package math.gemm;

import java.util.concurrent.Callable;

final class Dgemm4x4Worker implements Callable<Integer> {

    // controls for loop-execution
    private final int m_nblo;
    private final int m_nbhi;
    // private final int m_kblo;
    private final int m_kbhi;
    private final int m_mblo;
    private final int m_mbhi;
    private final boolean m_nbfull;
    // private final boolean m_kbfull;
    private final boolean m_mbfull;

    @SuppressWarnings("unused")
    private final String m_loop;
    @SuppressWarnings("unused")
    private final int m_tasks;

    // the whole rest
    private final int m_nc;
    private final int m_kc;
    private final int m_mc;
    private final double m_alpha;
    private final int m_offA;
    private final double[] m_A;
    private final int m_incRowA;
    private final int m_incColA;
    private final int m_offB;
    private final double[] m_B;
    private final int m_incRowB;
    private final int m_incColB;
    private final double m_beta;
    private final int m_offC;
    private final double[] m_C;
    private final int m_incRowC;
    private final int m_incColC;

    Dgemm4x4Worker(DgemmTasks.TaskConfig cfg, int _nc, int _kc, int _mc, double alpha, int offA, double[] A,
            int incRowA, int incColA, int offB, double[] B, int incRowB, int incColB, double beta, int offC, double[] C,
            int incRowC, int incColC) {

        // unpack loop controls
        m_nblo = cfg.nb_lo;
        m_nbhi = cfg.nb_hi;
        m_nbfull = !cfg.nb_hi_is_last;
        // m_kblo = cfg.kb_lo;
        m_kbhi = cfg.kb_hi;
        // m_kbfull = !cfg.kb_hi_is_last;
        m_mblo = cfg.mb_lo;
        m_mbhi = cfg.mb_hi;
        m_mbfull = !cfg.mb_hi_is_last;

        m_loop = cfg.loop;
        m_tasks = cfg.tasks;

        // set the rest
        m_nc = _nc;
        m_kc = _kc;
        m_mc = _mc;
        m_alpha = alpha;
        m_offA = offA;
        m_A = A;
        m_incRowA = incRowA;
        m_incColA = incColA;
        m_offB = offB;
        m_B = B;
        m_incRowB = incRowB;
        m_incColB = incColB;
        m_beta = beta;
        m_offC = offC;
        m_C = C;
        m_incRowC = incRowC;
        m_incColC = incColC;
    }

    @Override
    public Integer call() throws Exception {
        //
        // hoist local variables
        //
        final int _nc = m_nc;
        final int _kc = m_kc;
        final int _mc = m_mc;
        final double alpha = m_alpha;
        final int offA = m_offA;
        final double[] A = m_A;
        final int incRowA = m_incRowA;
        final int incColA = m_incColA;
        final int offB = m_offB;
        final double[] B = m_B;
        final int incRowB = m_incRowB;
        final int incColB = m_incColB;
        final double beta = m_beta;
        final int offC = m_offC;
        final double[] C = m_C;
        final int incRowC = m_incRowC;
        final int incColC = m_incColC;

        //
        // Local buffers for storing panels from A, B and C
        //
        final double[] _A = new double[Dgemm4x4.MC * Dgemm4x4.KC];
        final double[] _B = new double[Dgemm4x4.KC * Dgemm4x4.NC];
        final double[] _C = new double[Dgemm4x4.MR_Height * Dgemm4x4.NR_Width];
        final double[] AB = new double[Dgemm4x4.MR_Height * Dgemm4x4.NR_Width];

        // counter for microkernel calls
        int micro_kernel_calls = 0;

        // hoist loop bounds
        final int nblo = m_nblo;
        final int nbhi = m_nbhi;
        final boolean nbfull = m_nbfull;

        // final int kblo = m_kblo;
        final int kbhi = m_kbhi;
        // final boolean kbfull = m_kbfull;

        final int mblo = m_mblo;
        final int mbhi = m_mbhi;
        final boolean mbfull = m_mbfull;

//        System.out.println(
//                "computation on " + Thread.currentThread().getName() + " - loop : " + m_loop + " , tasks: " + m_tasks);

        for (int j = nblo; j < nbhi; ++j) {
            int nc = (_nc == 0 || nbfull || (!nbfull && j != nbhi - 1)) ? Dgemm4x4.NC : _nc;

            // kb loop can't be parallelized!
            for (int l = 0; l < kbhi; ++l) {
                int kc = (l != kbhi - 1 || _kc == 0) ? Dgemm4x4.KC : _kc;
                double _beta = (l == 0) ? beta : 1.0;

                Dgemm4x4.pack_B(kc, nc, (offB + l * Dgemm4x4.KC * incRowB + j * Dgemm4x4.NC * incColB), B, incRowB,
                        incColB, _B);

                for (int i = mblo; i < mbhi; ++i) {
                    int mc = (_mc == 0 || mbfull || (!mbfull && i != mbhi - 1)) ? Dgemm4x4.MC : _mc;

                    Dgemm4x4.pack_A(mc, kc, (offA + i * Dgemm4x4.MC * incRowA + l * Dgemm4x4.KC * incColA), A, incRowA,
                            incColA, _A);

                    micro_kernel_calls += Dgemm4x4.dgemm_macro_kernel(mc, nc, kc, alpha, _beta,
                            (offC + i * Dgemm4x4.MC * incRowC + j * Dgemm4x4.NC * incColC), C, incRowC, incColC, _A, _B,
                            AB, _C);
                }
            }
        }
        return micro_kernel_calls;
    }
}
