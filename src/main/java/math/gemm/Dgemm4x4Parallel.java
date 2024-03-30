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

import java.util.ArrayList;
import java.util.concurrent.Future;

import math.gemm.DgemmTasks.TaskConfig;

final class Dgemm4x4Parallel {

    //
    // Compute C <- beta*C + alpha*A*B
    //
    static int dgemm(int nb, int kb, int mb, int _nc, int _kc, int _mc, double alpha, int offA, double[] A, int incRowA,
            int incColA, int offB, double[] B, int incRowB, int incColB, double beta, int offC, double[] C, int incRowC,
            int incColC) {

        if (nb <= 1 && kb <= 1 && mb <= 1) {
            throw new AssertionError("nb=" + nb + ", kb=" + kb + ", mb=" + mb);
        }
        if (DgemmTasks.availableCores() < 2) {
            throw new AssertionError("can't parallelize: cores = " + DgemmTasks.availableCores());
        }

        DgemmTasks.TaskConfig[] cfgs = DgemmTasks.split(nb, kb, mb);
        ArrayList<Future<Integer>> results = new ArrayList<>(cfgs.length);
        for (TaskConfig cfg : cfgs) {
            Dgemm4x4Worker task = new Dgemm4x4Worker(cfg, _nc, _kc, _mc, alpha, offA, A, incRowA, incColA, offB, B,
                    incRowB, incColB, beta, offC, C, incRowC, incColC);
            Future<Integer> calls = ThreadPool.submit(task);
            results.add(calls);
        }
        int micro_kernel_calls = 0;
        for (Future<Integer> res : results) {
            try {
                micro_kernel_calls += res.get();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        return micro_kernel_calls;
    }
}
