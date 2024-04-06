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

//
// Compute A *= alpha
//
final class Gescal {

    static void gescal(int m, int n, double alpha, int A_start, double[] A, int incRowA, int incColA) {
        if (alpha == 1.0 || m <= 0 || n <= 0) {
            return;
        }

        if (incRowA == 1) {
            for (int j = 0; j < n; ++j) {
                Scal.scal(m, alpha, A, A_start + j * incColA, 1);
            }
        } else if (incColA == 1) {
            for (int i = 0; i < m; ++i) {
                Scal.scal(n, alpha, A, A_start + i * incRowA, 1);
            }
        } else {
            if (alpha != 0.0) {
                for (int j = 0; j < n; ++j) {
                    int base_A = A_start + j * incColA;
                    for (int i = 0; i < m; ++i) {
                        A[base_A + i * incRowA] *= alpha;
                    }
                }
            } else {
                for (int j = 0; j < n; ++j) {
                    int base_A = A_start + j * incColA;
                    for (int i = 0; i < m; ++i) {
                        A[base_A + i * incRowA] = 0.0;
                    }
                }
            }
        }
    }

    private Gescal() {
        throw new AssertionError();
    }
}
