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
package math.trmm;

/**
 * Compute X *= alpha.
 */
final class Scal {

    static void scal(int n, double alpha, double[] x, int xOff, int incX) {
        if (alpha != 1.0 && alpha != 0.0) {
            for (int i = 0; i < n; ++i) {
                x[xOff + i * incX] *= alpha;
            }
        } else if (alpha == 0.0) {
            for (int i = 0; i < n; ++i) {
                x[xOff + i * incX] *= 0.0;
            }
        }
    }

    private Scal() {
        throw new AssertionError();
    }
}
