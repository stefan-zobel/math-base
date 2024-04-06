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

//
// Compute Y += alpha*X
//
final class Geaxpy {

    static void geaxpy(int m, int n, double alpha, int X_start, double[] X, int incRowX, int incColX,
            int Y_start, double[] Y, int incRowY, int incColY) {
        if (m <= 0 || n <= 0 || alpha == 0.0) {
            return;
        }

        if (incRowX == 1 && incRowY == 1) {
            //
            // X and Y are both column major
            //
            for (int j = 0; j < n; ++j) {
                Axpy.axpy(m, alpha, X, X_start + j * incColX, 1, Y, Y_start + j * incColY, 1);
            }
        } else if (incColX == 1 && incColY == 1) {
            //
            // X and Y are both row major
            //
            for (int i = 0; i < m; ++i) {
                Axpy.axpy(n, alpha, X, X_start + i * incRowX, 1, Y, Y_start + i * incRowY, 1);
            }
        } else {
            //
            // General case
            //
            for (int j = 0; j < n; ++j) {
                int base_Y = Y_start + j * incColY;
                int _incColX = j * incColX;
                for (int i = 0; i < m; ++i) {
                    Y[base_Y + i * incRowY] += alpha * X[i * incRowX + _incColX];
                }
            }
        }
    }

    private Geaxpy() {
        throw new AssertionError();
    }
}
