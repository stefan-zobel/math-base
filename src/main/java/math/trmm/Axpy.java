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
 * Compute Y += alpha*X.
 */
final class Axpy {

    static void axpy(int n, double alpha, double[] x, int xOff, int incX, double[] y, int yOff, int incY) {
        if (n <= 0 || alpha == 0.0) {
            return;
        }

        if (incX == 1 && incY == 1) {
            // code for both increments equal to 1

            for (int i = 0; i < n; ++i) {
                y[i + yOff] += alpha * x[i + xOff];
            }

        } else {
            // code for unequal increments or
            // equal increments greater than 1
            if (incX <= 0 || incY <= 0) {
                throw new IllegalArgumentException("illegal non-positive incrememts: incX=" + incX + ", incY=" + incY);
            }
            for (int i = 0; i < n; ++i) {
                y[i * incY + yOff] += alpha * x[i * incX + xOff];
            }
        }
    }

    private Axpy() {
        throw new AssertionError();
    }
}
