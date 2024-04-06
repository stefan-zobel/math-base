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

final class Trupack {

    static void trupack(int mc, boolean unit, int U_start, double[] U, int incRowU, int incColU,
            double[] buffer) {
        final int mp = mc / BlockSizes.MR;
        final int mr_ = (mc & (BlockSizes.MR - 1));

        int buffer_start = 0;

        for (int i = 0; i < mp; ++i) {
            trupack_MRxk(mc - i * BlockSizes.MR, unit, U_start, U, incRowU, incColU, buffer, buffer_start);
            buffer_start += (mc - i * BlockSizes.MR) * BlockSizes.MR;
            U_start += BlockSizes.MR * (incRowU + incColU);
        }
        if (mr_ > 0) {
            for (int j = 0; j < mr_; ++j) {
                for (int i = 0; i < j; ++i) {
                    buffer[buffer_start + i] = U[U_start + i * incRowU];
                }
                buffer[buffer_start + j] = (unit) ? 1.0 : U[U_start + j * incRowU];
                for (int i = j + 1; i < BlockSizes.MR; ++i) {
                    buffer[buffer_start + i] = 0.0;
                }
                buffer_start += BlockSizes.MR;
                U_start += incColU;
            }
        }
    }

    private static void trupack_MRxk(int k, boolean unit, int U_start, double[] U, int incRowU, int incColU,
            double[] buffer, int buffer_start) {

        for (int j = 0; j < BlockSizes.MR; ++j) {
            for (int i = 0; i < j; ++i) {
                buffer[buffer_start + i] = U[U_start + i * incRowU];
            }
            buffer[buffer_start + j] = (unit) ? 1.0 : U[U_start + j * incRowU];
            for (int i = j + 1; i < BlockSizes.MR; ++i) {
                buffer[buffer_start + i] = 0.0;
            }
            buffer_start += BlockSizes.MR;
            U_start += incColU;
        }
        for (int j = 0; j < k - BlockSizes.MR; ++j) {
            for (int i = 0; i < BlockSizes.MR; ++i) {
                buffer[buffer_start + i] = U[U_start + i * incRowU];
            }
            buffer_start += BlockSizes.MR;
            U_start += incColU;
        }
    }

    private Trupack() {
        throw new AssertionError();
    }
}
