/*
 * Copyright 2023 Stefan Zobel
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
package math.imsl;

import math.MathConsts;

/**
 * The GELU activation function in deep learning.
 */
public final class GELU {

    /**
     * Computes the value of the GELU function at {@code x}.
     * 
     * @param x
     *            input value
     * @return GELU value at x
     */
    public static double gelu(double x) {
        if (x < -4.860999999999954) {
            return -0.0;
        }
        if (x > 4.861000000000067) {
            return x;
        }
        return 0.5 * x * (1.0 + Trig.tanh(MathConsts.SQRT_TWO_BY_PI * (x + 0.044715 * (x * x * x))));
    }

    /**
     * Computes the derivative of the GELU function with respect to {@code x}.
     * 
     * @param x
     *            input value
     * @return GELU derivative at x
     */
    public static double dgelu_dx(double x) {
        if (x < -7.446000000000483) {
            return -0.0;
        }
        if (x > 7.446000000000483) {
            return 1.0;
        }
        double x2 = x * x;
        double x3 = x2 * x;
        double y = MathConsts.SQRT_TWO_BY_PI * (x + 0.044715 * x3);
        return 0.5 * (1.0 + Trig.tanh(y)) + 0.5 * x * Trig.sech2(y) * MathConsts.SQRT_TWO_BY_PI * (1.0 + 0.134145 * x2);
    }

    private GELU() {
        throw new AssertionError();
    }
}
