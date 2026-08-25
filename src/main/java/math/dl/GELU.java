/*
 * Copyright 2023, 2024 Stefan Zobel
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
package math.dl;

import math.MathConsts;
import math.imsl.Trig;

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
        // measured: below -30.0 the result underflows to -0.0 anyway, and
        // from 7.09 upwards it is exactly x. The cutoffs used to sit at
        // +/-4.861, which is where the tangent gave up rather than where a
        // double does, so gelu(-4.861) answered -0.0 for a true -5.73e-07.
        // They also keep x * x * x, which overflows past 5.6e102, out of the
        // way
        if (x < -30.0) {
            return -0.0;
        }
        if (x > 7.09) {
            return x;
        }
        return x * halfOnePlusTanh(MathConsts.SQRT_TWO_BY_PI * (x + 0.044715 * (x * x * x)));
    }

    /**
     * Computes the value of the GELU function at {@code x}.
     * 
     * @param x
     *            input value
     * @return GELU value at x
     */
    public static float geluF(float x) {
        return (float) gelu(x);
    }

    /**
     * Computes the derivative of the GELU function with respect to {@code x}.
     * 
     * @param x
     *            input value
     * @return GELU derivative at x
     */
    public static double dgelu_dx(double x) {
        // measured: below -21.5 the result underflows to 0.0, and from 7.45
        // upwards it is exactly 1.0. The upper cutoff was already right; the
        // lower one used to mirror it and cut the tail off 14 units early.
        // They also keep x * x, which overflows past 1.3e154, out of the way
        if (x < -21.5) {
            return -0.0;
        }
        if (x > 7.45) {
            return 1.0;
        }
        double x2 = x * x;
        double x3 = x2 * x;
        double y = MathConsts.SQRT_TWO_BY_PI * (x + 0.044715 * x3);
        // sech2 goes through cosh and does not cancel, so only the first
        // term needs the careful form
        return halfOnePlusTanh(y) + 0.5 * x * Trig.sech2(y) * MathConsts.SQRT_TWO_BY_PI * (1.0 + 0.134145 * x2);
    }

    /**
     * Computes the derivative of the GELU function with respect to {@code x}.
     * 
     * @param x
     *            input value
     * @return GELU derivative at x
     */
    public static float dgeluF_dx(float x) {
        return (float) dgelu_dx(x);
    }

    // 0.5 * (1 + tanh(y)) written so that it does not cancel. As tanh(y)
    // approaches -1 the sum loses every digit it has -- at y = -19 the
    // straightforward form answers 5.551e-17 where the truth is 2.784e-17, and
    // from y = -38 on it answers zero. The logistic identity
    // 0.5 * (1 + tanh(y)) = 1 / (1 + e^(-2y)) has no subtraction in it. It is
    // only used below y = -1, where the sum starts to lose bits; above that the
    // two agree to the last bit and the tangent is the cheaper of the two,
    // because it needs no exp for |y| <= 1.
    private static double halfOnePlusTanh(double y) {
        if (y < -1.0) {
            return 1.0 / (1.0 + Math.exp(-2.0 * y));
        }
        return 0.5 * (1.0 + Trig.tanh(y));
    }

    private GELU() {
        throw new AssertionError();
    }
}
