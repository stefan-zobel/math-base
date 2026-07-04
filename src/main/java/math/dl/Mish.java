/*
 * Copyright 2026 Stefan Zobel
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

import math.imsl.Trig;

/**
 * The Mish activation function in deep learning. Mish is a self-regularized
 * non-monotonic activation function defined as: f(x) = x * tanh(softplus(x)) =
 * x * tanh(ln(1 + e^x))
 */
public final class Mish {

    /**
     * Computes the value of the Mish function at {@code x}.
     * 
     * @param x input value
     * @return Mish value at x
     */
    public static double mish(double x) {
        // Fast cutoff for highly negative inputs where Mish approaches 0
        if (x < -8.0) {
            return -0.0;
        }
        // Fast cutoff for highly positive inputs where Mish approaches x
        if (x > 16.0) {
            return x;
        }

        // softplus(x) = ln(1 + e^x)
        double softplus = Math.log(1.0 + Math.exp(x));
        return x * Trig.tanh(softplus);
    }

    /**
     * Computes the value of the Mish function at {@code x}.
     * 
     * @param x input value
     * @return Mish value at x
     */
    public static float mishF(float x) {
        return (float) mish(x);
    }

    /**
     * Computes the derivative of the Mish function with respect to {@code x}.
     * 
     * @param x input value
     * @return Mish derivative at x
     */
    public static double dmish_dx(double x) {
        // Derivative saturates to 0.0 for highly negative inputs
        if (x < -8.0) {
            return 0.0;
        }
        // Derivative saturates to 1.0 for highly positive inputs
        if (x > 16.0) {
            return 1.0;
        }

        double exp = Math.exp(x);
        double softplus = Math.log(1.0 + exp);
        double tanh = Trig.tanh(softplus);

        // Pre-calculating components to optimize performance
        // sigmoid(x) = exp(x) / (1 + exp(x))
        double sigmoid = exp / (1.0 + exp);
        double sech2 = 1.0 - (tanh * tanh);

        // Mathematical derivative: dMish/dx = tanh(softplus(x)) + x * sigmoid(x) *
        // sech^2(softplus(x))
        return tanh + x * sigmoid * sech2;
    }

    /**
     * Computes the derivative of the Mish function with respect to {@code x}.
     * 
     * @param x input value
     * @return Mish derivative at x
     */
    public static float dmishF_dx(float x) {
        return (float) dmish_dx(x);
    }

    private Mish() {
        throw new AssertionError();
    }
}
