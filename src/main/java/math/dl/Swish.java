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

/**
 * The Swish activation function in deep learning (with beta = 1).
 */
public final class Swish {

    /**
     * Computes the value of the Swish function at {@code x}.
     * 
     * @param x input value
     * @return Swish value at x
     */
    public static double swish(double x) {
        // Asymptotic limit for extremely negative inputs
        if (x < -17.0) {
            return -0.0;
        }
        // For very large positive inputs, Swish approaches exactly x
        if (x > 17.0) {
            return x;
        }
        return x / (1.0 + Math.exp(-x));
    }

    /**
     * Computes the value of the Swish function at {@code x}.
     * 
     * @param x input value
     * @return Swish value at x
     */
    public static float swishF(float x) {
        return (float) swish(x);
    }

    /**
     * Computes the derivative of the Swish function with respect to {@code x}.
     * 
     * @param x input value
     * @return Swish derivative at x
     */
    public static double dswish_dx(double x) {
        // Fast path for saturated derivative
        if (x < -17.0) {
            return 0.0;
        }
        if (x > 17.0) {
            return 1.0;
        }

        double sigmoid = 1.0 / (1.0 + Math.exp(-x));
        // Mathematical derivative: f'(x) = f(x) + sigmoid(x) * (1 - f(x))
        // Alternatively: f'(x) = swish(x) + sigmoid(x) * (1.0 - swish(x))
        return x * sigmoid + sigmoid * (1.0 - x * sigmoid);
    }

    /**
     * Computes the derivative of the Swish function with respect to {@code x}.
     * 
     * @param x input value
     * @return Swish derivative at x
     */
    public static float dswishF_dx(float x) {
        return (float) dswish_dx(x);
    }

    private Swish() {
        throw new AssertionError();
    }
}
