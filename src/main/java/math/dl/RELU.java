/*
 * Copyright 2024, 2026 Stefan Zobel
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
 * The RELU activation function in deep learning.
 */
public final class RELU {

    /**
     * Computes the value of the RELU function at {@code x}.
     * 
     * @param x
     *            input value
     * @return RELU value at x
     */
    public static double relu(double x) {
        // Math.max rather than a comparison: x > 0.0 is false for a NaN, so
        // the plain form answered 0.0 for one and hid whatever produced it.
        // Math.max propagates it, and still answers 0.0 for -0.0
        return Math.max(x, 0.0);
    }

    /**
     * Computes the value of the RELU function at {@code x}.
     * 
     * @param x
     *            input value
     * @return RELU value at x
     */
    public static float reluF(float x) {
        return Math.max(x, 0.0f);
    }

    /**
     * Computes the derivative of the RELU function with respect to {@code x}.
     * 
     * @param x
     *            input value
     * @return RELU derivative at x
     */
    public static double drelu_dx(double x) {
        if (Double.isNaN(x)) {
            return Double.NaN;
        }
        return (x > 0.0) ? 1.0 : 0.0;
    }

    /**
     * Computes the derivative of the RELU function with respect to {@code x}.
     * 
     * @param x
     *            input value
     * @return RELU derivative at x
     */
    public static float dreluF_dx(float x) {
        if (Float.isNaN(x)) {
            return Float.NaN;
        }
        return (x > 0.0f) ? 1.0f : 0.0f;
    }

    private RELU() {
        throw new AssertionError();
    }
}
