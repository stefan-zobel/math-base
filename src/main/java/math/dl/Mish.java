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
        // measured by bisection: the result reaches -0.0 at -745.133219,
        // where Math.exp(x) underflows, and is exactly x from
        // 19.061547465398494 upwards -- which is 27.5 * ln(2), the point at
        // which the tangent saturates, since softplus(x) approaches x. The
        // cutoffs used to sit at -8 and 16, where mish(-8.0000001) answered
        // -0.0 for a true -2.683e-03. The lower one is also what keeps
        // -Infinity * 0.0, which is NaN rather than zero, out of the result
        if (x < -745.13) {
            return -0.0;
        }
        if (x > 19.061547465398494) {
            return x;
        }

        double softplus = softplusOf(x, Math.exp(x));
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
        // measured by bisection: 0.0 from -745.133219 down, exactly 1.0 from
        // 20.573551 up. The cutoffs used to be -8 and 16, and the step at -8
        // was 2.347e-03
        if (x < -745.13) {
            return 0.0;
        }
        if (x > 20.573551) {
            return 1.0;
        }

        double exp = Math.exp(x);
        double softplus = softplusOf(x, exp);
        double tanh = Trig.tanh(softplus);

        // Pre-calculating components to optimize performance
        // sigmoid(x) = exp(x) / (1 + exp(x))
        double sigmoid = exp / (1.0 + exp);
        // 1 - tanh^2 erodes as tanh approaches 1 and collapses to zero at
        // softplus = 18.715, which the upper cutoff above now passes. Kept
        // anyway, measured: over the whole range its worst absolute error is
        // 8.438e-15. Trig.sech2 would give 4.441e-16 but costs a third exp
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

    // softplus(x) = ln(1 + e^x). The straightforward form loses the bits that
    // matter once e^x is small: 1 + e^x rounds away its low end, and the
    // relative error of the result is about eps / (2 * e^x) -- 1.020e-03 at
    // x = -30 and 4.271e-02 at x = -36, so with the plain form the widened
    // cutoff above would answer 3.333e-01 instead of 4.939e-16. Math.log1p has
    // no such loss but is the more elaborate routine: using it everywhere
    // costs 20 % on N(0,1) inputs, where the plain form is exact anyway.
    // Measured with the branch at -4: 1 %, and a worst relative error of
    // 6.266e-15 over the whole range
    private static double softplusOf(double x, double exp) {
        if (x < -4.0) {
            return Math.log1p(exp);
        }
        return Math.log(1.0 + exp);
    }

    private Mish() {
        throw new AssertionError();
    }
}
