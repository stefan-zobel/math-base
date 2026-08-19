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
package math.fun;

/**
 * Interface for the Jacobian matrix of once-differentiable vector-valued
 * functions over double[] arrays, the counterpart of {@link DGradient} for a
 * function with more than one result.
 *
 * @since 1.5.2
 */
public interface DJacobian {

    /**
     * The matrix of first derivatives (a.k.a. Jacobian) of a vector-valued
     * function over a double[] array evaluated at the input location {@code x}
     * gets stored into the output array {@code jacobian}.
     * <p>
     * The matrix is <em>flat and column-major</em>, as everywhere else in this
     * library: entry {@code (i, j)}, the derivative of result {@code i} with
     * respect to argument {@code j}, sits at {@code jacobian[j * m + i]}, where
     * {@code m} is the number of results. Column {@code j} therefore occupies
     * the contiguous range {@code [j * m, (j + 1) * m)}, which is the order a
     * derivative is usually cheapest to compute in.
     *
     * @param x
     *            a <code>double[]</code> input vector of length {@code n} (not
     *            modified)
     * @param jacobian
     *            a <code>double[]</code> output array of length {@code m * n}
     *            receiving the Jacobian at location {@code x} in column-major
     *            order (modified)
     */
    void jacobianAt(double[] x, double[] jacobian);
}
