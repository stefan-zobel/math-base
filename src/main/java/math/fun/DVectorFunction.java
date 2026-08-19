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
 * A vector-valued multivariate function over double[] arrays, taking {@code n}
 * arguments to {@code m} results.
 * <p>
 * This is the shape a nonlinear least squares problem has: {@code m} residuals
 * over {@code n} parameters, with {@code m >= n}. The results are written into
 * an output array rather than returned, so that a solver evaluating the
 * function thousands of times allocates nothing.
 *
 * @since 1.5.2
 */
public interface DVectorFunction {

    /**
     * The value of this function at the input location {@code x} gets stored
     * into the output vector {@code values}.
     *
     * @param x
     *            a <code>double[]</code> input vector of length {@code n} (not
     *            modified)
     * @param values
     *            a <code>double[]</code> output vector of length {@code m}
     *            receiving the value at location {@code x} (modified)
     */
    void valueAt(double[] x, double[] values);
}
