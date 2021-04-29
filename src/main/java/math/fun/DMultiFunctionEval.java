/*
 * Copyright 2013 Stefan Zobel
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
 * Describes the final outcome of an iterative evaluation procedure of a
 * {@link DMultiFunction}.
 */
public final class DMultiFunctionEval {

    /**
     * To which point has the evaluation settled?
     */
    public final double[] point;

    /**
     * What is the function's value at {@link #point}?
     */
    public final double value;

    /**
     * How many iterations were performed during the evaluation?
     */
    public final int iterations;

    /**
     * Has the evaluation converged?
     */
    public final boolean converged;

    /**
     * Creates a new {@code DMultiFunctionEval} instance from the given
     * arguments.
     * 
     * @param point
     *            the point at which the evaluation has settled
     * @param value
     *            the function's value at {@code point}
     * @param iterations
     *            the number of iterations that were performed during the
     *            evaluation
     * @param converged
     *            whether the evaluation has converged
     */
    public DMultiFunctionEval(double[] point, double value, int iterations, boolean converged) {
        this.point = point;
        this.value = value;
        this.iterations = iterations;
        this.converged = converged;
    }

    /**
     * Creates a new {@code DMultiFunctionEval} instance from the given
     * arguments.
     * 
     * @param point
     *            the point at which the evaluation has settled
     * @param value
     *            the function's value at {@code point}
     * @param iterations
     *            the number of iterations that were performed during the
     *            evaluation
     */
    public DMultiFunctionEval(double[] point, double value, int iterations) {
        this(point, value, iterations, true);
    }

    /**
     * Creates a new {@code DMultiFunctionEval} instance from the given
     * arguments.
     * 
     * @param point
     *            the point at which the evaluation has settled
     * @param value
     *            the function's value at {@code point}
     */
    public DMultiFunctionEval(double[] point, double value) {
        this(point, value, 0);
    }
}
