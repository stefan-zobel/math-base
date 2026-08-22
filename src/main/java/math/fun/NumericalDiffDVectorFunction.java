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

import math.MathConsts;

/**
 * Finite difference numerical Jacobian calculation using the forward difference
 * approximation {@code df/dx = (f(x + h) - f(x)) / h}, the counterpart of
 * {@link NumericalDiffDMultiFunction} for a function with more than one result.
 * <p>
 * Scaling of {@code h} is taken into account by each individual {@code h} being
 * based upon the absolute magnitude of the corresponding element in the vector
 * {@code x}, <em>floored at one</em>. One column of the Jacobian costs one
 * evaluation of {@link #valueAt(double[], double[])}, so a Jacobian costs
 * {@code n + 1} of them.
 * <p>
 * The floor is the one place where this class deviates from
 * {@link NumericalDiffDMultiFunction}, and it is not cosmetic. A purely
 * relative step collapses along with the component it is taken from: for
 * {@code F(x) = (x0^2 + x1^2 - 1, x0 + x1 - 1)} at {@code x0 = 1e-4} it gets
 * the first column wrong by 40 per cent, and by {@code x0 = 1e-9} the whole
 * column comes back as exactly zero, because the change it induces in
 * {@code F} is far below the last bit of {@code F} itself. A singular Jacobian
 * out of a perfectly ordinary point is worse than a coarse one, and the floor
 * costs an order of magnitude of accuracy where the relative step was working
 * to buy that back. It assumes the arguments are scaled to be of order one,
 * the same assumption a solver's step test makes; where they are not, say so
 * through {@code diffScale}.
 * <p>
 * The forward difference is accurate to about the square root of the machine
 * epsilon. Where the derivative can be written down, implement
 * {@link DJacobian} directly instead: a solver given the exact Jacobian both
 * converges faster and reaches a smaller residual.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Numerical_differentiation">
 *      Numerical differentiation</a>
 * @since 1.5.2
 */
public abstract class NumericalDiffDVectorFunction implements DiffDVectorFunction {

    /**
     * The number {@code m} of results this function produces
     */
    protected final int m;

    /**
     * The scaling factor to use for scaling of the individual {@code h}
     */
    protected final double diffScale;

    /**
     * Creates a new {@code NumericalDiffDVectorFunction} with a default scaling
     * factor.
     *
     * @param resultCount
     *            the number {@code m} of results
     *            {@link #valueAt(double[], double[])} produces, {@code 1} or
     *            greater
     * @throws IllegalArgumentException
     *             if {@code resultCount} is smaller than {@code 1}
     */
    public NumericalDiffDVectorFunction(int resultCount) {
        this(resultCount, 1.5 * Math.sqrt(MathConsts.MACH_EPS_DBL));
    }

    /**
     * Creates a new {@code NumericalDiffDVectorFunction} with the provided
     * scaling factor.
     *
     * @param resultCount
     *            the number {@code m} of results
     *            {@link #valueAt(double[], double[])} produces, {@code 1} or
     *            greater
     * @param diffScale
     *            scaling factor to use for {@code h}, finite and positive. A
     *            function whose values are themselves computed to less than
     *            machine precision -- out of a simulation, a quadrature or a
     *            table -- needs a larger one, or the difference measures the
     *            error in the values rather than the derivative
     * @throws IllegalArgumentException
     *             if {@code resultCount} is smaller than {@code 1} or if
     *             {@code diffScale} is not finite and positive
     */
    public NumericalDiffDVectorFunction(int resultCount, double diffScale) {
        if (resultCount < 1) {
            throw new IllegalArgumentException("resultCount must be >= 1 : " + resultCount);
        }
        if (!(diffScale > 0.0) || Double.isInfinite(diffScale)) {
            throw new IllegalArgumentException("diffScale must be finite and positive : " + diffScale);
        }
        this.m = resultCount;
        this.diffScale = diffScale;
    }

    /**
     * The number {@code m} of results {@link #valueAt(double[], double[])}
     * produces, i.e. the number of rows of the Jacobian.
     *
     * @return the number of results
     */
    public final int resultCount() {
        return m;
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code x} is written to and restored again, once per component, so this
     * method must not be entered twice on the same array at the same time. The
     * two vectors of length {@code m} it needs are allocated per call rather
     * than held as fields for that same reason: shared scratch would turn a
     * concurrent second caller from unsupported into silently wrong, and the
     * allocation is negligible beside the {@code n + 1} function evaluations.
     */
    @Override
    public final void jacobianAt(double[] x, double[] jacobian) {
        int n = x.length;
        double[] fx = new double[m];
        double[] fx_plus_hj = new double[m];
        this.valueAt(x, fx);

        for (int j = 0; j < n; ++j) {
            double xj = x[j];
            double hj = diffScale * Math.max(Math.abs(xj), 1.0);

            double xj_plus_hj = xj + hj;

            // account for potential round-off errors
            hj = xj_plus_hj - xj;

            x[j] = xj_plus_hj;
            // new function values for advance in variable j
            this.valueAt(x, fx_plus_hj);
            // estimated Jacobian column for variable j, column-major
            int col = j * m;
            for (int i = 0; i < m; ++i) {
                jacobian[col + i] = (fx_plus_hj[i] - fx[i]) / hj;
            }

            // restore the old value for variable j
            x[j] = xj;
        }
    }
}
