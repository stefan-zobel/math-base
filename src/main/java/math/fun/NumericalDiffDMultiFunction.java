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

import math.MathConsts;

/**
 * Finite difference numerical gradient calculation using the forward difference
 * approximation {@code f'(x) = (f(x + h) - f(x)) / h}.
 * <p>
 * Each individual {@code h} is scaled by the absolute magnitude of the
 * corresponding element of the vector {@code x}, <em>floored at one</em> --
 * the same step {@link NumericalDiffDVectorFunction} takes per column of a
 * Jacobian. One component of the gradient costs one evaluation of
 * {@link #apply(double[])}, so a gradient costs {@code n + 1} of them.
 * <p>
 * The floor is what keeps a component near zero from returning something that
 * is not a derivative at all. A purely relative step collapses along with the
 * component it is taken from: for {@code f(x) = sin(x0) + x1 * x1}, whose
 * derivative in {@code x0} is one to nine digits everywhere below, it comes
 * out 40 per cent wrong at {@code x0 = 1e-8} and exactly zero from
 * {@code x0 = 1e-9} downwards, because the change it induces in {@code f} is
 * below the last bit of {@code f} itself. Nor is the floor a concession where
 * the relative step works: for an {@code f} whose curvature is of order one it
 * is the more accurate of the two by up to three orders of magnitude, and from
 * {@code |x| >= 1} the two are the same number.
 * <p>
 * What it assumes is that the arguments are scaled to be of order one. Where
 * they are not -- where {@code f} varies on the scale of {@code x} itself and
 * that scale is far below one -- the floor steps clean over the feature, and
 * the way to say so is {@code diffScale}: multiplied by the magnitude the
 * arguments actually have, it buys the accuracy of the relative step back.
 * That covers arguments sharing one scale, not a vector whose components have
 * several.
 * <p>
 * The forward difference is accurate to about the square root of the machine
 * epsilon. Where the derivative can be written down, implement
 * {@link DGradient} directly instead.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Numerical_differentiation">
 *      Numerical differentiation</a>
 */
public abstract class NumericalDiffDMultiFunction implements DiffDMultiFunction {

    /**
     * The scaling factor to use for scaling of the individual {@code h}
     */
    protected final double diffScale;

    /**
     * Creates a new {@code NumericalDiffDMultiFunction} with a default scaling
     * factor.
     */
    public NumericalDiffDMultiFunction() {
        this(1.5 * Math.sqrt(MathConsts.MACH_EPS_DBL));
    }

    /**
     * Creates a new {@code NumericalDiffDMultiFunction} with the provided
     * scaling factor.
     * 
     * @param diffScale
     *            scaling factor to use for {@code h}, finite and positive. A
     *            function whose values are themselves computed to less than
     *            machine precision -- out of a simulation, a quadrature or a
     *            table -- needs a larger one, or the difference measures the
     *            error in the values rather than the derivative
     * @throws IllegalArgumentException
     *             if {@code diffScale} is not finite and positive
     */
    public NumericalDiffDMultiFunction(double diffScale) {
        if (!(diffScale > 0.0) || Double.isInfinite(diffScale)) {
            throw new IllegalArgumentException("diffScale must be finite and positive : " + diffScale);
        }
        this.diffScale = diffScale;
    }

    /**
     * {@inheritDoc}
     * <p>
     * {@code x} is written to and restored again, once per component, so this
     * method must not be entered twice on the same array at the same time.
     */
    @Override
    public final void derivativeAt(double[] x, double[] grad) {
        double fx = this.apply(x);

        for (int i = 0; i < x.length; ++i) {
            double xi = x[i];
            double hi = diffScale * Math.max(Math.abs(xi), 1.0);

            double xi_plus_hi = xi + hi;

            // account for potential round-off errors
            hi = xi_plus_hi - xi;

            x[i] = xi_plus_hi;
            // new function value for advance in variable i
            double fx_plus_hi = this.apply(x);
            // estimated gradient component for variable i
            grad[i] = (fx_plus_hi - fx) / hi;

            // restore the old value for variable i
            x[i] = xi;
        }
    }
}
