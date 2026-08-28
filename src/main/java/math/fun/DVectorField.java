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
 * A time dependent vector field: the right hand side of a first order system of
 * ordinary differential equations,
 * <p>
 * <code>y'(t) = f(t, y(t))</code>
 * <p>
 * with {@code y} and the derivative both of length {@code n}. The dimension is
 * not part of this interface because the initial value hands it to the solver
 * already, which keeps a one line equation a lambda rather than an anonymous
 * class.
 * <p>
 * As in {@link DVectorFunction}, the derivative is written into an output array
 * rather than returned, so that a solver evaluating the field several times per
 * step over thousands of steps allocates nothing. An implementation must not
 * write into {@code y}.
 * <p>
 * A system that does not depend on {@code t} simply ignores its first argument;
 * a higher order equation becomes a first order one by carrying the derivatives
 * as further components of {@code y}, and {@link DSecondOrderField} is the case
 * where that is worth avoiding.
 *
 * @since 1.5.3
 */
public interface DVectorField {

    /**
     * The derivative of the state {@code y} at time {@code t} gets stored into
     * the output vector {@code dydt}.
     *
     * @param t
     *            the time at which the field is evaluated
     * @param y
     *            a <code>double[]</code> state vector of length {@code n} (not
     *            modified)
     * @param dydt
     *            a <code>double[]</code> output vector of length {@code n}
     *            receiving the derivative at {@code (t, y)} (modified)
     */
    void valueAt(double t, double[] y, double[] dydt);
}
