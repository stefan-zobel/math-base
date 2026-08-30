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
 * A {@link DVectorField} that can state its own first derivatives, the
 * counterpart of {@link DiffDVectorFunction} for the right hand side of a
 * differential equation.
 * <p>
 * An implicit method has to linearize the field at the point it is stepping
 * from, and there are only two ways to come by that linearization: ask for it
 * or difference it. Differencing costs {@code n + 1} evaluations of the field
 * per step and is accurate to about the square root of the machine epsilon;
 * where the derivative can be written down, saying so here is both cheaper and
 * more accurate, and on a stiff problem it is also what keeps the step size up.
 * <p>
 * <b>Both derivatives arrive in one call</b>, because a method wants them at
 * the same {@code (t, y)} and nowhere else. A field that does not depend on the
 * time explicitly fills {@code dfdt} with zeros, which is one statement it has
 * to make rather than a default it could forget: a default that assumed the
 * field autonomous would be silently wrong for one that is not, and nothing
 * downstream could detect it.
 *
 * @since 1.5.3
 */
public interface DiffDVectorField extends DVectorField {

    /**
     * The two first derivatives of the field at {@code (t, y)} get stored into
     * the output arrays.
     * <p>
     * The Jacobian is <em>flat and column-major</em>, as {@link DJacobian}
     * prescribes and as everywhere else in this library: the derivative of
     * component {@code i} of the field with respect to component {@code j} of
     * the state sits at {@code dfdy[j * n + i]}, so column {@code j} occupies
     * the contiguous range {@code [j * n, (j + 1) * n)}. An implementation must
     * write every entry, including the zeros, and must not write into
     * {@code y}.
     *
     * @param t
     *            the time at which the derivatives are taken
     * @param y
     *            a <code>double[]</code> state vector of length {@code n} (not
     *            modified)
     * @param dfdy
     *            a <code>double[]</code> output array of length {@code n * n}
     *            receiving the Jacobian in column-major order (modified)
     * @param dfdt
     *            a <code>double[]</code> output array of length {@code n}
     *            receiving the partial derivative with respect to the time,
     *            all zeros for a field that does not depend on it (modified)
     */
    void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt);
}
