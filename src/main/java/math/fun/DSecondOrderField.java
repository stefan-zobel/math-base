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
 * The right hand side of a second order system of ordinary differential
 * equations,
 * <p>
 * <code>q''(t) = f(t, q(t), q'(t))</code>
 * <p>
 * with the position {@code q}, the velocity {@code v} and the resulting
 * acceleration all of length {@code n}. As in {@link DVectorField}, the result
 * is written into an output array and the dimension is not part of the
 * interface.
 * <p>
 * Every such system can be written as a first order one of twice the length,
 * and this interface exists for the methods where doing so throws away what
 * they exploit. Runge-Kutta-Nystroem methods use the second order form
 * directly; a symplectic method needs it, because the property it is chosen for
 * lives in the splitting of position and velocity.
 * <p>
 * <b>A method may require {@code f} not to depend on {@code v}.</b> The
 * separable case <code>q'' = f(t, q)</code> is the mechanical one -- an
 * acceleration that is a force over a mass and no damping -- and it is what a
 * symplectic method assumes. Such a method ignores the {@code v} it is handed;
 * passing it a velocity dependent field does not fail, it silently stops being
 * symplectic, and the class documents which of the two it is.
 *
 * @since 1.5.3
 */
public interface DSecondOrderField {

    /**
     * The acceleration at time {@code t}, position {@code q} and velocity
     * {@code v} gets stored into the output vector {@code acceleration}.
     *
     * @param t
     *            the time at which the field is evaluated
     * @param q
     *            a <code>double[]</code> position vector of length {@code n}
     *            (not modified)
     * @param v
     *            a <code>double[]</code> velocity vector of length {@code n}
     *            (not modified), which a method for the separable case ignores
     * @param acceleration
     *            a <code>double[]</code> output vector of length {@code n}
     *            receiving the second derivative (modified)
     */
    void valueAt(double t, double[] q, double[] v, double[] acceleration);
}
