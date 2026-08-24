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
 * Represents an operation that accepts three input arguments and returns no
 * result. This is the three-arity specialization of {@link DConsumer}. Unlike
 * most other functional interfaces, {@code DTriConsumer} is expected to operate
 * via side-effects.
 *
 * @since 1.5.2
 */
public interface DTriConsumer {
    /**
     * Performs this operation on the given arguments.
     *
     * @param x1
     *            the first input argument
     * @param x2
     *            the second input argument
     * @param x3
     *            the third input argument
     */
    void accept(double x1, double x2, double x3);
}
