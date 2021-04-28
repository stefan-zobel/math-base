/*
 * Copyright 2018 Stefan Zobel
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
 * Represents an operation that accepts two input arguments and returns no
 * result. This is the two-arity specialization of {@link DConsumer}. Unlike
 * most other functional interfaces, {@code DBiConsumer} is expected to operate
 * via side-effects.
 */
public interface DBiConsumer {
    /**
     * Performs this operation on the given arguments.
     *
     * @param x
     *            the first input argument
     * @param y
     *            the second input argument
     */
    void accept(double x, double y);
}
