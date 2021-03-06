/*
 * Copyright 2020 Stefan Zobel
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
 * Represents a function that accepts three arguments and produces a result. This
 * is the three-arity specialization of {@link DFunction}.
 */
public interface DTriFunction {
    /**
     * Applies this function to the given arguments.
     *
     * @param x1
     *            the first function argument
     * @param x2
     *            the second function argument
     * @param x3
     *            the third function argument            
     * @return the function result
     */
    double apply(double x1, double x2, double x3);
}
