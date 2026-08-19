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
 * A once-differentiable vector-valued multivariate function over double[]
 * arrays, the counterpart of {@link DiffDMultiFunction} for a function with
 * more than one result.
 *
 * @since 1.5.2
 */
public interface DiffDVectorFunction extends DVectorFunction, DJacobian {
}
