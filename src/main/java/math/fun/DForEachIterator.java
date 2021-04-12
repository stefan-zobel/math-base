/*
 * Copyright 2021 Stefan Zobel
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

import java.util.Objects;
import java.util.function.DoubleConsumer;

/**
 * An extension of {@link DIterator} adding the two (default) methods
 * {@link #remove()} and {@link #forEachRemaining(DoubleConsumer)}.
 */
public interface DForEachIterator extends DIterator {

    /**
     * Removes from the underlying collection the last element returned by this
     * iterator. This method can be called only once per call to {@link #next}.
     * The behavior of an iterator is unspecified if the underlying collection
     * is modified while the iteration is in progress in any way other than by
     * calling this method. The default implementation throws an instance of
     * {@link UnsupportedOperationException} and performs no other action.
     *
     * @throws UnsupportedOperationException
     *             if the {@code remove} operation is not supported by this
     *             iterator
     *
     * @throws IllegalStateException
     *             if the {@code next} method has not yet been called, or the
     *             {@code remove} method has already been called after the last
     *             call to the {@code next} method
     */
    default void remove() {
        throw new UnsupportedOperationException("remove");
    }

    /**
     * Performs the given action for each remaining element until all elements
     * have been processed or the action throws an exception. Actions are
     * performed in the order of iteration. Exceptions thrown by the action are
     * relayed to the caller.
     *
     * <p>
     * The default implementation behaves as if:
     * 
     * <pre>
     * {@code
     *     while (hasNext())
     *         action.accept(next());
     * }
     * </pre>
     *
     * @param action
     *            The action to be performed for each element
     * @throws NullPointerException
     *             if the specified action is null
     */
    default void forEachRemaining(DoubleConsumer action) {
        Objects.requireNonNull(action);
        while (hasNext()) {
            action.accept(next());
        }
    }
}
