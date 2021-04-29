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
 * A device for traversing elements of a source, either sequentially in bulk or
 * individually.
 */
public interface DForEachBi {

    /**
     * Performs the given action for each remaining element sequentially until
     * all elements have been processed or the action throws an exception.
     * Exceptions thrown by the action are relayed to the caller.
     *
     * @param action
     *            The action
     * @throws NullPointerException
     *             if the specified action is null
     */
    void forEachRemaining(DBiConsumer action);

    /**
     * If a remaining element exists, performs the given action on it, returning
     * {@code true}; else returns {@code false}. Exceptions thrown by the action
     * are relayed to the caller.
     *
     * @param action
     *            The action
     * @return {@code false} if no remaining elements existed upon entry to this
     *         method, else {@code true}.
     * @throws NullPointerException
     *             if the specified action is null
     */
    boolean tryAdvance(DBiConsumer action);
}
