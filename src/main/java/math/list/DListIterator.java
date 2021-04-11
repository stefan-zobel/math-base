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
package math.list;

import java.util.NoSuchElementException;

import math.fun.DForEachIterator;

/**
 * An iterator for lists that allows the programmer to traverse the list in
 * either direction, modify the list during iteration, and obtain the iterator's
 * current position in the list. A {@code DListIterator} has no current element;
 * its <i>cursor position</i> always lies between the element that would be
 * returned by a call to {@code previous()} and the element that would be
 * returned by a call to {@code next()}. An iterator for a list of length
 * {@code n} has {@code n+1} possible cursor positions, as illustrated by the
 * carets ({@code ^}) below:
 * 
 * <pre>
 *                      Element(0)   Element(1)   Element(2)   ... Element(n-1)
 * cursor positions:  ^            ^            ^            ^                  ^
 * </pre>
 * 
 * Note that the {@link #remove} and {@link #set(double)} methods are <i>not</i>
 * defined in terms of the cursor position; they are defined to operate on the
 * last element returned by a call to {@link #next} or {@link #previous()}.
 */
public interface DListIterator extends DForEachIterator {

    /**
     * Returns {@code true} if this list iterator has more elements when
     * traversing the list in the reverse direction.
     * 
     * @return {@code true} if the list iterator has more elements when
     *         traversing the list in the reverse direction
     */
    boolean hasPrevious();

    /**
     * Returns the previous element in the list and moves the cursor position
     * backwards. This method may be called repeatedly to iterate through the
     * list backwards, or intermixed with calls to {@link #next} to go back and
     * forth. (Note that alternating calls to {@code next} and {@code previous}
     * will return the same element repeatedly.)
     * 
     * @return the previous element in the list
     * @throws NoSuchElementException
     *             if the iteration has no previous element
     */
    double previous();

    /**
     * Returns the index of the element that would be returned by a subsequent
     * call to {@link #next}. (Returns list size if the list iterator is at the
     * end of the list.)
     *
     * @return the index of the element that would be returned by a subsequent
     *         call to {@code next}, or list size if the list iterator is at the
     *         end of the list
     */
    int nextIndex();

    /**
     * Returns the index of the element that would be returned by a subsequent
     * call to {@link #previous}. (Returns -1 if the list iterator is at the
     * beginning of the list.)
     *
     * @return the index of the element that would be returned by a subsequent
     *         call to {@code previous}, or -1 if the list iterator is at the
     *         beginning of the list
     */
    int previousIndex();

    /**
     * Replaces the last element returned by {@link #next} or {@link #previous}
     * with the specified element. This call can be made only if neither
     * {@link #remove} nor {@link #add} have been called after the last call to
     * {@code next} or {@code previous}.
     *
     * @param e
     *            the element with which to replace the last element returned by
     *            {@code next} or {@code previous}
     * @throws IllegalStateException
     *             if neither {@code next} nor {@code previous} have been
     *             called, or {@code remove} or {@code add} have been called
     *             after the last call to {@code next} or {@code previous}
     */
    void set(double e);

    /**
     * Inserts the specified element into the list. The element is inserted
     * immediately before the element that would be returned by {@link #next},
     * if any, and after the element that would be returned by
     * {@link #previous}, if any. (If the list contains no elements, the new
     * element becomes the sole element on the list.) The new element is
     * inserted before the implicit cursor: a subsequent call to {@code next}
     * would be unaffected, and a subsequent call to {@code previous} would
     * return the new element. (This call increases by one the value that would
     * be returned by a call to {@code nextIndex} or {@code previousIndex}.)
     *
     * @param e
     *            the element to insert
     */
    void add(double e);
}
