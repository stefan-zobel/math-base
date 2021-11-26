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

import java.util.Arrays;
import java.util.Objects;
import java.util.Spliterator;
import java.util.function.DoubleConsumer;
import java.util.stream.DoubleStream;
import java.util.stream.StreamSupport;

import math.fun.DForEachIterator;

/**
 * A resizable primitive double[] array with vector operations that behaves
 * similarly to a {@code java.util.List}.
 */
public interface DoubleList {

    /**
     * Returns a sequential {@code Stream} with this list as its source.
     * 
     * @return a sequential {@code Stream} over the elements in this list
     */
    default DoubleStream stream() {
        return StreamSupport.doubleStream(spliterator(), false);
    }

    /**
     * Returns a parallel {@code Stream} with this list as its source.
     * 
     * @return a parallel {@code Stream} over the elements in this list
     */
    default DoubleStream parallelStream() {
        return StreamSupport.doubleStream(spliterator(), true);
    }

    /**
     * Performs the given action for each element of this list in the order of
     * iteration until all elements have been processed or the action throws an
     * exception. Exceptions thrown by the action are relayed to the caller.
     * 
     * @param action
     *            The action to be performed for each element
     * @throws NullPointerException
     *             if the specified action is null
     */
    default void forEach(DoubleConsumer action) {
        Objects.requireNonNull(action);
        for (int i = 0; i < this.size(); ++i) {
            action.accept(this.get(i));
        }
    }

    /**
     * Sorts this list in ascending order.
     */
    default void sort() {
        double[] a = this.toArray();
        Arrays.sort(a);
        for (int i = 0; i < this.size(); ++i) {
            this.set(i, a[i]);
        }
    }

    /**
     * Constructs a list containing the passed values in the order of their
     * appearance.
     * 
     * @param values
     *            the values which are to be placed into the list
     * @return a new {@code DoubleList} constructed from the passed values
     */
    static DoubleList of(double... values) {
        return new DoubleArrayList(values, false);
    }

    /**
     * Tests for approximate equality (or "closeness") of the two DoubleLists
     * {@code a} and {@code b} where {@code a} and {@code b} must have the same
     * length and each element {@code xa} of {@code a} is tested for approximate
     * equality against the corresponding element {@code xb} of {@code b} as
     * defined in {@link #approxEqual(DoubleList, DoubleList, double, double)}.
     * This method is equivalent to a call of
     * {@linkplain #approxEqual(DoubleList, DoubleList, double, double)} with an
     * {@code relTol} argument of value {@code 1.0e-8} and an {@code absTol}
     * argument equal to {@code 0.0}:
     * 
     * <pre>
     * {@code approxEqual(a, b, 1.0e-8, 0.0)}
     * </pre>
     * 
     * @param a
     *            one of the two DoubleLists to test for approximate equality
     *            (it doesn't matter which one since the test is symmetric)
     * @param b
     *            the other one of the two DoubleLists to test for approximate
     *            equality (it doesn't matter which one since the test is
     *            symmetric)
     * @return {@code true} if {@code a} and {@code b} are approximately equal
     *         according to the criterion defined in
     *         {@link #approxEqual(DoubleList, DoubleList, double, double)}
     */
    static boolean approxEqual(DoubleList a, DoubleList b) {
        return approxEqual(a, b, 1.0e-8, 0.0);
    }

    /**
     * Tests for approximate equality (or "closeness") of the two DoubleLists
     * {@code a} and {@code b} where {@code a} and {@code b} must have the same
     * length and each element {@code xa} of {@code a} is tested for approximate
     * equality against the corresponding element {@code xb} of {@code b}.
     * <p>
     * If, for all pairs {@code (xa, xb)},
     * 
     * <pre>
     * {@code abs(xa - xb) <= max( relTol * max(abs(xa), abs(xb)), absTol )}
     * </pre>
     * 
     * the lists {@code a} and {@code b} are considered approximately equal,
     * otherwise they are not. This test is symmetric, so interchanging
     * {@code a} and {@code b} doesn't change the result.
     * <p>
     * <b>Implementation Note:</b><br>
     * The definition of approximate equality used here is the one employed in
     * Python's {@code math.isclose()} function defined in <a
     * href=https://www.python.org/dev/peps/pep-0485/>PEP 485 - A Function for
     * testing approximate equality</a>. This document gives a nice discussion
     * of the rationale for this approach, how to use it, and the alternatives
     * they had considered.
     * 
     * @param a
     *            one of the two DoubleLists to test for approximate equality
     *            (it doesn't matter which one since the test is symmetric)
     * @param b
     *            the other one of the two DoubleLists to test for approximate
     *            equality (it doesn't matter which one since the test is
     *            symmetric)
     * @param relTol
     *            relative tolerance, must be {@code >= 0.0}
     * @param absTol
     *            absolute tolerance, must be {@code >= 0.0}
     * @return {@code true} if {@code a} and {@code b} are approximately equal
     *         according to the criterion defined above, otherwise {@code false}
     * @throws IllegalArgumentException
     *             if {@code relTol < 0.0} or {@code absTol < 0.0}
     */
    static boolean approxEqual(DoubleList a, DoubleList b, double relTol, double absTol) {
        if (!DoubleArrayList.checkApproxEqualArgs(a, b, relTol, absTol)) {
            return false;
        }
        if (a == b) {
            return true;
        }
        double[] _a = a.getArrayUnsafe();
        double[] _b = b.getArrayUnsafe();
        int aoff = a.offset();
        int boff = b.offset();
        // a.size() == b.size() guaranteed here
        int len = a.size();
        for (int i = 0; i < len; ++i) {
            double xa = _a[aoff + i];
            double xb = _b[boff + i];
            if (xa != xb) {
                double diff = Math.abs(xa - xb);
                if (!((diff <= relTol * Math.max(Math.abs(xa), Math.abs(xb))) || (diff <= absTol))) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Constructs a random list of length {@code size} with random values
     * uniformly distributed between {@code min} and {@code max}.
     * 
     * @param min
     *            lower bound of the uniform distribution
     * @param max
     *            upper bound of the uniform distribution
     * @param size
     *            length of the list
     * @return a random list of length {@code size}
     */
    static DoubleList uniformRandom(double min, double max, int size) {
        return DoubleArrayList.uniformRandom(min, max, size);
    }

    /**
     * Constructs a random list of length {@code size} with random values
     * normally distributed with mean {@code mu} and standard deviation
     * {@code sigma}.
     * 
     * @param mu
     *            mean (expectation) of the normal distribution
     * @param sigma
     *            standard deviation of the normal distribution, must be
     *            {@code > 0.0}
     * @param size
     *            length of the list
     * @return a random list of length {@code size}
     */
    static DoubleList uniformNormal(double mu, double sigma, int size) {
        return DoubleArrayList.uniformNormal(mu, sigma, size);
    }

    /**
     * Returns the number of elements in this list.
     *
     * @return the number of elements in this list
     */
    int size();

    /**
     * Returns {@code true} if this list contains no elements.
     *
     * @return {@code true} if this list contains no elements
     */
    boolean isEmpty();

    /**
     * Returns {@code true} if this list contains the specified element. More
     * formally, returns {@code true} if and only if this list contains at least
     * one element {@code e} such that {@code o == e}.
     *
     * @param o
     *            element whose presence in this list is to be tested
     * @return {@code true} if this list contains the specified element
     */
    boolean contains(double o);

    /**
     * Returns {@code true} if this DoubleList contains all of the elements in
     * the specified DoubleList.
     *
     * @param c
     *            DoubleList to be checked for containment in this DoubleList
     * @return {@code true} if this DoubleList contains all of the elements in
     *         the specified DoubleList
     */
    boolean containsAll(DoubleList c);

    /**
     * Returns the index of the first occurrence of the specified element in
     * this list, or -1 if this list does not contain the element. More
     * formally, returns the lowest index {@code i} such that
     * {@code o == get(i)}, or -1 if there is no such index.
     * 
     * @param o
     *            element to search for
     * @return index of the first occurrence, or {@code -1} if not found
     */
    int indexOf(double o);

    /**
     * Returns the index of the last occurrence of the specified element in this
     * list, or -1 if this list does not contain the element. More formally,
     * returns the highest index {@code i} such that {@code o == get(i)}, or -1
     * if there is no such index.
     * 
     * @param o
     *            element to search for
     * @return index of the last occurrence, or {@code -1} if not found
     */
    int lastIndexOf(double o);

    /**
     * Returns an array containing all of the elements in this list in proper
     * sequence (from first to last element).
     *
     * <p>
     * The returned array will be "safe" in that no references to it are
     * maintained by this list. (In other words, this method must allocate a new
     * array). The caller is thus free to modify the returned array.
     *
     * @return an array containing all of the elements in this list in proper
     *         sequence
     */
    double[] toArray();

    /**
     * Returns the element at the specified position in this list.
     *
     * @param index
     *            index of the element to return
     * @return the element at the specified position in this list
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     */
    double get(int index);

    /**
     * Replaces the element at the specified position in this list with the
     * specified element.
     *
     * @param index
     *            index of the element to replace
     * @param element
     *            element to be stored at the specified position
     * @return the element previously at the specified position
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     */
    double set(int index, double element);

    /**
     * Appends the specified element to the end of this list.
     *
     * @param e
     *            element to be appended to this list
     * @return {@code true}
     */
    boolean add(double e);

    /**
     * Inserts the specified element at the specified position in this list.
     * Shifts the element currently at that position (if any) and any subsequent
     * elements to the right (adds one to their indices).
     *
     * @param index
     *            index at which the specified element is to be inserted
     * @param element
     *            element to be inserted
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     */
    void add(int index, double element);

    /**
     * Removes the element at the specified position in this list. Shifts any
     * subsequent elements to the left (subtracts one from their indices).
     *
     * @param index
     *            the index of the element to be removed
     * @return the element that was removed from the list
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     */
    double remove(int index);

    /**
     * Compares the specified object with this DoubleList for equality. Returns
     * {@code true} if and only if the specified object is also a DoubleList,
     * both lists have the same size, and all corresponding pairs of elements in
     * the two lists are <i>equal</i>. In other words, two DoubleList are
     * defined to be equal if they contain the same elements in the same order.
     *
     * @param o
     *            the object to be compared for equality with this DoubleList
     * @return {@code true} if the specified object is equal to this DoubleList
     */
    boolean equals(Object o);

    /**
     * Returns the hash code value for this DoubleList. The hash code of a
     * DoubleList is defined to be the result of the following calculation:
     * 
     * <pre>{@code
     *  int hashCode = 1;
     *  for (int i = 0; i < list.size(); i++) {
     *      hashCode = 31 * hashCode + Double.hashCode(list.get(i));
     *  }
     * }
     * </pre>
     * 
     * This ensures that {@code list1.equals(list2)} implies that
     * {@code list1.hashCode()==list2.hashCode()} for any two DoubleLists,
     * {@code list1} and {@code list2}, as required by the general contract of
     * {@link Object#hashCode}.
     *
     * @return the hash code value for this DoubleList
     */
    int hashCode();

    /**
     * Removes the first occurrence of the specified element from this list, if
     * it is present. If the list does not contain the element, it is unchanged.
     * More formally, removes the element with the lowest index {@code i} such
     * that {@code o == get(i)} (if such an element exists). Returns
     * {@code true} if this list contained the specified element (or
     * equivalently, if this list changed as a result of the call).
     *
     * @param o
     *            element to be removed from this list, if present
     * @return {@code true} if this list contained the specified element
     */
    boolean remove(double o);

    /**
     * Removes all of the elements from this list. The list will be empty after
     * this call returns.
     */
    void clear();

    /**
     * Appends all of the elements in the specified list to the end of this
     * list, in the order that they are returned by the specified list's
     * DIterator. The behavior of this operation is undefined if the specified
     * list is modified while the operation is in progress. (This implies that
     * the behavior of this call is undefined if the specified list is this
     * list, and this list is nonempty.)
     *
     * @param c
     *            list containing elements to be added to this list
     * @return {@code true} if this list changed as a result of the call
     * @throws NullPointerException
     *             if the specified list is null
     */
    boolean addAll(DoubleList c);

    /**
     * Inserts all of the elements in the specified list into this list,
     * starting at the specified position. Shifts the element currently at that
     * position (if any) and any subsequent elements to the right (increases
     * their indices). The new elements will appear in the list in the order
     * that they are returned by the specified list's iterator.
     *
     * @param index
     *            index at which to insert the first element from the specified
     *            list
     * @param c
     *            list containing elements to be added to this list
     * @return {@code true} if this list changed as a result of the call
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     * @throws NullPointerException
     *             if the specified list is null
     */
    boolean addAll(int index, DoubleList c);

    /**
     * Removes from this list all of its elements that are contained in the
     * specified list.
     *
     * @param c
     *            list containing elements to be removed from this list
     * @return {@code true} if this list changed as a result of the call
     */
    boolean removeAll(DoubleList c);

    /**
     * Retains only the elements in this list that are contained in the
     * specified list. In other words, removes from this list all of its
     * elements that are not contained in the specified list.
     *
     * @param c
     *            list containing elements to be retained in this list
     * @return {@code true} if this list changed as a result of the call
     */
    boolean retainAll(DoubleList c);

    /**
     * Returns a list iterator over the elements in this list (in proper
     * sequence), starting at the specified position in the list. The specified
     * index indicates the first element that would be returned by an initial
     * call to {@link DListIterator#next next}. An initial call to
     * {@link DListIterator#previous previous} would return the element with the
     * specified index minus one.
     *
     * <p>
     * The returned list iterator is <i>fail-fast</i>.
     *
     * @param index
     *            starting position
     * @return a list iterator
     * @throws IndexOutOfBoundsException
     *             if the index is out of range
     *             ({@code index < 0 || index >= size()})
     */
    DListIterator listIterator(int index);

    /**
     * Returns a list iterator over the elements in this list (in proper
     * sequence).
     *
     * <p>
     * The returned list iterator is <i>fail-fast</i>.
     * 
     * @return a list iterator
     */
    DListIterator listIterator();

    /**
     * Returns an iterator over the elements in this list in proper sequence.
     *
     * <p>
     * The returned iterator is <i>fail-fast</i>.
     *
     * @return an iterator over the elements in this list in proper sequence
     */
    DForEachIterator iterator();

    /**
     * Returns a view of the portion of this list between the specified
     * {@code fromIndex}, inclusive, and {@code toIndex}, exclusive. (If
     * {@code fromIndex} and {@code toIndex} are equal, the returned list is
     * empty.) The returned list is backed by this list, so non-structural
     * changes in the returned list are reflected in this list, and vice-versa.
     * The returned list supports all of the optional list operations.
     *
     * <p>
     * This method eliminates the need for explicit range operations (of the
     * sort that commonly exist for arrays). Any operation that expects a list
     * can be used as a range operation by passing a subList view instead of a
     * whole list. For example, the following idiom removes a range of elements
     * from a list:
     * 
     * <pre>
     * list.subList(from, to).clear();
     * </pre>
     * 
     * Similar idioms may be constructed for {@link #indexOf(double)} and
     * {@link #lastIndexOf(double)}, and all of the algorithms in a
     * {@code DoubleList} can be applied to a subList.
     *
     * <p>
     * The semantics of the list returned by this method become undefined if the
     * backing list (i.e., this list) is <i>structurally modified</i> in any way
     * other than via the returned list. (Structural modifications are those
     * that change the size of this list, or otherwise perturb it in such a
     * fashion that iterations in progress may yield incorrect results.)
     *
     * @param fromIndex
     *            start index, inclusive
     * @param toIndex
     *            end index, exclusive
     * @return a view of a portion of this list
     * @throws IndexOutOfBoundsException
     *             if an endpoint index value is out of range
     *             {@code (fromIndex < 0 || toIndex > size)}
     * @throws IllegalArgumentException
     *             if the endpoint indices are out of order
     *             {@code (fromIndex > toIndex)}
     */
    DoubleList subList(int fromIndex, int toIndex);

    /**
     * Creates a <em>late-binding</em> and <em>fail-fast</em>
     * {@link Spliterator} over the elements in this list.
     *
     * <p>
     * The {@code Spliterator.OfDouble} reports {@link Spliterator#SIZED},
     * {@link Spliterator#SUBSIZED}, and {@link Spliterator#ORDERED}. Overriding
     * implementations should document the reporting of additional
     * characteristic values.
     *
     * @return a {@code Spliterator.OfDouble} over the elements in this list
     */
    Spliterator.OfDouble spliterator();

    /**
     * Returns the offset of this list.
     * 
     * @return the offset of this list
     */
    int offset();

    /**
     * Returns the underlying array as an escape hatch. Note that the underlying
     * array may be the root array in the case of SubLists.
     * 
     * @return the underlying array
     */
    double[] getArrayUnsafe();

    /**
     * Assign the constant {@code val} to each element in this list.
     * 
     * @param val
     *            the constant to assign
     * @return this {@code DoubleList}
     */
    DoubleList assignConst(double val);

    /**
     * Add the scalar {@code val} to each element in this list.
     * 
     * @param val
     *            the scalar to add
     * @return this {@code DoubleList}
     */
    DoubleList plus(double val);

    /**
     * Multiply each element in this list by the scalar {@code val}.
     * 
     * @param val
     *            the scalar to multiply by
     * @return this {@code DoubleList}
     */
    DoubleList mul(double val);

    /**
     * Returns the minimum of the elements in this list if this list has at
     * least one element, otherwise throws {@code NoSuchElementException}.
     * 
     * @return the minimum of the elements in this list
     * @throws java.util.NoSuchElementException
     *             if this list is empty
     */
    double min();

    /**
     * Returns the maximum of the elements in this list if this list has at
     * least one element, otherwise throws {@code NoSuchElementException}.
     * 
     * @return the maximum of the elements in this list
     * @throws java.util.NoSuchElementException
     *             if this list is empty
     */
    double max();

    /**
     * Returns the average value of the elements in this list if this list has
     * at least one element, otherwise throws {@code NoSuchElementException}.
     * 
     * @return the average value of the elements in this list
     * @throws java.util.NoSuchElementException
     *             if this list is empty
     */
    double avg();

    /**
     * Returns the median value of the elements in this list if this list has at
     * least one element, otherwise throws {@code NoSuchElementException}.
     * 
     * @return the median value of the elements in this list
     * @throws java.util.NoSuchElementException
     *             if this list is empty
     */
    double median();

    /**
     * Returns the sum of the elements in this list.
     * 
     * @return the sum of the elements in this list
     */
    double sum();

    /**
     * Add the scalar {@code val} to the element at {@code index}.
     * 
     * @param index
     *            index of the element to increment by {@code val}
     * @param val
     *            the scalar to add to the element
     * @return the new value of the element after addition
     */
    double plusi(int index, double val);

    /**
     * Multiply the element at {@code index} by the scalar {@code val}.
     * 
     * @param index
     *            index of the element to multiply by {@code val}
     * @param val
     *            the scalar to multiply by
     * @return the new value of the element after multiplication
     */
    double muli(int index, double val);

    /**
     * Computes the dot product of this list and the given {@code list}.
     * 
     * @param list
     *            a DoubleList which must have the same {@link #size()} as this
     *            list
     * @return the dot product of this list and and the given {@code list}
     */
    double dot(DoubleList list);

    /**
     * Computes the dot product of this list and the given {@code array}.
     * 
     * @param array
     *            an array which must have the same length as this list's
     *            {@link #size()}
     * @return the dot product of this list and and the given {@code array}.
     */
    double dot(double[] array);

    /**
     * Computes the <a href="https://en.wikipedia.org/wiki/Cross_product">cross
     * product</a> of this DoubleList and the given {@code list} provided that
     * both, this DoubleList and the given {@code list}, have a dimension of
     * <b>3</b>. If this DoubleList is not of {@link #size()} {@code == 3} an
     * {@code UnsupportedOperationException} is thrown. Also, if the given
     * {@code list} doesn't satisfy {@code list.size() == 3} an
     * {@code IllegalArgumentException} is thrown.
     * 
     * @param list
     *            a DoubleList which must have the {@code size() == 3}
     * @return the the cross product of this DoubleList and the given
     *         {@code list}
     * @throws java.lang.UnsupportedOperationException
     *             if this DoubleList is not of {@code size() == 3}
     * @throws java.lang.IllegalArgumentException
     *             if the given {@code list} doesn't satisfy
     *             {@code list.size() == 3}
     */
    DoubleList cross(DoubleList list);

    /**
     * Computes the <a href="https://en.wikipedia.org/wiki/Cross_product">cross
     * product</a> of this DoubleList and the given {@code array} provided that
     * both, this DoubleList and the given {@code array}, have a dimension of
     * <b>3</b>. If this DoubleList is not of {@link #size()} {@code == 3} an
     * {@code UnsupportedOperationException} is thrown. Also, if the given
     * {@code array} doesn't satisfy {@code array.length == 3} an
     * {@code IllegalArgumentException} is thrown.
     * 
     * @param array
     *            an array which must have {@code length == 3}
     * @return the cross product of this DoubleList and the given {@code array}
     * @throws java.lang.UnsupportedOperationException
     *             if this DoubleList is not of {@code size() == 3}
     * @throws java.lang.IllegalArgumentException
     *             if the given {@code array} doesn't satisfy
     *             {@code array.length == 3}
     */
    DoubleList cross(double[] array);

    /**
     * Returns a new list that contains the contents of this list normalized to
     * a probability distribution according to the
     * <a href="https://en.wikipedia.org/wiki/Softmax_function">softmax</a>
     * function which sums up to {@code 1.0}.
     * 
     * @return a new list that contains the contents of this list normalized
     *         according to the softmax function
     */
    DoubleList softmax();

    /**
     * Returns a new list that contains the contents of this list normalized to
     * a distribution according to the
     * <a href="https://en.wikipedia.org/wiki/Softmax_function">softmax</a>
     * function which sums up to {@code normalizedSum}.
     * 
     * @param normalizedSum
     *            the overall sum of the desired distribution
     * @return a new list that contains the contents of this list normalized
     *         according to the softmax function which sums up to
     *         {@code normalizedSum}
     */
    DoubleList normalizedSoftmax(double normalizedSum);

    /**
     * Computes the euclidean norm (magnitude) of this list.
     * 
     * @return the euclidean norm of this list
     */
    double norm2();

    /**
     * Returns the logarithm of the sum of the exponentials of this list (see
     * <a href="https://en.wikipedia.org/wiki/LogSumExp">LogSumExp</a>) if it
     * has at least one element, otherwise {@code NoSuchElementException} is
     * thrown.
     * 
     * @return the logarithm of the sum of the exponentials of this list
     * @throws java.util.NoSuchElementException
     *             if this list is empty
     */
    double logSumExp();

    /**
     * Add the scalars from the given {@code list} to the elements in this list.
     * Note that {@code list} doesn't have to have the same {@link #size()} as
     * this list as the minimum of both list sizes is used to determine the
     * number of scalars to add.
     * 
     * @param list
     *            a DoubleList of scalars to add to this list
     * @return this {@code DoubleList}
     */
    DoubleList plusn(DoubleList list);

    /**
     * Add the scalars from the given {@code array} to the elements in this
     * list. Note that {@code array} doesn't have to have the same length as
     * this list as the minimum of both lengths is used to determine the number
     * of scalars to add.
     * 
     * @param array
     *            an array of scalars to add to this list
     * @return this {@code DoubleList}
     */
    DoubleList plusn(double[] array);

    /**
     * Multiplies the scalars from this list by the scalars from the given
     * {@code list}. Note that {@code list} doesn't have to have the same
     * {@link #size()} as this list as the minimum of both list sizes is used to
     * determine the number of scalars to multiply.
     * 
     * @param list
     *            a DoubleList of scalars to multiply by the elements in this
     *            list
     * @return this {@code DoubleList}
     */
    DoubleList muln(DoubleList list);

    /**
     * Multiplies the scalars from this list by the scalars from the given
     * {@code array}. Note that {@code array} doesn't have to have the same
     * length as this list as the minimum of both lengths is used to determine
     * the number of scalars to multiply.
     * 
     * @param array
     *            an array of scalars to multiply by the elements in this list
     * @return this {@code DoubleList}
     */
    DoubleList muln(double[] array);

    /**
     * Round the scalars in this list to {@code scale}.
     * 
     * @param scale
     *            if zero or positive, the scale is the number of digits to the
     *            right of the decimal point. If negative, the unscaled value of
     *            the number is multiplied by ten to the power of the negation
     *            of the scale.
     * @return this {@code DoubleList}
     */
    DoubleList round(int scale);

    /**
     * Set all scalars that are either NaN, positive or negative infinity to the
     * respective ersatz value provided by the {@code nanSurrogate},
     * {@code posInfSurrogate} and {@code negInfSurrogate} arguments.
     * 
     * @param nanSurrogate
     *            the substitution value to use for NaN values, NaNs won't get
     *            substituted if you pass {@link Double#NaN}
     * @param posInfSurrogate
     *            the substitution value to use for positive infinity values, no
     *            substitution happens if you pass
     *            {@link Double#POSITIVE_INFINITY}
     * @param negInfSurrogate
     *            the substitution value to use for negative infinity values, no
     *            substitution happens if you pass
     *            {@link Double#NEGATIVE_INFINITY}
     * @return this {@code DoubleList}
     */
    DoubleList sanitizeNonFinite(double nanSurrogate, double posInfSurrogate, double negInfSurrogate);
}
