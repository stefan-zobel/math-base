/*
 * Copyright 2021, 2023 Stefan Zobel
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Spliterator;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.DoubleConsumer;
import java.util.function.DoublePredicate;

import math.fun.DForEachIterator;

/**
 * Resizable primitive double[] array implementation. This is essentially a
 * shameless stripped-down copy of {@link ArrayList} specialized for primitive
 * doubles and augmented with mathematical vector operations.
 */
public class DoubleArrayList implements DoubleList, Cloneable, Externalizable {
    /**
     * Default initial capacity.
     */
    private static final int DEFAULT_CAPACITY = 10;

    /**
     * Shared empty array instance used for empty instances.
     */
    private static final double[] EMPTY_ELEMENTDATA = {};

    /**
     * Shared empty array instance used for default sized empty instances. We
     * distinguish this from EMPTY_ELEMENTDATA to know how much to inflate when
     * first element is added.
     */
    private static final double[] DEFAULTCAPACITY_EMPTY_ELEMENTDATA = {};

    /**
     * The array buffer into which the elements of the DoubleArrayList are
     * stored. The capacity of the DoubleArrayList is the length of this array
     * buffer. Any empty DoubleArrayList with elementData ==
     * DEFAULTCAPACITY_EMPTY_ELEMENTDATA will be expanded to DEFAULT_CAPACITY
     * when the first element is added.
     */
    double[] elementData; // non-private to simplify nested class access

    /**
     * The size of the DoubleArrayList (the number of elements it contains).
     */
    int size;

    int modCount = 0;

    /**
     * Constructs an empty list with the specified initial capacity.
     *
     * @param initialCapacity
     *            the initial capacity of the list
     * @throws IllegalArgumentException
     *             if the specified initial capacity is negative
     */
    public DoubleArrayList(int initialCapacity) {
        if (initialCapacity > 0) {
            elementData = new double[initialCapacity];
        } else if (initialCapacity == 0) {
            elementData = EMPTY_ELEMENTDATA;
        } else {
            throw new IllegalArgumentException("Illegal Capacity: " + initialCapacity);
        }
    }

    /**
     * Constructs a list with length {@code initialLength} where all elements
     * have {@code initialValue} as initial value.
     * 
     * @param initialLength
     *            the initial length of the list
     * @param initialValue
     *            the initial value of all elements in the list
     * @throws IllegalArgumentException
     *             if the specified initial length is negative
     */
    public DoubleArrayList(int initialLength, double initialValue) {
        if (initialLength > 0) {
            double[] es = new double[initialLength];
            Arrays.fill(es, initialValue);
            elementData = es;
            size = es.length;
        } else if (initialLength == 0) {
            elementData = EMPTY_ELEMENTDATA;
        } else {
            throw new IllegalArgumentException("Illegal Length: " + initialLength);
        }
    }

    /**
     * Constructs an empty list with an initial capacity of ten.
     */
    public DoubleArrayList() {
        elementData = DEFAULTCAPACITY_EMPTY_ELEMENTDATA;
    }

    /**
     * Constructs a list containing the elements of the specified array, in the
     * order they appear in the array.
     *
     * @param c
     *            the array whose elements are to be placed into this list
     * @throws NullPointerException
     *             if the specified array is null
     */
    public DoubleArrayList(double[] c) {
        this(c, true);
    }

    DoubleArrayList(double[] c, boolean doCopy) {
        if (doCopy) {
            if ((size = c.length) != 0) {
                elementData = Arrays.copyOf(c, size);
            } else {
                // replace with empty array.
                elementData = EMPTY_ELEMENTDATA;
            }
        } else {
            size = c.length;
            elementData = c;
        }
    }

    /**
     * Constructs a list containing the passed values in the order of their
     * appearance.
     * 
     * @param values
     *            the values which are to be placed into the list
     * @return a new {@code DoubleArrayList} constructed from the passed values
     */
    public static DoubleArrayList of(double... values) {
        return new DoubleArrayList(values, false);
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
    public static DoubleArrayList randomUniform(double min, double max, int size) {
        double[] values = new double[size];
        ThreadLocalRandom rng = ThreadLocalRandom.current();
        double spread = max - min;
        for (int i = 0; i < values.length; ++i) {
            values[i] = min + spread * rng.nextDouble();
        }
        return new DoubleArrayList(values, false);
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
    public static DoubleArrayList randomNormal(double mu, double sigma, int size) {
        if (sigma <= 0.0) {
            throw new IllegalArgumentException("Standard deviation must be positive (" + sigma + ")");
        }
        double[] values = new double[size];
        ThreadLocalRandom rng = ThreadLocalRandom.current();
        for (int i = 0; i < values.length; ++i) {
            values[i] = mu + sigma * rng.nextGaussian();
        }
        return new DoubleArrayList(values, false);
    }

    /**
     * Constructs a list containing the elements of the specified list, in the
     * order they are returned by the list's {@link DoubleList#toArray()}
     * method.
     *
     * @param c
     *            the list whose elements are to be placed into this list
     * @throws NullPointerException
     *             if the specified list is null
     */
    public DoubleArrayList(DoubleList c) {
        this(c.toArray());
    }

    /**
     * Trims the capacity of this {@code DoubleArrayList} instance to be the
     * list's current size. An application can use this operation to minimize
     * the storage of an {@code DoubleArrayList} instance.
     */
    public void trimToSize() {
        modCount++;
        if (size < elementData.length) {
            elementData = (size == 0) ? EMPTY_ELEMENTDATA : Arrays.copyOf(elementData, size);
        }
    }

    /**
     * Increases the capacity of this {@code DoubleArrayList} instance, if
     * necessary, to ensure that it can hold at least the number of elements
     * specified by the minimum capacity argument.
     *
     * @param minCapacity
     *            the desired minimum capacity
     */
    public void ensureCapacity(int minCapacity) {
        if (minCapacity > elementData.length
                && !(elementData == DEFAULTCAPACITY_EMPTY_ELEMENTDATA && minCapacity <= DEFAULT_CAPACITY)) {
            modCount++;
            grow(minCapacity);
        }
    }

    /**
     * The maximum size of array to allocate (unless necessary). Some VMs
     * reserve some header words in an array. Attempts to allocate larger arrays
     * may result in OutOfMemoryError: Requested array size exceeds VM limit
     */
    private static final int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 8;

    /**
     * Increases the capacity to ensure that it can hold at least the number of
     * elements specified by the minimum capacity argument.
     *
     * @param minCapacity
     *            the desired minimum capacity
     * @throws OutOfMemoryError
     *             if minCapacity is less than zero
     */
    private double[] grow(int minCapacity) {
        return elementData = Arrays.copyOf(elementData, newCapacity(minCapacity));
    }

    private double[] grow() {
        return grow(size + 1);
    }

    /**
     * Returns a capacity at least as large as the given minimum capacity.
     * Returns the current capacity increased by 50% if that suffices. Will not
     * return a capacity greater than MAX_ARRAY_SIZE unless the given minimum
     * capacity is greater than MAX_ARRAY_SIZE.
     *
     * @param minCapacity
     *            the desired minimum capacity
     * @throws OutOfMemoryError
     *             if minCapacity is less than zero
     */
    private int newCapacity(int minCapacity) {
        // overflow-conscious code
        int oldCapacity = elementData.length;
        int newCapacity = oldCapacity + (oldCapacity >> 1);
        if (newCapacity - minCapacity <= 0) {
            if (elementData == DEFAULTCAPACITY_EMPTY_ELEMENTDATA) {
                return Math.max(DEFAULT_CAPACITY, minCapacity);
            }
            if (minCapacity < 0) { // overflow // "lgtm[java/constant-comparison]"
                throw new OutOfMemoryError();
            }
            return minCapacity;
        }
        return (newCapacity - MAX_ARRAY_SIZE <= 0) ? newCapacity : hugeCapacity(minCapacity);
    }

    private static int hugeCapacity(int minCapacity) {
        if (minCapacity < 0) { // overflow
            throw new OutOfMemoryError();
        }
        return (minCapacity > MAX_ARRAY_SIZE) ? Integer.MAX_VALUE : MAX_ARRAY_SIZE;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int size() {
        return size;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean contains(double o) {
        return indexOf(o) >= 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean containsAll(DoubleList c) {
        for (int i = 0; i < c.size(); ++i) {
            if (!contains(c.get(i))) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int indexOf(double o) {
        return indexOfRange(o, 0, size);
    }

    int indexOfRange(double o, int start, int end) {
        double[] es = elementData;
        for (int i = start; i < end; i++) {
            if (o == es[i]) {
                return i;
            }
        }
        return -1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int lastIndexOf(double o) {
        return lastIndexOfRange(o, 0, size);
    }

    int lastIndexOfRange(double o, int start, int end) {
        double[] es = elementData;
        for (int i = end - 1; i >= start; i--) {
            if (o == es[i]) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Returns a deep copy of this {@code DoubleArrayList} instance.
     *
     * @return a clone of this {@code DoubleArrayList} instance
     */
    @Override
    public Object clone() {
        try {
            DoubleArrayList v = (DoubleArrayList) super.clone();
            v.elementData = Arrays.copyOf(elementData, size);
            v.modCount = 0;
            return v;
        } catch (CloneNotSupportedException e) {
            // this shouldn't happen, since we are Cloneable
            throw new InternalError(e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] toArray() {
        return Arrays.copyOf(elementData, size);
    }

    // Positional Access Operations

    double elementData(int index) {
        return elementData[index];
    }

    static double elementAt(double[] es, int index) {
        return es[index];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double get(int index) {
        checkIndex(index, size);
        return elementData(index);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double set(int index, double element) {
        checkIndex(index, size);
        double oldValue = elementData(index);
        elementData[index] = element;
        return oldValue;
    }

    /**
     * {@inheritDoc}
     */
    private void add(double e, double[] elementData, int s) {
        if (s == elementData.length) {
            elementData = grow();
        }
        elementData[s] = e;
        size = s + 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean add(double e) {
        modCount++;
        add(e, elementData, size);
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void add(int index, double element) {
        rangeCheckForAdd(index);
        modCount++;
        final int s;
        double[] elementData;
        if ((s = size) == (elementData = this.elementData).length) {
            elementData = grow();
        }
        System.arraycopy(elementData, index, elementData, index + 1, s - index);
        elementData[index] = element;
        size = s + 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double remove(int index) {
        checkIndex(index, size);
        final double[] es = elementData;

        double oldValue = es[index];
        fastRemove(es, index);

        return oldValue;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (!(o instanceof DoubleList)) {
            return false;
        }

        final int expectedModCount = modCount;
        // DoubleArrayList can be subclassed and given arbitrary behavior, but
        // we can still deal with the common case where o is DoubleArrayList
        // precisely
        boolean equal = (o.getClass() == DoubleArrayList.class) ? equalsDoubleList((DoubleArrayList) o)
                : equalsRange((DoubleList) o, 0, size);

        checkForComodification(expectedModCount);
        return equal;
    }

    boolean equalsRange(DoubleList other, int from, int to) {
        final double[] es = elementData;
        if (to > es.length) {
            throw new ConcurrentModificationException();
        }
        DForEachIterator oit = other.iterator();
        for (; from < to; from++) {
            if (!oit.hasNext() || es[from] != oit.next()) {
                return false;
            }
        }
        return !oit.hasNext();
    }

    private boolean equalsDoubleList(DoubleArrayList other) {
        final int otherModCount = other.modCount;
        final int s = size;
        boolean equal;
        if (equal = (s == other.size)) {
            final double[] otherEs = other.elementData;
            final double[] es = elementData;
            if (s > es.length || s > otherEs.length) {
                throw new ConcurrentModificationException();
            }
            for (int i = 0; i < s; i++) {
                if (es[i] != otherEs[i]) {
                    equal = false;
                    break;
                }
            }
        }
        other.checkForComodification(otherModCount);
        return equal;
    }

    private void checkForComodification(final int expectedModCount) {
        if (modCount != expectedModCount) {
            throw new ConcurrentModificationException();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int expectedModCount = modCount;
        int hash = hashCodeRange(0, size);
        checkForComodification(expectedModCount);
        return hash;
    }

    int hashCodeRange(int from, int to) {
        final double[] es = elementData;
        if (to > es.length) {
            throw new ConcurrentModificationException();
        }
        int hashCode = 1;
        for (int i = from; i < to; i++) {
            hashCode = 31 * hashCode + Double.hashCode(es[i]);
        }
        return hashCode;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean remove(double o) {
        final double[] es = elementData;
        final int size = this.size;
        int i = 0;
        found: {
            for (; i < size; i++) {
                if (o == es[i]) {
                    break found;
                }
            }
            return false;
        }
        fastRemove(es, i);
        return true;
    }

    /**
     * Private remove method that skips bounds checking and does not return the
     * value removed.
     */
    private void fastRemove(double[] es, int i) {
        modCount++;
        final int newSize;
        if ((newSize = size - 1) > i) {
            System.arraycopy(es, i + 1, es, i, newSize - i);
        }
        size = newSize;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void clear() {
        modCount++;
        size = 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean addAll(DoubleList c) {
        double[] a = c.toArray();
        modCount++;
        int numNew = a.length;
        if (numNew == 0) {
            return false;
        }
        double[] elementData;
        final int s;
        if (numNew > (elementData = this.elementData).length - (s = size)) {
            elementData = grow(s + numNew);
        }
        System.arraycopy(a, 0, elementData, s, numNew);
        size = s + numNew;
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean addAll(int index, DoubleList c) {
        rangeCheckForAdd(index);

        double[] a = c.toArray();
        modCount++;
        int numNew = a.length;
        if (numNew == 0) {
            return false;
        }
        double[] elementData;
        final int s;
        if (numNew > (elementData = this.elementData).length - (s = size)) {
            elementData = grow(s + numNew);
        }

        int numMoved = s - index;
        if (numMoved > 0) {
            System.arraycopy(elementData, index, elementData, index + numNew, numMoved);
        }
        System.arraycopy(a, 0, elementData, index, numNew);
        size = s + numNew;
        return true;
    }

    /**
     * Removes from this list all of the elements whose index is between
     * {@code fromIndex}, inclusive, and {@code toIndex}, exclusive. Shifts any
     * succeeding elements to the left (reduces their index). This call shortens
     * the list by {@code (toIndex - fromIndex)} elements. (If
     * {@code toIndex==fromIndex}, this operation has no effect.)
     *
     * @param fromIndex
     *            start index, inclusive
     * @param toIndex
     *            end index, exclusive
     * @throws IndexOutOfBoundsException
     *             if {@code fromIndex} or {@code toIndex} is out of range
     *             ({@code fromIndex < 0 ||
     *          toIndex > size() ||
     *          toIndex < fromIndex})
     */
    protected void removeRange(int fromIndex, int toIndex) {
        if (fromIndex > toIndex) {
            throw new IndexOutOfBoundsException(outOfBoundsMsg(fromIndex, toIndex));
        }
        modCount++;
        shiftTailOverGap(elementData, fromIndex, toIndex);
    }

    /** Erases the gap from lo to hi, by sliding down following elements. */
    private void shiftTailOverGap(double[] es, int lo, int hi) {
        System.arraycopy(es, hi, es, lo, size - hi);
        size -= hi - lo;
    }

    /**
     * A version of rangeCheck used by add and addAll.
     */
    private void rangeCheckForAdd(int index) {
        if (index > size || index < 0) {
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
        }
    }

    /**
     * Constructs an IndexOutOfBoundsException detail message. Of the many
     * possible refactorings of the error handling code, this "outlining"
     * performs best with both server and client VMs.
     */
    private String outOfBoundsMsg(int index) {
        return "Index: " + index + ", Size: " + size;
    }

    /**
     * A version used in checking (fromIndex > toIndex) condition
     */
    private static String outOfBoundsMsg(int fromIndex, int toIndex) {
        return "From Index: " + fromIndex + " > To Index: " + toIndex;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList filter(DoublePredicate predicate) {
        return filter(size, 0, elementData, predicate);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList shuffle() {
        shuffle(size, 0, elementData);
        return this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean removeAll(DoubleList c) {
        return batchRemove(c, false, 0, size);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean retainAll(DoubleList c) {
        return batchRemove(c, true, 0, size);
    }

    boolean batchRemove(DoubleList c, boolean complement, final int from, final int end) {
        Objects.requireNonNull(c);
        final double[] es = elementData;
        int r;
        // Optimize for initial run of survivors
        for (r = from;; r++) {
            if (r == end) {
                return false;
            }
            if (c.contains(es[r]) != complement) {
                break;
            }
        }
        int w = r++;
        try {
            for (double e; r < end; r++) {
                if (c.contains(e = es[r]) == complement) {
                    es[w++] = e;
                }
            }
        } catch (Throwable ex) {
            // Preserve behavioral compatibility with AbstractCollection,
            // even if c.contains() throws.
            System.arraycopy(es, r, es, w, end - r);
            w += end - r;
            throw ex;
        } finally {
            modCount += end - w;
            shiftTailOverGap(es, w, end);
        }
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DListIterator listIterator(int index) {
        rangeCheckForAdd(index);
        return new ListItr(index);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DListIterator listIterator() {
        return new ListItr(0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DForEachIterator iterator() {
        return new Itr();
    }

    /**
     * An optimized version of AbstractDoubleList.Itr
     */
    private class Itr implements DForEachIterator {
        int cursor; // index of next element to return
        int lastRet = -1; // index of last element returned; -1 if no such
        int expectedModCount = modCount;

        // prevent creating a synthetic constructor
        Itr() {
        }

        public boolean hasNext() {
            return cursor != size;
        }

        public double next() {
            checkForComodification();
            int i = cursor;
            if (i >= size) {
                throw new NoSuchElementException();
            }
            double[] elementData = DoubleArrayList.this.elementData;
            if (i >= elementData.length) {
                throw new ConcurrentModificationException();
            }
            cursor = i + 1;
            return elementData[lastRet = i];
        }

        public void remove() {
            if (lastRet < 0) {
                throw new IllegalStateException();
            }
            checkForComodification();

            try {
                DoubleArrayList.this.remove(lastRet);
                cursor = lastRet;
                lastRet = -1;
                expectedModCount = modCount;
            } catch (IndexOutOfBoundsException ex) {
                throw new ConcurrentModificationException();
            }
        }

        @Override
        public void forEachRemaining(DoubleConsumer action) {
            Objects.requireNonNull(action);
            final int size = DoubleArrayList.this.size;
            int i = cursor;
            if (i < size) {
                final double[] es = elementData;
                if (i >= es.length) {
                    throw new ConcurrentModificationException();
                }
                for (; i < size && modCount == expectedModCount; i++) {
                    action.accept(elementAt(es, i));
                }
                // update once at end to reduce heap write traffic
                cursor = i;
                lastRet = i - 1;
                checkForComodification();
            }
        }

        final void checkForComodification() {
            if (modCount != expectedModCount) {
                throw new ConcurrentModificationException();
            }
        }
    }

    /**
     * An optimized version of AbstractDoubleList.ListItr
     */
    private class ListItr extends Itr implements DListIterator {
        ListItr(int index) {
            super();
            cursor = index;
        }

        public boolean hasPrevious() {
            return cursor != 0;
        }

        public int nextIndex() {
            return cursor;
        }

        public int previousIndex() {
            return cursor - 1;
        }

        public double previous() {
            checkForComodification();
            int i = cursor - 1;
            if (i < 0) {
                throw new NoSuchElementException();
            }
            double[] elementData = DoubleArrayList.this.elementData;
            if (i >= elementData.length) {
                throw new ConcurrentModificationException();
            }
            cursor = i;
            return elementData[lastRet = i];
        }

        public void set(double e) {
            if (lastRet < 0) {
                throw new IllegalStateException();
            }
            checkForComodification();

            try {
                DoubleArrayList.this.set(lastRet, e);
            } catch (IndexOutOfBoundsException ex) {
                throw new ConcurrentModificationException();
            }
        }

        public void add(double e) {
            checkForComodification();

            try {
                int i = cursor;
                DoubleArrayList.this.add(i, e);
                cursor = i + 1;
                lastRet = -1;
                expectedModCount = modCount;
            } catch (IndexOutOfBoundsException ex) {
                throw new ConcurrentModificationException();
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList subList(int fromIndex, int toIndex) {
        subListRangeCheck(fromIndex, toIndex, size);
        return new SubList(this, fromIndex, toIndex);
    }

    private abstract static class AbstractDoubleList {
        int modCount;

        public abstract int size();

        public abstract double get(int index);

        public abstract double set(int index, double element);

        public abstract void add(int index, double element);

        public abstract DForEachIterator iterator();

        AbstractDoubleList() {
        }

        private class Itr implements DForEachIterator {
            /**
             * Index of element to be returned by subsequent call to next.
             */
            int cursor = 0;

            /**
             * Index of element returned by most recent call to next or
             * previous. Reset to -1 if this element is deleted by a call to
             * remove.
             */
            int lastRet = -1;

            /**
             * The modCount value that the iterator believes that the backing
             * List should have. If this expectation is violated, the iterator
             * has detected concurrent modification.
             */
            int expectedModCount = modCount;

            Itr() {
            }

            public boolean hasNext() {
                return cursor != size();
            }

            public double next() {
                checkForComodification();
                try {
                    int i = cursor;
                    double next = get(i);
                    lastRet = i;
                    cursor = i + 1;
                    return next;
                } catch (IndexOutOfBoundsException e) {
                    checkForComodification();
                    throw new NoSuchElementException();
                }
            }

            public void remove() {
                if (lastRet < 0) {
                    throw new IllegalStateException();
                }
                checkForComodification();

                try {
                    AbstractDoubleList.this.remove(lastRet);
                    if (lastRet < cursor) {
                        cursor--;
                    }
                    lastRet = -1;
                    expectedModCount = modCount;
                } catch (IndexOutOfBoundsException e) {
                    throw new ConcurrentModificationException();
                }
            }

            final void checkForComodification() {
                if (modCount != expectedModCount) {
                    throw new ConcurrentModificationException();
                }
            }
        }

        private class ListItr extends Itr implements DListIterator {
            ListItr(int index) {
                cursor = index;
            }

            public boolean hasPrevious() {
                return cursor != 0;
            }

            public double previous() {
                checkForComodification();
                try {
                    int i = cursor - 1;
                    double previous = get(i);
                    lastRet = cursor = i;
                    return previous;
                } catch (IndexOutOfBoundsException e) {
                    checkForComodification();
                    throw new NoSuchElementException();
                }
            }

            public int nextIndex() {
                return cursor;
            }

            public int previousIndex() {
                return cursor - 1;
            }

            public void set(double e) {
                if (lastRet < 0) {
                    throw new IllegalStateException();
                }
                checkForComodification();

                try {
                    AbstractDoubleList.this.set(lastRet, e);
                    expectedModCount = modCount;
                } catch (IndexOutOfBoundsException ex) {
                    throw new ConcurrentModificationException();
                }
            }

            public void add(double e) {
                checkForComodification();
                try {
                    int i = cursor;
                    AbstractDoubleList.this.add(i, e);
                    lastRet = -1;
                    cursor = i + 1;
                    expectedModCount = modCount;
                } catch (IndexOutOfBoundsException ex) {
                    throw new ConcurrentModificationException();
                }
            }
        }

        private String outOfBoundsMsg(int index) {
            return "Index: " + index + ", Size: " + size();
        }

        private void rangeCheckForAdd(int index) {
            if (index < 0 || index > size()) {
                throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
            }
        }

        public DListIterator listIterator() {
            return listIterator(0);
        }

        public DListIterator listIterator(int index) {
            rangeCheckForAdd(index);

            return new ListItr(index);
        }

        public boolean remove(double o) {
            DForEachIterator it = iterator();
            while (it.hasNext()) {
                if (o == it.next()) {
                    it.remove();
                    return true;
                }
            }
            return false;
        }

        public boolean add(double e) {
            add(size(), e);
            return true;
        }

        public boolean isEmpty() {
            return size() == 0;
        }

        protected void removeRange(int fromIndex, int toIndex) {
            DListIterator it = listIterator(fromIndex);
            for (int i = 0, n = toIndex - fromIndex; i < n; i++) {
                it.next();
                it.remove();
            }
        }

        public void clear() {
            removeRange(0, size());
        }

        public boolean containsAll(DoubleList c) {
            for (int i = 0; i < c.size(); ++i) {
                if (!contains(c.get(i))) {
                    return false;
                }
            }
            return true;
        }

        public boolean contains(double o) {
            for (int i = 0; i < size(); ++i) {
                if (o == get(i)) {
                    return true;
                }
            }
            return false;
        }
    }

    private static class SubList extends AbstractDoubleList implements DoubleList {
        final DoubleArrayList root;
        private final SubList parent;
        final int offset;
        int size;

        /**
         * Constructs a sublist of an arbitrary DoubleArrayList.
         */
        SubList(DoubleArrayList root, int fromIndex, int toIndex) {
            this.root = root;
            this.parent = null;
            this.offset = fromIndex;
            this.size = toIndex - fromIndex;
            this.modCount = root.modCount;
        }

        /**
         * Constructs a sublist of another SubList.
         */
        private SubList(SubList parent, int fromIndex, int toIndex) {
            this.root = parent.root;
            this.parent = parent;
            this.offset = parent.offset + fromIndex;
            this.size = toIndex - fromIndex;
            this.modCount = parent.modCount;
        }

        public double set(int index, double element) {
            checkIndex(index, size);
            checkForComodification();
            double oldValue = root.elementData(offset + index);
            root.elementData[offset + index] = element;
            return oldValue;
        }

        public double get(int index) {
            checkIndex(index, size);
            checkForComodification();
            return root.elementData(offset + index);
        }

        public int size() {
            checkForComodification();
            return size;
        }

        public void add(int index, double element) {
            rangeCheckForAdd(index);
            checkForComodification();
            root.add(offset + index, element);
            updateSizeAndModCount(1);
        }

        public double remove(int index) {
            checkIndex(index, size);
            checkForComodification();
            double result = root.remove(offset + index);
            updateSizeAndModCount(-1);
            return result;
        }

        protected void removeRange(int fromIndex, int toIndex) {
            checkForComodification();
            root.removeRange(offset + fromIndex, offset + toIndex);
            updateSizeAndModCount(fromIndex - toIndex);
        }

        public boolean addAll(DoubleList c) {
            return addAll(this.size, c);
        }

        public boolean addAll(int index, DoubleList c) {
            rangeCheckForAdd(index);
            int cSize = c.size();
            if (cSize == 0) {
                return false;
            }
            checkForComodification();
            root.addAll(offset + index, c);
            updateSizeAndModCount(cSize);
            return true;
        }

        public DoubleList filter(DoublePredicate predicate) {
            return DoubleArrayList.filter(size, offset, root.elementData, predicate);
        }

        public DoubleList shuffle() {
            DoubleArrayList.shuffle(size, offset, root.elementData);
            return this;
        }

        public boolean removeAll(DoubleList c) {
            return batchRemove(c, false);
        }

        public boolean retainAll(DoubleList c) {
            return batchRemove(c, true);
        }

        private boolean batchRemove(DoubleList c, boolean complement) {
            checkForComodification();
            int oldSize = root.size;
            boolean modified = root.batchRemove(c, complement, offset, offset + size);
            if (modified) {
                updateSizeAndModCount(root.size - oldSize);
            }
            return modified;
        }

        public double[] toArray() {
            checkForComodification();
            return Arrays.copyOfRange(root.elementData, offset, offset + size);
        }

        public void sort() {
            Arrays.sort(root.elementData, offset, offset + size);
            checkForComodification();
        }

        public boolean equals(Object o) {
            if (o == this) {
                return true;
            }
            if (!(o instanceof DoubleList)) {
                return false;
            }
            boolean equal = root.equalsRange((DoubleList) o, offset, offset + size);
            checkForComodification();
            return equal;
        }

        public int hashCode() {
            int hash = root.hashCodeRange(offset, offset + size);
            checkForComodification();
            return hash;
        }

        public int indexOf(double o) {
            int index = root.indexOfRange(o, offset, offset + size);
            checkForComodification();
            return index >= 0 ? index - offset : -1;
        }

        public int lastIndexOf(double o) {
            int index = root.lastIndexOfRange(o, offset, offset + size);
            checkForComodification();
            return index >= 0 ? index - offset : -1;
        }

        public boolean contains(double o) {
            return indexOf(o) >= 0;
        }

        public DForEachIterator iterator() {
            return listIterator();
        }

        public DListIterator listIterator(int index) {
            checkForComodification();
            rangeCheckForAdd(index);

            return new DListIterator() {
                int cursor = index;
                int lastRet = -1;
                int expectedModCount = SubList.this.modCount;

                public boolean hasNext() {
                    return cursor != SubList.this.size;
                }

                public double next() {
                    checkForComodification();
                    int i = cursor;
                    if (i >= SubList.this.size) {
                        throw new NoSuchElementException();
                    }
                    double[] elementData = root.elementData;
                    if (offset + i >= elementData.length) {
                        throw new ConcurrentModificationException();
                    }
                    cursor = i + 1;
                    return elementData[offset + (lastRet = i)];
                }

                public boolean hasPrevious() {
                    return cursor != 0;
                }

                public double previous() {
                    checkForComodification();
                    int i = cursor - 1;
                    if (i < 0) {
                        throw new NoSuchElementException();
                    }
                    double[] elementData = root.elementData;
                    if (offset + i >= elementData.length) {
                        throw new ConcurrentModificationException();
                    }
                    cursor = i;
                    return elementData[offset + (lastRet = i)];
                }

                public void forEachRemaining(DoubleConsumer action) {
                    Objects.requireNonNull(action);
                    final int size = SubList.this.size;
                    int i = cursor;
                    if (i < size) {
                        final double[] es = root.elementData;
                        if (offset + i >= es.length) {
                            throw new ConcurrentModificationException();
                        }
                        for (; i < size && root.modCount == expectedModCount; i++) {
                            action.accept(elementAt(es, offset + i));
                        }
                        // update once at end to reduce heap write traffic
                        cursor = i;
                        lastRet = i - 1;
                        checkForComodification();
                    }
                }

                public int nextIndex() {
                    return cursor;
                }

                public int previousIndex() {
                    return cursor - 1;
                }

                public void remove() {
                    if (lastRet < 0) {
                        throw new IllegalStateException();
                    }
                    checkForComodification();

                    try {
                        SubList.this.remove(lastRet);
                        cursor = lastRet;
                        lastRet = -1;
                        expectedModCount = SubList.this.modCount;
                    } catch (IndexOutOfBoundsException ex) {
                        throw new ConcurrentModificationException();
                    }
                }

                public void set(double e) {
                    if (lastRet < 0) {
                        throw new IllegalStateException();
                    }
                    checkForComodification();

                    try {
                        root.set(offset + lastRet, e);
                    } catch (IndexOutOfBoundsException ex) {
                        throw new ConcurrentModificationException();
                    }
                }

                public void add(double e) {
                    checkForComodification();

                    try {
                        int i = cursor;
                        SubList.this.add(i, e);
                        cursor = i + 1;
                        lastRet = -1;
                        expectedModCount = SubList.this.modCount;
                    } catch (IndexOutOfBoundsException ex) {
                        throw new ConcurrentModificationException();
                    }
                }

                final void checkForComodification() {
                    if (root.modCount != expectedModCount) {
                        throw new ConcurrentModificationException();
                    }
                }
            };
        }

        public DoubleList subList(int fromIndex, int toIndex) {
            subListRangeCheck(fromIndex, toIndex, size);
            return new SubList(this, fromIndex, toIndex);
        }

        private void rangeCheckForAdd(int index) {
            if (index < 0 || index > this.size) {
                throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
            }
        }

        private String outOfBoundsMsg(int index) {
            return "Index: " + index + ", Size: " + this.size;
        }

        private void checkForComodification() {
            if (root.modCount != modCount) {
                throw new ConcurrentModificationException();
            }
        }

        private void updateSizeAndModCount(int sizeChange) {
            SubList slist = this;
            do {
                slist.size += sizeChange;
                slist.modCount = root.modCount;
                slist = slist.parent;
            } while (slist != null);
        }

        public Spliterator.OfDouble spliterator() {
            checkForComodification();

            return new Spliterator.OfDouble() {
                private int index = offset; // current index, modified on advance/split
                private int fence = -1; // -1 until used; then one past last index
                private int expectedModCount; // initialized when fence set

                private int getFence() { // initialize fence to size on first use
                    int hi; // (a specialized variant appears in method forEach)
                    if ((hi = fence) < 0) {
                        expectedModCount = modCount;
                        hi = fence = offset + size;
                    }
                    return hi;
                }

                public Spliterator.OfDouble trySplit() {
                    int hi = getFence(), lo = index, mid = (lo + hi) >>> 1;
                    // ArrayListSpliterator can be used here as the source is
                    // already bound
                    return (lo >= mid) ? null : // divide range in half unless too small
                    root.new ArrayListSpliterator(lo, index = mid, expectedModCount);
                }

                public boolean tryAdvance(DoubleConsumer action) {
                    Objects.requireNonNull(action);
                    int hi = getFence(), i = index;
                    if (i < hi) {
                        index = i + 1;
                        action.accept(root.elementData[i]);
                        if (root.modCount != expectedModCount) {
                            throw new ConcurrentModificationException();
                        }
                        return true;
                    }
                    return false;
                }

                public void forEachRemaining(DoubleConsumer action) {
                    Objects.requireNonNull(action);
                    int i, hi, mc; // hoist accesses and checks from loop
                    DoubleArrayList lst = root;
                    double[] a;
                    if ((a = lst.elementData) != null) {
                        if ((hi = fence) < 0) {
                            mc = modCount;
                            hi = offset + size;
                        } else {
                            mc = expectedModCount;
                        }
                        if ((i = index) >= 0 && (index = hi) <= a.length) {
                            for (; i < hi; ++i) {
                                action.accept(a[i]);
                            }
                            if (lst.modCount == mc) {
                                return;
                            }
                        }
                    }
                    throw new ConcurrentModificationException();
                }

                public long estimateSize() {
                    return getFence() - index;
                }

                public int characteristics() {
                    return Spliterator.ORDERED | Spliterator.SIZED | Spliterator.SUBSIZED;
                }
            };
        }

        @Override
        public int offset() {
            return offset;
        }

        @Override
        public double[] getArrayUnsafe() {
            return root.getArrayUnsafe();
        }

        @Override
        public DoubleList assignConst(double val) {
            final int start = offset;
            final int end = start + size;
            Arrays.fill(root.elementData, start, end, val);
            checkForComodification();
            return this;
        }

        @Override
        public DoubleList round(int scale) {
            final double[] es = root.elementData;
            final int start = offset;
            final int end = start + size;
            for (int i = start; i < end; ++i) {
                es[i] = BigDecimal.valueOf(es[i]).setScale(scale, RoundingMode.HALF_EVEN).doubleValue();
            }
            checkForComodification();
            return this;
        }

        @Override
        public DoubleList plus(double val) {
            final double[] es = root.elementData;
            final int start = offset;
            final int end = start + size;
            for (int i = start; i < end; ++i) {
                es[i] += val;
            }
            checkForComodification();
            return this;
        }

        @Override
        public DoubleList mul(double val) {
            final double[] es = root.elementData;
            final int start = offset;
            final int end = start + size;
            for (int i = start; i < end; ++i) {
                es[i] *= val;
            }
            checkForComodification();
            return this;
        }

        @Override
        public double min() {
            if (isEmpty()) {
                throw new NoSuchElementException();
            }
            return DoubleArrayList.min(size, offset, root.elementData);
        }

        @Override
        public double max() {
            if (isEmpty()) {
                throw new NoSuchElementException();
            }
            return DoubleArrayList.max(size, offset, root.elementData);
        }

        @Override
        public double avg() {
            if (isEmpty()) {
                throw new NoSuchElementException();
            }
            return DoubleArrayList.avg(size, offset, root.elementData);
        }

        @Override
        public double stddev() {
            return DoubleArrayList.stddev(size, offset, root.elementData);
        }

        @Override
        public double median() {
            if (isEmpty()) {
                throw new NoSuchElementException();
            }
            return DoubleArrayList.median(size, offset, root.elementData);
        }

        @Override
        public double iqr() {
            int length = checkLengthGeq2(size());
            double[] sorted = DoubleArrayList.sorted(length, offset, root.elementData);
            if (length == 2) {
                return sorted[1] - sorted[0];
            }
            int size = length / 2;
            double q3;
            if (length % 2 == 0) {
                q3 = DoubleArrayList.medianOnSorted(size, size, sorted);
            } else {
                q3 = DoubleArrayList.medianOnSorted(size, size + 1, sorted);
            }
            return q3 - DoubleArrayList.medianOnSorted(size, 0, sorted);
        }

        @Override
        public double[] getOutlierFences() {
            double[] sorted = sorted(checkLengthGeq2(size()), offset, root.elementData);
            double q1 = lowerQuartile(sorted);
            double q3 = upperQuartile(sorted);
            double iqr = q3 - q1;
            double lowerF = q1 - 1.5 * iqr;
            double upperF = q3 + 1.5 * iqr;
            double[] minmax = minmax(size, offset, root.elementData);
            lowerF = Math.max(minmax[0], lowerF);
            upperF = Math.min(minmax[1], upperF);
            return new double[] { lowerF, upperF };
        }

        @Override
        public double lowerQuartile() {
            return lowerQuartile(sorted(checkLengthGeq2(size()), offset, root.elementData));
        }

        private double lowerQuartile(double[] sorted) {
            int length = sorted.length;
            if (length == 2) {
                return sorted[0];
            }
            return DoubleArrayList.medianOnSorted(length / 2, 0, sorted);
        }

        @Override
        public double upperQuartile() {
            return upperQuartile(sorted(checkLengthGeq2(size()), offset, root.elementData));
        }

        private double upperQuartile(double[] sorted) {
            int length = sorted.length;
            if (length == 2) {
                return sorted[1];
            }
            int size = length / 2;
            if (length % 2 == 0) {
                return DoubleArrayList.medianOnSorted(size, size, sorted);
            } else {
                return DoubleArrayList.medianOnSorted(size, size + 1, sorted);
            }
        }

        @Override
        public double sum() {
            return DoubleArrayList.sum(size, offset, root.elementData);
        }

        @Override
        public double plusi(int index, double val) {
            checkIndex(index, size);
            final double[] es = root.elementData;
            return (es[offset + index] += val);
        }

        @Override
        public double muli(int index, double val) {
            checkIndex(index, size);
            final double[] es = root.elementData;
            return (es[offset + index] *= val);
        }

        @Override
        public double dot(DoubleList list) {
            if (size != Objects.requireNonNull(list, "list").size()) {
                throw new IllegalArgumentException("incommensurate vector dimensions");
            }
            return dot(size, offset, root.elementData, list.offset(), list.getArrayUnsafe());
        }

        @Override
        public double dot(double[] array) {
            if (size != Objects.requireNonNull(array, "array").length) {
                throw new IllegalArgumentException("incommensurate vector dimensions");
            }
            return dot(size, offset, root.elementData, 0, array);
        }

        private static double dot(int length, int aoff, double[] a, int boff, double[] b) {
            double product = 0.0;
            for (int i = 0; i < length; ++i) {
                product += a[aoff + i] * b[boff + i];
            }
            return product;
        }

        @Override
        public DoubleList cross(DoubleList list) {
            checkDim3(this, false);
            checkDim3(list, true);
            return new DoubleArrayList(
                    DoubleArrayList.cross(offset, root.elementData, list.offset(), list.getArrayUnsafe()), false);
        }

        @Override
        public DoubleList cross(double[] array) {
            checkDim3(this, false);
            checkDim3(array);
            return new DoubleArrayList(DoubleArrayList.cross(offset, root.elementData, 0, array), false);
        }

        @Override
        public double norm2() {
            return Math.sqrt(dot(size, offset, root.elementData, offset, root.elementData));
        }

        @Override
        public DoubleList softmax() {
            return normalizedSoftmax(1.0);
        }

        @Override
        public DoubleList normalizedSoftmax(double normalizedSum) {
            if (isEmpty()) {
                return new DoubleArrayList(0);
            }
            return new DoubleArrayList(
                    DoubleArrayList.normalizedSoftmax(normalizedSum, size, offset, root.elementData, new double[size]),
                    false);
        }

        @Override
        public double logSumExp() {
            if (isEmpty()) {
                throw new NoSuchElementException();
            }
            return DoubleArrayList.logSumExp(size, offset, root.elementData);
        }

        @Override
        public DoubleList plusn(DoubleList list) {
            int length = Math.min(size, Objects.requireNonNull(list, "list").size());
            plusn(length, offset, root.elementData, list.offset(), list.getArrayUnsafe());
            return this;
        }

        @Override
        public DoubleList plusn(double[] array) {
            int length = Math.min(size, Objects.requireNonNull(array, "array").length);
            plusn(length, offset, root.elementData, 0, array);
            return this;
        }

        @Override
        public DoubleList minusn(DoubleList list) {
            int length = Math.min(size, Objects.requireNonNull(list, "list").size());
            minusn(length, offset, root.elementData, list.offset(), list.getArrayUnsafe());
            return this;
        }

        @Override
        public DoubleList minusn(double[] array) {
            int length = Math.min(size, Objects.requireNonNull(array, "array").length);
            minusn(length, offset, root.elementData, 0, array);
            return this;
        }

        private static void plusn(int length, int aoff, double[] a, int boff, double[] b) {
            for (int i = 0; i < length; ++i) {
                a[aoff + i] += b[boff + i];
            }
        }

        private static void minusn(int length, int aoff, double[] a, int boff, double[] b) {
            for (int i = 0; i < length; ++i) {
                a[aoff + i] -= b[boff + i];
            }
        }

        @Override
        public DoubleList muln(DoubleList list) {
            int length = Math.min(size, Objects.requireNonNull(list, "list").size());
            muln(length, offset, root.elementData, list.offset(), list.getArrayUnsafe());
            return this;
        }

        @Override
        public DoubleList muln(double[] array) {
            int length = Math.min(size, Objects.requireNonNull(array, "array").length);
            muln(length, offset, root.elementData, 0, array);
            return this;
        }

        private static void muln(int length, int aoff, double[] a, int boff, double[] b) {
            for (int i = 0; i < length; ++i) {
                a[aoff + i] *= b[boff + i];
            }
        }

        @Override
        public DoubleList sanitizeNonFinite(double nanSurrogate, double posInfSurrogate, double negInfSurrogate) {
            DoubleArrayList.sanitize(size, offset, root.elementData, nanSurrogate, posInfSurrogate, negInfSurrogate);
            return this;
        }

        @Override
        public String toString() {
            StringBuilder buf = new StringBuilder();
            final double[] es = root.elementData;
            final int start = offset;
            final int end = start + size;
            for (int i = start; i < end; ++i) {
                buf.append(String.format(FORMAT_D, es[i]));
                if (i < end - 1) {
                    buf.append(", ");
                }
            }
            return buf.toString();
        }

        @Override
        public void forEach(DoubleConsumer action) {
            Objects.requireNonNull(action);
            final int expectedModCount = modCount;
            final double[] es = root.elementData;
            final int start = offset;
            final int end = start + size;
            for (int i = start; modCount == expectedModCount && i < end; ++i) {
                action.accept(es[i]);
            }
            checkForComodification();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void forEach(DoubleConsumer action) {
        Objects.requireNonNull(action);
        final int expectedModCount = modCount;
        final double[] es = elementData;
        final int size = this.size;
        for (int i = 0; modCount == expectedModCount && i < size; i++) { // "lgtm[java/constant-comparison]"
            action.accept(elementAt(es, i));
        }
        if (modCount != expectedModCount) {
            throw new ConcurrentModificationException();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Spliterator.OfDouble spliterator() {
        return new ArrayListSpliterator(0, -1, 0);
    }

    /** Index-based split-by-two, lazily initialized Spliterator */
    final class ArrayListSpliterator implements Spliterator.OfDouble {

        private int index; // current index, modified on advance/split
        private int fence; // -1 until used; then one past last index
        private int expectedModCount; // initialized when fence set

        /** Creates new spliterator covering the given range. */
        ArrayListSpliterator(int origin, int fence, int expectedModCount) {
            this.index = origin;
            this.fence = fence;
            this.expectedModCount = expectedModCount;
        }

        private int getFence() { // initialize fence to size on first use
            int hi; // (a specialized variant appears in method forEach)
            if ((hi = fence) < 0) {
                expectedModCount = modCount;
                hi = fence = size;
            }
            return hi;
        }

        public ArrayListSpliterator trySplit() {
            int hi = getFence(), lo = index, mid = (lo + hi) >>> 1;
            return (lo >= mid) ? null : // divide range in half unless too small
                    new ArrayListSpliterator(lo, index = mid, expectedModCount);
        }

        public boolean tryAdvance(DoubleConsumer action) {
            if (action == null) {
                throw new NullPointerException();
            }
            int hi = getFence(), i = index;
            if (i < hi) {
                index = i + 1;
                action.accept(elementData[i]);
                if (modCount != expectedModCount) {
                    throw new ConcurrentModificationException();
                }
                return true;
            }
            return false;
        }

        public void forEachRemaining(DoubleConsumer action) {
            int i, hi, mc; // hoist accesses and checks from loop
            double[] a;
            if (action == null) {
                throw new NullPointerException();
            }
            if ((a = elementData) != null) {
                if ((hi = fence) < 0) {
                    mc = modCount;
                    hi = size;
                } else {
                    mc = expectedModCount;
                }
                if ((i = index) >= 0 && (index = hi) <= a.length) {
                    for (; i < hi; ++i) {
                        action.accept(a[i]);
                    }
                    if (modCount == mc) {
                        return;
                    }
                }
            }
            throw new ConcurrentModificationException();
        }

        public long estimateSize() {
            return getFence() - index;
        }

        public int characteristics() {
            return Spliterator.ORDERED | Spliterator.SIZED | Spliterator.SUBSIZED;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void sort() {
        final int expectedModCount = modCount;
        Arrays.sort(elementData, 0, size);
        if (modCount != expectedModCount) {
            throw new ConcurrentModificationException();
        }
        modCount++;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int offset() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getArrayUnsafe() {
        double[] es = elementData;
        if (es == DEFAULTCAPACITY_EMPTY_ELEMENTDATA || es == EMPTY_ELEMENTDATA) {
            es = new double[0];
            elementData = es;
        }
        return es;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList assignConst(double val) {
        final int expectedModCount = modCount;
        Arrays.fill(elementData, 0, size, val);
        if (modCount == expectedModCount) { // "lgtm[java/constant-comparison]"
            return this;
        }
        throw new ConcurrentModificationException();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList round(int scale) {
        final int expectedModCount = modCount;
        final double[] es = elementData;
        final int end = size;
        for (int i = 0; i < end; ++i) {
            es[i] = BigDecimal.valueOf(es[i]).setScale(scale, RoundingMode.HALF_EVEN).doubleValue();
        }
        if (modCount == expectedModCount) { // "lgtm[java/constant-comparison]"
            return this;
        }
        throw new ConcurrentModificationException();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList plus(double val) {
        final int expectedModCount = modCount;
        final double[] es = elementData;
        final int end = size;
        for (int i = 0; i < end; ++i) {
            es[i] += val;
        }
        if (modCount == expectedModCount) { // "lgtm[java/constant-comparison]"
            return this;
        }
        throw new ConcurrentModificationException();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList mul(double val) {
        final int expectedModCount = modCount;
        final double[] es = elementData;
        final int end = size;
        for (int i = 0; i < end; ++i) {
            es[i] *= val;
        }
        if (modCount == expectedModCount) { // "lgtm[java/constant-comparison]"
            return this;
        }
        throw new ConcurrentModificationException();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double sum() {
        return sum(size, 0, elementData);
    }

    static double sum(int length, int aoff, double[] a) {
        double sum = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            sum += a[i];
        }
        return sum;
    }

    /* check for list length greater or equal 2 */
    static int checkLengthGeq2(int length) {
        if (length < 2) {
            throw new IllegalArgumentException("length is : " + length);
        }
        return length;
    }

    static double stddev(int length, int aoff, double[] a) {
        length = checkLengthGeq2(length);
        double sum = 0.0;
        double sumSqr = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            double x = a[i];
            sum += x;
            sumSqr += x * x;
        }
        sumSqr = sumSqr / length;
        double mean = sum / length;
        double var = sumSqr - (mean * mean);
        return Math.sqrt(var);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double plusi(int index, double val) {
        checkIndex(index, size);
        final double[] es = elementData;
        return (es[index] += val);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double muli(int index, double val) {
        checkIndex(index, size);
        final double[] es = elementData;
        return (es[index] *= val);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double dot(DoubleList list) {
        if (size != Objects.requireNonNull(list, "list").size()) {
            throw new IllegalArgumentException("incommensurate vector dimensions");
        }
        return dot(size, elementData, list.offset(), list.getArrayUnsafe());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double dot(double[] array) {
        if (size != Objects.requireNonNull(array, "array").length) {
            throw new IllegalArgumentException("incommensurate vector dimensions");
        }
        return dot(size, elementData, 0, array);
    }

    private static double dot(int length, double[] a, int boff, double[] b) {
        double product = 0.0;
        for (int i = 0; i < length; ++i) {
            product += a[i] * b[boff + i];
        }
        return product;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList cross(DoubleList list) {
        checkDim3(this, false);
        checkDim3(list, true);
        return new DoubleArrayList(cross(0, elementData, list.offset(), list.getArrayUnsafe()), false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList cross(double[] array) {
        checkDim3(this, false);
        checkDim3(array);
        return new DoubleArrayList(cross(0, elementData, 0, array), false);
    }

    static void checkDim3(DoubleList list, boolean isArg) {
        if (Objects.requireNonNull(list, "list").size() != 3) {
            String msg = "vector dimension needs to be 3 but is ";
            if (isArg) {
                throw new IllegalArgumentException(msg + list.size());
            } else {
                throw new UnsupportedOperationException(msg + list.size());
            }
        }
    }

    static void checkDim3(double[] array) {
        if (Objects.requireNonNull(array, "array").length != 3) {
            throw new IllegalArgumentException("vector dimension needs to be 3 but is " + array.length);
        }
    }

    static double[] cross(int aoff, double[] a, int boff, double[] b) {
        double a1 = a[aoff];
        double a2 = a[aoff + 1];
        double a3 = a[aoff + 2];
        double b1 = b[boff];
        double b2 = b[boff + 1];
        double b3 = b[boff + 2];
        double[] c = new double[3];
        c[0] = a2 * b3 - a3 * b2;
        c[1] = a3 * b1 - a1 * b3;
        c[2] = a1 * b2 - a2 * b1;
        return c;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double norm2() {
        return Math.sqrt(dot(size, elementData, 0, elementData));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList softmax() {
        return normalizedSoftmax(1.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList normalizedSoftmax(double normalizedSum) {
        if (isEmpty()) {
            return new DoubleArrayList(0);
        }
        return new DoubleArrayList(normalizedSoftmax(normalizedSum, size, 0, elementData, new double[size]), false);
    }

    static double[] normalizedSoftmax(double normalizedSum, int length, int aoff, double[] a, double[] out) {
        double max = max(length, aoff, a);
        double s = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            double q = Math.exp(a[i] - max);
            s += q;
            out[i - aoff] = q;
        }
        s = normalizedSum / s;
        for (int i = 0; i < length; ++i) {
            out[i] *= s;
        }
        return out;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double logSumExp() {
        if (isEmpty()) {
            throw new NoSuchElementException();
        }
        return logSumExp(size, 0, elementData);
    }

    static double logSumExp(int length, int aoff, double[] a) {
        double max = max(length, aoff, a);
        double s = 0.0;
        for (int i = aoff; i < aoff + length; ++i) {
            s += Math.exp(a[i] - max);
        }
        return max + Math.log(s);
    }

    static double max(int length, int aoff, double[] a) {
        double max = a[aoff];
        for (int i = aoff + 1; i < aoff + length; ++i) {
            max = Math.max(max, a[i]);
        }
        return max;
    }

    static double min(int length, int aoff, double[] a) {
        double min = a[aoff];
        for (int i = aoff + 1; i < aoff + length; ++i) {
            min = Math.min(min, a[i]);
        }
        return min;
    }

    static double[] minmax(int length, int aoff, double[] a) {
        double min = a[aoff];
        double max = min;
        for (int i = aoff + 1; i < aoff + length; ++i) {
            double x = a[i];
            min = Math.min(min, x);
            max = Math.max(max, x);
        }
        return new double[] { min, max };
    }

    static double avg(int length, int aoff, double[] a) {
        return sum(length, aoff, a) / length;
    }

    static double median(int length, int aoff, double[] a) {
        if (length == 1) {
            return a[aoff];
        }
        return medianOnSorted(length, aoff, sorted(length, aoff, a));
    }

    static double medianOnSorted(int length, int aoff, double[] a) {
        if (length == 1) {
            return a[aoff];
        }
        int mid = aoff + length / 2;
        if (length % 2 == 0) {
            return (a[mid - 1] + a[mid]) / 2.0;
        }
        return a[mid];
    }

    static double[] sorted(int length, int aoff, double[] a) {
        double[] b = new double[length];
        System.arraycopy(a, aoff, b, 0, length);
        Arrays.sort(b);
        return b;
    }

    static void shuffle(int length, int aoff, double[] a) {
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        for (int i = length; i > 1; --i) {
            swap(aoff, a, i - 1, rnd.nextInt(i));
        }
    }

    static void swap(int aoff, double[] a, int i, int j) {
        int ii = aoff + i;
        int jj = aoff + j;
        double tmp = a[ii];
        a[ii] = a[jj];
        a[jj] = tmp;
    }

    static void sanitize(int length, int aoff, double[] a, double nanSurrogate, double posInfSurrogate,
            double negInfSurrogate) {
        boolean subNan = (nanSurrogate == nanSurrogate); // "lgtm[java/comparison-of-identical-expressions]"
        boolean subPInf = (posInfSurrogate != Double.POSITIVE_INFINITY);
        boolean subNInf = (negInfSurrogate != Double.NEGATIVE_INFINITY);
        if (!subNan && !subPInf && !subNInf) {
            return;
        }
        for (int i = aoff; i < aoff + length; ++i) {
            double x = a[i];
            if (x != x && subNan) { // "lgtm[java/comparison-of-identical-expressions]"
                a[i] = nanSurrogate;
            } else if (x == Double.POSITIVE_INFINITY && subPInf) {
                a[i] = posInfSurrogate;
            } else if (x == Double.NEGATIVE_INFINITY && subNInf) {
                a[i] = negInfSurrogate;
            }
        }
    }

    static DoubleList filter(int length, int aoff, double[] a, DoublePredicate predicate) {
        predicate = Objects.requireNonNull(predicate, "predicate");
        DoubleArrayList result = new DoubleArrayList();
        for (int i = aoff; i < aoff + length; ++i) {
            double x = a[i];
            if (predicate.test(x)) {
                result.add(x);
            }
        }
        return result;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double min() {
        if (isEmpty()) {
            throw new NoSuchElementException();
        }
        return min(size, 0, elementData);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double max() {
        if (isEmpty()) {
            throw new NoSuchElementException();
        }
        return max(size, 0, elementData);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double avg() {
        if (isEmpty()) {
            throw new NoSuchElementException();
        }
        return avg(size, 0, elementData);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double stddev() {
        return stddev(size, 0, elementData);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double median() {
        if (isEmpty()) {
            throw new NoSuchElementException();
        }
        return median(size, 0, elementData);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double iqr() {
        int length = checkLengthGeq2(size());
        double[] sorted = sorted(length, 0, elementData);
        if (length == 2) {
            return sorted[1] - sorted[0];
        }
        int size = length / 2;
        double q3;
        if (length % 2 == 0) {
            q3 = medianOnSorted(size, size, sorted);
        } else {
            q3 = medianOnSorted(size, size + 1, sorted);
        }
        return q3 - medianOnSorted(size, 0, sorted);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getOutlierFences() {
        double[] sorted = sorted(checkLengthGeq2(size()), 0, elementData);
        double q1 = lowerQuartile(sorted);
        double q3 = upperQuartile(sorted);
        double iqr = q3 - q1;
        double lowerF = q1 - 1.5 * iqr;
        double upperF = q3 + 1.5 * iqr;
        double[] minmax = minmax(size, 0, elementData);
        lowerF = Math.max(minmax[0], lowerF);
        upperF = Math.min(minmax[1], upperF);
        return new double[] { lowerF, upperF };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double lowerQuartile() {
        return lowerQuartile(sorted(checkLengthGeq2(size()), 0, elementData));
    }

    private double lowerQuartile(double[] sorted) {
        int length = sorted.length;
        if (length == 2) {
            return sorted[0];
        }
        return medianOnSorted(length / 2, 0, sorted);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double upperQuartile() {
        return upperQuartile(sorted(checkLengthGeq2(size()), 0, elementData));
    }

    private double upperQuartile(double[] sorted) {
        int length = sorted.length;
        if (length == 2) {
            return sorted[1];
        }
        int size = length / 2;
        if (length % 2 == 0) {
            return medianOnSorted(size, size, sorted);
        } else {
            return medianOnSorted(size, size + 1, sorted);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList plusn(DoubleList list) {
        int length = Math.min(size, Objects.requireNonNull(list, "list").size());
        plusn(length, elementData, list.offset(), list.getArrayUnsafe());
        return this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList plusn(double[] array) {
        int length = Math.min(size, Objects.requireNonNull(array, "array").length);
        plusn(length, elementData, 0, array);
        return this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList minusn(DoubleList list) {
        int length = Math.min(size, Objects.requireNonNull(list, "list").size());
        minusn(length, elementData, list.offset(), list.getArrayUnsafe());
        return this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList minusn(double[] array) {
        int length = Math.min(size, Objects.requireNonNull(array, "array").length);
        minusn(length, elementData, 0, array);
        return this;
    }

    private static void plusn(int length, double[] a, int boff, double[] b) {
        for (int i = 0; i < length; ++i) {
            a[i] += b[boff + i];
        }
    }

    private static void minusn(int length, double[] a, int boff, double[] b) {
        for (int i = 0; i < length; ++i) {
            a[i] -= b[boff + i];
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList muln(DoubleList list) {
        int length = Math.min(size, Objects.requireNonNull(list, "list").size());
        muln(length, elementData, list.offset(), list.getArrayUnsafe());
        return this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList muln(double[] array) {
        int length = Math.min(size, Objects.requireNonNull(array, "array").length);
        muln(length, elementData, 0, array);
        return this;
    }

    private static void muln(int length, double[] a, int boff, double[] b) {
        for (int i = 0; i < length; ++i) {
            a[i] *= b[boff + i];
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public DoubleList sanitizeNonFinite(double nanSurrogate, double posInfSurrogate, double negInfSurrogate) {
        sanitize(size, 0, elementData, nanSurrogate, posInfSurrogate, negInfSurrogate);
        return this;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        final double[] es = elementData;
        final int end = size;
        for (int i = 0; i < end; ++i) {
            buf.append(String.format(FORMAT_D, es[i]));
            if (i < end - 1) {
                buf.append(", ");
            }
        }
        return buf.toString();
    }

    static void checkIndex(int index, int length) {
        if (index < 0 || index >= length) {
            throw new IndexOutOfBoundsException(String.format("Index %d out-of-bounds for length %d", index, length));
        }
    }

    static void subListRangeCheck(int fromIndex, int toIndex, int size) {
        if (fromIndex < 0) {
            throw new IndexOutOfBoundsException("fromIndex = " + fromIndex);
        }
        if (toIndex > size) {
            throw new IndexOutOfBoundsException("toIndex = " + toIndex);
        }
        if (fromIndex > toIndex) {
            throw new IllegalArgumentException("fromIndex(" + fromIndex + ") > toIndex(" + toIndex + ")");
        }
    }

    static boolean checkApproxEqualArgs(DoubleList a, DoubleList b, double relTol, double absTol) {
        if (relTol < 0.0 || Double.isNaN(relTol) || Double.isInfinite(relTol)) {
            throw new IllegalArgumentException("illegal relTol : " + relTol);
        }
        if (absTol < 0.0 || Double.isNaN(absTol) || Double.isInfinite(absTol)) {
            throw new IllegalArgumentException("illegal absTol : " + absTol);
        }
        if (a.size() != b.size()) {
            return false;
        }
        return true;
    }

    private static final String FORMAT_D = "%.12E";
    private static final long serialVersionUID = 1L;

    @Override
    public void writeExternal(ObjectOutput out) throws IOException {
        int expectedModCount = modCount;
        out.writeInt(size);
        out.writeObject(elementData);
        if (modCount != expectedModCount) {
            throw new ConcurrentModificationException();
        }
    }

    @Override
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        size = in.readInt();
        elementData = (double[]) in.readObject();
    }
}
