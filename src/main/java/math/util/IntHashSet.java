/*
 * Copyright 2022 Stefan Zobel
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
package math.util;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

/**
 * A set that holds primitive integer values. An {@code int} value that cannot
 * occur as a value in the set must be passed as {@code sentinel} on
 * construction. For efficiency, the {@code int} taking methods
 * {@link #addInt(int)}, {@link #containsInt(int)}, {@link #removeInt(int)}
 * should be used.
 */
public class IntHashSet implements Set<Integer>, Cloneable {

    private IntIntHashMap map;

    private final int PRESENT = Integer.MAX_VALUE;

    private final int sentinel;

    public IntHashSet(int sentinel) {
        this.map = new IntIntHashMap(sentinel);
        this.sentinel = sentinel;
    }

    public IntHashSet(Collection<Integer> c, int sentinel) {
        this.map = new IntIntHashMap(Math.max((int) (c.size() / .75f) + 1, 16), sentinel);
        this.sentinel = sentinel;
        addAll(c);
    }

    public IntHashSet(int initialCapacity, float loadFactor, int sentinel) {
        this.map = new IntIntHashMap(initialCapacity, loadFactor, sentinel);
        this.sentinel = sentinel;
    }

    public IntHashSet(int initialCapacity, int sentinel) {
        this.map = new IntIntHashMap(initialCapacity, sentinel);
        this.sentinel = sentinel;
    }

    @Override
    public Iterator<Integer> iterator() {
        return map.keySet().iterator();
    }

    @Override
    public int size() {
        return map.size();
    }

    @Override
    public boolean isEmpty() {
        return map.isEmpty();
    }

    public boolean containsInt(int value) {
        return map.containsKeyInt(value);
    }

    @Override
    public boolean contains(Object o) {
        return map.containsKey(o);
    }

    public boolean addInt(int value) {
        return map.putInt(value, PRESENT) == sentinel;
    }

    @Override
    public boolean add(Integer e) {
        return map.put(e, PRESENT) == null;
    }

    public boolean removeInt(int value) {
        return map.removeInt(value) != sentinel;
    }

    @Override
    public boolean remove(Object o) {
        return map.remove(o) == PRESENT;
    }

    @Override
    public void clear() {
        map.clear();
    }

    @Override
    public Object[] toArray() {
        Object[] r = new Object[size()];
        Iterator<Integer> it = iterator();
        for (int i = 0; i < r.length; i++) {
            if (!it.hasNext()) {
                return Arrays.copyOf(r, i);
            }
            r[i] = it.next();
        }
        return it.hasNext() ? finishToArray(r, it) : r;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(T[] a) {
        int size = size();
        T[] r = a.length >= size ? a : (T[]) java.lang.reflect.Array.newInstance(a.getClass().getComponentType(), size);
        Iterator<Integer> it = iterator();

        for (int i = 0; i < r.length; i++) {
            if (!it.hasNext()) {
                if (a == r) {
                    r[i] = null;
                } else if (a.length < i) {
                    return Arrays.copyOf(r, i);
                } else {
                    System.arraycopy(r, 0, a, 0, i);
                    if (a.length > i) {
                        a[i] = null;
                    }
                }
                return a;
            }
            r[i] = (T) it.next();
        }
        return it.hasNext() ? finishToArray(r, it) : r;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        for (Object e : c) {
            if (!contains(e)) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean addAll(Collection<? extends Integer> c) {
        boolean modified = false;
        for (Integer e : c) {
            if (add(e)) {
                modified = true;
            }
        }
        return modified;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        boolean modified = false;
        Iterator<Integer> it = iterator();
        while (it.hasNext()) {
            if (!c.contains(it.next())) {
                it.remove();
                modified = true;
            }
        }
        return modified;
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        boolean modified = false;
        if (size() > c.size()) {
            for (Iterator<?> i = c.iterator(); i.hasNext();) {
                modified |= remove(i.next());
            }
        } else {
            for (Iterator<?> i = iterator(); i.hasNext();) {
                if (c.contains(i.next())) {
                    i.remove();
                    modified = true;
                }
            }
        }
        return modified;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof Set)) {
            return false;
        }
        @SuppressWarnings("rawtypes")
        Collection c = (Collection) o;
        if (c.size() != size()) {
            return false;
        }
        try {
            return containsAll(c);
        } catch (ClassCastException unused) {
            return false;
        } catch (NullPointerException unused) {
            return false;
        }
    }

    @Override
    public int hashCode() {
        int h = 0;
        Iterator<Integer> it = iterator();
        while (it.hasNext()) {
            Integer i = it.next();
            if (i != null) {
                h += i.hashCode();
            }
        }
        return h;
    }

    @Override
    public String toString() {
        Iterator<Integer> it = iterator();
        if (!it.hasNext()) {
            return "[]";
        }
        StringBuilder sb = new StringBuilder();
        sb.append('[');
        for (;;) {
            Integer e = it.next();
            sb.append(e);
            if (!it.hasNext()) {
                return sb.append(']').toString();
            }
            sb.append(',').append(' ');
        }
    }

    @Override
    public Object clone() {
        try {
            IntHashSet newSet = (IntHashSet) super.clone();
            newSet.map = (IntIntHashMap) map.clone();
            return newSet;
        } catch (CloneNotSupportedException e) {
            throw new InternalError();
        }
    }

    private static final int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 8;

    @SuppressWarnings("unchecked")
    private static <T> T[] finishToArray(T[] r, Iterator<Integer> it) {
        int i = r.length;
        while (it.hasNext()) {
            int cap = r.length;
            if (i == cap) {
                int newCap = cap + (cap >> 1) + 1;
                if (newCap - MAX_ARRAY_SIZE > 0) {
                    newCap = hugeCapacity(cap + 1);
                }
                r = Arrays.copyOf(r, newCap);
            }
            r[i++] = (T) it.next();
        }
        return (i == r.length) ? r : Arrays.copyOf(r, i);
    }

    private static int hugeCapacity(int minCapacity) {
        if (minCapacity < 0) { // overflow
            throw new OutOfMemoryError("Required array size too large");
        }
        return (minCapacity > MAX_ARRAY_SIZE) ? Integer.MAX_VALUE : MAX_ARRAY_SIZE;
    }
}
