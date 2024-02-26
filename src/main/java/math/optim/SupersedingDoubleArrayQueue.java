/*
 * Copyright (C) 2016 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

final class SupersedingDoubleArrayQueue {

    private int size;
    private final double[][] values;

    SupersedingDoubleArrayQueue(int capacity) {
        values = new double[capacity][];
    }

    int size() {
        return size;
    }

    double[] get(int index) {
        return values[index];
    }

    double[] removeFirst() {
        double[] fst = values[0];
        System.arraycopy(values, 1, values, 0, size - 1);
        --size;
        return fst;
    }

    void addLast(double[] value) {
        if (size < values.length) {
            double[] copy = new double[value.length];
            System.arraycopy(value, 0, copy, 0, copy.length);
            values[size] = copy;
            ++size;
        } else {
            double[] fst = values[0];
            System.arraycopy(value, 0, fst, 0, value.length);
            System.arraycopy(values, 1, values, 0, size - 1);
            values[values.length - 1] = fst;
        }
    }

    void addLastNoCopy(double[] value) {
        if (size < values.length) {
            values[size] = value;
            ++size;
        } else {
            System.arraycopy(values, 1, values, 0, size - 1);
            values[values.length - 1] = value;
        }
    }
}
