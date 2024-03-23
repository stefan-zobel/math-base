/*
 * Copyright (C) 2003 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math;

/**
 * A class of static utility functions for manipulating arrays of double.
 */
public final class VectorOps {
    /**
     * Sets every element of a double array to a given value.
     * 
     * @param m
     *            The array to modify
     * @param v
     *            The value
     */
    public static void setAll(double[] m, double v) {
        java.util.Arrays.fill(m, v);
    }

    public static void set(double[] dest, double[] source) {
        if (source.length != dest.length) {
            throw new IllegalArgumentException("source.length != dest.length");
        }
        System.arraycopy(source, 0, dest, 0, source.length);
    }

    /**
     * Multiplies every element in an array by a scalar.
     * 
     * @param m
     *            The array
     * @param factor
     *            The scalar
     */
    public static void timesEquals(double[] m, double factor) {
        for (int i = 0; i < m.length; i++) {
            m[i] *= factor;
        }
    }

    public static void plusEquals(double[] m1, double[] m2) {
        for (int i = 0; i < m1.length; i++) {
            if (Double.isInfinite(m1[i]) && Double.isInfinite(m2[i]) && (m1[i] * m2[i] < 0.0)) {
                m1[i] = 0.0;
            } else {
                m1[i] += m2[i];
            }
        }
    }

    public static void plusEquals(double[] m1, double[] m2, double factor) {
        for (int i = 0; i < m1.length; i++) {
            double m1i = m1[i];
            double m2i = m2[i];
            if (Double.isInfinite(m1i) && Double.isInfinite(m2i) && (m1[i] * m2[i] < 0.0)) {
                m1[i] = 0.0;
            } else {
                m1[i] += m2[i] * factor;
            }
        }
    }

    public static double dotProduct(double[] m1, double[] m2) {
        double ret = 0.0;
        for (int i = 0; i < m1.length; i++) {
            ret += m1[i] * m2[i];
        }
        return ret;
    }

    public static double absNorm(double[] m) {
        double ret = 0.0;
        for (int i = 0; i < m.length; i++) {
            ret += Math.abs(m[i]);
        }
        return ret;
    }

    public static double twoNorm(double[] m) {
        double ret = 0.0;
        for (int i = 0; i < m.length; i++) {
            ret += m[i] * m[i];
        }
        return Math.sqrt(ret);
    }

    public static double oneNorm(double[] m) {
        double ret = 0.0;
        for (int i = 0; i < m.length; i++) {
            ret += m[i];
        }
        return ret;
    }

    public static double infinityNorm(double[] m) {
        double ret = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < m.length; i++) {
            if (Math.abs(m[i]) > ret) {
                ret = Math.abs(m[i]);
            }
        }
        return ret;
    }

    public static double absNormalize(double[] m) {
        double norm = absNorm(m);
        if (norm > 0.0) {
            for (int i = 0; i < m.length; i++) {
                m[i] /= norm;
            }
        }
        return norm;
    }

    public static boolean isNaN(double[] m) {
        for (int i = 0; i < m.length; i++) {
            if (Double.isNaN(m[i])) {
                return true;
            }
        }
        return false;
    }

    public static double mean(double[] m) {
        double sum = 0.0;
        for (int i = 0; i < m.length; i++) {
            sum += m[i];
        }
        return sum / m.length;
    }

    public static double max(double[] elems) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < elems.length; i++) {
            double elem = elems[i];
            if (elem > max) {
                max = elem;
            }
        }
        return max;
    }

    public static double min(double[] elems) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < elems.length; i++) {
            double elem = elems[i];
            if (elem < min) {
                min = elem;
            }
        }
        return min;
    }
}
