/*
 * Copyright (C) 2003 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.linalg;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import jdk.incubator.vector.VectorOperators;

/**
 * A class of static utility functions for manipulating arrays of double.
 */
public final class VectorOps {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    /**
     * {@code sqrt(Double.MAX_VALUE)}. Squares of magnitudes at or below this
     * cannot overflow, so {@link #twoNorm(double[])} may sum them directly.
     */
    private static final double SQRT_MAX_VALUE = Math.sqrt(Double.MAX_VALUE);

    /**
     * {@code sqrt(Double.MIN_NORMAL)}. Squares of magnitudes at or above this
     * stay normal, so summing them directly loses no precision to gradual
     * underflow.
     */
    private static final double SQRT_MIN_NORMAL = Math.sqrt(Double.MIN_NORMAL);

    public static boolean isVectorized() {
        return true;
    }

    /**
     * Sets every element of a double array to a given value.
     * 
     * @param m The array to modify
     * @param v The value
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
     * @param m      The array
     * @param factor The scalar
     */
    public static void timesEquals(double[] m, double factor) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m.length);

        for (; i < upperBound; i += SPECIES.length()) {
            var va = DoubleVector.fromArray(SPECIES, m, i);
            va.mul(factor).intoArray(m, i);
        }
        for (; i < m.length; i++) {
            m[i] *= factor;
        }
    }

    public static void plusEquals(double[] m1, double[] m2) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m1.length);

        for (; i < upperBound; i += SPECIES.length()) {
            var v1 = DoubleVector.fromArray(SPECIES, m1, i);
            var v2 = DoubleVector.fromArray(SPECIES, m2, i);

            VectorMask<Double> m1Inf = v1.test(VectorOperators.IS_INFINITE);
            VectorMask<Double> m2Inf = v2.test(VectorOperators.IS_INFINITE);
            VectorMask<Double> bothInf = m1Inf.and(m2Inf);

            VectorMask<Double> diffSign = v1.mul(v2).compare(VectorOperators.LT, 0.0);

            VectorMask<Double> specialCase = bothInf.and(diffSign);

            var res = v1.add(v2);

            res = res.blend(0.0, specialCase);

            res.intoArray(m1, i);
        }

        for (; i < m1.length; i++) {
            if (Double.isInfinite(m1[i]) && Double.isInfinite(m2[i]) && (m1[i] * m2[i] < 0.0)) {
                m1[i] = 0.0;
            } else {
                m1[i] += m2[i];
            }
        }
    }

    public static void plusEquals(double[] m1, double[] m2, double factor) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m1.length);

        var vFactor = DoubleVector.broadcast(SPECIES, factor);

        for (; i < upperBound; i += SPECIES.length()) {
            var v1 = DoubleVector.fromArray(SPECIES, m1, i);
            var v2 = DoubleVector.fromArray(SPECIES, m2, i);

            var specialCase = v1.test(VectorOperators.IS_INFINITE)
                                .and(v2.test(VectorOperators.IS_INFINITE))
                                .and(v1.mul(v2).compare(VectorOperators.LT, 0.0));

            var res = v2.fma(vFactor, v1);

            res = res.blend(0.0, specialCase);

            res.intoArray(m1, i);
        }

        for (; i < m1.length; i++) {
            if (Double.isInfinite(m1[i]) && Double.isInfinite(m2[i]) && (m1[i] * m2[i] < 0.0)) {
                m1[i] = 0.0;
            } else {
                m1[i] += m2[i] * factor;
            }
        }
    }

    /**
     * The inner product {@code sum m1_i * m2_i}, over as many elements as
     * {@code m1} holds.
     * <p>
     * The products are accumulated as they come, so unlike
     * {@link #twoNorm(double[])} this can overflow for operands whose product
     * leaves the normal range.
     *
     * @param m1 The first array, and the one whose length is used
     * @param m2 The second array, at least as long as {@code m1}
     * @return the inner product, {@code 0.0} for empty arrays
     */
    public static double dotProduct(double[] m1, double[] m2) {
        double ret = 0.0;
        int i = 0;
        int upperBound = SPECIES.loopBound(m1.length);
        var acc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            var va = DoubleVector.fromArray(SPECIES, m1, i);
            var vb = DoubleVector.fromArray(SPECIES, m2, i);
            acc = va.fma(vb, acc); // Fused Multiply-Add: (va * vb) + acc
        }

        ret = acc.reduceLanes(VectorOperators.ADD);

        for (; i < m1.length; i++) {
            ret += m1[i] * m2[i];
        }
        return ret;
    }

    /**
     * The one-norm {@code sum |m_i|}, the sum of the absolute values.
     *
     * @param m The array
     * @return the one-norm, {@code 0.0} for an empty array
     */
    public static double absNorm(double[] m) {
        double ret = 0.0;
        for (int i = 0; i < m.length; i++) {
            ret += Math.abs(m[i]);
        }
        return ret;
    }

    /**
     * The Euclidean norm {@code sqrt(sum m_i^2)}, computed so that it neither
     * overflows for a large vector nor underflows for a small one.
     * <p>
     * The magnitude of the largest element decides. Where its square can be
     * formed and summed without leaving the normal range -- which is every
     * vector this library is likely to see -- the squares are summed directly
     * and the answer is the one the straightforward loop gives. Otherwise the
     * elements are scaled by a power of two, which is exact, and the result is
     * scaled back. Accumulating the squares unscaled instead returns
     * {@code Infinity} from about {@code 8.4e152} upwards and {@code 0.0}
     * below about {@code 1.1e-162}, with the accuracy already degrading an
     * order of magnitude before that.
     *
     * @param m The array
     * @return the Euclidean norm, {@code 0.0} for an empty array,
     *         {@code Double.NaN} if any element is {@code NaN}
     */
    public static double twoNorm(double[] m) {
        double max = maxAbs(m);
        if (max == 0.0) {
            return 0.0;
        }
        if (max >= SQRT_MIN_NORMAL && max * Math.sqrt(m.length) <= SQRT_MAX_VALUE) {
            return Math.sqrt(sumOfSquares(m));
        }
        // NaN and infinite input reach this branch as well and carry through
        // it unchanged: the scaling is a multiplication by a power of two
        int exponent = Math.getExponent(max);
        double scale = Math.scalb(1.0, -exponent);
        return Math.scalb(Math.sqrt(sumOfScaledSquares(m, scale)), exponent);
    }

    /**
     * The largest absolute value in the array, {@code 0.0} for an empty one and
     * {@code NaN} if any element is {@code NaN}. A maximum does not depend on
     * the order it is taken in, so unlike the sums around it this is bit for
     * bit what the scalar tree computes.
     */
    private static double maxAbs(double[] m) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m.length);
        var acc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            acc = acc.lanewise(VectorOperators.MAX, DoubleVector.fromArray(SPECIES, m, i).abs());
        }

        double max = acc.reduceLanes(VectorOperators.MAX);

        for (; i < m.length; i++) {
            max = Math.max(max, Math.abs(m[i]));
        }
        return max;
    }

    /** {@code sum m_i^2}, the loop {@link #twoNorm(double[])} used to be. */
    private static double sumOfSquares(double[] m) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m.length);
        var acc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            var va = DoubleVector.fromArray(SPECIES, m, i);
            acc = va.fma(va, acc);
        }

        double sum = acc.reduceLanes(VectorOperators.ADD);

        for (; i < m.length; i++) {
            sum += m[i] * m[i];
        }
        return sum;
    }

    /** {@code sum (scale * m_i)^2}, for a {@code scale} that is a power of two. */
    private static double sumOfScaledSquares(double[] m, double scale) {
        int i = 0;
        int upperBound = SPECIES.loopBound(m.length);
        var acc = DoubleVector.zero(SPECIES);
        var vScale = DoubleVector.broadcast(SPECIES, scale);

        for (; i < upperBound; i += SPECIES.length()) {
            var va = DoubleVector.fromArray(SPECIES, m, i).mul(vScale);
            acc = va.fma(va, acc);
        }

        double sum = acc.reduceLanes(VectorOperators.ADD);

        for (; i < m.length; i++) {
            double scaled = m[i] * scale;
            sum += scaled * scaled;
        }
        return sum;
    }

    /**
     * The signed sum {@code sum m_i}.
     *
     * @param m The array
     * @return the sum of the elements, {@code 0.0} for an empty array
     * @deprecated the name says norm but the absolute values are not taken, so
     *             this returns {@code -1.0} for {@code (3, -4, 5, -5)} where
     *             the one-norm is {@code 17}. Use {@link #absNorm(double[])}
     *             for the one-norm, or {@link #sum(double[])} for the signed
     *             sum this actually computes.
     */
    @Deprecated
    public static double oneNorm(double[] m) {
        return sum(m);
    }

    /**
     * The signed sum {@code sum m_i}.
     * <p>
     * This can overflow to an infinity, and nothing can be done about that: a
     * sum that leaves the range has no representable answer.
     * {@link #mean(double[])} does repair the case where the mean itself is
     * still representable.
     *
     * @param m The array
     * @return the sum of the elements, {@code 0.0} for an empty array
     * @since 1.5.2
     */
    public static double sum(double[] m) {
        double ret = 0.0;
        for (int i = 0; i < m.length; i++) {
            ret += m[i];
        }
        return ret;
    }

    /**
     * The maximum norm {@code max |m_i|}.
     *
     * @param m The array
     * @return the maximum norm, {@code 0.0} for an empty array,
     *         {@code Double.NaN} if any element is {@code NaN}
     */
    public static double infinityNorm(double[] m) {
        return maxAbs(m);
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

    /**
     * The arithmetic mean {@code (sum m_i) / n}, computed so that it does not
     * overflow for a vector whose sum leaves the range while its mean does
     * not.
     * <p>
     * The elements are summed directly, which for every vector that has no
     * problem is bit for bit what the straightforward loop gives. Only where
     * that sum comes back infinite -- a thousand elements of {@code 1e306}
     * reach it, and so does {@code (1e308, 1e308, -1e308)}, whose overflow
     * cancels away again -- are the elements scaled by a power of two, which
     * is exact, and the result scaled back. {@link #sum(double[])} has no such
     * repair available to it: a sum that overflows has no representable
     * answer.
     * <p>
     * Both loops are scalar in the Java 8 and the Java 25 source tree, so
     * unlike {@link #twoNorm(double[])} and
     * {@link #dotProduct(double[], double[])} this method returns the same
     * bits in either.
     *
     * @param m The array
     * @return the arithmetic mean, {@code Double.NaN} for an empty array
     */
    public static double mean(double[] m) {
        double sum = 0.0;
        for (int i = 0; i < m.length; i++) {
            sum += m[i];
        }
        if (Math.abs(sum) <= Double.MAX_VALUE) {
            return sum / m.length;
        }
        // NaN and infinite input reach this branch as well and carry through
        // it unchanged: the scaling is a multiplication by a power of two
        int exponent = Math.getExponent(maxAbs(m));
        double scale = Math.scalb(1.0, -exponent);
        double scaled = 0.0;
        for (int i = 0; i < m.length; i++) {
            scaled += m[i] * scale;
        }
        return Math.scalb(scaled / m.length, exponent);
    }

    public static double max(double[] elems) {
        int i = 0;
        int upperBound = SPECIES.loopBound(elems.length);
        var res = DoubleVector.broadcast(SPECIES, Double.NEGATIVE_INFINITY);

        for (; i < upperBound; i += SPECIES.length()) {
            var va = DoubleVector.fromArray(SPECIES, elems, i);
            res = res.lanewise(VectorOperators.MAX, va);
        }

        double max = res.reduceLanes(VectorOperators.MAX);

        for (; i < elems.length; i++) {
            if (elems[i] > max) {
                max = elems[i];
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
