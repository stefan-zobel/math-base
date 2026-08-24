/*
 * Copyright 2018, 2026 Stefan Zobel
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
package math.coord;

import java.util.Arrays;
import java.util.NoSuchElementException;

import math.fun.DBiConsumer;
import math.fun.DConsumer;
import math.fun.DForEach;
import math.fun.DForEachBi;
import math.fun.DFunction;
import math.fun.DIndexIterator;

/**
 * Evenly spaced points between {@code start} and {@code end}. The spacing
 * between the points is {@code (end-start)/(n-1)} where {@code n} is the number
 * of points. A {@code LinSpace} always includes the endpoints. If {@code end}
 * is smaller than {@code start}, then the {@code LinSpace} describes descending
 * values.
 * <p>
 * The number of points is always the number that was asked for. Equal
 * endpoints are not a special case: {@code linspace(a, a, n)} is {@code n}
 * copies of {@code a} with a spacing of zero. A {@code LinSpace} of a single
 * point sits on {@code start}, and {@link #end()} reports that same point.
 * <p>
 * Every accessor reports the same abscissa for the same position:
 * {@link #point(int)}, {@link #points()}, {@link #iterator()},
 * {@link #forEach()}, {@link #forEachBi()} and {@link #eval(DFunction)} agree
 * bit for bit, and both endpoints are exact. A slice begins and ends exactly
 * on the points of the range it was cut from; its interior points are its own,
 * because it carries its own spacing.
 * <p>
 * Note that indexes are 1-based!
 */
public final class LinSpace {

    private final double start;
    private final double stop;
    private final int numberOfPoints;
    private final double step;
    private double[] vec;

    private LinSpace(double start, double end, int numberOfPoints, double[] data) {
        checkArg(start, "start");
        checkArg(end, "end");
        if (numberOfPoints <= 0) {
            throw new IllegalArgumentException("numberOfPoints must be strictly positive : " + numberOfPoints);
        }
        if (numberOfPoints == 1) {
            // one point has no interval, so it sits on start the way
            // numpy's linspace(a, b, 1) does. stop follows start rather than
            // end, otherwise pos == 1 and pos == numberOfPoints coincide and
            // point(1) would answer start where points() answers stop
            this.start = start;
            this.stop = start;
            this.numberOfPoints = 1;
        } else {
            this.start = start;
            this.stop = end;
            this.numberOfPoints = numberOfPoints;
        }
        step = (this.numberOfPoints == 1) ? 0.0 : (this.stop - this.start) / (this.numberOfPoints - 1);
        if (data != null) {
            if (data.length != this.numberOfPoints) {
                throw new IllegalStateException(
                        "inconsistent vector dimension : " + this.numberOfPoints + " != " + data.length);
            }
            vec = data;
        }
    }

    private LinSpace(LinSpace other) {
        start = other.start;
        stop = other.stop;
        numberOfPoints = other.numberOfPoints;
        step = other.step;
        if (other.vec != null) {
            vec = other.vec.clone();
        }
    }

    public double spacing() {
        return Math.abs(step);
    }

    public int size() {
        return numberOfPoints;
    }

    public double start() {
        return start;
    }

    public double end() {
        return stop;
    }

    public DIndexIterator iterator() {
        return new DblIt(this);
    }

    /**
     * The abscissa at the 1-based position {@code pos}, unchecked. Every
     * accessor of this class goes through here, so that the same position
     * cannot have two values. Both endpoints come back exactly; the points in
     * between are {@code start + (pos - 1) * step}, which does not accumulate
     * rounding along the way the way a running sum does.
     * 
     * @param pos
     *            1-based position, in {@code [1, numberOfPoints]}
     * @return the abscissa at that position
     */
    private double abscissa(int pos) {
        if (pos == 1) {
            return start;
        }
        if (pos == numberOfPoints) {
            return stop;
        }
        return start + ((pos - 1) * step);
    }

    public LinSpace slice(int from, int to) {
        checkPosition(from, "from");
        checkPosition(to, "to");
        if (from > to) {
            int tmp = to;
            to = from;
            from = tmp;
        }
        if (from == 1 && to == numberOfPoints) {
            return new LinSpace(this);
        }
        boolean c = vec != null;
        // from >= 1
        // to <= n
        if (from == to) {
            int count = 1;
            double point = abscissa(from);
            return new LinSpace(point, point, count, c ? new double[] { vec[from - 1] } : null);
        }
        int count = 1 + to - from;
        double begin = abscissa(from);
        double end = abscissa(to);
        double[] d = c ? Arrays.copyOfRange(vec, from - 1, to) : null;
        return new LinSpace(begin, end, count, d);
    }

    public LinSpace sliceTo(int to) {
        checkPosition(to, "to");
        boolean c = vec != null;
        if (to == 1) {
            return new LinSpace(start, start, 1, c ? new double[] { vec[0] } : null);
        }
        if (to == numberOfPoints) {
            return new LinSpace(this);
        }
        // to is in [2..n-1]
        int count = to;
        double[] d = c ? Arrays.copyOfRange(vec, 0, to) : null;
        return new LinSpace(start, abscissa(to), count, d);
    }

    public LinSpace sliceFrom(int from) {
        checkPosition(from, "from");
        boolean c = vec != null;
        if (from == 1) {
            return new LinSpace(this);
        }
        if (from == numberOfPoints) {
            return new LinSpace(stop, stop, 1, c ? new double[] { vec[numberOfPoints - 1] } : null);
        }
        // from is in [2..n-1]
        int count = 1 + numberOfPoints - from;
        double[] d = c ? Arrays.copyOfRange(vec, from - 1, from - 1 + count) : null;
        return new LinSpace(abscissa(from), stop, count, d);
    }

    /**
     * The abscissa at the 1-based position {@code pos}. Both endpoints are
     * exact: {@code point(1)} is {@code start} and {@code point(size())} is
     * {@code end()}.
     * 
     * @param pos
     *            1-based position of the point
     * @return the abscissa at that position
     * @throws IndexOutOfBoundsException
     *             if {@code pos} is not in {@code [1, size()]}
     */
    public double point(int pos) {
        checkPosition(pos, "pos");
        return abscissa(pos);
    }

    public LinSpace allocate() {
        vec = new double[numberOfPoints];
        return this;
    }

    /**
     * All abscissas, in order, as a fresh array. Entry {@code i} is
     * {@code point(i + 1)}, exactly.
     * 
     * @return a new array of {@link #size()} abscissas
     */
    // escape hatch
    public double[] points() {
        double[] points = new double[numberOfPoints];
        for (int i = 0; i < points.length; ++i) {
            points[i] = abscissa(i + 1);
        }
        return points;
    }

    public double value(int pos) {
        checkPosition(pos, "pos");
        if (vec == null) {
            throw new NoSuchElementException("no data");
        }
        return vec[pos - 1];
    }

    // escape hatch
    public double[] values() {
        if (vec == null) {
            throw new NoSuchElementException("no data");
        }
        return vec;
    }

    public LinSpace setValue(int pos, double x) {
        checkPosition(pos, "pos");
        if (vec == null) {
            throw new NoSuchElementException("no data");
        }
        vec[pos - 1] = x;
        return this;
    }

    public DForEach forEach() {
        return new DblForEach(this);
    }

    public DForEachBi forEachBi() {
        if (!hasValues()) {
            throw new NoSuchElementException("no data");
        }
        return new DblForEachBi(this, vec);
    }

    /**
     * A new {@code LinSpace} over the same points, holding the value of
     * {@code fun} at each of them. {@code fun} is evaluated at
     * {@code point(i)} for every {@code i}, ascending and descending ranges
     * alike, so the value at the last point is the value at {@link #end()}.
     * 
     * @param fun
     *            the function to evaluate at each of the points
     * @return a new {@code LinSpace} carrying the values of {@code fun}
     */
    public LinSpace eval(DFunction fun) {
        LinSpace result = new LinSpace(start, stop, numberOfPoints, new double[numberOfPoints]);
        double[] y = result.vec;
        for (int i = 0; i < y.length; ++i) {
            y[i] = fun.apply(abscissa(i + 1));
        }
        return result;
    }

    private static final class DblForEach implements DForEach {
        private final LinSpace lsp;
        private final int total;
        private int remaining;

        DblForEach(LinSpace lsp) {
            this.lsp = lsp;
            total = lsp.numberOfPoints;
            remaining = total;
        }

        @Override
        public void forEachRemaining(DConsumer action) {
            while (remaining > 0) {
                int pos = total - remaining + 1;
                --remaining;
                action.accept(lsp.abscissa(pos));
            }
        }

        @Override
        public boolean tryAdvance(DConsumer action) {
            if (remaining > 0) {
                int pos = total - remaining + 1;
                --remaining;
                action.accept(lsp.abscissa(pos));
                return true;
            }
            return false;
        }
    }

    private static final class DblForEachBi implements DForEachBi {

        private final LinSpace lsp;
        private final int total;
        private final double[] data;
        private int remaining;

        DblForEachBi(LinSpace lsp, double[] vec) {
            this.lsp = lsp;
            total = lsp.numberOfPoints;
            remaining = total;
            data = vec;
        }

        @Override
        public void forEachRemaining(DBiConsumer action) {
            while (remaining > 0) {
                int pos = total - remaining + 1;
                --remaining;
                action.accept(lsp.abscissa(pos), data[pos - 1]);
            }
        }

        @Override
        public boolean tryAdvance(DBiConsumer action) {
            if (remaining > 0) {
                int pos = total - remaining + 1;
                --remaining;
                action.accept(lsp.abscissa(pos), data[pos - 1]);
                return true;
            }
            return false;
        }
    }

    private static final class DblIt implements DIndexIterator {
        private final LinSpace lsp;
        private final int total;
        private int remaining;

        DblIt(LinSpace lsp) {
            this.lsp = lsp;
            total = lsp.numberOfPoints;
            remaining = total;
        }

        @Override
        public boolean hasNext() {
            return remaining > 0;
        }

        @Override
        public int nextIndex() {
            if (remaining > 0) {
                return total - remaining + 1;
            }
            throw new NoSuchElementException("exhausted");
        }

        @Override
        public double next() {
            if (remaining > 0) {
                int pos = total - remaining + 1;
                --remaining;
                return lsp.abscissa(pos);
            }
            throw new NoSuchElementException("exhausted");
        }
    }

    /**
     * Returns {@code 128} evenly spaced points between {@code start} and
     * {@code end} (including the interval endpoints).
     * 
     * @param start
     *            start point of interval (included)
     * @param end
     *            endpoint of interval (included)
     * @return sample interval containing {@code 128} points
     */
    public static LinSpace linspace(double start, double end) {
        return new LinSpace(start, end, 128, null);
    }

    /**
     * Returns {@code numberOfPoints} evenly spaced points between
     * {@code start} and {@code end} (including the interval endpoints).
     * <p>
     * The count is honored as given. For {@code start == end} the result is
     * {@code numberOfPoints} copies of that value, and for
     * {@code numberOfPoints == 1} it is the single point {@code start}.
     * 
     * @param start
     *            start point of interval (included)
     * @param end
     *            endpoint of interval (included)
     * @param numberOfPoints
     *            the number of points, must be strictly positive
     * @return sample interval containing {@code numberOfPoints} points
     * @throws IllegalArgumentException
     *             if {@code numberOfPoints} is not strictly positive, or if
     *             {@code start} or {@code end} is infinite or {@code NaN}
     */
    public static LinSpace linspace(double start, double end, int numberOfPoints) {
        return new LinSpace(start, end, numberOfPoints, null);
    }

    /**
     * Returns {@code numberOfPoints} evenly spaced points between
     * {@code start} and {@code end}, each carrying the value of {@code fun} at
     * that point. {@code fun} is evaluated once per point, so a degenerate
     * interval yields {@code numberOfPoints} copies of the same value.
     * 
     * @param start
     *            start point of interval (included)
     * @param end
     *            endpoint of interval (included)
     * @param numberOfPoints
     *            the number of points, must be strictly positive
     * @param fun
     *            the function to evaluate at each of the points
     * @return sample interval of {@code numberOfPoints} points holding the
     *         values of {@code fun}
     * @throws IllegalArgumentException
     *             if {@code numberOfPoints} is not strictly positive, or if
     *             {@code start} or {@code end} is infinite or {@code NaN}
     */
    public static LinSpace compute(double start, double end, int numberOfPoints, DFunction fun) {
        return linspace(start, end, numberOfPoints).eval(fun);
    }

    public static LinSpace centeredIntIndexed(double[] data) {
        final int length = data.length;
        if (length < 1) {
            throw new IllegalArgumentException("data.length must be strictly positive : 0");
        }
        if (length == 1) {
            return new LinSpace(0.0, 0.0, 1, new double[] {data[0]});
        }
        double sym = (((double) length) - 1.0) / 2.0;
        double start = (length % 2 != 0) ? -Math.floor(sym) : -Math.floor(sym) - 1.0;
        double end = start + length - 1.0;
        return new LinSpace(start, end, length, data.clone());
    }

    public static LinSpace centeredDoubleIndexed(double[] data) {
        final int length = data.length;
        if (length < 1) {
            throw new IllegalArgumentException("data.length must be strictly positive : 0");
        }
        if (length == 1) {
            return new LinSpace(0.0, 0.0, 1, new double[] {data[0]});
        }
        double sym = (((double) length) - 1.0) / 2.0;
        return new LinSpace(-sym, sym, length, data.clone());
    }

    public boolean hasValues() {
        return vec != null;
    }

    @Override
    public String toString() {
        return "1x" + numberOfPoints + " :  [" + start + "  ...  " + stop + "]";
    }

    private void checkPosition(int pos, String name) {
        if (pos < 1 || pos > numberOfPoints) {
            throw new IndexOutOfBoundsException(
                    name + " = " + pos + " (indexes are 1-based, size is " + numberOfPoints + ")");
        }
    }

    private static void checkArg(double a, String name) {
        if (isBadNum(a)) {
            throw new IllegalArgumentException("Bad argument : " + name + " (Inf or NaN)");
        }
    }

    private static boolean isBadNum(double v) {
        if (Double.isInfinite(v) || Double.isNaN(v)) {
            return true;
        }
        return false;
    }
}
