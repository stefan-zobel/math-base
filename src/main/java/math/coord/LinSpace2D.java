package math.coord;

import java.util.NoSuchElementException;

import math.fun.DBiConsumer;
import math.fun.DBiFunction;
import math.fun.DTriConsumer;

/**
 * A rectangular grid of evenly spaced points, the tensor product of two
 * {@link LinSpace} axes: the point at the 1-based position {@code (i, j)} is
 * {@code (x().point(i), y().point(j))}. Optionally the grid carries one value
 * per point.
 * <p>
 * All the abscissa arithmetic lives in {@link LinSpace} and none of it is
 * repeated here, so a grid inherits its properties: both endpoints of either
 * axis are exact, a degenerate axis is not folded onto one point, and a span
 * beyond {@code Double.MAX_VALUE} stays finite.
 * <p>
 * Values are held in one flat array in column-major order with the first
 * variable running fastest: {@code values()[(j - 1) * sizeX() + (i - 1)]} is
 * the value at {@code (i, j)}. That is the layout of {@code math.linalg}, so
 * the array is directly an {@code sizeX() x sizeY()} matrix buffer, and the
 * indexing is the one numpy calls {@code meshgrid(..., indexing='ij')}.
 * <p>
 * Note that indexes are 1-based!
 *
 * @see <a href="https://en.wikipedia.org/wiki/Regular_grid">Regular grid</a>
 * @since 1.5.2
 */
public final class LinSpace2D {

    private final LinSpace x;
    private final LinSpace y;
    private final int nx;
    private final int ny;
    private double[] vec;

    private LinSpace2D(LinSpace x, LinSpace y, double[] data) {
        this.x = x;
        this.y = y;
        nx = x.size();
        ny = y.size();
        int total = checkTotal(nx, ny);
        if (data != null) {
            if (data.length != total) {
                throw new IllegalStateException(
                        "inconsistent grid dimension : " + total + " != " + data.length);
            }
            vec = data;
        }
    }

    /**
     * A copy of {@code axis} without its values. A {@code LinSpace} is
     * determined by its two endpoints and its size, so this reproduces the
     * abscissas bit for bit.
     */
    private static LinSpace axis(LinSpace ls) {
        return LinSpace.linspace(ls.start(), ls.end(), ls.size());
    }

    /**
     * A grid of {@code nx} times {@code ny} evenly spaced points over the
     * rectangle spanned by the two intervals, both endpoints of either
     * interval included.
     *
     * @param xStart
     *            start point of the interval in the first variable (included)
     * @param xEnd
     *            endpoint of the interval in the first variable (included)
     * @param nx
     *            the number of points in the first variable, must be strictly
     *            positive
     * @param yStart
     *            start point of the interval in the second variable (included)
     * @param yEnd
     *            endpoint of the interval in the second variable (included)
     * @param ny
     *            the number of points in the second variable, must be strictly
     *            positive
     * @return a new grid of {@code nx * ny} points
     * @throws IllegalArgumentException
     *             if either count is not strictly positive, if any of the four
     *             bounds is infinite or {@code NaN}, or if {@code nx * ny}
     *             exceeds {@code Integer.MAX_VALUE}
     */
    public static LinSpace2D linspace2D(double xStart, double xEnd, int nx, double yStart, double yEnd, int ny) {
        return new LinSpace2D(LinSpace.linspace(xStart, xEnd, nx), LinSpace.linspace(yStart, yEnd, ny), null);
    }

    /**
     * A grid over two axes that already exist. Both are copied without their
     * values, so the arguments stay independent of the grid.
     *
     * @param x
     *            the axis of the first variable
     * @param y
     *            the axis of the second variable
     * @return a new grid of {@code x.size() * y.size()} points
     * @throws IllegalArgumentException
     *             if {@code x.size() * y.size()} exceeds
     *             {@code Integer.MAX_VALUE}
     */
    public static LinSpace2D over(LinSpace x, LinSpace y) {
        if (x == null || y == null) {
            throw new NullPointerException("x or y");
        }
        return new LinSpace2D(axis(x), axis(y), null);
    }

    /**
     * A grid of {@code nx} times {@code ny} evenly spaced points, each
     * carrying the value of {@code fun} at that point. {@code fun} is
     * evaluated once per point.
     *
     * @param xStart
     *            start point of the interval in the first variable (included)
     * @param xEnd
     *            endpoint of the interval in the first variable (included)
     * @param nx
     *            the number of points in the first variable, must be strictly
     *            positive
     * @param yStart
     *            start point of the interval in the second variable (included)
     * @param yEnd
     *            endpoint of the interval in the second variable (included)
     * @param ny
     *            the number of points in the second variable, must be strictly
     *            positive
     * @param fun
     *            the function to evaluate at each of the points
     * @return a new grid of {@code nx * ny} points holding the values of
     *         {@code fun}
     * @throws IllegalArgumentException
     *             if either count is not strictly positive, if any of the four
     *             bounds is infinite or {@code NaN}, or if {@code nx * ny}
     *             exceeds {@code Integer.MAX_VALUE}
     */
    public static LinSpace2D compute(double xStart, double xEnd, int nx, double yStart, double yEnd, int ny,
            DBiFunction fun) {
        return linspace2D(xStart, xEnd, nx, yStart, yEnd, ny).eval(fun);
    }

    /**
     * A grid over integer abscissas centered on zero in both variables,
     * holding a copy of {@code data}. Entry {@code data[i - 1][j - 1]} is the
     * value at the 1-based position {@code (i, j)}, so the <em>first</em>
     * index runs along the first variable, as it does in {@link #toArrays()}.
     * A kernel written out as an array literal is therefore read column by
     * column rather than row by row:
     * {@code {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}}} differentiates in
     * {@code x}, not in {@code y}.
     * <p>
     * Each variable is centered on its own: for an odd length its abscissas
     * are symmetric, {@code -2, -1, 0, 1, 2}, and for an even one there is one
     * more below zero than above, {@code -2, -1, 0, 1}.
     *
     * @param data
     *            the values, copied into the result, rectangular and not
     *            empty in either direction
     * @return a new grid of {@code data.length * data[0].length} points
     * @throws IllegalArgumentException
     *             if {@code data} is empty in either direction, if its rows
     *             differ in length, or if it holds more than
     *             {@code Integer.MAX_VALUE} values
     */
    public static LinSpace2D centeredIntIndexed(double[][] data) {
        int ny = checkRectangular(data);
        int nx = data.length;
        return new LinSpace2D(axis(LinSpace.centeredIntIndexed(new double[nx])),
                axis(LinSpace.centeredIntIndexed(new double[ny])), flatten(data, nx, ny));
    }

    /**
     * A grid over abscissas centered on zero at a spacing of one in both
     * variables, holding a copy of {@code data}. Entry
     * {@code data[i - 1][j - 1]} is the value at the 1-based position
     * {@code (i, j)}, the first index running along the first variable, as in
     * {@link #centeredIntIndexed(double[][])}.
     * <p>
     * Each variable is centered on its own, and its abscissas are symmetric
     * for either length: {@code -2, -1, 0, 1, 2} for an odd one and
     * {@code -1.5, -0.5, 0.5, 1.5} for an even one.
     *
     * @param data
     *            the values, copied into the result, rectangular and not
     *            empty in either direction
     * @return a new grid of {@code data.length * data[0].length} points
     * @throws IllegalArgumentException
     *             if {@code data} is empty in either direction, if its rows
     *             differ in length, or if it holds more than
     *             {@code Integer.MAX_VALUE} values
     */
    public static LinSpace2D centeredDoubleIndexed(double[][] data) {
        int ny = checkRectangular(data);
        int nx = data.length;
        return new LinSpace2D(axis(LinSpace.centeredDoubleIndexed(new double[nx])),
                axis(LinSpace.centeredDoubleIndexed(new double[ny])), flatten(data, nx, ny));
    }

    /**
     * The axis of the first variable itself, not a copy. It carries no values,
     * and giving it any has no effect on this grid, which reads nothing but
     * its abscissas.
     *
     * @return the axis of the first variable
     */
    // escape hatch
    public LinSpace x() {
        return x;
    }

    /**
     * The axis of the second variable itself, not a copy. It carries no
     * values, and giving it any has no effect on this grid, which reads
     * nothing but its abscissas.
     *
     * @return the axis of the second variable
     */
    // escape hatch
    public LinSpace y() {
        return y;
    }

    /**
     * The number of points along the first variable.
     *
     * @return the number of points in {@code x}, at least 1
     */
    public int sizeX() {
        return nx;
    }

    /**
     * The number of points along the second variable.
     *
     * @return the number of points in {@code y}, at least 1
     */
    public int sizeY() {
        return ny;
    }

    /**
     * The number of points of the whole grid, {@code sizeX() * sizeY()}.
     *
     * @return the number of grid points, at least 1
     */
    public int size() {
        return nx * ny;
    }

    /**
     * The abscissa in the first variable at the 1-based position {@code i},
     * exactly {@code x().point(i)}.
     *
     * @param i
     *            1-based position along the first variable
     * @return the abscissa at that position
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not in {@code [1, sizeX()]}
     */
    public double pointX(int i) {
        return x.point(i);
    }

    /**
     * The abscissa in the second variable at the 1-based position {@code j},
     * exactly {@code y().point(j)}.
     *
     * @param j
     *            1-based position along the second variable
     * @return the abscissa at that position
     * @throws IndexOutOfBoundsException
     *             if {@code j} is not in {@code [1, sizeY()]}
     */
    public double pointY(int j) {
        return y.point(j);
    }

    /**
     * Gives this grid a value array of {@link #size()} zeros, discarding any
     * values it already held.
     *
     * @return this grid
     */
    public LinSpace2D allocate() {
        vec = new double[nx * ny];
        return this;
    }

    /**
     * Whether this grid holds values, i.e. whether {@link #values()} and
     * {@link #value(int, int)} can be called.
     *
     * @return {@code true} if a value array is present
     */
    public boolean hasValues() {
        return vec != null;
    }

    /**
     * The value at the 1-based position {@code (i, j)}.
     *
     * @param i
     *            1-based position along the first variable
     * @param j
     *            1-based position along the second variable
     * @return the value at that position
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not in {@code [1, sizeX()]} or {@code j} is
     *             not in {@code [1, sizeY()]}
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public double value(int i, int j) {
        checkX(i, "i");
        checkY(j, "j");
        requireValues();
        return vec[(j - 1) * nx + (i - 1)];
    }

    /**
     * Sets the value at the 1-based position {@code (i, j)}.
     *
     * @param i
     *            1-based position along the first variable
     * @param j
     *            1-based position along the second variable
     * @param v
     *            the value to store
     * @return this grid
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not in {@code [1, sizeX()]} or {@code j} is
     *             not in {@code [1, sizeY()]}
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public LinSpace2D setValue(int i, int j, double v) {
        checkX(i, "i");
        checkY(j, "j");
        requireValues();
        vec[(j - 1) * nx + (i - 1)] = v;
        return this;
    }

    /**
     * The value array itself, not a copy: writing into the returned array
     * changes this grid. Entry {@code (j - 1) * sizeX() + (i - 1)} is the
     * value at {@code (i, j)}, so the array is an
     * {@code sizeX() x sizeY()} matrix in column-major order.
     *
     * @return the internal array of {@link #size()} values
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    // escape hatch
    public double[] values() {
        requireValues();
        return vec;
    }

    /**
     * A new grid over the same points, holding the value of {@code fun} at
     * each of them. {@code fun} is evaluated at
     * {@code (pointX(i), pointY(j))} for every position, ascending and
     * descending axes alike.
     *
     * @param fun
     *            the function to evaluate at each of the points
     * @return a new grid carrying the values of {@code fun}
     */
    public LinSpace2D eval(DBiFunction fun) {
        double[] xs = x.points();
        double[] ys = y.points();
        double[] z = new double[nx * ny];
        for (int j = 0; j < ny; ++j) {
            double yj = ys[j];
            int off = j * nx;
            for (int i = 0; i < nx; ++i) {
                z[off + i] = fun.apply(xs[i], yj);
            }
        }
        return new LinSpace2D(x, y, z);
    }

    /**
     * The sub-grid between the 1-based positions {@code fromX} and
     * {@code toX} along the first variable and {@code fromY} and {@code toY}
     * along the second, all four included, as a new grid. Either pair may be
     * given in either order. The sub-grid begins and ends exactly on the
     * points of this grid; its interior points are its own, because each of
     * its axes carries its own spacing. Values, if this grid has any, are
     * copied along.
     *
     * @param fromX
     *            1-based position of one end along the first variable
     * @param toX
     *            1-based position of the other end along the first variable
     * @param fromY
     *            1-based position of one end along the second variable
     * @param toY
     *            1-based position of the other end along the second variable
     * @return a new grid of {@code (|toX - fromX| + 1) * (|toY - fromY| + 1)}
     *         points
     * @throws IndexOutOfBoundsException
     *             if any of the four positions is out of range
     */
    public LinSpace2D slice(int fromX, int toX, int fromY, int toY) {
        checkX(fromX, "fromX");
        checkX(toX, "toX");
        checkY(fromY, "fromY");
        checkY(toY, "toY");
        if (fromX > toX) {
            int tmp = toX;
            toX = fromX;
            fromX = tmp;
        }
        if (fromY > toY) {
            int tmp = toY;
            toY = fromY;
            fromY = tmp;
        }
        int mx = 1 + toX - fromX;
        int my = 1 + toY - fromY;
        double[] d = null;
        if (vec != null) {
            d = new double[mx * my];
            for (int j = 0; j < my; ++j) {
                System.arraycopy(vec, (fromY - 1 + j) * nx + (fromX - 1), d, j * mx, mx);
            }
        }
        return new LinSpace2D(axis(x.slice(fromX, toX)), axis(y.slice(fromY, toY)), d);
    }

    /**
     * The cut through this grid along the first variable at the 1-based
     * position {@code j} of the second: a {@code LinSpace} over the axis
     * {@link #x()}, carrying the {@link #sizeX()} values at
     * {@code (1, j) ... (sizeX(), j)}.
     *
     * @param j
     *            1-based position along the second variable
     * @return a new {@code LinSpace} of {@code sizeX()} points with values
     * @throws IndexOutOfBoundsException
     *             if {@code j} is not in {@code [1, sizeY()]}
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public LinSpace alongX(int j) {
        checkY(j, "j");
        requireValues();
        LinSpace line = axis(x).allocate();
        System.arraycopy(vec, (j - 1) * nx, line.values(), 0, nx);
        return line;
    }

    /**
     * The cut through this grid along the second variable at the 1-based
     * position {@code i} of the first: a {@code LinSpace} over the axis
     * {@link #y()}, carrying the {@link #sizeY()} values at
     * {@code (i, 1) ... (i, sizeY())}.
     *
     * @param i
     *            1-based position along the first variable
     * @return a new {@code LinSpace} of {@code sizeY()} points with values
     * @throws IndexOutOfBoundsException
     *             if {@code i} is not in {@code [1, sizeX()]}
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public LinSpace alongY(int i) {
        checkX(i, "i");
        requireValues();
        LinSpace line = axis(y).allocate();
        double[] d = line.values();
        for (int j = 0; j < ny; ++j) {
            d[j] = vec[j * nx + (i - 1)];
        }
        return line;
    }

    /**
     * A new grid with the two variables exchanged: its value at
     * {@code (j, i)} is the value of this grid at {@code (i, j)}.
     *
     * @return a new grid over the swapped axes
     */
    public LinSpace2D transpose() {
        double[] d = null;
        if (vec != null) {
            d = new double[nx * ny];
            for (int j = 0; j < ny; ++j) {
                int off = j * nx;
                for (int i = 0; i < nx; ++i) {
                    d[i * ny + j] = vec[off + i];
                }
            }
        }
        return new LinSpace2D(y, x, d);
    }

    /**
     * The 1-based index {@code i} of the cell along the first variable that
     * contains {@code px}, so that {@code px} lies between {@code pointX(i)}
     * and {@code pointX(i + 1)}. The last point belongs to the last cell, and
     * a {@code px} outside the interval answers the first or the last cell.
     * <p>
     * Because the points are evenly spaced this is arithmetic, not a search.
     *
     * @param px
     *            the abscissa in the first variable
     * @return the 1-based cell index, in {@code [1, sizeX() - 1]}
     * @throws IllegalArgumentException
     *             if {@code px} is infinite or {@code NaN}
     * @throws IllegalStateException
     *             if {@code sizeX() == 1}, as a single point has no cell
     */
    public int cellX(double px) {
        return cell(x, px, "px");
    }

    /**
     * The 1-based index {@code j} of the cell along the second variable that
     * contains {@code py}, so that {@code py} lies between {@code pointY(j)}
     * and {@code pointY(j + 1)}. The last point belongs to the last cell, and
     * a {@code py} outside the interval answers the first or the last cell.
     * <p>
     * Because the points are evenly spaced this is arithmetic, not a search.
     *
     * @param py
     *            the abscissa in the second variable
     * @return the 1-based cell index, in {@code [1, sizeY() - 1]}
     * @throws IllegalArgumentException
     *             if {@code py} is infinite or {@code NaN}
     * @throws IllegalStateException
     *             if {@code sizeY() == 1}, as a single point has no cell
     */
    public int cellY(double py) {
        return cell(y, py, "py");
    }

    private static int cell(LinSpace ax, double q, String name) {
        checkArg(q, name);
        int n = ax.size();
        if (n == 1) {
            throw new IllegalStateException("an axis of one point has no cell");
        }
        double h = ax.spacing();
        if (h == 0.0) {
            // a degenerate axis: every one of its cells has width zero
            return 1;
        }
        boolean up = ax.end() >= ax.start();
        if (!up) {
            h = -h;
        }
        // measure from whichever endpoint is the nearer one, the way the
        // abscissas themselves are formed. The forward form alone is exact on
        // an ordinary axis but walks up to 499 positions on a span wider than
        // the double range, because there the far abscissas are measured from
        // the far end
        double mid = ax.start() / 2.0 + ax.end() / 2.0;
        double t;
        if (up ? (q <= mid) : (q >= mid)) {
            t = (q - ax.start()) / h;
        } else {
            t = (n - 1) - ((ax.end() - q) / h);
        }
        int i;
        if (!(t > 0.0)) {
            i = 1;
        } else if (t >= n - 1) {
            i = n - 1;
        } else {
            i = 1 + (int) t;
            if (i > n - 1) {
                i = n - 1;
            }
        }
        // a guard, not a search: measured over 130200 lookups on ordinary,
        // degenerate, descending and wide-span axes it never took a step
        while (i > 1 && (up ? q < ax.point(i) : q > ax.point(i))) {
            --i;
        }
        while (i < n - 1 && (up ? q > ax.point(i + 1) : q < ax.point(i + 1))) {
            ++i;
        }
        return i;
    }

    /**
     * Performs {@code action} on every grid point, the first variable running
     * fastest, which is the order the values are stored in.
     *
     * @param action
     *            the action, receiving {@code (pointX(i), pointY(j))}
     */
    public void forEachPoint(DBiConsumer action) {
        double[] xs = x.points();
        double[] ys = y.points();
        for (int j = 0; j < ny; ++j) {
            double yj = ys[j];
            for (int i = 0; i < nx; ++i) {
                action.accept(xs[i], yj);
            }
        }
    }

    /**
     * Performs {@code action} on every grid point paired with its value, the
     * first variable running fastest, which is the order the values are stored
     * in.
     *
     * @param action
     *            the action, receiving
     *            {@code (pointX(i), pointY(j), value(i, j))}
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public void forEachValue(DTriConsumer action) {
        requireValues();
        double[] xs = x.points();
        double[] ys = y.points();
        for (int j = 0; j < ny; ++j) {
            double yj = ys[j];
            int off = j * nx;
            for (int i = 0; i < nx; ++i) {
                action.accept(xs[i], yj, vec[off + i]);
            }
        }
    }

    /**
     * The values as a fresh array of arrays, {@code z[i - 1][j - 1]} being the
     * value at {@code (i, j)}. This is the layout the interpolators of
     * {@code math.fit} take, next to {@code x().points()} and
     * {@code y().points()}.
     *
     * @return a new {@code sizeX() x sizeY()} array of values
     * @throws NoSuchElementException
     *             if this grid holds no values
     */
    public double[][] toArrays() {
        requireValues();
        double[][] z = new double[nx][ny];
        for (int j = 0; j < ny; ++j) {
            int off = j * nx;
            for (int i = 0; i < nx; ++i) {
                z[i][j] = vec[off + i];
            }
        }
        return z;
    }

    @Override
    public String toString() {
        return nx + "x" + ny + " :  [" + x.start() + "  ...  " + x.end() + "] x [" + y.start() + "  ...  "
                + y.end() + "]";
    }

    private void requireValues() {
        if (vec == null) {
            throw new NoSuchElementException("no data");
        }
    }

    private void checkX(int i, String name) {
        if (i < 1 || i > nx) {
            throw new IndexOutOfBoundsException(
                    name + " = " + i + " (indexes are 1-based, size is " + nx + ")");
        }
    }

    private void checkY(int j, String name) {
        if (j < 1 || j > ny) {
            throw new IndexOutOfBoundsException(
                    name + " = " + j + " (indexes are 1-based, size is " + ny + ")");
        }
    }

    private static void checkArg(double a, String name) {
        if (Double.isInfinite(a) || Double.isNaN(a)) {
            throw new IllegalArgumentException("Bad argument : " + name + " (Inf or NaN)");
        }
    }

    private static int checkTotal(int nx, int ny) {
        long total = (long) nx * (long) ny;
        if (total > Integer.MAX_VALUE) {
            // the grid is addressed by a single int, so this has to be caught
            // whether values are allocated or not: size() must not lie
            throw new IllegalArgumentException(
                    "grid too large : " + nx + " x " + ny + " = " + total + " points");
        }
        return (int) total;
    }

    /**
     * Checks that {@code data} is a non-empty rectangle and answers the
     * length of its rows. A jagged argument is the one way of getting the
     * shape wrong that the one-dimensional class cannot have.
     */
    private static int checkRectangular(double[][] data) {
        if (data == null) {
            throw new NullPointerException("data");
        }
        if (data.length < 1) {
            throw new IllegalArgumentException("data.length must be strictly positive : 0");
        }
        if (data[0] == null) {
            throw new NullPointerException("data[0]");
        }
        int ny = data[0].length;
        if (ny < 1) {
            throw new IllegalArgumentException("data[0].length must be strictly positive : 0");
        }
        for (int i = 1; i < data.length; ++i) {
            if (data[i] == null) {
                throw new NullPointerException("data[" + i + "]");
            }
            if (data[i].length != ny) {
                throw new IllegalArgumentException("data must be rectangular : data[" + i + "].length = "
                        + data[i].length + " != " + ny);
            }
        }
        return ny;
    }

    private static double[] flatten(double[][] data, int nx, int ny) {
        double[] d = new double[checkTotal(nx, ny)];
        for (int i = 0; i < nx; ++i) {
            double[] row = data[i];
            for (int j = 0; j < ny; ++j) {
                d[j * nx + i] = row[j];
            }
        }
        return d;
    }
}
