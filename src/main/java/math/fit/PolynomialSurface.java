package math.fit;

import java.util.Arrays;

import math.fun.DBiFunction;

/**
 * A piecewise polynomial surface of degree at most three in each variable over
 * a rectangular grid, evaluated in the local coordinates
 * {@code t = x - knotX(i)} and {@code v = y - knotY(j)} of the cell that
 * contains the point.
 * <p>
 * Instances come from {@link BilinearInterpolator#interpolate(double[], double[], double[][])}
 * and {@link BicubicInterpolator#interpolate(double[], double[], double[][])},
 * which differ in how they choose the sixteen coefficients per cell. This class
 * only evaluates them, but it evaluates both first partial derivatives and the
 * exact double integral over a rectangle as well as the value, all in closed
 * form -- the integral separates per cell because the cell polynomial is a
 * tensor product. Since it is a {@link math.fun.DBiFunction}, it can be handed
 * to {@code math.solve.AdaptiveGaussKronrod} directly.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Bicubic_interpolation">Bicubic
 * interpolation</a>
 *
 * @since 1.5.2
 */
public final class PolynomialSurface implements DBiFunction {

    private final double[] x;
    private final double[] y;
    /** coefficients[i][j][a * 4 + b] multiplies {@code t^a v^b} on cell (i, j). */
    private final double[][][] coefficients;

    /**
     * @param x
     *            the abscissae in the first variable, strictly increasing
     * @param y
     *            the abscissae in the second variable, strictly increasing
     * @param coefficients
     *            {@code (x.length - 1) x (y.length - 1)} arrays of 16
     */
    PolynomialSurface(double[] x, double[] y, double[][][] coefficients) {
        this.x = x.clone();
        this.y = y.clone();
        this.coefficients = new double[coefficients.length][][];
        for (int i = 0; i < coefficients.length; i++) {
            this.coefficients[i] = new double[coefficients[i].length][];
            for (int j = 0; j < coefficients[i].length; j++) {
                this.coefficients[i][j] = coefficients[i][j].clone();
            }
        }
    }

    /**
     * The number of cells in the first variable, one less than the number of
     * abscissae.
     *
     * @return the cell count along {@code x}
     */
    public int cellsX() {
        return x.length - 1;
    }

    /**
     * The number of cells in the second variable.
     *
     * @return the cell count along {@code y}
     */
    public int cellsY() {
        return y.length - 1;
    }

    /**
     * The {@code i}-th abscissa in the first variable.
     *
     * @param i
     *            the index, {@code 0} to {@link #cellsX()} inclusive
     * @return that abscissa
     */
    public double knotX(int i) {
        if (i < 0 || i >= x.length) {
            throw new IllegalArgumentException("x index out of range [0, " + (x.length - 1) + "] : " + i);
        }
        return x[i];
    }

    /**
     * The {@code j}-th abscissa in the second variable.
     *
     * @param j
     *            the index, {@code 0} to {@link #cellsY()} inclusive
     * @return that abscissa
     */
    public double knotY(int j) {
        if (j < 0 || j >= y.length) {
            throw new IllegalArgumentException("y index out of range [0, " + (y.length - 1) + "] : " + j);
        }
        return y[j];
    }

    /**
     * The lower end of the grid in the first variable.
     *
     * @return the smallest {@code x} this surface is defined at
     */
    public double lowerBoundX() {
        return x[0];
    }

    /**
     * The upper end of the grid in the first variable.
     *
     * @return the largest {@code x} this surface is defined at
     */
    public double upperBoundX() {
        return x[x.length - 1];
    }

    /**
     * The lower end of the grid in the second variable.
     *
     * @return the smallest {@code y} this surface is defined at
     */
    public double lowerBoundY() {
        return y[0];
    }

    /**
     * The upper end of the grid in the second variable.
     *
     * @return the largest {@code y} this surface is defined at
     */
    public double upperBoundY() {
        return y[y.length - 1];
    }

    /**
     * A copy of the sixteen coefficients of one cell, where entry
     * {@code a * 4 + b} multiplies {@code t^a v^b} in the local coordinates
     * {@code t = x - knotX(cellX)} and {@code v = y - knotY(cellY)}.
     *
     * @param cellX
     *            the cell index along {@code x}, below {@link #cellsX()}
     * @param cellY
     *            the cell index along {@code y}, below {@link #cellsY()}
     * @return the sixteen coefficients of that cell
     */
    public double[] coefficients(int cellX, int cellY) {
        if (cellX < 0 || cellX >= coefficients.length) {
            throw new IllegalArgumentException(
                    "x cell index out of range [0, " + (coefficients.length - 1) + "] : " + cellX);
        }
        if (cellY < 0 || cellY >= coefficients[cellX].length) {
            throw new IllegalArgumentException(
                    "y cell index out of range [0, " + (coefficients[cellX].length - 1) + "] : " + cellY);
        }
        return coefficients[cellX][cellY].clone();
    }

    /**
     * The value of the surface at {@code (px, py)}.
     *
     * @param px
     *            the abscissa in the first variable, within the grid
     * @param py
     *            the abscissa in the second variable, within the grid
     * @return the interpolated value
     */
    public double value(double px, double py) {
        int i = cellOfX(px);
        int j = cellOfY(py);
        double t = px - x[i];
        double v = py - y[j];
        double[] c = coefficients[i][j];
        double sum = 0.0;
        for (int a = 3; a >= 0; a--) {
            sum = sum * t + inV(c, a, v);
        }
        return sum;
    }

    @Override
    public double apply(double px, double py) {
        return value(px, py);
    }

    /**
     * The partial derivative with respect to the first variable.
     *
     * @param px
     *            the abscissa in the first variable, within the grid
     * @param py
     *            the abscissa in the second variable, within the grid
     * @return {@code df/dx} at that point
     */
    public double dx(double px, double py) {
        int i = cellOfX(px);
        int j = cellOfY(py);
        double t = px - x[i];
        double v = py - y[j];
        double[] c = coefficients[i][j];
        double sum = 0.0;
        for (int a = 3; a >= 1; a--) {
            sum = sum * t + a * inV(c, a, v);
        }
        return sum;
    }

    /**
     * The partial derivative with respect to the second variable.
     *
     * @param px
     *            the abscissa in the first variable, within the grid
     * @param py
     *            the abscissa in the second variable, within the grid
     * @return {@code df/dy} at that point
     */
    public double dy(double px, double py) {
        int i = cellOfX(px);
        int j = cellOfY(py);
        double t = px - x[i];
        double v = py - y[j];
        double[] c = coefficients[i][j];
        double sum = 0.0;
        for (int a = 3; a >= 0; a--) {
            double inner = 0.0;
            for (int b = 3; b >= 1; b--) {
                inner = inner * v + b * c[a * 4 + b];
            }
            sum = sum * t + inner;
        }
        return sum;
    }

    /**
     * The integral of the surface over the rectangle, in closed form. It is
     * antisymmetric in each pair of limits, so exchanging {@code ax} and
     * {@code bx} negates the result.
     *
     * @param ax
     *            one limit in the first variable, within the grid
     * @param bx
     *            the other limit in the first variable, within the grid
     * @param ay
     *            one limit in the second variable, within the grid
     * @param by
     *            the other limit in the second variable, within the grid
     * @return the definite double integral
     */
    public double integrate(double ax, double bx, double ay, double by) {
        checkX(ax, "ax");
        checkX(bx, "bx");
        checkY(ay, "ay");
        checkY(by, "by");
        double sign = 1.0;
        double xLo = ax;
        double xHi = bx;
        if (xLo > xHi) {
            xLo = bx;
            xHi = ax;
            sign = -sign;
        }
        double yLo = ay;
        double yHi = by;
        if (yLo > yHi) {
            yLo = by;
            yHi = ay;
            sign = -sign;
        }
        if (xLo == xHi || yLo == yHi) {
            return 0.0;
        }
        double sum = 0.0;
        double[] ia = new double[4];
        double[] ib = new double[4];
        for (int i = 0; i < coefficients.length; i++) {
            double cellXLo = Math.max(xLo, x[i]);
            double cellXHi = Math.min(xHi, x[i + 1]);
            if (!(cellXHi > cellXLo)) {
                continue;
            }
            antiderivatives(cellXLo - x[i], cellXHi - x[i], ia);
            for (int j = 0; j < coefficients[i].length; j++) {
                double cellYLo = Math.max(yLo, y[j]);
                double cellYHi = Math.min(yHi, y[j + 1]);
                if (!(cellYHi > cellYLo)) {
                    continue;
                }
                antiderivatives(cellYLo - y[j], cellYHi - y[j], ib);
                double[] c = coefficients[i][j];
                double cell = 0.0;
                for (int a = 0; a < 4; a++) {
                    double inner = 0.0;
                    for (int b = 0; b < 4; b++) {
                        inner += c[a * 4 + b] * ib[b];
                    }
                    cell += inner * ia[a];
                }
                sum += cell;
            }
        }
        return sign * sum;
    }

    @Override
    public String toString() {
        return "PolynomialSurface[" + cellsX() + " x " + cellsY() + " cells over [" + lowerBoundX() + ", "
                + upperBoundX() + "] x [" + lowerBoundY() + ", " + upperBoundY() + "]]";
    }

    /** The inner Horner in {@code v} of the coefficients of one power of {@code t}. */
    private static double inV(double[] c, int a, double v) {
        double inner = 0.0;
        for (int b = 3; b >= 0; b--) {
            inner = inner * v + c[a * 4 + b];
        }
        return inner;
    }

    /** {@code out[a]} becomes the integral of {@code s^a} from {@code lo} to {@code hi}. */
    private static void antiderivatives(double lo, double hi, double[] out) {
        double powLo = lo;
        double powHi = hi;
        for (int a = 0; a < 4; a++) {
            out[a] = (powHi - powLo) / (a + 1);
            powLo *= lo;
            powHi *= hi;
        }
    }

    private int cellOfX(double px) {
        checkX(px, "x");
        return cellOf(x, px, coefficients.length);
    }

    private int cellOfY(double py) {
        checkY(py, "y");
        return cellOf(y, py, coefficients[0].length);
    }

    private static int cellOf(double[] knots, double q, int cells) {
        int i = Arrays.binarySearch(knots, q);
        if (i < 0) {
            i = -i - 2;
        }
        // the last abscissa belongs to the last cell
        if (i >= cells) {
            i = cells - 1;
        }
        return i;
    }

    private void checkX(double px, String name) {
        // NaN fails both comparisons, so it is rejected explicitly
        if (Double.isNaN(px)) {
            throw new IllegalArgumentException(name + " is NaN");
        }
        if (px < x[0] || px > x[x.length - 1]) {
            throw new IllegalArgumentException(
                    name + " out of range [" + x[0] + ", " + x[x.length - 1] + "] : " + px);
        }
    }

    private void checkY(double py, String name) {
        if (Double.isNaN(py)) {
            throw new IllegalArgumentException(name + " is NaN");
        }
        if (py < y[0] || py > y[y.length - 1]) {
            throw new IllegalArgumentException(
                    name + " out of range [" + y[0] + ", " + y[y.length - 1] + "] : " + py);
        }
    }

    /**
     * The precondition every grid interpolator in this package shares: two
     * strictly increasing abscissa vectors of at least two entries each, a
     * value array shaped to match them, and everything finite.
     */
    static void checkGrid(double[] x, double[] y, double[][] z) {
        if (x == null || y == null || z == null) {
            throw new IllegalArgumentException("null argument");
        }
        checkAxis(x, "x");
        checkAxis(y, "y");
        if (z.length != x.length) {
            throw new IllegalArgumentException("z.length != x.length : " + z.length + " != " + x.length);
        }
        for (int i = 0; i < z.length; i++) {
            if (z[i] == null) {
                throw new IllegalArgumentException("z[" + i + "] is null");
            }
            if (z[i].length != y.length) {
                throw new IllegalArgumentException(
                        "z[" + i + "].length != y.length : " + z[i].length + " != " + y.length);
            }
            for (int j = 0; j < z[i].length; j++) {
                if (Double.isNaN(z[i][j]) || Double.isInfinite(z[i][j])) {
                    throw new IllegalArgumentException("z[" + i + "][" + j + "] is not finite : " + z[i][j]);
                }
            }
        }
    }

    private static void checkAxis(double[] t, String name) {
        if (t.length < 2) {
            throw new IllegalArgumentException("at least 2 " + name + " values are needed, got " + t.length);
        }
        for (int i = 0; i < t.length; i++) {
            if (Double.isNaN(t[i]) || Double.isInfinite(t[i])) {
                throw new IllegalArgumentException(name + "[" + i + "] is not finite : " + t[i]);
            }
        }
        for (int i = 1; i < t.length; i++) {
            if (!(t[i - 1] < t[i])) {
                throw new IllegalArgumentException(name + " must be strictly increasing, but " + name + "["
                        + (i - 1) + "] = " + t[i - 1] + " and " + name + "[" + i + "] = " + t[i]);
            }
        }
    }
}
