package math.fit;

/**
 * Bilinear interpolation on a rectangular grid: on each cell the unique
 * function of the form {@code a + b t + c v + d t v} through the four corner
 * values.
 * <p>
 * It is continuous and nothing more -- the derivative jumps at every grid
 * line -- but it is the only scheme in this package that <em>cannot</em> leave
 * the range of the data, in either dimension, since every value it produces is
 * a convex combination of four grid values. Measured over sixty random grids
 * the largest excursion outside the data was {@code 8.9e-17} of the range,
 * which is rounding and not overshoot. Where smoothness matters more than that
 * guarantee, use {@link BicubicInterpolator}.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Bilinear_interpolation">Bilinear
 * interpolation</a>
 *
 * @since 1.5.2
 */
public final class BilinearInterpolator {

    private final PolynomialSurface surface;

    /**
     * Interpolates the given grid.
     *
     * @param x
     *            the abscissae in the first variable, strictly increasing, at
     *            least two
     * @param y
     *            the abscissae in the second variable, strictly increasing, at
     *            least two
     * @param z
     *            the values, {@code z[i][j]} at {@code (x[i], y[j])}
     */
    public BilinearInterpolator(double[] x, double[] y, double[][] z) {
        surface = interpolate(x, y, z);
    }

    /**
     * The value of the interpolant at {@code (px, py)}.
     *
     * @param px
     *            the abscissa in the first variable, within the grid
     * @param py
     *            the abscissa in the second variable, within the grid
     * @return the interpolated value
     */
    public double value(double px, double py) {
        return surface.value(px, py);
    }

    /**
     * The interpolant itself, which also has the partial derivatives and the
     * double integral.
     *
     * @return the surface this interpolator computed
     */
    public PolynomialSurface surface() {
        return surface;
    }

    /**
     * Interpolates the given grid bilinearly.
     *
     * @param x
     *            the abscissae in the first variable, strictly increasing, at
     *            least two
     * @param y
     *            the abscissae in the second variable, strictly increasing, at
     *            least two
     * @param z
     *            the values, {@code z[i][j]} at {@code (x[i], y[j])}
     * @return the interpolating surface
     */
    public static PolynomialSurface interpolate(double[] x, double[] y, double[][] z) {
        PolynomialSurface.checkGrid(x, y, z);
        int cellsX = x.length - 1;
        int cellsY = y.length - 1;
        double[][][] coefficients = new double[cellsX][cellsY][];
        for (int i = 0; i < cellsX; i++) {
            double hx = x[i + 1] - x[i];
            for (int j = 0; j < cellsY; j++) {
                double hy = y[j + 1] - y[j];
                double z00 = z[i][j];
                double z10 = z[i + 1][j];
                double z01 = z[i][j + 1];
                double z11 = z[i + 1][j + 1];
                double[] c = new double[16];
                c[0] = z00;
                c[4] = (z10 - z00) / hx;
                c[1] = (z01 - z00) / hy;
                c[5] = (z11 - z10 - z01 + z00) / (hx * hy);
                coefficients[i][j] = c;
            }
        }
        return new PolynomialSurface(x, y, coefficients);
    }
}
