package math.fit;

/**
 * Two-dimensional interpolation by applying a one-dimensional scheme twice:
 * along {@code x} through every grid column, then across {@code y} through the
 * results.
 * <p>
 * This is how the schemes that {@link BicubicInterpolator} cannot use get into
 * two dimensions. That class builds a genuine tensor product, which works only
 * because the natural spline is linear in the data; the rules of
 * {@link KrugerInterpolator} and {@link AkimaInterpolator} are not, so for them
 * the sweep is the construction rather than a way of computing one. Two things
 * follow, and both are measured rather than assumed:
 * <ul>
 * <li>The result <em>depends on which direction is swept first</em> for every
 * scheme except {@link Scheme#NATURAL}. Sweeping {@code y} before {@code x}
 * changed the surface by up to 15% of its own size on an unevenly spaced grid.
 * This class always sweeps {@code x} first.</li>
 * <li>Range preservation carries over from one dimension but shape does not
 * carry over as far as one might hope. Over two hundred random grids
 * {@link Scheme#KRUGER} never left the range of the data at all, while
 * {@link Scheme#AKIMA} left it by up to 184% -- far worse than the 1.1% the
 * same rule gives in one dimension -- and {@link Scheme#NATURAL} by 490%.</li>
 * </ul>
 * With {@link Scheme#NATURAL} this produces the same surface as
 * {@link BicubicInterpolator}, to about {@code 2.4e-15} of it, by an entirely
 * separate route. Prefer that class when the scheme is the natural spline: it
 * evaluates in constant time, and it can differentiate and integrate.
 *
 * @since 1.5.2
 */
public final class SuccessiveInterpolator {

    /** The one-dimensional rule applied in both directions. */
    public enum Scheme {

        /** The natural cubic spline of {@link SplineInterpolator}. */
        NATURAL,

        /** The constrained cubic spline of {@link KrugerInterpolator}. */
        KRUGER,

        /** Akima's rule, {@link AkimaInterpolator.Variant#CLASSIC}. */
        AKIMA,

        /** Akima's rule, {@link AkimaInterpolator.Variant#MODIFIED}. */
        AKIMA_MODIFIED
    }

    private SuccessiveInterpolator() {
        throw new AssertionError();
    }

    /**
     * Interpolates the given grid by sweeping {@code x} and then {@code y}.
     *
     * @param x
     *            the abscissae in the first variable, strictly increasing, at
     *            least two
     * @param y
     *            the abscissae in the second variable, strictly increasing, at
     *            least two
     * @param z
     *            the values, {@code z[i][j]} at {@code (x[i], y[j])}
     * @param scheme
     *            the one-dimensional rule to apply in both directions
     * @return the interpolating surface
     */
    public static SuccessiveSurface interpolate(double[] x, double[] y, double[][] z, Scheme scheme) {
        if (scheme == null) {
            throw new IllegalArgumentException("scheme must not be null");
        }
        PolynomialSurface.checkGrid(x, y, z);
        int nx = x.length;
        int ny = y.length;
        CubicSpline[] alongX = new CubicSpline[ny];
        double[] column = new double[nx];
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                column[i] = z[i][j];
            }
            alongX[j] = oneDimensional(x, column, scheme);
        }
        return new SuccessiveSurface(x, y, alongX, scheme);
    }

    /** Dispatches to the one-dimensional interpolator the scheme names. */
    static CubicSpline oneDimensional(double[] t, double[] w, Scheme scheme) {
        switch (scheme) {
        case NATURAL:
            return SplineInterpolator.interpolate(t, w);
        case KRUGER:
            return KrugerInterpolator.interpolate(t, w);
        case AKIMA:
            return AkimaInterpolator.interpolate(t, w, AkimaInterpolator.Variant.CLASSIC);
        default:
            return AkimaInterpolator.interpolate(t, w, AkimaInterpolator.Variant.MODIFIED);
        }
    }
}
