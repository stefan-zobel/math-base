package math.linalg;

/**
 * Centered and scaled copy of a design matrix and its response, the form every
 * penalized estimator in this package fits on.
 * <p>
 * Neither the ridge nor the lasso penalty is scale invariant, so the columns
 * have to be brought to a common scale before the penalty means the same thing
 * for each of them. Each column is centered on its mean and divided by its root
 * mean square; the response is centered. The intercept is therefore not part of
 * the fit and is recovered afterwards as {@code yBar - sum(beta_j * xBar_j)},
 * which keeps it out of the penalty where it belongs.
 * <p>
 * After the transform every column satisfies {@code sum_i x_ij^2 == n}, which
 * is what lets the coordinate update in {@link CoordinateDescent} divide by
 * {@code n} instead of by a per column norm.
 * <p>
 * That is also why this class is not the general one and is not public. It
 * divides by the root mean square {@code sqrt(ss/n)} rather than by the sample
 * standard deviation {@code sqrt(ss/(n-1))}, and it does so because
 * {@link CoordinateDescent} requires it, not because that is the customary
 * scale: the divisor is coupled to an inner loop one class away. For
 * standardizing a table of observations, use {@link Standardization}, which
 * divides by the sample standard deviation, keeps no response beside the
 * design, and can be applied again to data it was not fitted on.
 */
final class ScaledDesign {

    /** Number of rows actually used. */
    final int n;
    /** Number of columns. */
    final int p;
    /** Centered and scaled design, {@code n x p}, column-major. */
    final double[] x;
    /** Centered response, length {@code n}. */
    final double[] y;
    /** Column means of the original design, length {@code p}. */
    final double[] xBar;
    /** Column scales, {@code sqrt(sum (x_ij - xBar_j)^2 / n)}, length {@code p}. */
    final double[] scale;
    /** Mean of the original response. */
    final double yBar;

    private ScaledDesign(int n, int p, double[] x, double[] y, double[] xBar, double[] scale, double yBar) {
        this.n = n;
        this.p = p;
        this.x = x;
        this.y = y;
        this.xBar = xBar;
        this.scale = scale;
        this.yBar = yBar;
    }

    /**
     * Transforms the given design and response.
     *
     * @param xs
     *            the design, {@code rows x p}, column-major; not modified
     * @param ys
     *            the response, length {@code rows}; not modified
     * @param rows
     *            the number of rows {@code xs} and {@code ys} hold
     * @param p
     *            the number of columns of {@code xs}
     * @param use
     *            the row indices to use, in the order given, or {@code null}
     *            for all of them
     * @return the transform, holding its own arrays
     * @throws IllegalArgumentException
     *             if a column is constant over the rows used
     */
    static ScaledDesign of(double[] xs, double[] ys, int rows, int p, int[] use) {
        int n = (use == null) ? rows : use.length;
        double[] x = new double[n * p];
        double[] y = new double[n];
        double[] xBar = new double[p];
        double[] scale = new double[p];

        double yBar = 0.0;
        for (int i = 0; i < n; i++) {
            yBar += ys[(use == null) ? i : use[i]];
        }
        yBar /= n;
        for (int i = 0; i < n; i++) {
            y[i] = ys[(use == null) ? i : use[i]] - yBar;
        }

        for (int j = 0; j < p; j++) {
            int src = j * rows;
            int dst = j * n;
            double mean = 0.0;
            for (int i = 0; i < n; i++) {
                mean += xs[src + ((use == null) ? i : use[i])];
            }
            mean /= n;
            double ss = 0.0;
            for (int i = 0; i < n; i++) {
                double centered = xs[src + ((use == null) ? i : use[i])] - mean;
                x[dst + i] = centered;
                ss += centered * centered;
            }
            double s = Math.sqrt(ss / n);
            if (s == 0.0) {
                throw new IllegalArgumentException(
                        "column " + j + " of X is constant; the intercept is fitted separately, "
                                + "so a constant column carries no information");
            }
            xBar[j] = mean;
            scale[j] = s;
            for (int i = 0; i < n; i++) {
                x[dst + i] /= s;
            }
        }

        return new ScaledDesign(n, p, x, y, xBar, scale, yBar);
    }

    /**
     * Turns coefficients of the standardized problem back into coefficients of
     * the original one.
     *
     * @param betaScaled
     *            coefficients in the standardized scale, length {@code p}
     * @return a new array of coefficients in the original scale
     */
    double[] unscale(double[] betaScaled) {
        double[] beta = new double[p];
        for (int j = 0; j < p; j++) {
            beta[j] = betaScaled[j] / scale[j];
        }
        return beta;
    }

    /**
     * The intercept belonging to coefficients in the original scale.
     *
     * @param beta
     *            coefficients in the original scale, length {@code p}
     * @return {@code yBar - sum(beta_j * xBar_j)}
     */
    double intercept(double[] beta) {
        double intercept = yBar;
        for (int j = 0; j < p; j++) {
            intercept -= beta[j] * xBar[j];
        }
        return intercept;
    }
}
