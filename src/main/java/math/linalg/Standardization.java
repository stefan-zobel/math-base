package math.linalg;

/**
 * Centers each column of a table of observations on its mean and divides it by
 * its standard deviation, the transform usually called the z-score.
 * <p>
 * The rows are observations and the columns are variables, the layout
 * {@link CovariancePCA#transform(double[][])} and
 * {@link math.probe.CovarianceAccumulator#addAll(double[][])} already use. It
 * is a <em>fitted</em> transform: the mean and the scale are taken from the
 * data handed to {@link #of(double[][])} and kept, so the same transform can be
 * applied to rows it was not fitted on -- a hold-out set, or one observation
 * arriving later -- and undone again with {@link #inverse(double[][])}.
 * Re-standardizing each batch on its own would give each of them a different
 * scale, which is the mistake this class exists to make hard.
 * <p>
 * Standardizing is what makes a scale-dependent method answer a question about
 * shape rather than about units. A principal component analysis of a table
 * whose columns are millimeters and kilograms is dominated by whichever column
 * happens to carry the larger numbers; after this transform every column
 * contributes its variance rather than its unit.
 * <p>
 * {@link #of(double[][])} divides by the sample standard deviation
 * {@code sqrt(sum (x - xBar)^2 / (n - 1))} and {@link #ofPopulation(double[][])}
 * by {@code sqrt(sum (x - xBar)^2 / n)}. Which of the two is wanted is a real
 * question rather than a detail, so it is spelled in the method name; the same
 * pair is named {@code variance()} and {@code populationCovariance()} on
 * {@code CovarianceAccumulator}.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Standard_score">Standard
 *      score</a>
 * @since 1.5.2
 */
public final class Standardization {

    /** Column means of the data this was fitted on. */
    private final double[] mean;

    /** Column standard deviations, all strictly positive. */
    private final double[] scale;

    private Standardization(double[] mean, double[] scale) {
        this.mean = mean;
        this.scale = scale;
    }

    /**
     * Fits the transform, dividing by the sample standard deviation.
     *
     * @param samples
     *            the observations, one per row, all rows of the same length;
     *            not modified and not retained
     * @return the fitted transform
     * @throws IllegalArgumentException
     *             if {@code samples} is {@code null}, holds fewer than two
     *             rows, has no columns, has rows of differing length, has a
     *             column that is constant, or holds a value that is not finite
     */
    public static Standardization of(double[][] samples) {
        return fit(samples, false);
    }

    /**
     * Fits the transform, dividing by the population standard deviation.
     * <p>
     * This is the right divisor when the rows are the whole population rather
     * than a sample drawn from one. It differs from {@link #of(double[][])} by
     * the factor {@code sqrt((n - 1) / n)} in every column, so the two agree in
     * the limit and differ by about half a per cent at {@code n = 100}.
     *
     * @param samples
     *            the observations, one per row, all rows of the same length;
     *            not modified and not retained
     * @return the fitted transform
     * @throws IllegalArgumentException
     *             if {@code samples} is {@code null}, holds no rows, has no
     *             columns, has rows of differing length, has a column that is
     *             constant, or holds a value that is not finite
     */
    public static Standardization ofPopulation(double[][] samples) {
        return fit(samples, true);
    }

    private static Standardization fit(double[][] samples, boolean population) {
        if (samples == null) {
            throw new IllegalArgumentException("samples is null");
        }
        int n = samples.length;
        int minimum = population ? 1 : 2;
        if (n < minimum) {
            throw new IllegalArgumentException("need at least " + minimum + " row"
                    + (minimum == 1 ? "" : "s") + " to fit "
                    + (population ? "a population" : "a sample") + " standard deviation, have " + n);
        }
        if (samples[0] == null) {
            throw new IllegalArgumentException("row 0 is null");
        }
        int p = samples[0].length;
        if (p == 0) {
            throw new IllegalArgumentException("samples has no columns");
        }
        for (int i = 1; i < n; i++) {
            if (samples[i] == null) {
                throw new IllegalArgumentException("row " + i + " is null");
            }
            if (samples[i].length != p) {
                throw new IllegalArgumentException(
                        "row " + i + " has " + samples[i].length + " columns, row 0 has " + p);
            }
        }

        double[] mean = new double[p];
        double[] scale = new double[p];
        double divisor = population ? n : (n - 1);
        // two passes per column, the mean first and the deviations from it
        // afterwards, rather than one pass over the squares: the one-pass form
        // subtracts two large numbers to get a small one and loses the answer
        // when the mean is large beside the spread
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += samples[i][j];
            }
            double m = sum / n;
            double ss = 0.0;
            for (int i = 0; i < n; i++) {
                double deviation = samples[i][j] - m;
                ss += deviation * deviation;
            }
            double s = Math.sqrt(ss / divisor);
            // a non-finite entry poisons the mean and then the deviations, so
            // it arrives here as a NaN scale rather than as a zero one; saying
            // "constant" about it would name the wrong problem
            if (Double.isNaN(s)) {
                throw new IllegalArgumentException("column " + j + " holds a value that is not"
                        + " finite, so it has neither a mean nor a scale");
            }
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("column " + j + " is constant, so it has no scale"
                        + " to standardize to; every value in it is " + samples[0][j]);
            }
            mean[j] = m;
            scale[j] = s;
        }
        return new Standardization(mean, scale);
    }

    /**
     * Centers and scales one variable, the whole of it, and hands back the
     * z-scores.
     * <p>
     * The one-column case of {@link #of(double[][])}, for a caller who has a
     * single variable rather than a table and does not need the transform
     * afterwards. It divides by the sample standard deviation.
     *
     * @param v
     *            the values; not modified
     * @return a new array of the same length holding
     *         {@code (v_i - mean) / sd}
     * @throws IllegalArgumentException
     *             if {@code v} is {@code null}, holds fewer than two values, is
     *             constant, or holds a value that is not finite
     */
    public static double[] standardize(double[] v) {
        if (v == null) {
            throw new IllegalArgumentException("v is null");
        }
        if (v.length < 2) {
            throw new IllegalArgumentException(
                    "need at least 2 values to fit a sample standard deviation, have " + v.length);
        }
        int n = v.length;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += v[i];
        }
        double mean = sum / n;
        double ss = 0.0;
        for (int i = 0; i < n; i++) {
            double deviation = v[i] - mean;
            ss += deviation * deviation;
        }
        double sd = Math.sqrt(ss / (n - 1));
        if (Double.isNaN(sd)) {
            throw new IllegalArgumentException(
                    "the values hold one that is not finite, so they have neither a mean nor a scale");
        }
        if (!(sd > 0.0)) {
            throw new IllegalArgumentException(
                    "the values are constant, so they have no scale to standardize to; every one of"
                            + " them is " + v[0]);
        }
        double[] z = new double[n];
        for (int i = 0; i < n; i++) {
            z[i] = (v[i] - mean) / sd;
        }
        return z;
    }

    /**
     * The number of columns this transform was fitted on.
     *
     * @return the number of variables
     */
    public int dimension() {
        return mean.length;
    }

    /**
     * The column means that are subtracted.
     *
     * @return a new array of length {@link #dimension()}
     */
    public double[] mean() {
        return mean.clone();
    }

    /**
     * The column standard deviations that are divided by, at whichever divisor
     * the transform was fitted with.
     *
     * @return a new array of length {@link #dimension()}, all entries strictly
     *         positive
     */
    public double[] scale() {
        return scale.clone();
    }

    /**
     * Applies the transform to a block of observations.
     *
     * @param samples
     *            the observations, one per row, each of length
     *            {@link #dimension()}; not modified
     * @return a new block of the same shape
     * @throws IllegalArgumentException
     *             if {@code samples} is {@code null} or holds a row that is
     *             {@code null} or of the wrong length
     */
    public double[][] transform(double[][] samples) {
        return map(samples, true);
    }

    /**
     * Undoes the transform on a block of standardized observations.
     *
     * @param standardized
     *            the standardized observations, one per row, each of length
     *            {@link #dimension()}; not modified
     * @return a new block of the same shape, in the original units
     * @throws IllegalArgumentException
     *             if {@code standardized} is {@code null} or holds a row that
     *             is {@code null} or of the wrong length
     */
    public double[][] inverse(double[][] standardized) {
        return map(standardized, false);
    }

    /**
     * Applies the transform to one observation.
     *
     * @param sample
     *            one observation of length {@link #dimension()}; not modified
     * @return a new array of the same length
     * @throws IllegalArgumentException
     *             if {@code sample} is {@code null} or of the wrong length
     */
    public double[] transformRow(double[] sample) {
        return mapRow(sample, true);
    }

    /**
     * Undoes the transform on one standardized observation.
     *
     * @param standardized
     *            one standardized observation of length {@link #dimension()};
     *            not modified
     * @return a new array of the same length, in the original units
     * @throws IllegalArgumentException
     *             if {@code standardized} is {@code null} or of the wrong
     *             length
     */
    public double[] inverseRow(double[] standardized) {
        return mapRow(standardized, false);
    }

    private double[][] map(double[][] rows, boolean forward) {
        if (rows == null) {
            throw new IllegalArgumentException("samples is null");
        }
        double[][] out = new double[rows.length][];
        for (int i = 0; i < rows.length; i++) {
            if (rows[i] == null) {
                throw new IllegalArgumentException("row " + i + " is null");
            }
            if (rows[i].length != mean.length) {
                throw new IllegalArgumentException("row " + i + " has " + rows[i].length
                        + " columns, this transform was fitted on " + mean.length);
            }
            out[i] = mapRow(rows[i], forward);
        }
        return out;
    }

    private double[] mapRow(double[] row, boolean forward) {
        if (row == null) {
            throw new IllegalArgumentException("sample is null");
        }
        int p = mean.length;
        if (row.length != p) {
            throw new IllegalArgumentException("sample has " + row.length
                    + " columns, this transform was fitted on " + p);
        }
        double[] out = new double[p];
        for (int j = 0; j < p; j++) {
            out[j] = forward ? (row[j] - mean[j]) / scale[j] : row[j] * scale[j] + mean[j];
        }
        return out;
    }
}
