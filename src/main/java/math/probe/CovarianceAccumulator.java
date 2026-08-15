package math.probe;

/**
 * Streaming accumulator for the mean vector and the covariance (or correlation)
 * matrix of a multivariate sample, using Welford's online algorithm.
 * <p>
 * Samples are pushed in one at a time with {@link #add(double[])} and are never
 * retained, so the memory footprint is {@code O(n^2)} in the number of features
 * regardless of how many samples pass through. This is the situation the
 * accumulator exists for: the covariance matrix of a data set that does not fit
 * in memory, or that arrives as a stream and is discarded as it goes.
 * <p>
 * <b>Why Welford and not the textbook formula.</b> The naive route accumulates
 * {@code sum(x_i * x_j)} and subtracts {@code count * mean_i * mean_j} at the
 * end. Both terms are of order {@code count * mean^2} while their difference is
 * of order {@code count * variance}, so for data whose mean is large against
 * its spread the subtraction cancels away most of the significant digits - with
 * values around {@code 1e9} and unit spread, all of them. Welford's recurrence
 * updates the co-moment directly and never forms the two large intermediates,
 * which costs one extra pass over the mean vector per sample and buys back the
 * lost accuracy.
 * <p>
 * Matrices are handed out flat and column-major ({@code C[j*n + i]}), the
 * layout used throughout {@code math.linalg}; for a symmetric matrix that is
 * the same thing as row-major. Only one triangle is updated internally and
 * mirrored on output, so the emitted matrix is <em>exactly</em> symmetric
 * rather than symmetric up to rounding - which is what the symmetric
 * eigensolver downstream expects.
 * <p>
 * Instances are not thread-safe. Independent accumulators can be filled on
 * separate threads and combined afterwards with {@link #merge(CovarianceAccumulator)},
 * which uses Chan's pairwise update and is exact in the same sense as the
 * sequential path.
 *
 * @since 1.5.1
 */
public final class CovarianceAccumulator {

    private final int n;
    private final double[] mean;    // running per-feature mean, length n
    private final double[] cm;      // co-moments, column-major, lower triangle (i >= j)
    private final double[] delta;   // scratch: x - mean before the update
    private long count;

    /**
     * Creates an empty accumulator for samples with {@code n} features.
     *
     * @param n number of features per sample, at least 1
     * @throws IllegalArgumentException if {@code n < 1}
     */
    public CovarianceAccumulator(int n) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1, got " + n);
        }
        this.n = n;
        this.mean = new double[n];
        this.cm = new double[n * n];
        this.delta = new double[n];
    }

    /**
     * Incorporates one sample.
     *
     * @param x sample with {@code n} features; not retained and not modified
     * @throws IllegalArgumentException if {@code x.length != n}
     */
    public void add(double[] x) {
        if (x.length != n) {
            throw new IllegalArgumentException("expected " + n + " features, got " + x.length);
        }
        count++;
        final double inv = 1.0 / count;
        // delta = x - mean_old, then advance the mean
        for (int j = 0; j < n; j++) {
            double d = x[j] - mean[j];
            delta[j] = d;
            mean[j] += d * inv;
        }
        // cm_ij += (x_i - mean_i_old) * (x_j - mean_j_new); lower triangle only.
        // The two triangles differ in floating point but agree in exact
        // arithmetic, so mirroring one of them (see covariance()) is as good an
        // estimator as either and is exactly symmetric.
        for (int j = 0; j < n; j++) {
            final double dNew = x[j] - mean[j];
            if (dNew == 0.0) {
                continue;
            }
            final int off = j * n;
            for (int i = j; i < n; i++) {
                cm[off + i] += delta[i] * dNew;
            }
        }
    }

    /**
     * Incorporates every row of a row-major block of samples.
     *
     * @param samples row-major data, one sample per row, {@code n} features each
     * @throws IllegalArgumentException if any row has the wrong length
     */
    public void addAll(double[][] samples) {
        for (int i = 0; i < samples.length; i++) {
            add(samples[i]);
        }
    }

    /**
     * Folds another accumulator over the same number of features into this one,
     * as if its samples had been added here (Chan's pairwise update). The other
     * accumulator is left unchanged.
     *
     * @param other accumulator to absorb
     * @throws IllegalArgumentException if {@code other} has a different feature count
     */
    public void merge(CovarianceAccumulator other) {
        if (other.n != n) {
            throw new IllegalArgumentException(
                    "feature count mismatch: " + n + " vs " + other.n);
        }
        if (other.count == 0L) {
            return;
        }
        if (count == 0L) {
            System.arraycopy(other.mean, 0, mean, 0, n);
            System.arraycopy(other.cm, 0, cm, 0, cm.length);
            count = other.count;
            return;
        }
        final double na = (double) count;
        final double nb = (double) other.count;
        final double total = na + nb;
        final double w = na * nb / total;               // weight of the mean gap
        for (int j = 0; j < n; j++) {
            delta[j] = other.mean[j] - mean[j];
        }
        for (int j = 0; j < n; j++) {
            final int off = j * n;
            final double dj = delta[j];
            for (int i = j; i < n; i++) {
                cm[off + i] += other.cm[off + i] + delta[i] * dj * w;
            }
        }
        for (int j = 0; j < n; j++) {
            mean[j] += delta[j] * (nb / total);
        }
        count += other.count;
    }

    /**
     * Number of features per sample.
     *
     * @return the feature count {@code n}
     */
    public int dimension() {
        return n;
    }

    /**
     * Number of samples seen so far.
     *
     * @return the sample count
     */
    public long count() {
        return count;
    }

    /**
     * Per-feature mean.
     *
     * @return a fresh copy of the mean vector, length {@code n} (all zero while
     *         no sample has been added)
     */
    public double[] mean() {
        return mean.clone();
    }

    /**
     * Unbiased sample covariance matrix, normalized by {@code count - 1}.
     *
     * @return a fresh, exactly symmetric {@code n x n} matrix in column-major
     *         layout
     * @throws IllegalStateException if fewer than two samples have been added
     */
    public double[] covariance() {
        if (count < 2L) {
            throw new IllegalStateException(
                    "need at least 2 samples for an unbiased covariance, have " + count);
        }
        return scaled(1.0 / (count - 1L));
    }

    /**
     * Population covariance matrix, normalized by {@code count}.
     *
     * @return a fresh, exactly symmetric {@code n x n} matrix in column-major
     *         layout
     * @throws IllegalStateException if no sample has been added
     */
    public double[] populationCovariance() {
        if (count < 1L) {
            throw new IllegalStateException("no samples");
        }
        return scaled(1.0 / count);
    }

    /**
     * Per-feature sample variance, i.e. the diagonal of {@link #covariance()}.
     *
     * @return a fresh vector of length {@code n}
     * @throws IllegalStateException if fewer than two samples have been added
     */
    public double[] variance() {
        if (count < 2L) {
            throw new IllegalStateException(
                    "need at least 2 samples for an unbiased variance, have " + count);
        }
        final double inv = 1.0 / (count - 1L);
        double[] v = new double[n];
        for (int j = 0; j < n; j++) {
            v[j] = cm[j * n + j] * inv;
        }
        return v;
    }

    /**
     * Pearson correlation matrix, i.e. {@link #covariance()} scaled to unit
     * diagonal. A feature with zero variance has no defined correlation with
     * anything; its row and column are set to zero and its diagonal entry to
     * one.
     *
     * @return a fresh, exactly symmetric {@code n x n} matrix in column-major
     *         layout
     * @throws IllegalStateException if fewer than two samples have been added
     */
    public double[] correlation() {
        double[] c = covariance();
        double[] inv = new double[n];
        for (int j = 0; j < n; j++) {
            double var = c[j * n + j];
            inv[j] = (var > 0.0) ? 1.0 / Math.sqrt(var) : 0.0;
        }
        for (int j = 0; j < n; j++) {
            final int off = j * n;
            for (int i = j; i < n; i++) {
                double r = (i == j) ? 1.0 : c[off + i] * inv[i] * inv[j];
                c[off + i] = r;
                c[i * n + j] = r;
            }
        }
        return c;
    }

    /** Mirrors the accumulated lower triangle into a full matrix, scaled by {@code f}. */
    private double[] scaled(double f) {
        double[] c = new double[n * n];
        for (int j = 0; j < n; j++) {
            final int off = j * n;
            for (int i = j; i < n; i++) {
                double v = cm[off + i] * f;
                c[off + i] = v;
                c[i * n + j] = v;
            }
        }
        return c;
    }

    @Override
    public String toString() {
        return "CovarianceAccumulator{n=" + n + ", count=" + count + "}";
    }
}
