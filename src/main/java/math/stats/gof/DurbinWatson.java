package math.stats.gof;

import math.MathConsts;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.SymmetricJacobiEigen;

/**
 * The null distribution of the Durbin-Watson statistic
 * {@code d = sum (e_i - e_(i-1))^2 / sum e_i^2} for the residuals {@code e} of a
 * linear least squares fit with normal errors.
 * <p>
 * There is no table for this test, and there cannot be one. Writing
 * {@code e = M y} with {@code M = I - X (X'X)^-1 X'} and
 * {@code A = tridiag(-1, [1 2 ... 2 1], -1)} makes {@code d} the ratio
 * {@code z' B z / z' z} of two quadratic forms in the same standard normal
 * vector, where {@code B} is {@code A} restricted to the residual space. Its
 * eigenvalues {@code nu_1 .. nu_m}, {@code m = n - k}, depend on the design
 * matrix, so the distribution does too -- which is why the classical answer was
 * a pair of bounds that hold for every design rather than one p-value.
 * <p>
 * With the eigenvalues in hand the tail is exact:
 * {@code P[d <= t] = P[sum (nu_i - t) z_i^2 <= 0]}, which is a weighted sum of
 * chi-squares and therefore {@link WeightedChiSquare} at {@code x = 0}.
 * <p>
 * https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic
 *
 * @since 1.5.3
 */
public final class DurbinWatson {

    /**
     * The eigenvalues that determine the null distribution for a given design
     * matrix.
     * <p>
     * They are the eigenvalues of {@code A} restricted to the orthogonal
     * complement of the column space of {@code X}, in descending order, and
     * there are {@code n - k} of them. Since {@code A} is positive
     * semi-definite they are never negative, and they lie in {@code [0, 4]}.
     * <p>
     * The work is an {@code n x k} singular value decomposition, an
     * {@code O(n^2 k)} assembly and a symmetric eigendecomposition of an
     * {@code n x n} matrix. That last term is cubic in the number of
     * observations, so this is a computation to do once per design and keep,
     * not once per test -- which is why it is a separate method.
     *
     * @param designColumnMajor
     *            the design matrix {@code X}, column-major and of length
     *            {@code n * k}, element {@code (i, j)} at {@code [j * n + i]}
     * @param n
     *            the number of observations
     * @param k
     *            the number of columns of the design matrix, at least
     *            {@code 1} and less than {@code n}
     * @return the {@code n - k} eigenvalues, in descending order
     * @throws IllegalArgumentException
     *             if {@code designColumnMajor} is {@code null} or has the wrong
     *             length, if {@code k} is not in {@code [1, n)}, or if the
     *             design matrix does not have full column rank
     * @throws ArithmeticException
     *             if either decomposition failed to converge
     */
    public static double[] nullEigenvalues(double[] designColumnMajor, int n, int k) {
        if (designColumnMajor == null) {
            throw new IllegalArgumentException("designColumnMajor must not be null");
        }
        if (k < 1 || k >= n) {
            throw new IllegalArgumentException("k must lie in [1, n) : k = " + k + ", n = " + n);
        }
        if (designColumnMajor.length != (long) n * k) {
            throw new IllegalArgumentException(
                    "designColumnMajor must hold n * k = " + ((long) n * k) + " values, got "
                            + designColumnMajor.length);
        }

        // an orthonormal basis U of the column space of X. The singular values
        // come with it, so the rank question is answered by the same call
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(designColumnMajor, n, k);
        if (!svd.converged) {
            throw new ArithmeticException("the singular value decomposition of the design matrix did not converge");
        }
        double tolerance = Math.max(n, k) * MathConsts.MACH_EPS_DBL;
        if (svd.sigma[0] <= 0.0 || svd.sigma[k - 1] <= tolerance * svd.sigma[0]) {
            throw new IllegalArgumentException("the design matrix does not have full column rank: its smallest "
                    + "singular value is " + svd.sigma[k - 1] + " against a largest of " + svd.sigma[0]
                    + ", below the relative tolerance " + tolerance);
        }

        // M A M = A - U W' - W U' + U S U' with W = A U and S = U' W, which
        // costs O(n^2 k) instead of the two n x n products it is written as
        double[] u = svd.U;
        double[] w = new double[n * k];
        for (int j = 0; j < k; j++) {
            applyA(u, j * n, w, j * n, n);
        }
        double[] s = new double[k * k];
        for (int a = 0; a < k; a++) {
            for (int b = 0; b < k; b++) {
                double sum = 0.0;
                for (int i = 0; i < n; i++) {
                    sum += u[a * n + i] * w[b * n + i];
                }
                s[b * k + a] = sum;
            }
        }

        double[] mam = new double[n * n];
        fillA(mam, n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                double correction = 0.0;
                for (int a = 0; a < k; a++) {
                    correction -= u[a * n + i] * w[a * n + j];
                    correction -= w[a * n + i] * u[a * n + j];
                    double sub = 0.0;
                    for (int b = 0; b < k; b++) {
                        sub += s[b * k + a] * u[b * n + j];
                    }
                    correction += u[a * n + i] * sub;
                }
                mam[j * n + i] += correction;
            }
        }

        SymmetricJacobiEigen.Result eigen = new SymmetricJacobiEigen().decomposeInPlace(mam, n);
        if (!eigen.converged) {
            throw new ArithmeticException("the eigendecomposition of the restricted difference matrix did not "
                    + "converge");
        }

        // the k smallest eigenvalues belong to the column space of X, where the
        // restriction is zero. A is positive semi-definite, so dropping the k
        // smallest keeps exactly the m eigenvalues of the restriction -- zeros
        // of the restriction itself included
        double[] nu = new double[n - k];
        System.arraycopy(eigen.lambda, 0, nu, 0, n - k);
        for (int i = 0; i < nu.length; i++) {
            if (nu[i] < 0.0) {
                nu[i] = 0.0;
            }
        }
        return nu;
    }

    /**
     * The lower tail {@code P[d <= t]} of the null distribution, which is the
     * p-value against positive autocorrelation.
     * <p>
     * See {@link WeightedChiSquare#barF} for how far down a small p-value here
     * keeps its digits: the two tails are the same integral read from either
     * side of one half, and the limit that puts on them is the same.
     *
     * @param eigenvalues
     *            the eigenvalues from
     *            {@link #nullEigenvalues(double[], int, int)}
     * @param t
     *            the value of the statistic
     * @return {@code P[d <= t]}
     * @throws IllegalArgumentException
     *             if {@code eigenvalues} is {@code null} or empty, or if
     *             {@code t} is not a number
     */
    public static double cdf(double[] eigenvalues, double t) {
        return WeightedChiSquare.cdf(shifted(eigenvalues, t), 0.0);
    }

    /**
     * The upper tail {@code P[d >= t]} of the null distribution, which is the
     * p-value against negative autocorrelation.
     * <p>
     * <b>Both tails are a number near one half plus or minus the same
     * integral</b>, so a small p-value here is a difference of two numbers near
     * one half -- the one thing the rest of this package is written to avoid,
     * and here it is Imhof's formula rather than a choice.
     * {@link WeightedChiSquare#barF} carries the measurement of how far down
     * that leaves a p-value its digits; every conventional decision level is
     * far above the floor.
     *
     * @param eigenvalues
     *            the eigenvalues from
     *            {@link #nullEigenvalues(double[], int, int)}
     * @param t
     *            the value of the statistic
     * @return {@code P[d >= t]}
     * @throws IllegalArgumentException
     *             if {@code eigenvalues} is {@code null} or empty, or if
     *             {@code t} is not a number
     */
    public static double barF(double[] eigenvalues, double t) {
        return WeightedChiSquare.barF(shifted(eigenvalues, t), 0.0);
    }

    /**
     * The shifted eigenvalues {@code nu_i - t}, which turn the ratio into a
     * plain sum: {@code d} is at most {@code t} exactly when
     * {@code sum (nu_i - t) z_i^2} is at most zero, the denominator
     * {@code sum z_i^2} being positive with probability one.
     * <p>
     * An eigenvalue sitting exactly at {@code t} shifts to zero and is dropped
     * by {@link WeightedChiSquare}, which is why nothing is filtered here.
     */
    private static double[] shifted(double[] eigenvalues, double t) {
        if (eigenvalues == null || eigenvalues.length == 0) {
            throw new IllegalArgumentException("eigenvalues must not be null or empty");
        }
        if (Double.isNaN(t)) {
            throw new IllegalArgumentException("t must not be NaN");
        }
        double[] lambda = new double[eigenvalues.length];
        for (int i = 0; i < eigenvalues.length; i++) {
            lambda[i] = eigenvalues[i] - t;
        }
        return lambda;
    }

    /** {@code y = A x}, where {@code A} is the second difference matrix. */
    private static void applyA(double[] x, int xFrom, double[] y, int yFrom, int n) {
        if (n == 1) {
            y[yFrom] = 0.0;
            return;
        }
        y[yFrom] = x[xFrom] - x[xFrom + 1];
        for (int i = 1; i < n - 1; i++) {
            y[yFrom + i] = 2.0 * x[xFrom + i] - x[xFrom + i - 1] - x[xFrom + i + 1];
        }
        y[yFrom + n - 1] = x[xFrom + n - 1] - x[xFrom + n - 2];
    }

    /** {@code A} itself, column-major and dense. */
    private static void fillA(double[] a, int n) {
        for (int i = 0; i < n; i++) {
            a[i * n + i] = (i == 0 || i == n - 1) ? 1.0 : 2.0;
        }
        for (int i = 0; i < n - 1; i++) {
            a[(i + 1) * n + i] = -1.0;
            a[i * n + i + 1] = -1.0;
        }
    }

    private DurbinWatson() {
        throw new AssertionError();
    }
}
