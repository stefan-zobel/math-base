package math.stats.gof;

import math.MathConsts;
import math.fun.DFunction;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.SymmetricJacobiEigen;
import math.solve.AdaptiveGaussKronrod;

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
 * {@code P[d <= t] = P[sum (nu_i - t) z_i^2 <= 0]}, and that is the probability
 * Imhof's inversion of the characteristic function computes.
 * <p>
 * https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic
 *
 * @since 1.5.3
 */
public final class DurbinWatson {

    /**
     * The accuracy the quadrature behind {@link #cdf} and {@link #barF} is
     * asked for, per panel.
     * <p>
     * It is an <em>absolute</em> tolerance, and both tails read the integral
     * as {@code 1/2 - I/pi} and {@code 1/2 + I/pi}, so it also bounds how far
     * down a p-value keeps its digits. Asking for less does not help: measured
     * against a closed form, tightening it from {@code 1e-11} to {@code 1e-16}
     * changed nothing at all, because what limits the answer is the
     * subtraction and not the rule.
     */
    public static final double QUADRATURE_TOLERANCE = 1.0e-13;

    /** The rule the panels are integrated with. */
    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;

    /** How deep a single panel may be bisected before its value is taken as it stands. */
    private static final int MAX_DEPTH = 14;

    /** The size of the tail of the integral that is dropped rather than integrated. */
    private static final double TRUNCATION_TOLERANCE = 1.0e-15;

    /** Bisections of the first panel taken before the error estimate is consulted. */
    private static final int FORCED_SUBDIVISIONS = 3;

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
     * See {@link #barF} for how far down a small p-value here keeps its digits:
     * the two tails are the same integral read from either side of one half, and
     * the limit that puts on them is the same.
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
        return tail(eigenvalues, t, false);
    }

    /**
     * The upper tail {@code P[d >= t]} of the null distribution, which is the
     * p-value against negative autocorrelation.
     * <p>
     * <b>Both tails are a number near one half plus or minus the same
     * integral</b>, so a small p-value here is a difference of two numbers near
     * one half -- the one thing the rest of this package is written to avoid,
     * and here it is Imhof's formula rather than a choice. Measured against an
     * exact closed form, the relative error was below {@code 1e-11} for a
     * p-value of {@code 1e-3}, below {@code 1e-9} at {@code 1e-4} and below
     * {@code 1e-5} at {@code 1e-6}; from about {@code 1e-8} downwards only the
     * order of magnitude survives, and further down the answer is clamped to
     * zero. Every conventional decision level is far above that.
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
        return tail(eigenvalues, t, true);
    }

    /**
     * Either tail, from Imhof's inversion of the characteristic function of
     * {@code Q = sum (nu_i - t) z_i^2}: {@code P[Q > 0]} is {@code 1/2 + I/pi}
     * and {@code P[Q < 0]} is {@code 1/2 - I/pi}.
     * <p>
     * The integral runs over a half line and is cut where Imhof's own bound on
     * the remainder falls below {@link #TRUNCATION_TOLERANCE}. That point can be
     * very far out when an eigenvalue sits close to {@code t}, so the range is
     * covered by panels that double in width rather than by one interval.
     */
    private static double tail(double[] eigenvalues, double t, boolean upper) {
        if (eigenvalues == null || eigenvalues.length == 0) {
            throw new IllegalArgumentException("eigenvalues must not be null or empty");
        }
        if (Double.isNaN(t)) {
            throw new IllegalArgumentException("t must not be NaN");
        }
        final double[] lambda = shifted(eigenvalues, t);
        if (lambda.length == 0) {
            // every eigenvalue sits exactly at t, so the statistic is t with
            // no spread at all and neither tail excludes anything
            return 1.0;
        }
        if (lambda.length == 1) {
            // one term left, so the sign of that eigenvalue decides everything
            return upper == (lambda[0] > 0.0) ? 1.0 : 0.0;
        }
        DFunction integrand = new DFunction() {
            @Override
            public double apply(double v) {
                if (v <= 0.0) {
                    // the limit at the origin, where sin(theta)/v is not yet
                    // an expression a machine can evaluate
                    double sum = 0.0;
                    for (int i = 0; i < lambda.length; i++) {
                        sum += lambda[i];
                    }
                    return 0.5 * sum;
                }
                double theta = 0.0;
                for (int i = 0; i < lambda.length; i++) {
                    theta += Math.atan(lambda[i] * v);
                }
                return Math.sin(0.5 * theta) * Math.exp(-Math.log(v) - logRho(lambda, v));
            }
        };

        // one panel for the interval the integrand is O(1) on, then panels
        // that double in width. A single panel over the whole range would
        // not do: the rule would place its nodes so far apart that it never
        // sees the part that carries the integral, and its own error
        // estimate would agree that there was nothing there
        double integral = AdaptiveGaussKronrod.integrate1DAdaptive(RULE, integrand, 0.0, 1.0,
                QUADRATURE_TOLERANCE, MAX_DEPTH, FORCED_SUBDIVISIONS);
        double upperLimit = truncationLimit(lambda);
        double from = 1.0;
        while (from < upperLimit) {
            double to = 2.0 * from;
            integral += AdaptiveGaussKronrod.integrate1DAdaptive(RULE, integrand, from, to,
                    QUADRATURE_TOLERANCE, MAX_DEPTH, 0);
            from = to;
        }
        return clamp(upper ? 0.5 + integral / Math.PI : 0.5 - integral / Math.PI);
    }

    /**
     * The logarithm of Imhof's {@code rho(v) = prod (1 + lambda_i^2 v^2)^(1/4)}.
     * <p>
     * Written as a sum because the product passes the largest {@code double}
     * for a few dozen eigenvalues, and each term avoids {@code log1p} of a
     * square that would overflow before its logarithm would.
     */
    private static double logRho(double[] lambda, double v) {
        double sum = 0.0;
        for (int i = 0; i < lambda.length; i++) {
            double lv = lambda[i] * v;
            double magnitude = Math.abs(lv);
            sum += magnitude > 1.0e8 ? 2.0 * Math.log(magnitude) : Math.log1p(lv * lv);
        }
        return 0.25 * sum;
    }

    /**
     * Where the rest of the integral is smaller than
     * {@link #TRUNCATION_TOLERANCE}.
     * <p>
     * Beyond {@code U} each factor of {@code rho} grows at least as
     * {@code (v/U)^(1/2)}, so the remainder is bounded by
     * {@code 2 / (pi m rho(U))} -- Imhof's own bound, which is why the half
     * line does not have to be substituted away.
     */
    private static double truncationLimit(double[] lambda) {
        double u = 1.0;
        for (int doublings = 0; doublings < 200; doublings++) {
            if (2.0 / (Math.PI * lambda.length * Math.exp(logRho(lambda, u))) < TRUNCATION_TOLERANCE) {
                return u;
            }
            u *= 2.0;
        }
        return u;
    }

    /** The eigenvalues of the shifted form, with the exact zeros left out. */
    private static double[] shifted(double[] eigenvalues, double t) {
        int kept = 0;
        for (int i = 0; i < eigenvalues.length; i++) {
            if (eigenvalues[i] != t) {
                kept++;
            }
        }
        double[] lambda = new double[kept];
        int at = 0;
        for (int i = 0; i < eigenvalues.length; i++) {
            if (eigenvalues[i] != t) {
                lambda[at++] = eigenvalues[i] - t;
            }
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

    private static double clamp(double p) {
        if (p <= 0.0) {
            return 0.0;
        }
        if (p >= 1.0) {
            return 1.0;
        }
        return p;
    }

    private DurbinWatson() {
        throw new AssertionError();
    }
}
