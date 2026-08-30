package math.ode;

import java.util.Arrays;

/**
 * The coefficients of a Rosenbrock method, and no code that uses them.
 * <p>
 * A Rosenbrock method is what is left of a diagonally implicit Runge-Kutta
 * method when the nonlinear solve in each stage is replaced by a single
 * Newton step taken with the Jacobian of the field. What remains is
 * <em>linearly</em> implicit: every stage is one linear system in the same
 * matrix, so a step costs one factorization and one back substitution per
 * stage, and no iteration converges or fails to. With {@code J = df/dy} and
 * {@code dT = df/dt} taken once per step at the point the step starts from,
 *
 * <pre>
 * W = I / (gamma h) - J
 *
 * W k[1] = f(t, y) + h d[1] dT
 * W k[i] = f(t + alpha[i] h, y + sum a[i][j] k[j])
 *          + sum (c[i][j] / h) k[j] + h d[i] dT
 *
 * y[n+1] = y[n] + sum b[i] k[i]
 * error  = sum bError[i] k[i]
 * </pre>
 *
 * <b>The abscissae are called {@code alpha} here</b> because {@code c} is
 * needed for the coupling matrix, which is the name it carries in the
 * literature and in every implementation. They are genuinely independent
 * coefficients: {@code alpha[i]} is not the row sum of {@code a[i]}, because
 * the stage vectors {@code k} are not increments of the solution.
 * <p>
 * <b>What is not checked here.</b> The constructor validates shapes and
 * finiteness and nothing else. The order conditions are awkward to write down
 * and easy to get wrong, and there is a better check: on {@code y' = lambda y}
 * one step is exactly {@code R(z) y} with {@code z = h lambda} and {@code R}
 * rational, so the order, the error constant, the order of the embedded
 * estimate and L-stability all fall out of one number the machinery computes
 * itself. That is what the tests do, and what polices a transcribed table.
 * <p>
 * Instances are immutable and shareable; every accessor that hands out an
 * array hands out a copy.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Rosenbrock_methods">Wikipedia
 * Rosenbrock methods</a>.
 *
 * @since 1.5.3
 */
public final class RosenbrockTableau {

    /**
     * The two stage L-stable method of order two, with an embedded estimate of
     * order one.
     * <p>
     * Every coefficient is a closed form in the square root of two --
     * {@code gamma = 1 + 1/sqrt(2)}, {@code a[2][1] = 2 - sqrt(2)},
     * {@code c[2][1] = -2 (2 - sqrt(2))} -- so this table carries no
     * transcribed digit at all. It is the cheapest stiff method here and the
     * one the machinery is proved on before a table is copied from anywhere.
     */
    public static final RosenbrockTableau ROS2 = ros2();

    /**
     * The six stage L-stable method of order four with an embedded estimate of
     * order three, stiffly accurate, and with a continuous extension of its
     * own. This is the workhorse of the package, and the method a stiff problem
     * should be given unless there is a reason not to.
     * <p>
     * Its coefficients are transcribed, from Hairer and Wanner by way of the
     * table in {@code OrdinaryDiffEq.jl}; only the weights are derived here,
     * because a stiffly accurate method's are the last row of its stage matrix
     * with a one appended and its error estimate is the last stage vector
     * alone.
     */
    public static final RosenbrockTableau RODAS4 = rodas4();

    /**
     * The same shape as {@link #RODAS4} with coefficients chosen against the
     * order reduction that afflicts it on a parabolic problem or a differential
     * algebraic one of index one.
     * <p>
     * On an ordinary stiff system there is little between them, so this is the
     * one to reach for when a problem comes out of a discretized partial
     * differential equation and the measured order falls short of four.
     */
    public static final RosenbrockTableau RODAS4P = rodas4p();

    private final String name;
    private final int order;
    private final int embeddedOrder;
    private final double gamma;
    private final double[] alpha;
    private final double[][] a;
    private final double[][] c;
    private final double[] d;
    private final double[] b;
    private final double[] bError;
    private final double[][] dense;

    /**
     * A tableau from its coefficients.
     *
     * @param name
     *            how the method is known, for {@link #toString()}
     * @param order
     *            the order of the propagating solution, at least one
     * @param embeddedOrder
     *            the order of the embedded solution, or zero if there is none
     * @param gamma
     *            the diagonal coefficient, finite and not zero
     * @param alpha
     *            the abscissae, one per stage, the first of them zero (copied)
     * @param a
     *            the stage matrix, strictly lower triangular and ragged: row
     *            {@code i} holds {@code i} entries (copied)
     * @param c
     *            the coupling matrix, of the same shape as {@code a} (copied)
     * @param d
     *            the coefficients of the time derivative, one per stage
     *            (copied)
     * @param b
     *            the weights of the propagating solution, one per stage
     *            (copied)
     * @param bError
     *            the weights of the embedded solution, one per stage, or
     *            {@code null} if the method has no error estimate (copied)
     * @param dense
     *            the rows of the interpolation matrix, each holding one entry
     *            per stage, or {@code null} if the method has no continuous
     *            extension of its own (copied)
     * @throws IllegalArgumentException
     *             if an argument is {@code null} where it may not be, a shape
     *             does not fit, or a coefficient is not finite
     */
    public RosenbrockTableau(String name, int order, int embeddedOrder, double gamma, double[] alpha,
            double[][] a, double[][] c, double[] d, double[] b, double[] bError, double[][] dense) {
        if (name == null || name.isEmpty()) {
            throw new IllegalArgumentException("name must not be null or empty");
        }
        if (order < 1) {
            throw new IllegalArgumentException("order must be at least 1, not " + order);
        }
        if (embeddedOrder < 0) {
            throw new IllegalArgumentException("embeddedOrder must not be negative, but is " + embeddedOrder);
        }
        if (!isFinite(gamma) || gamma == 0.0) {
            throw new IllegalArgumentException("gamma must be finite and not zero, but is " + gamma);
        }
        if (alpha == null || alpha.length < 1) {
            throw new IllegalArgumentException("alpha must hold at least one abscissa");
        }
        int s = alpha.length;
        checkFinite(alpha, "alpha");
        if (alpha[0] != 0.0) {
            throw new IllegalArgumentException("alpha[0] must be zero, but is " + alpha[0]);
        }
        checkTriangular(a, s, "a");
        checkTriangular(c, s, "c");
        checkVector(d, s, "d");
        checkVector(b, s, "b");
        if (bError != null) {
            checkVector(bError, s, "bError");
        }
        if (bError == null && embeddedOrder != 0) {
            throw new IllegalArgumentException("embeddedOrder is " + embeddedOrder + " without any bError");
        }
        if (bError != null && embeddedOrder == 0) {
            throw new IllegalArgumentException("bError was given without an embeddedOrder");
        }
        if (dense != null) {
            if (dense.length < 1) {
                throw new IllegalArgumentException("dense must hold at least one row or be null");
            }
            for (int j = 0; j < dense.length; ++j) {
                checkVector(dense[j], s, "dense[" + j + "]");
            }
        }
        this.name = name;
        this.order = order;
        this.embeddedOrder = embeddedOrder;
        this.gamma = gamma;
        this.alpha = alpha.clone();
        this.a = copy(a);
        this.c = copy(c);
        this.d = d.clone();
        this.b = b.clone();
        this.bError = (bError == null) ? null : bError.clone();
        this.dense = (dense == null) ? null : copy(dense);
    }

    /**
     * How the method is known.
     *
     * @return the name given to the constructor
     */
    public String name() {
        return name;
    }

    /**
     * The order of the propagating solution.
     *
     * @return the order, at least one
     */
    public int order() {
        return order;
    }

    /**
     * The order of the embedded solution, one below the propagating one for
     * every method here.
     *
     * @return the embedded order, or zero if there is no error estimate
     */
    public int embeddedOrder() {
        return embeddedOrder;
    }

    /**
     * The number of stages, which is the number of linear systems solved per
     * step and, for every method here, also the number of times the field is
     * evaluated.
     *
     * @return the number of stages, at least one
     */
    public int stages() {
        return alpha.length;
    }

    /**
     * The diagonal coefficient, the {@code gamma} of
     * {@code W = I / (gamma h) - J}.
     *
     * @return the diagonal coefficient, finite and not zero
     */
    public double gamma() {
        return gamma;
    }

    /**
     * Whether the method carries an embedded solution to estimate its own
     * error against, which is what an adaptive step size needs.
     *
     * @return {@code true} if {@link #bError()} is not {@code null}
     */
    public boolean hasErrorEstimate() {
        return bError != null;
    }

    /**
     * Whether the method carries a continuous extension of its own, built out
     * of the stage vectors it has already computed and therefore free.
     *
     * @return {@code true} if {@link #dense()} is not {@code null}
     */
    public boolean hasDenseOutput() {
        return dense != null;
    }

    /**
     * The number of rows of the interpolation matrix, which is the degree of
     * the continuous extension beyond the straight line through the two ends.
     *
     * @return the number of dense rows, or zero if there is no continuous
     *         extension
     */
    public int denseRows() {
        return (dense == null) ? 0 : dense.length;
    }

    /**
     * Whether the propagating solution is the last stage plus its own stage
     * vector, that is whether {@code b} is the last row of {@code a} with a one
     * appended.
     * <p>
     * A stiffly accurate method reproduces the exact solution of a very stiff
     * component in one step instead of merely damping it, which is what makes
     * it usable on a problem whose fast modes are already at rest.
     *
     * @return {@code true} if the method is stiffly accurate
     */
    public boolean isStifflyAccurate() {
        int s = alpha.length;
        if (b[s - 1] != 1.0) {
            return false;
        }
        for (int j = 0; j < s - 1; ++j) {
            if (b[j] != a[s - 1][j]) {
                return false;
            }
        }
        return true;
    }

    /**
     * The abscissae, one per stage.
     *
     * @return a fresh <code>double[]</code> of length {@link #stages()}
     */
    public double[] alpha() {
        return alpha.clone();
    }

    /**
     * The stage matrix, ragged and strictly lower triangular.
     *
     * @return a fresh ragged <code>double[][]</code> whose row {@code i} has
     *         length {@code i}
     */
    public double[][] a() {
        return copy(a);
    }

    /**
     * The coupling matrix, ragged and strictly lower triangular. Its entries
     * are divided by the step size where they are used.
     *
     * @return a fresh ragged <code>double[][]</code> whose row {@code i} has
     *         length {@code i}
     */
    public double[][] c() {
        return copy(c);
    }

    /**
     * The coefficients of the time derivative, one per stage. They are all
     * zero past the point where the method stops needing it.
     *
     * @return a fresh <code>double[]</code> of length {@link #stages()}
     */
    public double[] d() {
        return d.clone();
    }

    /**
     * The weights of the propagating solution.
     *
     * @return a fresh <code>double[]</code> of length {@link #stages()}
     */
    public double[] b() {
        return b.clone();
    }

    /**
     * The weights of the embedded solution, which combine into the error
     * estimate directly rather than as a difference of two solutions.
     *
     * @return a fresh <code>double[]</code> of length {@link #stages()}, or
     *         {@code null} if the method has no error estimate
     */
    public double[] bError() {
        return (bError == null) ? null : bError.clone();
    }

    /**
     * The rows of the interpolation matrix, each combining the stage vectors
     * into one coefficient of the continuous extension.
     *
     * @return a fresh <code>double[][]</code> with {@link #denseRows()} rows of
     *         length {@link #stages()}, or {@code null} if the method has no
     *         continuous extension of its own
     */
    public double[][] dense() {
        return (dense == null) ? null : copy(dense);
    }

    /**
     * The name, the order and the size of the method.
     */
    @Override
    public String toString() {
        return name + "[order " + order + (hasErrorEstimate() ? "(" + embeddedOrder + ")" : "") + ", "
                + alpha.length + " stages]";
    }

    /**
     * Two tableaus are equal when every coefficient is, the name included.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (!(other instanceof RosenbrockTableau)) {
            return false;
        }
        RosenbrockTableau that = (RosenbrockTableau) other;
        return name.equals(that.name) && order == that.order && embeddedOrder == that.embeddedOrder
                && Double.compare(gamma, that.gamma) == 0 && Arrays.equals(alpha, that.alpha)
                && Arrays.deepEquals(a, that.a) && Arrays.deepEquals(c, that.c) && Arrays.equals(d, that.d)
                && Arrays.equals(b, that.b) && Arrays.equals(bError, that.bError)
                && Arrays.deepEquals(dense, that.dense);
    }

    /**
     * A hash over the same coefficients {@link #equals(Object)} compares.
     */
    @Override
    public int hashCode() {
        int hash = name.hashCode();
        hash = 31 * hash + order;
        hash = 31 * hash + embeddedOrder;
        hash = 31 * hash + Double.valueOf(gamma).hashCode();
        hash = 31 * hash + Arrays.hashCode(alpha);
        hash = 31 * hash + Arrays.deepHashCode(a);
        hash = 31 * hash + Arrays.deepHashCode(c);
        hash = 31 * hash + Arrays.hashCode(d);
        hash = 31 * hash + Arrays.hashCode(b);
        hash = 31 * hash + Arrays.hashCode(bError);
        hash = 31 * hash + Arrays.deepHashCode(dense);
        return hash;
    }

    /**
     * The two stage method, every entry written as what it is rather than as
     * the decimal it comes to.
     */
    private static RosenbrockTableau ros2() {
        double root = Math.sqrt(2.0);
        double gamma = 1.0 + 1.0 / root;
        double step = 2.0 - root;
        return new RosenbrockTableau("ROS2", 2, 1, gamma, new double[] { 0.0, 1.0 },
                new double[][] { {}, { step } }, new double[][] { {}, { -2.0 * step } },
                new double[] { gamma, -gamma }, new double[] { 1.5 * step, 0.5 * step },
                new double[] { 0.5 * step, 0.5 * step }, null);
    }

    /**
     * Hairer and Wanner's method, six stages of order four. The published
     * table is the abscissae, the two matrices, the coefficients of the time
     * derivative and the two rows of the interpolation matrix; the weights
     * follow from the structure and are not copied.
     */
    private static RosenbrockTableau rodas4() {
        double[] alpha = { 0.0, 0.386, 0.21, 0.63, 1.0, 1.0 };
        double[][] a = { {}, { 1.544 }, { 0.9466785280815826, 0.2557011698983284 },
                { 3.314825187068521, 2.896124015972201, 0.9986419139977817 },
                { 1.221224509226641, 6.019134481288629, 12.53708332932087, -0.687886036105895 },
                { 1.221224509226641, 6.019134481288629, 12.53708332932087, -0.687886036105895, 1.0 } };
        double[][] c = { {}, { -5.6688 }, { -2.430093356833875, -0.2063599157091915 },
                { -0.1073529058151375, -9.594562251023355, -20.47028614809616 },
                { 7.496443313967647, -10.24680431464352, -33.99990352819905, 11.7089089320616 },
                { 8.083246795921522, -7.981132988064893, -31.52159432874371, 16.31930543123136,
                        -6.058818238834054 } };
        double[] d = { 0.25, -0.1043, 0.1035, -0.0362, 0.0, 0.0 };
        double[][] dense = {
                { 10.12623508344586, -7.487995877610167, -34.80091861555747, -7.992771707568823,
                        1.025137723295662, 0.0 },
                { -0.6762803392801253, 6.087714651680015, 16.43084320892478, 24.76722511418386,
                        -6.594389125716872, 0.0 } };
        return stifflyAccurate("RODAS4", 4, 3, 0.25, alpha, a, c, d, dense);
    }

    /** Steinebach's coefficients in the same six stage shape. */
    private static RosenbrockTableau rodas4p() {
        double[] alpha = { 0.0, 0.75, 0.21, 0.63, 1.0, 1.0 };
        double[][] a = { {}, { 3.0 }, { 1.831036793486759, 0.4955183967433795 },
                { 2.304376582692669, -0.05249275245743001, -1.176798761832782 },
                { -7.170454962423024, -4.741636671481785, -16.31002631330971, -1.062004044111401 },
                { -7.170454962423024, -4.741636671481785, -16.31002631330971, -1.062004044111401,
                        1.0 } };
        double[][] c = { {}, { -12.0 }, { -8.791795173947035, -2.207865586973518 },
                { 10.81793056857153, 6.780270611428266, 19.5348594464241 },
                { 34.19095006749676, 15.49671153725963, 54.7476087596413, 14.16005392148534 },
                { 34.62605830930532, 15.30084976114473, 56.99955578662667, 18.40807009793095,
                        -5.714285714285717 } };
        double[] d = { 0.25, -0.5, -0.023504, -0.0362, 0.0, 0.0 };
        double[][] dense = {
                { 25.09876703708589, 11.62013104361867, 28.49148307714626, -5.664021568594133, 0.0,
                        0.0 },
                { 1.638054557396973, -0.7373619806678748, 8.47791821923899, 15.9925314877952,
                        -1.882352941176471, 0.0 } };
        return stifflyAccurate("RODAS4P", 4, 3, 0.25, alpha, a, c, d, dense);
    }

    /**
     * A tableau whose solution is its last stage plus that stage's own vector,
     * which fixes both sets of weights and leaves nothing to transcribe.
     */
    private static RosenbrockTableau stifflyAccurate(String name, int order, int embeddedOrder,
            double gamma, double[] alpha, double[][] a, double[][] c, double[] d, double[][] dense) {
        int s = alpha.length;
        double[] b = new double[s];
        System.arraycopy(a[s - 1], 0, b, 0, s - 1);
        b[s - 1] = 1.0;
        double[] bError = new double[s];
        bError[s - 1] = 1.0;
        return new RosenbrockTableau(name, order, embeddedOrder, gamma, alpha, a, c, d, b, bError, dense);
    }

    private static void checkTriangular(double[][] m, int s, String name) {
        if (m == null || m.length != s) {
            throw new IllegalArgumentException(name + " must have " + s + " rows");
        }
        for (int i = 0; i < s; ++i) {
            if (m[i] == null || m[i].length != i) {
                throw new IllegalArgumentException(name + "[" + i + "] must hold " + i + " entries");
            }
            checkFinite(m[i], name + "[" + i + "]");
        }
    }

    private static void checkVector(double[] v, int s, String name) {
        if (v == null || v.length != s) {
            throw new IllegalArgumentException(name + " must hold " + s + " entries");
        }
        checkFinite(v, name);
    }

    private static void checkFinite(double[] v, String name) {
        for (int i = 0; i < v.length; ++i) {
            if (!isFinite(v[i])) {
                throw new IllegalArgumentException(name + "[" + i + "] is not finite: " + v[i]);
            }
        }
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static double[][] copy(double[][] m) {
        double[][] out = new double[m.length][];
        for (int i = 0; i < m.length; ++i) {
            out[i] = m[i].clone();
        }
        return out;
    }
}
