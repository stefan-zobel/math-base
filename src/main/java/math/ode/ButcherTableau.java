package math.ode;

import java.util.Arrays;

/**
 * The coefficients of an explicit Runge-Kutta method, as an immutable value.
 * <p>
 * A step of length {@code h} from <code>(t, y)</code> evaluates the field at
 * {@code s} stages,
 * <p>
 * <code>k[i] = f(t + c[i] h, y + h sum(a[i][j] k[j], j &lt; i))</code>
 * <p>
 * and advances the solution by <code>y + h sum(b[i] k[i])</code>. An embedded
 * pair carries a second weight vector <code>b*</code> of a lower order over the
 * same stages, and the difference of the two is an estimate of the error of the
 * step -- which is what makes adaptive step size control possible without
 * halving the step and comparing.
 * <p>
 * <b>This is data, not an algorithm.</b> {@link ExplicitRungeKutta} executes
 * whichever tableau it is handed, so a further method is a set of numbers
 * rather than a class, and the order conditions that decide whether those
 * numbers are the ones intended can be checked generically. That check is worth
 * naming, because the failure it prevents is silent: a mistyped coefficient
 * does not throw, it lowers the order of the method, and the symptom appears
 * much later as an accuracy that does not improve the way it should.
 * <p>
 * The two tableaux the library ships are {@link #CLASSIC_RK4} and
 * {@link #DORMAND_PRINCE_45}.
 * <p>
 * <b>Stage counts come in two.</b> {@link #stages()} is how many evaluations
 * advance the solution; {@link #denseStages()} is how many are needed when a
 * value strictly inside the step is wanted. For the tableaux here the two are
 * equal, since the continuous extension of Dormand-Prince costs no extra
 * evaluation, but a higher order interpolant does need further stages and the
 * distinction is what lets one be added without changing this class.
 * <p>
 * The arrays handed out by the accessors are copies, so a caller may keep and
 * modify them; the constants are therefore immutable in fact and not only by
 * convention.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">Wikipedia
 * Runge-Kutta methods</a>.
 *
 * @since 1.5.3
 */
public final class ButcherTableau {

    /**
     * The row sums of {@code a} may deviate from {@code c} by this much before
     * the constructor refuses the tableau. The coefficients are rationals
     * evaluated in binary, so exact agreement cannot be required.
     */
    private static final double CONSISTENCY_TOL = 1.0e-12;

    /**
     * The classical fourth order Runge-Kutta method: four stages, no error
     * estimate and no continuous extension.
     * <p>
     * It is here for the fixed step case, where its cost per step is the lowest
     * of any fourth order method, and as the reference the adaptive methods are
     * measured against. Having no embedded pair, it cannot control its own step
     * size, and a caller that wants a value between two steps gets the cubic
     * Hermite interpolant {@link ExplicitRungeKutta} falls back on.
     */
    public static final ButcherTableau CLASSIC_RK4 = classicRk4();

    /**
     * Dormand-Prince 5(4), the standard explicit method for a non-stiff problem
     * at moderate accuracy.
     * <p>
     * Seven stages carrying a fifth order solution and an embedded fourth order
     * one, and the last stage of a step is the first stage of the next -- the
     * property called FSAL -- so an accepted step costs six evaluations rather
     * than seven. Its continuous extension is of order four and uses no stage
     * beyond the seven already computed.
     */
    public static final ButcherTableau DORMAND_PRINCE_45 = dormandPrince45();

    private final String name;
    private final int order;
    private final int embeddedOrder;
    private final double[] c;
    private final double[][] a;
    private final double[] b;
    private final double[] bStar;
    private final double[][] dense;
    private final boolean fsal;

    /**
     * A tableau without an embedded pair and without a continuous extension.
     *
     * @param name
     *            how the method is known, used in messages
     * @param order
     *            the order of the solution {@code b} advances
     * @param c
     *            the {@code s} stage times, {@code c[0]} being zero
     * @param a
     *            the strictly lower triangular stage coefficients,
     *            {@code a[i]} of length {@code i}
     * @param b
     *            the {@code s} weights of the advancing solution
     * @throws IllegalArgumentException
     *             if the shapes disagree or a row of {@code a} does not sum to
     *             its stage time
     */
    public ButcherTableau(String name, int order, double[] c, double[][] a, double[] b) {
        this(name, order, c, a, b, null, 0, null);
    }

    /**
     * A tableau in full.
     *
     * @param name
     *            how the method is known, used in messages
     * @param order
     *            the order of the solution {@code b} advances
     * @param c
     *            the stage times, {@code c[0]} being zero, one entry per stage
     *            including any that only the continuous extension needs
     * @param a
     *            the strictly lower triangular stage coefficients,
     *            {@code a[i]} of length {@code i}
     * @param b
     *            the weights of the advancing solution, one per propagating
     *            stage, at most as many as there are stage times
     * @param bStar
     *            the weights of the embedded solution of the lower order, of
     *            the same length as {@code b}, or {@code null} for a method
     *            that estimates no error
     * @param embeddedOrder
     *            the order of the embedded solution, ignored when
     *            {@code bStar} is {@code null}
     * @param dense
     *            per stage, the coefficients of the polynomial in {@code theta}
     *            that weighs that stage inside the step, lowest power first and
     *            the constant term omitted because it is always zero, or
     *            {@code null} for a method with no continuous extension
     * @throws IllegalArgumentException
     *             if the shapes disagree, a row of {@code a} does not sum to
     *             its stage time, or the continuous extension does not reproduce
     *             {@code b} at the end of the step
     */
    public ButcherTableau(String name, int order, double[] c, double[][] a, double[] b, double[] bStar,
            int embeddedOrder, double[][] dense) {
        if (name == null) {
            throw new IllegalArgumentException("name must not be null");
        }
        if (order < 1) {
            throw new IllegalArgumentException("order must be at least 1, not " + order);
        }
        if (c == null || a == null || b == null) {
            throw new IllegalArgumentException("c, a and b must not be null");
        }
        if (a.length != c.length) {
            throw new IllegalArgumentException("a has " + a.length + " rows but there are " + c.length
                    + " stage times");
        }
        if (b.length < 1 || b.length > c.length) {
            throw new IllegalArgumentException("b is of length " + b.length + ", which is not between 1 and "
                    + c.length);
        }
        for (int i = 0; i < a.length; ++i) {
            if (a[i] == null || a[i].length != i) {
                throw new IllegalArgumentException("row " + i + " of a must be of length " + i);
            }
            double rowSum = 0.0;
            for (int j = 0; j < i; ++j) {
                rowSum += a[i][j];
            }
            if (Math.abs(rowSum - c[i]) > CONSISTENCY_TOL) {
                throw new IllegalArgumentException("row " + i + " of a sums to " + rowSum + " but c[" + i
                        + "] is " + c[i]);
            }
        }
        if (bStar != null) {
            if (bStar.length != b.length) {
                throw new IllegalArgumentException("bStar is of length " + bStar.length + ", not " + b.length);
            }
            if (embeddedOrder < 1) {
                throw new IllegalArgumentException("embeddedOrder must be at least 1, not " + embeddedOrder);
            }
        }
        if (dense != null) {
            if (dense.length != c.length) {
                throw new IllegalArgumentException("dense has " + dense.length + " rows but there are "
                        + c.length + " stage times");
            }
            int degree = -1;
            for (int i = 0; i < dense.length; ++i) {
                if (dense[i] == null || dense[i].length < 1) {
                    throw new IllegalArgumentException("row " + i + " of dense must not be empty");
                }
                if (degree < 0) {
                    degree = dense[i].length;
                } else if (dense[i].length != degree) {
                    throw new IllegalArgumentException("row " + i + " of dense is of length " + dense[i].length
                            + ", not " + degree);
                }
            }
            for (int i = 0; i < dense.length; ++i) {
                double atEnd = 0.0;
                for (int j = 0; j < dense[i].length; ++j) {
                    atEnd += dense[i][j];
                }
                double weight = (i < b.length) ? b[i] : 0.0;
                if (Math.abs(atEnd - weight) > CONSISTENCY_TOL) {
                    throw new IllegalArgumentException("the continuous extension of stage " + i
                            + " reaches " + atEnd + " at the end of the step but b[" + i + "] is " + weight);
                }
            }
        }

        this.name = name;
        this.order = order;
        this.embeddedOrder = (bStar == null) ? 0 : embeddedOrder;
        this.c = c.clone();
        this.a = copy(a);
        this.b = b.clone();
        this.bStar = (bStar == null) ? null : bStar.clone();
        this.dense = (dense == null) ? null : copy(dense);
        this.fsal = isFirstSameAsLast(this.c, this.a, this.b);
    }

    /**
     * How the method is known.
     *
     * @return the name given to the constructor, never {@code null}
     */
    public String name() {
        return name;
    }

    /**
     * The order of the solution the method advances: halving the step size
     * divides the error accumulated over a fixed interval by two to this power.
     *
     * @return the order, at least one
     */
    public int order() {
        return order;
    }

    /**
     * The number of field evaluations that advance the solution by one step.
     *
     * @return the number of propagating stages, at least one
     */
    public int stages() {
        return b.length;
    }

    /**
     * The number of field evaluations a step costs when a value strictly inside
     * it is wanted, which is {@link #stages()} for a method whose continuous
     * extension is free.
     *
     * @return the total number of stages, at least {@link #stages()}
     */
    public int denseStages() {
        return c.length;
    }

    /**
     * Whether the method carries an embedded solution of a lower order and can
     * therefore estimate the error of its own step.
     *
     * @return {@code true} if a second weight vector is present
     */
    public boolean hasEmbedded() {
        return bStar != null;
    }

    /**
     * The order of the embedded solution, one below the order of the method for
     * the pairs here.
     *
     * @return the embedded order, or zero if there is no embedded pair
     */
    public int embeddedOrder() {
        return embeddedOrder;
    }

    /**
     * Whether the method carries a continuous extension, that is a polynomial
     * in the fraction of the step that is accurate between the endpoints rather
     * than only at them.
     *
     * @return {@code true} if dense output coefficients are present
     */
    public boolean hasDenseOutput() {
        return dense != null;
    }

    /**
     * The degree in {@code theta} of the continuous extension.
     *
     * @return the degree, or zero if there is no continuous extension
     */
    public int denseDegree() {
        return (dense == null) ? 0 : dense[0].length;
    }

    /**
     * Whether the last propagating stage of a step is the first stage of the
     * next one, so that an accepted step costs one evaluation less than it has
     * stages. This holds when the last stage time is one and the last row of
     * {@code a} is {@code b}.
     *
     * @return {@code true} if the method is first same as last
     */
    public boolean isFsal() {
        return fsal;
    }

    /**
     * The stage times, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #denseStages()}
     */
    public double[] c() {
        return c.clone();
    }

    /**
     * The stage coefficients, as a fresh copy, strictly lower triangular with
     * row {@code i} of length {@code i}.
     *
     * @return a <code>double[][]</code> of {@link #denseStages()} rows
     */
    public double[][] a() {
        return copy(a);
    }

    /**
     * The weights of the advancing solution, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #stages()}
     */
    public double[] b() {
        return b.clone();
    }

    /**
     * The weights of the embedded solution, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #stages()}, or
     *         {@code null} if there is no embedded pair
     */
    public double[] bStar() {
        return (bStar == null) ? null : bStar.clone();
    }

    /**
     * The continuous extension, as a fresh copy: row {@code i} holds the
     * coefficients of the polynomial in {@code theta} weighing stage {@code i},
     * lowest power first, starting at the linear term because the constant one
     * is zero for every stage.
     *
     * @return a <code>double[][]</code> of {@link #denseStages()} rows, or
     *         {@code null} if there is no continuous extension
     */
    public double[][] dense() {
        return (dense == null) ? null : copy(dense);
    }

    /**
     * The name of the method and its order.
     */
    @Override
    public String toString() {
        return name + " (order " + order + ", " + stages() + " stages"
                + (hasEmbedded() ? ", embedded order " + embeddedOrder : "") + (fsal ? ", FSAL" : "") + ")";
    }

    private static boolean isFirstSameAsLast(double[] c, double[][] a, double[] b) {
        int last = b.length - 1;
        if (last < 1 || Math.abs(c[last] - 1.0) > CONSISTENCY_TOL || a[last].length != last) {
            return false;
        }
        for (int j = 0; j < last; ++j) {
            if (Math.abs(a[last][j] - b[j]) > CONSISTENCY_TOL) {
                return false;
            }
        }
        return Math.abs(b[last]) <= CONSISTENCY_TOL;
    }

    private static double[][] copy(double[][] x) {
        double[][] out = new double[x.length][];
        for (int i = 0; i < x.length; ++i) {
            out[i] = x[i].clone();
        }
        return out;
    }

    private static ButcherTableau classicRk4() {
        double[] c = { 0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 };
        double[][] a = { {}, { 1.0 / 2.0 }, { 0.0, 1.0 / 2.0 }, { 0.0, 0.0, 1.0 } };
        double[] b = { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 };
        return new ButcherTableau("classical Runge-Kutta", 4, c, a, b);
    }

    private static ButcherTableau dormandPrince45() {
        double[] c = { 0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0 };
        double[][] a = {
                {},
                { 1.0 / 5.0 },
                { 3.0 / 40.0, 9.0 / 40.0 },
                { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0 },
                { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0 },
                { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0 },
                { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0 } };
        double[] b = { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0 };
        double[] bStar = { 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
                187.0 / 2100.0, 1.0 / 40.0 };
        double[] d = { -12715105075.0 / 11282082432.0, 0.0, 87487479700.0 / 32700410799.0,
                -10690763975.0 / 1880347072.0, 701980252875.0 / 199316789632.0, -1453857185.0 / 822651844.0,
                69997945.0 / 29380423.0 };
        return new ButcherTableau("Dormand-Prince 5(4)", 5, c, a, b, bStar, 4, quarticFsalDense(b, d));
    }

    /**
     * Turns the continuous extension of Dormand-Prince from the form it is
     * published in into one polynomial per stage.
     * <p>
     * The published form evaluates, with <code>u = 1 - theta</code>,
     * <p>
     * <code>y0 + theta (r2 + u (r3 + theta (r4 + u r5)))</code>
     * <p>
     * where <code>r2 = y1 - y0</code>, <code>r3 = h k[0] - r2</code>,
     * <code>r4 = r2 - h k[s-1] - r3</code> and
     * <code>r5 = h sum(d[i] k[i])</code>. Substituting
     * <code>r2 = h sum(b[i] k[i])</code> and collecting the powers of
     * {@code theta} leaves, for stage {@code i} with {@code e} one at the first
     * stage and {@code l} one at the last,
     * <p>
     * <code>e theta + (3 b - 2 e - l + d) theta^2 + (-2 b + e + l - 2 d)
     * theta^3 + d theta^4</code>
     * <p>
     * which is what this returns. Written this way the four conditions that
     * make it an interpolant -- value and derivative at each end of the step --
     * hold whatever the {@code d} are, so the coefficients that had to be
     * copied from a table are exactly the ones the fourth order accuracy
     * between the endpoints rests on, and nothing else.
     */
    private static double[][] quarticFsalDense(double[] b, double[] d) {
        int s = b.length;
        double[][] dense = new double[s][4];
        for (int i = 0; i < s; ++i) {
            double e = (i == 0) ? 1.0 : 0.0;
            double l = (i == s - 1) ? 1.0 : 0.0;
            dense[i][0] = e;
            dense[i][1] = 3.0 * b[i] - 2.0 * e - l + d[i];
            dense[i][2] = -2.0 * b[i] + e + l - 2.0 * d[i];
            dense[i][3] = d[i];
        }
        return dense;
    }

    /**
     * Whether the given object is a tableau with the same coefficients.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (!(other instanceof ButcherTableau)) {
            return false;
        }
        ButcherTableau that = (ButcherTableau) other;
        return order == that.order && embeddedOrder == that.embeddedOrder && Arrays.equals(c, that.c)
                && Arrays.deepEquals(a, that.a) && Arrays.equals(b, that.b)
                && Arrays.equals(bStar, that.bStar) && Arrays.deepEquals(dense, that.dense);
    }

    /**
     * A hash consistent with {@link #equals(Object)}.
     */
    @Override
    public int hashCode() {
        int h = 31 * order + embeddedOrder;
        h = 31 * h + Arrays.hashCode(c);
        h = 31 * h + Arrays.deepHashCode(a);
        h = 31 * h + Arrays.hashCode(b);
        h = 31 * h + Arrays.hashCode(bStar);
        return 31 * h + Arrays.deepHashCode(dense);
    }
}
