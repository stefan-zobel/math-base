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
 * The three tableaux the library ships are {@link #CLASSIC_RK4},
 * {@link #DORMAND_PRINCE_45} and {@link #DOP853}.
 * <p>
 * <b>Stage counts come in two.</b> {@link #stages()} is how many evaluations
 * advance the solution; {@link #denseStages()} is how many are needed when a
 * value strictly inside the step is wanted. They are equal for the first two,
 * whose continuous extensions cost no extra evaluation, and differ by three for
 * DOP853, whose seventh order interpolant does not come free -- which is why
 * {@link ExplicitRungeKutta} evaluates those three only when an interior value
 * is actually asked for.
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
     * How stiff a step may look before it counts against the equation, for a
     * method that does not say. It is Hairer's value for a fifth order pair,
     * and since it is a property of the stability region a method with a longer
     * one gives its own.
     */
    private static final double DEFAULT_STIFFNESS_THRESHOLD = 3.25;

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

    /**
     * DOP853, Hairer and Wanner's eighth order method, for a non-stiff problem
     * at an accuracy Dormand-Prince reaches only by taking very many steps.
     * <p>
     * Twelve stages carry the eighth order solution, and the thirteenth is the
     * field at the new point, which is the first stage of the next step -- so
     * an accepted step costs twelve evaluations. Three further stages exist
     * only for the continuous extension, which is of order seven, and are
     * evaluated when a value strictly inside a step is asked for and not
     * otherwise.
     * <p>
     * <b>It carries two error estimates and needs both.</b> One embedded
     * solution is of order five and one of order three; either alone would
     * control the step on the accuracy of a solution that is discarded, and the
     * combination behaves like the local error of the eighth order solution
     * that is kept.
     */
    public static final ButcherTableau DOP853 = dop853();

    private final String name;
    private final int order;
    private final int embeddedOrder;
    private final int secondaryEmbeddedOrder;
    private final double[] c;
    private final double[][] a;
    private final double[] b;
    private final double[] bStar;
    private final double[] bStarSecondary;
    private final double[][] dense;
    private final double stiffnessThreshold;
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
        this(name, order, c, a, b, bStar, embeddedOrder, null, 0, dense, DEFAULT_STIFFNESS_THRESHOLD);
    }

    /**
     * A tableau in full, with a second error estimate and a stiffness threshold
     * of its own.
     * <p>
     * <b>Two estimates rather than one</b> is what an eighth order method needs.
     * A single embedded solution of order {@code p} makes the difference behave
     * like <code>h^(p+1)</code>, and controlling the step on that is
     * controlling the accuracy of the solution that is thrown away rather than
     * of the one that is kept. Combining a high and a low estimate --
     * {@link StepController#errorNorm(double[], double[], double[], double[])}
     * does the combining, since it owns the tolerances -- gives back a number
     * that falls at the order of the solution actually propagated.
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
     * @param bStarSecondary
     *            the weights of a second embedded solution, of a lower order
     *            again and of the same length as {@code b}, or {@code null} for
     *            a method with only one estimate
     * @param secondaryEmbeddedOrder
     *            the order of that second solution, ignored when
     *            {@code bStarSecondary} is {@code null}
     * @param dense
     *            per stage, the coefficients of the polynomial in {@code theta}
     *            that weighs that stage inside the step, lowest power first and
     *            the constant term omitted because it is always zero, or
     *            {@code null} for a method with no continuous extension
     * @param stiffnessThreshold
     *            how large {@link OdeStepper#stiffnessMeasure()} may grow before
     *            a step counts as stiff, which is a property of how far the
     *            stability region reaches along the negative real axis
     * @throws IllegalArgumentException
     *             if the shapes disagree, a row of {@code a} does not sum to
     *             its stage time, the continuous extension does not reproduce
     *             {@code b} at the end of the step, a second estimate is given
     *             without a first, or the threshold is not positive and finite
     */
    public ButcherTableau(String name, int order, double[] c, double[][] a, double[] b, double[] bStar,
            int embeddedOrder, double[] bStarSecondary, int secondaryEmbeddedOrder, double[][] dense,
            double stiffnessThreshold) {
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
        if (bStarSecondary != null) {
            if (bStar == null) {
                throw new IllegalArgumentException("a second embedded solution needs a first one");
            }
            if (bStarSecondary.length != b.length) {
                throw new IllegalArgumentException("bStarSecondary is of length " + bStarSecondary.length
                        + ", not " + b.length);
            }
            if (secondaryEmbeddedOrder < 1) {
                throw new IllegalArgumentException("secondaryEmbeddedOrder must be at least 1, not "
                        + secondaryEmbeddedOrder);
            }
        }
        if (!(stiffnessThreshold > 0.0) || Double.isInfinite(stiffnessThreshold)) {
            throw new IllegalArgumentException("stiffnessThreshold must be positive and finite, not "
                    + stiffnessThreshold);
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
        this.secondaryEmbeddedOrder = (bStarSecondary == null) ? 0 : secondaryEmbeddedOrder;
        this.c = c.clone();
        this.a = copy(a);
        this.b = b.clone();
        this.bStar = (bStar == null) ? null : bStar.clone();
        this.bStarSecondary = (bStarSecondary == null) ? null : bStarSecondary.clone();
        this.dense = (dense == null) ? null : copy(dense);
        this.stiffnessThreshold = stiffnessThreshold;
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
     * Whether the method carries a second embedded solution beside the first,
     * of a lower order again, so that the two together can be combined into an
     * error estimate that falls at the order of the solution the method
     * actually advances.
     *
     * @return {@code true} if a third weight vector is present
     */
    public boolean hasSecondaryEmbedded() {
        return bStarSecondary != null;
    }

    /**
     * The order of the second embedded solution.
     *
     * @return that order, or zero if there is no second embedded solution
     */
    public int secondaryEmbeddedOrder() {
        return secondaryEmbeddedOrder;
    }

    /**
     * The weights of the second embedded solution, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #stages()}, or
     *         {@code null} if there is no second embedded solution
     */
    public double[] bStarSecondary() {
        return (bStarSecondary == null) ? null : bStarSecondary.clone();
    }

    /**
     * How large a stiffness measure may grow over a step of this method before
     * the step counts as stiff. It is how far the stability region reaches
     * along the negative real axis, so a method of a higher order, whose region
     * is longer, tolerates more.
     *
     * @return the threshold, positive and finite
     */
    public double stiffnessThreshold() {
        return stiffnessThreshold;
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
                + (denseStages() > stages() ? " and " + (denseStages() - stages()) + " for the interpolant"
                        : "")
                + (hasEmbedded() ? ", embedded order " + embeddedOrder : "")
                + (hasSecondaryEmbedded() ? " and " + secondaryEmbeddedOrder : "")
                + (fsal ? ", FSAL" : "") + ")";
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
     * The coefficients of DOP853, in the layout this class expects.
     * <p>
     * Sixteen rows: the twelve stages that carry the eighth order solution, then
     * the field at the new point -- whose {@code a} row is {@code b} and whose
     * stage time is one, which is what makes the method first same as last --
     * and then the three the continuous extension needs. Only thirteen of them
     * carry weight in the solution, so {@link #stages()} is thirteen where
     * {@link #denseStages()} is sixteen.
     * <p>
     * <b>The two error rows are derived rather than copied.</b> What is
     * published is the <em>difference</em> between the eighth order weights and
     * the fifth order ones, so the embedded weights are that subtracted from
     * {@code b}; and the third order row is published as weights already. That
     * the first sums to zero and the second to one is the cheapest check there
     * is that both were read the right way round, and the order conditions are
     * the thorough one.
     */
    private static ButcherTableau dop853() {
        double[] c = {
                0.0, 0.526001519587677318785587544488e-01, 0.789002279381515978178381316732e-01,
                0.118350341907227396726757197510e+00, 0.281649658092772603273242802490e+00,
                0.333333333333333333333333333333e+00, 0.25e+00, 0.307692307692307692307692307692e+00,
                0.651282051282051282051282051282e+00, 0.6e+00, 0.857142857142857142857142857142e+00, 1.0,
                1.0, 0.1e+00, 0.2e+00, 0.777777777777777777777777777778e+00 };

        double[] b = { 5.42937341165687622380535766363e-2, 0.0, 0.0, 0.0, 0.0,
                4.45031289275240888144113950566e0, 1.89151789931450038304281599044e0,
                -5.8012039600105847814672114227e0, 3.1116436695781989440891606237e-1,
                -1.52160949662516078556178806805e-1, 2.01365400804030348374776537501e-1,
                4.47106157277725905176885569043e-2, 0.0 };

        double[] er = { 0.1312004499419488073250102996e-01, 0.0, 0.0, 0.0, 0.0,
                -0.1225156446376204440720569753e+01, -0.4957589496572501915214079952e+00,
                0.1664377182454986536961530415e+01, -0.3503288487499736816886487290e+00,
                0.3341791187130174790297318841e+00, 0.8192320648511571246570742613e-01,
                -0.2235530786388629525884427845e-01, 0.0 };

        double[] bhh = { 0.244094488188976377952755905512e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.733846688281611857341361741547e+00, 0.0, 0.0, 0.220588235294117647058823529412e-01, 0.0 };

        double[][] a = {
                {}, {5.26001519587677318785587544488e-2},
                {1.97250569845378994544595329183e-2, 5.91751709536136983633785987549e-2},
                {2.95875854768068491816892993775e-2, 0.0, 8.87627564304205475450678981324e-2},
                {2.41365134159266685502369798665e-1, 0.0, -8.84549479328286085344864962717e-1,
                9.24834003261792003115737966543e-1},
                {3.7037037037037037037037037037e-2, 0.0, 0.0, 1.70828608729473871279604482173e-1,
                1.25467687566822425016691814123e-1},
                {3.7109375e-2, 0.0, 0.0, 1.70252211019544039314978060272e-1,
                6.02165389804559606850219397283e-2, -1.7578125e-2},
                {3.70920001185047927108779319836e-2, 0.0, 0.0, 1.70383925712239993810214054705e-1,
                1.07262030446373284651809199168e-1, -1.53194377486244017527936158236e-2,
                8.27378916381402288758473766002e-3},
                {6.24110958716075717114429577812e-1, 0.0, 0.0, -3.36089262944694129406857109825e0,
                -8.68219346841726006818189891453e-1, 2.75920996994467083049415600797e1,
                2.01540675504778934086186788979e1, -4.34898841810699588477366255144e1},
                {4.77662536438264365890433908527e-1, 0.0, 0.0, -2.48811461997166764192642586468e0,
                -5.90290826836842996371446475743e-1, 2.12300514481811942347288949897e1,
                1.52792336328824235832596922938e1, -3.32882109689848629194453265587e1,
                -2.03312017085086261358222928593e-2},
                {-9.3714243008598732571704021658e-1, 0.0, 0.0, 5.18637242884406370830023853209e0,
                1.09143734899672957818500254654e0, -8.14978701074692612513997267357e0,
                -1.85200656599969598641566180701e1, 2.27394870993505042818970056734e1,
                2.49360555267965238987089396762e0, -3.0467644718982195003823669022e0},
                {2.27331014751653820792359768449e0, 0.0, 0.0, -1.05344954667372501984066689879e1,
                -2.00087205822486249909675718444e0, -1.79589318631187989172765950534e1,
                2.79488845294199600508499808837e1, -2.85899827713502369474065508674e0,
                -8.87285693353062954433549289258e0, 1.23605671757943030647266201528e1,
                6.43392746015763530355970484046e-1},
                Arrays.copyOf(b, 12),
                {5.61675022830479523392909219681e-2, 0.0, 0.0, 0.0, 0.0, 0.0,
                2.53500210216624811088794765333e-1, -2.46239037470802489917441475441e-1,
                -1.24191423263816360469010140626e-1, 1.5329179827876569731206322685e-1,
                8.20105229563468988491666602057e-3, 7.56789766054569976138603589584e-3, -8.298e-3},
                {3.18346481635021405060768473261e-2, 0.0, 0.0, 0.0, 0.0, 2.83009096723667755288322961402e-2,
                5.35419883074385676223797384372e-2, -5.49237485713909884646569340306e-2, 0.0, 0.0,
                -1.08347328697249322858509316994e-4, 3.82571090835658412954920192323e-4,
                -3.40465008687404560802977114492e-4, 1.41312443674632500278074618366e-1},
                {-4.28896301583791923408573538692e-1, 0.0, 0.0, 0.0, 0.0,
                -4.69762141536116384314449447206e0, 7.68342119606259904184240953878e0,
                4.06898981839711007970213554331e0, 3.56727187455281109270669543021e-1, 0.0, 0.0, 0.0,
                -1.39902416515901462129418009734e-3, 2.9475147891527723389556272149e0,
                -9.15095847217987001081870187138e0} };

        double[][] d = {
                {-0.84289382761090128651353491142e+01, 0.0, 0.0, 0.0, 0.0,
                0.56671495351937776962531783590e+00, -0.30689499459498916912797304727e+01,
                0.23846676565120698287728149680e+01, 0.21170345824450282767155149946e+01,
                -0.87139158377797299206789907490e+00, 0.22404374302607882758541771650e+01,
                0.63157877876946881815570249290e+00, -0.88990336451333310820698117400e-01,
                0.18148505520854727256656404962e+02, -0.91946323924783554000451984436e+01,
                -0.44360363875948939664310572000e+01},
                {0.10427508642579134603413151009e+02, 0.0, 0.0, 0.0, 0.0,
                0.24228349177525818288430175319e+03, 0.16520045171727028198505394887e+03,
                -0.37454675472269020279518312152e+03, -0.22113666853125306036270938578e+02,
                0.77334326684722638389603898808e+01, -0.30674084731089398182061213626e+02,
                -0.93321305264302278729567221706e+01, 0.15697238121770843886131091075e+02,
                -0.31139403219565177677282850411e+02, -0.93529243588444783865713862664e+01,
                0.35816841486394083752465898540e+02},
                {0.19985053242002433820987653617e+02, 0.0, 0.0, 0.0, 0.0,
                -0.38703730874935176555105901742e+03, -0.18917813819516756882830838328e+03,
                0.52780815920542364900561016686e+03, -0.11573902539959630126141871134e+02,
                0.68812326946963000169666922661e+01, -0.10006050966910838403183860980e+01,
                0.77771377980534432092869265740e+00, -0.27782057523535084065932004339e+01,
                -0.60196695231264120758267380846e+02, 0.84320405506677161018159903784e+02,
                0.11992291136182789328035130030e+02},
                {-0.25693933462703749003312586129e+02, 0.0, 0.0, 0.0, 0.0,
                -0.15418974869023643374053993627e+03, -0.23152937917604549567536039109e+03,
                0.35763911791061412378285349910e+03, 0.93405324183624310003907691704e+02,
                -0.37458323136451633156875139351e+02, 0.10409964950896230045147246184e+03,
                0.29840293426660503123344363579e+02, -0.43533456590011143754432175058e+02,
                0.96324553959188282948394950600e+02, -0.39177261675615439165231486172e+02,
                -0.14972683625798562581422125276e+03} };

        double[] bStar = new double[b.length];
        for (int i = 0; i < b.length; ++i) {
            bStar[i] = b[i] - er[i];
        }
        return new ButcherTableau("DOP853", 8, c, a, b, bStar, 5, bhh, 3, septicDense(b, d), 6.1);
    }

    /**
     * Turns the continuous extension of DOP853 from the form it is published in
     * into one polynomial per stage.
     * <p>
     * The published form evaluates, with <code>u = 1 - theta</code>,
     * <p>
     * <code>y0 + theta r2 + theta u r3 + theta^2 u r4 + theta^2 u^2 r5
     * + theta^3 u^2 r6 + theta^3 u^3 r7 + theta^4 u^3 r8</code>
     * <p>
     * where <code>r2 = y1 - y0</code>, <code>r3 = h k[0] - r2</code>,
     * <code>r4 = r2 - h k[s-1] - r3</code> and <code>r5</code> through
     * <code>r8</code> are the four published combinations of every stage.
     * Substituting <code>r2 = h sum(b[i] k[i])</code> and expanding the powers
     * of {@code u} leaves one polynomial of degree seven per stage, which is
     * what this returns.
     * <p>
     * As with {@link #quarticFsalDense(double[], double[])}, writing it this way
     * makes the four conditions that turn a polynomial into an interpolant --
     * the value and the derivative at each end of the step -- hold whatever the
     * four combinations are: at {@code theta = 0} only the first two terms
     * survive and leave {@code k[0]}, and at {@code theta = 1} every term
     * carrying {@code u^2} or more drops out and leaves the last stage. So the
     * numbers that had to be copied from a table are exactly the ones the
     * seventh order between the endpoints rests on, and nothing else.
     */
    private static double[][] septicDense(double[] b, double[][] d) {
        int s = d[0].length;
        double[][] dense = new double[s][];
        for (int i = 0; i < s; ++i) {
            double weight = (i < b.length) ? b[i] : 0.0;
            double first = (i == 0) ? 1.0 : 0.0;
            double last = (i == b.length - 1) ? 1.0 : 0.0;
            double[] poly = new double[8];
            addTerm(poly, weight, 1, 0);
            addTerm(poly, first - weight, 1, 1);
            addTerm(poly, 2.0 * weight - last - first, 2, 1);
            addTerm(poly, d[0][i], 2, 2);
            addTerm(poly, d[1][i], 3, 2);
            addTerm(poly, d[2][i], 3, 3);
            addTerm(poly, d[3][i], 4, 3);
            dense[i] = Arrays.copyOfRange(poly, 1, poly.length);
        }
        return dense;
    }

    /**
     * Adds <code>w theta^p (1 - theta)^q</code> to a polynomial held as its
     * coefficients, lowest power first.
     */
    private static void addTerm(double[] poly, double w, int p, int q) {
        double binomial = 1.0;
        for (int j = 0; j <= q; ++j) {
            poly[p + j] += (((j & 1) == 0) ? w : -w) * binomial;
            binomial = binomial * (q - j) / (j + 1);
        }
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
        return order == that.order && embeddedOrder == that.embeddedOrder
                && secondaryEmbeddedOrder == that.secondaryEmbeddedOrder
                && stiffnessThreshold == that.stiffnessThreshold && Arrays.equals(c, that.c)
                && Arrays.deepEquals(a, that.a) && Arrays.equals(b, that.b)
                && Arrays.equals(bStar, that.bStar) && Arrays.equals(bStarSecondary, that.bStarSecondary)
                && Arrays.deepEquals(dense, that.dense);
    }

    /**
     * A hash consistent with {@link #equals(Object)}.
     */
    @Override
    public int hashCode() {
        int h = 31 * order + embeddedOrder;
        h = 31 * h + secondaryEmbeddedOrder;
        h = 31 * h + Double.valueOf(stiffnessThreshold).hashCode();
        h = 31 * h + Arrays.hashCode(c);
        h = 31 * h + Arrays.deepHashCode(a);
        h = 31 * h + Arrays.hashCode(b);
        h = 31 * h + Arrays.hashCode(bStar);
        h = 31 * h + Arrays.hashCode(bStarSecondary);
        return 31 * h + Arrays.deepHashCode(dense);
    }
}
