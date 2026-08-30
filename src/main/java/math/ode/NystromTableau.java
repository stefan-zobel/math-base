package math.ode;

import java.util.Arrays;

/**
 * The coefficients of an explicit Runge-Kutta-Nystroem method, as an immutable
 * value.
 * <p>
 * A step of length {@code h} from <code>(t, q, v)</code> evaluates the
 * acceleration at {@code s} stages,
 * <p>
 * <code>Q[i] = q + c[i] h v + h^2 sum(a[i][j] g[j], j &lt; i)</code>,
 * <code>g[i] = f(t + c[i] h, Q[i])</code>
 * <p>
 * and advances the state by
 * <code>q + h v + h^2 sum(bbar[i] g[i])</code> and
 * <code>v + h sum(b[i] g[i])</code>. Two weight vectors rather than one is what
 * separates this from a {@link ButcherTableau}, and the stage matrix multiplies
 * <code>h^2</code> rather than {@code h}.
 * <p>
 * <b>This is for the second order form and not for a mechanical problem.</b>
 * The equation is <code>q'' = f(t, q)</code>, which is as much a
 * semi-discretized wave equation or a shooting reduction of a boundary value
 * problem as it is a force law. What the form buys is stages: because the
 * acceleration does not read the velocity, order conditions drop out and a
 * given order needs fewer evaluations than the same order on the flattened
 * system.
 * <p>
 * <b>This is data, not an algorithm</b>, as {@link ButcherTableau} is;
 * {@link NystromRungeKutta} executes whichever tableau it is handed. The
 * constructor checks the conditions every such method must satisfy whatever its
 * order -- that the weights sum to <code>1/2</code> and to {@code 1}, and that
 * each row of the stage matrix sums to <code>c[i]^2 / 2</code> -- because a
 * mistyped coefficient does not throw, it silently lowers the order, and the
 * embedded error estimate then understates the damage in the same proportion.
 * <p>
 * The arrays handed out by the accessors are copies, so the constants are
 * immutable in fact and not only by convention.
 * <p>
 * <b>See</b> <a href=
 * "https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Nystr%C3%B6m_methods">
 * Wikipedia Nystroem methods</a>.
 *
 * @since 1.5.3
 */
public final class NystromTableau {

    /**
     * How far a sum may sit from the value the order conditions demand. The
     * coefficients of {@link #RKN6_4} are algebraic rather than rational, so
     * they are formed in floating point and cannot be expected to hit the mark
     * exactly; the conditions still held to {@code 4e-16} when they were
     * measured, and a transcription error misses by a factor rather than by an
     * ulp -- the one this class caught while it was being written was out by
     * 37 %.
     */
    private static final double CONDITION_TOLERANCE = 1.0e-12;

    private static final double ROOT_8581 = Math.sqrt(8581.0);

    /**
     * RKN6(4)6FM of Dormand, El-Mikkawy and Prince: order six with an embedded
     * four over six stages, first-same-as-last, and a continuous extension of
     * degree four.
     * <p>
     * Its coefficients are algebraic in <code>sqrt(8581)</code> and are written
     * here as the expressions that generate them rather than as decimals. The
     * observed order was measured at {@code 6.00} for the advancing solution
     * and {@code 4.00} for the embedded one, by halving the step on problems
     * with an exact solution, and the tests keep that check rather than trusting
     * the transcription.
     */
    public static final NystromTableau RKN6_4 = rkn64();

    private final String name;
    private final int order;
    private final int embeddedOrder;
    private final double[] c;
    private final double[][] a;
    private final double[] bbar;
    private final double[] b;
    private final double[] bbarStar;
    private final double[] bStar;
    private final double[][] densePosition;
    private final double[][] denseVelocity;
    private final boolean fsal;

    /**
     * A tableau in full.
     *
     * @param name
     *            how the method is known, used in messages
     * @param order
     *            the order of the solution {@code bbar} and {@code b} advance
     * @param embeddedOrder
     *            the order of the embedded solution, ignored when
     *            {@code bbarStar} is {@code null}
     * @param c
     *            the stage times, {@code c[0]} being zero, one per stage
     * @param a
     *            the strictly lower triangular stage matrix, {@code a[i]} of
     *            length {@code i}, multiplying <code>h^2</code>
     * @param bbar
     *            the weights of the position update, one per stage
     * @param b
     *            the weights of the velocity update, one per stage
     * @param bbarStar
     *            the position weights of the embedded solution, or {@code null}
     *            for a method that estimates no error
     * @param bStar
     *            its velocity weights, of the same length, or {@code null}
     * @param densePosition
     *            per stage, the coefficients of the polynomial in
     *            <code>theta</code> that weighs that stage in the position
     *            interpolant, lowest power first, or {@code null} for a method
     *            with no continuous extension
     * @param denseVelocity
     *            the same for the velocity interpolant
     * @throws IllegalArgumentException
     *             if the shapes disagree, a coefficient is not finite, the
     *             weights do not sum to <code>1/2</code> and {@code 1}, a row of
     *             {@code a} does not sum to <code>c[i]^2 / 2</code>, or the
     *             continuous extension does not reproduce the weights at the end
     *             of the step
     */
    public NystromTableau(String name, int order, int embeddedOrder, double[] c, double[][] a,
            double[] bbar, double[] b, double[] bbarStar, double[] bStar,
            double[][] densePosition, double[][] denseVelocity) {
        if (name == null) {
            throw new IllegalArgumentException("name must not be null");
        }
        if (order < 1) {
            throw new IllegalArgumentException("order must be at least 1, not " + order);
        }
        if (c == null || a == null || bbar == null || b == null) {
            throw new IllegalArgumentException("c, a, bbar and b must not be null");
        }
        int s = c.length;
        if (s < 1) {
            throw new IllegalArgumentException("a method needs at least one stage");
        }
        if (a.length != s || bbar.length != s || b.length != s) {
            throw new IllegalArgumentException("a, bbar and b must have one entry per stage, "
                    + s + ", not " + a.length + ", " + bbar.length + " and " + b.length);
        }
        if (c[0] != 0.0) {
            throw new IllegalArgumentException("the first stage must sit at the start of the step,"
                    + " so c[0] must be zero, not " + c[0]);
        }
        finite("c", c);
        finite("bbar", bbar);
        finite("b", b);
        for (int i = 0; i < s; ++i) {
            if (a[i] == null || a[i].length != i) {
                throw new IllegalArgumentException("a[" + i + "] must be of length " + i + ", not "
                        + (a[i] == null ? "null" : Integer.toString(a[i].length)));
            }
            finite("a[" + i + "]", a[i]);
            double rowSum = 0.0;
            for (int j = 0; j < i; ++j) {
                rowSum += a[i][j];
            }
            near("row " + i + " of a", rowSum, 0.5 * c[i] * c[i]);
        }
        near("the sum of bbar", sum(bbar), 0.5);
        near("the sum of b", sum(b), 1.0);

        if ((bbarStar == null) != (bStar == null)) {
            throw new IllegalArgumentException(
                    "an embedded solution needs both of its weight vectors or neither");
        }
        if (bbarStar != null) {
            if (bbarStar.length != s || bStar.length != s) {
                throw new IllegalArgumentException("the embedded weights must have one entry per"
                        + " stage, " + s + ", not " + bbarStar.length + " and " + bStar.length);
            }
            finite("bbarStar", bbarStar);
            finite("bStar", bStar);
            near("the sum of bbarStar", sum(bbarStar), 0.5);
            near("the sum of bStar", sum(bStar), 1.0);
            if (embeddedOrder < 1) {
                throw new IllegalArgumentException(
                        "embeddedOrder must be at least 1, not " + embeddedOrder);
            }
        }

        if ((densePosition == null) != (denseVelocity == null)) {
            throw new IllegalArgumentException(
                    "a continuous extension needs both of its coefficient sets or neither");
        }
        if (densePosition != null) {
            checkDense("densePosition", densePosition, bbar, s);
            checkDense("denseVelocity", denseVelocity, b, s);
        }

        this.name = name;
        this.order = order;
        this.embeddedOrder = (bbarStar == null) ? 0 : embeddedOrder;
        this.c = c.clone();
        this.a = copy(a);
        this.bbar = bbar.clone();
        this.b = b.clone();
        this.bbarStar = (bbarStar == null) ? null : bbarStar.clone();
        this.bStar = (bStar == null) ? null : bStar.clone();
        this.densePosition = (densePosition == null) ? null : copy(densePosition);
        this.denseVelocity = (denseVelocity == null) ? null : copy(denseVelocity);
        this.fsal = isFirstSameAsLast(this.c, this.a, this.bbar);
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
     * The order of the solution the method advances.
     *
     * @return the order, which is what a step controller raises the error to
     */
    public int order() {
        return order;
    }

    /**
     * The order of the embedded solution, or zero where there is none.
     *
     * @return the embedded order
     */
    public int embeddedOrder() {
        return embeddedOrder;
    }

    /**
     * How many times a step evaluates the acceleration.
     *
     * @return the number of stages
     */
    public int stages() {
        return c.length;
    }

    /**
     * Whether the method carries a second solution to judge its steps by.
     *
     * @return {@code true} where an adaptive step size is possible
     */
    public boolean hasErrorEstimate() {
        return bbarStar != null;
    }

    /**
     * Whether the method can be asked for the state strictly inside a step at
     * its own accuracy rather than by a Hermite interpolant.
     *
     * @return {@code true} where a continuous extension is present
     */
    public boolean hasDenseOutput() {
        return densePosition != null;
    }

    /**
     * Whether the last stage of a step is the first stage of the next, which
     * saves one evaluation per step. Derived from the coefficients rather than
     * declared: it holds when the last stage sits at the end of the step and
     * its row equals the position weights.
     *
     * @return {@code true} for a first-same-as-last method
     */
    public boolean isFsal() {
        return fsal;
    }

    /**
     * The stage times.
     *
     * @return a copy of {@code c}
     */
    public double[] c() {
        return c.clone();
    }

    /**
     * The stage matrix, which multiplies the square of the step size.
     *
     * @return a copy of {@code a}
     */
    public double[][] a() {
        return copy(a);
    }

    /**
     * The weights of the position update.
     *
     * @return a copy of {@code bbar}
     */
    public double[] bbar() {
        return bbar.clone();
    }

    /**
     * The weights of the velocity update.
     *
     * @return a copy of {@code b}
     */
    public double[] b() {
        return b.clone();
    }

    /**
     * The position weights of the embedded solution.
     *
     * @return a copy, or {@code null} where there is no error estimate
     */
    public double[] bbarStar() {
        return (bbarStar == null) ? null : bbarStar.clone();
    }

    /**
     * The velocity weights of the embedded solution.
     *
     * @return a copy, or {@code null} where there is no error estimate
     */
    public double[] bStar() {
        return (bStar == null) ? null : bStar.clone();
    }

    /**
     * The position interpolant, per stage the coefficients of a polynomial in
     * <code>theta</code>, lowest power first. The interpolant it belongs to is
     * <code>q + theta h v + h^2 theta^2 sum(P[i](theta) g[i])</code>.
     *
     * @return a copy, or {@code null} where there is no continuous extension
     */
    public double[][] densePosition() {
        return (densePosition == null) ? null : copy(densePosition);
    }

    /**
     * The velocity interpolant, in the same shape, belonging to
     * <code>v + h theta sum(Q[i](theta) g[i])</code>.
     *
     * @return a copy, or {@code null} where there is no continuous extension
     */
    public double[][] denseVelocity() {
        return (denseVelocity == null) ? null : copy(denseVelocity);
    }

    /**
     * The name, the two orders and the stage count.
     */
    @Override
    public String toString() {
        return "NystromTableau[" + name + ", order " + order
                + (hasErrorEstimate() ? "(" + embeddedOrder + ")" : "") + ", " + stages()
                + " stages" + (fsal ? ", FSAL" : "") + (hasDenseOutput() ? ", dense" : "") + "]";
    }

    /**
     * Two tableaux are equal when every coefficient is.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (!(other instanceof NystromTableau)) {
            return false;
        }
        NystromTableau that = (NystromTableau) other;
        return order == that.order && embeddedOrder == that.embeddedOrder
                && Arrays.equals(c, that.c) && Arrays.deepEquals(a, that.a)
                && Arrays.equals(bbar, that.bbar) && Arrays.equals(b, that.b)
                && Arrays.equals(bbarStar, that.bbarStar) && Arrays.equals(bStar, that.bStar)
                && Arrays.deepEquals(densePosition, that.densePosition)
                && Arrays.deepEquals(denseVelocity, that.denseVelocity);
    }

    @Override
    public int hashCode() {
        int h = order;
        h = 31 * h + embeddedOrder;
        h = 31 * h + Arrays.hashCode(c);
        h = 31 * h + Arrays.deepHashCode(a);
        h = 31 * h + Arrays.hashCode(bbar);
        h = 31 * h + Arrays.hashCode(b);
        h = 31 * h + Arrays.hashCode(bbarStar);
        h = 31 * h + Arrays.hashCode(bStar);
        h = 31 * h + Arrays.deepHashCode(densePosition);
        h = 31 * h + Arrays.deepHashCode(denseVelocity);
        return h;
    }

    private static double sum(double[] x) {
        double s = 0.0;
        for (int i = 0; i < x.length; ++i) {
            s += x[i];
        }
        return s;
    }

    private static void finite(String what, double[] x) {
        for (int i = 0; i < x.length; ++i) {
            if (Double.isNaN(x[i]) || Double.isInfinite(x[i])) {
                throw new IllegalArgumentException(what + "[" + i + "] is not finite: " + x[i]);
            }
        }
    }

    private static void near(String what, double got, double want) {
        double scale = Math.max(Math.abs(want), 1.0);
        if (Math.abs(got - want) > CONDITION_TOLERANCE * scale) {
            throw new IllegalArgumentException(what + " must be " + want + ", not " + got
                    + ", which is off by " + Math.abs(got - want)
                    + ". A Nystroem method of any order satisfies this, so the coefficients are"
                    + " not the ones intended");
        }
    }

    private static void checkDense(String what, double[][] dense, double[] weights, int s) {
        if (dense.length != s) {
            throw new IllegalArgumentException(what + " must have one polynomial per stage, " + s
                    + ", not " + dense.length);
        }
        for (int i = 0; i < s; ++i) {
            if (dense[i] == null || dense[i].length < 1) {
                throw new IllegalArgumentException(what + "[" + i + "] must hold at least one"
                        + " coefficient");
            }
            finite(what + "[" + i + "]", dense[i]);
            double atOne = sum(dense[i]);
            if (Math.abs(atOne - weights[i]) > CONDITION_TOLERANCE) {
                throw new IllegalArgumentException(what + "[" + i + "] comes to " + atOne
                        + " at the end of the step where the weight it belongs to is " + weights[i]
                        + ", so the continuous extension does not meet the step it extends");
            }
        }
    }

    private static boolean isFirstSameAsLast(double[] c, double[][] a, double[] bbar) {
        int s = c.length;
        if (s < 2 || c[s - 1] != 1.0) {
            return false;
        }
        for (int j = 0; j < s - 1; ++j) {
            if (a[s - 1][j] != bbar[j]) {
                return false;
            }
        }
        return bbar[s - 1] == 0.0;
    }

    private static double[][] copy(double[][] x) {
        double[][] out = new double[x.length][];
        for (int i = 0; i < x.length; ++i) {
            out[i] = x[i].clone();
        }
        return out;
    }

    private static NystromTableau rkn64() {
        double r = ROOT_8581;
        double[] c = { 0.0, (209.0 - r) / 900.0, (209.0 - r) / 450.0, (209.0 + r) / 450.0,
                9.0 / 10.0, 1.0 };
        double[][] a = { {}, { (26131.0 - 209.0 * r) / 810000.0 },
                { (26131.0 - 209.0 * r) / 607500.0, (26131.0 - 209.0 * r) / 303750.0 },
                { (980403512254.0 + 7781688431.0 * r) / 11694469921875.0,
                        -(1262884486208.0 + 15385481287.0 * r) / 11694469921875.0,
                        (7166233891441.0 + 78694563299.0 * r) / 46777879687500.0 },
                { -9.0 * (329260.0 + 3181.0 * r) / 27040000.0,
                        27.0 * (35129.0 + 3331.0 * r) / 13520000.0,
                        -27.0 * (554358343.0 + 31040327.0 * r) / 464060480000.0,
                        153.0 * (8555257.0 - 67973.0 * r) / 2745920000.0 },
                { 329.0 / 4212.0, 0.0, (84119543.0 + 366727.0 * r) / 409622616.0,
                        (84119543.0 - 366727.0 * r) / 409622616.0, 200.0 / 17901.0 } };
        double[] bbar = { 329.0 / 4212.0, 0.0, (84119543.0 + 366727.0 * r) / 409622616.0,
                (84119543.0 - 366727.0 * r) / 409622616.0, 200.0 / 17901.0, 0.0 };
        double[] b = { 329.0 / 4212.0, 0.0, (389225579.0 + 96856.0 * r) / 1024056540.0,
                (389225579.0 - 96856.0 * r) / 1024056540.0, 2000.0 / 17901.0, 1.0 / 20.0 };
        // the embedded solution, as the advancing one less the published difference
        double[] bbarStar = { bbar[0] - (329.0 / 4212.0 - (2701.0 + 23.0 * r) / 4563.0),
                bbar[1] - (9829.0 + 131.0 * r) / 9126.0,
                bbar[2] - ((84119543.0 + 366727.0 * r) / 409622616.0
                        - 5.0 * (1798.0 + 17.0 * r) / 9126.0),
                bbar[3] - bbar[3], bbar[4] - bbar[4], bbar[5] };
        double[] bStar = { b[0] - (329.0 / 4212.0 - 115.0 / 2106.0), b[1],
                b[2] - ((389225579.0 + 96856.0 * r) / 1024056540.0
                        - (84119543.0 + 366727.0 * r) / 256014135.0),
                b[3] - ((389225579.0 - 96856.0 * r) / 1024056540.0
                        - (84119543.0 - 366727.0 * r) / 256014135.0),
                b[4] - (2000.0 / 17901.0 - 6950.0 / 17901.0),
                b[5] - (1.0 / 20.0 + 1.0 / 10.0) };
        double[][] densePosition = {
                { 2106.0 / 4212.0, -5244.0 / 4212.0, 6386.0 / 4212.0, -3819.0 / 4212.0,
                        900.0 / 4212.0 },
                { 0.0, 0.0, 0.0, 0.0, 0.0 },
                { 0.0, 18.0 * (81356461.0 + 25954829.0 * r) / 22529243880.0,
                        (22190560391.0 - 1109665151.0 * r) / 22529243880.0,
                        -6.0 * (4929647204.0 - 156109769.0 * r) / 22529243880.0,
                        1800.0 * (5860823.0 - 152228.0 * r) / 22529243880.0 },
                { 0.0, 18.0 * (81356461.0 - 25954829.0 * r) / 22529243880.0,
                        (22190560391.0 + 1109665151.0 * r) / 22529243880.0,
                        -6.0 * (4929647204.0 + 156109769.0 * r) / 22529243880.0,
                        1800.0 * (5860823.0 + 152228.0 * r) / 22529243880.0 },
                { 0.0, 200.0 * 195.0 / 17901.0, -200.0 * 620.0 / 17901.0,
                        200.0 * 651.0 / 17901.0, -200.0 * 225.0 / 17901.0 },
                { 0.0, -117.0 / 110.0, 757.0 / 220.0, -823.0 / 220.0, 15.0 / 11.0 } };
        double[][] denseVelocity = {
                { 1.0, -15732.0 / 4212.0, 25544.0 / 4212.0, -19095.0 / 4212.0, 5400.0 / 4212.0 },
                { 0.0, 0.0, 0.0, 0.0, 0.0 },
                { 0.0, 27.0 * (81356461.0 + 25954829.0 * r) / 11264621940.0,
                        2.0 * (22190560391.0 - 1109665151.0 * r) / 11264621940.0,
                        -15.0 * (4929647204.0 - 156109769.0 * r) / 11264621940.0,
                        5400.0 * (5860823.0 - 152228.0 * r) / 11264621940.0 },
                { 0.0, 27.0 * (81356461.0 - 25954829.0 * r) / 11264621940.0,
                        2.0 * (22190560391.0 + 1109665151.0 * r) / 11264621940.0,
                        -15.0 * (4929647204.0 + 156109769.0 * r) / 11264621940.0,
                        5400.0 * (5860823.0 + 152228.0 * r) / 11264621940.0 },
                { 0.0, 1000.0 * 117.0 / 17901.0, -1000.0 * 496.0 / 17901.0,
                        1000.0 * 651.0 / 17901.0, -1000.0 * 270.0 / 17901.0 },
                { 0.0, -702.0 / 220.0, 3028.0 / 220.0, -4115.0 / 220.0, 1800.0 / 220.0 } };
        return new NystromTableau("RKN6(4)6FM", 6, 4, c, a, bbar, b, bbarStar, bStar,
                densePosition, denseVelocity);
    }
}
