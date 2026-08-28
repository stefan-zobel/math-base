package math.ode;

import java.util.Arrays;

/**
 * The coefficients of an explicit symplectic method for a separable second
 * order system, as an immutable value.
 * <p>
 * A step of length {@code h} is a sequence of drifts and kicks:
 * <p>
 * <code>q += a[i] h v</code>, then <code>v += b[i] h f(q)</code>, for
 * {@code i} from zero
 * <p>
 * with <code>sum(a) = 1</code> and <code>sum(b) = 1</code>. A drift moves the
 * position at constant velocity and a kick changes the velocity at constant
 * position; each is the exact flow of one half of the Hamiltonian, and each is
 * therefore symplectic on its own. A composition of symplectic maps is
 * symplectic, which is the whole argument -- the property is structural and
 * does not depend on the coefficients being right, only the <em>order</em>
 * does.
 * <p>
 * <b>This is the same thing as a symplectic Runge-Kutta-Nystroem method.</b>
 * For a separable system the two constructions coincide, so a set of published
 * SRKN coefficients drops straight in here.
 * <p>
 * <b>A method begins either with a drift or with a kick, and it matters to what
 * a step costs.</b> One that begins with a drift ends with one -- its last kick
 * coefficient is zero -- and costs one evaluation per non-zero kick. One that
 * begins with a kick ends with a kick too, and the two that meet at a step
 * boundary sit at the same time and the same position, so the field is
 * evaluated there once and used twice: such a method costs one evaluation
 * <em>less</em> than it has kicks, and only its very first step pays the extra.
 * That is the saving Dormand-Prince gets from being first-same-as-last, reached
 * from the other side.
 * <p>
 * <b>What separates methods of the same order is the size of their
 * sub-steps.</b> Every method of order four or higher has at least one negative
 * coefficient, so parts of the step run backwards in time and the leading error
 * term grows with them. {@link #oneNorm()} measures that: it is one for a
 * method whose sub-steps all go forward, {@code 4.40} for {@link #YOSHIDA_4},
 * whose middle step is {@code -1.70} times the step it is embedded in, and only
 * {@code 1.16} for {@link #BLANES_MOAN_4}, which barely runs backwards at all.
 * The leading error constants that follow -- {@code 1.0e-01},
 * {@code 4.5e-03}, {@code 1.6e-05} for the three fourth order methods here --
 * are four orders of magnitude apart at the same order, which is why order
 * alone is not an argument for anything.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Symplectic_integrator">Wikipedia
 * symplectic integrator</a>.
 *
 * @since 1.5.3
 */
public final class SplittingCoefficients {

    /** How far the coefficient sums may stray from one. */
    private static final double CONSISTENCY_TOL = 1.0e-13;

    /**
     * The Stoermer-Verlet method, also called leapfrog: drift a half step, kick
     * a whole one, drift the other half.
     * <p>
     * Order two at one evaluation per step, and the base every composition here
     * is built from. Its sub-steps are both positive and its
     * {@link #oneNorm()} is one, which no method of higher order can manage.
     */
    public static final SplittingCoefficients VERLET = new SplittingCoefficients("Stoermer-Verlet", 2,
            new double[] { 0.5, 0.5 }, new double[] { 1.0, 0.0 });

    /**
     * Yoshida's triple jump: three Verlet steps of lengths
     * <code>w1 h, w0 h, w1 h</code>, with
     * <code>w1 = 1 / (2 - 2^(1/3))</code> and
     * <code>w0 = -2^(1/3) / (2 - 2^(1/3))</code>.
     * <p>
     * Order four at three evaluations per step, and the cheapest fourth order
     * method there is. It pays for that with the largest error constant in the
     * family: the middle step is {@code -1.70} times {@code h}, so the method
     * covers {@code 4.40} units of time to advance one, and the leading error
     * term grows with the sub-steps rather than with the step.
     */
    public static final SplittingCoefficients YOSHIDA_4 = yoshida4();

    /**
     * Suzuki's five-fold composition: five Verlet steps of lengths
     * <code>w h, w h, (1 - 4w) h, w h, w h</code> with
     * <code>w = 1 / (4 - 4^(1/3))</code>.
     * <p>
     * Order four at five evaluations per step. It buys two more evaluations
     * than {@link #YOSHIDA_4} and spends them on smaller sub-steps -- the
     * negative one is {@code -0.66} rather than {@code -1.70} -- which is the
     * trade every optimized method of this kind makes. Whether it is worth it
     * is a question about accuracy per evaluation and not about order, and
     * {@code SymplecticNystromTest} measures it.
     */
    public static final SplittingCoefficients SUZUKI_4 = suzuki4();

    /**
     * Blanes and Moan's {@code SRKN 6b}: order four at six evaluations, in the
     * form that begins and ends with a kick.
     * <p>
     * Six evaluations against Yoshida's three, and a leading error constant of
     * {@code 1.6e-05} against Yoshida's {@code 1.0e-01} -- a factor of six
     * thousand at the same order. After paying for the extra evaluations it
     * reaches a given accuracy for about <b>a fifth</b> of the cost of the
     * triple jump and a third of Suzuki's, whatever the accuracy is, since the
     * ratio between two methods of the same order does not depend on it.
     * <p>
     * The coefficients are from Blanes and Moan, <i>Practical symplectic
     * partitioned Runge-Kutta and Runge-Kutta-Nystrom methods</i>, J. Comput.
     * Appl. Math. 142 (2002) 313-330; the two that close the sequence follow
     * from the sums rather than being copied, and the order is verified rather
     * than trusted.
     */
    public static final SplittingCoefficients BLANES_MOAN_4 = blanesMoan4();

    /**
     * Blanes and Moan's {@code SRKN 11b}: order six at eleven evaluations, in
     * the form that begins and ends with a kick.
     * <p>
     * Two orders higher than anything else here, at a leading error constant of
     * {@code 2.2e-07}. Being of a different order it cannot be compared to the
     * others by a single number: it is cheaper than {@link #BLANES_MOAN_4} at
     * every accuracy tried from {@code 1e-6} downwards, and the margin widens
     * as the accuracy demanded tightens, which is what a higher order buys.
     * <p>
     * Same source as {@link #BLANES_MOAN_4}.
     */
    public static final SplittingCoefficients BLANES_MOAN_6 = blanesMoan6();

    private final String name;
    private final int order;
    private final double[] a;
    private final double[] b;
    private final boolean velocityFirst;
    private final int evaluations;

    /**
     * Coefficients given directly.
     *
     * @param name
     *            how the method is known, used in messages
     * @param order
     *            the order of the method
     * @param a
     *            the drift coefficients, summing to one
     * @param b
     *            the kick coefficients, of the same length, summing to one; a
     *            zero entry is a kick that does not happen
     * @throws IllegalArgumentException
     *             if the shapes disagree, a value is not finite, or a sum is
     *             not one
     */
    public SplittingCoefficients(String name, int order, double[] a, double[] b) {
        if (name == null) {
            throw new IllegalArgumentException("name must not be null");
        }
        if (order < 1) {
            throw new IllegalArgumentException("order must be at least 1, not " + order);
        }
        if (a == null || b == null) {
            throw new IllegalArgumentException("a and b must not be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a is of length " + a.length + " and b of length " + b.length
                    + ", which must agree");
        }
        if (a.length < 1) {
            throw new IllegalArgumentException("there must be at least one stage");
        }
        double driftSum = 0.0;
        double kickSum = 0.0;
        int kicks = 0;
        for (int i = 0; i < a.length; ++i) {
            if (Double.isNaN(a[i]) || Double.isInfinite(a[i]) || Double.isNaN(b[i])
                    || Double.isInfinite(b[i])) {
                throw new IllegalArgumentException("stage " + i + " is not finite");
            }
            driftSum += a[i];
            kickSum += b[i];
            if (b[i] != 0.0) {
                ++kicks;
            }
        }
        if (Math.abs(driftSum - 1.0) > CONSISTENCY_TOL) {
            throw new IllegalArgumentException("the drift coefficients sum to " + driftSum + ", not one");
        }
        if (Math.abs(kickSum - 1.0) > CONSISTENCY_TOL) {
            throw new IllegalArgumentException("the kick coefficients sum to " + kickSum + ", not one");
        }
        this.name = name;
        this.order = order;
        this.a = a.clone();
        this.b = b.clone();
        this.velocityFirst = (a[0] == 0.0) && (b[b.length - 1] != 0.0);
        this.evaluations = velocityFirst ? kicks - 1 : kicks;
    }

    /**
     * A composition of Verlet steps, which is how every method here beyond
     * Verlet itself is built.
     * <p>
     * Each fraction is one Verlet step of that share of {@code h}. Adjacent
     * drifts merge -- the trailing half step of one and the leading half step
     * of the next are one drift -- so a composition of {@code s} Verlet steps
     * costs {@code s} evaluations rather than {@code 2s}, and the fractions
     * become the kick coefficients unchanged.
     *
     * @param name
     *            how the composition is known
     * @param order
     *            the order it claims, which the tests verify rather than trust
     * @param fractions
     *            the share of the step each Verlet step takes, summing to one;
     *            beyond order two at least one of them is negative
     * @return the drift and kick coefficients of the composition
     * @throws IllegalArgumentException
     *             if the fractions are missing, empty, not finite or do not sum
     *             to one
     */
    public static SplittingCoefficients composition(String name, int order, double[] fractions) {
        if (fractions == null || fractions.length < 1) {
            throw new IllegalArgumentException("there must be at least one fraction");
        }
        int s = fractions.length;
        double[] a = new double[s + 1];
        double[] b = new double[s + 1];
        for (int i = 0; i < s; ++i) {
            b[i] = fractions[i];
        }
        b[s] = 0.0;
        a[0] = 0.5 * fractions[0];
        for (int i = 1; i < s; ++i) {
            a[i] = 0.5 * (fractions[i - 1] + fractions[i]);
        }
        a[s] = 0.5 * fractions[s - 1];
        return new SplittingCoefficients(name, order, a, b);
    }

    /**
     * A palindromic method that begins and ends with a kick, from the half of
     * its coefficients that gets published.
     * <p>
     * The whole sequence is
     * <code>B A B A ... B A B</code> with one more kick than drift, and it is
     * symmetric, so a table gives only the first half of each. Which half has
     * the repeated middle element follows from the two lengths and does not
     * have to be said: as many kicks as drifts means the middle drift stands
     * alone, one more kick than drift means the middle kick does.
     *
     * @param name
     *            how the method is known
     * @param order
     *            the order it claims, which the tests verify rather than trust
     * @param halfDrifts
     *            the drift coefficients up to and including the middle
     * @param halfKicks
     *            the kick coefficients up to and including the middle, as many
     *            as there are drifts or one more
     * @return the full drift and kick sequences, with a leading zero drift
     * @throws IllegalArgumentException
     *             if either half is missing or empty, the two lengths are not
     *             equal or one apart, or the mirrored coefficients do not sum
     *             to one
     */
    public static SplittingCoefficients velocityFirst(String name, int order, double[] halfDrifts,
            double[] halfKicks) {
        if (halfDrifts == null || halfKicks == null || halfDrifts.length < 1 || halfKicks.length < 1) {
            throw new IllegalArgumentException("there must be at least one drift and one kick");
        }
        double[] drifts;
        double[] kicks;
        if (halfKicks.length == halfDrifts.length + 1) {
            drifts = mirror(halfDrifts, false);
            kicks = mirror(halfKicks, true);
        } else if (halfKicks.length == halfDrifts.length) {
            drifts = mirror(halfDrifts, true);
            kicks = mirror(halfKicks, false);
        } else {
            throw new IllegalArgumentException("there must be as many kicks as drifts or one more, not "
                    + halfKicks.length + " against " + halfDrifts.length);
        }
        double[] a = new double[kicks.length];
        System.arraycopy(drifts, 0, a, 1, drifts.length);
        return new SplittingCoefficients(name, order, a, kicks);
    }

    private static double[] mirror(double[] half, boolean singleMiddle) {
        int n = singleMiddle ? (2 * half.length - 1) : (2 * half.length);
        double[] out = new double[n];
        for (int i = 0; i < half.length; ++i) {
            out[i] = half[i];
            out[n - 1 - i] = half[i];
        }
        return out;
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
     * The order of the method: halving the step size divides the error
     * accumulated over a fixed interval by two to this power.
     *
     * @return the order, at least one
     */
    public int order() {
        return order;
    }

    /**
     * The number of drift and kick pairs in one step.
     *
     * @return the number of stages, at least one
     */
    public int stages() {
        return a.length;
    }

    /**
     * The number of field evaluations one step costs in a run: the number of
     * kicks that are not zero, less one for a method that begins with a kick,
     * because that kick falls at the same time and the same position as the
     * last kick of the step before and is not computed twice. The very first
     * step of a run pays the one extra.
     *
     * @return the cost of a step, at least one
     */
    public int evaluations() {
        return evaluations;
    }

    /**
     * Whether the method begins with a kick rather than with a drift, which is
     * what lets one evaluation be shared between consecutive steps.
     *
     * @return {@code true} if the first drift is zero and the last kick is not
     */
    public boolean isVelocityFirst() {
        return velocityFirst;
    }

    /**
     * The sum of the absolute kick coefficients: the amount of time the method
     * travels through in order to advance by one step.
     * <p>
     * One for a method whose sub-steps all go forward, which above order two is
     * impossible, and larger the further a method has to run backwards. It is
     * not the error constant, but it is what the error constant is driven by,
     * and it separates two methods of the same order at a glance.
     *
     * @return the one-norm of the kick coefficients, at least one
     */
    public double oneNorm() {
        double sum = 0.0;
        for (int i = 0; i < b.length; ++i) {
            sum += Math.abs(b[i]);
        }
        return sum;
    }

    /**
     * Whether the drifts and the kicks read the same forwards and backwards,
     * which is what makes the method time reversible and forces its order to be
     * even.
     *
     * @return {@code true} if the sequence is palindromic
     */
    public boolean isPalindromic() {
        int s = a.length;
        if (velocityFirst) {
            // B A B ... A B, so the kicks mirror about their own middle and the
            // drifts about theirs, one place along
            for (int i = 0; i < s; ++i) {
                if (b[i] != b[s - 1 - i]) {
                    return false;
                }
            }
            for (int i = 1; i < s; ++i) {
                if (a[i] != a[s - i]) {
                    return false;
                }
            }
            return true;
        }
        if (b[s - 1] != 0.0) {
            return false;
        }
        for (int i = 0; i < s; ++i) {
            if (a[i] != a[s - 1 - i]) {
                return false;
            }
        }
        for (int i = 0; i <= s - 2; ++i) {
            if (b[i] != b[s - 2 - i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * The drift coefficients, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #stages()}
     */
    public double[] a() {
        return a.clone();
    }

    /**
     * The kick coefficients, as a fresh copy.
     *
     * @return a <code>double[]</code> of length {@link #stages()}
     */
    public double[] b() {
        return b.clone();
    }

    /**
     * The name of the method, its order and what a step costs.
     */
    @Override
    public String toString() {
        return name + " (order " + order + ", " + evaluations + " evaluations per step, one-norm "
                + oneNorm() + ")";
    }

    /**
     * Whether the given object holds the same coefficients.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (!(other instanceof SplittingCoefficients)) {
            return false;
        }
        SplittingCoefficients that = (SplittingCoefficients) other;
        return order == that.order && Arrays.equals(a, that.a) && Arrays.equals(b, that.b);
    }

    /**
     * A hash consistent with {@link #equals(Object)}.
     */
    @Override
    public int hashCode() {
        return 31 * (31 * order + Arrays.hashCode(a)) + Arrays.hashCode(b);
    }

    /**
     * The triple jump, whose two conditions are that the fractions sum to one
     * and that their cubes cancel -- the second being what kills the third
     * order term of a symmetric method and so lifts it from two to four.
     */
    private static SplittingCoefficients yoshida4() {
        double cubeRootOfTwo = Math.cbrt(2.0);
        double outer = 1.0 / (2.0 - cubeRootOfTwo);
        double middle = -cubeRootOfTwo / (2.0 - cubeRootOfTwo);
        return composition("Yoshida 4", 4, new double[] { outer, middle, outer });
    }

    /**
     * The five-fold composition, which solves the same two conditions with four
     * equal outer steps instead of two, buying smaller sub-steps for two more
     * evaluations.
     */
    private static SplittingCoefficients suzuki4() {
        double outer = 1.0 / (4.0 - Math.cbrt(4.0));
        double middle = 1.0 - 4.0 * outer;
        return composition("Suzuki 4", 4, new double[] { outer, outer, middle, outer, outer });
    }

    /**
     * {@code SRKN 6b} of Blanes and Moan, transcribed. The last coefficient of
     * each half is not in the table: it is what the sum has to be, and deriving
     * it rather than copying it is one fewer digit to get wrong.
     */
    private static SplittingCoefficients blanesMoan4() {
        double b1 = 0.0829844064174052;
        double b2 = 0.396309801498368;
        double b3 = -0.0390563049223486;
        double b4 = 1.0 - 2.0 * (b1 + b2 + b3);
        double a1 = 0.245298957184271;
        double a2 = 0.604872665711080;
        double a3 = 0.5 - (a1 + a2);
        return velocityFirst("Blanes-Moan 4", 4, new double[] { a1, a2, a3 },
                new double[] { b1, b2, b3, b4 });
    }

    /**
     * {@code SRKN 11b} of Blanes and Moan, transcribed. Here it is the middle
     * <em>kick</em> that stands alone rather than the middle drift, which is
     * why the two halves are of the same length.
     */
    private static SplittingCoefficients blanesMoan6() {
        double b1 = 0.0414649985182624;
        double b2 = 0.198128671918067;
        double b3 = -0.0400061921041533;
        double b4 = 0.0752539843015807;
        double b5 = -0.0115113874206879;
        double b6 = 0.5 - (b1 + b2 + b3 + b4 + b5);
        double a1 = 0.123229775946271;
        double a2 = 0.290553797799558;
        double a3 = -0.127049212625417;
        double a4 = -0.246331761062075;
        double a5 = 0.357208872795928;
        double a6 = 1.0 - 2.0 * (a1 + a2 + a3 + a4 + a5);
        return velocityFirst("Blanes-Moan 6", 6, new double[] { a1, a2, a3, a4, a5, a6 },
                new double[] { b1, b2, b3, b4, b5, b6 });
    }
}
