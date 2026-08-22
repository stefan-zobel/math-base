package math.solve;

import java.util.Locale;

import math.fun.DFunction;

/**
 * Double-exponential quadrature, the tanh-sinh family, over a finite, a half
 * infinite or a doubly infinite interval. See
 * <a href="https://en.wikipedia.org/wiki/Tanh-sinh_quadrature">tanh-sinh
 * quadrature</a>.
 * <p>
 * The substitution turns the integral into a trapezoidal rule in {@code t},
 * {@code x = tanh(pi/2 * sinh t)} for {@code [a, b]},
 * {@code x = a + exp(pi/2 * sinh t)} for {@code [a, inf)} and
 * {@code x = sinh(pi/2 * sinh t)} for the whole line. The step is halved until
 * two levels agree.
 *
 * <h2>Why the nodes crowd towards the ends</h2>
 * The transformed integrand decays doubly exponentially, so the truncated
 * trapezoidal rule converges faster than any fixed order, and the nodes
 * approach an endpoint fast enough that it is never evaluated. That is what
 * makes an integrable singularity there a non-event rather than the hard case.
 * Against the adaptive Gauss-Kronrod rule of {@link AdaptiveGaussKronrod}, on
 * the same interval at the same tolerance:
 * <ul>
 *   <li>{@code x^-0.9} on {@code [0, 1]}: 17 digits for 74 evaluations against
 *       1.0 digit for 5,880.</li>
 *   <li>{@code 1/sqrt(x)} on {@code [0, 1]}: 17 digits for 74 against 5.0 for
 *       1,620.</li>
 *   <li>{@code exp(-x)/sqrt(x)} on {@code [0, inf)}, a half line and a pole at
 *       once: 15.6 digits for 437 evaluations.</li>
 * </ul>
 * It is not a replacement for that rule. On a polynomial or an oscillation
 * Gauss-Kronrod is the more accurate of the two - it is exact to degree 22 and
 * subdivides where the feature is - and it adapts to the interior, which a
 * fixed node set cannot. This rule is the specialist for what happens at the
 * ends.
 *
 * <h2>Halving the step reuses every evaluation</h2>
 * The nodes of one level are {@code t = j*h}, so the previous level is the even
 * {@code j} and only the odd ones are new. The weights depend on {@code t}
 * alone and every new node has a new {@code t}, so nothing is reused except the
 * sum itself - which is carried forward, making the whole ladder cost as much
 * as its finest level and need no value cache.
 *
 * <h2>When it declines to answer</h2>
 * A level is accepted when it differs from its predecessor by no more than
 * {@code epsTol} times the integral of {@code |f|} over the same nodes, which
 * is the scale QUADPACK uses. Measuring against the integral of {@code |f|}
 * rather than against the value, or against a fixed floor, is what makes the
 * rule refuse a particular kind of wrong answer: when only a single node lands
 * on the mass, halving the step halves the value, and a succession of levels
 * reading {@code 9.8e-154}, {@code 4.9e-154}, {@code 2.4e-154} is the signature
 * of having missed the integrand rather than of convergence. A fixed floor
 * accepts that at once. Two further consequences:
 * <ul>
 *   <li>An integrand that is zero at every node yields
 *       {@code converged == false}, because the rule cannot tell a function
 *       that is zero everywhere from one that is zero everywhere it looked. The
 *       value returned is {@code 0.0} and the caller decides.</li>
 *   <li>An integral that is genuinely zero, or genuinely {@code 1e-200},
 *       converges normally, since the scale comes from {@code |f|} and not from
 *       the result.</li>
 * </ul>
 *
 * <h2>A pole at a non-zero endpoint costs half the digits</h2>
 * The abscissa is built as a distance from the nearer endpoint, so that the
 * node next to a pole keeps its accuracy. That distance survives down to
 * {@code 1e-307} at an endpoint of zero, but at any other endpoint {@code b}
 * the difference {@code b - delta} collapses onto {@code b} as soon as
 * {@code delta} falls under one ulp, and the rule stops there. Measured on
 * {@code 1/sqrt(x(1-x))}: over {@code [0, 0.5]}, where the only pole sits at
 * zero, 15.8 digits for 148 evaluations; over {@code [0.5, 1]}, the same
 * function with its pole at one, 8.0 digits and no convergence reported. This
 * is a property of binary
 * floating point at the endpoint rather than of the rule, and no integrand can
 * work around it as long as it is handed a position instead of a distance.
 *
 * @since 1.5.2
 */
public final class DoubleExponential {

    /** tanh-sinh, a finite interval. */
    private static final int FINITE = 0;
    /** exp-sinh, {@code [a, inf)}. */
    private static final int UPPER = 1;
    /** sinh-sinh, the whole line. */
    private static final int WHOLE = 2;

    private static final double HALF_PI = 0.5 * Math.PI;

    /** Step of level 0, giving between 10 and 13 nodes. */
    private static final double H0 = 1.0;

    /** No level below this may be accepted, so that a claim rests on real work. */
    private static final int MIN_LEVEL = 2;

    /** Default finest level: about 38,000 evaluations, and only if it gets there. */
    private static final int MAX_LEVEL = 12;

    /** Ceiling on the explicit overload; level 16 is already 2.4 million evaluations. */
    private static final int LEVEL_CEILING = 16;

    /** The result of one integration, with the estimate that decided when to stop. */
    public static final class IntegralResult {

        /** The approximated integral, always taken from the finest level computed. */
        public final double value;
        /** The change from the second finest to the finest level. */
        public final double approximatedErrorEstimate;
        /** Integrand evaluations spent, over the whole ladder. */
        public final long evaluations;
        /** The finest level computed; level {@code k} has a step of {@code 2^-k}. */
        public final int level;
        /** Whether {@link #approximatedErrorEstimate} met the tolerance. */
        public final boolean converged;

        IntegralResult(double value, double approximatedErrorEstimate, long evaluations, int level,
                boolean converged) {
            this.value = value;
            this.approximatedErrorEstimate = approximatedErrorEstimate;
            this.evaluations = evaluations;
            this.level = level;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT,
                    "Value: %.8f (approx. Error: %.2e, %d evaluations, level %d, converged: %b)", value,
                    approximatedErrorEstimate, evaluations, level, converged);
        }
    }

    /** The two running sums of the ladder and what they have cost. */
    private static final class Accumulator {
        double sum;
        double absSum;
        long calls;
    }

    private DoubleExponential() {
        throw new AssertionError();
    }

    // =========================================================================
    // ENTRY POINTS
    // =========================================================================

    /**
     * Integrates {@code f} over {@code [a, b]}, where either limit may be
     * infinite, halving the step up to at most level 12.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
     * @param b
     *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
     * @param epsTol
     *            error tolerance, relative to the integral of {@code |f|}
     * @return the approximated integral, its error estimate, what it cost and
     *         whether the tolerance was met
     * @throws IllegalArgumentException
     *             if a limit is {@code NaN}, if {@code a} is
     *             {@code Double.POSITIVE_INFINITY} or if {@code b} is
     *             {@code Double.NEGATIVE_INFINITY}
     */
    public static IntegralResult integrate1D(DFunction f, double a, double b, double epsTol) {
        return integrate1D(f, a, b, epsTol, MIN_LEVEL, MAX_LEVEL);
    }

    /**
     * Integrates {@code f} over {@code [a, b]} with an explicit level range, for
     * an integrand known to be harder or easier than the default assumes.
     * {@code minLevel} is clamped to at least 2, {@code maxLevel} to at most 16,
     * and {@code maxLevel} is raised to {@code minLevel} if it is smaller.
     * <p>
     * A generous ceiling is close to free, because the ladder stops as soon as
     * two levels agree and each level costs as much as everything before it put
     * together. Level 12 is about 38,000 evaluations for a finite interval and
     * 56,000 for an infinite one, and is reached only by an integrand that never
     * settles.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
     * @param b
     *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
     * @param epsTol
     *            error tolerance, relative to the integral of {@code |f|}
     * @param minLevel
     *            the coarsest level that may be accepted
     * @param maxLevel
     *            the finest level to compute
     * @return the approximated integral, its error estimate, what it cost and
     *         whether the tolerance was met
     * @throws IllegalArgumentException
     *             if a limit is {@code NaN}, if {@code a} is
     *             {@code Double.POSITIVE_INFINITY} or if {@code b} is
     *             {@code Double.NEGATIVE_INFINITY}
     */
    public static IntegralResult integrate1D(DFunction f, double a, double b, double epsTol, int minLevel,
            int maxLevel) {
        if (Double.isNaN(a) || Double.isNaN(b)) {
            throw new IllegalArgumentException("limits must not be NaN: [" + a + ", " + b + "]");
        }
        if (a == Double.POSITIVE_INFINITY || b == Double.NEGATIVE_INFINITY) {
            throw new IllegalArgumentException("limits are the wrong way round: [" + a + ", " + b + "]");
        }
        int lo = Math.max(MIN_LEVEL, minLevel);
        int hi = Math.min(LEVEL_CEILING, Math.max(lo, maxLevel));

        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        if (aInf && bInf) {
            return ladder(f, WHOLE, 0.0, 0.0, epsTol, lo, hi);
        }
        if (bInf) {
            return ladder(f, UPPER, a, 0.0, epsTol, lo, hi);
        }
        if (aInf) {
            // (-inf, b] is [-b, inf) of the reflected integrand
            return ladder(x -> f.apply(-x), UPPER, -b, 0.0, epsTol, lo, hi);
        }
        if (a == b) {
            return new IntegralResult(0.0, 0.0, 0L, lo, true);
        }
        return ladder(f, FINITE, a, b, epsTol, lo, hi);
    }

    // =========================================================================
    // THE LADDER
    // =========================================================================

    private static IntegralResult ladder(DFunction f, int kind, double a, double b, double epsTol, int minLevel,
            int maxLevel) {
        Accumulator acc = new Accumulator();
        double scale = (kind == FINITE) ? 0.5 * (b - a) : 1.0;

        sweep(f, kind, a, b, H0, false, acc);
        double previous = acc.sum * H0 * scale;
        double err = Double.POSITIVE_INFINITY;

        for (int level = 1; level <= maxLevel; level++) {
            double h = H0 * Math.pow(2.0, -level);
            sweep(f, kind, a, b, h, true, acc);
            double value = acc.sum * h * scale;
            double absValue = acc.absSum * h * Math.abs(scale);
            err = Math.abs(value - previous);
            previous = value;
            // nothing may be certified before the rule has seen a non-zero
            // integrand, and nothing that agrees only on the scale of the
            // answer rather than on the scale of the integrand
            if (level >= minLevel && absValue > 0.0 && err <= epsTol * absValue) {
                return new IntegralResult(value, err, acc.calls, level, true);
            }
        }
        return new IntegralResult(previous, err, acc.calls, maxLevel, false);
    }

    /**
     * Adds the nodes {@code t = j*h} to the running sums: all {@code j} when
     * {@code odd} is false, which is a fresh level, only the odd {@code j} when
     * it is true, which is the refinement of the level before.
     */
    private static void sweep(DFunction f, int kind, double a, double b, double h, boolean odd, Accumulator acc) {
        int step = odd ? 2 : 1;
        for (int j = odd ? 1 : 0; ; j += step) {
            if (!node(f, kind, a, b, j * h, j == 0, acc)) {
                break;
            }
        }
        if (kind == UPPER) {
            // exp-sinh is not symmetric in t, so the negative half is its own loop
            for (int j = -1; ; j -= step) {
                if (!node(f, kind, a, b, j * h, false, acc)) {
                    break;
                }
            }
        }
    }

    /**
     * Adds the node at {@code t}, and reports whether the loop may go on. It may
     * not once the rule has run out of representable abscissae, which turns on
     * {@code t} and the limits alone and never on the integrand - one that
     * happens to vanish out there must not be allowed to cut the sweep short.
     */
    private static boolean node(DFunction f, int kind, double a, double b, double t, boolean center,
            Accumulator acc) {
        double u = HALF_PI * Math.sinh(t);
        double coshT = Math.cosh(t);

        if (kind == FINITE) {
            double ch = Math.cosh(u);
            if (Double.isInfinite(ch)) {
                return false;
            }
            double w = HALF_PI * coshT / (ch * ch);
            if (!(w > 0.0)) {
                return false;
            }
            // 1 - tanh(u) as a distance from the endpoint. Never as a
            // subtraction: mid +- half*tanh(u) cancels away exactly the
            // accuracy near the endpoint that this rule exists to obtain.
            double delta = 0.5 * (b - a) * (2.0 / (Math.exp(2.0 * u) + 1.0));
            double right = b - delta;
            double left = a + delta;
            boolean rightOk = right < b;
            boolean leftOk = left > a;
            if (!rightOk && !leftOk) {
                return false;
            }
            // at t = 0 both sides land on the midpoint, so the left one is
            // skipped there rather than counted twice
            double s = 0.0;
            double as = 0.0;
            if (rightOk) {
                double v = value(f, right, acc);
                s += v;
                as += Math.abs(v);
            }
            if (leftOk && !center) {
                double v = value(f, left, acc);
                s += v;
                as += Math.abs(v);
            }
            acc.sum += w * s;
            acc.absSum += w * as;
            return true;
        }

        if (kind == WHOLE) {
            double x = Math.sinh(u);
            double w = HALF_PI * coshT * Math.cosh(u);
            if (Double.isInfinite(x) || Double.isInfinite(w)) {
                return false;
            }
            if (center) {
                double v = value(f, x, acc);
                acc.sum += w * v;
                acc.absSum += w * Math.abs(v);
                return true;
            }
            double vp = value(f, x, acc);
            double vm = value(f, -x, acc);
            acc.sum += w * (vp + vm);
            acc.absSum += w * (Math.abs(vp) + Math.abs(vm));
            return true;
        }

        double e = Math.exp(u);
        double w = HALF_PI * coshT * e;
        double x = a + e;
        if (Double.isInfinite(x) || Double.isInfinite(w) || !(w > 0.0) || x <= a) {
            return false;
        }
        double v = value(f, x, acc);
        acc.sum += w * v;
        acc.absSum += w * Math.abs(v);
        return true;
    }

    /**
     * A node the integrand cannot evaluate contributes nothing. Out in the tails
     * an overflow times an underflow is a routine {@code NaN} that says nothing
     * about the integral, and letting it through would poison the whole sum.
     */
    private static double value(DFunction f, double x, Accumulator acc) {
        acc.calls++;
        double v = f.apply(x);
        return (Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v;
    }

    // =========================================================================
    // SELF CHECK
    // =========================================================================

    /**
     * Integrates a handful of integrands whose values are known in closed form
     * and reports the digits reached and what they cost.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        double tol = 1.0e-13;
        double inf = Double.POSITIVE_INFINITY;
        double nInf = Double.NEGATIVE_INFINITY;
        Object[][] cases = {
                { "1/sqrt(x) on [0,1]", (DFunction) (x -> 1.0 / Math.sqrt(x)), 0.0, 1.0, 2.0 },
                { "log(x) on [0,1]", (DFunction) (x -> Math.log(x)), 0.0, 1.0, -1.0 },
                { "x^-0.9 on [0,1]", (DFunction) (x -> Math.pow(x, -0.9)), 0.0, 1.0, 10.0 },
                { "exp(x) on [-1,1]", (DFunction) (x -> Math.exp(x)), -1.0, 1.0, Math.E - 1.0 / Math.E },
                { "exp(-x)/sqrt(x) on [0,inf)", (DFunction) (x -> Math.exp(-x) / Math.sqrt(x)), 0.0, inf,
                        Math.sqrt(Math.PI) },
                { "1/(1+x*x) on [0,inf)", (DFunction) (x -> 1.0 / (1.0 + x * x)), 0.0, inf, 0.5 * Math.PI },
                { "1/(1+x*x) on (-inf,0]", (DFunction) (x -> 1.0 / (1.0 + x * x)), nInf, 0.0, 0.5 * Math.PI },
                { "exp(-x*x) on R", (DFunction) (x -> Math.exp(-x * x)), nInf, inf, Math.sqrt(Math.PI) } };

        boolean ok = true;
        System.out.println(String.format(Locale.ROOT, "%-28s %22s %10s %8s %6s", "integrand", "value", "rel err",
                "evals", "level"));
        for (Object[] c : cases) {
            DFunction f = (DFunction) c[1];
            double a = ((Double) c[2]).doubleValue();
            double b = ((Double) c[3]).doubleValue();
            double want = ((Double) c[4]).doubleValue();
            IntegralResult r = integrate1D(f, a, b, tol);
            double rel = Math.abs(r.value - want) / Math.abs(want);
            ok &= r.converged && rel < 1.0e-12;
            System.out.println(String.format(Locale.ROOT, "%-28s %22.15f %10.2e %8d %6d", c[0], r.value, rel,
                    r.evaluations, r.level));
        }

        // an integrand that is zero wherever the rule looks must not be certified
        IntegralResult zero = integrate1D(x -> 0.0, 0.0, 1.0, tol);
        System.out.println(String.format(Locale.ROOT, "%-28s %22.15f %10s %8d %6s", "the zero function on [0,1]",
                zero.value, "-", zero.evaluations, zero.converged ? "certified" : "refused"));
        ok &= !zero.converged && zero.value == 0.0;

        System.out.println(ok ? ">>> OK" : ">>> FAILED");
    }
}
