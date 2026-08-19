package math.optim;

import math.MathConsts;
import math.fun.DFunction;
import math.fun.DiffDFunction;

/**
 * Minimization of a function of one variable by Brent's method, with the golden
 * section search that brackets the minimum first. Both are taken from Numerical
 * Recipes, where they are {@code mnbrak}, {@code brent} and {@code dbrent}. See
 * <a href="https://en.wikipedia.org/wiki/Brent%27s_method">Brent's method</a>.
 * <p>
 * This is the minimizer counterpart to
 * {@link math.solve.RootFinder#brentDekker(double, double, DFunction, double)},
 * which is Brent's method for <em>roots</em>. Minimization is used in two
 * steps: {@link #bracket(DFunction, double, double)} finds a triple
 * {@code a < b < c} with {@code f(b) <= min(f(a), f(c))}, which is what
 * guarantees a minimum in between, and one of the {@code minimize} methods then
 * narrows it down.
 * <p>
 * The default tolerance is {@code sqrt(MACH_EPS_DBL)}, about {@code 1.05e-8}.
 * Nothing smaller is worth asking for: a function is quadratically flat at its
 * minimum, so the abscissa cannot be resolved more finely than the square root
 * of the machine precision no matter which method is used.
 *
 * @since 1.5.2
 */
public final class BrentMinimizer {

    /** Accuracy goal if none is given, the square root of the machine epsilon. */
    private static final double DEFAULT_TOL = Math.sqrt(MathConsts.MACH_EPS_DBL);
    /** Iteration budget of one minimization if none is given. */
    private static final int DEFAULT_MAX_ITERATIONS = 100;
    /** Number of expansion steps the bracketing may take if none is given. */
    private static final int DEFAULT_MAX_BRACKET_STEPS = 200;

    /** Default ratio by which the bracketing interval grows. */
    private static final double GOLD = 1.618034;
    /** Largest factor by which a parabolic step may leave the interval. */
    private static final double GLIMIT = 100.0;
    /** Guards the parabolic step against a division by zero. */
    private static final double TINY = 1.0e-20;
    /** The golden ratio complement, {@code (3 - sqrt(5)) / 2}. */
    private static final double CGOLD = 0.3819660112501051;
    /**
     * Absolute floor under the convergence tolerance. Without it the tolerance
     * is purely relative and collapses to zero for a minimum that happens to
     * lie at {@code x == 0}, where the convergence test can then never fire.
     */
    private static final double ZEPS = 1.0e-10;

    private final double tol;
    private final int maxIterations;
    private final int maxBracketSteps;

    /** An interval known to contain a minimum, and the values at its ends. */
    public static final class Bracket {

        /** Left end. */
        public final double a;
        /** Interior point, {@code fb <= min(fa, fc)}. */
        public final double b;
        /** Right end. */
        public final double c;
        /** {@code f(a)}. */
        public final double fa;
        /** {@code f(b)}. */
        public final double fb;
        /** {@code f(c)}. */
        public final double fc;
        /**
         * Whether a minimum was actually bracketed. {@code false} means the
         * search left the finite range or ran out of expansion steps, which is
         * what a function without a minimum in this direction looks like; the
         * other fields then hold the last point reached and must not be handed
         * to {@code minimize}.
         */
        public final boolean bracketed;

        Bracket(double a, double b, double c, double fa, double fb, double fc, boolean bracketed) {
            this.a = a;
            this.b = b;
            this.c = c;
            this.fa = fa;
            this.fb = fb;
            this.fc = fc;
            this.bracketed = bracketed;
        }

        @Override
        public String toString() {
            return String.format("Bracket(%.6g, %.6g, %.6g)%s", a, b, c, bracketed ? "" : " NOT BRACKETED");
        }
    }

    /** Where a minimization settled. */
    public static final class Result {

        /** The abscissa of the minimum. */
        public final double x;
        /** {@code f(x)}. */
        public final double value;
        /** Iterations performed. */
        public final int iterations;
        /**
         * Whether the tolerance was met. {@code false} means the iteration
         * budget ran out first, and {@link #x} is only the best point seen.
         */
        public final boolean converged;

        Result(double x, double value, int iterations, boolean converged) {
            this.x = x;
            this.value = value;
            this.iterations = iterations;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return String.format("Result(x=%.10g, f=%.10g, %d iterations)%s", x, value, iterations,
                    converged ? "" : " NOT CONVERGED");
        }
    }

    /**
     * Creates a minimizer with the default accuracy goal and budgets.
     */
    public BrentMinimizer() {
        this(DEFAULT_TOL, DEFAULT_MAX_ITERATIONS, DEFAULT_MAX_BRACKET_STEPS);
    }

    /**
     * Creates a minimizer with explicit settings.
     *
     * @param tol
     *            relative accuracy goal for the abscissa of the minimum,
     *            greater than {@code 0}; asking for less than
     *            {@code sqrt(MACH_EPS_DBL)} does not make the answer better
     * @param maxIterations
     *            iteration budget of one minimization, {@code 1} or greater;
     *            exhausting it ends the search without reporting convergence
     * @param maxBracketSteps
     *            number of expansion steps the bracketing may take, {@code 1}
     *            or greater
     */
    public BrentMinimizer(double tol, int maxIterations, int maxBracketSteps) {
        if (!(tol > 0.0) || Double.isInfinite(tol)) {
            throw new IllegalArgumentException("tol must be finite and positive : " + tol);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException("maxIterations must be >= 1 : " + maxIterations);
        }
        if (maxBracketSteps < 1) {
            throw new IllegalArgumentException("maxBracketSteps must be >= 1 : " + maxBracketSteps);
        }
        this.tol = tol;
        this.maxIterations = maxIterations;
        this.maxBracketSteps = maxBracketSteps;
    }

    /**
     * Searches downhill from {@code a} towards {@code b} until the function
     * turns back up, and returns the triple that straddles the minimum.
     *
     * @param f
     *            the function to bracket a minimum of
     * @param a
     *            one end of the initial interval
     * @param b
     *            the other end, which sets the direction and the initial step
     * @return the bracket; check {@link Bracket#bracketed} before using it
     */
    public Bracket bracket(DFunction f, double a, double b) {
        if (f == null) {
            throw new IllegalArgumentException("f is null");
        }
        if (!isFinite(a)) {
            throw new IllegalArgumentException("a is not finite : " + a);
        }
        if (!isFinite(b)) {
            throw new IllegalArgumentException("b is not finite : " + b);
        }
        if (a == b) {
            throw new IllegalArgumentException("a and b must differ, both are " + a);
        }

        double ax = a;
        double fa = f.apply(ax);
        double bx = b;
        double fb = f.apply(bx);

        if (fb > fa) {
            // walk downhill, so swap if we picked the wrong direction
            double tmp = fa;
            fa = fb;
            fb = tmp;
            tmp = ax;
            ax = bx;
            bx = tmp;
        }

        double cx = bx + GOLD * (bx - ax);
        double fc = f.apply(cx);

        for (int step = 0; step < maxBracketSteps; step++) {
            if (!(fb > fc)) {
                return bracketOf(ax, bx, cx, fa, fb, fc);
            }
            // Leaving the finite range means the function keeps falling and no
            // minimum exists in this direction. Continuing would iterate on
            // infinities until the comparison above happens to be false for a
            // NaN, which looks exactly like success.
            if (!isFinite(cx) || !isFinite(fc)) {
                return new Bracket(ax, bx, cx, fa, fb, fc, false);
            }

            double r = (bx - ax) * (fb - fc);
            double q = (bx - cx) * (fb - fa);
            double u = bx - ((bx - cx) * q - (bx - ax) * r)
                    / (2.0 * sign(Math.max(Math.abs(q - r), TINY), q - r));
            double fu;
            double ulim = bx + GLIMIT * (cx - bx);

            if ((bx - u) * (u - cx) > 0.0) {
                // parabolic minimum between b and c
                fu = f.apply(u);
                if (fu < fc) {
                    return bracketOf(bx, u, cx, fb, fu, fc);
                } else if (fu > fb) {
                    return bracketOf(ax, bx, u, fa, fb, fu);
                }
                u = cx + GOLD * (cx - bx);
                fu = f.apply(u);
            } else if ((cx - u) * (u - ulim) > 0.0) {
                // parabolic minimum beyond c, within the allowed limit
                fu = f.apply(u);
                if (fu < fc) {
                    bx = cx;
                    cx = u;
                    u = cx + GOLD * (cx - bx);
                    fb = fc;
                    fc = fu;
                    fu = f.apply(u);
                }
            } else if ((u - ulim) * (ulim - cx) >= 0.0) {
                // do not take a parabolic step further than the limit
                u = ulim;
                fu = f.apply(u);
            } else {
                // reject the parabolic step and grow by the golden ratio
                u = cx + GOLD * (cx - bx);
                fu = f.apply(u);
            }

            ax = bx;
            bx = cx;
            cx = u;
            fa = fb;
            fb = fc;
            fc = fu;
        }

        return new Bracket(ax, bx, cx, fa, fb, fc, false);
    }

    /**
     * Narrows a bracket down to the minimum using the derivative, which is
     * Numerical Recipes' {@code dbrent}. Slightly fewer function evaluations
     * than the derivative-free variant, at the cost of one derivative
     * evaluation per iteration.
     *
     * @param f
     *            the function to minimize
     * @param bracket
     *            a bracket from {@link #bracket(DFunction, double, double)},
     *            with {@link Bracket#bracketed} set
     * @return the minimum found, see {@link Result}
     */
    public Result minimize(DiffDFunction f, Bracket bracket) {
        checkArguments(f, bracket);

        double ax = bracket.a;
        double cx = bracket.c;
        double a = (ax < cx) ? ax : cx;
        double b = (ax > cx) ? ax : cx;

        double d = 0.0;
        double e = 0.0;

        double x = bracket.b;
        double v = x;
        double w = x;
        double fx = bracket.fb;
        double fv = fx;
        double fw = fx;
        double dx = f.derivativeAt(x);
        double dv = dx;
        double dw = dx;

        for (int iteration = 1; iteration <= maxIterations; iteration++) {
            double xm = 0.5 * (a + b);
            double tol1 = tol * Math.abs(x) + ZEPS;
            double tol2 = 2.0 * tol1;
            if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                return new Result(x, fx, iteration, true);
            }

            double u;
            if (Math.abs(e) > tol1) {
                // secant step from each of the two stored derivatives
                double d1 = 2.0 * (b - a);
                double d2 = d1;
                if (dw != dx) {
                    d1 = (w - x) * dx / (dx - dw);
                }
                if (dv != dx) {
                    d2 = (v - x) * dx / (dx - dv);
                }
                double u1 = x + d1;
                double u2 = x + d2;
                boolean ok1 = ((a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0);
                boolean ok2 = ((a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0);
                double olde = e;
                e = d;
                if (ok1 || ok2) {
                    if (ok1 && ok2) {
                        d = (Math.abs(d1) < Math.abs(d2)) ? d1 : d2;
                    } else if (ok1) {
                        d = d1;
                    } else {
                        d = d2;
                    }
                    if (Math.abs(d) <= Math.abs(0.5 * olde)) {
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2) {
                            d = sign(tol1, xm - x);
                        }
                    } else {
                        e = (dx >= 0.0) ? (a - x) : (b - x);
                        d = 0.5 * e;
                    }
                } else {
                    e = (dx >= 0.0) ? (a - x) : (b - x);
                    d = 0.5 * e;
                }
            } else {
                e = (dx >= 0.0) ? (a - x) : (b - x);
                d = 0.5 * e;
            }

            double fu;
            if (Math.abs(d) >= tol1) {
                u = x + d;
                fu = f.apply(u);
            } else {
                u = x + sign(tol1, d);
                fu = f.apply(u);
                if (fu > fx) {
                    // a step of the smallest size already goes uphill
                    return new Result(x, fx, iteration, true);
                }
            }

            double du = f.derivativeAt(u);
            if (fu <= fx) {
                if (u >= x) {
                    a = x;
                } else {
                    b = x;
                }
                v = w;
                fv = fw;
                dv = dw;
                w = x;
                fw = fx;
                dw = dx;
                x = u;
                fx = fu;
                dx = du;
            } else {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
                if (fu <= fw || w == x) {
                    v = w;
                    fv = fw;
                    dv = dw;
                    w = u;
                    fw = fu;
                    dw = du;
                } else if (fu < fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                    dv = du;
                }
            }
        }

        return new Result(x, fx, maxIterations, false);
    }

    /**
     * Narrows a bracket down to the minimum without derivatives, by golden
     * section search with parabolic interpolation, which is Numerical Recipes'
     * {@code brent}.
     * <p>
     * Note that a {@link DiffDFunction} passed through a variable of this
     * static type selects this method, not the one that uses the derivative.
     *
     * @param f
     *            the function to minimize
     * @param bracket
     *            a bracket from {@link #bracket(DFunction, double, double)},
     *            with {@link Bracket#bracketed} set
     * @return the minimum found, see {@link Result}
     */
    public Result minimize(DFunction f, Bracket bracket) {
        checkArguments(f, bracket);

        double ax = bracket.a;
        double cx = bracket.c;
        double a = (ax < cx) ? ax : cx;
        double b = (ax > cx) ? ax : cx;

        double d = 0.0;
        double e = 0.0;

        double x = bracket.b;
        double v = x;
        double w = x;
        double fx = bracket.fb;
        double fv = fx;
        double fw = fx;

        for (int iteration = 1; iteration <= maxIterations; iteration++) {
            double xm = 0.5 * (a + b);
            double tol1 = tol * Math.abs(x) + ZEPS;
            double tol2 = 2.0 * tol1;
            if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                return new Result(x, fx, iteration, true);
            }

            if (Math.abs(e) > tol1) {
                // fit a parabola through x, v and w
                double r = (x - w) * (fx - fv);
                double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) {
                    p = -p;
                }
                q = Math.abs(q);
                double olde = e;
                e = d;
                if (Math.abs(p) >= Math.abs(0.5 * q * olde) || p <= q * (a - x) || p >= q * (b - x)) {
                    // the parabolic step is unacceptable, take a golden one
                    e = (x >= xm) ? (a - x) : (b - x);
                    d = CGOLD * e;
                } else {
                    d = p / q;
                    double u = x + d;
                    if (u - a < tol2 || b - u < tol2) {
                        d = sign(tol1, xm - x);
                    }
                }
            } else {
                e = (x >= xm) ? (a - x) : (b - x);
                d = CGOLD * e;
            }

            double u = (Math.abs(d) >= tol1) ? (x + d) : (x + sign(tol1, d));
            double fu = f.apply(u);

            if (fu <= fx) {
                if (u >= x) {
                    a = x;
                } else {
                    b = x;
                }
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            } else {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
                if (fu <= fw || w == x) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }

        return new Result(x, fx, maxIterations, false);
    }

    private static void checkArguments(DFunction f, Bracket bracket) {
        if (f == null) {
            throw new IllegalArgumentException("f is null");
        }
        if (bracket == null) {
            throw new IllegalArgumentException("bracket is null");
        }
        if (!bracket.bracketed) {
            throw new IllegalArgumentException(
                    "no minimum was bracketed, so there is nothing to narrow down : " + bracket);
        }
    }

    /**
     * A valid bracket, or one flagged as invalid if the values do not actually
     * straddle a minimum.
     */
    private static Bracket bracketOf(double a, double b, double c, double fa, double fb, double fc) {
        boolean ok = isFinite(a) && isFinite(b) && isFinite(c) && isFinite(fa) && isFinite(fb) && isFinite(fc)
                && fb <= fa && fb <= fc;
        return new Bracket(a, b, c, fa, fb, fc, ok);
    }

    /** The magnitude of {@code x} with the sign of {@code y}. */
    private static double sign(double x, double y) {
        return (y >= 0.0) ? Math.abs(x) : -Math.abs(x);
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
