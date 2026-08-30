package math.solve;

import java.util.Locale;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * One entry point for numerical integration in one, two and three dimensions,
 * over finite and over infinite domains alike.
 * <p>
 * {@link #integrate(DFunction, double, double, double)} and its siblings hand
 * the work to {@link InfiniteIntegrator}, which substitutes away whichever
 * limits are infinite and passes a domain that is finite throughout on to
 * {@link MetaIntegrator}, where a Clenshaw-Curtis rule and an adaptive
 * Gauss-Kronrod subdivision compete for it. The dimension follows from the
 * integrand, and the rule parameter those classes take - an enum with one
 * constant - does not appear here.
 * <p>
 * The rule that is <em>not</em> reached by any of that is the
 * double-exponential one, because deciding by heuristic whether an integrand
 * has an endpoint singularity is guesswork.
 * {@link #integrateSingular(DFunction, double, double, double)} lets a caller
 * who knows say so instead.
 * <p>
 * Every method here returns a bare {@code double}. {@link WithError} carries
 * the same seven signatures and returns what the rule thought of its own
 * answer.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Numerical_integration">Numerical
 *      integration</a>
 * @since 1.5.2
 */
public final class Quadrature {

    /** The one rule the layers underneath accept. */
    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;

    /**
     * What one integration returned and what the rule thought of it.
     *
     * @since 1.5.3
     */
    public static final class Result {

        /** The approximated integral, the same value the bare form returns. */
        public final double value;
        /**
         * The error estimate of whichever rule answered: the sum of the panel
         * estimates for a Gauss-Kronrod subdivision, the change from the second
         * finest to the finest grid for Clenshaw-Curtis and for the
         * double-exponential rule. Read {@link WithError} before believing it.
         */
        public final double approximatedErrorEstimate;
        /**
         * Whether the rule stopped because it met the tolerance it was given,
         * rather than because it ran out of budget or levels.
         * <p>
         * The name is deliberately narrow. {@code false} is a reliable warning
         * that something did not resolve. {@code true} says only that the
         * algorithm terminated on its own success criterion, which is not the
         * same as the answer being right - see {@link WithError} for how far
         * apart those two can be. The rules underneath call this same flag
         * {@code converged}; here it is named for what it actually tests,
         * because this is where a caller outside the package reads it.
         */
        public final boolean toleranceMet;

        Result(double value, double approximatedErrorEstimate, boolean toleranceMet) {
            this.value = value;
            this.approximatedErrorEstimate = approximatedErrorEstimate;
            this.toleranceMet = toleranceMet;
        }

        @Override
        public String toString() {
            return String.format(Locale.ROOT, "Value: %.8f (approx. Error: %.2e, tolerance met: %b)", value,
                    approximatedErrorEstimate, toleranceMet);
        }
    }

    /**
     * The seven methods of {@link Quadrature} again, with the same parameter
     * lists, returning a {@link Result} instead of a bare {@code double}.
     *
     * <h2>What the estimate is</h2>
     * For the adaptive Gauss-Kronrod path it is the sum of the {@code |K - G|}
     * of every panel, which is what QUADPACK calls {@code abserr}. Because the
     * tolerance is halved at every level of the subdivision, the panel
     * tolerances of a full tree sum to the tolerance the caller asked for, so a
     * run that met its tolerance is bounded by that tolerance.
     *
     * <h2>What it is not</h2>
     * It is not a bound on the true error, and it is not tight. Measured over
     * 42 integrands with known values, it exceeded the true error in 40 of them
     * and by as much as six orders of magnitude, so it is a conservative
     * indication rather than an error bar.
     * <p>
     * {@link Result#toleranceMet} is the more useful half. It said
     * {@code false} on {@code 1/sqrt(x)} and {@code log(x)} over
     * {@code [0, 1]}, where the value really is wrong - the endpoint
     * singularity is what
     * {@link #integrateSingular(DFunction, double, double, double)} is for -
     * and it says {@code false} whenever a panel was handed out without meeting
     * the tolerance it was given.
     * <p>
     * <strong>It does not follow that {@code toleranceMet} means the answer is
     * right.</strong> A feature narrower than the spacing of the forced panels
     * is missed by the Kronrod rule and by the Gauss rule alike; the two then
     * agree, the estimate is tiny, and nothing notices. See the class comment
     * of {@link AdaptiveGaussKronrod}. Measured through this facade, a
     * Gaussian of width {@code 0.005} in the unit cube at a tolerance of
     * {@code 1e-6} integrates to {@code 2.80e-09} where the answer is
     * {@code 6.96e-07} - a relative error of {@code 0.996} - and reports
     * {@code toleranceMet} with an estimate of {@code 2.80e-09}. At a width of
     * {@code 0.002} the peak is missed altogether. The same Gaussian in one
     * dimension is exact to {@code 4e-16}: the forced bisections cover that
     * case per axis, and in three dimensions the price of covering it is not
     * paid by default.
     * <p>
     * The same blindness covers a second case, and it is the more common one: a
     * rule asked for the wrong domain integrates it correctly and reports so.
     * Nothing measured here can catch that, because {@code |K - G|} only ever
     * compares what was already sampled. On an infinite domain
     * {@link InfiniteIntegrator} answers it with a probe rather than with an
     * estimate, and throws when the substitution never sampled the mass; on a
     * finite one nothing probes.
     *
     * @see <a href="https://en.wikipedia.org/wiki/QUADPACK">QUADPACK</a>
     * @since 1.5.3
     */
    public static final class WithError {

        /**
         * Integrates {@code f} over {@code [a, b]}, where either limit may be
         * infinite. See
         * {@link Quadrature#integrate(DFunction, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param a
         *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
         * @param b
         *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
         * @param epsTol
         *            error tolerance
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws ArithmeticException
         *             if an infinite domain was substituted away without the
         *             substitution ever sampling the mass of the integrand
         */
        public static Result integrate(DFunction f, double a, double b, double epsTol) {
            return of(InfiniteIntegrator.integrate1DInfiniteWithError(RULE, f, a, b, epsTol));
        }

        /**
         * Integrates {@code f} over {@code [a, b]} with the substitution
         * centered on {@code center}. See
         * {@link Quadrature#integrate(DFunction, double, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param a
         *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
         * @param b
         *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
         * @param epsTol
         *            error tolerance
         * @param center
         *            where the mass of the integrand lies, or
         *            {@code Double.NaN} for none
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws IllegalArgumentException
         *             if {@code center} is infinite, or does not lie strictly
         *             inside an interval that has a finite end
         * @throws ArithmeticException
         *             if the substitution never sampled the mass of the
         *             integrand
         */
        public static Result integrate(DFunction f, double a, double b, double epsTol, double center) {
            return of(InfiniteIntegrator.integrate1DInfiniteWithError(RULE, f, a, b, epsTol, center));
        }

        /**
         * Integrates {@code f} over a rectangle whose sides may be infinite.
         * See {@link Quadrature#integrate(DBiFunction, double, double, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param ax
         *            lower limit in x, possibly infinite
         * @param bx
         *            upper limit in x, possibly infinite
         * @param ay
         *            lower limit in y, possibly infinite
         * @param by
         *            upper limit in y, possibly infinite
         * @param epsTol
         *            error tolerance
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws ArithmeticException
         *             if the substitution never sampled the mass of the
         *             integrand, or sampled it but disagrees with the same
         *             integral split around it
         */
        public static Result integrate(DBiFunction f, double ax, double bx, double ay, double by, double epsTol) {
            return of(InfiniteIntegrator.integrate2DInfiniteWithError(RULE, f, ax, bx, ay, by, epsTol));
        }

        /**
         * Integrates {@code f} over a rectangle whose sides may be infinite,
         * with each infinite axis centered where the caller says the mass lies.
         * See {@link Quadrature#integrate(DBiFunction, double, double, double, double, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param ax
         *            lower limit in x, possibly infinite
         * @param bx
         *            upper limit in x, possibly infinite
         * @param ay
         *            lower limit in y, possibly infinite
         * @param by
         *            upper limit in y, possibly infinite
         * @param epsTol
         *            error tolerance
         * @param centerX
         *            where the mass lies in x, or {@code Double.NaN} for none
         * @param centerY
         *            where the mass lies in y, or {@code Double.NaN} for none
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws IllegalArgumentException
         *             if a center is infinite, or does not lie strictly inside
         *             an axis that has a finite end
         * @throws ArithmeticException
         *             if the substitution never sampled the mass of the
         *             integrand, or sampled it but disagrees with the same
         *             integral split around it
         */
        public static Result integrate(DBiFunction f, double ax, double bx, double ay, double by, double epsTol,
                double centerX, double centerY) {
            return of(InfiniteIntegrator.integrate2DInfiniteWithError(RULE, f, ax, bx, ay, by, epsTol, centerX,
                    centerY));
        }

        /**
         * Integrates {@code f} over a box whose sides may be infinite. See
         * {@link Quadrature#integrate(DTriFunction, double, double, double, double, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param ax
         *            lower limit in x, possibly infinite
         * @param bx
         *            upper limit in x, possibly infinite
         * @param ay
         *            lower limit in y, possibly infinite
         * @param by
         *            upper limit in y, possibly infinite
         * @param az
         *            lower limit in z, possibly infinite
         * @param bz
         *            upper limit in z, possibly infinite
         * @param epsTol
         *            error tolerance
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws ArithmeticException
         *             if the substitution never sampled the mass of the
         *             integrand, or sampled it but disagrees with the same
         *             integral split around it
         */
        public static Result integrate(DTriFunction f, double ax, double bx, double ay, double by, double az,
                double bz, double epsTol) {
            return of(InfiniteIntegrator.integrate3DInfiniteWithError(RULE, f, ax, bx, ay, by, az, bz, epsTol));
        }

        /**
         * Integrates {@code f} over a box whose sides may be infinite, with
         * each infinite axis centered where the caller says the mass lies. See
         * {@link Quadrature#integrate(DTriFunction, double, double, double, double, double, double, double, double, double, double)}.
         *
         * @param f
         *            the integrand
         * @param ax
         *            lower limit in x, possibly infinite
         * @param bx
         *            upper limit in x, possibly infinite
         * @param ay
         *            lower limit in y, possibly infinite
         * @param by
         *            upper limit in y, possibly infinite
         * @param az
         *            lower limit in z, possibly infinite
         * @param bz
         *            upper limit in z, possibly infinite
         * @param epsTol
         *            error tolerance
         * @param centerX
         *            where the mass lies in x, or {@code Double.NaN} for none
         * @param centerY
         *            where the mass lies in y, or {@code Double.NaN} for none
         * @param centerZ
         *            where the mass lies in z, or {@code Double.NaN} for none
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol}
         * @throws IllegalArgumentException
         *             if a center is infinite, or does not lie strictly inside
         *             an axis that has a finite end
         * @throws ArithmeticException
         *             if the substitution never sampled the mass of the
         *             integrand, or sampled it but disagrees with the same
         *             integral split around it
         */
        public static Result integrate(DTriFunction f, double ax, double bx, double ay, double by, double az,
                double bz, double epsTol, double centerX, double centerY, double centerZ) {
            return of(InfiniteIntegrator.integrate3DInfiniteWithError(RULE, f, ax, bx, ay, by, az, bz, epsTol, centerX,
                    centerY, centerZ));
        }

        /**
         * Integrates {@code f} over {@code [a, b]} with the double-exponential
         * rule of {@link DoubleExponential}, for an integrand that is singular
         * at one of the limits. Either limit may be infinite.
         * <p>
         * Where {@link Quadrature#integrateSingular(DFunction, double, double, double)}
         * refuses a run that did not converge, because it can only hand out a
         * bare number, this one returns it with
         * {@link Result#toleranceMet} set to {@code false} and lets the caller
         * decide.
         *
         * @param f
         *            the integrand
         * @param a
         *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
         * @param b
         *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
         * @param epsTol
         *            error tolerance
         * @return the approximated integral, its estimate and whether the rule
         *         reached {@code epsTol} before it ran out of levels
         */
        public static Result integrateSingular(DFunction f, double a, double b, double epsTol) {
            DoubleExponential.IntegralResult r = DoubleExponential.integrate1D(f, a, b, epsTol);
            return new Result(r.value, r.approximatedErrorEstimate, r.converged);
        }

        private static Result of(AdaptiveGaussKronrod.IntegralResult r) {
            return new Result(r.value, r.approximatedErrorEstimate, r.converged);
        }

        private WithError() {
            throw new AssertionError();
        }
    }

    /**
     * Integrates {@code f} over {@code [a, b]}, where either limit may be
     * infinite.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
     * @param b
     *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
     * @param epsTol
     *            error tolerance
     * @return the approximated integral
     * @throws ArithmeticException
     *             if an infinite domain was substituted away without the
     *             substitution ever sampling the mass of the integrand
     */
    public static double integrate(DFunction f, double a, double b, double epsTol) {
        return InfiniteIntegrator.integrate1DInfinite(RULE, f, a, b, epsTol);
    }

    /**
     * Integrates {@code f} over {@code [a, b]} with the substitution centered
     * on {@code center}, for mass so far from the origin that no probe can find
     * it. See
     * {@link InfiniteIntegrator#integrate1DInfinite(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double, double)}.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
     * @param b
     *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
     * @param epsTol
     *            error tolerance
     * @param center
     *            where the mass of the integrand lies, or {@code Double.NaN}
     *            for none
     * @return the approximated integral
     * @throws IllegalArgumentException
     *             if {@code center} is infinite, or does not lie strictly
     *             inside an interval that has a finite end
     * @throws ArithmeticException
     *             if the substitution never sampled the mass of the integrand
     */
    public static double integrate(DFunction f, double a, double b, double epsTol, double center) {
        return InfiniteIntegrator.integrate1DInfinite(RULE, f, a, b, epsTol, center);
    }

    /**
     * Integrates {@code f} over a rectangle whose sides may be infinite.
     *
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x, possibly infinite
     * @param bx
     *            upper limit in x, possibly infinite
     * @param ay
     *            lower limit in y, possibly infinite
     * @param by
     *            upper limit in y, possibly infinite
     * @param epsTol
     *            error tolerance
     * @return the approximated integral
     * @throws ArithmeticException
     *             if the substitution never sampled the mass of the integrand,
     *             or sampled it but disagrees with the same integral split
     *             around it
     */
    public static double integrate(DBiFunction f, double ax, double bx, double ay, double by, double epsTol) {
        return InfiniteIntegrator.integrate2DInfinite(RULE, f, ax, bx, ay, by, epsTol);
    }

    /**
     * Integrates {@code f} over a rectangle whose sides may be infinite, with
     * each infinite axis centered where the caller says the mass lies.
     *
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x, possibly infinite
     * @param bx
     *            upper limit in x, possibly infinite
     * @param ay
     *            lower limit in y, possibly infinite
     * @param by
     *            upper limit in y, possibly infinite
     * @param epsTol
     *            error tolerance
     * @param centerX
     *            where the mass lies in x, or {@code Double.NaN} for none
     * @param centerY
     *            where the mass lies in y, or {@code Double.NaN} for none
     * @return the approximated integral
     * @throws IllegalArgumentException
     *             if a center is infinite, or does not lie strictly inside an
     *             axis that has a finite end
     * @throws ArithmeticException
     *             if the substitution never sampled the mass of the integrand,
     *             or sampled it but disagrees with the same integral split
     *             around it
     */
    public static double integrate(DBiFunction f, double ax, double bx, double ay, double by, double epsTol,
            double centerX, double centerY) {
        return InfiniteIntegrator.integrate2DInfinite(RULE, f, ax, bx, ay, by, epsTol, centerX, centerY);
    }

    /**
     * Integrates {@code f} over a box whose sides may be infinite.
     *
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x, possibly infinite
     * @param bx
     *            upper limit in x, possibly infinite
     * @param ay
     *            lower limit in y, possibly infinite
     * @param by
     *            upper limit in y, possibly infinite
     * @param az
     *            lower limit in z, possibly infinite
     * @param bz
     *            upper limit in z, possibly infinite
     * @param epsTol
     *            error tolerance
     * @return the approximated integral
     * @throws ArithmeticException
     *             if the substitution never sampled the mass of the integrand,
     *             or sampled it but disagrees with the same integral split
     *             around it
     */
    public static double integrate(DTriFunction f, double ax, double bx, double ay, double by, double az, double bz,
            double epsTol) {
        return InfiniteIntegrator.integrate3DInfinite(RULE, f, ax, bx, ay, by, az, bz, epsTol);
    }

    /**
     * Integrates {@code f} over a box whose sides may be infinite, with each
     * infinite axis centered where the caller says the mass lies.
     *
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x, possibly infinite
     * @param bx
     *            upper limit in x, possibly infinite
     * @param ay
     *            lower limit in y, possibly infinite
     * @param by
     *            upper limit in y, possibly infinite
     * @param az
     *            lower limit in z, possibly infinite
     * @param bz
     *            upper limit in z, possibly infinite
     * @param epsTol
     *            error tolerance
     * @param centerX
     *            where the mass lies in x, or {@code Double.NaN} for none
     * @param centerY
     *            where the mass lies in y, or {@code Double.NaN} for none
     * @param centerZ
     *            where the mass lies in z, or {@code Double.NaN} for none
     * @return the approximated integral
     * @throws IllegalArgumentException
     *             if a center is infinite, or does not lie strictly inside an
     *             axis that has a finite end
     * @throws ArithmeticException
     *             if the substitution never sampled the mass of the integrand,
     *             or sampled it but disagrees with the same integral split
     *             around it
     */
    public static double integrate(DTriFunction f, double ax, double bx, double ay, double by, double az, double bz,
            double epsTol, double centerX, double centerY, double centerZ) {
        return InfiniteIntegrator.integrate3DInfinite(RULE, f, ax, bx, ay, by, az, bz, epsTol, centerX, centerY,
                centerZ);
    }

    /**
     * Integrates {@code f} over {@code [a, b]} with the double-exponential
     * rule of {@link DoubleExponential}, whose nodes crowd towards the ends of
     * the interval, for an integrand that is singular at one of them. Either
     * limit may be infinite.
     * <p>
     * That rule reports whether it reached the tolerance and this method
     * returns a bare number, so a run that did not converge is refused rather
     * than handed out without its flag.
     * {@link WithError#integrateSingular(DFunction, double, double, double)}
     * hands it out with the flag instead.
     *
     * @param f
     *            the integrand
     * @param a
     *            lower limit, possibly {@code Double.NEGATIVE_INFINITY}
     * @param b
     *            upper limit, possibly {@code Double.POSITIVE_INFINITY}
     * @param epsTol
     *            error tolerance
     * @return the approximated integral
     * @throws ArithmeticException
     *             if the rule exhausted its levels without reaching
     *             {@code epsTol}
     */
    public static double integrateSingular(DFunction f, double a, double b, double epsTol) {
        DoubleExponential.IntegralResult r = DoubleExponential.integrate1D(f, a, b, epsTol);
        if (!r.converged) {
            throw new ArithmeticException(String.format(Locale.ROOT,
                    "the double-exponential rule did not reach %.3e: it stopped at level %d after %d evaluations "
                            + "with %.10g and an estimated error of %.3e. Ask for a tolerance it can reach, or "
                            + "integrate the piece that is singular separately.",
                    Double.valueOf(epsTol), Integer.valueOf(r.level), Long.valueOf(r.evaluations),
                    Double.valueOf(r.value), Double.valueOf(r.approximatedErrorEstimate)));
        }
        return r.value;
    }

    private Quadrature() {
        throw new AssertionError();
    }
}
