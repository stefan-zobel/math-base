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
 *
 * @see <a href="https://en.wikipedia.org/wiki/Numerical_integration">Numerical
 *      integration</a>
 * @since 1.5.2
 */
public final class Quadrature {

    /** The one rule the layers underneath accept. */
    private static final AdaptiveGaussKronrod.G7_K15 RULE = AdaptiveGaussKronrod.G7_K15.POINTS_15;

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
