package math.solve;

import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLong;

import math.fun.DBiFunction;
import math.fun.DFunction;
import math.fun.DTriFunction;

/**
 * Integration over an infinite or semi-infinite interval, by an algebraic
 * substitution that maps it onto a finite one.
 *
 * <h2>The substitution can miss the integrand entirely</h2>
 * {@code x = t / (1 - t^2)} carries the whole line onto {@code [-1, 1]}, and
 * mass that sits far from the origin arrives as a spike pressed against
 * {@code |t| = 1}. Whether the subdivision finds that spike turns on how heavy
 * the tail is and on how far away the mass is <em>measured in its own width</em>,
 * not on the distance alone. A Cauchy centered at a thousand is integrated
 * exactly. A normal of unit width centered at a hundred, left to itself, would
 * come back as {@code 6e-101}, and one of width ten centered at three hundred as
 * {@code 3e-60}, where both ought to be one -- those are the values the bare
 * substitution produces, and neither is handed out any more: in one dimension
 * both now come back as one, from the fallback described below, and where that
 * fallback cannot answer either the call is refused. The one-sided forms
 * substitute {@code x = a + t / (1 - t)} and fail the same way when the mass
 * lies far from the finite end.
 * <p>
 * The multidimensional forms fail sooner, because they are allowed less
 * recursion per axis: {@code integrate2DSmart} spends a {@code maxDepth} of ten
 * on both axes together and {@code integrate3DSmart} fourteen on three. Left to
 * itself a product of two unit normals centered at {@code (10, 10)} yields
 * {@code 1.0042} and at {@code (30, 30)} {@code 1.07e-33}, where the 1D form is
 * still exact at thirty.
 *
 * <h2>None of them returns such an answer quietly any more</h2>
 * Whenever at least one limit is infinite -- a call with finite limits
 * throughout involves no substitution and is passed straight through unchecked
 * -- the largest {@code |f|} the substitution actually evaluated is recorded and
 * compared against a direct probe of the integrand. If the substitution never
 * came within a millionth of what the probe found, it looked in the wrong
 * place.
 * <p>
 * In one dimension that is not the end of it. The probe knows <em>where</em>
 * the mass is, and splitting the domain there is exactly what the refusal used
 * to tell the caller to do; so the split is made here instead and both halves
 * go to {@link DoubleExponential}, whose nodes crowd towards the finite end of
 * each half - which is now the peak itself. The result is used only if both
 * halves report convergence, which is also the check on the probe: handed a
 * split point a hundred widths off the mass, that rule declines rather than
 * agrees. When it declines, an {@link ArithmeticException} names where the mass
 * is and what to do instead, as before. In two and three dimensions there is no
 * such fallback, because the rule exists in one dimension only, and the refusal
 * is the whole of the answer there.
 * <p>
 * The probe walks a geometric ladder of magnitudes outwards from whichever end
 * of an axis is finite, or from the origin when neither is; an axis bounded on
 * both sides is scanned uniformly instead. In 2D and 3D the axes are combined
 * into a tensor grid, which only has to find <em>one</em> point where the
 * integrand is not zero, because coordinate sweeps then climb from there to the
 * peak. In 1D a local refinement around the best rung does the same job. That
 * step is what makes any of it work: the rungs are a fixed <em>fraction</em>
 * apart, so a narrow feature far out falls between two of them, and for a
 * unit-width normal at three hundred the nearest rung is fifteen widths off and
 * reports {@code 6e-51} where the density is {@code 0.4}.
 * <p>
 * In 2D and 3D a second check follows, because there the substitution can
 * sample the mass and still integrate it badly, which no probe can detect: when
 * the located peak sits more than three away from the point the substitution
 * is centered on, the domain is split there into orthants and the two results
 * are compared. That is what catches the {@code (10, 10)} case, whose error is
 * only four parts in a thousand. It costs four further integrations in 2D and
 * eight in 3D, and is reached only when the mass really does sit away from the
 * center.
 *
 * <h2>What it catches, measured</h2>
 * Over the batteries collected while this was written, and with no false alarm
 * anywhere -- an integrand that cancels to zero, the zero function, an integral
 * that genuinely is {@code 1e-200} and heavy tails all pass untouched:
 * <ul>
 *   <li><b>1D</b>: of 9 known failures, 8 are now <em>answered</em> to nine
 *       digits or better by the fallback and the ninth - an oscillation riding
 *       on mass that sits far out, which no double-exponential node set
 *       resolves - is refused. None is returned wrong. 0 of 9 legitimate
 *       integrals disturbed. The probe costs about 1,200 evaluations of the
 *       integrand per call, and the fallback a few hundred to a few thousand
 *       more, on the path that was about to refuse.</li>
 *   <li><b>2D</b>: 9 of 11 refused, 0 of 6 disturbed, about 50,000
 *       evaluations.</li>
 *   <li><b>3D</b>: 5 of 6 refused, 0 of 4 disturbed, about 410,000
 *       evaluations. A 3D integration over an infinite domain costs at least
 *       216,000 by itself, so this is of the same order as the work it
 *       guards.</li>
 * </ul>
 * What escapes is a feature narrower than the ladder can resolve. A unit-width
 * normal at {@code 1e4} in 1D, or at {@code (1000, 1000)} in 2D, underflows to
 * exactly zero at every rung, so there is nothing to refine around and the old
 * silent zero remains. Detection by sampling has that limit in principle, and
 * it bites sooner in more dimensions because the integrand is a product: the
 * rung nearest the mass is off by the same relative amount on every axis at
 * once. No probe can help that, since the fallback is reached only through a
 * probe that found something: what the probe cannot see, the fallback is never
 * asked about.
 *
 * <h2>Which is why the caller may name the center</h2>
 * The cure for that residue is to center the substitution on the integrand
 * rather than to guard it. In one dimension splitting at a peak the probe found
 * is that cure, and it is taken automatically -- but only for a peak the probe
 * <em>did</em> find. Where it found nothing, the overloads taking a
 * {@code center} let a caller say where the mass is, in one, two and three
 * dimensions, per axis. The substitution then becomes
 * {@code x = center + t/(1 - t^2)} whatever the ends of the interval are, whose
 * derivative at the center is one: a feature of unit width there keeps its
 * width, which is what makes it resolvable. A finite end is not dropped but
 * carried over into the limit of {@code t} that maps onto it. Named, the normal
 * at {@code 1e4} above comes back to thirteen digits.
 * <p>
 * Centering is never <em>inferred</em> in two or three dimensions: an axis with
 * no center is substituted exactly as it was. And a center is a promise rather
 * than a licence -- the probe samples it in addition to everything it sampled
 * before, one evaluation per axis, so a center in the wrong place leaves every
 * guard here standing. In one dimension it does better than refuse: where the
 * substitution sampled the mass without resolving it, which only a center the
 * caller chose can bring about, the domain is split where the mass really is
 * and the answer comes from there.
 */
public class InfiniteIntegrator {

    /**
     * How far below the direct probe the substitution's own view of the
     * integrand may fall before the result is refused. A factor of a million
     * apart is not a matter of accuracy; it means the two never looked at the
     * same part of the domain.
     */
    private static final double MISS_FACTOR = 1.0e-6;

    /**
     * Mass this close to the point a substitution is centered on is always
     * resolved, so the multidimensional split check is not worth its cost
     * there.
     */
    private static final double SAFE_RADIUS = 3.0;

    /** Rungs of the 1D probing ladder above and below unit magnitude. */
    private static final int LADDER_STEPS = 380;

    /** Ratio between neighbouring rungs of the 1D probing ladder. */
    private static final double LADDER_RATIO = 1.2;

    /** The rung of the 1D ladder whose magnitude is one. */
    private static final int LADDER_UNIT = 200;

    /** Ratio between neighbouring rungs of the 2D probing grid. */
    private static final double GRID_RATIO_2D = 1.3;

    /** Rungs of the 2D probing grid above and below unit magnitude. */
    private static final int GRID_RUNGS_2D = 53;

    /** Ratio between neighbouring rungs of the 3D probing grid. */
    private static final double GRID_RATIO_3D = 1.9;

    /** Rungs of the 3D probing grid above and below unit magnitude. */
    private static final int GRID_RUNGS_3D = 18;

    /** Rounds of local refinement around the best rung of a 1D ladder. */
    private static final int REFINEMENTS = 2;

    /** Samples taken on each side of the best point in a refinement round. */
    private static final int REFINEMENT_SAMPLES = 100;

    /** How much each 1D refinement round narrows the window it scans. */
    private static final double REFINEMENT_NARROWING = 50.0;

    /** Sweeps over all axes after the multidimensional grid. */
    private static final int SWEEP_ROUNDS = 4;

    /** Narrowing rounds within one axis of a sweep. */
    private static final int SWEEP_REFINEMENTS = 3;

    /** How much each sweep round narrows the window it scans. */
    private static final double SWEEP_NARROWING = 20.0;

    /**
     * The largest {@code |f|} seen so far and the point it was seen at. Used by
     * the probes, which are sequential.
     */
    private static final class Peak {
        double value;
        double x;
        double y;
        double z;

        void record(double px, double v) {
            double magnitude = Math.abs(v);
            if (magnitude > value) {
                value = magnitude;
                x = px;
            }
        }

        void record(double px, double py, double v) {
            double magnitude = Math.abs(v);
            if (magnitude > value) {
                value = magnitude;
                x = px;
                y = py;
            }
        }

        void record(double px, double py, double pz, double v) {
            double magnitude = Math.abs(v);
            if (magnitude > value) {
                value = magnitude;
                x = px;
                y = py;
                z = pz;
            }
        }
    }

    /**
     * The largest {@code |f|} the substitution evaluated. The 2D and 3D rules
     * sum in parallel, so this has to be safe against concurrent updates; the
     * common path is a single read that fails the comparison.
     */
    private static final class Watch {
        private final AtomicLong bits = new AtomicLong();

        double record(double v) {
            double magnitude = Math.abs(v);
            long candidate = Double.doubleToRawLongBits(magnitude);
            long current = bits.get();
            while (magnitude > Double.longBitsToDouble(current)) {
                if (bits.compareAndSet(current, candidate)) {
                    break;
                }
                current = bits.get();
            }
            return v;
        }

        double value() {
            return Double.longBitsToDouble(bits.get());
        }
    }

    // =========================================================================
    // ONEDIMENSIONAL INFINITE INTEGRATOR (1D)
    // =========================================================================
    /**
     * Integrates {@code f} over {@code [a, b]}, where either limit may be
     * infinite. Supports {@code [-inf, +inf]}, {@code [a, +inf]},
     * {@code [-inf, b]} and, as a fallback, an ordinary finite interval.
     * <p>
     * For the three infinite cases the result is checked against a direct probe
     * of the integrand. If the substitution never sampled the region that
     * carries the mass, the domain is split where the probe found it and both
     * halves are integrated by {@link DoubleExponential} instead; only if that
     * declines to converge is the call refused. See the class comment.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand and the double-exponential
     *             fallback could not resolve it either
     */
    public static double integrate1DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                               double a, double b, double epsTol) {
        return integrate1D(setup, f, a, b, epsTol, Double.NaN);
    }

    /**
     * Integrates {@code f} over {@code [a, b]} with the substitution centered
     * on {@code center} instead of on the origin or on the finite end.
     * <p>
     * This is for the case the class comment ends on: mass that no probe can
     * find, because it underflows to zero at every rung of the ladder. A caller
     * who knows where the mass lies says so here and the substitution is built
     * around that point: {@code x = center + t/(1 - t^2)}, whose derivative at
     * the center is one, so a feature of unit width there stays unit width in
     * {@code t}. A semi-infinite interval uses the same substitution, its
     * finite end carried over into the limit of {@code t} that maps onto it,
     * because merely scaling the semi-infinite substitution would move the mass
     * to the middle and compress it by the same factor, which is not centering.
     * <p>
     * A center is a promise, and it is still checked. The probe samples it in
     * addition to everything it sampled before, so a center that is not where
     * the mass is leaves the existing guard intact rather than disarming it.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *            for none, which is what the overload without it passes. It has
     *            no effect on a finite interval, where nothing is substituted
     * @return the approximated integral
     * @throws IllegalArgumentException
     *             if {@code center} is infinite, or does not lie strictly
     *             inside an interval that has a finite end
     * @throws ArithmeticException
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand and the double-exponential
     *             fallback could not resolve it either
     * @since 1.5.2
     */
    public static double integrate1DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                               double a, double b, double epsTol, double center) {
        return integrate1D(setup, f, a, b, epsTol, center);
    }

    private static double integrate1D(AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                               double a, double b, double epsTol, double rawCenter) {

        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        final double center = checkCenter(rawCenter, a, b, "center");

        // CASE Z: any infinite interval, centered on a point the caller named
        if (!Double.isNaN(center)) {
            final Watch seen = new Watch();
            // Transformation: x = c + t / (1 - t^2), the doubly infinite one
            // moved onto the center. Its derivative at t = 0 is one, so a
            // feature of unit width at the center stays unit width in t --
            // which is what centering has to mean and what scaling the
            // semi-infinite substitution instead would not give. A finite end
            // is not dropped but carried into the limit of t that maps onto it.
            DFunction transformed = t -> {
                double t2 = t * t;
                double divisor = 1.0 - t2;
                if (Math.abs(divisor) < 1e-15) {
                    return 0.0;
                }
                double x = center + t / divisor;
                double derivative = (1.0 + t2) / (divisor * divisor);
                return seen.record(f.apply(x)) * derivative;
            };
            double lo = aInf ? -1.0 : preimage(a - center);
            double hi = bInf ? 1.0 : preimage(b - center);
            double result = MetaIntegrator.integrate1DSmart(setup, transformed, lo, hi, epsTol);
            return verify1D(f, seen, a, b, result, epsTol, center);
        }

        // CASE A: Doubly infinite [-inf, +inf]
        if (aInf && bInf) {
            final Watch seen = new Watch();
            // Transformation: x = t / (1 - t^2) on the interval [-1, 1]
            DFunction transformed = t -> {
                double t2 = t * t;
                double divisor = 1.0 - t2;
                // Guard against exact division by zero at the boundaries -1 and 1
                if (Math.abs(divisor) < 1e-15) return 0.0;

                double x = t / divisor;
                double derivative = (1.0 + t2) / (divisor * divisor);
                return seen.record(f.apply(x)) * derivative;
            };
            // We integrate the transformed function strictly from -1 to 1
            double result = MetaIntegrator.integrate1DSmart(setup, transformed, -1.0, 1.0, epsTol);
            return verify1D(f, seen, a, b, result, epsTol, center);
        }

        // CASE B: Semi-infinite upward [a, +inf]
        if (!aInf && bInf) {
            final Watch seen = new Watch();
            // Transformation: x = a + t / (1 - t) on the interval [0, 1]
            DFunction transformed = t -> {
                double divisor = 1.0 - t;
                if (Math.abs(divisor) < 1e-15) return 0.0;

                double x = a + t / divisor;
                double derivative = 1.0 / (divisor * divisor);
                return seen.record(f.apply(x)) * derivative;
            };
            double result = MetaIntegrator.integrate1DSmart(setup, transformed, 0.0, 1.0, epsTol);
            return verify1D(f, seen, a, b, result, epsTol, center);
        }

        // CASE C: Semi-infinite downward [-inf, b]
        if (aInf && !bInf) {
            final Watch seen = new Watch();
            // Transformation: x = b - t / (1 - t) on the interval [0, 1]
            DFunction transformed = t -> {
                double divisor = 1.0 - t;
                if (Math.abs(divisor) < 1e-15) return 0.0;

                double x = b - t / divisor;
                double derivative = 1.0 / (divisor * divisor);
                return seen.record(f.apply(x)) * derivative;
            };
            double result = MetaIntegrator.integrate1DSmart(setup, transformed, 0.0, 1.0, epsTol);
            return verify1D(f, seen, a, b, result, epsTol, center);
        }

        // CASE D: Ordinary finite integral [a, b]
        return MetaIntegrator.integrate1DSmart(setup, f, a, b, epsTol);
    }

    /**
     * The point {@code t} that {@code x = c + t/(1 - t^2)} carries onto a
     * finite end at distance {@code d} from the center, which is the root of
     * {@code d*t^2 + t - d} that lies in {@code (-1, 1)}. Written as
     * {@code 2d / (1 + sqrt(1 + 4d^2))} rather than as the quadratic formula:
     * that form neither divides by zero at {@code d = 0} nor overflows when
     * {@code d} is large, and {@code hypot} keeps the square root safe too.
     */
    private static double preimage(double d) {
        return (2.0 * d) / (1.0 + Math.hypot(1.0, 2.0 * d));
    }

    // =========================================================================
    // MULTIDIMENSIONAL INFINITE INTEGRATOR (2D)
    // =========================================================================
    /**
     * Integrates {@code f} over a rectangle whose sides may be infinite.
     * <p>
     * The result is checked against a probe of the integrand and, when the mass
     * sits away from the origin, against the same integral split into quadrants
     * around it; see the class comment.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand, or sampled it but disagrees with
     *             the same integral split around it
     */
    public static double integrate2DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
                                               double ax, double bx, double ay, double by, double epsTol) {
        return integrate2D(setup, f, ax, bx, ay, by, epsTol, true, Double.NaN, Double.NaN);
    }

    /**
     * Integrates {@code f} over a rectangle whose sides may be infinite, with
     * each infinite axis of the substitution centered where the caller says the
     * mass lies instead of on the origin or on the finite end.
     * <p>
     * See the single-argument sibling
     * {@link #integrate1DInfinite(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double, double)}
     * for what a center is and why supplying one does not disarm the guard.
     * Centering is never inferred here: an axis whose center is
     * {@code Double.NaN} is substituted exactly as it was before, so a caller
     * may center one axis and leave the other alone.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand, or sampled it but disagrees with
     *             the same integral split around it
     * @since 1.5.2
     */
    public static double integrate2DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
                                               double ax, double bx, double ay, double by, double epsTol,
                                               double centerX, double centerY) {
        return integrate2D(setup, f, ax, bx, ay, by, epsTol, true, centerX, centerY);
    }

    private static double integrate2D(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
                                        double ax, double bx, double ay, double by, double epsTol,
                                        boolean verify, double rawCx, double rawCy) {

        // 1. Analyze bounds
        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);

        final double cx = checkCenter(rawCx, ax, bx, "centerX");
        final double cy = checkCenter(rawCy, ay, by, "centerY");
        final boolean onX = !Double.isNaN(cx);
        final boolean onY = !Double.isNaN(cy);

        // 2. Set target integration bounds for the transformed variables
        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }

        // A centered axis uses the doubly infinite substitution about its
        // center whatever its ends are, so its finite end becomes a limit in t
        if (onX) { transAx = axInf ? -1.0 : preimage(ax - cx); transBx = bxInf ? 1.0 : preimage(bx - cx); }
        if (onY) { transAy = ayInf ? -1.0 : preimage(ay - cy); transBy = byInf ? 1.0 : preimage(by - cy); }

        final Watch seen = new Watch();

        // 3. Build dynamic wrapper function
        DBiFunction transformed = (tX, tY) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (onX) {
                double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = cx + tX / div; jX = (1.0 + tX * tX) / (div * div);
            } else if (axInf && bxInf) {
                double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = tX / div; jX = (1.0 + tX * tX) / (div * div);
            } else if (!axInf && bxInf) {
                double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = ax + tX / div; jX = 1.0 / (div * div);
            } else if (axInf && !bxInf) {
                double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0;
                x = bx - tX / div; jX = 1.0 / (div * div);
            }

            // Transform Y-axis
            double y = tY, jY = 1.0;
            if (onY) {
                double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = cy + tY / div; jY = (1.0 + tY * tY) / (div * div);
            } else if (ayInf && byInf) {
                double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = tY / div; jY = (1.0 + tY * tY) / (div * div);
            } else if (!ayInf && byInf) {
                double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = ay + tY / div; jY = 1.0 / (div * div);
            } else if (ayInf && !byInf) {
                double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0;
                y = by - tY / div; jY = 1.0 / (div * div);
            }

            return seen.record(f.apply(x, y)) * jX * jY;
        };

        // Delegate to the parallel 2D meta-integrator
        double result = MetaIntegrator.integrate2DSmart(setup, transformed, transAx, transBx, transAy, transBy,
                epsTol);

        boolean anyInfinite = axInf || bxInf || ayInf || byInf;
        if (!verify || !anyInfinite) {
            return result;
        }

        Peak probed = probe2(f, ax, bx, ay, by, cx, cy);
        if (probed.value <= 0.0) {
            return result;
        }
        if (seen.value() < MISS_FACTOR * probed.value) {
            throw missed("[" + limit(ax) + ", " + limit(bx) + "] x [" + limit(ay) + ", " + limit(by) + "]",
                    seen.value(), probed.value, point(probed.x, probed.y));
        }
        if (!farFromAnchor(probed.x, ax, bx, cx) && !farFromAnchor(probed.y, ay, by, cy)) {
            return result;
        }

        // The comparison has to be the one this check always made, so the
        // orthants are integrated uncentered: a center of the whole domain need
        // not even lie inside the piece being integrated here.
        double split = 0.0;
        for (int quadrant = 0; quadrant < 4; ++quadrant) {
            double[] rx = orthant(probed.x, ax, bx, (quadrant & 1) == 0);
            double[] ry = orthant(probed.y, ay, by, (quadrant & 2) == 0);
            split += integrate2D(setup, f, rx[0], rx[1], ry[0], ry[1], epsTol, false, Double.NaN, Double.NaN);
        }
        if (disagree(split, result)) {
            throw disagreement(point(probed.x, probed.y), result, split);
        }
        return result;
    }

    // =========================================================================
    // MULTIDIMENSIONAL INFINITE INTEGRATOR (3D)
    // =========================================================================
    /**
     * Integrates {@code f} over a box whose sides may be infinite.
     * <p>
     * The result is checked against a probe of the integrand and, when the mass
     * sits away from the origin, against the same integral split into octants
     * around it; see the class comment.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand, or sampled it but disagrees with
     *             the same integral split around it
     */
    public static double integrate3DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f,
                                               double ax, double bx, double ay, double by, double az, double bz, double epsTol) {
        return integrate3D(setup, f, ax, bx, ay, by, az, bz, epsTol, true, Double.NaN, Double.NaN, Double.NaN);
    }

    /**
     * Integrates {@code f} over a box whose sides may be infinite, with each
     * infinite axis of the substitution centered where the caller says the mass
     * lies instead of on the origin or on the finite end.
     * <p>
     * See the single-argument sibling
     * {@link #integrate1DInfinite(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double, double)}
     * for what a center is and why supplying one does not disarm the guard.
     * Centering is never inferred here: an axis whose center is
     * {@code Double.NaN} is substituted exactly as it was before.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
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
     *             if the substitution demonstrably never sampled the region
     *             that carries the integrand, or sampled it but disagrees with
     *             the same integral split around it
     * @since 1.5.2
     */
    public static double integrate3DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f,
                                               double ax, double bx, double ay, double by, double az, double bz,
                                               double epsTol, double centerX, double centerY, double centerZ) {
        return integrate3D(setup, f, ax, bx, ay, by, az, bz, epsTol, true, centerX, centerY, centerZ);
    }

    private static double integrate3D(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f,
                                        double ax, double bx, double ay, double by, double az, double bz,
                                        double epsTol, boolean verify, double rawCx, double rawCy, double rawCz) {

        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);
        boolean azInf = (az == Double.NEGATIVE_INFINITY); boolean bzInf = (bz == Double.POSITIVE_INFINITY);

        final double cx = checkCenter(rawCx, ax, bx, "centerX");
        final double cy = checkCenter(rawCy, ay, by, "centerY");
        final double cz = checkCenter(rawCz, az, bz, "centerZ");
        final boolean onX = !Double.isNaN(cx);
        final boolean onY = !Double.isNaN(cy);
        final boolean onZ = !Double.isNaN(cz);

        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;
        double transAz = (azInf && bzInf) ? -1.0 : 0.0; double transBz = (azInf && bzInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }
        if (!azInf && !bzInf) { transAz = az; transBz = bz; }

        // See the 2D form: a centered axis is substituted about its center
        if (onX) { transAx = axInf ? -1.0 : preimage(ax - cx); transBx = bxInf ? 1.0 : preimage(bx - cx); }
        if (onY) { transAy = ayInf ? -1.0 : preimage(ay - cy); transBy = byInf ? 1.0 : preimage(by - cy); }
        if (onZ) { transAz = azInf ? -1.0 : preimage(az - cz); transBz = bzInf ? 1.0 : preimage(bz - cz); }

        final Watch seen = new Watch();

        DTriFunction transformed = (tX, tY, tZ) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (onX) { double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0; x = cx + tX / div; jX = (1.0 + tX * tX) / (div * div); }
            else if (axInf && bxInf) { double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0; x = tX / div; jX = (1.0 + tX * tX) / (div * div); }
            else if (!axInf && bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = ax + tX / div; jX = 1.0 / (div * div); }
            else if (axInf && !bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = bx - tX / div; jX = 1.0 / (div * div); }

            // Transform Y-axis
            double y = tY, jY = 1.0;
            if (onY) { double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0; y = cy + tY / div; jY = (1.0 + tY * tY) / (div * div); }
            else if (ayInf && byInf) { double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0; y = tY / div; jY = (1.0 + tY * tY) / (div * div); }
            else if (!ayInf && byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = ay + tY / div; jY = 1.0 / (div * div); }
            else if (ayInf && !byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = by - tY / div; jY = 1.0 / (div * div); }

            // Transform Z-axis
            double z = tZ, jZ = 1.0;
            if (onZ) { double div = 1.0 - tZ * tZ; if (Math.abs(div) < 1e-15) return 0.0; z = cz + tZ / div; jZ = (1.0 + tZ * tZ) / (div * div); }
            else if (azInf && bzInf) { double div = 1.0 - tZ * tZ; if (Math.abs(div) < 1e-15) return 0.0; z = tZ / div; jZ = (1.0 + tZ * tZ) / (div * div); }
            else if (!azInf && bzInf) { double div = 1.0 - tZ; if (Math.abs(div) < 1e-15) return 0.0; z = az + tZ / div; jZ = 1.0 / (div * div); }
            else if (azInf && !bzInf) { double div = 1.0 - tZ; if (Math.abs(div) < 1e-15) return 0.0; z = bz - tZ / div; jZ = 1.0 / (div * div); }

            return seen.record(f.apply(x, y, z)) * jX * jY * jZ;
        };

        double result = MetaIntegrator.integrate3DSmart(setup, transformed, transAx, transBx, transAy, transBy,
                transAz, transBz, epsTol);

        boolean anyInfinite = axInf || bxInf || ayInf || byInf || azInf || bzInf;
        if (!verify || !anyInfinite) {
            return result;
        }

        Peak probed = probe3(f, ax, bx, ay, by, az, bz, cx, cy, cz);
        if (probed.value <= 0.0) {
            return result;
        }
        if (seen.value() < MISS_FACTOR * probed.value) {
            throw missed("[" + limit(ax) + ", " + limit(bx) + "] x [" + limit(ay) + ", " + limit(by) + "] x ["
                    + limit(az) + ", " + limit(bz) + "]", seen.value(), probed.value,
                    point(probed.x, probed.y, probed.z));
        }
        if (!farFromAnchor(probed.x, ax, bx, cx) && !farFromAnchor(probed.y, ay, by, cy)
                && !farFromAnchor(probed.z, az, bz, cz)) {
            return result;
        }

        // Uncentered for the same reason as in 2D: the comparison has to stay
        // the one this check always made.
        double split = 0.0;
        for (int octant = 0; octant < 8; ++octant) {
            double[] rx = orthant(probed.x, ax, bx, (octant & 1) == 0);
            double[] ry = orthant(probed.y, ay, by, (octant & 2) == 0);
            double[] rz = orthant(probed.z, az, bz, (octant & 4) == 0);
            split += integrate3D(setup, f, rx[0], rx[1], ry[0], ry[1], rz[0], rz[1], epsTol, false, Double.NaN,
                    Double.NaN, Double.NaN);
        }
        if (disagree(split, result)) {
            throw disagreement(point(probed.x, probed.y, probed.z), result, split);
        }
        return result;
    }

    // =========================================================================
    // PROBES
    // =========================================================================

    /**
     * Probes the integrand along a geometric ladder of magnitudes running
     * outwards from the finite end of the domain, or from the origin when there
     * is none, then refines around the best rung.
     */
    private static Peak probe(DFunction f, double a, double b, double center) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        double anchor = (aInf && bInf) ? 0.0 : (aInf ? b : a);

        Peak peak = new Peak();
        peak.record(anchor, f.apply(anchor));
        // A center is added to what the ladder samples, never substituted for
        // it: one evaluation buys the caller's point without giving up the
        // guard against a caller who named the wrong one.
        if (!Double.isNaN(center)) {
            peak.record(center, f.apply(center));
        }
        for (int k = 0; k <= LADDER_STEPS; ++k) {
            double magnitude = Math.pow(LADDER_RATIO, k - LADDER_UNIT);
            if (aInf && bInf) {
                peak.record(magnitude, f.apply(magnitude));
                peak.record(-magnitude, f.apply(-magnitude));
            } else if (bInf) {
                double x = a + magnitude;
                peak.record(x, f.apply(x));
            } else {
                double x = b - magnitude;
                peak.record(x, f.apply(x));
            }
        }

        double span = Math.abs(peak.x) * (LADDER_RATIO - 1.0 / LADDER_RATIO);
        for (int round = 0; round < REFINEMENTS && span > 0.0; ++round) {
            double centre = peak.x;
            for (int i = -REFINEMENT_SAMPLES; i <= REFINEMENT_SAMPLES; ++i) {
                double x = centre + (span * i) / REFINEMENT_SAMPLES;
                if (x < a || x > b) {
                    continue;
                }
                peak.record(x, f.apply(x));
            }
            span /= REFINEMENT_NARROWING;
        }
        return peak;
    }

    /**
     * The sample positions of one axis: a geometric ladder outwards from
     * whichever end is finite, both signs when neither is, and a uniform scan
     * when the axis is bounded on both sides. A center is appended to that
     * list, which costs one position per axis rather than a second grid.
     */
    private static double[] axisSamples(double a, double b, double center, double ratio, int rungs) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        int n = 2 * rungs + 1;
        double[] samples;
        if (!aInf && !bInf) {
            samples = new double[n];
            for (int i = 0; i < n; ++i) {
                samples[i] = a + ((b - a) * i) / (n - 1.0);
            }
        } else if (aInf && bInf) {
            samples = new double[2 * n + 1];
            samples[0] = 0.0;
            for (int k = 0; k < n; ++k) {
                double magnitude = Math.pow(ratio, k - rungs);
                samples[1 + 2 * k] = magnitude;
                samples[2 + 2 * k] = -magnitude;
            }
        } else {
            double anchor = aInf ? b : a;
            samples = new double[n + 1];
            samples[0] = anchor;
            for (int k = 0; k < n; ++k) {
                double magnitude = Math.pow(ratio, k - rungs);
                samples[1 + k] = aInf ? (anchor - magnitude) : (anchor + magnitude);
            }
        }
        if (Double.isNaN(center)) {
            return samples;
        }
        double[] withCenter = Arrays.copyOf(samples, samples.length + 1);
        withCenter[samples.length] = center;
        return withCenter;
    }

    private static Peak probe2(DBiFunction f, double ax, double bx, double ay, double by, double cx, double cy) {
        double[] xs = axisSamples(ax, bx, cx, GRID_RATIO_2D, GRID_RUNGS_2D);
        double[] ys = axisSamples(ay, by, cy, GRID_RATIO_2D, GRID_RUNGS_2D);
        Peak peak = new Peak();
        for (int i = 0; i < xs.length; ++i) {
            for (int j = 0; j < ys.length; ++j) {
                peak.record(xs[i], ys[j], f.apply(xs[i], ys[j]));
            }
        }
        for (int round = 0; round < SWEEP_ROUNDS; ++round) {
            for (int axis = 0; axis < 2; ++axis) {
                double span = Math.max(Math.abs(axis == 0 ? peak.x : peak.y), 1.0);
                for (int r = 0; r < SWEEP_REFINEMENTS; ++r) {
                    double centre = (axis == 0) ? peak.x : peak.y;
                    for (int i = -REFINEMENT_SAMPLES; i <= REFINEMENT_SAMPLES; ++i) {
                        double c = centre + (span * i) / REFINEMENT_SAMPLES;
                        double x = (axis == 0) ? c : peak.x;
                        double y = (axis == 1) ? c : peak.y;
                        if (x < ax || x > bx || y < ay || y > by) {
                            continue;
                        }
                        peak.record(x, y, f.apply(x, y));
                    }
                    span /= SWEEP_NARROWING;
                }
            }
        }
        return peak;
    }

    private static Peak probe3(DTriFunction f, double ax, double bx, double ay, double by, double az, double bz,
            double cx, double cy, double cz) {
        double[] xs = axisSamples(ax, bx, cx, GRID_RATIO_3D, GRID_RUNGS_3D);
        double[] ys = axisSamples(ay, by, cy, GRID_RATIO_3D, GRID_RUNGS_3D);
        double[] zs = axisSamples(az, bz, cz, GRID_RATIO_3D, GRID_RUNGS_3D);
        Peak peak = new Peak();
        for (int i = 0; i < xs.length; ++i) {
            for (int j = 0; j < ys.length; ++j) {
                for (int k = 0; k < zs.length; ++k) {
                    peak.record(xs[i], ys[j], zs[k], f.apply(xs[i], ys[j], zs[k]));
                }
            }
        }
        for (int round = 0; round < SWEEP_ROUNDS; ++round) {
            for (int axis = 0; axis < 3; ++axis) {
                double start = (axis == 0) ? peak.x : (axis == 1) ? peak.y : peak.z;
                double span = Math.max(Math.abs(start), 1.0);
                for (int r = 0; r < SWEEP_REFINEMENTS; ++r) {
                    double centre = (axis == 0) ? peak.x : (axis == 1) ? peak.y : peak.z;
                    for (int i = -REFINEMENT_SAMPLES; i <= REFINEMENT_SAMPLES; ++i) {
                        double c = centre + (span * i) / REFINEMENT_SAMPLES;
                        double x = (axis == 0) ? c : peak.x;
                        double y = (axis == 1) ? c : peak.y;
                        double z = (axis == 2) ? c : peak.z;
                        if (x < ax || x > bx || y < ay || y > by || z < az || z > bz) {
                            continue;
                        }
                        peak.record(x, y, z, f.apply(x, y, z));
                    }
                    span /= SWEEP_NARROWING;
                }
            }
        }
        return peak;
    }

    // =========================================================================
    // VERDICTS
    // =========================================================================

    /**
     * Returns the substitution's own result when the probe agrees that it
     * looked in the right place, the double-exponential rule's result when it
     * did not but that rule can supply one, and throws otherwise. A probe that
     * finds nothing at all does not accuse: an integrand negligible everywhere
     * the ladder reaches is entitled to integrate to zero.
     */
    private static double verify1D(DFunction f, Watch seen, double a, double b, double result, double epsTol,
            double center) {
        Peak probed = probe(f, a, b, center);
        if (probed.value <= 0.0) {
            // The ladder found nothing anywhere, and an integrand negligible
            // everywhere it reached is entitled to integrate to zero -- unless
            // the caller named a place to look, which is worth trying before
            // that zero is handed out.
            if (!Double.isNaN(center)) {
                double named = rescue1D(f, a, b, center, epsTol);
                if (!Double.isNaN(named)) {
                    return named;
                }
            }
            return result;
        }
        if (seen.value() >= MISS_FACTOR * probed.value) {
            // The substitution did sample the mass, which used to settle it:
            // without a center the mass can only be resolved or missed, never
            // sampled badly, because the substitution sits on the anchor. A
            // center the caller chose can be far from where the mass turned out
            // to be, and then sampling it is not the same as resolving it --
            // the failure the 2D and 3D forms answer with their split check. In
            // one dimension there is a better answer than a refusal: split
            // where the mass is and let the rule that crowds its nodes there
            // settle the disagreement.
            if (!Double.isNaN(center) && farFromAnchor(probed.x, a, b, center)) {
                double split = rescue1D(f, a, b, probed.x, epsTol);
                if (!Double.isNaN(split) && disagree(split, result)) {
                    return split;
                }
            }
            return result;
        }
        // The substitution looked in the wrong place, and the message below
        // would tell the caller to split the range where the mass is. The
        // probe knows where that is, so do it here instead: split there and
        // hand both halves to DoubleExponential, whose nodes crowd towards the
        // finite end of each - which is now the peak itself.
        double rescued = rescue1D(f, a, b, probed.x, epsTol);
        if (!Double.isNaN(rescued)) {
            return rescued;
        }
        throw missed("[" + limit(a) + ", " + limit(b) + "]", seen.value(), probed.value, point(probed.x));
    }

    /**
     * The integral over {@code [a, cut]} plus the one over {@code [cut, b]}, or
     * {@code NaN} when either half declines to converge - which is what happens
     * when {@code cut} is not in fact where the mass is, so the two flags
     * together are the check on the probe as well as on the rule.
     */
    private static double rescue1D(DFunction f, double a, double b, double cut, double epsTol) {
        if (!(cut > a) || !(cut < b)) {
            return Double.NaN;
        }
        DoubleExponential.IntegralResult lower = DoubleExponential.integrate1D(f, a, cut, epsTol);
        if (!lower.converged) {
            return Double.NaN;
        }
        DoubleExponential.IntegralResult upper = DoubleExponential.integrate1D(f, cut, b, epsTol);
        if (!upper.converged) {
            return Double.NaN;
        }
        return lower.value + upper.value;
    }

    /** Is this coordinate far from the point its substitution is centered on? */
    private static boolean farFromAnchor(double coordinate, double a, double b, double center) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        if (!aInf && !bInf) {
            return false;
        }
        double anchor = !Double.isNaN(center) ? center : ((aInf && bInf) ? 0.0 : (aInf ? b : a));
        return Math.abs(coordinate - anchor) > SAFE_RADIUS;
    }

    /**
     * Validates one axis's center and returns it, {@code NaN} standing for
     * none. An axis with two finite ends is not substituted at all, so a center
     * given for one is dropped here rather than silently changing the rule.
     */
    private static double checkCenter(double center, double a, double b, String name) {
        if (Double.isNaN(center)) {
            return Double.NaN;
        }
        if (Double.isInfinite(center)) {
            throw new IllegalArgumentException(name + " must be finite, was " + center);
        }
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        if (!aInf && !bInf) {
            return Double.NaN;
        }
        if (!aInf && !(center > a)) {
            throw new IllegalArgumentException(String.format(Locale.ROOT,
                    "%s must lie strictly above the finite end %.6g of [%s, %s], was %.6g", name, Double.valueOf(a),
                    limit(a), limit(b), Double.valueOf(center)));
        }
        if (!bInf && !(center < b)) {
            throw new IllegalArgumentException(String.format(Locale.ROOT,
                    "%s must lie strictly below the finite end %.6g of [%s, %s], was %.6g", name, Double.valueOf(b),
                    limit(a), limit(b), Double.valueOf(center)));
        }
        return center;
    }

    /** One side of a split of {@code [a, b]} at {@code cut}. */
    private static double[] orthant(double cut, double a, double b, boolean upper) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        if (!aInf && !bInf) {
            return new double[] { a, b };
        }
        return upper ? new double[] { cut, b } : new double[] { a, cut };
    }

    private static boolean disagree(double split, double result) {
        double scale = Math.max(Math.abs(split), Math.abs(result));
        return Math.abs(split - result) > 1.0e-6 * scale;
    }

    private static ArithmeticException missed(String domain, double seen, double probed, String where) {
        return new ArithmeticException(String.format(Locale.ROOT,
                "the substitution for %s never sampled the integrand where its mass lies: the largest |f| it "
                        + "evaluated was %.3e, while |f%s| is %.3e. Integrate over a finite domain containing "
                        + "%s instead, or split the range there, which resolves it.",
                domain, Double.valueOf(seen), where, Double.valueOf(probed), where));
    }

    private static ArithmeticException disagreement(String where, double result, double split) {
        return new ArithmeticException(String.format(Locale.ROOT,
                "the substitution sampled the integrand but did not resolve it: over the whole domain it gives "
                        + "%.10g, while splitting at %s -- where the mass is -- gives %.10g. Split the domain "
                        + "there yourself, or integrate over a finite one around it.",
                Double.valueOf(result), where, Double.valueOf(split)));
    }

    private static String point(double x) {
        return "(" + coordinate(x) + ")";
    }

    private static String point(double x, double y) {
        return "(" + coordinate(x) + ", " + coordinate(y) + ")";
    }

    private static String point(double x, double y, double z) {
        return "(" + coordinate(x) + ", " + coordinate(y) + ", " + coordinate(z) + ")";
    }

    private static String coordinate(double x) {
        return String.format(Locale.ROOT, "%.6g", Double.valueOf(x));
    }

    /** Renders an integration limit for a refusal message. */
    private static String limit(double x) {
        if (x == Double.NEGATIVE_INFINITY) {
            return "-inf";
        }
        if (x == Double.POSITIVE_INFINITY) {
            return "+inf";
        }
        return coordinate(x);
    }
}
