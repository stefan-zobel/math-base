package math.solve;

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
 * substitution produces, and what the guard below now refuses to hand out. The
 * one-sided forms substitute {@code x = a + t / (1 - t)} and fail the same way
 * when the mass lies far from the finite end.
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
 * place, and an {@link ArithmeticException} names where the mass is and what to
 * do instead.
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
 *   <li><b>1D</b>: 7 of 7 known failures refused, 0 of 9 legitimate integrals
 *       disturbed, about 1,200 evaluations of the integrand per call.</li>
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
 * once. The real cure is to center the substitution on the integrand, which is
 * a change to the method rather than a guard on it.
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
     * of the integrand and refused if the substitution never sampled the region
     * that carries the mass; see the class comment.
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
     *             that carries the integrand
     */
    public static double integrate1DInfinite(AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                               double a, double b, double epsTol) {

        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);

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
            refuseIfMissed1D(f, seen, a, b);
            return result;
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
            refuseIfMissed1D(f, seen, a, b);
            return result;
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
            refuseIfMissed1D(f, seen, a, b);
            return result;
        }

        // CASE D: Ordinary finite integral [a, b]
        return MetaIntegrator.integrate1DSmart(setup, f, a, b, epsTol);
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
        return integrate2D(setup, f, ax, bx, ay, by, epsTol, true);
    }

    private static double integrate2D(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
                                        double ax, double bx, double ay, double by, double epsTol,
                                        boolean verify) {

        // 1. Analyze bounds
        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);

        // 2. Set target integration bounds for the transformed variables
        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }

        final Watch seen = new Watch();

        // 3. Build dynamic wrapper function
        DBiFunction transformed = (tX, tY) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (axInf && bxInf) {
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
            if (ayInf && byInf) {
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

        Peak probed = probe2(f, ax, bx, ay, by);
        if (probed.value <= 0.0) {
            return result;
        }
        if (seen.value() < MISS_FACTOR * probed.value) {
            throw missed("[" + limit(ax) + ", " + limit(bx) + "] x [" + limit(ay) + ", " + limit(by) + "]",
                    seen.value(), probed.value, point(probed.x, probed.y));
        }
        if (!farFromAnchor(probed.x, ax, bx) && !farFromAnchor(probed.y, ay, by)) {
            return result;
        }

        double split = 0.0;
        for (int quadrant = 0; quadrant < 4; ++quadrant) {
            double[] rx = orthant(probed.x, ax, bx, (quadrant & 1) == 0);
            double[] ry = orthant(probed.y, ay, by, (quadrant & 2) == 0);
            split += integrate2D(setup, f, rx[0], rx[1], ry[0], ry[1], epsTol, false);
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
        return integrate3D(setup, f, ax, bx, ay, by, az, bz, epsTol, true);
    }

    private static double integrate3D(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f,
                                        double ax, double bx, double ay, double by, double az, double bz,
                                        double epsTol, boolean verify) {

        boolean axInf = (ax == Double.NEGATIVE_INFINITY); boolean bxInf = (bx == Double.POSITIVE_INFINITY);
        boolean ayInf = (ay == Double.NEGATIVE_INFINITY); boolean byInf = (by == Double.POSITIVE_INFINITY);
        boolean azInf = (az == Double.NEGATIVE_INFINITY); boolean bzInf = (bz == Double.POSITIVE_INFINITY);

        double transAx = (axInf && bxInf) ? -1.0 : 0.0; double transBx = (axInf && bxInf) ? 1.0 : 1.0;
        double transAy = (ayInf && byInf) ? -1.0 : 0.0; double transBy = (ayInf && byInf) ? 1.0 : 1.0;
        double transAz = (azInf && bzInf) ? -1.0 : 0.0; double transBz = (azInf && bzInf) ? 1.0 : 1.0;

        if (!axInf && !bxInf) { transAx = ax; transBx = bx; }
        if (!ayInf && !byInf) { transAy = ay; transBy = by; }
        if (!azInf && !bzInf) { transAz = az; transBz = bz; }

        final Watch seen = new Watch();

        DTriFunction transformed = (tX, tY, tZ) -> {
            // Transform X-axis
            double x = tX, jX = 1.0;
            if (axInf && bxInf) { double div = 1.0 - tX * tX; if (Math.abs(div) < 1e-15) return 0.0; x = tX / div; jX = (1.0 + tX * tX) / (div * div); }
            else if (!axInf && bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = ax + tX / div; jX = 1.0 / (div * div); }
            else if (axInf && !bxInf) { double div = 1.0 - tX; if (Math.abs(div) < 1e-15) return 0.0; x = bx - tX / div; jX = 1.0 / (div * div); }

            // Transform Y-axis
            double y = tY, jY = 1.0;
            if (ayInf && byInf) { double div = 1.0 - tY * tY; if (Math.abs(div) < 1e-15) return 0.0; y = tY / div; jY = (1.0 + tY * tY) / (div * div); }
            else if (!ayInf && byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = ay + tY / div; jY = 1.0 / (div * div); }
            else if (ayInf && !byInf) { double div = 1.0 - tY; if (Math.abs(div) < 1e-15) return 0.0; y = by - tY / div; jY = 1.0 / (div * div); }

            // Transform Z-axis
            double z = tZ, jZ = 1.0;
            if (azInf && bzInf) { double div = 1.0 - tZ * tZ; if (Math.abs(div) < 1e-15) return 0.0; z = tZ / div; jZ = (1.0 + tZ * tZ) / (div * div); }
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

        Peak probed = probe3(f, ax, bx, ay, by, az, bz);
        if (probed.value <= 0.0) {
            return result;
        }
        if (seen.value() < MISS_FACTOR * probed.value) {
            throw missed("[" + limit(ax) + ", " + limit(bx) + "] x [" + limit(ay) + ", " + limit(by) + "] x ["
                    + limit(az) + ", " + limit(bz) + "]", seen.value(), probed.value,
                    point(probed.x, probed.y, probed.z));
        }
        if (!farFromAnchor(probed.x, ax, bx) && !farFromAnchor(probed.y, ay, by)
                && !farFromAnchor(probed.z, az, bz)) {
            return result;
        }

        double split = 0.0;
        for (int octant = 0; octant < 8; ++octant) {
            double[] rx = orthant(probed.x, ax, bx, (octant & 1) == 0);
            double[] ry = orthant(probed.y, ay, by, (octant & 2) == 0);
            double[] rz = orthant(probed.z, az, bz, (octant & 4) == 0);
            split += integrate3D(setup, f, rx[0], rx[1], ry[0], ry[1], rz[0], rz[1], epsTol, false);
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
    private static Peak probe(DFunction f, double a, double b) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        double anchor = (aInf && bInf) ? 0.0 : (aInf ? b : a);

        Peak peak = new Peak();
        peak.record(anchor, f.apply(anchor));
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
     * when the axis is bounded on both sides.
     */
    private static double[] axisSamples(double a, double b, double ratio, int rungs) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        int n = 2 * rungs + 1;
        if (!aInf && !bInf) {
            double[] samples = new double[n];
            for (int i = 0; i < n; ++i) {
                samples[i] = a + ((b - a) * i) / (n - 1.0);
            }
            return samples;
        }
        if (aInf && bInf) {
            double[] samples = new double[2 * n + 1];
            samples[0] = 0.0;
            for (int k = 0; k < n; ++k) {
                double magnitude = Math.pow(ratio, k - rungs);
                samples[1 + 2 * k] = magnitude;
                samples[2 + 2 * k] = -magnitude;
            }
            return samples;
        }
        double anchor = aInf ? b : a;
        double[] samples = new double[n + 1];
        samples[0] = anchor;
        for (int k = 0; k < n; ++k) {
            double magnitude = Math.pow(ratio, k - rungs);
            samples[1 + k] = aInf ? (anchor - magnitude) : (anchor + magnitude);
        }
        return samples;
    }

    private static Peak probe2(DBiFunction f, double ax, double bx, double ay, double by) {
        double[] xs = axisSamples(ax, bx, GRID_RATIO_2D, GRID_RUNGS_2D);
        double[] ys = axisSamples(ay, by, GRID_RATIO_2D, GRID_RUNGS_2D);
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

    private static Peak probe3(DTriFunction f, double ax, double bx, double ay, double by, double az, double bz) {
        double[] xs = axisSamples(ax, bx, GRID_RATIO_3D, GRID_RUNGS_3D);
        double[] ys = axisSamples(ay, by, GRID_RATIO_3D, GRID_RUNGS_3D);
        double[] zs = axisSamples(az, bz, GRID_RATIO_3D, GRID_RUNGS_3D);
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
     * Throws when the substitution never came near the part of the domain that
     * carries the integrand. A probe that finds nothing at all does not accuse:
     * an integrand negligible everywhere the ladder reaches is entitled to
     * integrate to zero.
     */
    private static void refuseIfMissed1D(DFunction f, Watch seen, double a, double b) {
        Peak probed = probe(f, a, b);
        if (probed.value <= 0.0 || seen.value() >= MISS_FACTOR * probed.value) {
            return;
        }
        throw missed("[" + limit(a) + ", " + limit(b) + "]", seen.value(), probed.value, point(probed.x));
    }

    /** Is this coordinate far from the point its substitution is centered on? */
    private static boolean farFromAnchor(double coordinate, double a, double b) {
        boolean aInf = (a == Double.NEGATIVE_INFINITY);
        boolean bInf = (b == Double.POSITIVE_INFINITY);
        if (!aInf && !bInf) {
            return false;
        }
        double anchor = (aInf && bInf) ? 0.0 : (aInf ? b : a);
        return Math.abs(coordinate - anchor) > SAFE_RADIUS;
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
