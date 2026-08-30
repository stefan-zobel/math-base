package math.solve;

import java.util.logging.Level;
import java.util.logging.Logger;

import math.fun.DFunction;
import math.fun.DBiFunction;
import math.fun.DTriFunction;

/**
 * Picks an integration strategy per call: a cheap oscillation heuristic decides
 * between the refining Clenshaw-Curtis rule of {@link ClenshawCurtis} for an
 * oscillating integrand and the adaptive Gauss-Kronrod subdivision of
 * {@link AdaptiveGaussKronrod} for everything else.
 * <p>
 * Clenshaw-Curtis is asked for the tolerance and answers whether it reached it.
 * If it did not - its grid tops out at 2049 points in 1D, 257 per axis in 2D
 * and 65 per axis in 3D - the result is discarded and the adaptive
 * Gauss-Kronrod path runs instead. That fallback is expensive, in 3D by more
 * than an order of magnitude, but a slow right answer beats a fast wrong one,
 * and it is reached only by an integrand that oscillates faster than a few
 * hundred half waves per axis.
 * <p>
 * The heuristic that routes into the Clenshaw-Curtis branch rests on the error
 * estimate of the <em>undivided</em> domain and is therefore no more reliable
 * than that estimate; see the class comment of {@link AdaptiveGaussKronrod}.
 * Both branches are safe against a wrong routing decision, so this affects the
 * cost of a call rather than its result.
 * <p>
 * Nothing routes to {@link DoubleExponential} from here, and nothing is meant
 * to: a heuristic for "has an endpoint singularity" would be guesswork where
 * this one is merely unreliable. {@link Quadrature} offers that rule by name
 * instead, so the choice is the caller's and explicit.
 */
public class MetaIntegrator {

    private static final Logger LOG = Logger.getLogger(MetaIntegrator.class.getName());

    // ==========================================
    // SMART INTEGRATOR 1D
    // ==========================================
    public static double integrate1DSmart(AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                           double a, double b, double epsTol) {
        return integrate1DSmartWithError(setup, f, a, b, epsTol).value;
    }

    /**
     * The same integral as
     * {@link #integrate1DSmart(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double)},
     * with the error estimate of whichever rule answered and whether that rule
     * reached the tolerance. The two probe integrations that route the call
     * belong to the heuristic and contribute nothing to the reported error.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param a
     *            lower limit
     * @param b
     *            upper limit
     * @param epsTol
     *            error tolerance
     * @return the approximated integral, its estimate and whether the rule that
     *         supplied it reached {@code epsTol}
     * @since 1.5.3
     */
    public static AdaptiveGaussKronrod.IntegralResult integrate1DSmartWithError(
                                           AdaptiveGaussKronrod.G7_K15 setup, DFunction f,
                                           double a, double b, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate1D(setup, f, a, b);

        // There is deliberately no early return on the value of firstStep here.
        // On the undivided domain the error estimate is worthless whenever the
        // integrand has a feature narrower than the node spacing - see the
        // class comment of AdaptiveGaussKronrod. firstStep is used for the
        // oscillation heuristic only; the tolerance test happens in
        // integrate1DAdaptive, after its forced subdivisions.
        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate1D(setup, x -> Math.abs(f.apply(x)), a, b);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Aliasing triggers an oscillation switch ONLY if the values genuinely cancel each other out
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value)) && (oscillationIndex < 0.1);

        // case 2: Genuine Oscillation (Low oscillation index OR heavy aliasing noise)
        if ((oscillationIndex < 0.05 || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) { // case 2: Oscillation
            if (LOG.isLoggable(Level.FINE)) {
                LOG.fine("[Meta1D] Oscillation detected (Index: " + oscillationIndex + ", Noise: " + isAliasedOscillation + ") -> Switch to Clenshaw-Curtis.");
            }
            ClenshawCurtis.IntegralResult cc = ClenshawCurtis.integrate1D(f, a, b, epsTol);
            if (cc.converged) {
                return new AdaptiveGaussKronrod.IntegralResult(cc.value, cc.approximatedErrorEstimate, true);
            }
            // the grid could not resolve the oscillation and says so; fall
            // through to the subdivision rather than return what it disowns
        }

        // case 3: Subdivision
        return AdaptiveGaussKronrod.integrate1DAdaptiveWithError(setup, f, a, b, epsTol, 20);
    }

    // ==========================================
    // SMART INTEGRATOR 2D
    // ==========================================
    public static double integrate2DSmart(AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
            double ax, double bx, double ay, double by, double epsTol) {
        return integrate2DSmartWithError(setup, f, ax, bx, ay, by, epsTol).value;
    }

    /**
     * The same integral as
     * {@link #integrate2DSmart(AdaptiveGaussKronrod.G7_K15, DBiFunction, double, double, double, double, double)},
     * with the error estimate of whichever rule answered and whether that rule
     * reached the tolerance. See
     * {@link #integrate1DSmartWithError(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double)}.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x
     * @param bx
     *            upper limit in x
     * @param ay
     *            lower limit in y
     * @param by
     *            upper limit in y
     * @param epsTol
     *            error tolerance
     * @return the approximated integral, its estimate and whether the rule that
     *         supplied it reached {@code epsTol}
     * @since 1.5.3
     */
    public static AdaptiveGaussKronrod.IntegralResult integrate2DSmartWithError(
            AdaptiveGaussKronrod.G7_K15 setup, DBiFunction f,
            double ax, double bx, double ay, double by, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate2DParallel(setup, f, ax, bx, ay, by);

        // no early return on firstStep, see integrate1DSmart
        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate2DParallel(setup, (x, y) -> Math.abs(f.apply(x, y)), ax, bx, ay, by);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Aliasing criteria: high error combined with clean cancellation signature
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value)) && (oscillationIndex < 0.1);

        if (((oscillationIndex < 0.05) || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) {
            if (LOG.isLoggable(Level.FINE)) {
                LOG.fine("[Meta2D] Oscillation detected -> Route to efficient Clenshaw-Curtis.");
            }
            ClenshawCurtis.IntegralResult cc = ClenshawCurtis.integrate2D(f, ax, bx, ay, by, epsTol);
            if (cc.converged) {
                return new AdaptiveGaussKronrod.IntegralResult(cc.value, cc.approximatedErrorEstimate, true);
            }
            // no convergence, see integrate1DSmart
        }

        return AdaptiveGaussKronrod.integrate2DAdaptiveWithError(setup, f, ax, bx, ay, by, epsTol, 10);
    }

    // ==========================================
    // SMART INTEGRATOR 3D
    // ==========================================

    public static double integrate3DSmart(AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f, double ax, double bx,
            double ay, double by, double az, double bz, double epsTol) {
        return integrate3DSmartWithError(setup, f, ax, bx, ay, by, az, bz, epsTol).value;
    }

    /**
     * The same integral as
     * {@link #integrate3DSmart(AdaptiveGaussKronrod.G7_K15, DTriFunction, double, double, double, double, double, double, double)},
     * with the error estimate of whichever rule answered and whether that rule
     * reached the tolerance. See
     * {@link #integrate1DSmartWithError(AdaptiveGaussKronrod.G7_K15, DFunction, double, double, double)}.
     *
     * @param setup
     *            the Gauss-Kronrod rule to use
     * @param f
     *            the integrand
     * @param ax
     *            lower limit in x
     * @param bx
     *            upper limit in x
     * @param ay
     *            lower limit in y
     * @param by
     *            upper limit in y
     * @param az
     *            lower limit in z
     * @param bz
     *            upper limit in z
     * @param epsTol
     *            error tolerance
     * @return the approximated integral, its estimate and whether the rule that
     *         supplied it reached {@code epsTol}
     * @since 1.5.3
     */
    public static AdaptiveGaussKronrod.IntegralResult integrate3DSmartWithError(
            AdaptiveGaussKronrod.G7_K15 setup, DTriFunction f, double ax, double bx,
            double ay, double by, double az, double bz, double epsTol) {
        AdaptiveGaussKronrod.IntegralResult firstStep = AdaptiveGaussKronrod.integrate3DParallel(setup, f, ax, bx, ay,
                by, az, bz);

        // no early return on firstStep, see integrate1DSmart
        AdaptiveGaussKronrod.IntegralResult absStep = AdaptiveGaussKronrod.integrate3DParallel(setup,
                (x, y, z) -> Math.abs(f.apply(x, y, z)), ax, bx, ay, by, az, bz);
        double oscillationIndex = Math.abs(firstStep.value) / (absStep.value + 1e-15);

        // Block aliasing triggers on sharp 3D spikes by enforcing the cancellation index check (< 0.1)
        boolean isAliasedOscillation = (firstStep.approximatedErrorEstimate > Math.abs(firstStep.value))
                && (oscillationIndex < 0.1);

        if (((oscillationIndex < 0.05) || isAliasedOscillation) && firstStep.approximatedErrorEstimate > epsTol) {
            if (LOG.isLoggable(Level.FINE)) {
                LOG.fine("[Meta3D] Oscillation detected -> Route to efficient Clenshaw-Curtis.");
            }
            ClenshawCurtis.IntegralResult cc = ClenshawCurtis.integrate3D(f, ax, bx, ay, by, az, bz, epsTol);
            if (cc.converged) {
                return new AdaptiveGaussKronrod.IntegralResult(cc.value, cc.approximatedErrorEstimate, true);
            }
            // no convergence, see integrate1DSmart. The subdivision below costs
            // more than an order of magnitude over the Clenshaw-Curtis ladder
            // in 3D, which is the price of not returning an aliased value
        }

        // Two bisections per axis and a budget that leaves room for the
        // adaptive part on top of them. The 8 that used to stand here was too
        // small to resolve a narrow peak even once the premature accept above
        // was gone: 22 of 31 peak positions came out wrong by 100 percent,
        // and raising only the budget did not help, because three forced
        // levels do not sample finely enough to notice the peak at all.
        return AdaptiveGaussKronrod.integrate3DAdaptiveWithError(setup, f, ax, bx, ay, by, az, bz, epsTol, 14, 6);
    }
}
