package math.stats.bayes;

import math.MathConsts;
import math.fun.DFunction;
import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.optim.NelderMead;
import math.solve.Quadrature;

/**
 * The posterior of two parameters, normalized by quadrature.
 * <p>
 * The two-parameter counterpart of {@link ScalarPosterior}: a caller hands over
 * the <b>logarithm</b> of an unnormalized posterior as a {@link DMultiFunction}
 * together with the box the parameters live in, and gets back the evidence, the
 * mode, the density, both moments and the marginal densities. Nothing is
 * sampled and there is no chain to diagnose.
 * <p>
 * <b>Why this is not just {@code ScalarPosterior} with one more axis.</b> The
 * substitution that carries an infinite domain onto a finite one is applied per
 * axis, so it can be centered but not <em>rotated</em> -- and a two-parameter
 * posterior is almost always a diagonal ridge, since the slope and the
 * intercept of a fit are correlated by construction. Integrated on the caller's
 * axes, a correlated normal pair whose true integral is one comes back as
 * {@code 1.0015} at a correlation of {@code 0.999}, as {@code 1.08} at
 * {@code 0.9999}, and as <b>{@code 2.9e-198}</b> when the two parameters also
 * sit at scales {@code 100} and {@code 0.001} apart -- with no exception and no
 * warning, because the substitution does sample the mass and only fails to
 * resolve it.
 * <p>
 * So the parameters are <b>whitened</b> first: the curvature {@code H} of the
 * negative log posterior at the mode is factored as {@code H = L L^T}, and the
 * integration runs over {@code u = L^T (x - mode)}, where the posterior is
 * approximately a standard normal pair. All four cases above then come back to
 * within {@code 1.5e-13}.
 * <p>
 * <b>The whitening is the Laplace approximation</b>, so
 * {@link #logEvidenceLaplace()} costs nothing extra and the gap between it and
 * {@link #logEvidence()} says how far the posterior is from Gaussian.
 * <p>
 * <b>What construction costs.</b> Seven integrations -- the normalizer, a
 * second one with the axes exchanged as a check, two for the mean and three for
 * the covariance -- at roughly 120 000 evaluations of the caller's function
 * each. An instance is expensive to build, cheap to ask, immutable, and
 * shareable between threads; the integrand is evaluated in parallel, so the
 * function handed in must be thread safe.
 * <p>
 * <b>What it assumes.</b> That the posterior is proper, that it has one peak,
 * and that the peak is a genuine maximum -- a saddle gives a curvature that is
 * not positive definite and is refused rather than whitened. Multimodality is
 * not detected; {@link #errorEstimate()} is the check, since its second route
 * subdivides the same integral differently.
 *
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Bayesian_inference">Wikipedia Bayesian
 *      inference</a>
 * @since 1.5.3
 */
public final class BivariatePosterior {

    /** Error tolerance handed to {@link Quadrature} if none is given. */
    private static final double DEFAULT_EPS_TOL = 1.0e-10;

    /**
     * The fourth root of the machine epsilon, about {@code 1.2e-4}. The step of
     * a central second difference: smaller and the cancellation of the three
     * function values dominates, larger and the difference stops being local.
     */
    private static final double EPS_FOURTH = Math.pow(MathConsts.MACH_EPS_DBL, 0.25);

    /** How far the search for an axis scale may double before giving up. */
    private static final int MAX_SCALE_STEPS = 200;

    private final DMultiFunction logPosterior;
    private final double[] lo;
    private final double[] hi;
    private final double epsTol;

    private final double[] mode;
    /** {@code logPosterior(mode)}, the shift that keeps the integrand alive. */
    private final double logMax;
    /** The curvature at the mode, in the caller's axis order. */
    private final DMatrix curvature;

    /**
     * Which of the caller's axes is the first and which the second of the
     * whitened frame. A bounded axis is put second, for the reason
     * {@link #uLimits()} gives.
     */
    private final int first;
    private final int second;

    /** The Cholesky factor of the permuted curvature, lower triangular. */
    private final double l00;
    private final double l10;
    private final double l11;
    private final double detL;
    /** {@code (L^T)^-1}, upper triangular: {@code z = mode + M u}. */
    private final double m00;
    private final double m01;
    private final double m11;

    private final double normalizer;
    private final double logEvidence;
    private final double logEvidenceLaplace;
    private final double[] mean;
    private final DMatrix covariance;
    private final double errorEstimate;

    /**
     * The posterior of {@code logPosterior} over the box {@code [lo, hi]}, with
     * a starting point and an error tolerance derived from the box.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density of
     *            two parameters. It is never evaluated outside the box, it may
     *            return {@link Double#NEGATIVE_INFINITY} where the posterior is
     *            zero, and it is evaluated from several threads at once, so it
     *            must not keep mutable state
     * @param lo
     *            the lower end of each axis, two entries, either possibly
     *            {@link Double#NEGATIVE_INFINITY}. Not modified
     * @param hi
     *            the upper end of each axis, two entries, either possibly
     *            {@link Double#POSITIVE_INFINITY}. Not modified
     * @return the normalized posterior
     * @throws IllegalArgumentException
     *             if {@code logPosterior}, {@code lo} or {@code hi} is
     *             {@code null}, if either does not hold two entries, if
     *             {@code lo[i]} is not below {@code hi[i]}, if any of them is
     *             {@code NaN}, if no peak can be found, or if the peak is not a
     *             maximum
     */
    public static BivariatePosterior of(DMultiFunction logPosterior, double[] lo, double[] hi) {
        checkBox(lo, hi);
        return new BivariatePosterior(logPosterior, lo, hi, defaultStart(lo, hi), DEFAULT_EPS_TOL);
    }

    /**
     * The posterior of {@code logPosterior} over the box {@code [lo, hi]},
     * starting the search for the peak at {@code start}.
     * <p>
     * Naming {@code start} is what rescues a posterior whose mass sits far from
     * anything the short form would guess at.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density of
     *            two parameters
     * @param lo
     *            the lower end of each axis, two entries. Not modified
     * @param hi
     *            the upper end of each axis, two entries. Not modified
     * @param start
     *            two finite coordinates inside the box, near the mass. Not
     *            modified
     * @param epsTol
     *            the error tolerance of each integration, strictly positive
     * @return the normalized posterior
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #of(DMultiFunction, double[], double[])} states, and
     *             additionally if {@code start} is {@code null}, does not hold
     *             two finite entries inside the box, or if {@code epsTol} is
     *             not strictly positive
     */
    public static BivariatePosterior of(DMultiFunction logPosterior, double[] lo, double[] hi, double[] start,
            double epsTol) {
        checkBox(lo, hi);
        return new BivariatePosterior(logPosterior, lo, hi, start, epsTol);
    }

    private BivariatePosterior(DMultiFunction logPosterior, double[] lo, double[] hi, double[] start,
            double epsTol) {
        if (logPosterior == null) {
            throw new IllegalArgumentException("logPosterior must not be null");
        }
        if (start == null || start.length != 2) {
            throw new IllegalArgumentException("start must hold two coordinates");
        }
        for (int i = 0; i < 2; i++) {
            if (!isFinite(start[i]) || start[i] < lo[i] || start[i] > hi[i]) {
                throw new IllegalArgumentException("start[" + i + "] must be finite and inside [" + lo[i] + ", "
                        + hi[i] + "] : " + start[i]);
            }
        }
        if (!(epsTol > 0.0)) {
            throw new IllegalArgumentException("epsTol must be strictly positive : " + epsTol);
        }
        this.logPosterior = logPosterior;
        this.lo = lo.clone();
        this.hi = hi.clone();
        this.epsTol = epsTol;

        // ---- the peak -------------------------------------------------
        // walled off outside the box with a *finite* value, not an infinity: a
        // simplex step forms differences of function values just as a parabolic
        // one does, and Infinity - Infinity is NaN. One nat above the value at
        // the starting point is enough, since the search only ever descends
        // from there
        final double reference = logPosterior.apply(start.clone());
        if (!isFinite(reference)) {
            throw new IllegalArgumentException("the log posterior is not finite at the starting point ("
                    + start[0] + ", " + start[1] + ") : " + reference
                    + "; name a start where the posterior is positive");
        }
        final double wall = -reference + 1.0;
        DMultiFunction fenced = new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double outside = 0.0;
                for (int i = 0; i < 2; i++) {
                    if (p[i] < BivariatePosterior.this.lo[i]) {
                        outside += BivariatePosterior.this.lo[i] - p[i];
                    } else if (p[i] > BivariatePosterior.this.hi[i]) {
                        outside += p[i] - BivariatePosterior.this.hi[i];
                    }
                }
                if (outside > 0.0) {
                    return wall + outside;
                }
                double value = BivariatePosterior.this.logPosterior.apply(p);
                if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
                    return wall;
                }
                return -value;
            }
        };
        DMultiFunctionEval peak = new NelderMead().minimize(fenced, start.clone(), initialStep(lo, hi, start));
        double[] at = peak.point.clone();
        for (int i = 0; i < 2; i++) {
            if (at[i] < lo[i]) {
                at[i] = lo[i];
            }
            if (at[i] > hi[i]) {
                at[i] = hi[i];
            }
        }
        this.mode = at;
        this.logMax = logPosterior.apply(at.clone());
        if (!isFinite(logMax)) {
            throw new IllegalArgumentException("the log posterior is not finite at its peak (" + mode[0] + ", "
                    + mode[1] + ") : " + logMax);
        }

        // ---- the curvature --------------------------------------------
        // the step of a second difference wants the width of the peak, not the
        // magnitude of its position: a posterior over a few thousand
        // observations is a needle whose width has nothing to do with where it
        // sits, and a step taken from |mode| would sample where the density is
        // already zero
        double[] scale = { axisScale(0), axisScale(1) };
        this.curvature = hessian(scale);

        // ---- the whitening --------------------------------------------
        // a bounded axis goes second, because L is lower triangular: the second
        // whitened coordinate depends on the second parameter alone, so a box
        // constraint on it stays a box constraint and the region of integration
        // keeps its corners square
        int bounded0 = isBounded(0) ? 1 : 0;
        int bounded1 = isBounded(1) ? 1 : 0;
        if (bounded0 == 1 && bounded1 == 0) {
            this.first = 1;
            this.second = 0;
        } else {
            this.first = 0;
            this.second = 1;
        }
        DMatrix permuted = new DMatrix(2, 2);
        permuted.set(0, 0, curvature.get(first, first));
        permuted.set(1, 1, curvature.get(second, second));
        double off = curvature.get(first, second);
        permuted.set(0, 1, off);
        permuted.set(1, 0, off);
        DMatrix L;
        try {
            L = CholeskyDecomp.cholesky(permuted);
        } catch (RuntimeException e) {
            throw new IllegalArgumentException("the curvature at (" + mode[0] + ", " + mode[1]
                    + ") is not positive definite, so the peak found is not a maximum: " + permuted);
        }
        this.l00 = L.get(0, 0);
        this.l10 = L.get(1, 0);
        this.l11 = L.get(1, 1);
        if (!(l00 > 0.0) || !(l11 > 0.0) || !isFinite(l10)) {
            throw new IllegalArgumentException(
                    "the curvature at (" + mode[0] + ", " + mode[1] + ") could not be factored: " + permuted);
        }
        this.detL = l00 * l11;
        this.m00 = 1.0 / l00;
        this.m01 = -l10 / (l00 * l11);
        this.m11 = 1.0 / l11;

        // ---- the normalizer -------------------------------------------
        double total = integrate(Moment.ZEROTH, 0, 0.0, 0.0, first, second, l00, l10, l11, detL, m00, m01, m11);
        if (!(total > 0.0) || Double.isInfinite(total)) {
            throw new IllegalArgumentException("the posterior does not integrate to a finite positive number over ["
                    + lo[0] + ", " + hi[0] + "] x [" + lo[1] + ", " + hi[1] + "]");
        }
        this.normalizer = total;
        this.logEvidence = logMax + Math.log(total);
        // the Gaussian whose curvature this is integrates to 2 pi / det(L)
        this.logEvidenceLaplace = logMax + Math.log(2.0 * Math.PI) - Math.log(detL);

        // ---- the moments ----------------------------------------------
        this.mean = new double[2];
        for (int i = 0; i < 2; i++) {
            mean[i] = integrate(Moment.FIRST, i, 0.0, 0.0, first, second, l00, l10, l11, detL, m00, m01, m11)
                    / total;
        }
        // two passes rather than E[XY] - E[X]E[Y], which cancels for a
        // posterior that sits far from the origin
        this.covariance = new DMatrix(2, 2);
        double c00 = integrate(Moment.CENTRAL, 0, mean[0], mean[0], first, second, l00, l10, l11, detL, m00, m01,
                m11) / total;
        double c01 = integrate(Moment.CROSS, 0, mean[0], mean[1], first, second, l00, l10, l11, detL, m00, m01, m11)
                / total;
        double c11 = integrate(Moment.CENTRAL, 1, mean[1], mean[1], first, second, l00, l10, l11, detL, m00, m01,
                m11) / total;
        covariance.set(0, 0, c00);
        covariance.set(1, 1, c11);
        covariance.set(0, 1, c01);
        covariance.set(1, 0, c01);

        // ---- the second route -----------------------------------------
        // the same integral with the axes exchanged: a different Cholesky
        // factor, a different triangular substitution, a different
        // subdivision. The *unwhitened* route cannot serve as a check here,
        // because it is the thing that returns 2.9e-198
        double swapped;
        try {
            swapped = swappedNormalizer();
        } catch (RuntimeException e) {
            swapped = Double.NaN;
        }
        this.errorEstimate = Double.isNaN(swapped) ? Double.NaN : Math.abs(swapped - total) / total;
    }

    // ------------------------------------------------------------------
    // the density
    // ------------------------------------------------------------------

    /**
     * The posterior density at {@code x}.
     * <p>
     * The exponential of {@link #logPdf(double[])}, which is where the accuracy
     * is: a posterior over a few thousand observations has a density below the
     * smallest {@code double} everywhere except within a few widths of its
     * mode.
     *
     * @param x
     *            the two coordinates. Not modified
     * @return the density there, {@code 0.0} outside the box
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or does not hold two entries
     */
    public double pdf(double[] x) {
        return Math.exp(logPdf(x));
    }

    /**
     * The natural logarithm of the posterior density at {@code x}.
     * <p>
     * The log posterior as it was given, less {@link #logEvidence()}. The one
     * member that costs a single evaluation of the caller's function and no
     * integration at all.
     *
     * @param x
     *            the two coordinates. Not modified
     * @return the log density there, {@link Double#NEGATIVE_INFINITY} outside
     *         the box
     * @throws IllegalArgumentException
     *             if {@code x} is {@code null} or does not hold two entries
     */
    public double logPdf(double[] x) {
        if (x == null || x.length != 2) {
            throw new IllegalArgumentException("x must hold two coordinates");
        }
        for (int i = 0; i < 2; i++) {
            if (Double.isNaN(x[i]) || x[i] < lo[i] || x[i] > hi[i]) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return logPosterior.apply(x.clone()) - logEvidence;
    }

    /**
     * The natural logarithm of the normalizing constant, which for a posterior
     * built as prior times likelihood is the <b>log marginal likelihood</b> of
     * the data -- the evidence, the thing a Bayes factor is a ratio of.
     *
     * @return the natural logarithm of the integral of
     *         {@code exp(logPosterior)} over the box
     */
    public double logEvidence() {
        return logEvidence;
    }

    /**
     * The Laplace approximation to {@link #logEvidence()}:
     * {@code logMax + log(2 pi) - log(det L)}, where {@code L L^T} is the
     * curvature at the mode.
     * <p>
     * It costs nothing, because the whitening had to factor that curvature
     * anyway. <b>The gap between the two is a diagnostic</b>: measured, it is
     * around {@code 1e-8} for a posterior that really is Gaussian -- which is
     * the accuracy of the difference-quotient curvature and not a real
     * difference -- and {@code 4.7e-4} for a two-parameter posterior over two
     * thousand observations, whose scale parameter is genuinely skew. A gap
     * far above {@code 1e-7} means the approximation is missing something the
     * quadrature saw.
     *
     * @return the Laplace approximation to the log evidence
     */
    public double logEvidenceLaplace() {
        return logEvidenceLaplace;
    }

    /**
     * Writes the mode, where the posterior peaks, into {@code out}.
     *
     * @param out
     *            where the result is written, two entries. Its previous
     *            contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or does not hold two entries
     */
    public void mode(double[] out) {
        checkOut(out);
        out[0] = mode[0];
        out[1] = mode[1];
    }

    /**
     * The curvature the parameters were whitened by: the Hessian of the
     * <em>negative</em> log posterior at the mode, in the caller's axis order,
     * by central differences.
     * <p>
     * Its inverse is the covariance the Laplace approximation would report, so
     * comparing it with {@link #covariance()} is the same diagnostic
     * {@link #logEvidenceLaplace()} gives for the evidence.
     *
     * @return a fresh {@code 2 x 2} curvature matrix
     */
    public DMatrix curvature() {
        return curvature.copy();
    }

    /**
     * The relative gap between the two routes to the normalizing constant: the
     * whitened integral as it stands, and the same integral with the two axes
     * exchanged.
     * <p>
     * Exchanging the axes gives a different Cholesky factor and therefore a
     * different substitution and a different subdivision, so agreement is
     * evidence about the answer rather than a restatement of the tolerance that
     * was asked for. It is what would notice a second peak the whitening was
     * not laid around. <b>The unwhitened integral cannot serve as this check</b>
     * -- it is the route that fails silently, which is why the class exists.
     *
     * @return the relative gap, or {@link Double#NaN} if the second route could
     *         not be computed, which is a failure of the check and not of the
     *         answer
     */
    public double errorEstimate() {
        return errorEstimate;
    }

    // ------------------------------------------------------------------
    // the moments
    // ------------------------------------------------------------------

    /**
     * Writes the posterior mean of each parameter into {@code out}.
     *
     * @param out
     *            where the result is written, two entries. Its previous
     *            contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or does not hold two entries
     */
    public void mean(double[] out) {
        checkOut(out);
        out[0] = mean[0];
        out[1] = mean[1];
    }

    /**
     * The posterior covariance matrix of the two parameters.
     * <p>
     * Summed in two passes rather than as {@code E[XY] - E[X]E[Y]}, which
     * cancels for a posterior far from the origin. Where the posterior is
     * Gaussian this is the inverse of {@link #curvature()}; where it is not,
     * the difference between the two is what the Laplace approximation misses.
     * <p>
     * A moment that does not exist is <b>not detected</b>, for the reason
     * {@link ScalarPosterior#variance()} sets out at length: the rule truncates
     * a divergent tail somewhere and reports what it accumulated.
     *
     * @return a fresh {@code 2 x 2} covariance matrix, exactly symmetric
     */
    public DMatrix covariance() {
        return covariance.copy();
    }

    /**
     * The marginal density of parameter {@code i} at {@code x}: the posterior
     * with the other parameter integrated out.
     * <p>
     * <b>This is the expensive member.</b> One call is a whole one-dimensional
     * integration -- measured at about 1400 evaluations of the caller's
     * function and 16 ms -- because a marginal density is an integral and there
     * is nothing precomputed to look it up in. Wrapping it in a
     * {@link ScalarPosterior} to get a marginal <em>law</em> costs on the order
     * of a minute, and the library does not hide that by doing it for you.
     *
     * @param i
     *            which parameter, {@code 0} or {@code 1}
     * @param x
     *            where to evaluate its marginal density
     * @return the marginal density, {@code 0.0} outside that parameter's range
     * @throws IllegalArgumentException
     *             if {@code i} is neither {@code 0} nor {@code 1}
     */
    public double marginalDensity(int i, final double x) {
        if (i != 0 && i != 1) {
            throw new IllegalArgumentException("there is no parameter " + i);
        }
        if (Double.isNaN(x) || x < lo[i] || x > hi[i]) {
            return 0.0;
        }
        final int free = 1 - i;
        final int fixed = i;
        // the conditional mode of the free axis under the Gaussian the
        // whitening fits, which is where the mass of the inner integral is
        final double center = conditionalMode(fixed, free, x);
        // and its conditional width, which the same curvature gives. Naming the
        // center is not enough here: the substitution has derivative one there,
        // so a feature far narrower than one unit stays far narrower than one
        // unit and no subdivision finds it. Measured, a conditional width of
        // 4.5e-5 returned zero. Scaling the axis first makes it unit width,
        // which is the whitening of the one dimension that is left
        final double width = conditionalWidth(free);
        DFunction inner = new DFunction() {
            @Override
            public double apply(double s) {
                double t = center + width * s;
                if (t < lo[free] || t > hi[free]) {
                    return 0.0;
                }
                double[] p = new double[2];
                p[fixed] = x;
                p[free] = t;
                double value = logPosterior.apply(p);
                if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
                    return 0.0;
                }
                return Math.exp(value - logMax) * width;
            }
        };
        double sLo = isFinite(lo[free]) ? (lo[free] - center) / width : Double.NEGATIVE_INFINITY;
        double sHi = isFinite(hi[free]) ? (hi[free] - center) / width : Double.POSITIVE_INFINITY;
        double integral = Quadrature.integrate(inner, sLo, sHi, epsTol, centerInside(sLo, sHi));
        return integral / normalizer;
    }

    // ------------------------------------------------------------------
    // the integration
    // ------------------------------------------------------------------

    /** Which quantity a whitened integration is to take. */
    private enum Moment {
        ZEROTH, FIRST, CENTRAL, CROSS
    }

    /**
     * One whitened integration. Two routes, and the choice between them is
     * structural rather than a heuristic:
     * <ul>
     * <li>at most one axis bounded -- the constraint is a box in {@code u},
     * because {@code L} is triangular, so the native two-dimensional rule
     * integrates the exact region;</li>
     * <li>both axes bounded -- the first parameter's constraint cuts
     * diagonally, so the region is a parallelogram and the integration is
     * nested, {@code u1} outside and {@code u0} between two lines whose
     * equations are closed form.</li>
     * </ul>
     * An indicator inside the integrand was measured and rejected: it costs 5
     * to 7 decimal digits, because an adaptive rule accepts a panel straddling
     * a jump without noticing it.
     */
    private double integrate(final Moment moment, final int which, final double centerA, final double centerB,
            final int fi, final int se, final double a00, final double a10, final double a11, final double dL,
            final double b00, final double b01, final double b11) {
        final double modeFirst = mode[fi];
        final double modeSecond = mode[se];
        final double loFirst = lo[fi];
        final double hiFirst = hi[fi];
        final double loSecond = lo[se];
        final double hiSecond = hi[se];

        final double u1lo = isFinite(loSecond) ? a11 * (loSecond - modeSecond) : Double.NEGATIVE_INFINITY;
        final double u1hi = isFinite(hiSecond) ? a11 * (hiSecond - modeSecond) : Double.POSITIVE_INFINITY;

        if (isFinite(loFirst) || isFinite(hiFirst)) {
            // The first parameter is constrained, and its constraint is the one
            // that cuts diagonally in u -- so the region is a parallelogram and
            // the limits have to be carried, not tested for. Either end alone
            // is enough: falling back to the native rule with a membership test
            // was measured at 6.2e-06 where this reaches 1e-15, because an
            // adaptive rule accepts a panel straddling a jump without noticing
            return Quadrature.integrate(new DFunction() {
                @Override
                public double apply(double u1) {
                    final double uu = u1;
                    double a0 = isFinite(loFirst) ? (loFirst - modeFirst - b01 * uu) / b00
                            : Double.NEGATIVE_INFINITY;
                    double b0 = isFinite(hiFirst) ? (hiFirst - modeFirst - b01 * uu) / b00
                            : Double.POSITIVE_INFINITY;
                    if (!(a0 < b0)) {
                        return 0.0;
                    }
                    return Quadrature.integrate(new DFunction() {
                        @Override
                        public double apply(double u0) {
                            return value(u0, uu, moment, which, centerA, centerB, fi, se, modeFirst, modeSecond, dL,
                                    b00, b01, b11);
                        }
                    }, a0, b0, epsTol, centerInside(a0, b0));
                }
            }, u1lo, u1hi, epsTol, centerInside(u1lo, u1hi));
        }
        return Quadrature.integrate((u0, u1) -> value(u0, u1, moment, which, centerA, centerB, fi, se, modeFirst,
                modeSecond, dL, b00, b01, b11), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, u1lo, u1hi,
                epsTol, 0.0, centerInside(u1lo, u1hi));
    }

    /**
     * The integrand: the shifted density at the point {@code u} maps to, times
     * whatever the moment asks for. Allocates its own array, because the rule
     * evaluates this from several threads at once.
     */
    private double value(double u0, double u1, Moment moment, int which, double centerA, double centerB, int fi,
            int se, double modeFirst, double modeSecond, double dL, double b00, double b01, double b11) {
        double zFirst = modeFirst + b00 * u0 + b01 * u1;
        double zSecond = modeSecond + b11 * u1;
        if (zFirst < lo[fi] || zFirst > hi[fi] || zSecond < lo[se] || zSecond > hi[se]) {
            return 0.0;
        }
        double[] p = new double[2];
        p[fi] = zFirst;
        p[se] = zSecond;
        double logValue = logPosterior.apply(p);
        if (Double.isNaN(logValue) || logValue == Double.NEGATIVE_INFINITY) {
            return 0.0;
        }
        double density = Math.exp(logValue - logMax) / dL;
        switch (moment) {
        case FIRST:
            return p[which] * density;
        case CENTRAL: {
            double d = p[which] - centerA;
            return d * d * density;
        }
        case CROSS:
            return (p[0] - centerA) * (p[1] - centerB) * density;
        default:
            return density;
        }
    }

    /** The normalizer again, with the two axes exchanged. */
    private double swappedNormalizer() {
        int fi = second;
        int se = first;
        DMatrix permuted = new DMatrix(2, 2);
        permuted.set(0, 0, curvature.get(fi, fi));
        permuted.set(1, 1, curvature.get(se, se));
        double off = curvature.get(fi, se);
        permuted.set(0, 1, off);
        permuted.set(1, 0, off);
        DMatrix L = CholeskyDecomp.cholesky(permuted);
        double a00 = L.get(0, 0);
        double a10 = L.get(1, 0);
        double a11 = L.get(1, 1);
        double dL = a00 * a11;
        return integrate(Moment.ZEROTH, 0, 0.0, 0.0, fi, se, a00, a10, a11, dL, 1.0 / a00, -a10 / (a00 * a11),
                1.0 / a11);
    }

    // ------------------------------------------------------------------
    // construction helpers
    // ------------------------------------------------------------------

    /**
     * The width of the peak along one axis: the distance over which the log
     * posterior falls by half a nat, found by doubling outward. Used only to
     * scale a difference quotient, so a factor of two either way is harmless.
     */
    private double axisScale(int axis) {
        double target = logMax - 0.5;
        double reachUp = isFinite(hi[axis]) ? hi[axis] - mode[axis] : Double.POSITIVE_INFINITY;
        double reachDown = isFinite(lo[axis]) ? mode[axis] - lo[axis] : Double.POSITIVE_INFINITY;
        double reach = Math.max(reachUp, reachDown);
        double direction = reachUp >= reachDown ? 1.0 : -1.0;
        if (!(reach > 0.0)) {
            return Math.max(Math.abs(mode[axis]), 1.0) * EPS_FOURTH;
        }
        double d = Math.min(1.0, reach) / 64.0;
        if (!(d > 0.0)) {
            d = 1.0;
        }
        double[] p = new double[2];
        for (int step = 0; step < MAX_SCALE_STEPS; step++) {
            double far = Math.min(d, reach);
            p[0] = mode[0];
            p[1] = mode[1];
            p[axis] = mode[axis] + direction * far;
            double value = logPosterior.apply(p);
            if (Double.isNaN(value) || value <= target) {
                return far;
            }
            if (far == reach) {
                return reach;
            }
            d += d;
        }
        throw new IllegalArgumentException("the log posterior does not fall away from its peak along axis " + axis
                + "; it cannot be normalized");
    }

    /**
     * The Hessian of the negative log posterior at the mode, by central
     * differences.
     * <p>
     * Each off-diagonal is computed once and written into both places. That is
     * not tidiness: {@link CholeskyDecomp} tests symmetry with an
     * <em>absolute</em> tolerance of {@code 1e-10}, so a curvature whose
     * entries are of order {@code 1e6} -- an ordinary posterior over a few
     * thousand observations -- would be rejected over its own round-off if the
     * two halves were computed separately.
     */
    private DMatrix hessian(double[] scale) {
        double h0 = EPS_FOURTH * scale[0];
        double h1 = EPS_FOURTH * scale[1];
        double f0 = -logMax;
        double dxx = (negAt(h0, 0.0) - 2.0 * f0 + negAt(-h0, 0.0)) / (h0 * h0);
        double dyy = (negAt(0.0, h1) - 2.0 * f0 + negAt(0.0, -h1)) / (h1 * h1);
        double dxy = (negAt(h0, h1) - negAt(h0, -h1) - negAt(-h0, h1) + negAt(-h0, -h1)) / (4.0 * h0 * h1);
        DMatrix h = new DMatrix(2, 2);
        h.set(0, 0, dxx);
        h.set(1, 1, dyy);
        h.set(0, 1, dxy);
        h.set(1, 0, dxy);
        return h;
    }

    /** The negative log posterior a step away from the mode, clamped into the box. */
    private double negAt(double d0, double d1) {
        double[] p = { mode[0] + d0, mode[1] + d1 };
        for (int i = 0; i < 2; i++) {
            if (p[i] < lo[i]) {
                p[i] = lo[i];
            }
            if (p[i] > hi[i]) {
                p[i] = hi[i];
            }
        }
        double value = logPosterior.apply(p);
        if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
            return Double.MAX_VALUE / 8.0;
        }
        return -value;
    }

    /**
     * The conditional mode of the free axis given the fixed one, under the
     * Gaussian the whitening fits. One multiplication, because the covariance
     * of that Gaussian is the inverse of a triangular product.
     */
    private double conditionalMode(int fixed, int free, double x) {
        double hff = curvature.get(free, free);
        double hfx = curvature.get(free, fixed);
        if (!(hff > 0.0)) {
            return mode[free];
        }
        double candidate = mode[free] - hfx * (x - mode[fixed]) / hff;
        if (candidate <= lo[free] || candidate >= hi[free] || !isFinite(candidate)) {
            return centerInside(lo[free], hi[free]);
        }
        return candidate;
    }

    /**
     * The conditional width of the free axis under the same Gaussian, which is
     * what the inner integral of a marginal has to be scaled by.
     */
    private double conditionalWidth(int free) {
        double hff = curvature.get(free, free);
        if (!(hff > 0.0) || !isFinite(hff)) {
            return 1.0;
        }
        double width = 1.0 / Math.sqrt(hff);
        if (!(width > 0.0) || !isFinite(width)) {
            return 1.0;
        }
        return width;
    }

    /** A point strictly inside the interval, for the centered substitution. */
    private static double centerInside(double a, double b) {
        if (a < 0.0 && b > 0.0) {
            return 0.0;
        }
        if (isFinite(a) && isFinite(b)) {
            return 0.5 * (a + b);
        }
        if (isFinite(a)) {
            return a + 1.0;
        }
        if (isFinite(b)) {
            return b - 1.0;
        }
        return 0.0;
    }

    private boolean isBounded(int axis) {
        return isFinite(lo[axis]) || isFinite(hi[axis]);
    }

    private static void checkBox(double[] lo, double[] hi) {
        if (lo == null || hi == null) {
            throw new IllegalArgumentException("lo and hi must not be null");
        }
        if (lo.length != 2 || hi.length != 2) {
            throw new IllegalArgumentException("lo and hi must hold two entries each");
        }
        for (int i = 0; i < 2; i++) {
            if (Double.isNaN(lo[i]) || Double.isNaN(hi[i])) {
                throw new IllegalArgumentException("the range of axis " + i + " must not be NaN");
            }
            if (!(lo[i] < hi[i])) {
                throw new IllegalArgumentException(
                        "lo[" + i + "] must be below hi[" + i + "] : [" + lo[i] + ", " + hi[i] + "]");
            }
        }
    }

    private static void checkOut(double[] out) {
        if (out == null || out.length != 2) {
            throw new IllegalArgumentException("out must hold two coordinates");
        }
    }

    /** Where to start looking for the peak when the caller does not say. */
    private static double[] defaultStart(double[] lo, double[] hi) {
        double[] start = new double[2];
        for (int i = 0; i < 2; i++) {
            if (isFinite(lo[i]) && isFinite(hi[i])) {
                start[i] = lo[i] + (hi[i] - lo[i]) / 2.0;
            } else if (isFinite(lo[i])) {
                start[i] = lo[i] + 1.0;
            } else if (isFinite(hi[i])) {
                start[i] = hi[i] - 1.0;
            } else {
                start[i] = 0.0;
            }
        }
        return start;
    }

    /** The first simplex edge, small enough to stay in the box. */
    private static double[] initialStep(double[] lo, double[] hi, double[] start) {
        double[] step = new double[2];
        for (int i = 0; i < 2; i++) {
            double span = hi[i] - lo[i];
            double s = isFinite(span) ? span / 100.0 : 1.0;
            if (!(s > 0.0)) {
                s = 1.0;
            }
            if (start[i] + s > hi[i]) {
                s = -s;
            }
            step[i] = s;
        }
        return step;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
