package math.stats.bayes;

import math.MathConsts;
import math.distribution.ContinuousDistribution;
import math.fun.DFunction;
import math.optim.BrentMinimizer;
import math.solve.Quadrature;
import math.solve.RootFinder;

/**
 * The posterior of a single parameter, normalized by quadrature.
 * <p>
 * A caller hands over the <b>logarithm</b> of an unnormalized posterior --
 * typically a log prior plus a log likelihood, both of which every distribution
 * in {@code math.distribution} can supply -- together with the range the
 * parameter lives in, and gets back an ordinary
 * {@link ContinuousDistribution}: density, distribution function, quantile,
 * mean and variance, and with them credible intervals. Nothing is sampled and
 * there is no chain to diagnose; the normalizing constant is an integral, and
 * an integral comes with an estimate of its own error.
 * <p>
 * <b>Why the logarithm.</b> Quadrature needs a density, but a posterior over
 * more than a few hundred observations has an unnormalized density that is
 * zero in every {@code double} at every point a rule will evaluate. So the peak
 * is located first and the integrand is {@code exp(logPosterior - logMax)},
 * whose maximum is exactly one whatever the scale of the problem;
 * {@code logMax} goes back into {@link #logEvidence()} afterwards. This is the
 * log-sum-exp trick applied to an integral, and without it everything here
 * would be {@code 0/0}.
 * <p>
 * <b>What construction costs.</b> Everything is computed once, in
 * {@code of}: the mode, the width of the posterior on either side of it, one
 * integration per panel of a fixed ladder around the mode, and both moments.
 * That is a few dozen integrations, so an instance is expensive to build and
 * cheap to ask. It is immutable and can be shared between threads.
 * <p>
 * <b>What it assumes.</b> That the posterior is proper -- integrable over the
 * range given -- and that it has one peak. A multimodal posterior is not
 * detected: the mode search finds one of the peaks and the panel ladder is
 * laid around it, so the others are covered only by the outer panels. The
 * agreement reported by {@link #errorEstimate()} is the check on that, since
 * the second route to the normalizer subdivides the range differently.
 *
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Bayesian_inference">Wikipedia Bayesian
 *      inference</a>
 * @since 1.5.3
 */
public final class ScalarPosterior implements ContinuousDistribution {

    /**
     * Panel edges are placed at these multiples of the posterior's half-width
     * on either side of the mode. Dense where the mass is, sparse where the
     * tail only has to be caught.
     */
    private static final double[] LADDER = { 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0,
            15.0, 20.0 };

    /** Error tolerance handed to {@link Quadrature} if none is given. */
    private static final double DEFAULT_EPS_TOL = 1.0e-10;

    /**
     * How far the search for a half-width may step out from the mode, in
     * multiples of its starting step, before the posterior is called improper.
     */
    private static final int MAX_WIDTH_STEPS = 200;

    /**
     * The square root of the machine epsilon, about {@code 1.05e-8}. Nothing
     * finer is worth asking of a search for the width of a peak, for the reason
     * {@link BrentMinimizer} gives about its own tolerance: a function is
     * quadratically flat at its extremum.
     */
    private static final double SQRT_EPS = Math.sqrt(MathConsts.MACH_EPS_DBL);

    private final DFunction logPosterior;
    private final double lo;
    private final double hi;
    private final double epsTol;

    private final double mode;
    /** {@code logPosterior(mode)}, the shift that keeps the integrand alive. */
    private final double logMax;
    private final double leftScale;
    private final double rightScale;

    /** The panel edges, ascending, {@code edges[0] == lo} and last {@code == hi}. */
    private final double[] edges;
    /** The distribution function at each edge, from {@code 0.0} to {@code 1.0}. */
    private final double[] cumulative;
    /** The integral of the shifted integrand, before the shift is undone. */
    private final double normalizer;

    private final double logEvidence;
    private final double mean;
    private final double variance;
    private final double errorEstimate;

    /** An interval and the posterior mass it holds. */
    public static final class Interval {

        /** The lower end. */
        public final double lower;
        /** The upper end. */
        public final double upper;
        /** The mass between them, which is what was asked for up to rounding. */
        public final double mass;

        Interval(double lower, double upper, double mass) {
            this.lower = lower;
            this.upper = upper;
            this.mass = mass;
        }

        /**
         * The width of the interval.
         *
         * @return {@code upper - lower}
         */
        public double width() {
            return upper - lower;
        }

        @Override
        public String toString() {
            return "[" + lower + ", " + upper + "] holding " + mass;
        }
    }

    /**
     * The posterior of {@code logPosterior} over {@code [lo, hi]}, with a
     * starting point and an error tolerance derived from the range.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density. It
     *            is never evaluated outside {@code [lo, hi]}, and it may return
     *            {@link Double#NEGATIVE_INFINITY} where the posterior is zero
     * @param lo
     *            the lower end of the range, possibly
     *            {@link Double#NEGATIVE_INFINITY}
     * @param hi
     *            the upper end of the range, possibly
     *            {@link Double#POSITIVE_INFINITY}
     * @return the normalized posterior
     * @throws IllegalArgumentException
     *             if {@code logPosterior} is {@code null}, if {@code lo} is not
     *             below {@code hi}, if either is {@code NaN}, if no peak can be
     *             bracketed, or if the posterior does not fall away from its
     *             peak and so cannot be normalized
     */
    public static ScalarPosterior of(DFunction logPosterior, double lo, double hi) {
        return new ScalarPosterior(logPosterior, lo, hi, defaultStart(lo, hi), DEFAULT_EPS_TOL);
    }

    /**
     * The posterior of {@code logPosterior} over {@code [lo, hi]}, starting the
     * search for the peak at {@code start}.
     * <p>
     * Naming {@code start} is what rescues a posterior whose mass sits far from
     * anything the short form would guess at -- the search walks outward from
     * it and has to find the peak before the peak underflows.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density
     * @param lo
     *            the lower end of the range, possibly
     *            {@link Double#NEGATIVE_INFINITY}
     * @param hi
     *            the upper end of the range, possibly
     *            {@link Double#POSITIVE_INFINITY}
     * @param start
     *            a finite point inside {@code [lo, hi]} near the mass
     * @param epsTol
     *            the error tolerance of each integration, strictly positive
     * @return the normalized posterior
     * @throws IllegalArgumentException
     *             under the conditions {@link #of(DFunction, double, double)}
     *             states, and additionally if {@code start} is not finite or
     *             lies outside {@code [lo, hi]}, or if {@code epsTol} is not
     *             strictly positive
     */
    public static ScalarPosterior of(DFunction logPosterior, double lo, double hi, double start, double epsTol) {
        return new ScalarPosterior(logPosterior, lo, hi, start, epsTol);
    }

    private ScalarPosterior(DFunction logPosterior, double lo, double hi, double start, double epsTol) {
        if (logPosterior == null) {
            throw new IllegalArgumentException("logPosterior must not be null");
        }
        if (Double.isNaN(lo) || Double.isNaN(hi)) {
            throw new IllegalArgumentException("the range must not be NaN : [" + lo + ", " + hi + "]");
        }
        if (!(lo < hi)) {
            throw new IllegalArgumentException("lo must be below hi : [" + lo + ", " + hi + "]");
        }
        if (!isFinite(start) || start < lo || start > hi) {
            throw new IllegalArgumentException("start must be finite and inside [" + lo + ", " + hi + "] : " + start);
        }
        if (!(epsTol > 0.0)) {
            throw new IllegalArgumentException("epsTol must be strictly positive : " + epsTol);
        }
        this.logPosterior = logPosterior;
        this.lo = lo;
        this.hi = hi;
        this.epsTol = epsTol;

        // ---- the peak -------------------------------------------------
        // minimized rather than maximized because that is what BrentMinimizer
        // does, and walled off outside the range so that the downhill walk
        // turns back at the ends instead of evaluating the caller's function
        // where it was never promised to be defined.
        //
        // The wall is finite, and that is not a detail: an infinity here poisons
        // the parabolic extrapolation of BrentMinimizer.bracket, whose step
        // forms differences of function values, so Infinity - Infinity makes it
        // NaN and no peak is ever bracketed. A wall one nat above the value at
        // the starting point is enough, since the walk only ever moves downhill
        // from there and so cannot mistake the wall for a descent
        final double reference = logPosterior.apply(start);
        if (!isFinite(reference)) {
            throw new IllegalArgumentException("the log posterior is not finite at the starting point " + start
                    + " : " + reference + "; name a start where the posterior is positive");
        }
        final double wall = -reference + 1.0;
        DFunction fenced = new DFunction() {
            @Override
            public double apply(double x) {
                if (x < ScalarPosterior.this.lo) {
                    return wall + (ScalarPosterior.this.lo - x);
                }
                if (x > ScalarPosterior.this.hi) {
                    return wall + (x - ScalarPosterior.this.hi);
                }
                double value = ScalarPosterior.this.logPosterior.apply(x);
                if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
                    return wall;
                }
                return -value;
            }
        };
        BrentMinimizer minimizer = new BrentMinimizer();
        double step = initialStep(lo, hi, start);
        BrentMinimizer.Bracket bracket = minimizer.bracket(fenced, start, start + step);
        if (!bracket.bracketed) {
            throw new IllegalArgumentException("no peak could be bracketed from start = " + start
                    + "; name a start point near the mass, or check that the posterior is proper");
        }
        BrentMinimizer.Result peak = minimizer.minimize(fenced, bracket);
        double peakAt = peak.x;
        if (peakAt < lo) {
            peakAt = lo;
        }
        if (peakAt > hi) {
            peakAt = hi;
        }
        this.mode = peakAt;
        this.logMax = logPosterior.apply(peakAt);
        if (!isFinite(logMax)) {
            throw new IllegalArgumentException(
                    "the log posterior is not finite at its peak " + peakAt + " : " + logMax);
        }

        // ---- the width on either side ---------------------------------
        // half a nat below the peak, which for a normal is exactly one standard
        // deviation. Taken separately left and right, so that a skewed
        // posterior does not get a grid built for a symmetric one
        double probe = Math.abs(bracket.c - bracket.a);
        if (!isFinite(probe) || probe <= 0.0) {
            probe = 1.0;
        }
        this.leftScale = halfWidth(-1.0, probe);
        this.rightScale = halfWidth(1.0, probe);

        // ---- the panels -----------------------------------------------
        this.edges = buildEdges();
        int panels = edges.length - 1;
        double[] mass = new double[panels];
        double total = 0.0;
        for (int i = 0; i < panels; i++) {
            mass[i] = panelIntegral(edges[i], edges[i + 1], Moment.ZEROTH, 0.0);
            total += mass[i];
        }
        if (!(total > 0.0) || Double.isInfinite(total)) {
            throw new IllegalArgumentException(
                    "the posterior does not integrate to a finite positive number over [" + lo + ", " + hi + "]");
        }
        this.normalizer = total;
        this.cumulative = new double[edges.length];
        double running = 0.0;
        for (int i = 0; i < panels; i++) {
            running += mass[i];
            cumulative[i + 1] = running / total;
        }
        // assigned rather than accumulated, so that the last panel closes the
        // distribution function exactly
        cumulative[panels] = 1.0;
        this.logEvidence = logMax + Math.log(total);

        // ---- the moments ----------------------------------------------
        double first = 0.0;
        for (int i = 0; i < panels; i++) {
            first += panelIntegral(edges[i], edges[i + 1], Moment.FIRST, 0.0);
        }
        this.mean = first / total;
        // two passes rather than E[X*X] - E[X]*E[X], which cancels for a
        // posterior that sits far from the origin -- and one located at 1e4 is
        // exactly the case this class is built to survive
        double second = 0.0;
        for (int i = 0; i < panels; i++) {
            second += panelIntegral(edges[i], edges[i + 1], Moment.CENTRAL, mean);
        }
        this.variance = second / total;

        // ---- the second route to the same number ----------------------
        // one adaptive integration over the whole range, subdivided by the
        // rule rather than by the ladder above. The gap between the two is an
        // error estimate that owes nothing to epsTol
        double whole;
        try {
            whole = Quadrature.integrate(shifted(Moment.ZEROTH, 0.0), lo, hi, epsTol, mode);
        } catch (ArithmeticException e) {
            // the substitution could not find the mass even told where it is,
            // which the panel sweep did not have to do; report it rather than
            // fail the construction over a cross-check
            whole = Double.NaN;
        }
        this.errorEstimate = Double.isNaN(whole) ? Double.NaN : Math.abs(whole - total) / total;
    }

    // ------------------------------------------------------------------
    // the density
    // ------------------------------------------------------------------

    /**
     * {@inheritDoc}
     * <p>
     * The exponential of {@link #logPdf(double)}, which is where the accuracy
     * is: a posterior over a few thousand observations has a density below the
     * smallest {@code double} everywhere except within a few widths of its
     * mode.
     */
    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    /**
     * {@inheritDoc}
     * <p>
     * The log posterior as it was given, less {@link #logEvidence()}. This is
     * the one member that costs a single evaluation of the caller's function
     * and no integration at all.
     */
    @Override
    public double logPdf(double x) {
        if (x < lo || x > hi || Double.isNaN(x)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logPosterior.apply(x) - logEvidence;
    }

    /**
     * The natural logarithm of the normalizing constant, which for a posterior
     * built as prior times likelihood is the <b>log marginal likelihood</b> of
     * the data -- the evidence, the thing a Bayes factor is a ratio of.
     * <p>
     * It is {@code 0} for a density that was already normalized, and it comes
     * back exactly as given for a constant added to the log posterior, which is
     * what makes it usable for model comparison.
     *
     * @return the natural logarithm of the integral of {@code exp(logPosterior)}
     */
    public double logEvidence() {
        return logEvidence;
    }

    /**
     * Where the posterior peaks.
     *
     * @return the mode
     */
    public double mode() {
        return mode;
    }

    /**
     * The relative gap between the two routes to the normalizing constant: the
     * panel sweep this class lays out, and one adaptive integration over the
     * whole range.
     * <p>
     * The two subdivide the same integral differently, so agreement is evidence
     * about the answer rather than a restatement of the tolerance that was
     * asked for. It is the one number here that would notice a second peak the
     * ladder was not laid around.
     *
     * @return the relative gap, or {@link Double#NaN} if the single integration
     *         could not find the mass at all -- which is not a failure of the
     *         answer, only of the check
     */
    public double errorEstimate() {
        return errorEstimate;
    }

    /**
     * The half-width of the posterior to the left of its mode: the distance
     * over which the log posterior falls by half a nat, which for a normal is
     * exactly one standard deviation.
     *
     * @return the left half-width, strictly positive
     */
    public double leftScale() {
        return leftScale;
    }

    /**
     * The half-width of the posterior to the right of its mode. It differs from
     * {@link #leftScale()} exactly to the extent the posterior is skewed.
     *
     * @return the right half-width, strictly positive
     */
    public double rightScale() {
        return rightScale;
    }

    @Override
    public double supportLowerBound() {
        return lo;
    }

    @Override
    public double supportUpperBound() {
        return hi;
    }

    // ------------------------------------------------------------------
    // the distribution function
    // ------------------------------------------------------------------

    /**
     * {@inheritDoc}
     * <p>
     * The cumulative mass of the panels below {@code x} plus one integration
     * over the single panel {@code x} falls inside. The panel sums are formed
     * once and are non-decreasing by construction, and the partial panel is
     * clamped into the gap it belongs to, so this is monotone <b>exactly</b>
     * and not merely up to the tolerance of a quadrature -- which is what makes
     * {@link #inverseCdf(double)} a well posed question in the tails.
     */
    @Override
    public double cdf(double x) {
        if (Double.isNaN(x)) {
            return Double.NaN;
        }
        if (x <= lo) {
            return 0.0;
        }
        if (x >= hi) {
            return 1.0;
        }
        int j = panelOf(x);
        if (x == edges[j]) {
            return cumulative[j];
        }
        double partial = panelIntegral(edges[j], x, Moment.ZEROTH, 0.0) / normalizer;
        double value = cumulative[j] + partial;
        if (value < cumulative[j]) {
            return cumulative[j];
        }
        if (value > cumulative[j + 1]) {
            return cumulative[j + 1];
        }
        return value;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The panel table locates the answer to within one panel by a binary
     * search, and {@link RootFinder#brentDekker} finds it inside that panel.
     * The Newton iteration of
     * {@link ContinuousDistribution#findRoot(double, double, double, double)}
     * is not used: every one of its steps would cost an integration, and it has
     * no bracket to fall back on.
     */
    @Override
    public double inverseCdf(double probability) {
        double p = probability;
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("probability must be in [0, 1] : " + p);
        }
        if (p == 0.0) {
            return lo;
        }
        if (p == 1.0) {
            return hi;
        }
        int j = 0;
        int high = cumulative.length - 1;
        while (j < high) {
            int mid = (j + high) >>> 1;
            if (cumulative[mid] >= p) {
                high = mid;
            } else {
                j = mid + 1;
            }
        }
        // j is the first edge whose mass reaches p, so the answer is in the
        // panel below it
        int panel = j - 1;
        if (panel < 0) {
            panel = 0;
        }
        double a = edges[panel];
        double b = edges[panel + 1];
        if (!isFinite(a)) {
            a = walkOut(b, -1.0, p);
        }
        if (!isFinite(b)) {
            b = walkOut(a, 1.0, p);
        }
        final double target = p;
        // the tightest the root search can usefully be asked for, rather than
        // epsTol: the distribution function is good to a few ulps here because
        // the panel sums are formed once, so the quantile can be too
        return RootFinder.brentDekker(a, b, new DFunction() {
            @Override
            public double apply(double x) {
                return cdf(x) - target;
            }
        }, MathConsts.MIN_TOL);
    }

    // ------------------------------------------------------------------
    // the moments
    // ------------------------------------------------------------------

    @Override
    public double mean() {
        return mean;
    }

    /**
     * {@inheritDoc}
     * <p>
     * <b>A variance that does not exist is not detected, and this is the one
     * place where a number comes back that should not.</b> The second moment is
     * summed over the same panels as everything else, and where its integral
     * diverges the rule truncates it somewhere in the outermost panel and hands
     * back whatever it accumulated -- a finite number that says more about
     * where the tail was cut than about the posterior. Measured on a Student t
     * with two degrees of freedom, whose variance is infinite, the value over
     * {@code [-R, R]} grows by {@code log(100)} for every <em>decade</em> of
     * {@code R}, which is the signature of a logarithmic divergence being cut
     * off rather than integrated: {@code 7.90} at {@code 1e2}, {@code 12.51} at
     * {@code 1e3}, {@code 17.11} at {@code 1e4}, {@code 26.32} at {@code 1e6}.
     * Over the whole line it reports {@code 38.2}, which is simply where the
     * rule happened to stop.
     * <p>
     * {@link #errorEstimate()} does not catch it either, since that watches the
     * normalizing constant, which converges perfectly well for such a
     * posterior. For a tail heavy enough to be in doubt, read a credible
     * interval instead -- {@link #credibleInterval(double)} and
     * {@link #highestDensityInterval(double)} are defined whether a second
     * moment exists or not.
     *
     * @return the variance of the posterior, or a truncation of it where the
     *         posterior has none
     */
    @Override
    public double variance() {
        return variance;
    }

    // ------------------------------------------------------------------
    // credible intervals
    // ------------------------------------------------------------------

    /**
     * The equal-tailed credible interval at the given level: the two quantiles
     * that leave the same mass outside on either side.
     * <p>
     * Cheap -- two calls to {@link #inverseCdf(double)} -- and the usual
     * choice, but for a skewed posterior it is not the shortest interval
     * holding that mass; {@link #highestDensityInterval(double)} is.
     *
     * @param level
     *            the mass the interval is to hold, in {@code (0, 1)}
     * @return the interval
     * @throws IllegalArgumentException
     *             if {@code level} is {@code NaN} or is not in {@code (0, 1)}
     */
    public Interval credibleInterval(double level) {
        checkLevel(level);
        double tail = (1.0 - level) / 2.0;
        double lower = inverseCdf(tail);
        double upper = inverseCdf(1.0 - tail);
        return new Interval(lower, upper, cdf(upper) - cdf(lower));
    }

    /**
     * The highest posterior density interval at the given level: the
     * <b>shortest</b> interval holding that mass, which for a posterior with
     * one peak is the set of points whose density exceeds a common threshold.
     * <p>
     * Found through the condition that defines it -- the density is equal at
     * the two ends -- by moving the lower quantile until it is. For a symmetric
     * posterior it is the equal-tailed interval; for a skewed one it is
     * strictly shorter, and it is the one that does not put its lower end where
     * the density is higher than at points it excludes.
     * <p>
     * <b>This is the expensive member.</b> Every step evaluates
     * {@link #inverseCdf(double)} twice and each of those integrates, so a call
     * costs on the order of a hundred times what an equal-tailed interval
     * costs. <b>It assumes one peak</b>, as the class does throughout; against
     * a multimodal posterior it returns an interval that holds the right mass
     * but need not be the shortest.
     * <p>
     * Where the mode sits at an end of the range the interval is one-sided and
     * is returned as such, since there is nothing on the other side of the peak
     * to balance against.
     *
     * @param level
     *            the mass the interval is to hold, in {@code (0, 1)}
     * @return the interval
     * @throws IllegalArgumentException
     *             if {@code level} is {@code NaN} or is not in {@code (0, 1)}
     */
    public Interval highestDensityInterval(double level) {
        checkLevel(level);
        double atMode = cdf(mode);
        // a peak at an end leaves the density monotone, so the shortest
        // interval starts at that end
        if (atMode <= 1.0 - level) {
            if (atMode == 0.0) {
                double upper = inverseCdf(level);
                return new Interval(lo, upper, cdf(upper));
            }
        }
        if (atMode >= level) {
            if (atMode == 1.0) {
                double lower = inverseCdf(1.0 - level);
                return new Interval(lower, hi, 1.0 - cdf(lower));
            }
        }
        final double mass = level;
        // g is the density at the lower end less the density at the upper end,
        // as a function of how much mass is left below the interval. It falls
        // through zero exactly once for a posterior with one peak
        DFunction g = new DFunction() {
            @Override
            public double apply(double q) {
                return pdf(inverseCdf(q)) - pdf(inverseCdf(q + mass));
            }
        };
        double span = 1.0 - level;
        double a = span * 1.0e-6;
        double b = span * (1.0 - 1.0e-6);
        double ga = g.apply(a);
        double gb = g.apply(b);
        if (ga * gb > 0.0) {
            // the density never balances inside the range, which is what a peak
            // pressed against an end looks like; take the end the density
            // favours
            if (ga < 0.0) {
                double upper = inverseCdf(level);
                return new Interval(lo, upper, cdf(upper) - cdf(lo));
            }
            double lower = inverseCdf(1.0 - level);
            return new Interval(lower, hi, cdf(hi) - cdf(lower));
        }
        double q = RootFinder.brentDekker(a, b, g, Math.max(epsTol, MathConsts.MIN_TOL));
        double lower = inverseCdf(q);
        double upper = inverseCdf(q + level);
        return new Interval(lower, upper, cdf(upper) - cdf(lower));
    }

    // ------------------------------------------------------------------
    // construction helpers
    // ------------------------------------------------------------------

    /** Which moment of the shifted integrand a panel integration is to take. */
    private enum Moment {
        ZEROTH, FIRST, CENTRAL
    }

    /**
     * The integrand, {@code exp(logPosterior(x) - logMax)} times a power of
     * {@code x}. Its maximum is exactly one at the mode for
     * {@link Moment#ZEROTH}, whatever the scale of the problem.
     */
    private DFunction shifted(final Moment moment, final double center) {
        return new DFunction() {
            @Override
            public double apply(double x) {
                if (x < lo || x > hi) {
                    return 0.0;
                }
                double value = logPosterior.apply(x);
                if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
                    return 0.0;
                }
                double density = Math.exp(value - logMax);
                switch (moment) {
                case FIRST:
                    return x * density;
                case CENTRAL:
                    double d = x - center;
                    return d * d * density;
                default:
                    return density;
                }
            }
        };
    }

    /**
     * Integrates one panel. The centered substitution is used only where a
     * limit is infinite, which is where it is the difference between finding
     * the mass and not; a panel with two finite ends involves no substitution
     * at all.
     */
    private double panelIntegral(double a, double b, Moment moment, double center) {
        if (a == b) {
            return 0.0;
        }
        DFunction f = shifted(moment, center);
        if (isFinite(a) && isFinite(b)) {
            return Quadrature.integrate(f, a, b, epsTol);
        }
        // The centered substitution insists the center lie strictly inside the
        // panel, so the mode itself will not do here: an unbounded panel is a
        // tail and the mode is on the far side of its finite end. What the
        // substitution needs to resolve is the shoulder just inside that end,
        // which is where whatever mass this panel holds actually sits
        double tailCenter;
        if (!isFinite(a)) {
            tailCenter = b - leftScale;
        } else {
            tailCenter = a + rightScale;
        }
        return Quadrature.integrate(f, a, b, epsTol, tailCenter);
    }

    /**
     * The distance from the mode over which the log posterior falls by half a
     * nat, in the given direction. Stepped out geometrically until the drop is
     * found and then narrowed by {@link RootFinder#brentDekker}.
     */
    private double halfWidth(double direction, double probe) {
        double target = logMax - 0.5;
        double bound = direction < 0.0 ? lo : hi;
        double reach = isFinite(bound) ? Math.abs(bound - mode) : Double.POSITIVE_INFINITY;
        if (reach == 0.0) {
            // the mode sits on this end, so there is no room on this side; the
            // other side sets the scale and this one is never stepped into
            return Math.max(Math.abs(mode), 1.0) * SQRT_EPS;
        }
        double d = Math.min(probe, reach) / 8.0;
        if (!(d > 0.0)) {
            d = Math.min(1.0, reach);
        }
        double far = 0.0;
        boolean dropped = false;
        for (int step = 0; step < MAX_WIDTH_STEPS; step++) {
            far = Math.min(d, reach);
            double at = mode + direction * far;
            double value = logPosterior.apply(at);
            if (Double.isNaN(value) || value <= target) {
                dropped = true;
                break;
            }
            if (far == reach) {
                // the whole side was walked and the posterior never fell far
                // enough, so this side is bounded and short
                return reach;
            }
            d += d;
        }
        if (!dropped) {
            throw new IllegalArgumentException("the log posterior does not fall away from its peak at " + mode
                    + " towards " + (direction < 0.0 ? lo : hi) + "; it cannot be normalized");
        }
        final double level = target;
        final double sign = direction;
        double root = RootFinder.brentDekker(0.0, far, new DFunction() {
            @Override
            public double apply(double t) {
                double value = logPosterior.apply(mode + sign * t);
                if (Double.isNaN(value)) {
                    return -1.0;
                }
                if (value == Double.NEGATIVE_INFINITY) {
                    return -1.0;
                }
                return value - level;
            }
        }, SQRT_EPS);
        double width = Math.abs(root);
        if (!(width > 0.0) || !isFinite(width)) {
            width = Math.max(Math.abs(mode), 1.0) * SQRT_EPS;
        }
        return width;
    }

    /** The panel edges: the range ends, the mode, and the ladder around it. */
    private double[] buildEdges() {
        double[] tmp = new double[2 * LADDER.length + 3];
        int n = 0;
        tmp[n++] = lo;
        for (int i = LADDER.length - 1; i >= 0; i--) {
            double x = mode - leftScale * LADDER[i];
            if (x > lo && x < mode && x > tmp[n - 1]) {
                tmp[n++] = x;
            }
        }
        if (mode > tmp[n - 1]) {
            tmp[n++] = mode;
        }
        for (int i = 0; i < LADDER.length; i++) {
            double x = mode + rightScale * LADDER[i];
            if (x > tmp[n - 1] && x < hi) {
                tmp[n++] = x;
            }
        }
        if (hi > tmp[n - 1]) {
            tmp[n++] = hi;
        }
        double[] result = new double[n];
        System.arraycopy(tmp, 0, result, 0, n);
        return result;
    }

    /** The panel {@code x} falls in, the largest {@code j} with {@code edges[j] <= x}. */
    private int panelOf(double x) {
        int low = 0;
        int high = edges.length - 1;
        while (low < high - 1) {
            int mid = (low + high) >>> 1;
            if (edges[mid] <= x) {
                low = mid;
            } else {
                high = mid;
            }
        }
        return low;
    }

    /**
     * A finite end for a root search in an unbounded outer panel: steps out
     * from {@code from} until the distribution function has passed {@code p}.
     */
    private double walkOut(double from, double direction, double p) {
        double scale = direction < 0.0 ? leftScale : rightScale;
        double d = Math.max(scale, Math.abs(from) * SQRT_EPS);
        if (!(d > 0.0)) {
            d = 1.0;
        }
        double at = from;
        for (int step = 0; step < MAX_WIDTH_STEPS; step++) {
            at = from + direction * d;
            double c = cdf(at);
            if (direction < 0.0 ? c <= p : c >= p) {
                return at;
            }
            d += d;
        }
        return at;
    }

    private static void checkLevel(double level) {
        if (Double.isNaN(level) || !(level > 0.0) || !(level < 1.0)) {
            throw new IllegalArgumentException("level must be in (0, 1) : " + level);
        }
    }

    /** Where to start looking for the peak when the caller does not say. */
    private static double defaultStart(double lo, double hi) {
        if (isFinite(lo) && isFinite(hi)) {
            return lo + (hi - lo) / 2.0;
        }
        if (isFinite(lo)) {
            return lo + 1.0;
        }
        if (isFinite(hi)) {
            return hi - 1.0;
        }
        return 0.0;
    }

    /** The first step of the downhill walk, small enough to stay in the range. */
    private static double initialStep(double lo, double hi, double start) {
        double span = hi - lo;
        double step = isFinite(span) ? span / 100.0 : 1.0;
        if (!(step > 0.0)) {
            step = 1.0;
        }
        if (start + step > hi) {
            step = -step;
        }
        return step;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
