package math.stats.bayes;

import math.MathConsts;
import math.distribution.MultivariateNormal;
import math.distribution.Normal;
import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.optim.NelderMead;

/**
 * The Laplace approximation of a posterior: the Gaussian that matches it at its
 * mode.
 * <p>
 * A caller hands over the <b>logarithm</b> of an unnormalized posterior of any
 * number of parameters, and gets back the mode, the curvature there, the
 * approximate log evidence, and the whole thing as a
 * {@link MultivariateNormal}. It is the cheap answer to the question
 * {@link ScalarPosterior} and {@link BivariatePosterior} answer exactly:
 * measured, three parameters cost about 520 evaluations of the caller's
 * function here against four to seven million for one quadrature integral, and
 * a full quadrature construction is ten of those.
 * <p>
 * <b>What it costs and how far it reaches.</b> The mode comes from
 * {@link NelderMead}, which is derivative free, and the curvature from central
 * differences, which is {@code 2 d^2 + 1} evaluations. Measured against an
 * equicorrelated normal in {@code d} dimensions, where every part has a closed
 * form:
 *
 * <pre>
 *   d    correlation   worst log evidence (exact answer 0)   worst mode error
 *   3    0 to 0.9      -7.4e-08                              1.8e-08
 *   8    0 to 0.9      -2.5e-07                              1.8e-08
 *   20   0 to 0.9      -1.8e-06                              3.1e-08
 *   40   0 to 0.9      -6.6e-06                              3.7e-08
 *   60   0 to 0.9      -1.9e-05                              4.2e-08
 * </pre>
 *
 * The evidence error grows with the dimension and the mode error does not,
 * because the two are not the same quantity: the log determinant is a sum of
 * {@code d} logarithms, each carrying the accuracy of a difference quotient, so
 * it accumulates where an abscissa does not. At {@code d = 60} a construction
 * is a few hundred milliseconds and at {@code d = 20} about ten.
 * <p>
 * <b>When it holds.</b> Exactly for a Gaussian posterior, and to the extent the
 * posterior is Gaussian otherwise -- which for a posterior over {@code n}
 * observations improves as {@code n} grows. Measured against the exact answer
 * from {@link BivariatePosterior} on the posterior of a mean and a log scale,
 * the error in the log evidence runs {@code 9.2e-03} at a hundred observations,
 * {@code 2.3e-03} at four hundred and {@code 4.7e-04} at two thousand. For one
 * and two parameters the exact answer is available and there is no reason to
 * approximate; this class is for three and above, where it is the only answer
 * that is affordable.
 * <p>
 * Instances are immutable and can be shared between threads.
 *
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Laplace%27s_approximation">Wikipedia
 *      Laplace's approximation</a>
 * @since 1.5.3
 */
public final class LaplaceApproximation {

    /**
     * The fourth root of the machine epsilon, about {@code 1.2e-4}: the step of
     * a central second difference. Smaller and the cancellation of the three
     * function values dominates, larger and the difference stops being local.
     */
    private static final double EPS_FOURTH = Math.pow(MathConsts.MACH_EPS_DBL, 0.25);

    /** How far the search for an axis scale may double before giving up. */
    private static final int MAX_SCALE_STEPS = 200;

    private final DMultiFunction logPosterior;
    private final double[] lo;
    private final double[] hi;

    private final double[] mode;
    private final double logMax;
    private final boolean modeConverged;
    private final DMatrix curvature;
    private final double logEvidence;
    private final MultivariateNormal normal;
    private final double massOutsideSupport;

    /**
     * The Laplace approximation of {@code logPosterior} over the box
     * {@code [lo, hi]}, starting the search for the peak at the middle of
     * whatever the box bounds.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density. It
     *            is never evaluated outside the box and may return
     *            {@link Double#NEGATIVE_INFINITY} where the posterior is zero
     * @param lo
     *            the lower end of each axis, entries possibly
     *            {@link Double#NEGATIVE_INFINITY}. Not modified
     * @param hi
     *            the upper end of each axis, entries possibly
     *            {@link Double#POSITIVE_INFINITY}. Not modified
     * @return the approximation
     * @throws IllegalArgumentException
     *             if any argument is {@code null}, if {@code lo} and {@code hi}
     *             are not of the same non-zero length, if {@code lo[i]} is not
     *             below {@code hi[i]}, if either is {@code NaN}, if no peak can
     *             be found, or if the peak is not a maximum
     */
    public static LaplaceApproximation of(DMultiFunction logPosterior, double[] lo, double[] hi) {
        checkBox(lo, hi);
        return new LaplaceApproximation(logPosterior, lo, hi, defaultStart(lo, hi));
    }

    /**
     * The Laplace approximation of {@code logPosterior} over the box
     * {@code [lo, hi]}, starting the search for the peak at {@code start}.
     *
     * @param logPosterior
     *            the natural logarithm of an unnormalized posterior density
     * @param lo
     *            the lower end of each axis. Not modified
     * @param hi
     *            the upper end of each axis. Not modified
     * @param start
     *            finite coordinates inside the box, near the mass. Not modified
     * @return the approximation
     * @throws IllegalArgumentException
     *             under the conditions
     *             {@link #of(DMultiFunction, double[], double[])} states, and
     *             additionally if {@code start} is {@code null} or does not
     *             hold finite coordinates inside the box
     */
    public static LaplaceApproximation of(DMultiFunction logPosterior, double[] lo, double[] hi, double[] start) {
        checkBox(lo, hi);
        return new LaplaceApproximation(logPosterior, lo, hi, start);
    }

    private LaplaceApproximation(DMultiFunction logPosterior, double[] lo, double[] hi, double[] start) {
        if (logPosterior == null) {
            throw new IllegalArgumentException("logPosterior must not be null");
        }
        int d = lo.length;
        if (start == null || start.length != d) {
            throw new IllegalArgumentException("start must hold " + d + " coordinates");
        }
        for (int i = 0; i < d; i++) {
            if (!isFinite(start[i]) || start[i] < lo[i] || start[i] > hi[i]) {
                throw new IllegalArgumentException("start[" + i + "] must be finite and inside [" + lo[i] + ", "
                        + hi[i] + "] : " + start[i]);
            }
        }
        this.logPosterior = logPosterior;
        this.lo = lo.clone();
        this.hi = hi.clone();

        // ---- the peak -------------------------------------------------
        // walled off outside the box with a *finite* value: a simplex step
        // forms differences of function values, and Infinity - Infinity is
        // NaN, after which no peak is ever bracketed
        final double reference = logPosterior.apply(start.clone());
        if (!isFinite(reference)) {
            throw new IllegalArgumentException("the log posterior is not finite at the starting point : "
                    + reference + "; name a start where the posterior is positive");
        }
        final double wall = -reference + 1.0;
        DMultiFunction fenced = new DMultiFunction() {
            @Override
            public double apply(double[] p) {
                double outside = 0.0;
                for (int i = 0; i < p.length; i++) {
                    if (p[i] < LaplaceApproximation.this.lo[i]) {
                        outside += LaplaceApproximation.this.lo[i] - p[i];
                    } else if (p[i] > LaplaceApproximation.this.hi[i]) {
                        outside += p[i] - LaplaceApproximation.this.hi[i];
                    }
                }
                if (outside > 0.0) {
                    return wall + outside;
                }
                double value = LaplaceApproximation.this.logPosterior.apply(p);
                if (Double.isNaN(value) || value == Double.NEGATIVE_INFINITY) {
                    return wall;
                }
                return -value;
            }
        };
        DMultiFunctionEval peak = new NelderMead().minimize(fenced, start.clone(), initialStep(lo, hi, start));
        double[] at = peak.point.clone();
        for (int i = 0; i < d; i++) {
            if (at[i] < lo[i]) {
                at[i] = lo[i];
            }
            if (at[i] > hi[i]) {
                at[i] = hi[i];
            }
        }
        this.mode = at;
        this.modeConverged = peak.converged;
        this.logMax = logPosterior.apply(at.clone());
        if (!isFinite(logMax)) {
            throw new IllegalArgumentException("the log posterior is not finite at its peak : " + logMax);
        }

        // ---- the curvature --------------------------------------------
        // the step of a second difference wants the width of the peak and not
        // the magnitude of its position: a posterior over a few thousand
        // observations is a needle whose width has nothing to do with where it
        // sits, and a step taken from |mode| samples where the density is zero
        double[] scale = new double[d];
        for (int i = 0; i < d; i++) {
            scale[i] = axisScale(i);
        }
        this.curvature = hessian(scale);

        // ---- the factor, and everything that comes out of it ----------
        DMatrix factor;
        try {
            factor = CholeskyDecomp.cholesky(curvature);
        } catch (RuntimeException e) {
            throw new IllegalArgumentException(
                    "the curvature at the peak is not positive definite, so what was found is not a maximum: "
                            + curvature);
        }
        for (int i = 0; i < d; i++) {
            if (!(factor.get(i, i) > 0.0)) {
                throw new IllegalArgumentException("the curvature at the peak could not be factored: " + curvature);
            }
        }
        this.logEvidence = logMax + 0.5 * d * Math.log(2.0 * Math.PI)
                - 0.5 * CholeskyDecomp.logDeterminant(factor);

        // the covariance is the inverse of the curvature, column by column
        // through the factor rather than by forming an inverse: that is what
        // CholeskyDecomp.solve is for and it does not square the conditioning
        DMatrix covariance = new DMatrix(d, d);
        double[] unit = new double[d];
        double[] column = new double[d];
        for (int j = 0; j < d; j++) {
            java.util.Arrays.fill(unit, 0.0);
            unit[j] = 1.0;
            CholeskyDecomp.solve(factor, unit, column);
            for (int i = 0; i < d; i++) {
                covariance.set(i, j, column[i]);
            }
        }
        // an inverse computed column by column is symmetric only to round-off,
        // and MultivariateNormal compares relatively -- but averaging here is
        // free and makes what it is handed exact
        for (int i = 0; i < d; i++) {
            for (int j = i + 1; j < d; j++) {
                double m = 0.5 * (covariance.get(i, j) + covariance.get(j, i));
                covariance.set(i, j, m);
                covariance.set(j, i, m);
            }
        }
        this.normal = new MultivariateNormal(mode, covariance);

        // ---- how much of it falls outside the box ---------------------
        double leak = 0.0;
        for (int i = 0; i < d; i++) {
            double sd = Math.sqrt(covariance.get(i, i));
            if (!(sd > 0.0) || !isFinite(sd)) {
                continue;
            }
            Normal axis = new Normal(mode[i], sd);
            if (isFinite(lo[i])) {
                leak += axis.cdf(lo[i]);
            }
            if (isFinite(hi[i])) {
                leak += 1.0 - axis.cdf(hi[i]);
            }
        }
        this.massOutsideSupport = Math.min(1.0, leak);
    }

    /**
     * The number of parameters.
     *
     * @return the dimension, one or more
     */
    public int dimension() {
        return mode.length;
    }

    /**
     * Writes the mode, where the posterior peaks, into {@code out}.
     *
     * @param out
     *            where the result is written, one entry per parameter. Its
     *            previous contents are overwritten
     * @throws IllegalArgumentException
     *             if {@code out} is {@code null} or is not as long as there are
     *             parameters
     */
    public void mode(double[] out) {
        if (out == null || out.length != mode.length) {
            throw new IllegalArgumentException("out must hold " + mode.length + " coordinates");
        }
        System.arraycopy(mode, 0, out, 0, mode.length);
    }

    /**
     * The value of the log posterior at the mode, before any normalization.
     *
     * @return {@code logPosterior(mode)}
     */
    public double logMax() {
        return logMax;
    }

    /**
     * Whether the mode search reported convergence.
     * <p>
     * <b>It is a measurement and not a verdict, which is why it is a flag
     * rather than an exception.</b> Handed {@link NelderMead}'s own default
     * initial simplex, forty parameters at a correlation of {@code 0.9} run out
     * of the iteration budget after 368 000 evaluations and report
     * {@code false} -- with a mode that is nonetheless right to {@code 9.3e-09}.
     * <b>The size of that first simplex is what decides it, not the
     * dimension</b>: with the step this class sets instead -- a hundredth of a
     * bounded span, one unit where the axis is unbounded -- the flag came back
     * {@code true} in every case measured, up to sixty parameters, correlations
     * to {@code 0.9} and starting points fifty units from the mass.
     * <p>
     * So a {@code false} here is a reason to check the answer rather than to
     * discard it, and naming a {@code start} nearer the mass is what usually
     * clears it.
     *
     * @return whether the simplex search converged
     */
    public boolean modeConverged() {
        return modeConverged;
    }

    /**
     * The curvature at the mode: the Hessian of the <em>negative</em> log
     * posterior, by central differences, exactly symmetric.
     * <p>
     * Its inverse is the covariance of {@link #asNormal()}.
     *
     * @return a fresh curvature matrix
     */
    public DMatrix curvature() {
        return curvature.copy();
    }

    /**
     * The approximate log marginal likelihood,
     * {@code logMax + (d/2) log(2 pi) - 0.5 log det(H)}.
     * <p>
     * Exact when the posterior is Gaussian, and otherwise as good as the
     * posterior is Gaussian; the class comment carries the measured figures.
     * The ratio of two of these is a Bayes factor.
     *
     * @return the approximate natural logarithm of the evidence
     */
    public double logEvidence() {
        return logEvidence;
    }

    /**
     * The approximation itself: the normal centered at the mode whose
     * covariance is the inverse curvature there.
     *
     * @return the Gaussian this approximation consists of
     */
    public MultivariateNormal asNormal() {
        return normal;
    }

    /**
     * How much of {@link #asNormal()} falls outside the box the posterior was
     * defined on.
     * <p>
     * A Gaussian has all of {@code R^d} for its support, so an approximation to
     * a posterior over a positive parameter always puts some mass where the
     * posterior has none. This reports how much rather than keeping quiet about
     * it: near zero means the bound does not matter, and a value that is not
     * small means the approximation is describing a shape the posterior does
     * not have.
     * <p>
     * <b>It is a union bound.</b> The tails are summed axis by axis, so a point
     * outside on two axes at once is counted twice; the number is therefore an
     * upper bound, exact in one dimension and conservative above it. It is
     * clamped at one.
     *
     * @return an upper bound on the mass outside the box, in {@code [0, 1]}
     */
    public double massOutsideSupport() {
        return massOutsideSupport;
    }

    // ------------------------------------------------------------------
    // construction helpers
    // ------------------------------------------------------------------

    /**
     * The width of the peak along one axis: the distance over which the log
     * posterior falls half a nat, found by doubling outward. Used only to scale
     * a difference quotient, so a factor of two either way is harmless.
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
        double step = Math.min(1.0, reach) / 64.0;
        if (!(step > 0.0)) {
            step = 1.0;
        }
        double[] p = new double[mode.length];
        for (int k = 0; k < MAX_SCALE_STEPS; k++) {
            double far = Math.min(step, reach);
            System.arraycopy(mode, 0, p, 0, mode.length);
            p[axis] = mode[axis] + direction * far;
            double value = logPosterior.apply(p);
            if (Double.isNaN(value) || value <= target) {
                return far;
            }
            if (far == reach) {
                return reach;
            }
            step += step;
        }
        throw new IllegalArgumentException("the log posterior does not fall away from its peak along axis " + axis
                + "; it cannot be approximated");
    }

    /**
     * The Hessian of the negative log posterior at the mode.
     * <p>
     * Each off-diagonal is computed once and written into both places. That is
     * not tidiness: {@link CholeskyDecomp} tests symmetry with an
     * <em>absolute</em> tolerance, so a curvature whose entries are of order
     * {@code 1e6} -- an ordinary posterior over a few thousand observations --
     * would be rejected over its own round-off if the two halves were computed
     * separately.
     */
    private DMatrix hessian(double[] scale) {
        int d = mode.length;
        double[] h = new double[d];
        for (int i = 0; i < d; i++) {
            h[i] = EPS_FOURTH * scale[i];
        }
        double f0 = -logMax;
        DMatrix result = new DMatrix(d, d);
        double[] delta = new double[d];
        for (int i = 0; i < d; i++) {
            java.util.Arrays.fill(delta, 0.0);
            delta[i] = h[i];
            double plus = negAt(delta);
            delta[i] = -h[i];
            double minus = negAt(delta);
            result.set(i, i, (plus - 2.0 * f0 + minus) / (h[i] * h[i]));
            for (int j = i + 1; j < d; j++) {
                java.util.Arrays.fill(delta, 0.0);
                delta[i] = h[i];
                delta[j] = h[j];
                double pp = negAt(delta);
                delta[j] = -h[j];
                double pm = negAt(delta);
                delta[i] = -h[i];
                double mm = negAt(delta);
                delta[j] = h[j];
                double mp = negAt(delta);
                double value = (pp - pm - mp + mm) / (4.0 * h[i] * h[j]);
                result.set(i, j, value);
                result.set(j, i, value);
            }
        }
        return result;
    }

    /** The negative log posterior a step away from the mode, clamped into the box. */
    private double negAt(double[] delta) {
        double[] p = new double[mode.length];
        for (int i = 0; i < p.length; i++) {
            p[i] = mode[i] + delta[i];
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

    private static void checkBox(double[] lo, double[] hi) {
        if (lo == null || hi == null) {
            throw new IllegalArgumentException("lo and hi must not be null");
        }
        if (lo.length == 0 || lo.length != hi.length) {
            throw new IllegalArgumentException("lo and hi must hold the same non-zero number of entries");
        }
        for (int i = 0; i < lo.length; i++) {
            if (Double.isNaN(lo[i]) || Double.isNaN(hi[i])) {
                throw new IllegalArgumentException("the range of axis " + i + " must not be NaN");
            }
            if (!(lo[i] < hi[i])) {
                throw new IllegalArgumentException(
                        "lo[" + i + "] must be below hi[" + i + "] : [" + lo[i] + ", " + hi[i] + "]");
            }
        }
    }

    /** Where to start looking for the peak when the caller does not say. */
    private static double[] defaultStart(double[] lo, double[] hi) {
        double[] start = new double[lo.length];
        for (int i = 0; i < start.length; i++) {
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
        double[] step = new double[lo.length];
        for (int i = 0; i < step.length; i++) {
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
