package math.ode;

/**
 * Decides how long the next step should be, from the error the last one
 * estimated for itself.
 * <p>
 * The estimate is a vector, and the first thing that happens to it is that it
 * becomes a number: each component is divided by what an acceptable error at
 * that component would be, <code>atol + rtol max(|y|, |yNext|)</code>, and the
 * root mean square of those is taken. A step is kept when that number is at
 * most one, so the two tolerances are not a bound on the error of the answer
 * but on the error one step is allowed to add -- the difference being what
 * accumulates over the run, and why halving {@code rtol} does not halve the
 * error at the end.
 * <p>
 * <b>The absolute tolerance is what keeps a component that passes through zero
 * from asking for an infinitely small step</b>, and it is therefore not
 * optional and not a formality: it carries the scale below which a component
 * stops mattering, which only the caller knows.
 * <p>
 * The new step follows <code>safety err^-a errPrevious^b</code>, clamped so
 * that no single step can grow or shrink without limit. With {@code b} at zero
 * that is the textbook controller, which reacts to the last step alone and
 * tends to oscillate -- a step that was slightly too long is followed by one
 * that is too short, and the sequence rings instead of settling. The
 * proportional-integral form, {@code b} at {@code 0.04} as Hairer's
 * {@code dopri5} sets it, damps that by remembering the step before.
 * <p>
 * <b>The damped form is the default, and it is not free.</b> Measured against
 * the textbook one, it costs between two and eight percent more evaluations on
 * every smooth problem tried -- the harmonic oscillator, Kepler orbits at
 * eccentricities from {@code 0.6} to {@code 0.99}, Lotka-Volterra, the
 * Brusselator -- because on those the step size moves slowly and there is
 * nothing to damp. What it buys shows up where the step size has to move fast:
 * on van der Pol at {@code mu = 1000} it takes the rejected steps from
 * {@code 42444} down to {@code 72}, and saves evaluations doing it. A bounded
 * loss on the easy problems against an unbounded gain on the hard ones is the
 * right way round for a default, but a caller who knows the equation is smooth
 * can set {@code beta} to zero and keep the few percent.
 * <p>
 * <b>The exponent is one over the order of the method</b>, because the estimate
 * measures the error of the <em>lower</em> half of the embedded pair and that
 * one falls with the order of the method rather than one below it. It is the
 * exponent Hairer uses for both of his pairs.
 * <p>
 * Instances are immutable and can be shared between threads and between
 * integrators.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Adaptive_step_size">Wikipedia adaptive
 * step size</a>.
 *
 * @since 1.5.3
 */
public final class StepController {

    /** What Hairer starts the remembered error at, and its floor thereafter. */
    static final double ERROR_FLOOR = 1.0e-4;

    private final double rtol;
    private final double atol;
    private final double safety;
    private final double minScale;
    private final double maxScale;
    private final double beta;
    private final double maxStepSize;
    private final int maxSteps;

    /**
     * A controller with the settings that are right until measurement says
     * otherwise: safety {@code 0.9}, a step that may shrink to a fifth or grow
     * tenfold, the proportional-integral term at {@code 0.04}, no bound on the
     * step size, and at most {@code 100000} steps.
     *
     * @param rtol
     *            the relative tolerance, the share of the state a step may err
     *            by, typically between {@code 1e-3} and {@code 1e-12}
     * @param atol
     *            the absolute tolerance, the size below which a component is
     *            not worth resolving; it must be positive
     * @throws IllegalArgumentException
     *             if a tolerance is not positive or not finite
     */
    public StepController(double rtol, double atol) {
        this(rtol, atol, 0.9, 0.2, 10.0, 0.04, Double.POSITIVE_INFINITY, 100000);
    }

    /**
     * A controller with the default settings but a step budget of its own.
     *
     * @param rtol
     *            the relative tolerance
     * @param atol
     *            the absolute tolerance
     * @param maxSteps
     *            how many steps, accepted and rejected together, a run may take
     *            before it gives up
     * @throws IllegalArgumentException
     *             if a tolerance is not positive or {@code maxSteps} is not
     */
    public StepController(double rtol, double atol, int maxSteps) {
        this(rtol, atol, 0.9, 0.2, 10.0, 0.04, Double.POSITIVE_INFINITY, maxSteps);
    }

    /**
     * A controller in full.
     *
     * @param rtol
     *            the relative tolerance
     * @param atol
     *            the absolute tolerance
     * @param safety
     *            the factor the proposed step is shortened by, between zero and
     *            one, so that a step aimed exactly at the tolerance is not
     *            rejected half the time
     * @param minScale
     *            the most one step may shrink by, between zero and one
     * @param maxScale
     *            the most one step may grow by, greater than one
     * @param beta
     *            the weight of the previous error, zero for the textbook
     *            controller and around {@code 0.04} for the damped one; it must
     *            be at least zero and less than one over the order
     * @param maxStepSize
     *            an upper bound on the length of a step, or
     *            {@link Double#POSITIVE_INFINITY} for none
     * @param maxSteps
     *            how many steps a run may take before it gives up
     * @throws IllegalArgumentException
     *             if a setting is outside the range described
     */
    public StepController(double rtol, double atol, double safety, double minScale, double maxScale,
            double beta, double maxStepSize, int maxSteps) {
        if (!(rtol > 0.0) || Double.isInfinite(rtol)) {
            throw new IllegalArgumentException("rtol must be positive and finite, not " + rtol);
        }
        if (!(atol > 0.0) || Double.isInfinite(atol)) {
            throw new IllegalArgumentException("atol must be positive and finite, not " + atol);
        }
        if (!(safety > 0.0) || !(safety <= 1.0)) {
            throw new IllegalArgumentException("safety must be in (0, 1], not " + safety);
        }
        if (!(minScale > 0.0) || !(minScale < 1.0)) {
            throw new IllegalArgumentException("minScale must be in (0, 1), not " + minScale);
        }
        if (!(maxScale > 1.0)) {
            throw new IllegalArgumentException("maxScale must be greater than 1, not " + maxScale);
        }
        if (!(beta >= 0.0) || !(beta < 0.5)) {
            throw new IllegalArgumentException("beta must be in [0, 0.5), not " + beta);
        }
        if (!(maxStepSize > 0.0)) {
            throw new IllegalArgumentException("maxStepSize must be positive, not " + maxStepSize);
        }
        if (maxSteps < 1) {
            throw new IllegalArgumentException("maxSteps must be at least 1, not " + maxSteps);
        }
        this.rtol = rtol;
        this.atol = atol;
        this.safety = safety;
        this.minScale = minScale;
        this.maxScale = maxScale;
        this.beta = beta;
        this.maxStepSize = maxStepSize;
        this.maxSteps = maxSteps;
    }

    /**
     * The relative tolerance.
     *
     * @return the share of the state one step may err by
     */
    public double relativeTolerance() {
        return rtol;
    }

    /**
     * The absolute tolerance.
     *
     * @return the size below which a component is not resolved
     */
    public double absoluteTolerance() {
        return atol;
    }

    /**
     * The factor a proposed step is shortened by.
     *
     * @return the safety factor
     */
    public double safety() {
        return safety;
    }

    /**
     * The most one step may shrink by.
     *
     * @return the smallest factor a step size is multiplied by
     */
    public double minScale() {
        return minScale;
    }

    /**
     * The most one step may grow by.
     *
     * @return the largest factor a step size is multiplied by
     */
    public double maxScale() {
        return maxScale;
    }

    /**
     * The weight of the previous error, zero for the textbook controller.
     *
     * @return the proportional-integral term
     */
    public double beta() {
        return beta;
    }

    /**
     * The upper bound on the length of a step.
     *
     * @return the bound, or {@link Double#POSITIVE_INFINITY} if there is none
     */
    public double maxStepSize() {
        return maxStepSize;
    }

    /**
     * The step budget of a run.
     *
     * @return how many steps, accepted and rejected together, may be taken
     */
    public int maxSteps() {
        return maxSteps;
    }

    /**
     * Turns a componentwise error estimate into the number that decides the
     * step: the root mean square of the components, each measured against what
     * an acceptable error there would be.
     *
     * @param error
     *            the estimate the step produced
     * @param y
     *            the state the step started from
     * @param yNext
     *            the state it reached
     * @return the scaled error, a step being acceptable at one or below
     */
    public double errorNorm(double[] error, double[] y, double[] yNext) {
        double sum = 0.0;
        for (int i = 0; i < error.length; ++i) {
            double scale = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(yNext[i]));
            double e = error[i] / scale;
            sum += e * e;
        }
        return Math.sqrt(sum / error.length);
    }

    /**
     * The same, for a method that produced two estimates of different orders.
     * <p>
     * <b>Why they are combined rather than one of them used.</b> Take a method
     * of order {@code q} carrying embedded solutions of orders {@code p1} and
     * {@code p2}, and write {@code s1} and {@code s2} for the sums of the
     * squared scaled components of the two estimates. Then
     * <p>
     * <code>err = s1 sqrt(1 / (n (s1 + eps s2)))</code>
     * <p>
     * with {@code eps} of {@code 0.01}. This is Hairer's expression, which
     * carries a leading <code>|h|</code> because there the two estimates are
     * combinations of derivatives; the ones arriving here are already
     * multiplied by the step size, and that factor comes out of the ratio by
     * itself. Where the first estimate dominates this
     * is the plain root mean square of it. Where the second does -- which is
     * what happens as the step shrinks, since the lower order estimate goes to
     * zero more slowly -- it becomes {@code |h| s1 / sqrt(eps n s2)}, and with
     * <code>s1 ~ h^(2 p1 + 2)</code> and <code>s2 ~ h^(2 p2 + 2)</code> that
     * falls like <code>h^(2 p1 - p2 + 1)</code>. For DOP853, where the orders
     * are 5 and 3, the exponent is 8: the number falls at the order of the
     * solution the method actually advances, and controlling the step on it
     * controls the accuracy that is kept rather than the accuracy that is
     * thrown away.
     * <p>
     * The formula is Hairer's, and so is the guard: a denominator that is not
     * positive -- both estimates exactly zero -- is replaced by one, which
     * leaves {@code err} at zero and the step accepted.
     *
     * @param error
     *            the estimate of the higher order, from
     *            {@link OdeStepper#step(double, double[], double, double[], double[])}
     * @param secondary
     *            the estimate of the lower order, from
     *            {@link OdeStepper#secondaryError()}
     * @param y
     *            the state the step started from
     * @param yNext
     *            the state it reached
     * @return the scaled error, a step being acceptable at one or below
     */
    public double errorNorm(double[] error, double[] secondary, double[] y, double[] yNext) {
        double sum = 0.0;
        double sumSecondary = 0.0;
        for (int i = 0; i < error.length; ++i) {
            double scale = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(yNext[i]));
            double e = error[i] / scale;
            double e2 = secondary[i] / scale;
            sum += e * e;
            sumSecondary += e2 * e2;
        }
        double denominator = sum + 0.01 * sumSecondary;
        if (!(denominator > 0.0)) {
            denominator = 1.0;
        }
        return sum * Math.sqrt(1.0 / (error.length * denominator));
    }

    /**
     * The factor to multiply the step size by after a step that was kept.
     *
     * @param error
     *            the scaled error of the step that was kept, at most one
     * @param previousError
     *            the scaled error of the step before it, which damps the
     *            response; {@link #ERROR_FLOOR} stands in before there is one
     * @param order
     *            the order of the method
     * @return a factor between {@link #minScale()} and {@link #maxScale()}
     */
    public double scale(double error, double previousError, int order) {
        if (error <= 0.0) {
            return maxScale;
        }
        double exponent = 1.0 / order - 0.75 * beta;
        double factor = safety * Math.pow(error, -exponent) * Math.pow(previousError, beta);
        return clamp(factor);
    }

    /**
     * The factor to multiply the step size by after a step that was thrown
     * away. It never grows the step and it does not look at the step before,
     * whose error says nothing about a proposal that has already failed.
     *
     * @param error
     *            the scaled error of the rejected step, greater than one, or
     *            not a number if the step produced a state that was not finite
     * @param order
     *            the order of the method
     * @return a factor between {@link #minScale()} and one
     */
    public double scaleAfterRejection(double error, int order) {
        if (Double.isNaN(error) || Double.isInfinite(error)) {
            return minScale;
        }
        double factor = safety * Math.pow(error, -1.0 / order);
        return Math.max(minScale, Math.min(1.0, factor));
    }

    /**
     * A first step size, from two evaluations of the field and nothing else.
     * <p>
     * The state and its rate of change give a length over which the solution
     * changes by about a hundredth of itself; a trial step of that length gives
     * the curvature, and the step that would meet the tolerance at that
     * curvature is the answer. This is the estimate from Hairer, Norsett and
     * Wanner, and it is worth the two evaluations: a first step guessed an
     * order of magnitude too long is rejected repeatedly before the controller
     * finds the scale, and one guessed too short is paid for over the whole
     * run.
     *
     * @param stepper
     *            the method, asked twice for the rate of change
     * @param t0
     *            the time to start from
     * @param y0
     *            the state to start from
     * @param t1
     *            the time to reach, which fixes the direction and bounds the
     *            first step
     * @return a signed step size, pointing from {@code t0} towards {@code t1}
     * @throws IllegalArgumentException
     *             if the arguments are out of shape or the interval is empty
     */
    public double initialStep(OdeStepper stepper, double t0, double[] y0, double t1) {
        if (stepper == null) {
            throw new IllegalArgumentException("stepper must not be null");
        }
        int n = stepper.dimension();
        if (y0 == null || y0.length != n) {
            throw new IllegalArgumentException("y0 must be of length " + n);
        }
        if (t1 == t0) {
            throw new IllegalArgumentException("t1 must differ from t0, both being " + t0);
        }
        double span = Math.abs(t1 - t0);
        double direction = (t1 > t0) ? 1.0 : -1.0;

        double[] f0 = new double[n];
        double[] trial = new double[n];
        double[] f1 = new double[n];
        stepper.derivative(t0, y0, f0);

        double d0 = 0.0;
        double d1 = 0.0;
        for (int i = 0; i < n; ++i) {
            double scale = atol + rtol * Math.abs(y0[i]);
            double a = y0[i] / scale;
            double b = f0[i] / scale;
            d0 += a * a;
            d1 += b * b;
        }
        d0 = Math.sqrt(d0 / n);
        d1 = Math.sqrt(d1 / n);

        double h0 = (d0 < 1.0e-5 || d1 < 1.0e-5) ? 1.0e-6 : 0.01 * (d0 / d1);
        h0 = Math.min(h0, span);
        for (int i = 0; i < n; ++i) {
            trial[i] = y0[i] + direction * h0 * f0[i];
        }
        stepper.derivative(t0 + direction * h0, trial, f1);

        double d2 = 0.0;
        for (int i = 0; i < n; ++i) {
            double scale = atol + rtol * Math.abs(y0[i]);
            double e = (f1[i] - f0[i]) / scale;
            d2 += e * e;
        }
        d2 = Math.sqrt(d2 / n) / h0;

        double worst = Math.max(d1, d2);
        double h1 = (worst <= 1.0e-15) ? Math.max(1.0e-6, h0 * 1.0e-3)
                : Math.pow(0.01 / worst, 1.0 / stepper.order());
        double h = Math.min(100.0 * h0, Math.min(h1, span));
        h = Math.min(h, maxStepSize);
        return direction * h;
    }

    /**
     * The settings, in the order the full constructor takes them.
     */
    @Override
    public String toString() {
        return "StepController[rtol " + rtol + ", atol " + atol + ", safety " + safety + ", scale ["
                + minScale + ", " + maxScale + "], beta " + beta + ", maxStepSize " + maxStepSize
                + ", maxSteps " + maxSteps + "]";
    }

    private double clamp(double factor) {
        if (factor < minScale) {
            return minScale;
        }
        if (factor > maxScale) {
            return maxScale;
        }
        return factor;
    }
}
