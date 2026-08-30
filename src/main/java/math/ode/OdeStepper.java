package math.ode;

/**
 * One step of a numerical method for an initial value problem, and what a
 * driver needs to know about it.
 * <p>
 * <b>A stepper owns its field.</b> It is constructed with the equation it
 * solves rather than handed one per call, which is what lets a method for a
 * second order system -- a Runge-Kutta-Nystroem or a symplectic one, working on
 * {@link math.fun.DSecondOrderField} with the position and the velocity stacked
 * into one state -- reach a driver through the same interface as a method for
 * {@link math.fun.DVectorField}. The driver never learns which of the two it is
 * advancing.
 * <p>
 * <b>Instances are stateful.</b> A step leaves behind what
 * {@link #interpolate(double, double[])} needs, and an implementation may keep
 * a stage across calls when the point it was evaluated at has not moved. Two
 * threads cannot share one.
 * <p>
 * <b>A step does not judge itself.</b> It fills in the new state and, if it can,
 * an estimate of its own error, and that is all: whether the step is accepted,
 * how the estimate becomes a scalar and what happens to a state that came back
 * non-finite are decisions of {@link StepController} and {@link OdeIntegrator}.
 * The same stepper therefore serves the fixed step and the adaptive path
 * without knowing which one is driving it.
 *
 * @since 1.5.3
 */
public interface OdeStepper {

    /**
     * The number of components of the state this stepper advances.
     *
     * @return {@code n}, one or more
     */
    int dimension();

    /**
     * The order of the method: over a fixed interval, halving the step size
     * divides the accumulated error by two to this power.
     *
     * @return the order, at least one
     */
    int order();

    /**
     * Whether {@link #step(double, double[], double, double[], double[])} fills
     * in the error estimate it is offered. A method without one can only be
     * driven at a fixed step size.
     *
     * @return {@code true} if an error estimate is produced
     */
    boolean hasErrorEstimate();

    /**
     * A second, independent estimate of the error of the step just taken, of a
     * lower order than the first, or {@code null} for a method that carries
     * only one.
     * <p>
     * <b>Why a method would want two.</b> An embedded solution of order
     * {@code p} makes the difference of the pair fall like
     * <code>h^(p+1)</code>, so a step controlled on it is a step controlled on
     * the accuracy of the solution that is discarded. For a pair whose two
     * orders are close that is near enough; for a method of order eight
     * carrying a fifth order estimate it is not, and the step would be held to
     * what the fifth order solution can do. Two estimates of different orders
     * can be combined into a number that falls at the order of the solution the
     * method actually advances, and
     * {@link StepController#errorNorm(double[], double[], double[], double[])}
     * is where that combining happens, because it takes the tolerances to do it.
     * <p>
     * The array belongs to the stepper and is valid until the next step. A
     * caller that keeps it must copy it.
     *
     * @return the second estimate, or {@code null}
     */
    default double[] secondaryError() {
        return null;
    }

    /**
     * Whether the method carries a continuous extension of its own: an
     * interpolant built over stages it has already computed, and one order
     * below the step it spans rather than three.
     * {@link #interpolate(double, double[])} answers either way, so a
     * {@code false} here is a statement about the accuracy and the cost of an
     * interior value, not about whether one can be had.
     *
     * @return {@code true} if the interpolant is the method's own
     */
    boolean hasDenseOutput();

    /**
     * The rate of change of the state at {@code (t, y)}, which for a method on
     * a second order system is the velocity stacked on the acceleration.
     * <p>
     * A driver needs this to guess a first step size before it has taken one,
     * and it is the only thing about the equation a stepper lets out. The call
     * counts towards {@link #evaluations()} and leaves nothing behind: it does
     * not disturb a step in progress and cannot be interpolated over.
     *
     * @param t
     *            the time at which to evaluate
     * @param y
     *            a <code>double[]</code> state of length {@link #dimension()}
     *            (not modified)
     * @param dydt
     *            a <code>double[]</code> of length {@link #dimension()}
     *            receiving the rate of change (modified)
     */
    void derivative(double t, double[] y, double[] dydt);

    /**
     * Advances the state by one step of length {@code h}, which may be
     * negative to integrate backwards.
     * <p>
     * Neither the arguments nor the result are checked for being finite. A
     * field that returns {@code NaN} for a state the step reached produces a
     * {@code NaN} here, and it is the driver that decides whether that is a
     * step to shrink or a problem to give up on.
     *
     * @param t
     *            the time the step starts at
     * @param y
     *            a <code>double[]</code> state of length {@link #dimension()}
     *            at time {@code t} (not modified)
     * @param h
     *            the length of the step, not zero
     * @param yOut
     *            a <code>double[]</code> of length {@link #dimension()}
     *            receiving the state at <code>t + h</code> (modified), which
     *            must not be the array passed as {@code y}
     * @param errOut
     *            a <code>double[]</code> of length {@link #dimension()}
     *            receiving a componentwise estimate of the error of this step
     *            (modified), or {@code null} to ask for none; ignored by a
     *            method that has no estimate to give
     */
    void step(double t, double[] y, double h, double[] yOut, double[] errOut);

    /**
     * The state inside the step that was taken last, at the fraction
     * {@code theta} of it.
     * <p>
     * <b>The endpoints are exact.</b> {@code theta} of zero reproduces the
     * state the last step started from and {@code theta} of one the state it
     * reached, both bit for bit, so an output point that happens to fall on a
     * step boundary cannot disagree with the step itself. In between, a
     * stepper reporting {@link #hasDenseOutput()} interpolates at an order the
     * implementation names; one that does not falls back on the cubic through
     * the two states and the two derivatives.
     *
     * @param theta
     *            the fraction of the last step, zero at its start and one at
     *            its end; a value outside that range extrapolates, at an
     *            accuracy nothing here claims
     * @param out
     *            a <code>double[]</code> of length {@link #dimension()}
     *            receiving the state (modified)
     * @throws IllegalStateException
     *             if no step has been taken since the last {@link #reset()}
     */
    void interpolate(double theta, double[] out);

    /**
     * How stiff the equation looked over the step that was taken last, as
     * <code>|h| |lambda|</code> with {@code lambda} the largest eigenvalue of
     * the Jacobian the step could see.
     * <p>
     * An explicit method is stable only while that product stays inside a
     * region of a few units across, so a value well past three means the step
     * size is being held down by stability and not by accuracy -- which is what
     * stiffness is, and what an explicit method cannot do anything about. It
     * costs no evaluation: it is a difference quotient between two stages the
     * step computed at the same time and at different states, and a method with
     * no such pair of stages simply has no estimate to give.
     *
     * @return the estimate, or {@link Double#NaN} if the method cannot form one
     *         or no step has been taken
     */
    double stiffnessMeasure();

    /**
     * How large {@link #stiffnessMeasure()} may grow over a step before the
     * step counts as one the equation is forcing rather than the tolerance.
     * <p>
     * It is not a universal number: it is how far the stability region of the
     * particular method reaches along the negative real axis, and a method of a
     * higher order has a longer one. The default is Hairer's value for a fifth
     * order pair, which is what a stepper whose measure is always
     * {@link Double#NaN} may as well answer, and a method that reaches further
     * says so.
     *
     * @return the threshold, positive and finite
     */
    default double stiffnessThreshold() {
        return 3.25;
    }

    /**
     * A controller this stepper has to be driven by, rather than any controller
     * with the same settings, or {@code null} for one that does not care --
     * which almost every method is.
     * <p>
     * <b>Why a stepper would care.</b> The division of labor described above
     * holds because a step is one thing and judging it is another. A stepper
     * that chooses <em>between</em> methods breaks that in one narrow place: to
     * decide whether the cheaper method could carry the next step it has to
     * know what an acceptable step is, and that is the tolerances, which live
     * in the controller. Answering here is how such a stepper says that the
     * controller it was built with and the one the driver holds must be the
     * same object, so that the decision and the acceptance cannot be made
     * against different tolerances.
     * <p>
     * {@link OdeIntegrator} checks this when it is constructed and refuses the
     * pair rather than letting the two disagree silently.
     *
     * @return the controller this stepper is bound to, or {@code null}
     */
    default StepController requiredController() {
        return null;
    }

    /**
     * How often the field has been evaluated since this stepper was created.
     * The count runs across {@link #reset()}, so a driver can report the cost
     * of a whole solve.
     *
     * @return the number of evaluations, never decreasing
     */
    long evaluations();

    /**
     * Forgets the last step and anything cached from it, leaving the evaluation
     * count alone. A driver calls this when the state it is advancing changes
     * behind the stepper's back, which is what happens when a solve begins and
     * when an event stops one short.
     */
    void reset();
}
