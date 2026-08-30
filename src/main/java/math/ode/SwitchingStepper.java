package math.ode;

import math.fun.DVectorField;
import math.fun.DiffDVectorField;

/**
 * Two methods behind one stepper, an explicit one and an implicit one, with the
 * step itself deciding which of them the equation currently wants.
 * <p>
 * <b>Why a switch is worth anything.</b> An explicit method is cheap and is
 * stable only while the step size times the largest eigenvalue stays inside a
 * region a few units across; an implicit one has no such limit and costs a
 * Jacobian, a factorization and a back substitution per step. What makes
 * switching pay is an equation that is stiff only for <em>part</em> of its run
 * -- a relaxation oscillator, a reaction with a fast transient, anything with a
 * dial in it.
 * <p>
 * <b>The decision is a price, not a diagnosis.</b> Every {@code probeEvery}
 * steps the method that is <em>not</em> stepping is given the step to try, and
 * what that trial costs per unit of time is compared with what the method in
 * charge costs per unit of time. Neither side can diagnose stiffness on its own
 * -- an implicit method damps the fast modes between its own stages and so
 * never sees what limits an explicit one, and an explicit method whose step has
 * settled at the stability limit sits just <em>below</em>
 * {@link OdeStepper#stiffnessThreshold()} rather than above it -- and a cost
 * comparison needs neither reading. The two directions ask the same question
 * but not by the same reading: {@link #theImplicitSettlesCheaper} says why one
 * trial step is enough coming back and not enough going out, and
 * {@link #stabilityVetoesTheReturn} is the one thing the price alone cannot say.
 * <p>
 * <b>A successful trial is not thrown away.</b> It is the step, so a trial that
 * says "switch" costs nothing at all and only a trial that says "stay" costs
 * the evaluations it took.
 * <p>
 * <b>What it is worth.</b> Measured against the better of the two methods run
 * alone, and against an oracle switching at exactly the right instants, it
 * closes <b>81 %</b> and <b>67 %</b> of the distance to that oracle on the dial
 * and on van der Pol at {@code mu = 1000}, and costs about 3 % on Robertson,
 * which is stiff throughout and has nothing to win. Where there is no stiffness
 * at all it takes no trial and is the pure explicit run <em>bit for bit</em>,
 * so the premium is paid only by an equation that really is stiff.
 * <p>
 * <b>The one regime where it loses is a large system without a written
 * Jacobian</b>, and by enough to matter: both factors of a trial's bill grow
 * with the dimension at once, since the settled step has to climb past
 * <code>(n + 7) / 6</code> and every climbing step costs
 * <code>stages + n + 1</code> evaluations. On a semi-discretized
 * diffusion-reaction PDE at {@code n = 100} the pair costs <b>2.20 times</b> the
 * better pure method with the Jacobian differenced and <b>1.11</b> with it
 * written. <b>So: with a large system and no analytic Jacobian, prefer
 * {@link Ode#solveStiff} outright.</b> Where a Jacobian can be written, the
 * {@link DiffDVectorField} constructor removes the problem rather than easing
 * it -- it drops the ratio a switch must clear to a constant, so it stops
 * depending on {@code n} at all.
 * <p>
 * <b>Where switching is not the answer</b> is an equation with no stiffness at
 * all, which {@link Ode#solve(math.fun.DVectorField, double, double[], double,
 * double)} does at the same cost and with none of this apparatus, and one that
 * is stiff from beginning to end, which
 * {@link Ode#solveStiff(math.fun.DVectorField, double, double[], double,
 * double)} does a few percent cheaper. This is for the equations that are
 * neither, and for the caller who does not know which of the three they have.
 * <p>
 * <b>What it reports about stiffness.</b> {@link #stiffnessMeasure()} is
 * {@link Double#NaN} while the implicit side is stepping, because there is
 * genuinely nothing there to measure, so
 * {@link OdeIntegrator.Result#seemsStiff} speaks for the explicit side alone
 * and is not the question to ask here; {@link #stiffSteps()} and
 * {@link #switches()} are.
 * <p>
 * Instances are stateful and cannot be shared between threads.
 * <p>
 * <b>See</b> <a href="https://en.wikipedia.org/wiki/Stiff_equation">Wikipedia
 * stiff equation</a>.
 *
 * @since 1.5.3
 */
public final class SwitchingStepper implements OdeStepper {

    /**
     * How many steps go by before the method that is not stepping is given one
     * to try. Swept over {@code 5}, {@code 10}, {@code 20} and {@code 40} on
     * four problems and six tolerances, where {@code 20} was best or equal best
     * throughout.
     * <p>
     * <b>A caller who already knows the equation is stiff from beginning to end
     * can pass {@code 50}</b> to the constructor that takes it and save about
     * two percent, since every trial such a run takes is one it throws away.
     * That is the whole of what the cadence is worth: adaptive policies -- a
     * budget on the trial overhead, or backing off while the answers keep coming
     * back the same -- were measured twice and the best of them beat this
     * constant by {@code 0.9} percent while adding a setting of its own.
     */
    public static final int DEFAULT_PROBE_EVERY = 20;

    /**
     * The share of {@link OdeStepper#stiffnessThreshold()} the free estimate has
     * to reach before an implicit trial is worth taking at all. Half leaves a
     * factor of two below the value a stability-limited step settles at, which
     * is about {@code 2.97} of {@code 3.25} for Dormand-Prince, and is still
     * some thirty times what the measure reads on a smooth problem.
     */
    private static final double VETO_FRACTION = 0.5;

    /**
     * How many times {@link OdeStepper#stiffnessThreshold()} the bound
     * <code>|h| ||J||_inf</code> has to reach before the equation is not given
     * back to the explicit method, whatever the trial step cost. It sits in the
     * middle of eight clear decades between a return that was kept and one that
     * killed the run; see {@link #stabilityVetoesTheReturn} for the populations
     * and for why the bound and not the free estimate.
     */
    private static final double RETURN_VETO_FACTOR = 1.0e6;

    /**
     * The most steps the implicit method may take under its own step control
     * before it is judged. It is a ceiling and not a target: the trial stops as
     * soon as the answer is settled either way, and reaches this only where the
     * step size is still climbing and has not yet climbed far enough. See
     * {@link #theImplicitSettlesCheaper}.
     */
    private static final int MAX_SETTLE_STEPS = 16;

    /**
     * The growth per accepted trial step under which the implicit method's step
     * size counts as having stopped climbing. See
     * {@link #theImplicitSettlesCheaper}.
     */
    private static final double STALL_FACTOR = 1.05;

    private final ExplicitRungeKutta nonStiff;
    private final Rosenbrock stiff;
    private final Rosenbrock scout;
    private final StepController controller;
    private final int n;
    private final int probeEvery;
    private final double threshold;
    private final double explicitCost;
    private final double implicitCost;

    private final double[] scoutY;
    private final double[] scoutNext;
    private final double[] scoutError;

    private boolean onStiff;
    private boolean lastStiff;
    private boolean stepped;
    private int sinceProbe;
    private boolean alarmed;
    private boolean aboveThreshold;

    private long switches;
    private long trials;
    private long stiffAttempts;
    private long nonStiffAttempts;

    /**
     * Dormand-Prince against RODAS4, with the Jacobian differenced out of the
     * field.
     *
     * @param field
     *            the right hand side of the equation
     * @param dimension
     *            the number of components of the state
     * @param controller
     *            the controller the {@link OdeIntegrator} will be given, which
     *            has to be that very object
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the dimension is not
     *             positive
     */
    public SwitchingStepper(DVectorField field, int dimension, StepController controller) {
        this(ButcherTableau.DORMAND_PRINCE_45, RosenbrockTableau.RODAS4, field, null, dimension,
                controller, DEFAULT_PROBE_EVERY);
    }

    /**
     * Dormand-Prince against RODAS4, with the Jacobian the field supplies.
     *
     * @param field
     *            the right hand side together with its two first derivatives
     * @param dimension
     *            the number of components of the state
     * @param controller
     *            the controller the {@link OdeIntegrator} will be given, which
     *            has to be that very object
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the dimension is not
     *             positive
     */
    public SwitchingStepper(DiffDVectorField field, int dimension, StepController controller) {
        this(ButcherTableau.DORMAND_PRINCE_45, RosenbrockTableau.RODAS4, field, field, dimension,
                controller, DEFAULT_PROBE_EVERY);
    }

    /**
     * Any pair of methods, at a cadence of its caller's choosing.
     *
     * @param nonStiff
     *            the explicit method, which must estimate the error of its own
     *            steps
     * @param stiff
     *            the implicit method, which must too
     * @param field
     *            the right hand side of the equation
     * @param dimension
     *            the number of components of the state
     * @param controller
     *            the controller the {@link OdeIntegrator} will be given, which
     *            has to be that very object
     * @param probeEvery
     *            how many steps go by before the method that is not stepping is
     *            given one to try, at least one
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, the dimension or the cadence
     *             is not positive, or a method cannot do what is asked of it
     */
    public SwitchingStepper(ButcherTableau nonStiff, RosenbrockTableau stiff, DVectorField field,
            int dimension, StepController controller, int probeEvery) {
        this(nonStiff, stiff, field, null, dimension, controller, probeEvery);
    }

    /**
     * Any pair of methods, with the Jacobian the field supplies.
     *
     * @param nonStiff
     *            the explicit method
     * @param stiff
     *            the implicit method
     * @param field
     *            the right hand side together with its two first derivatives
     * @param dimension
     *            the number of components of the state
     * @param controller
     *            the controller the {@link OdeIntegrator} will be given, which
     *            has to be that very object
     * @param probeEvery
     *            how many steps go by before the method that is not stepping is
     *            given one to try, at least one
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, the dimension or the cadence
     *             is not positive, or a method cannot do what is asked of it
     */
    public SwitchingStepper(ButcherTableau nonStiff, RosenbrockTableau stiff, DiffDVectorField field,
            int dimension, StepController controller, int probeEvery) {
        this(nonStiff, stiff, field, field, dimension, controller, probeEvery);
    }

    private SwitchingStepper(ButcherTableau explicitTableau, RosenbrockTableau implicitTableau,
            DVectorField field, DiffDVectorField analytic, int dimension, StepController controller,
            int probeEvery) {
        if (explicitTableau == null) {
            throw new IllegalArgumentException("the explicit tableau must not be null");
        }
        if (implicitTableau == null) {
            throw new IllegalArgumentException("the implicit tableau must not be null");
        }
        if (controller == null) {
            throw new IllegalArgumentException("controller must not be null");
        }
        if (probeEvery < 1) {
            throw new IllegalArgumentException("probeEvery must be at least 1, not " + probeEvery);
        }
        if (!explicitTableau.hasEmbedded()) {
            throw new IllegalArgumentException(explicitTableau
                    + " estimates no error of its own, so its steps cannot be judged");
        }
        if (!implicitTableau.hasErrorEstimate()) {
            throw new IllegalArgumentException(implicitTableau
                    + " estimates no error of its own, so its steps cannot be judged");
        }
        this.nonStiff = new ExplicitRungeKutta(explicitTableau, field, dimension);
        this.stiff = (analytic == null) ? new Rosenbrock(implicitTableau, field, dimension)
                : new Rosenbrock(implicitTableau, analytic, dimension);
        // a second copy of the implicit method, so that letting it run ahead to
        // see where its step size settles leaves the one carrying the solution
        // untouched and still able to interpolate over the step it last took
        this.scout = (analytic == null) ? new Rosenbrock(implicitTableau, field, dimension)
                : new Rosenbrock(implicitTableau, analytic, dimension);
        this.controller = controller;
        this.n = dimension;
        this.probeEvery = probeEvery;
        this.threshold = explicitTableau.stiffnessThreshold();
        // what one step of either method evaluates the field: the explicit one
        // once per stage, less the one a first-same-as-last method carries over;
        // the implicit one once per stage but for the first, which reuses the
        // evaluation the linearization took, and once per component and once
        // more for the time where the Jacobian has to be differenced out. The
        // factorization the implicit step also pays is not an evaluation and
        // does not appear here, so this understates it and the trial errs
        // towards staying where it is
        this.explicitCost = explicitTableau.stages() - (explicitTableau.isFsal() ? 1 : 0);
        this.implicitCost = implicitTableau.stages()
                + (this.stiff.hasAnalyticJacobian() ? 0 : dimension + 1);
        this.scoutY = new double[dimension];
        this.scoutNext = new double[dimension];
        this.scoutError = new double[dimension];
    }

    /**
     * The explicit method, for a caller who wants its counts on their own.
     *
     * @return the stepper that carries the non-stiff stretches
     */
    public ExplicitRungeKutta nonStiffStepper() {
        return nonStiff;
    }

    /**
     * The implicit method, for a caller who wants its counts on their own.
     *
     * @return the stepper that carries the stiff stretches
     */
    public Rosenbrock stiffStepper() {
        return stiff;
    }

    /**
     * The controller this stepper judges its own trials by, which is the one
     * the {@link OdeIntegrator} has to be given as well.
     *
     * @return the controller from the constructor
     */
    public StepController controller() {
        return controller;
    }

    /**
     * Whether the implicit method is the one that would take the next step.
     *
     * @return {@code true} while the equation is being treated as stiff
     */
    public boolean isStiffActive() {
        return onStiff;
    }

    /**
     * How often the two methods have changed places since the last
     * {@link #reset()}, counted in both directions together.
     *
     * @return the number of switches, zero for a run that never left the
     *         explicit method
     */
    public long switches() {
        return switches;
    }

    /**
     * How often the explicit method was given a step to try since the last
     * {@link #reset()}, which is the whole cost of looking for the way back.
     * <p>
     * A trial that succeeded became the step and cost nothing beyond it; the
     * difference between this and {@link #switches()} is how many were thrown
     * away.
     *
     * @return the number of trial steps taken
     */
    public long trials() {
        return trials;
    }

    /**
     * How many steps the implicit method was asked for since the last
     * {@link #reset()}.
     * <p>
     * These are steps <em>attempted</em>, since a stepper is never told which
     * of its steps the driver kept; {@link OdeIntegrator.Result#steps} and
     * {@link OdeIntegrator.Result#rejected} give that split for the run as a
     * whole. The steps taken to find out where the implicit method's step size
     * settles are not among them, since they advance nothing;
     * {@link #explorationEvaluations()} is what those cost.
     *
     * @return the number of implicit steps attempted
     */
    public long stiffSteps() {
        return stiffAttempts;
    }

    /**
     * How many steps the explicit method was asked for since the last
     * {@link #reset()}, attempted rather than kept, as {@link #stiffSteps()}
     * describes.
     *
     * @return the number of explicit steps attempted
     */
    public long nonStiffSteps() {
        return nonStiffAttempts;
    }

    @Override
    public int dimension() {
        return n;
    }

    /**
     * The order of the method that took the last step, which is what the
     * controller has to scale the next step size with. Before any step it is
     * the order of the method that would take one.
     */
    @Override
    public int order() {
        return active().order();
    }

    @Override
    public boolean hasErrorEstimate() {
        return nonStiff.hasErrorEstimate() && stiff.hasErrorEstimate();
    }

    @Override
    public boolean hasDenseOutput() {
        return active().hasDenseOutput();
    }

    /**
     * How far the stability region of the <em>explicit</em> method reaches,
     * which is the only side a measure is ever formed on.
     */
    @Override
    public double stiffnessThreshold() {
        return threshold;
    }

    /**
     * The controller given to the constructor, which the {@link OdeIntegrator}
     * driving this stepper has to be given as well: the trial step that decides
     * whether the explicit method can take over is judged by it, and judging it
     * by tolerances the driver does not accept by would let the two disagree.
     */
    @Override
    public StepController requiredController() {
        return controller;
    }

    @Override
    public void derivative(double t, double[] y, double[] dydt) {
        active().derivative(t, y, dydt);
    }

    @Override
    public void step(double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (y == null || y.length != n) {
            throw new IllegalArgumentException("y must be of length " + n);
        }
        // steps are counted as they are attempted and not as they are kept: a
        // run whose steps are being rejected is a run that may want the other
        // method, and one that never probed until a step was accepted could
        // shrink towards zero without ever asking
        if (errOut != null && dueForTrial()) {
            ++trials;
            sinceProbe = 0;
            alarmed = false;
            boolean cheaper;
            if (onStiff) {
                // the trial is taken into the arrays the answer goes into, so a
                // trial that says "switch" is the step and costs nothing at all;
                // only one that says "stay" is thrown away
                take(false, t, y, h, yOut, errOut);
                cheaper = theOtherLooksCheaper(y, yOut, errOut, h);
            } else {
                cheaper = theImplicitSettlesCheaper(t, y, h);
                if (cheaper) {
                    take(true, t, y, h, yOut, errOut);
                }
            }
            if (cheaper) {
                onStiff = !onStiff;
                ++switches;
                lastStiff = onStiff;
                stepped = true;
                noteAlarm();
                return;
            }
        }

        ++sinceProbe;
        take(onStiff, t, y, h, yOut, errOut);
        lastStiff = onStiff;
        stepped = true;
        noteAlarm();
    }

    @Override
    public void interpolate(double theta, double[] out) {
        if (!stepped) {
            throw new IllegalStateException("no step has been taken");
        }
        active().interpolate(theta, out);
    }

    @Override
    public double stiffnessMeasure() {
        if (!stepped || lastStiff) {
            return Double.NaN;
        }
        return nonStiff.stiffnessMeasure();
    }

    @Override
    public double[] secondaryError() {
        return active().secondaryError();
    }

    @Override
    public long evaluations() {
        return nonStiff.evaluations() + stiff.evaluations() + scout.evaluations();
    }

    /**
     * The part of {@link #evaluations()} spent letting the implicit method find
     * its own step size rather than advancing the solution.
     * <p>
     * It is the price of the trial being able to see past one step, and it is
     * paid only while the explicit method is stepping and the free estimate has
     * not vetoed the trial: nothing at all on an equation with no stiffness in
     * it, and between 3 % and 23 % of the run on one where the two methods
     * change places, growing with the dimension.
     *
     * @return the number of evaluations spent exploring, never decreasing
     */
    public long explorationEvaluations() {
        return scout.evaluations();
    }

    /**
     * Forgets the last step on both sides and returns to the explicit method,
     * along with the switch and step counts, which are therefore those of one
     * run. The evaluation counts are left alone, as the interface requires.
     */
    @Override
    public void reset() {
        nonStiff.reset();
        stiff.reset();
        scout.reset();
        onStiff = false;
        lastStiff = false;
        stepped = false;
        sinceProbe = 0;
        alarmed = false;
        aboveThreshold = false;
        switches = 0L;
        trials = 0L;
        stiffAttempts = 0L;
        nonStiffAttempts = 0L;
    }

    /**
     * The two methods, which of them is stepping, and how often they have
     * changed places.
     */
    @Override
    public String toString() {
        return "SwitchingStepper[" + nonStiff.tableau() + " against " + stiff.tableau() + ", on the "
                + (onStiff ? "stiff" : "non-stiff") + " side, " + switches + " switches, dimension " + n
                + "]";
    }

    /** The side that took the last step, or the one that would take the next. */
    private OdeStepper active() {
        if (stepped) {
            return lastStiff ? stiff : nonStiff;
        }
        return onStiff ? stiff : nonStiff;
    }

    /**
     * Whether this step should be given to the other method to try, which is
     * where the free estimate is used -- to time the trials, never to decide
     * one.
     * <p>
     * <b>The free estimate is no good at saying "stiff", very good at saying
     * "not stiff", and it does not lie when it does speak up.</b> That makes it
     * a schedule rather than a test, and while the explicit method is stepping
     * it reads three ways:
     * <ul>
     * <li><b>Rising past the threshold, ask now.</b> A step outside the
     * stability region that the local error test kept anyway, which happens at
     * a loose tolerance. This is Hairer's signal and its known failure is
     * silence rather than a false alarm, so acting at once is safe. <b>On the
     * rising edge and not on the level</b>, or a run that stays above the
     * threshold while the trial keeps saying "stay" asks on every step.</li>
     * <li><b>Far below the threshold, never ask.</b> An implicit method of a
     * <em>lower</em> order cannot take longer steps than an explicit one whose
     * step is short for accuracy rather than stability, so the trial could only
     * come back saying "stay". This is what keeps an equation with no stiffness
     * in it from paying anything at all.</li>
     * <li><b>In between, ask on the cadence.</b> A step size that has
     * <em>settled</em> at the stability limit sits just under the threshold,
     * and nothing but a trial sees into that band.</li>
     * </ul>
     * A measure that is not a number is no information, so the cadence decides.
     * On the implicit side there is no such estimate and the cadence is all
     * there is.
     */
    private boolean dueForTrial() {
        if (alarmed) {
            return true;
        }
        if (onStiff || !stepped) {
            return sinceProbe >= probeEvery;
        }
        double measure = nonStiff.stiffnessMeasure();
        if (!Double.isNaN(measure) && measure <= VETO_FRACTION * threshold) {
            return false;
        }
        return sinceProbe >= probeEvery;
    }

    /**
     * Notes whether the step just taken crossed <em>up</em> through the
     * threshold, which owes the next step a trial. See {@link #dueForTrial()}
     * for why the edge and not the level.
     */
    private void noteAlarm() {
        if (lastStiff) {
            aboveThreshold = false;
            return;
        }
        double measure = nonStiff.stiffnessMeasure();
        boolean above = !Double.isNaN(measure) && measure > threshold;
        if (above && !aboveThreshold) {
            alarmed = true;
        }
        aboveThreshold = above;
    }

    /**
     * Whether the implicit method, allowed to find its own step size, would now
     * be the cheaper of the two.
     * <p>
     * <b>Why this side cannot be asked the one-step question.</b>
     * {@link #theOtherLooksCheaper} reads the step size a single trial would
     * propose next, <code>h safety err^(-1/p)</code>, which is exact inside the
     * asymptotic regime and useless across it. An implicit method probed at the
     * explicit method's stability-limited step is deep inside that regime while
     * its own natural step is not, and under-predicts by three orders: probed at
     * {@code h = 3.3e-05} it proposed {@code 2.7 h} where a pure implicit run
     * was taking steps <b>1500 times longer</b>. One step cannot show a ramp,
     * and the ramp is the whole of what an implicit method has to offer. The
     * switch then happened only while {@code n < 9}, a bound with nothing in it
     * about the equation. So this side lets the implicit method run under its
     * own control and reads where the step size <em>settles</em>.
     * <p>
     * <b>How long it settles is not a fixed number, and a fixed number is what
     * failed next</b>, because the ratio a switch has to clear grows with
     * {@code n} while what a fixed count reveals does not -- four steps moved
     * the bound only to about {@code n < 21}. Instead the trial is given the
     * ratio it must clear and returns the moment it clears it, and gives up the
     * moment the climb stops -- under {@link #STALL_FACTOR} per accepted step --
     * because a climb that died short will not resume. {@link #MAX_SETTLE_STEPS}
     * is only the ceiling. The same exit that finds a long ramp early also
     * abandons a hopeless probe early, which is why it improves both ends at
     * once: exploration on runs that probe and never switch fell from 30 % to
     * 4.1 %.
     * <p>
     * <b>Only this side.</b> Coming back the other way the one-step reading was
     * never shown wrong, and settling there costs the most, since a trial on the
     * implicit side is nearly always answered "stay": doing it on both sides
     * took Robertson's premium from 2 % to 34 %.
     * <p>
     * <b>What is still not seen.</b> Where the stiff modes <em>rotate</em> as
     * well as decay, an implicit method has to resolve the rotation and cannot
     * ramp away from the explicit method's step at all, so no amount of settling
     * finds a ramp that is not there. There the trial is right to answer "stay"
     * and the cost of asking is what remains.
     */
    private boolean theImplicitSettlesCheaper(double t, double[] y, double h) {
        double settled = settledStep(t, y, h);
        if (!(settled > 0.0)) {
            return false;
        }
        if (Double.isInfinite(settled)) {
            return true;
        }
        return implicitCost / settled < explicitCost / Math.abs(h);
    }

    /**
     * The step size the implicit method reaches under its own step control,
     * starting from {@code h}, or a negative number if it cannot take even one
     * step. It stops early once the reading can no longer change the verdict --
     * once it is past what a switch needs, or once it has stopped climbing --
     * and at {@link #MAX_SETTLE_STEPS} steps in any case.
     * <p>
     * This is the loop {@link OdeIntegrator} runs, on a copy of the method so
     * that the one carrying the solution is not disturbed, and bounded in
     * attempts so that a point where the implicit method cannot make progress
     * either costs a fixed amount rather than an unbounded one.
     * <p>
     * <b>{@code enough} is not a truncation of the decision, it is the decision,
     * written as an early exit.</b> <code>|trial| &gt; enough</code> says exactly
     * that <code>implicitCost / trial &lt; explicitCost / |h|</code> holds --
     * which is why the loop is allowed to stop there, and why the two must be
     * changed together or not at all. Moving one alone opens a band in which the
     * probe stops before it can answer and the answer is then always no: with
     * {@code enough} lowered by a sixth on its own, the probe exits after a
     * single step and the unscaled comparison refuses the result, which on
     * OREGO turned 33 trials into 10006 and the run into <b>342 times</b> its
     * cost. Most probes end at the stall rather than here -- 19 of 22 and 42 of
     * 47 on the reference problems -- and the ceiling is reached at most twice
     * a run.
     */
    private double settledStep(double t, double[] y, double h) {
        scout.reset();
        System.arraycopy(y, 0, scoutY, 0, n);
        double[] current = scoutY;
        double[] next = scoutNext;
        double time = t;
        double trial = h;
        double previous = StepController.ERROR_FLOOR;
        int order = scout.order();
        int taken = 0;
        int attempts = 0;
        double enough = Math.abs(h) * implicitCost / explicitCost;
        while (taken < MAX_SETTLE_STEPS && attempts < 8 * MAX_SETTLE_STEPS + 8) {
            ++attempts;
            scout.step(time, current, trial, next, scoutError);
            boolean finite = allFinite(next) && allFinite(scoutError);
            double scaled = finite ? controller.errorNorm(scoutError, current, next) : Double.NaN;
            if (finite && scaled <= 1.0) {
                ++taken;
                double[] swap = current;
                current = next;
                next = swap;
                time += trial;
                double factor = controller.scale(scaled, previous, order);
                previous = Math.max(scaled, StepController.ERROR_FLOOR);
                trial *= factor;
                if (Math.abs(trial) > enough) {
                    // the answer is yes and no further step can unmake it
                    return trial;
                }
                if (factor < STALL_FACTOR) {
                    // the climb has died away short of what was needed, and
                    // waiting longer only spends more to hear the same no
                    return trial;
                }
            } else {
                trial *= controller.scaleAfterRejection(scaled, order);
            }
        }
        return (taken == 0) ? -1.0 : trial;
    }

    /** One step of whichever side is named, counted against that side. */
    private void take(boolean useStiff, double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (useStiff) {
            stiff.step(t, y, h, yOut, errOut);
            ++stiffAttempts;
        } else {
            nonStiff.step(t, y, h, yOut, errOut);
            ++nonStiffAttempts;
        }
    }

    /**
     * Whether the trial step just taken says the method that took it would now
     * be the cheaper of the two.
     * <p>
     * <b>The question is not whether the trial was acceptable.</b> An
     * acceptable step says nothing about how long the next one may be, and it
     * is the length of the next one that decides what a method costs. What the
     * trial gives is its own error, and from that the step size its method
     * would propose; the two cost densities -- evaluations per unit of time --
     * then compare directly. This is the same question the side that is
     * currently stepping has already answered for itself, since the step size
     * it is being driven at is the one it proposed.
     * <p>
     * <b>Never mentioning stability is what makes it work in both directions,
     * and it is also its own blind spot.</b> A method that is cheap per unit of
     * time and unstable is not a method that can carry the equation, and only
     * the way <em>back</em> can choose one -- which is what
     * {@link #stabilityVetoesTheReturn} is for.
     * <p>
     * The proposal is formed unclamped, without the bounds a controller puts on
     * a single step: those keep one step from jumping, and this is a comparison
     * rather than a step.
     */
    private boolean theOtherLooksCheaper(double[] y, double[] yOut, double[] errOut, double h) {
        if (!theOtherCostsLess(y, yOut, errOut, h)) {
            return false;
        }
        return !stabilityVetoesTheReturn(h);
    }

    /**
     * The cost comparison itself, which is the whole of the decision on the way
     * to the implicit side and all but the last word on the way back. See
     * {@link #theOtherLooksCheaper}.
     */
    private boolean theOtherCostsLess(double[] y, double[] yOut, double[] errOut, double h) {
        if (!allFinite(yOut) || !allFinite(errOut)) {
            return false;
        }
        OdeStepper trial = onStiff ? nonStiff : stiff;
        double[] second = trial.secondaryError();
        if (second != null && !allFinite(second)) {
            return false;
        }
        double scaled = (second == null) ? controller.errorNorm(errOut, y, yOut)
                : controller.errorNorm(errOut, second, y, yOut);
        if (Double.isNaN(scaled) || Double.isInfinite(scaled)) {
            return false;
        }
        double here = (onStiff ? implicitCost : explicitCost) / Math.abs(h);
        double thereCost = onStiff ? explicitCost : implicitCost;
        if (scaled <= 0.0) {
            // an exact step, which nothing about the accuracy is holding back
            return true;
        }
        double proposed = Math.abs(h) * controller.safety() * Math.pow(scaled, -1.0 / trial.order());
        if (!(proposed > 0.0)) {
            return false;
        }
        if (Double.isInfinite(proposed)) {
            return true;
        }
        return thereCost / proposed < here;
    }

    /**
     * Whether the explicit method would be so far outside its stability region
     * that it must not be handed the equation, whatever the trial step cost.
     * <p>
     * <b>This is the counterpart {@link #dueForTrial()} has and this direction
     * did not.</b> Going to the implicit side the free estimate vetoes a trial
     * that has nothing to find; coming back, nothing stopped a trial that
     * happened to be cheap from taking the equation off a method that was
     * holding it together, and that kills a run rather than slowing one.
     * <p>
     * <b>The reading is the implicit method's own Jacobian and not the free
     * estimate</b>, which is the part that had to be measured rather than
     * assumed: where {@link Rosenbrock#jacobianNorm()} is right,
     * {@link OdeStepper#stiffnessMeasure()} can be fourteen decades low, because
     * the difference quotient loses a fast mode whose own combination has
     * decayed away while the norm is carried by components far apart. Over 172
     * returns in a corpus of 64 runs, <code>|h| ||J||_inf</code> against the
     * threshold separated a kept return (up to {@code 2.1e+03}) from a run that
     * dies ({@code 4.6e+11} and above) by <b>eight clear decades</b>, and
     * {@link #RETURN_VETO_FACTOR} sits in the middle of them.
     * <p>
     * <b>What it is honestly.</b> At those decisions the explicit trial is not
     * malfunctioning -- it converges at its own order and returns a smooth
     * state, and every local signal endorses the handover. This vetoes it
     * anyway, because a stability margin of {@code 1e12} leaves nothing to
     * recover from when the fast mode comes back. It is a veto on the margin
     * and not on the reading.
     * <p>
     * <b>What it costs.</b> One pass over a matrix the implicit method already
     * holds and factors on every step anyway, taken only where the cost
     * comparison has already said "switch" -- 2.3 % of return trials. The
     * reading is as old as the last linearization, which at this point is the
     * start of the previous accepted step; measured, that staleness moves none
     * of the populations. Where there is no Jacobian to read the answer is
     * {@link Double#NaN} and nothing is vetoed, the same way
     * {@link #dueForTrial()} treats a measure it cannot form: no information is
     * not evidence.
     */
    private boolean stabilityVetoesTheReturn(double h) {
        if (!onStiff) {
            return false;
        }
        double norm = stiff.jacobianNorm();
        if (Double.isNaN(norm)) {
            return false;
        }
        return Math.abs(h) * norm > RETURN_VETO_FACTOR * threshold;
    }

    private static boolean allFinite(double[] x) {
        for (int i = 0; i < x.length; ++i) {
            if (Double.isNaN(x[i]) || Double.isInfinite(x[i])) {
                return false;
            }
        }
        return true;
    }
}
