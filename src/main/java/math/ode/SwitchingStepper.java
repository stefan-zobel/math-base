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
 * Jacobian, a factorization and a back substitution per step. On an equation
 * that is stiff throughout, the implicit method wins everywhere and there is
 * nothing to switch. What makes switching pay is an equation that is stiff
 * only for part of its run -- a relaxation oscillator, a reaction with a fast
 * transient, anything with a dial in it.
 * <p>
 * <b>The decision is a price, not a diagnosis.</b> Every {@code probeEvery}
 * steps the method that is <em>not</em> stepping is given the step to try, and
 * what that trial costs per unit of time is compared with what the method in
 * charge costs per unit of time. The same question decides both directions,
 * though not by the same reading: coming back to the explicit method one trial
 * step is enough, and going to the implicit one it is not, because a
 * stability-limited step size is no baseline to extrapolate from.
 * {@link #theImplicitSettlesCheaper} is where that is measured out.
 * <p>
 * <b>Neither side can diagnose stiffness on its own, which is why the price is
 * asked instead.</b> Both of the obvious tests were built first and measured,
 * and both have a blind spot:
 * <ul>
 * <li>an <b>implicit</b> method damps the fast modes <em>between its own
 * stages</em>, so no difference quotient it can form sees what would limit an
 * explicit method; on Robertson's reaction such an estimate disagreed with the
 * truth on half the steps;</li>
 * <li>an <b>explicit</b> method whose step size has settled at the stability
 * limit sits just <em>below</em> {@link OdeStepper#stiffnessThreshold()} rather
 * than above it, so a threshold test almost never fires. On a dial rising
 * smoothly from <code>lambda = 1</code> to <code>1e5</code>,
 * {@link OdeStepper#stiffnessMeasure()} read a mean of {@code 2.97} and a worst
 * of {@code 3.2501} against a threshold of {@code 3.25} over 29158 steps,
 * exceeding it on five. Hairer's test catches a run being pushed past the limit
 * and rejected; it does not catch one that has already adapted to it.</li>
 * </ul>
 * A cost comparison has neither blind spot, because it never mentions
 * stability: an explicit method held down by stability is simply an explicit
 * method whose steps are short, and short steps are expensive.
 * <p>
 * <b>The free estimate is still used, but to time the trials rather than to
 * decide one.</b> It is bad at saying "stiff", good at saying "not stiff", and
 * it does not lie when it does speak up, so it schedules: far below its
 * threshold no trial is taken at all, on the rising edge past it one is taken
 * at once, and in between the cadence decides. {@link #dueForTrial()} has the
 * three cases and what each was measured against. The first is what keeps an
 * equation with no stiffness in it from paying for any of this -- on the
 * harmonic oscillator no trial is taken and the run is the pure explicit run
 * <em>bit for bit</em>.
 * <p>
 * <b>A successful trial is not thrown away.</b> It is the step, so a trial that
 * says "switch" costs nothing at all and only a trial that says "stay" costs
 * the evaluations it took.
 * <p>
 * <b>What it costs, measured</b> against the better of the two methods run
 * alone, and against an oracle switching at exactly the right instants, which
 * bounds every switching solver from above:
 * <table border="1">
 * <caption>evaluations, at {@code probeEvery} of 20</caption>
 * <tr><th>problem</th><th>this</th><th>better pure</th><th>oracle</th></tr>
 * <tr><td>a dial from 1 to 1e5 and back, rtol 1e-6</td>
 * <td>1789</td><td>2928</td><td>1468</td></tr>
 * <tr><td>van der Pol at mu = 1000, rtol 1e-6</td>
 * <td>5017</td><td>8902</td><td>2747</td></tr>
 * <tr><td>Robertson to 1e5, rtol 1e-8</td>
 * <td>3250</td><td>3147</td><td>3073</td></tr>
 * <tr><td>a chain of 50 cells under the same dial, rtol 1e-6</td>
 * <td>8757</td><td>19948</td><td>--</td></tr>
 * </table>
 * <p>
 * The first two are what switching is for, and it closes 78 % and 63 % of the
 * distance to the oracle. Robertson is an equation that is stiff throughout,
 * where there is nothing to win and the trials are an insurance premium of a
 * few percent on a claim never made. Over tolerances from {@code 1e-4} to
 * {@code 1e-12} the two winners stay winners -- van der Pol between 41 % and
 * 85 % of the pure implicit run, the dial never above 99 % of it -- and
 * Robertson's premium stays between 2 % and 17 %.
 * <p>
 * <b>The dimension is not among the things that decide.</b> The last row is
 * there because it once was: a trial judged on one step could only hand over
 * while one implicit step cost less than a couple of explicit ones, which for a
 * differenced Jacobian is {@code n < 9}, and above that the stepper never
 * handed over at all. See {@link #theImplicitSettlesCheaper} for what that was
 * and why it is gone. A caller with a large system should still reach for the
 * {@link DiffDVectorField} constructor where they can: the differencing is
 * {@code n + 1} evaluations a step, and on that chain of 50 the written
 * Jacobian costs 1414 evaluations against 8757.
 * <p>
 * <b>Which explicit method to pair.</b> {@link ButcherTableau#DOP853} works
 * here unchanged and is worth reaching for under two conditions at once: the
 * tolerance is below about {@code 1e-6}, where an eighth order method beats a
 * fifth order one at all, <em>and</em> the explicit side carries a large share
 * of the evaluations. Only that share is exposed to the improvement, so the
 * whole effect is
 * <p>
 * <code>cost(DOP853 pair) / cost(DP45 pair) = (1 - share) + share * q</code>
 * <p>
 * with {@code q} what DOP853 costs on the explicit part alone, and the measured
 * ratios follow that expression to a few percent. What it comes to:
 * <table border="1">
 * <caption>DOP853 pair against the Dormand-Prince pair, at the best tolerance
 * for each</caption>
 * <tr><th>problem</th><th>explicit share of the bill</th><th>best ratio</th></tr>
 * <tr><td>Robertson to 1e5</td><td>4 % to 6 %</td><td>1.03 -- never pays</td></tr>
 * <tr><td>a dial with the stiff window a third of the run</td>
 * <td>9 % to 37 %</td><td>0.98 at 1e-8 -- not worth it</td></tr>
 * <tr><td>van der Pol at mu = 1000</td><td>about 50 %</td>
 * <td>0.86 at 1e-9</td></tr>
 * <tr><td>a dial with the stiff window a twentieth of the run</td>
 * <td>62 % to 83 %</td><td>0.73 at 1e-6</td></tr>
 * </table>
 * <p>
 * <b>The share is of the evaluations and not of the time axis</b>, and the two
 * are nothing alike: on van der Pol the explicit method holds under 1 % of the
 * time and half the bill, because what it carries is the brief violent
 * transitions. <b>And the window closes again at a tight enough tolerance</b>,
 * because the implicit half of the bill grows like <code>rtol^(-1/4)</code>
 * against DOP853's <code>rtol^(-1/9)</code>: on the dial the explicit share
 * falls from 80 % at {@code 1e-6} to 11 % at {@code 1e-11}, and with it the
 * advantage. Only van der Pol, whose share stays near half at every tolerance,
 * keeps the gain all the way down.
 * <p>
 * <b>Where switching is not the answer</b> is an equation with no stiffness at
 * all, which {@link Ode#solve(math.fun.DVectorField, double, double[], double,
 * double)} does at exactly the same cost and with none of this apparatus, and
 * one that is stiff from beginning to end, which
 * {@link Ode#solveStiff(math.fun.DVectorField, double, double[], double,
 * double)} does a few percent cheaper. This is for the equations that are
 * neither, and for the caller who does not know which of the three they have.
 * <p>
 * <b>What it reports about stiffness.</b> {@link #stiffnessMeasure()} is
 * {@link Double#NaN} while the implicit side is stepping, because there is
 * genuinely nothing there to measure, and
 * {@link OdeIntegrator.Result#stiffness} therefore stays exactly what it is
 * documented to be -- the largest such measure actually seen. Nothing here is
 * decided by it. {@link OdeIntegrator.Result#seemsStiff} likewise speaks for the
 * explicit side alone and is not the question to ask here; {@link #stiffSteps()}
 * and {@link #switches()} are.
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
     * to try. Measured over {@code 5}, {@code 10}, {@code 20} and {@code 40} on
     * the three problems above, {@code 20} was the best or equal best on all of
     * them.
     * <p>
     * <b>A caller who already knows the equation is stiff from beginning to
     * end can pass {@code 50}</b> to the constructor that takes it and save
     * about two percent, since every trial such a run takes is one it throws
     * away. That is the whole of what the cadence is worth: policies that vary
     * it -- holding the trial overhead to a share of the bill, or backing off
     * while the answers keep coming back the same -- were measured over four
     * problems and six tolerances, and the best of them beat this constant by
     * {@code 0.9} percent on average while adding a setting of its own.
     */
    public static final int DEFAULT_PROBE_EVERY = 20;

    /**
     * The share of {@link OdeStepper#stiffnessThreshold()} the free estimate has
     * to reach before an implicit trial is worth taking at all.
     */
    private static final double VETO_FRACTION = 0.5;

    /**
     * How many steps the implicit method is allowed to take under its own step
     * control before it is judged. See {@link #theImplicitSettlesCheaper}.
     */
    private static final int SETTLE_STEPS = 4;

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
     * <b>The free estimate is no good at saying "stiff" and is very good at
     * saying "not stiff", and it does not lie when it does speak up.</b> Those
     * three facts make it a schedule rather than a test, and while the explicit
     * method is stepping it reads three ways:
     * <ul>
     * <li><b>Rising past the threshold, ask now.</b> The step that was just
     * taken ran outside the stability region and was kept anyway, which the
     * local error test will do at a loose tolerance -- on Robertson at
     * {@code rtol = 1e-4} the very first step is accepted at a measure of
     * {@code 6.64} against a threshold of {@code 3.25}. Waiting the ordinary
     * cadence out costs the run: twenty attempts later the second concentration
     * stands at {@code -0.12} where it belongs near {@code 1e-05}, and no
     * method recovers from that. This is Hairer's signal, and its known failure
     * is silence rather than a false alarm, so acting on it at once is safe.
     * <b>On the rising edge and not on the level</b>, because a run that stays
     * above the threshold while the trial keeps answering "stay" would
     * otherwise ask on every step -- measured at 39 trials in 40 steps before
     * the edge was put in.</li>
     * <li><b>Far below the threshold, never ask.</b> An implicit method of a
     * <em>lower</em> order cannot take longer steps than an explicit one whose
     * step is short for accuracy rather than stability, so the trial could only
     * come back saying "stay". This is what keeps an equation with no stiffness
     * in it from paying anything at all.</li>
     * <li><b>In between, ask on the cadence.</b> This is the blind band the
     * measurements found -- a step size that has <em>settled</em> at the
     * stability limit sits just under the threshold, at a mean of {@code 2.97}
     * of {@code 3.25} over 29158 steps -- and nothing but a trial sees into it.
     * Half the threshold is the lower edge, which leaves a factor of two below
     * that settled value and is still some thirty times what the measure reads
     * on a smooth problem.</li>
     * </ul>
     * A measure that is not a number -- a method that cannot form one, or a
     * step that came back infinite -- is no information, so the cadence decides.
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
     * asymptotic regime and useless across it. The step size an explicit method
     * is being driven at while it is stability-limited is far below anything the
     * implicit method would choose, so the trial lands deep inside that regime
     * and its own natural step does not: measured on a chain of diffusing cells
     * under a dial of {@code 1e5}, RODAS4 probed at the explicit method's
     * {@code h = 3.3e-05} reported an error of {@code 0.0131} and therefore
     * proposed {@code 2.7 h}, where a pure implicit run at that instant was
     * taking steps <b>1500 times longer</b>. One step cannot show a ramp, and the
     * ramp is the whole of what an implicit method has to offer.
     * <p>
     * <b>What that cost.</b> The switch then happened only while
     * <code>implicitCost &lt; explicitCost * 2.7</code>, which for a differenced
     * Jacobian is <code>n + 7 &lt; 16.2</code>, or {@code n < 9}: a bound with
     * nothing in it about the equation. Above it the stepper never handed over
     * at all -- at {@code n = 50}, 2677596 evaluations against the pure implicit
     * method's 19948, a hundredfold loss. So this side lets the implicit method
     * run {@link #SETTLE_STEPS} steps under its own control and reads where the
     * step size <em>settles</em>.
     * <p>
     * <b>Only this side.</b> Coming back the other way the one-step reading was
     * never shown wrong, and settling there is what costs the most, since a
     * trial on the implicit side is nearly always answered "stay": doing it on
     * both sides took Robertson's premium from 2 % to 34 %, and doing it here
     * alone costs between 0.1 % and 6.8 % on everything that already worked.
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
     * The step size the implicit method reaches after {@link #SETTLE_STEPS} of
     * its own steps, starting from {@code h}, or a negative number if it cannot
     * take even one.
     * <p>
     * This is the loop {@link OdeIntegrator} runs, on a copy of the method so
     * that the one carrying the solution is not disturbed, and bounded in
     * attempts so that a point where the implicit method cannot make progress
     * either costs a fixed amount rather than an unbounded one.
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
        while (taken < SETTLE_STEPS && attempts < 8 * SETTLE_STEPS + 8) {
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
     * <b>It is also the whole of the switching logic, in both directions.</b>
     * The obvious alternative for the way <em>to</em> the implicit side is
     * {@link OdeStepper#stiffnessMeasure()} against its threshold, and that was
     * built first and measured: an explicit method whose step size has settled
     * at the stability limit sits just <em>below</em> the threshold rather than
     * above it, so the test almost never fires. On a dial rising smoothly from
     * <code>lambda = 1</code> to <code>1e5</code> the measure read a mean of
     * {@code 2.97} and a worst of {@code 3.2501} against a threshold of
     * {@code 3.25} over 29158 steps, exceeding it on five of them, and the run
     * cost 44 % more than the pure implicit method rather than 22 % less.
     * Hairer's test catches a run being pushed past the limit and rejected; it
     * does not catch one that has already adapted to it. A cost comparison has
     * no such blind spot, because it never mentions stability.
     * <p>
     * The proposal is formed unclamped, without the bounds a controller puts on
     * a single step: those keep one step from jumping, and this is a comparison
     * rather than a step.
     */
    private boolean theOtherLooksCheaper(double[] y, double[] yOut, double[] errOut, double h) {
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

    private static boolean allFinite(double[] x) {
        for (int i = 0; i < x.length; ++i) {
            if (Double.isNaN(x[i]) || Double.isInfinite(x[i])) {
                return false;
            }
        }
        return true;
    }
}
