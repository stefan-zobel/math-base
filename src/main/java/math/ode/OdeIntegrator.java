package math.ode;

import java.util.Arrays;

import math.fun.DFunction;
import math.solve.RootFinder;

/**
 * Drives an {@link OdeStepper} from one end of an interval to the other and
 * records what it passed through.
 * <p>
 * The division of labor is the one {@link OdeStepper} describes: a stepper
 * knows how to advance the state by one step and how to estimate the error of
 * doing so, and everything that is a decision -- how long the step should be,
 * whether to keep it, what a state that came back non-finite means, where to
 * write down an answer -- happens here. That is why the same driver runs an
 * explicit Runge-Kutta method and, in time, a symplectic one: it never learns
 * which of the two it is holding.
 * <p>
 * <b>Fixed steps are not a degenerate case of adaptive ones.</b> The forms that
 * take a step count make every step that length and consult no
 * {@link StepController} even if one was given, which is what a method with no
 * error estimate can be driven at, and what a symplectic method has to be
 * driven at: the long term behavior such a method is chosen for is a property
 * of a constant step size and does not survive one that varies. The forms
 * without a step count are the adaptive ones and need a controller.
 * <p>
 * <b>Time is accumulated, not multiplied out.</b> Each step starts at the time
 * the one before it ended at, rather than at {@code t0 + i h}, and the two are
 * not the same number: the second differs from the first in the last place or
 * two, and a stepper carrying a stage over from the previous step recognizes
 * the point it is handed by comparing it. Multiplying out costs four of the
 * nineteen carried stages over twenty steps of Dormand-Prince, for no gain --
 * the price is that a fixed step run reaches {@code t1} only to within the
 * rounding accumulated on the way, which is far below the error of the
 * integration itself. An adaptive run lands on {@code t1} exactly, because its
 * last step is cut to fit.
 * <p>
 * <b>There is no stepping mode here</b>, because {@link OdeStepper} is one.
 * Where {@link math.ts.KalmanFilter} offers a batch form and an instance that
 * takes one observation at a time, the two live in different classes here, and
 * the batch form is a recording loop over the other.
 * <p>
 * <b>Output does not have to be the steps.</b> Handed a list of times, the
 * driver puts the answer at those times and nowhere else, reaching inside the
 * steps through {@link OdeStepper#interpolate(double, double[])}. A time that
 * falls on a step boundary is the state at that boundary exactly, so a grid
 * that happens to line up with the steps cannot disagree with them -- which
 * matters more in the adaptive case, where the steps are not where a caller
 * would have put them and asking for a grid is the normal thing to do.
 * <p>
 * <b>And it does not have to be a time at all.</b> An {@link Event} is a
 * quantity along the solution whose sign change ends a step bracketing a root;
 * {@link RootFinder#brentDekker} then finds it on the interpolant, at a
 * precision that has nothing to do with the step size. That is the whole reason
 * the continuous extension is not an afterthought here: without it, "when did
 * the ball land" could only be answered by shortening every step until the
 * grid nearly contained the answer.
 * <p>
 * An integrator is as shareable between threads as the stepper it drives, which
 * is to say not at all.
 * <p>
 * <b>See</b> <a href=
 * "https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations">
 * Wikipedia numerical methods for ordinary differential equations</a>.
 *
 * @since 1.5.3
 */
public final class OdeIntegrator {

    /** How many steps in a row have to look stiff before the run says so. */
    private static final int STIFF_STREAK = 15;

    /** How many calm steps in a row wipe the suspicion out again. */
    private static final int CALM_STREAK = 6;

    private static final Event[] NO_EVENTS = new Event[0];
    private static final double[] NO_TIMES = new double[0];
    private static final int[] NO_INDICES = new int[0];
    private static final double[][] NO_STATES = new double[0][];

    private final OdeStepper stepper;
    private final StepController controller;
    private final Event[] events;

    /**
     * An integrator for fixed steps, whose length its caller decides.
     *
     * @param stepper
     *            the method to drive
     * @throws IllegalArgumentException
     *             if {@code stepper} is {@code null}
     */
    public OdeIntegrator(OdeStepper stepper) {
        this(stepper, null, null);
    }

    /**
     * An integrator that chooses its own step sizes.
     *
     * @param stepper
     *            the method to drive, which must estimate the error of its own
     *            steps
     * @param controller
     *            the policy that turns those estimates into step sizes, or
     *            {@code null} for an integrator that only takes fixed steps
     * @throws IllegalArgumentException
     *             if {@code stepper} is {@code null}, or a controller is given
     *             for a method with no error estimate
     */
    public OdeIntegrator(OdeStepper stepper, StepController controller) {
        this(stepper, controller, null);
    }

    /**
     * An integrator that watches for events, at a fixed step size if no
     * controller is given.
     *
     * @param stepper
     *            the method to drive
     * @param controller
     *            the policy that chooses the step sizes, or {@code null} for
     *            fixed steps only
     * @param events
     *            the quantities to watch along the solution, or {@code null}
     *            for none; the array is copied
     * @throws IllegalArgumentException
     *             if {@code stepper} is {@code null}, a controller is given for
     *             a method with no error estimate, or an entry of
     *             {@code events} is {@code null}
     */
    public OdeIntegrator(OdeStepper stepper, StepController controller, Event[] events) {
        if (stepper == null) {
            throw new IllegalArgumentException("stepper must not be null");
        }
        if (controller != null && !stepper.hasErrorEstimate()) {
            throw new IllegalArgumentException(
                    "a step controller needs a method that estimates the error of its own steps, and "
                            + stepper + " does not");
        }
        if (events != null) {
            for (int i = 0; i < events.length; ++i) {
                if (events[i] == null) {
                    throw new IllegalArgumentException("events[" + i + "] must not be null");
                }
            }
        }
        this.stepper = stepper;
        this.controller = controller;
        this.events = (events == null || events.length == 0) ? NO_EVENTS : events.clone();
    }

    /**
     * The method this integrator drives.
     *
     * @return the stepper given to the constructor
     */
    public OdeStepper stepper() {
        return stepper;
    }

    /**
     * The policy that chooses the step sizes.
     *
     * @return the controller, or {@code null} for a fixed step integrator
     */
    public StepController controller() {
        return controller;
    }

    /**
     * The quantities being watched along the solution.
     *
     * @return a fresh copy of the events, empty if there are none
     */
    public Event[] events() {
        return events.clone();
    }

    /**
     * Integrates from {@code t0} to {@code t1} in equally long steps, keeping
     * the state at the end of every one of them.
     *
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state, of the stepper's dimension (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param steps
     *            how many equal steps to take, at least one
     * @return the state at {@code steps + 1} times, the first being {@code t0},
     *         or fewer if a terminal event stopped the run
     * @throws IllegalArgumentException
     *             if an argument is out of shape or {@code t1} equals
     *             {@code t0}
     * @throws ArithmeticException
     *             if the state stops being finite
     */
    public Result solve(double t0, double[] y0, double t1, int steps) {
        return runFixed(t0, y0, t1, steps, null);
    }

    /**
     * Integrates from {@code t0} to {@code t1} in equally long steps, keeping
     * the state at the given times and at no others.
     *
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state, of the stepper's dimension (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param steps
     *            how many equal steps to take, at least one
     * @param times
     *            the times to report the state at, running from {@code t0}
     *            towards {@code t1} without turning back and not passing the
     *            far end; they need not line up with the steps
     * @return the state at the requested times, in the order they were given,
     *         cut short if a terminal event stopped the run
     * @throws IllegalArgumentException
     *             if an argument is out of shape, {@code t1} equals {@code t0},
     *             or a time leaves the interval or the direction
     * @throws ArithmeticException
     *             if the state stops being finite
     */
    public Result solve(double t0, double[] y0, double t1, int steps, double[] times) {
        if (times == null) {
            throw new IllegalArgumentException("times must not be null");
        }
        return runFixed(t0, y0, t1, steps, times);
    }

    /**
     * Integrates from {@code t0} to {@code t1}, choosing the step sizes, and
     * keeps the state at the end of every step that was accepted.
     * <p>
     * Where those steps fall is the controller's business and not something to
     * rely on; {@link #solve(double, double[], double, double[])} is the form
     * for a caller who needs the answer somewhere in particular.
     *
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state, of the stepper's dimension (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @return the state at {@code t0} and at the end of every accepted step,
     *         the last of which is {@code t1} exactly unless a terminal event
     *         stopped the run first
     * @throws IllegalStateException
     *             if this integrator has no step controller
     * @throws IllegalArgumentException
     *             if an argument is out of shape or {@code t1} equals
     *             {@code t0}
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public Result solve(double t0, double[] y0, double t1) {
        return runAdaptive(t0, y0, t1, null);
    }

    /**
     * Integrates from {@code t0} to {@code t1}, choosing the step sizes, and
     * keeps the state at the given times and at no others.
     *
     * @param t0
     *            the time the initial state belongs to
     * @param y0
     *            the initial state, of the stepper's dimension (not modified)
     * @param t1
     *            the time to reach, which may lie before {@code t0}
     * @param times
     *            the times to report the state at, running from {@code t0}
     *            towards {@code t1} without turning back and not passing
     *            {@code t1}
     * @return the state at the requested times, in the order they were given,
     *         cut short if a terminal event stopped the run
     * @throws IllegalStateException
     *             if this integrator has no step controller
     * @throws IllegalArgumentException
     *             if an argument is out of shape, {@code t1} equals {@code t0},
     *             or a time leaves the interval or the direction
     * @throws ArithmeticException
     *             if the step size collapses or the step budget runs out
     */
    public Result solve(double t0, double[] y0, double t1, double[] times) {
        if (times == null) {
            throw new IllegalArgumentException("times must not be null");
        }
        return runAdaptive(t0, y0, t1, times);
    }

    private Result runFixed(double t0, double[] y0, double t1, int steps, double[] times) {
        int n = check(t0, y0, t1);
        if (steps < 1) {
            throw new IllegalArgumentException("steps must be at least 1, not " + steps);
        }
        boolean forward = t1 > t0;
        double h = (t1 - t0) / steps;
        if (times != null) {
            // the accumulated end may overshoot t1 by a place or two, and the
            // times this run reports have to be times it would also accept
            double reached = t0;
            for (int i = 0; i < steps; ++i) {
                reached += h;
            }
            checkTimes(times, t0, forward ? Math.max(t1, reached) : Math.min(t1, reached), forward);
        }

        stepper.reset();
        long before = stepper.evaluations();
        Recorder recorder = new Recorder(times, forward, (times == null) ? steps + 1 : times.length);
        Watch watch = new Watch(events, n);
        Stiffness stiffness = new Stiffness(stepper);
        double[] y = y0.clone();
        double[] yNext = new double[n];
        recorder.initial(t0, y);
        watch.initial(t0, y);

        double ta = t0;
        long taken = 0L;
        for (int i = 0; i < steps; ++i) {
            stepper.step(ta, y, h, yNext, null);
            if (!allFinite(yNext)) {
                throw new ArithmeticException("the state at step " + (i + 1) + ", t = " + (ta + h)
                        + ", is not finite; the step size " + h + " is too large for this equation"
                        + stiffness.verdict());
            }
            ++taken;
            stiffness.afterStep(stepper);
            double tb = ta + h;
            if (watch.afterStep(stepper, ta, h, tb, yNext)) {
                recorder.stopAt(stepper, ta, h, watch.stopTime(), watch.stopState());
                return recorder.finish(taken, 0L, stepper.evaluations() - before, watch, stiffness);
            }
            recorder.afterStep(stepper, ta, h, tb, yNext, i == steps - 1);
            double[] swap = y;
            y = yNext;
            yNext = swap;
            ta = tb;
        }
        return recorder.finish(taken, 0L, stepper.evaluations() - before, watch, stiffness);
    }

    private Result runAdaptive(double t0, double[] y0, double t1, double[] times) {
        if (controller == null) {
            throw new IllegalStateException(
                    "this integrator has no step controller, so the number of steps has to be given");
        }
        int n = check(t0, y0, t1);
        boolean forward = t1 > t0;
        if (times != null) {
            checkTimes(times, t0, t1, forward);
        }

        stepper.reset();
        long before = stepper.evaluations();
        Recorder recorder = new Recorder(times, forward, (times == null) ? 64 : times.length);
        Watch watch = new Watch(events, n);
        Stiffness stiffness = new Stiffness(stepper);
        double[] y = y0.clone();
        double[] yNext = new double[n];
        double[] error = new double[n];
        recorder.initial(t0, y);
        watch.initial(t0, y);

        int order = stepper.order();
        int budget = controller.maxSteps();
        double maxStep = controller.maxStepSize();
        double h = controller.initialStep(stepper, t0, y, t1);
        if (!isFinite(h) || h == 0.0) {
            throw new ArithmeticException("no first step size could be found at t0 = " + t0
                    + "; the right hand side is not finite there");
        }
        double previousError = StepController.ERROR_FLOOR;
        long accepted = 0L;
        long rejected = 0L;
        double t = t0;

        while (forward ? (t < t1) : (t > t1)) {
            if (accepted + rejected >= budget) {
                throw new ArithmeticException("the step budget of " + budget + " ran out at t = " + t
                        + ", with " + rejected + " steps rejected" + stiffness.verdict());
            }
            // the resolution of the time the step starts from, and not of the
            // interval: a stiff run to t1 = 1e11 begins at t = 0 with a step of
            // 1e-04, which is far below the ulp of t1 and perfectly able to
            // make progress. At t = 0 itself the guard cannot fire and does not
            // need to, since the next step starts somewhere else; the step
            // budget above is what catches a run that crawls
            if (Math.abs(h) <= 16.0 * Math.ulp(t)) {
                throw new ArithmeticException("the step size collapsed to " + h + " at t = " + t
                        + ", which is the resolution of the time itself; the solution may end there,"
                        + " or the tolerance may be below what double precision can carry"
                        + stiffness.verdict());
            }
            boolean last = forward ? (t + h >= t1) : (t + h <= t1);
            double trial = last ? (t1 - t) : h;

            stepper.step(t, y, trial, yNext, error);
            // a method that estimates its error twice, at two orders, hands the
            // second one over here and the controller combines them, because
            // combining takes the tolerances and those belong to the controller
            double[] second = stepper.secondaryError();
            boolean finite = allFinite(yNext) && allFinite(error) && (second == null || allFinite(second));
            double scaled = Double.NaN;
            if (finite) {
                scaled = (second == null) ? controller.errorNorm(error, y, yNext)
                        : controller.errorNorm(error, second, y, yNext);
            }

            if (finite && scaled <= 1.0) {
                ++accepted;
                double tb = last ? t1 : (t + trial);
                stiffness.afterStep(stepper);
                if (watch.afterStep(stepper, t, trial, tb, yNext)) {
                    recorder.stopAt(stepper, t, trial, watch.stopTime(), watch.stopState());
                    return recorder.finish(accepted, rejected, stepper.evaluations() - before, watch,
                            stiffness);
                }
                recorder.afterStep(stepper, t, trial, tb, yNext, last);
                double[] swap = y;
                y = yNext;
                yNext = swap;
                t = tb;
                double factor = controller.scale(scaled, previousError, order);
                previousError = Math.max(scaled, StepController.ERROR_FLOOR);
                h = bound(trial * factor, maxStep);
            } else {
                ++rejected;
                h = bound(trial * controller.scaleAfterRejection(scaled, order), maxStep);
            }
        }
        return recorder.finish(accepted, rejected, stepper.evaluations() - before, watch, stiffness);
    }

    private int check(double t0, double[] y0, double t1) {
        int n = stepper.dimension();
        if (y0 == null || y0.length != n) {
            throw new IllegalArgumentException("y0 must be of length " + n);
        }
        if (!isFinite(t0) || !isFinite(t1)) {
            throw new IllegalArgumentException("t0 and t1 must be finite");
        }
        if (t1 == t0) {
            throw new IllegalArgumentException("t1 must differ from t0, both being " + t0);
        }
        if (!allFinite(y0)) {
            throw new IllegalArgumentException("y0 must be finite");
        }
        return n;
    }

    private static double bound(double h, double maxStep) {
        if (Math.abs(h) > maxStep) {
            return (h < 0.0) ? -maxStep : maxStep;
        }
        return h;
    }

    private static void checkTimes(double[] times, double t0, double t1, boolean forward) {
        double previous = t0;
        for (int i = 0; i < times.length; ++i) {
            double time = times[i];
            if (!isFinite(time)) {
                throw new IllegalArgumentException("times[" + i + "] is not finite");
            }
            if (forward ? (time < t0 || time > t1) : (time > t0 || time < t1)) {
                throw new IllegalArgumentException("times[" + i + "] is " + time
                        + ", which lies outside the interval from " + t0 + " to " + t1);
            }
            if (forward ? (time < previous) : (time > previous)) {
                throw new IllegalArgumentException("times[" + i + "] is " + time + ", which turns back from "
                        + previous);
            }
            previous = time;
        }
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static boolean allFinite(double[] x) {
        for (int i = 0; i < x.length; ++i) {
            if (!isFinite(x[i])) {
                return false;
            }
        }
        return true;
    }

    /**
     * Keeps Hairer's non-stiffness test over a run: the largest
     * <code>|h| |lambda|</code> seen, and whether enough steps in a row have
     * looked stiff to say so out loud.
     */
    private static final class Stiffness {

        private final double threshold;
        private double worst = Double.NaN;
        private int stiffStreak;
        private int calmStreak;
        private boolean seems;

        Stiffness(OdeStepper stepper) {
            // how far the method's own stability region reaches, which is not
            // the same number for a fifth order pair and for an eighth order one
            this.threshold = stepper.stiffnessThreshold();
        }

        void afterStep(OdeStepper stepper) {
            double measure = stepper.stiffnessMeasure();
            if (Double.isNaN(measure)) {
                return;
            }
            if (Double.isNaN(worst) || measure > worst) {
                worst = measure;
            }
            if (measure > threshold) {
                calmStreak = 0;
                ++stiffStreak;
                if (stiffStreak >= STIFF_STREAK) {
                    seems = true;
                }
            } else {
                ++calmStreak;
                if (calmStreak >= CALM_STREAK) {
                    stiffStreak = 0;
                    calmStreak = 0;
                }
            }
        }

        String verdict() {
            if (seems) {
                return "; the equation is stiff, at " + worst + " times the stability limit of an"
                        + " explicit method, and wants an implicit one";
            }
            if (Double.isNaN(worst)) {
                return "";
            }
            return "; stiffness was measured at " + worst + " against a threshold of " + threshold;
        }
    }

    /**
     * Watches the events over a run: the sign each was last seen with, how
     * often it has fired, and where a terminal one stopped things.
     * <p>
     * Only a change of sign counts, and the sign remembered is the last
     * <em>non-zero</em> one, which is what keeps a value that touches zero from
     * being counted once on the way in and again on the way out.
     */
    private static final class Watch {

        private final Event[] events;
        private final int[] lastSign;
        private final int[] count;
        private final double[] buffer;
        private final double[] thetas;
        private final int[] order;

        private double[] eventT;
        private int[] eventIndex;
        private double[][] eventY;
        private int written;

        private boolean stopped;
        private double stopTime;
        private double[] stopState;

        Watch(Event[] events, int n) {
            this.events = events;
            this.lastSign = new int[events.length];
            this.count = new int[events.length];
            this.buffer = (events.length == 0) ? null : new double[n];
            this.thetas = new double[events.length];
            this.order = new int[events.length];
            this.eventT = NO_TIMES;
            this.eventIndex = NO_INDICES;
            this.eventY = NO_STATES;
        }

        void initial(double t0, double[] y0) {
            for (int i = 0; i < events.length; ++i) {
                double value = events[i].function().valueAt(t0, y0);
                lastSign[i] = (value > 0.0) ? 1 : ((value < 0.0) ? -1 : 0);
            }
        }

        /**
         * Looks at the step that has just been accepted and records whatever
         * crossed. Returns whether a terminal event ended the run.
         */
        boolean afterStep(final OdeStepper stepper, final double ta, final double h, double tb,
                double[] yb) {
            int found = 0;
            for (int i = 0; i < events.length; ++i) {
                if (count[i] >= events[i].maxCount()) {
                    continue;
                }
                double value = events[i].function().valueAt(tb, yb);
                int before = lastSign[i];
                if (value == 0.0) {
                    lastSign[i] = 0;
                    if (before != 0 && events[i].accepts(before)) {
                        thetas[found] = 1.0;
                        order[found] = i;
                        ++found;
                    }
                    continue;
                }
                int now = (value > 0.0) ? 1 : -1;
                lastSign[i] = now;
                if (before == 0 || now == before || !events[i].accepts(before)) {
                    continue;
                }
                thetas[found] = root(stepper, ta, h, events[i]);
                order[found] = i;
                ++found;
            }
            if (found == 0) {
                return false;
            }
            sortByTheta(found);
            for (int k = 0; k < found; ++k) {
                int i = order[k];
                double theta = thetas[k];
                double time = (theta == 1.0) ? tb : (ta + theta * h);
                stepper.interpolate(theta, buffer);
                add(time, i, buffer);
                ++count[i];
                if (events[i].isTerminal()) {
                    stopped = true;
                    stopTime = time;
                    stopState = buffer.clone();
                    return true;
                }
            }
            return false;
        }

        private double root(final OdeStepper stepper, final double ta, final double h, Event event) {
            final OdeEvent function = event.function();
            final double[] scratch = buffer;
            DFunction along = new DFunction() {
                @Override
                public double apply(double theta) {
                    stepper.interpolate(theta, scratch);
                    return function.valueAt(ta + theta * h, scratch);
                }
            };
            double theta = RootFinder.brentDekker(0.0, 1.0, along, event.tolerance() / Math.abs(h));
            if (theta < 0.0) {
                return 0.0;
            }
            return (theta > 1.0) ? 1.0 : theta;
        }

        private void sortByTheta(int found) {
            for (int a = 1; a < found; ++a) {
                double theta = thetas[a];
                int index = order[a];
                int b = a - 1;
                while (b >= 0 && thetas[b] > theta) {
                    thetas[b + 1] = thetas[b];
                    order[b + 1] = order[b];
                    --b;
                }
                thetas[b + 1] = theta;
                order[b + 1] = index;
            }
        }

        private void add(double time, int index, double[] state) {
            if (written == eventT.length) {
                int room = Math.max(4, eventT.length * 2);
                eventT = Arrays.copyOf(eventT, room);
                eventIndex = Arrays.copyOf(eventIndex, room);
                eventY = Arrays.copyOf(eventY, room);
            }
            eventT[written] = time;
            eventIndex[written] = index;
            eventY[written] = state.clone();
            ++written;
        }

        boolean stopped() {
            return stopped;
        }

        double stopTime() {
            return stopTime;
        }

        double[] stopState() {
            return stopState;
        }
    }

    /**
     * Collects the answer, either at every step or at the times that were asked
     * for, and is the one place where a step becomes a row in the result.
     */
    private static final class Recorder {

        private final double[] times;
        private final boolean forward;
        private double[] t;
        private double[][] y;
        private double[] point;
        private int written;
        private int next;

        Recorder(double[] times, boolean forward, int capacity) {
            this.times = times;
            this.forward = forward;
            int room = Math.max(1, capacity);
            this.t = new double[room];
            this.y = new double[room][];
        }

        void initial(double t0, double[] y0) {
            if (times == null) {
                add(t0, y0);
            }
        }

        void afterStep(OdeStepper stepper, double ta, double h, double tb, double[] yNext, boolean last) {
            if (times == null) {
                add(tb, yNext);
                return;
            }
            while (next < times.length) {
                double want = times[next];
                if (!last && (forward ? (want > tb) : (want < tb))) {
                    return;
                }
                // the ends of a step are recognized rather than divided for:
                // (tb - ta) / h is not exactly one, and a requested time that
                // is a step boundary has to be that step
                double theta;
                if (want == tb) {
                    theta = 1.0;
                } else if (want == ta) {
                    theta = 0.0;
                } else {
                    theta = (want - ta) / h;
                    if (theta < 0.0) {
                        theta = 0.0;
                    } else if (theta > 1.0) {
                        theta = 1.0;
                    }
                }
                if (point == null) {
                    point = new double[yNext.length];
                }
                stepper.interpolate(theta, point);
                add(want, point);
                ++next;
            }
        }

        /**
         * A terminal event cut the step short: report whatever was asked for up
         * to the event, and then the event itself, so that the last state in
         * the result is the one the run stopped at.
         */
        void stopAt(OdeStepper stepper, double ta, double h, double time, double[] state) {
            afterStep(stepper, ta, h, time, state, false);
            if (written == 0 || t[written - 1] != time) {
                add(time, state);
            }
        }

        Result finish(long steps, long rejected, long evaluations, Watch watch, Stiffness stiffness) {
            return new Result(Arrays.copyOf(t, written), Arrays.copyOf(y, written), steps, rejected,
                    evaluations, Arrays.copyOf(watch.eventT, watch.written),
                    Arrays.copyOf(watch.eventIndex, watch.written),
                    Arrays.copyOf(watch.eventY, watch.written), watch.stopped(), stiffness.seems,
                    stiffness.worst);
        }

        private void add(double time, double[] state) {
            if (written == t.length) {
                int room = t.length * 2;
                t = Arrays.copyOf(t, room);
                y = Arrays.copyOf(y, room);
            }
            t[written] = time;
            y[written] = state.clone();
            ++written;
        }
    }

    /**
     * What one run of an {@link OdeIntegrator} produced, as public final fields
     * in the manner of {@link math.ts.KalmanFilter.Result}.
     * <p>
     * The arrays are the result's own and are handed out live rather than
     * copied; a caller that writes into them is changing what it was given.
     *
     * @since 1.5.3
     */
    public static final class Result {

        /** How many times the state was recorded at. */
        public final int length;

        /** The times, one per recorded state, in the order they were reached. */
        public final double[] t;

        /** The states, {@code y[i]} belonging to {@code t[i]}. */
        public final double[][] y;

        /** How many steps were accepted. */
        public final long steps;

        /**
         * How many steps were computed and thrown away because their estimated
         * error was too large, which is zero at a fixed step size.
         */
        public final long rejected;

        /**
         * How often this run evaluated the right hand side, which is the honest
         * measure of what it cost.
         */
        public final long evaluations;

        /** The times the events happened at, in the order they happened. */
        public final double[] eventTimes;

        /**
         * Which event each occurrence belongs to, as an index into the array
         * given to the integrator.
         */
        public final int[] eventIndices;

        /** The state at each occurrence, from the interpolant of its step. */
        public final double[][] eventStates;

        /** Whether a terminal event ended the run before {@code t1}. */
        public final boolean stoppedByEvent;

        /**
         * Whether the equation held the step size down by stability rather than
         * by accuracy for long enough to say so.
         * <p>
         * When this is set, an explicit method is the wrong tool and no
         * tolerance will make it the right one: the step is short because a
         * longer one would be unstable, not because a longer one would be
         * inaccurate. What it is telling a caller to reach for is an implicit
         * method, and
         * {@link Ode#solveStiff(math.fun.DVectorField, double, double[], double, double)}
         * is one.
         */
        public final boolean seemsStiff;

        /**
         * The largest <code>|h| |lambda|</code> seen over the run, against
         * {@link OdeStepper#stiffnessThreshold()} -- {@code 3.25} for
         * Dormand-Prince and {@code 6.1} for DOP853, since the two stability
         * regions do not reach equally far. {@link Double#NaN} if the method
         * cannot estimate it, as the classical one cannot.
         */
        public final double stiffness;

        Result(double[] t, double[][] y, long steps, long rejected, long evaluations, double[] eventTimes,
                int[] eventIndices, double[][] eventStates, boolean stoppedByEvent, boolean seemsStiff,
                double stiffness) {
            this.length = t.length;
            this.t = t;
            this.y = y;
            this.steps = steps;
            this.rejected = rejected;
            this.evaluations = evaluations;
            this.eventTimes = eventTimes;
            this.eventIndices = eventIndices;
            this.eventStates = eventStates;
            this.stoppedByEvent = stoppedByEvent;
            this.seemsStiff = seemsStiff;
            this.stiffness = stiffness;
        }

        /**
         * How many events happened.
         *
         * @return the number of occurrences, over all events together
         */
        public int eventCount() {
            return eventTimes.length;
        }

        /**
         * The state at the end of the run.
         *
         * @return {@code y[length - 1]}, live rather than copied
         */
        public double[] finalState() {
            return y[length - 1];
        }

        /**
         * The time at the end of the run.
         *
         * @return {@code t[length - 1]}
         */
        public double finalTime() {
            return t[length - 1];
        }

        /**
         * The share of the steps that were computed and thrown away, which says
         * whether the step size was being chosen well.
         *
         * @return the rejected steps over all of them, zero if there were none
         */
        public double rejectionRate() {
            long all = steps + rejected;
            return (all == 0L) ? 0.0 : rejected / (double) all;
        }

        /**
         * How far the run got, at what cost.
         */
        @Override
        public String toString() {
            return "Result[" + length + " points from " + t[0] + " to " + t[length - 1] + ", " + steps
                    + " steps, " + rejected + " rejected, " + evaluations + " evaluations, "
                    + eventTimes.length + " events" + (stoppedByEvent ? ", stopped by one" : "")
                    + (seemsStiff ? ", stiff" : "") + "]";
        }
    }
}
