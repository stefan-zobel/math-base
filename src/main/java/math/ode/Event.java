package math.ode;

/**
 * An {@link OdeEvent} together with what the driver should do about it: which
 * way the crossing has to go, how precisely the time is wanted, how often it
 * may happen, and whether it ends the run.
 * <p>
 * <b>Only sign changes count.</b> A value that starts at zero is not an event,
 * and a value that touches zero and comes back the way it came is one event and
 * not two. That is what makes the count meaningful and what keeps a condition
 * resting on the boundary from firing at every step.
 * <p>
 * <b>The direction is usually the point.</b> A ball whose height crosses zero
 * on the way down has landed; the same crossing on the way up is the throw, and
 * a model that stops at either has stopped at the wrong one. The default is
 * {@link Direction#EITHER} because it is the honest default, not because it is
 * usually right.
 * <p>
 * Instances are immutable and hold no state from a run, so one can be shared
 * between integrators and between threads.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Brent%27s_method">Wikipedia Brent's
 * method</a>, which is what finds the crossing once a step has bracketed it.
 *
 * @since 1.5.3
 */
public final class Event {

    /**
     * Which way through zero counts.
     *
     * @since 1.5.3
     */
    public enum Direction {

        /** From negative to positive. */
        INCREASING,

        /** From positive to negative. */
        DECREASING,

        /** Either way. */
        EITHER
    }

    /** The default precision of the reported time, in the units of {@code t}. */
    public static final double DEFAULT_TOLERANCE = 1.0e-12;

    private final OdeEvent function;
    private final Direction direction;
    private final boolean terminal;
    private final double tolerance;
    private final int maxCount;

    /**
     * An event that is recorded whichever way it is crossed and does not stop
     * the run.
     *
     * @param function
     *            the quantity whose sign changes are the events
     * @throws IllegalArgumentException
     *             if {@code function} is {@code null}
     */
    public Event(OdeEvent function) {
        this(function, Direction.EITHER, false, DEFAULT_TOLERANCE, Integer.MAX_VALUE);
    }

    /**
     * An event with a direction, which does not stop the run.
     *
     * @param function
     *            the quantity whose sign changes are the events
     * @param direction
     *            which way through zero counts
     * @throws IllegalArgumentException
     *             if an argument is {@code null}
     */
    public Event(OdeEvent function, Direction direction) {
        this(function, direction, false, DEFAULT_TOLERANCE, Integer.MAX_VALUE);
    }

    /**
     * An event in full.
     *
     * @param function
     *            the quantity whose sign changes are the events
     * @param direction
     *            which way through zero counts
     * @param terminal
     *            whether the first crossing that counts ends the run there
     * @param tolerance
     *            how precisely the time is wanted, in the units of {@code t};
     *            it must be positive
     * @param maxCount
     *            how often this event may be recorded before the driver stops
     *            watching it, at least one
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or outside the range described
     */
    public Event(OdeEvent function, Direction direction, boolean terminal, double tolerance, int maxCount) {
        if (function == null) {
            throw new IllegalArgumentException("function must not be null");
        }
        if (direction == null) {
            throw new IllegalArgumentException("direction must not be null");
        }
        if (!(tolerance > 0.0) || Double.isInfinite(tolerance)) {
            throw new IllegalArgumentException("tolerance must be positive and finite, not " + tolerance);
        }
        if (maxCount < 1) {
            throw new IllegalArgumentException("maxCount must be at least 1, not " + maxCount);
        }
        this.function = function;
        this.direction = direction;
        this.terminal = terminal;
        this.tolerance = tolerance;
        this.maxCount = maxCount;
    }

    /**
     * An event that ends the run at its first crossing in the given direction.
     *
     * @param function
     *            the quantity whose sign changes are the events
     * @param direction
     *            which way through zero counts
     * @return an event that is terminal and happens at most once
     * @throws IllegalArgumentException
     *             if an argument is {@code null}
     */
    public static Event terminal(OdeEvent function, Direction direction) {
        return new Event(function, direction, true, DEFAULT_TOLERANCE, 1);
    }

    /**
     * The quantity being watched.
     *
     * @return the function given to the constructor
     */
    public OdeEvent function() {
        return function;
    }

    /**
     * Which way through zero counts.
     *
     * @return the direction, never {@code null}
     */
    public Direction direction() {
        return direction;
    }

    /**
     * Whether the run ends where this happens.
     *
     * @return {@code true} if the event is terminal
     */
    public boolean isTerminal() {
        return terminal;
    }

    /**
     * How precisely the time of the crossing is wanted.
     *
     * @return the tolerance, in the units of {@code t}
     */
    public double tolerance() {
        return tolerance;
    }

    /**
     * How often this may be recorded before the driver stops watching.
     *
     * @return the limit, {@link Integer#MAX_VALUE} for none
     */
    public int maxCount() {
        return maxCount;
    }

    /**
     * Whether a crossing whose value went from {@code before} to {@code after}
     * is one this event counts.
     *
     * @param before
     *            the sign before the crossing, negative or positive
     * @return {@code true} if the direction matches
     */
    boolean accepts(int before) {
        if (direction == Direction.EITHER) {
            return true;
        }
        return (direction == Direction.INCREASING) ? (before < 0) : (before > 0);
    }

    /**
     * The direction and whether the event is terminal.
     */
    @Override
    public String toString() {
        return "Event[" + direction + (terminal ? ", terminal" : "") + ", tolerance " + tolerance
                + (maxCount == Integer.MAX_VALUE ? "" : ", at most " + maxCount) + "]";
    }
}
