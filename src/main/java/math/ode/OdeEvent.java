package math.ode;

/**
 * A quantity along the solution whose zeros are worth stopping for.
 * <p>
 * The ball hits the ground when its height is zero, the orbit reaches apoapsis
 * when the radial velocity is zero, the reaction has finished when the
 * concentration falls to a threshold. In every case the question is not "what
 * is the state at the times I picked" but "at what time did this happen", and
 * that is a different question: the answer is not on the grid, and refining the
 * grid until it nearly is wastes most of the work.
 * <p>
 * This lives in {@code math.ode} rather than in {@code math.fun} because it is
 * not a shape but a role. The same signature would describe any scalar function
 * of a time and a state; what makes it an event is that a driver is going to
 * look for its roots along a trajectory, which is a statement about what the
 * solver does and not about the function.
 *
 * @since 1.5.3
 */
public interface OdeEvent {

    /**
     * The value of the quantity at the given time and state, whose sign changes
     * are what the driver watches for.
     * <p>
     * It has to be continuous along the solution, or a sign change need not
     * bracket a root and the time that comes back means nothing. A condition
     * that is naturally a comparison becomes a difference: not "has the height
     * fallen below zero" but "the height".
     *
     * @param t
     *            the time
     * @param y
     *            the state at that time (not modified)
     * @return the value, whose zeros are the events
     */
    double valueAt(double t, double[] y);
}
