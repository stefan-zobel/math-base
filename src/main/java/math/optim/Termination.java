package math.optim;

/**
 * Why a search stopped.
 * <p>
 * A boolean convergence flag has to stand for eight different outcomes, and
 * the ones a caller most needs to tell apart -- a stopping rule that fired on
 * the gradient, one that fired on the objective, an exhausted budget, a line
 * search that could not move, an evaluator that intervened -- all collapse
 * into it. This enum keeps them apart. {@link #isConvergence()} reproduces the
 * old flag exactly, so nothing that read it changes meaning.
 * <p>
 * <b>What a stationary exit does and does not prove.</b>
 * {@link #GRADIENT_TOLERANCE} and {@link #ZERO_GRADIENT} report that the
 * gradient itself was small, which is a <em>necessary</em> condition for a
 * maximum and not a sufficient one. An objective whose supremum is not
 * attained approaches it along a ray on which the gradient vanishes too: an
 * unpenalized logistic likelihood on separable data reaches a gradient norm of
 * {@code 5e-04} at coefficients of size 24, and {@code 7e-13} at size 71, both
 * of them reported as {@code GRADIENT_TOLERANCE} and both of them artifacts of
 * where the tolerance happened to bite. No local rule can do better, because
 * locally there is nothing to see.
 * <p>
 * The test that does work is not a single run at all: <b>tighten the
 * tolerances and look at the parameters.</b> On a well posed problem they stay
 * where they were -- the same penalized fit moves from a norm of 3.592 to
 * 3.588 -- and on one whose maximum does not exist they grow without bound.
 * {@code LimitedMemoryBFGS.getGradientNorm()} supplies the other half of the
 * evidence, the value the rule fired on.
 *
 * @since 1.5.2
 */
public enum Termination {

    /** No search has been run yet. */
    NOT_STARTED(false, false),

    /** The gradient was exactly zero. */
    ZERO_GRADIENT(true, true),

    /**
     * The norm of the gradient fell below the gradient tolerance, so the point
     * reached is stationary to within that tolerance. Read the caution in the
     * class comment before taking it for a maximum.
     */
    GRADIENT_TOLERANCE(true, true),

    /**
     * The relative change of the objective between two iterations fell below
     * the value tolerance. This says the objective stopped moving; the
     * gradient may still be far from zero, which is the usual outcome on a
     * well conditioned problem whose optimum is a point rather than a
     * direction.
     */
    VALUE_TOLERANCE(true, false),

    /**
     * The line search could not take a step in the current direction. That is
     * the flat region near a maximum as often as it is a failure, and the
     * parameters are left at the best point found.
     */
    LINE_SEARCH_STALLED(false, false),

    /** The optimizer's own iteration budget was exhausted. */
    ITERATION_LIMIT(false, false),

    /**
     * {@code optimize(numIterations)} performed the iterations the caller
     * asked for and returned. The search is unfinished and can be continued by
     * calling {@code optimize} again.
     */
    PARTIAL_RUN(false, false),

    /** An {@link OptimizerEvaluator} stopped the search from outside. */
    EVALUATOR_ABORTED(false, false);

    private final boolean convergence;
    private final boolean stationary;

    private Termination(boolean convergence, boolean stationary) {
        this.convergence = convergence;
        this.stationary = stationary;
    }

    /**
     * Whether this outcome is what {@code isConverged()} reports as
     * convergence, that is whether one of the stopping tolerances was met
     * rather than the search being cut short.
     *
     * @return {@code true} for the three tolerance outcomes, {@code false} for
     *         a budget, a stalled line search, an aborted run and a partial one
     */
    public boolean isConvergence() {
        return convergence;
    }

    /**
     * Whether the rule that fired was about the gradient rather than about the
     * objective value. Necessary for a maximum, not sufficient -- see the
     * class comment.
     *
     * @return {@code true} for {@link #ZERO_GRADIENT} and
     *         {@link #GRADIENT_TOLERANCE}
     */
    public boolean isStationary() {
        return stationary;
    }
}
