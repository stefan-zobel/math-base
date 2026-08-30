/**
 * Initial value problems for ordinary differential equations.
 * <p>
 * Given a state at one time and a rule for how fast it is changing, work out
 * the state at another time. {@link math.ode.Ode} is the shortcut for a caller
 * who does not want to choose anything; everything below it is separable
 * because the choices are real ones.
 * <p>
 * <b>The package is four things that do not know about each other.</b> A
 * {@link math.ode.ButcherTableau} is the coefficients of a method and holds no
 * code; an {@link math.ode.OdeStepper} advances the state by one step and makes
 * no decisions; a {@link math.ode.StepController} decides how long the next
 * step should be and never sees the equation; and an
 * {@link math.ode.OdeIntegrator} drives the one with the other and writes down
 * the answer. Adding a method means adding coefficients, and adding a kind of
 * method means adding a stepper -- neither touches the driver.
 * <p>
 * <b>What is here</b> is three families. For a general system
 * {@code y' = f(t, y)} that is not stiff: Dormand-Prince 5(4) with a continuous
 * extension and an error estimate, DOP853 for an answer wanted to many digits,
 * and the classical fourth order method for a fixed step size. For a mechanical
 * one {@code q'' = f(t, q)}, where the position and the velocity are worth
 * keeping apart, two methods that answer opposite questions:
 * {@link math.ode.SymplecticNystrom} over five splitting methods, from Verlet
 * up to Blanes-Moan of order six, which keeps an invariant inside a band over a
 * long run at a constant step size; and {@link math.ode.NystromRungeKutta},
 * which adapts its step size instead and is two to five times cheaper than the
 * same equation flattened onto the first order form. And for a stiff one,
 * {@link math.ode.Rosenbrock} over
 * {@link math.ode.RosenbrockTableau#RODAS4} and two others -- linearly
 * implicit, so a step costs one Jacobian and one matrix factorization and no
 * iteration converges or fails to. {@link math.ode.SwitchingStepper} is the
 * first and third of those behind one stepper, for an equation that is stiff
 * over part of its run and not over the rest.
 * <p>
 * <b>Among the symplectic five the order is not the thing to choose on.</b>
 * Yoshida's triple jump, Suzuki's five-fold composition and Blanes-Moan's
 * {@code SRKN 6b} are all of order four and their leading error constants span
 * four orders of magnitude, because the first travels {@code 4.40} units of
 * time to advance one and the last only {@code 1.16}.
 *
 * <pre>
 * DVectorField decay = (t, y, dydt) -&gt; dydt[0] = -0.5 * y[0];
 *
 * OdeIntegrator.Result r = Ode.solve(decay, 0.0, new double[] { 1.0 }, 10.0, 1.0e-10);
 * // r.finalState()[0]  the state at t = 10
 * // r.evaluations      what it cost
 * // r.seemsStiff       whether the answer is worth anything
 * </pre>
 *
 * Three things follow from the design and are worth naming, because they are
 * what makes an integrator worth reaching for over a loop that adds
 * {@code h f(t, y)}:
 * <ul>
 * <li><b>The step size is not the resolution of the answer.</b> A method with a
 * continuous extension can be asked for the state anywhere inside a step, at
 * very nearly the accuracy of the step itself and at no further cost. So the
 * output grid, the events and the step sizes are three independent things.</li>
 * <li><b>An {@link math.ode.Event} is a question the grid cannot answer.</b>
 * When did the ball land, when did the orbit reach apoapsis -- the answer is a
 * root along the trajectory, found by
 * {@link math.solve.RootFinder#brentDekker} on the interpolant of the step that
 * brackets it, and measured against analytic zeros it comes back to
 * {@code 1.1e-12}.</li>
 * <li><b>The tolerance bounds what one step adds, not what the run
 * accumulates.</b> The two are proportional, which is what makes the setting
 * useful, but the factor is a property of the equation and the interval and it
 * is not one.</li>
 * </ul>
 * <p>
 * <b>Four questions, and which call answers each.</b>
 * <p>
 * <b>Is it stiff?</b> An explicit method is stable only while the step size
 * times the largest eigenvalue of the Jacobian stays inside a region a few
 * units across, so on a stiff equation the step is held down by stability and
 * not by accuracy, and no tolerance changes that. The failure does not look
 * like inaccuracy, it looks like a program that has stopped, which is why
 * {@link math.ode.OdeIntegrator.Result#seemsStiff} exists. The answer to a
 * {@code true} there is {@link math.ode.Ode#solveStiff(math.fun.DVectorField,
 * double, double[], double, double)}. <b>But a stiff method is not a better
 * method</b> -- it pays for a Jacobian, a factorization and a back substitution
 * per stage, and where there is no stiffness to pay for it loses. On van der
 * Pol's oscillator the two cross near {@code mu = 1000}; Robertson's reaction
 * is the other side outright, where no tolerance lets an explicit method reach
 * the end at all.
 * <p>
 * <b>Does that answer hold for the whole run?</b> Often not: a relaxation
 * oscillator is stiff on its slow manifold and not during its transitions, and
 * anything with a dial in it is stiff while the dial is turned up. Choosing one
 * method for such a run is choosing which part of it to be wrong about, and
 * {@link math.ode.SwitchingStepper} -- reached through
 * {@link math.ode.Ode#solveAuto(math.fun.DVectorField, double, double[], double,
 * double)} -- declines to choose. Neither side can diagnose stiffness, so the
 * two are compared on what they cost per unit of time instead; that class
 * documents where the comparison is blind and what is done about it.
 * <p>
 * <b>How tight a tolerance?</b> This decides between the two explicit methods
 * rather than between families. Dormand-Prince needs steps proportional to
 * <code>rtol^(-1/5)</code> and DOP853 to <code>rtol^(-1/9)</code>, and a step
 * of the second costs twice a step of the first, so <b>there is a crossing and
 * it is near {@code 1e-07}</b>. Past it the advantage keeps growing, because
 * the exponents differ and the factor of two does not, and DOP853 also
 * interpolates far better -- {@code 2^8} per halving against {@code 2^5}.
 * <p>
 * <b>Over what horizon?</b> The energy of a two body orbit does not change;
 * Dormand-Prince loses it steadily, by a factor of ten per factor of ten in
 * time, and tightening the tolerance lowers that line without changing its
 * slope. A symplectic method does not conserve energy either, but its error
 * stays inside a band, and the angular momentum it conserves outright.
 * <b>Bounded is not the same as accurate.</b> At matched cost after ten orbits
 * Dormand-Prince is nearly three orders of magnitude closer to the truth; what
 * changes is the rate, since the symplectic position error grows linearly with
 * time and the adaptive one grows as its square. They cross at some ten
 * thousand orbits. For a short accurate answer use the first family, for a long
 * qualitative one the second.
 */
package math.ode;
