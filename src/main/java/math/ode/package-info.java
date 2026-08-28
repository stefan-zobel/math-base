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
 * extension and an error estimate, and the classical fourth order method for a
 * fixed step size. For a mechanical one {@code q'' = f(t, q)}, where the
 * position and the velocity are worth keeping apart:
 * {@link math.ode.SymplecticNystrom} over five methods, from Verlet up to
 * Blanes-Moan of order six. And for a stiff one, {@link math.ode.Rosenbrock}
 * over {@link math.ode.RosenbrockTableau#RODAS4} and two others -- linearly
 * implicit, so a step costs one Jacobian and one matrix factorization and no
 * iteration converges or fails to.
 * <p>
 * <b>Among those five, the order is not the thing to choose on.</b> Yoshida's
 * triple jump, Suzuki's five-fold composition and Blanes-Moan's
 * {@code SRKN 6b} are all of order four, and their leading error constants are
 * {@code 1.0e-01}, {@code 4.5e-03} and {@code 1.6e-05} -- four orders of
 * magnitude at the same order of accuracy, because the first travels
 * {@code 4.40} units of time to advance one and the last only {@code 1.16}.
 * Blanes-Moan reaches any given accuracy for about a fifth of the triple jump's
 * cost.
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
 * output grid, the events and the step sizes are three independent things, and
 * a caller who wants a value every millisecond does not thereby ask for a step
 * every millisecond.</li>
 * <li><b>An {@link math.ode.Event} is a question the grid cannot answer.</b>
 * When did the ball land, when did the orbit reach apoapsis, when did the
 * concentration fall through the threshold -- the answer is a root along the
 * trajectory, found by {@link math.solve.RootFinder#brentDekker} on the
 * interpolant of the step that brackets it. Measured against the analytic
 * zeros of a cosine, the times come back to {@code 1.1e-12}.</li>
 * <li><b>The tolerance bounds what one step adds, not what the run
 * accumulates.</b> The two are proportional, which is what makes the setting
 * useful: over eight decades of it, the error at the end of a particular run
 * stays at about six times {@code rtol}. The factor is a property of the
 * equation and the interval, and it is not one.</li>
 * </ul>
 * <p>
 * <b>Two questions, and a measurement for each.</b>
 * <p>
 * The first is <b>stiffness</b>. An explicit method is stable only while the
 * step size times the largest eigenvalue of the Jacobian stays inside a region
 * a few units across, so on a stiff equation the step is held down by stability
 * and not by accuracy, and no tolerance changes that. The failure does not look
 * like inaccuracy, it looks like a program that has stopped -- which is why
 * {@link math.ode.OdeIntegrator.Result#seemsStiff} exists and why the
 * exceptions say what they suspect. The answer to a {@code true} there is
 * {@link math.ode.Ode#solveStiff(math.fun.DVectorField, double, double[],
 * double, double)}.
 * <p>
 * <b>But a stiff method is not a better method.</b> It pays for a Jacobian, a
 * factorization and a back substitution per stage, and where there is no
 * stiffness to pay for it loses. Measured on van der Pol's oscillator, where
 * the stiffness is a dial: at {@code mu = 10} Dormand-Prince costs 789 field
 * evaluations against {@code RODAS4}'s 1763; at {@code mu = 100}, 3087 against
 * 5204; at {@code mu = 1000} they have crossed; at {@code mu = 1e4} it is 74205
 * against 10034; and at {@code mu = 1e5} the explicit method does not finish at
 * all. The explicit cost is proportional to {@code mu} and the implicit one
 * grows like its logarithm, by about two thousand evaluations a decade, so
 * there is a crossing and it is worth knowing which side of it a problem is
 * on. Robertson's reaction is the other side outright: {@code RODAS4} reaches
 * {@code t = 1e11} in 395 steps, and no tolerance lets an explicit method reach
 * it at all.
 * <p>
 * The second question is <b>what happens to an invariant over a long
 * time</b>, and the answer depends on which family is used. The energy of a
 * two body orbit does
 * not change; Dormand-Prince loses it steadily, by {@code 1.28e-09},
 * {@code 1.18e-08} and {@code 1.16e-07} over ten, a hundred and a thousand
 * orbits -- a factor of ten per factor of ten, and tightening the tolerance
 * lowers that line without changing its slope. A symplectic method does not
 * conserve the energy either, but its error stays inside a band: measured over
 * two hundred orbits, the band over the last twenty is the band over the first
 * twenty to three decimal places. The angular momentum it conserves outright,
 * because neither a drift nor a kick can move it.
 * <p>
 * <b>Bounded is not the same as accurate, and the choice turns on the
 * horizon.</b> At matched cost after ten orbits, Dormand-Prince is nearly three
 * orders of magnitude closer to the truth than Suzuki's method. What changes is
 * the rate: the symplectic position error grows exactly linearly with time and
 * the adaptive one grows as its square, so they cross at some ten thousand
 * orbits and past that the symplectic method is ahead and stays ahead. For a
 * short accurate answer, use the first family; for a long qualitative one, the
 * second.
 */
package math.ode;
