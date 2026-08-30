package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DVectorField;
import math.fun.DiffDVectorField;

/**
 * The linearly implicit stepper, checked on the one problem where its answer
 * can be written down.
 * <p>
 * <b>The apparatus is the stability function.</b> On {@code y' = lambda y} the
 * Jacobian is a number, every stage is a division, and one step is exactly
 * {@code R(z) y} with {@code z = h lambda} and {@code R} rational. From that
 * single number, with no trajectory integrated anywhere: the order, as the
 * power of {@code z} at which {@code R} and {@code exp} part company; the order
 * of the embedded estimate, as the same for {@code R} minus the estimate;
 * L-stability, as {@code R(z)} going to zero as {@code z} goes to minus
 * infinity; and A-stability, from the determinant of the {@code 2 x 2} map of a
 * rotation, which is {@code |R(i w)|} squared. That is the round-two
 * exact-map idea again, and it is what will police a transcribed table.
 * <p>
 * <b>What it does not cover</b> is the nonlinear half of the order conditions,
 * because a linear problem cannot see them. So the order is measured a second
 * time, by halving the step on {@code y' = -y^2} against its closed form.
 */
public final class RosenbrockTest {

    private static final double SQRT_EPS = 1.5e-8;

    /** {@code y' = lambda y}, with the Jacobian written down. */
    private static DiffDVectorField linear(final double lambda) {
        return new DiffDVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = lambda * y[0];
            }

            @Override
            public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
                dfdy[0] = lambda;
                dfdt[0] = 0.0;
            }
        };
    }

    /** The rotation {@code y' = (0, -w; w, 0) y}. */
    private static DiffDVectorField rotation(final double w) {
        return new DiffDVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = -w * y[1];
                dydt[1] = w * y[0];
            }

            @Override
            public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
                dfdy[0] = 0.0;
                dfdy[1] = w;
                dfdy[2] = -w;
                dfdy[3] = 0.0;
                dfdt[0] = 0.0;
                dfdt[1] = 0.0;
            }
        };
    }

    /** {@code y' = -y^2}, {@code y(0) = 1}, exact {@code y(t) = 1 / (1 + t)}. */
    private static final DiffDVectorField SQUARE = new DiffDVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = -y[0] * y[0];
        }

        @Override
        public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
            dfdy[0] = -2.0 * y[0];
            dfdt[0] = 0.0;
        }
    };

    /**
     * {@code y' = -L (y - cos t) - sin t}, whose solution is {@code cos t} for
     * every {@code L}: the answer does not change and the difficulty of getting
     * it does. This is the family the stiffness detector was calibrated on.
     */
    private static DiffDVectorField damped(final double lambda) {
        return new DiffDVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                dydt[0] = -lambda * (y[0] - Math.cos(t)) - Math.sin(t);
            }

            @Override
            public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
                dfdy[0] = -lambda;
                dfdt[0] = -lambda * Math.sin(t) - Math.cos(t);
            }
        };
    }

    /**
     * Robertson's reaction, the canonical stiff problem: three concentrations
     * spanning ten orders of magnitude, and a linear invariant a Rosenbrock
     * method conserves outright.
     */
    private static final DiffDVectorField ROBERTSON = new DiffDVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] f) {
            f[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            f[2] = 3.0e7 * y[1] * y[1];
            f[1] = -f[0] - f[2];
        }

        @Override
        public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
            dfdy[0] = -0.04;
            dfdy[1] = 0.04;
            dfdy[2] = 0.0;
            dfdy[3] = 1.0e4 * y[2];
            dfdy[4] = -1.0e4 * y[2] - 6.0e7 * y[1];
            dfdy[5] = 6.0e7 * y[1];
            dfdy[6] = 1.0e4 * y[1];
            dfdy[7] = -1.0e4 * y[1];
            dfdy[8] = 0.0;
            dfdt[0] = 0.0;
            dfdt[1] = 0.0;
            dfdt[2] = 0.0;
        }
    };

    /** Van der Pol's oscillator, whose stiffness is the parameter. */
    private static DiffDVectorField vanDerPol(final double mu) {
        return new DiffDVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] f) {
                f[0] = y[1];
                f[1] = mu * ((1.0 - y[0] * y[0]) * y[1] - y[0]);
            }

            @Override
            public void jacobianAt(double t, double[] y, double[] dfdy, double[] dfdt) {
                dfdy[0] = 0.0;
                dfdy[1] = mu * (-2.0 * y[0] * y[1] - 1.0);
                dfdy[2] = 1.0;
                dfdy[3] = mu * (1.0 - y[0] * y[0]);
                dfdt[0] = 0.0;
                dfdt[1] = 0.0;
            }
        };
    }

    /** {@code y'' = -y} as a first order pair, exact {@code (cos t, -sin t)}. */
    private static final DVectorField OSCILLATOR = new DVectorField() {
        @Override
        public void valueAt(double t, double[] y, double[] dydt) {
            dydt[0] = y[1];
            dydt[1] = -y[0];
        }
    };

    /** One step of length one at {@code lambda = z}, so the step is {@code R(z)}. */
    private static double[] map(RosenbrockTableau tableau, double z) {
        Rosenbrock stepper = new Rosenbrock(tableau, linear(z), 1);
        double[] out = new double[1];
        double[] error = new double[1];
        stepper.step(0.0, new double[] { 1.0 }, 1.0, out, error);
        return new double[] { out[0], out[0] - error[0] };
    }

    /** Every method the package ships, all held to the same measurements. */
    private static final RosenbrockTableau[] ALL = { RosenbrockTableau.ROS2, RosenbrockTableau.RODAS4,
            RosenbrockTableau.RODAS4P };

    /**
     * The power of {@code z} at which the propagating or the embedded solution
     * parts company with {@code exp(z)}, halving until the gap reaches the
     * floor of double precision. A method with a small error constant reaches
     * that floor at a larger {@code z}, which is why the range cannot be fixed.
     */
    private static double measuredOrder(RosenbrockTableau tableau, boolean embedded) {
        double previous = Double.NaN;
        double order = Double.NaN;
        for (double z = 0.4; z > 1.0e-4; z *= 0.5) {
            double[] both = map(tableau, z);
            double gap = Math.abs((embedded ? both[1] : both[0]) - Math.exp(z));
            if (gap < 1.0e-14) {
                break;
            }
            if (!Double.isNaN(previous)) {
                order = Math.log(previous / gap) / Math.log(2.0) - 1.0;
            }
            previous = gap;
        }
        return order;
    }

    @Test
    public void testTheStabilityFunctionHasTheOrderTheTableauClaims() {
        for (RosenbrockTableau tableau : ALL) {
            assertEquals(tableau + ": R(z) parts company with exp(z) at z^(order+1)", tableau.order(),
                    measuredOrder(tableau, false), 0.05);
        }
    }

    @Test
    public void testTheEmbeddedSolutionIsExactlyOneOrderBelow() {
        for (RosenbrockTableau tableau : ALL) {
            assertEquals(tableau + ": an estimate of the same order would estimate nothing",
                    tableau.embeddedOrder(), measuredOrder(tableau, true), 0.05);
            assertEquals("and one order below is what the tableau says it is", tableau.order() - 1,
                    tableau.embeddedOrder());
        }
    }

    /**
     * The number the order alone does not give, and the one that decides
     * between two methods of the same order.
     */
    @Test
    public void testTheErrorConstantIsWhereItSettles() {
        assertEquals("ROS2", -1.375, errorConstant(RosenbrockTableau.ROS2, 1.0e-4), 0.005);
        assertEquals("RODAS4", -7.349e-04, errorConstant(RosenbrockTableau.RODAS4, 1.25e-2), 1.0e-7);
        assertEquals("RODAS4P", -8.564e-04, errorConstant(RosenbrockTableau.RODAS4P, 1.25e-2), 1.0e-7);
    }

    private static double errorConstant(RosenbrockTableau tableau, double smallest) {
        double constant = Double.NaN;
        for (double z = 0.05; z >= smallest; z *= 0.5) {
            constant = (map(tableau, z)[0] - Math.exp(z)) / Math.pow(z, tableau.order() + 1);
        }
        return constant;
    }

    /**
     * The property an explicit method cannot have at any step size, and the
     * whole reason this class exists: a component whose eigenvalue is far out
     * in the left half plane is not merely kept bounded, it is annihilated.
     */
    @Test
    public void testEveryMethodIsLStable() {
        for (RosenbrockTableau tableau : ALL) {
            double previous = 1.0;
            for (double z : new double[] { -1.0e1, -1.0e3, -1.0e6, -1.0e12 }) {
                double r = Math.abs(map(tableau, z)[0]);
                assertTrue(tableau + ": R(" + z + ") = " + r + " is not inside the unit circle", r < 1.0);
                assertTrue(tableau + ": R is not falling towards zero at z = " + z, r < previous);
                previous = r;
            }
            assertTrue(tableau + ": R(-1e12) should be at the size of its reciprocal",
                    Math.abs(map(tableau, -1.0e12)[0]) < 1.0e-10);
        }
    }

    @Test
    public void testEveryMethodIsAStable() {
        for (RosenbrockTableau tableau : ALL) {
            double worst = 0.0;
            double at = 0.0;
            for (double w = 0.01; w < 1.0e6; w *= 2.0) {
                Rosenbrock stepper = new Rosenbrock(tableau, rotation(w), 2);
                double[] first = new double[2];
                double[] second = new double[2];
                stepper.step(0.0, new double[] { 1.0, 0.0 }, 1.0, first, null);
                stepper.reset();
                stepper.step(0.0, new double[] { 0.0, 1.0 }, 1.0, second, null);
                double determinant = first[0] * second[1] - second[0] * first[1];
                if (determinant > worst) {
                    worst = determinant;
                    at = w;
                }
            }
            assertTrue(tableau + ": |R(i w)|^2 reached " + worst + " at w = " + at, worst <= 1.0);
        }
    }

    /**
     * The linear problem cannot see the nonlinear order conditions, so the
     * order is measured a second time where it can. Taken across the whole
     * range rather than between two neighbors, because a single ratio is
     * noisy where the error is already near the rounding.
     */
    @Test
    public void testTheOrderHoldsOnANonlinearProblemToo() {
        assertEquals("ROS2", 2.0, nonlinearOrder(RosenbrockTableau.ROS2, 20, 640), 0.15);
        assertEquals("RODAS4", 4.0, nonlinearOrder(RosenbrockTableau.RODAS4, 4, 128), 0.25);
        assertEquals("RODAS4P", 4.0, nonlinearOrder(RosenbrockTableau.RODAS4P, 4, 128), 0.25);
    }

    private static double nonlinearOrder(RosenbrockTableau tableau, int from, int to) {
        double coarse = fixedStepError(tableau, from);
        double fine = fixedStepError(tableau, to);
        return Math.log(coarse / fine) / Math.log((double) to / from);
    }

    private static double fixedStepError(RosenbrockTableau tableau, int steps) {
        OdeIntegrator.Result r = new OdeIntegrator(new Rosenbrock(tableau, SQUARE, 1)).solve(0.0,
                new double[] { 1.0 }, 1.0, steps);
        return Math.abs(r.finalState()[0] - 0.5);
    }

    /**
     * The weights of a stiffly accurate method are not transcribed: they are
     * the last row of the stage matrix with a one appended, and the error
     * estimate is the last stage vector by itself.
     */
    @Test
    public void testTheStifflyAccurateWeightsAreDerivedAndNotCopied() {
        for (RosenbrockTableau tableau : new RosenbrockTableau[] { RosenbrockTableau.RODAS4,
                RosenbrockTableau.RODAS4P }) {
            assertTrue(tableau + " should be stiffly accurate", tableau.isStifflyAccurate());
            double[] b = tableau.b();
            double[] lastRow = tableau.a()[tableau.stages() - 1];
            for (int j = 0; j < lastRow.length; ++j) {
                assertEquals(tableau + ": weight " + j, lastRow[j], b[j], 0.0);
            }
            assertEquals(1.0, b[tableau.stages() - 1], 0.0);
            double[] bError = tableau.bError();
            for (int j = 0; j < tableau.stages() - 1; ++j) {
                assertEquals(tableau + ": the estimate is the last stage alone", 0.0, bError[j], 0.0);
            }
            assertEquals(1.0, bError[tableau.stages() - 1], 0.0);
            assertEquals(0.25, tableau.gamma(), 0.0);
            assertEquals(6, tableau.stages());
            assertEquals(2, tableau.denseRows());
        }
    }

    /**
     * The continuous extension is built out of the stage vectors the step has
     * already computed, so it costs nothing, and it is exact where the step
     * is. Measured over five halvings, the error inside the step falls as
     * {@code h^4} where the error at its end falls as {@code h^5}: one order
     * below the step it spans, which is the same relation
     * {@link ExplicitRungeKutta} reports for Dormand-Prince.
     */
    @Test
    public void testTheDenseOutputCostsNothingAndIsOneOrderBelowTheStep() {
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.RODAS4, SQUARE, 1);
        assertTrue(stepper.hasDenseOutput());
        double[] start = { 1.0 };
        double[] end = new double[1];
        stepper.step(0.0, start, 0.25, end, null);
        long afterStep = stepper.evaluations();
        double[] out = new double[1];
        stepper.interpolate(0.0, out);
        assertEquals(start[0], out[0], 0.0);
        stepper.interpolate(1.0, out);
        assertEquals(end[0], out[0], 0.0);
        stepper.interpolate(0.5, out);
        assertEquals("no evaluation is spent on an interior value", afterStep, stepper.evaluations());
        // one step of 0.25 is coarse; measured, the middle of it is off by 8.7e-05
        assertEquals(1.0 / 1.125, out[0], 2.0e-4);

        // the two finest step sizes: the ratio is still climbing towards
        // sixteen at coarser ones, and reads 3.89 here
        double coarse = worstInsideAStep(0.05);
        double fine = worstInsideAStep(0.0125);
        double order = Math.log(coarse / fine) / Math.log(4.0);
        assertEquals("the interpolant is one order below the step", 4.0, order, 0.2);
    }

    private static double worstInsideAStep(double h) {
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.RODAS4, SQUARE, 1);
        double[] out = new double[1];
        stepper.step(0.0, new double[] { 1.0 }, h, new double[1], null);
        double worst = 0.0;
        for (double theta = 0.05; theta < 1.0; theta += 0.05) {
            stepper.interpolate(theta, out);
            worst = Math.max(worst, Math.abs(out[0] - 1.0 / (1.0 + theta * h)));
        }
        return worst;
    }

    /**
     * Six stages, six evaluations. The first of them is the field at the point
     * the step starts from, which the linearization needs anyway; the rest are
     * one per stage. No two stages of these tables fall at the same point, so
     * there is nothing here to carry over.
     */
    @Test
    public void testEachStageCostsOneEvaluation() {
        for (RosenbrockTableau tableau : ALL) {
            Rosenbrock stepper = new Rosenbrock(tableau, SQUARE, 1);
            OdeIntegrator.Result r = new OdeIntegrator(stepper).solve(0.0, new double[] { 1.0 }, 1.0, 10);
            assertEquals(tableau.toString(), 10L * tableau.stages(), r.evaluations);
            assertEquals(10L, stepper.jacobians());
            assertEquals(10L, stepper.factorizations());
        }
    }

    @Test
    public void testTheDifferencedJacobianReachesTheSameAnswerAsTheWrittenOne() {
        OdeIntegrator.Result written = new OdeIntegrator(new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1))
                .solve(0.0, new double[] { 1.0 }, 1.0, 200);
        OdeIntegrator.Result differenced = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.ROS2, (DVectorField) SQUARE, 1))
                        .solve(0.0, new double[] { 1.0 }, 1.0, 200);
        double gap = Math.abs(written.finalState()[0] - differenced.finalState()[0]);
        assertTrue("the two answers differ by " + gap, gap < SQRT_EPS);
        assertTrue("but not by nothing, the Jacobian being nonlinear here", gap > 0.0);
        assertEquals("and the differenced one pays one evaluation per component and one for the time",
                written.evaluations + 2L * 200L, differenced.evaluations);
    }

    /**
     * The Jacobian depends on the point and not on the step size, so a step
     * rejected and taken again from the same place reuses it. The count says so
     * exactly: every attempt costs one factorization, every distinct point one
     * Jacobian.
     */
    @Test
    public void testARetriedStepCostsNoSecondJacobian() {
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.ROS2, damped(1.0e4), 1);
        OdeIntegrator.Result r = new OdeIntegrator(stepper, new StepController(1.0e-8, 1.0e-8))
                .solve(0.0, new double[] { 1.0 }, 10.0);
        assertTrue("the run has to reject something for this to say anything", r.rejected > 0L);
        assertEquals("one factorization per attempt", r.steps + r.rejected, stepper.factorizations());
        assertEquals("but only one Jacobian per point", r.steps, stepper.jacobians());
        assertEquals("the first stage of a retry is free as well",
                2L + 2L * r.steps + r.rejected, r.evaluations);
    }

    @Test
    public void testTheDriverTakesItUnchangedOnBothPaths() {
        OdeIntegrator.Result fixed = new OdeIntegrator(new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1))
                .solve(0.0, new double[] { 1.0 }, 1.0, 500);
        assertEquals(501, fixed.length);
        assertEquals(500L, fixed.steps);
        assertEquals(0L, fixed.rejected);
        assertEquals("second order at a step of 0.002", 0.5, fixed.finalState()[0], 1.0e-5);

        OdeIntegrator.Result adaptive = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1), new StepController(1.0e-9, 1.0e-9))
                        .solve(0.0, new double[] { 1.0 }, 1.0);
        assertEquals(0.5, adaptive.finalState()[0], 1.0e-8);
        assertTrue(adaptive.steps > 0L);
    }

    /**
     * The reason for the whole round: on a problem whose answer is
     * {@code cos t} whatever {@code L} is, the explicit method's cost climbs
     * with {@code L} and this one's does not.
     */
    @Test
    public void testTheCostDoesNotClimbWithTheStiffness() {
        long cheapest = Long.MAX_VALUE;
        long dearest = 0L;
        long atFour = 0L;
        for (double lambda : new double[] { 1.0e2, 1.0e4, 1.0e6 }) {
            OdeIntegrator.Result r = new OdeIntegrator(
                    new Rosenbrock(RosenbrockTableau.ROS2, damped(lambda), 1),
                    new StepController(1.0e-6, 1.0e-6)).solve(0.0, new double[] { 1.0 }, 10.0);
            assertEquals("the answer is cos t whatever L is", Math.cos(10.0), r.finalState()[0], 1.0e-5);
            cheapest = Math.min(cheapest, r.evaluations);
            dearest = Math.max(dearest, r.evaluations);
            if (lambda == 1.0e4) {
                atFour = r.evaluations;
            }
        }
        // it is not monotone -- 11115, 7038 and 7263 evaluations -- and that is
        // the point: four decades of stiffness barely move the cost at all
        assertTrue("the cost ran from " + cheapest + " to " + dearest + " over four decades of L",
                dearest < 3L * cheapest);

        OdeIntegrator.Result explicit = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, damped(1.0e4), 1),
                new StepController(1.0e-6, 1.0e-6)).solve(0.0, new double[] { 1.0 }, 10.0);
        assertEquals("the explicit method gets the same answer", Math.cos(10.0),
                explicit.finalState()[0], 1.0e-5);
        assertTrue("but it should be the one paying for the stiffness: " + explicit.evaluations
                + " against " + atFour, explicit.evaluations > 10L * atFour);
        assertTrue("and it should say so", explicit.seemsStiff);
    }

    @Test
    public void testTheEndsOfTheInterpolationAreExact() {
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1);
        double[] start = { 1.0 };
        double[] end = new double[1];
        stepper.step(0.0, start, 0.25, end, null);
        double[] out = new double[1];
        stepper.interpolate(0.0, out);
        assertEquals(start[0], out[0], 0.0);
        stepper.interpolate(1.0, out);
        assertEquals(end[0], out[0], 0.0);
    }

    /**
     * A method without a continuous extension of its own falls back on the
     * cubic, whose far end derivative is the one evaluation it costs -- and
     * only if somebody asks for an interior value.
     */
    @Test
    public void testTheFallbackCostsOneEvaluationAndOnlyWhenAsked() {
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1);
        assertFalse(stepper.hasDenseOutput());
        stepper.step(0.0, new double[] { 1.0 }, 0.05, new double[1], null);
        long afterStep = stepper.evaluations();
        double[] out = new double[1];
        stepper.interpolate(0.5, out);
        assertEquals("one for the derivative at the far end", afterStep + 1L, stepper.evaluations());
        // 1 / (1 + t) at the middle of the step, to the accuracy of the step
        // itself: one step of an order two method is a coarse thing
        assertEquals(1.0 / 1.025, out[0], 1.0e-3);
        stepper.interpolate(0.75, out);
        assertEquals("and not a second time", afterStep + 1L, stepper.evaluations());
        assertEquals(1.0 / 1.0375, out[0], 1.0e-3);
    }

    @Test
    public void testAnOutputGridAndAnEventWorkWithNoChange() {
        double[] times = { 0.0, 1.0, 2.5, 7.0 };
        OdeIntegrator.Result grid = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.ROS2, OSCILLATOR, 2), new StepController(1.0e-7, 1.0e-7))
                        .solve(0.0, new double[] { 1.0, 0.0 }, 7.0, times);
        assertEquals(times.length, grid.length);
        for (int i = 0; i < times.length; ++i) {
            assertEquals(times[i], grid.t[i], 0.0);
            assertEquals(Math.cos(times[i]), grid.y[i][0], 1.0e-5);
        }

        Event crossing = new Event(new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return y[0];
            }
        });
        OdeIntegrator.Result found = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.ROS2, OSCILLATOR, 2), new StepController(1.0e-7, 1.0e-7),
                new Event[] { crossing }).solve(0.0, new double[] { 1.0, 0.0 }, 7.0);
        assertEquals(2, found.eventCount());
        assertEquals(0.5 * Math.PI, found.eventTimes[0], 1.0e-5);
        assertEquals(1.5 * Math.PI, found.eventTimes[1], 1.0e-5);
    }

    /**
     * Where the matrix is exactly singular the step gives up by handing back a
     * state that is not finite, which is what the driver already knows what to
     * do with. Nothing new was written for this.
     */
    @Test
    public void testASingularMatrixComesBackAsAStateThatIsNotFinite() {
        double h = 0.5;
        final double singular = 1.0 / (RosenbrockTableau.ROS2.gamma() * h);
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.ROS2, linear(singular), 1);
        double[] out = new double[1];
        double[] error = new double[1];
        stepper.step(0.0, new double[] { 1.0 }, h, out, error);
        assertTrue("the state should not be finite, but is " + out[0], Double.isNaN(out[0]));
        assertTrue(Double.isNaN(error[0]));
        try {
            stepper.interpolate(0.5, out);
            fail("expected a refusal: no step was taken");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("no step"));
        }
    }

    /**
     * The invariant is not conserved to the order of the method, it is
     * conserved outright. The three rates sum to zero identically, so the
     * columns of the Jacobian do too, so every stage vector sums to zero, and
     * so does the solution: over four hundred steps spanning eleven decades of
     * time the sum stays within a few times the rounding.
     */
    @Test
    public void testRobertsonConservesItsLinearInvariantOutright() {
        OdeIntegrator.Result r = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, ROBERTSON, 3),
                new StepController(1.0e-8, 1.0e-8)).solve(0.0, new double[] { 1.0, 0.0, 0.0 }, 1.0e11);
        double worst = 0.0;
        for (int i = 0; i < r.length; ++i) {
            worst = Math.max(worst, Math.abs(r.y[i][0] + r.y[i][1] + r.y[i][2] - 1.0));
        }
        assertTrue("the sum wandered by " + worst, worst < 1.0e-14);
        assertTrue("eleven decades of time in " + r.steps + " steps", r.steps < 600L);
        // everything has reacted through to the third species by then
        assertEquals(1.0, r.finalState()[2], 1.0e-6);
    }

    /**
     * The same problem, and the reason the round exists: the explicit method
     * does not merely take longer, it does not finish at all, and what it says
     * on the way out is what it should say.
     */
    @Test
    public void testTheExplicitMethodCannotDoRobertsonAtAll() {
        try {
            new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, ROBERTSON, 3),
                    new StepController(1.0e-8, 1.0e-8, 20000)).solve(0.0, new double[] { 1.0, 0.0, 0.0 },
                            1.0e5);
            fail("an explicit method has no business finishing Robertson");
        } catch (ArithmeticException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("step budget"));
            assertTrue(expected.getMessage(), expected.getMessage().contains("wants an implicit one"));
        }
    }

    /**
     * How the two families answer to stiffness, on the one problem where it is
     * a dial. The explicit method's cost is proportional to {@code mu}; this
     * one's grows like its logarithm, by about two thousand evaluations a
     * decade. They cross between {@code mu} of a hundred and a thousand, and
     * below that the explicit method is the right answer -- which is worth
     * asserting too, because a stiff solver is not a better solver.
     */
    @Test
    public void testTheCostGrowsWithTheLogarithmOfTheStiffness() {
        double[] start = { 2.0, 0.0 };
        long atTen = evaluations(new Rosenbrock(RosenbrockTableau.RODAS4, vanDerPol(10.0), 2), start);
        long atTenThousand = evaluations(
                new Rosenbrock(RosenbrockTableau.RODAS4, vanDerPol(1.0e4), 2), start);
        assertTrue("three decades of stiffness cost " + atTen + " against " + atTenThousand,
                atTenThousand < 8L * atTen);

        long explicitAtTen = evaluations(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, vanDerPol(10.0), 2), start);
        long explicitAtTenThousand = evaluations(
                new ExplicitRungeKutta(ButcherTableau.DORMAND_PRINCE_45, vanDerPol(1.0e4), 2), start);
        assertTrue("below the crossover the explicit method is the cheaper one: " + explicitAtTen
                + " against " + atTen, explicitAtTen < atTen);
        assertTrue("and above it the implicit one is, by a factor of at least five: "
                + explicitAtTenThousand + " against " + atTenThousand,
                explicitAtTenThousand > 5L * atTenThousand);
    }

    private static long evaluations(OdeStepper stepper, double[] start) {
        return new OdeIntegrator(stepper, new StepController(1.0e-7, 1.0e-7)).solve(0.0, start, 2.0).evaluations;
    }

    /**
     * The forward difference has to follow a component below one, which is what
     * a stiff system is full of. Measured on Robertson, whose second
     * concentration sits at {@code 7e-08}: with the relative step floored at
     * one this run takes 1652 steps, and with the step floored at {@code 1e-05}
     * it takes 314 against the written Jacobian's 305.
     */
    @Test
    public void testTheDifferencedJacobianFollowsAComponentFarBelowOne() {
        double[] start = { 1.0, 0.0, 0.0 };
        OdeIntegrator.Result written = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, ROBERTSON, 3),
                new StepController(1.0e-8, 1.0e-8)).solve(0.0, start, 1.0e5);
        OdeIntegrator.Result differenced = new OdeIntegrator(
                new Rosenbrock(RosenbrockTableau.RODAS4, (DVectorField) ROBERTSON, 3),
                new StepController(1.0e-8, 1.0e-8)).solve(0.0, start, 1.0e5);
        assertTrue("the differenced run took " + differenced.steps + " steps against " + written.steps,
                differenced.steps < 1.2 * written.steps);
        assertEquals("and reaches the same answer", written.finalState()[0], differenced.finalState()[0],
                1.0e-7);
    }

    /**
     * The reported norm is the largest absolute row sum of the Jacobian the
     * field wrote at the point the step was linearized at -- not at wherever the
     * caller last was, which is the part {@link SwitchingStepper} depends on.
     */
    @Test
    public void testTheJacobianNormIsTheRowSumBoundOfTheLastLinearization() {
        Rosenbrock written = new Rosenbrock(RosenbrockTableau.RODAS4, ROBERTSON, 3);
        assertTrue("nothing is linearized yet", Double.isNaN(written.jacobianNorm()));

        double[] start = { 1.0, 0.0, 0.0 };
        double[] later = { 0.5, 3.6e-5, 0.5 };
        double[] out = new double[3];
        double[] err = new double[3];

        written.step(0.0, start, 1.0e-6, out, err);
        assertEquals("the bound at the start, where two of the terms are zero",
                robertsonRowSumBound(start), written.jacobianNorm(), 0.0);

        written.step(1.0e-6, later, 1.0e-6, out, err);
        assertEquals("and it followed the linearization rather than staying behind",
                robertsonRowSumBound(later), written.jacobianNorm(), 0.0);

        written.reset();
        assertTrue("a reset forgets the linearization", Double.isNaN(written.jacobianNorm()));

        // the same reading out of a differenced Jacobian, where the components
        // that carry the norm are alive; at the start they are not, and a
        // difference quotient of zero against zero is no test of anything
        Rosenbrock differenced = new Rosenbrock(RosenbrockTableau.RODAS4, (DVectorField) ROBERTSON,
                3);
        differenced.step(0.0, later, 1.0e-6, out, err);
        double gap = Math.abs(differenced.jacobianNorm() / robertsonRowSumBound(later) - 1.0);
        assertTrue("the differenced norm is off by " + gap, gap < 1.0e-5);
    }

    /** The largest absolute row sum of Robertson's written Jacobian at {@code y}. */
    private static double robertsonRowSumBound(double[] y) {
        double[] dfdy = new double[9];
        double[] dfdt = new double[3];
        ROBERTSON.jacobianAt(0.0, y, dfdy, dfdt);
        double worst = 0.0;
        for (int i = 0; i < 3; ++i) {
            double row = 0.0;
            for (int j = 0; j < 3; ++j) {
                row += Math.abs(dfdy[j * 3 + i]);
            }
            worst = Math.max(worst, row);
        }
        return worst;
    }

    @Test
    public void testWhatTheStepperSaysAboutItself() {
        Rosenbrock written = new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 3);
        assertEquals(3, written.dimension());
        assertEquals(2, written.order());
        assertTrue(written.hasErrorEstimate());
        assertFalse(written.hasDenseOutput());
        assertTrue(Double.isNaN(written.stiffnessMeasure()));
        assertTrue(written.hasAnalyticJacobian());
        assertEquals(RosenbrockTableau.ROS2, written.tableau());
        assertTrue(written.toString(), written.toString().contains("analytic"));

        Rosenbrock differenced = new Rosenbrock(RosenbrockTableau.ROS2, (DVectorField) SQUARE, 1);
        assertFalse(differenced.hasAnalyticJacobian());
        assertTrue(differenced.toString(), differenced.toString().contains("differenced"));
    }

    @Test
    public void testTheArgumentsAreChecked() {
        try {
            new Rosenbrock(null, SQUARE, 1);
            fail("expected a refusal naming the tableau");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("tableau"));
        }
        try {
            new Rosenbrock(RosenbrockTableau.ROS2, (DVectorField) null, 1);
            fail("expected a refusal naming the field");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("field"));
        }
        try {
            new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 0);
            fail("expected a refusal naming the dimension");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("dimension"));
        }
        Rosenbrock stepper = new Rosenbrock(RosenbrockTableau.ROS2, SQUARE, 1);
        double[] y = { 1.0 };
        refusesStep(stepper, y, new double[2], 0.1, "yOut must be of length");
        refusesStep(stepper, new double[2], new double[1], 0.1, "y must be of length");
        refusesStep(stepper, y, y, 0.1, "must not be the array");
        refusesStep(stepper, y, new double[1], 0.0, "h must not be zero");
        try {
            stepper.step(0.0, y, 0.1, new double[1], new double[2]);
            fail("expected a refusal naming errOut");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("errOut"));
        }
        try {
            stepper.derivative(0.0, y, new double[2]);
            fail("expected a refusal naming dydt");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("dydt"));
        }
    }

    private static void refusesStep(Rosenbrock stepper, double[] y, double[] out, double h,
            String fragment) {
        try {
            stepper.step(0.0, y, h, out, null);
            fail("expected a refusal mentioning " + fragment);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(fragment));
        }
    }
}
