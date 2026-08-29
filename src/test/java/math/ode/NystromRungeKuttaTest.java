package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.fun.DSecondOrderField;

/**
 * The adaptive Nystroem stepper, checked on what it claims: the order of both
 * of its solutions, the order of its interpolant, and its refusal to integrate
 * an equation it is not a method for.
 * <p>
 * The order is measured rather than asserted from the tableau, by halving the
 * step until rounding takes over. That is what polices a transcribed set of
 * coefficients, since a wrong one lowers the order without any other symptom.
 */
public class NystromRungeKuttaTest {

    private static final DSecondOrderField OSCILLATOR = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acc) {
            acc[0] = -q[0];
        }
    };

    private static final DSecondOrderField KEPLER = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acc) {
            double r = Math.sqrt(q[0] * q[0] + q[1] * q[1]);
            double r3 = r * r * r;
            acc[0] = -q[0] / r3;
            acc[1] = -q[1] / r3;
        }
    };

    /** a force that reads the velocity, which this method has no order on */
    private static final DSecondOrderField DAMPED = new DSecondOrderField() {
        @Override
        public void valueAt(double t, double[] q, double[] v, double[] acc) {
            acc[0] = -q[0] - 0.1 * v[0];
        }
    };

    /**
     * Marches the given number of equal steps and returns the worst error at the
     * end, taking the embedded solution instead of the advancing one where
     * asked.
     */
    private static double march(DSecondOrderField f, double[] y0, double t1, int steps,
            double[] want, boolean embedded) {
        int dim = y0.length;
        int dof = dim / 2;
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, f, dof);
        double[] y = y0.clone();
        double[] yOut = new double[dim];
        double[] err = new double[dim];
        double h = t1 / steps;
        for (int i = 0; i < steps; ++i) {
            s.step(i * h, y, h, yOut, err);
            for (int k = 0; k < dim; ++k) {
                // the embedded solution is the advancing one less the estimate
                y[k] = embedded ? yOut[k] - err[k] : yOut[k];
            }
        }
        double worst = 0.0;
        for (int k = 0; k < dim; ++k) {
            worst = Math.max(worst, Math.abs(y[k] - want[k]));
        }
        return worst;
    }

    private static void assertOrder(String what, DSecondOrderField f, double[] y0, double t1,
            double[] want, boolean embedded, double expected) {
        double previous = Double.NaN;
        int checked = 0;
        for (int m = 0; m < 6; ++m) {
            double e = march(f, y0, t1, 8 << m, want, embedded);
            if (m > 0 && e > 1.0e-12 && previous > 0.0) {
                assertEquals(what + ", halving " + m, expected,
                        Math.log(previous / e) / Math.log(2.0), 0.35);
                ++checked;
            }
            previous = e;
        }
        assertTrue(what + " gave only " + checked + " usable halvings", checked >= 2);
    }

    @Test
    public void testTheAdvancingSolutionIsOfOrderSix() {
        double t1 = 8.0;
        assertOrder("oscillator", OSCILLATOR, new double[] { 1.0, 0.0 }, t1,
                new double[] { Math.cos(t1), -Math.sin(t1) }, false, 6.0);
        double tk = 6.0;
        assertOrder("Kepler", KEPLER, new double[] { 1.0, 0.0, 0.0, 1.0 }, tk,
                new double[] { Math.cos(tk), Math.sin(tk), -Math.sin(tk), Math.cos(tk) }, false,
                6.0);
    }

    @Test
    public void testTheEmbeddedSolutionIsOfOrderFour() {
        double t1 = 8.0;
        assertOrder("oscillator", OSCILLATOR, new double[] { 1.0, 0.0 }, t1,
                new double[] { Math.cos(t1), -Math.sin(t1) }, true, 4.0);
    }

    /**
     * A method of order six integrates a polynomial of degree six exactly, and
     * the first degree it misses says where the order stops.
     */
    @Test
    public void testAPolynomialIsIntegratedExactlyThroughTheSixthDegree() {
        for (int k = 2; k <= 7; ++k) {
            final int kk = k;
            DSecondOrderField poly = new DSecondOrderField() {
                @Override
                public void valueAt(double t, double[] q, double[] v, double[] acc) {
                    acc[0] = kk * (kk - 1) * Math.pow(t, kk - 2);
                }
            };
            double e = march(poly, new double[] { 0.0, 0.0 }, 1.0, 1, new double[] { 1.0, kk },
                    false);
            if (k <= 6) {
                assertTrue("t^" + k + " should be exact, off by " + e, e < 1.0e-14);
            } else {
                assertTrue("t^" + k + " should not be, off by only " + e, e > 1.0e-4);
            }
        }
    }

    /**
     * The continuous extension interpolates at the order of the method, which is
     * the whole reason for carrying it rather than a Hermite polynomial over the
     * state. Measured on the nonlinear problem: the oscillator flatters it.
     */
    @Test
    public void testTheInterpolantIsOfTheMethodsOwnOrder() {
        double previous = Double.NaN;
        int checked = 0;
        for (int m = 0; m < 6; ++m) {
            double h = 0.5 / (1 << m);
            NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
            double[] y = { 1.0, 0.0, 0.0, 1.0 };
            double[] yOut = new double[4];
            double[] err = new double[4];
            s.step(0.0, y, h, yOut, err);
            double[] out = new double[4];
            double worst = 0.0;
            for (int k = 1; k < 10; ++k) {
                double theta = k / 10.0;
                s.interpolate(theta, out);
                double t = theta * h;
                // over the whole state: the position alone is not yet in the
                // asymptotic regime at the first halving and reads 7.2 there
                worst = Math.max(worst, Math.abs(out[0] - Math.cos(t)));
                worst = Math.max(worst, Math.abs(out[1] - Math.sin(t)));
                worst = Math.max(worst, Math.abs(out[2] + Math.sin(t)));
                worst = Math.max(worst, Math.abs(out[3] - Math.cos(t)));
            }
            if (m > 0 && worst > 1.0e-13 && previous > 0.0) {
                assertEquals("interpolation, halving " + m, 6.0,
                        Math.log(previous / worst) / Math.log(2.0), 0.35);
                ++checked;
            }
            previous = worst;
        }
        assertTrue("only " + checked + " usable halvings", checked >= 2);
    }

    @Test
    public void testTheEndsOfTheInterpolantAreTheEndsOfTheStep() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        double[] y = { 1.0, 0.0, 0.0, 1.0 };
        double[] yOut = new double[4];
        double[] err = new double[4];
        s.step(0.0, y, 0.3, yOut, err);
        double[] out = new double[4];
        s.interpolate(0.0, out);
        for (int k = 0; k < 4; ++k) {
            assertEquals("theta = 0, component " + k, y[k], out[k], 0.0);
        }
        s.interpolate(1.0, out);
        for (int k = 0; k < 4; ++k) {
            assertEquals("theta = 1, component " + k, yOut[k], out[k], 0.0);
        }
    }

    /**
     * The velocity the interpolant returns is the derivative of the position it
     * returns, which is what makes it a state and not two curves.
     */
    @Test
    public void testTheInterpolatedVelocityIsTheDerivativeOfThePosition() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        double h = 0.4;
        double[] yOut = new double[4];
        double[] err = new double[4];
        s.step(0.0, new double[] { 1.0, 0.0, 0.0, 1.0 }, h, yOut, err);
        double[] lo = new double[4];
        double[] hi = new double[4];
        double[] mid = new double[4];
        double d = 1.0e-6;
        for (int k = 2; k <= 7; ++k) {
            double theta = k / 10.0;
            s.interpolate(theta - d, lo);
            s.interpolate(theta + d, hi);
            s.interpolate(theta, mid);
            for (int i = 0; i < 2; ++i) {
                double slope = (hi[i] - lo[i]) / (2.0 * d * h);
                assertEquals("dq/dt against v at theta = " + theta, mid[2 + i], slope, 1.0e-7);
            }
        }
    }

    /**
     * A field that reads the velocity is refused rather than integrated at an
     * order lower than the one it is asked for.
     */
    @Test
    public void testAVelocityDependentFieldIsRefused() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, DAMPED, 1);
        try {
            s.step(0.0, new double[] { 1.0, 1.0 }, 0.1, new double[2], new double[2]);
            fail("a velocity-dependent field was accepted");
        } catch (IllegalArgumentException e) {
            String m = e.getMessage();
            assertTrue("names the form it needs, not " + m, m.contains("q'' = f(t, q)"));
            assertTrue("and says what does accept it, not " + m,
                    m.contains("SymplecticNystrom"));
        }
        // and the symplectic stepper does accept it, which is the contrast
        SymplecticNystrom other = new SymplecticNystrom(SplittingCoefficients.BLANES_MOAN_6,
                DAMPED, 1);
        other.step(0.0, new double[] { 1.0, 1.0 }, 0.1, new double[2], null);
    }

    /**
     * The probe costs one evaluation beyond the stages of the first step,
     * because the other of its two is that step's first stage.
     */
    @Test
    public void testTheVelocityProbeCostsOneEvaluation() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        double[] y = { 1.0, 0.0, 0.0, 1.0 };
        double[] yOut = new double[4];
        double[] err = new double[4];
        s.step(0.0, y, 0.1, yOut, err);
        assertEquals("six stages and one probe", 7L, s.evaluations());
        // the second step reuses the last stage of the first, so it costs five
        System.arraycopy(yOut, 0, y, 0, 4);
        s.step(0.1, y, 0.1, yOut, err);
        assertEquals("and the next step reuses the carried stage", 12L, s.evaluations());
    }

    @Test
    public void testTheDriverAdaptsTheStepSize() {
        StepController c = new StepController(1.0e-9, 1.0e-9);
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        OdeIntegrator.Result r = new OdeIntegrator(s, c).solve(0.0,
                new double[] { 1.0, 0.0, 0.0, 1.0 }, 20.0);
        assertTrue("more than one step", r.steps > 1L);
        double first = r.t[1] - r.t[0];
        boolean varied = false;
        for (int i = 2; i < r.length; ++i) {
            if (Math.abs((r.t[i] - r.t[i - 1]) - first) > 1.0e-12) {
                varied = true;
            }
        }
        assertTrue("the step size actually varied", varied);
        double[] end = r.finalState();
        assertEquals("x against cos", Math.cos(20.0), end[0], 1.0e-7);
        assertEquals("y against sin", Math.sin(20.0), end[1], 1.0e-7);
    }

    @Test
    public void testWhatTheStepperSaysAboutItself() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        assertEquals("dimension", 4, s.dimension());
        assertEquals("degrees of freedom", 2, s.degreesOfFreedom());
        assertEquals("order", 6, s.order());
        assertTrue("has an error estimate", s.hasErrorEstimate());
        assertTrue("has dense output", s.hasDenseOutput());
        assertTrue("no stiffness measure", Double.isNaN(s.stiffnessMeasure()));
        assertEquals("the tableau it was given", NystromTableau.RKN6_4, s.tableau());
        assertTrue("names itself", s.toString().contains("NystromRungeKutta"));
    }

    @Test
    public void testResetForgetsTheStepAndKeepsTheEvaluations() {
        NystromRungeKutta s = new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 2);
        s.step(0.0, new double[] { 1.0, 0.0, 0.0, 1.0 }, 0.1, new double[4], new double[4]);
        long spent = s.evaluations();
        s.reset();
        assertEquals("the evaluations survive", spent, s.evaluations());
        try {
            s.interpolate(0.5, new double[4]);
            fail("interpolation after a reset was allowed");
        } catch (IllegalStateException expected) {
            // no step has been taken
        }
        // and the velocity probe is asked again, so the next step costs seven
        s.step(0.0, new double[] { 1.0, 0.0, 0.0, 1.0 }, 0.1, new double[4], new double[4]);
        assertEquals("the probe ran again", spent + 7L, s.evaluations());
    }

    @Test
    public void testTheConstructorChecksItsArguments() {
        try {
            new NystromRungeKutta(null, KEPLER, 2);
            fail("a null tableau was accepted");
        } catch (IllegalArgumentException expected) {
            // named
        }
        try {
            new NystromRungeKutta(NystromTableau.RKN6_4, null, 2);
            fail("a null field was accepted");
        } catch (IllegalArgumentException expected) {
            // named
        }
        try {
            new NystromRungeKutta(NystromTableau.RKN6_4, KEPLER, 0);
            fail("a system of no coordinates was accepted");
        } catch (IllegalArgumentException expected) {
            // named
        }
        NystromTableau t = NystromTableau.RKN6_4;
        NystromTableau plain = new NystromTableau("no estimate", 6, 0, t.c(), t.a(), t.bbar(),
                t.b(), null, null, null, null);
        try {
            new NystromRungeKutta(plain, KEPLER, 2);
            fail("a tableau with no embedded solution was accepted");
        } catch (IllegalArgumentException e) {
            assertTrue("says why, not " + e.getMessage(),
                    e.getMessage().contains("estimates no error"));
        }
    }

    @Test
    public void testAMethodWithoutADenseExtensionFallsBackToHermite() {
        NystromTableau t = NystromTableau.RKN6_4;
        NystromTableau bare = new NystromTableau("no dense", 6, 4, t.c(), t.a(), t.bbar(), t.b(),
                t.bbarStar(), t.bStar(), null, null);
        NystromRungeKutta s = new NystromRungeKutta(bare, KEPLER, 2);
        assertFalse("says it has none", s.hasDenseOutput());
        double[] yOut = new double[4];
        s.step(0.0, new double[] { 1.0, 0.0, 0.0, 1.0 }, 0.2, yOut, new double[4]);
        double[] out = new double[4];
        s.interpolate(0.5, out);
        // still a usable answer, just at the Hermite order rather than the method's
        assertEquals("x at the middle of the step", Math.cos(0.1), out[0], 1.0e-5);
        s.interpolate(1.0, out);
        assertEquals("and the end is still exact", yOut[0], out[0], 0.0);
    }
}
