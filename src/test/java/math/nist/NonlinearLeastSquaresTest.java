package math.nist;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import math.optim.BoundedLevenbergMarquardt;
import math.optim.LevenbergMarquardt;

/**
 * The nonlinear half of the certification: five reference sets, each fitted
 * from <b>both</b> starting points NIST prescribes. The far start is the
 * interesting one -- a solver that only works from near the answer is only half
 * a solver, and it is where packages are usually found out.
 * <p>
 * As on the linear side the floors were measured rather than chosen, and they
 * exist so that a later change cannot quietly lose digits.
 */
public class NonlinearLeastSquaresTest {

    /** Digits to reach from either start: parameters, standard deviations, residual sum of squares. */
    private static final double[][] FLOORS = {
            { 9.0, 7.0, 10.0 }, // Misra1a
            { 5.0, 4.0, 11.0 }, // Chwirut1
            { 4.0, 3.0, 8.0 },  // Thurber
            { 3.0, 3.0, 8.0 },  // MGH09
            { 8.0, 9.0, 11.0 }  // MGH10
    };

    @Test
    public void theCollectionIsIntact() {
        StRD.NonlinearSet[] sets = StRD.nonlinear();
        assertEquals(5, sets.length);
        assertEquals(FLOORS.length, sets.length);
        for (int k = 0; k < sets.length; ++k) {
            StRD.NonlinearSet set = sets[k];
            assertEquals(set.name + " predictors", set.observations, set.x().length);
            assertEquals(set.name + " responses", set.observations, set.y().length);
            assertEquals(set.name + " far start", set.parameters, set.start(1).length);
            assertEquals(set.name + " near start", set.parameters, set.start(2).length);
            assertEquals(set.name + " certified parameters", set.parameters, set.certifiedBeta().length);
            assertTrue(set.name + " residual sum of squares", set.residualSumOfSquares > 0.0);
            assertNotNull(set.name + " model", Models.of(set));
        }
        assertEquals("Misra1a", sets[0].name);
        assertEquals("MGH10", sets[4].name);
    }

    /** Every set from every prescribed start, against the digits it reached. */
    @Test
    public void everySetIsFittedFromBothStarts() {
        StRD.NonlinearSet[] sets = StRD.nonlinear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.NonlinearSet set = sets[k];
            for (int start = 1; start <= 2; ++start) {
                LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(Models.of(set), set.start(start),
                        set.observations);
                String where = set.name + " from start " + start;

                assertTrue(where + " did not converge: " + r.status, r.converged);
                assertTrue(where + " reports success but not a stopping rule", r.status.isSuccess());

                double beta = Digits.worstOf(r.parameters, set.certifiedBeta());
                double rss = Digits.of(r.sumOfSquares, set.residualSumOfSquares);
                assertTrue(where + ": parameters agree to " + beta + " digits", beta >= FLOORS[k][0]);
                assertTrue(where + ": residual sum of squares agrees to " + rss + " digits", rss >= FLOORS[k][2]);

                assertNotNull(where + ": no standard errors", r.standardErrors);
                double sd = Digits.worstOf(r.standardErrors, set.certifiedStandardDeviation());
                assertTrue(where + ": deviations agree to " + sd + " digits", sd >= FLOORS[k][1]);

                assertEquals(where + " degrees of freedom", set.observations - set.parameters, r.degreesOfFreedom);
            }
        }
    }

    /**
     * The residual sum of squares is always known better than the parameters
     * are, because the objective is flat where they sit: the value at the
     * minimum is determined long before the position of it is.
     */
    @Test
    public void theValueIsKnownBetterThanThePosition() {
        StRD.NonlinearSet[] sets = StRD.nonlinear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.NonlinearSet set = sets[k];
            LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(Models.of(set), set.start(2),
                    set.observations);
            double beta = Digits.worstOf(r.parameters, set.certifiedBeta());
            double rss = Digits.of(r.sumOfSquares, set.residualSumOfSquares);
            assertTrue(set.name + ": " + rss + " digits on the value against " + beta + " on the parameters",
                    rss >= beta);
        }
    }

    /** The far start costs more than the near one, and finds the same place. */
    @Test
    public void theFarStartCostsMoreAndArrives() {
        StRD.NonlinearSet[] sets = StRD.nonlinear();
        int dearer = 0;
        for (int k = 0; k < sets.length; ++k) {
            StRD.NonlinearSet set = sets[k];
            LevenbergMarquardt.Result far = new LevenbergMarquardt().solve(Models.of(set), set.start(1),
                    set.observations);
            LevenbergMarquardt.Result near = new LevenbergMarquardt().solve(Models.of(set), set.start(2),
                    set.observations);

            assertTrue(set.name + " from far", far.converged);
            assertTrue(set.name + " from near", near.converged);
            assertEquals(set.name + ": the two starts must reach the same minimum", far.sumOfSquares,
                    near.sumOfSquares, 1.0e-6 * Math.abs(near.sumOfSquares));
            if (far.functionEvaluations > near.functionEvaluations) {
                dearer++;
            }
        }
        assertTrue("the far start should be the dearer one on most sets: " + dearer + " of 5", dearer >= 4);
    }

    /**
     * The models must go through {@code StrictMath}, which is specified to
     * return the same value on every implementation, and not through
     * {@code Math}, which is allowed one unit in the last place and takes it:
     * at one of MGH10's sixteen model arguments the two disagree. That single
     * bit is worth two digits of the fitted parameters, because the objective
     * is flat where the minimum lies, and it made the fit differ between JDK 8
     * and JDK 25 until the models were changed.
     */
    @Test
    public void theModelsAreReproducibleAcrossImplementations() {
        StRD.NonlinearSet mgh10 = StRD.nonlinear()[4];
        assertEquals("MGH10", mgh10.name);
        double[] x = mgh10.x();
        double[] b = mgh10.certifiedBeta();

        int differing = 0;
        for (int i = 0; i < x.length; ++i) {
            double argument = b[1] / (x[i] + b[2]);
            if (Double.doubleToRawLongBits(Math.exp(argument)) != Double
                    .doubleToRawLongBits(StrictMath.exp(argument))) {
                differing++;
            }
            assertEquals("the model must use StrictMath at argument " + i,
                    Double.doubleToRawLongBits(b[0] * StrictMath.exp(argument)),
                    Double.doubleToRawLongBits(Models.value("MGH10", b, x[i])));
        }
        assertTrue("Math and StrictMath have to disagree somewhere here, or this guard proves nothing",
                differing > 0);
    }

    /**
     * With bounds too wide to bind, the bounded solver has the same problem to
     * solve. It arrives, and it arrives with fewer digits -- the projection is
     * a different algorithm, not the same one with a guard bolted on.
     */
    @Test
    public void boundsThatDoNotBindStillCostAccuracy() {
        StRD.NonlinearSet[] sets = StRD.nonlinear();
        int looser = 0;
        for (int k = 0; k < sets.length; ++k) {
            StRD.NonlinearSet set = sets[k];
            double[] lower = new double[set.parameters];
            double[] upper = new double[set.parameters];
            for (int j = 0; j < set.parameters; ++j) {
                lower[j] = -1.0e10;
                upper[j] = 1.0e10;
            }

            BoundedLevenbergMarquardt.Result bounded = new BoundedLevenbergMarquardt().solve(Models.of(set),
                    set.start(2), set.observations, lower, upper);
            LevenbergMarquardt.Result free = new LevenbergMarquardt().solve(Models.of(set), set.start(2),
                    set.observations);

            assertTrue(set.name + " bounded did not converge: " + bounded.status, bounded.converged);
            double boundedDigits = Digits.worstOf(bounded.parameters, set.certifiedBeta());
            double freeDigits = Digits.worstOf(free.parameters, set.certifiedBeta());
            assertTrue(set.name + ": bounded reached " + boundedDigits + " digits", boundedDigits >= 4.0);
            if (boundedDigits < freeDigits) {
                looser++;
            }
        }
        assertTrue("the bounded solver is the less accurate one on most sets: " + looser + " of 5", looser >= 3);
    }
}
