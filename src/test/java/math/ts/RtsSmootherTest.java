package math.ts;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.MultivariateNormal;
import math.linalg.DMatrix;
import math.linalg.SymmetricJacobiEigen;

/**
 * {@link RtsSmoother}, against the same oracle the forward pass was checked
 * with.
 * <p>
 * {@link StackedJoint} conditions the whole joint on <b>all</b> the
 * observations rather than the first {@code t} of them, which is the definition
 * of what the backward pass computes, and it does so by ordinary Gaussian
 * conditioning rather than by any recursion. The two agree to {@code 9.4e-15}
 * relative on the smoothed mean, {@code 1.3e-14} on the covariance and
 * {@code 2.5e-14} on the lag-one covariance, over five models with and without
 * gaps.
 * <p>
 * The lag-one covariance is worth a separate word: it is the one quantity here
 * that is not symmetric, so a transposed factor would pass every other test in
 * this file. It is checked against the joint's own cross-covariance, in the
 * stated orientation.
 */
public final class RtsSmootherTest {

    /** Measured worst case 2.5e-14 relative, on the local linear trend. */
    private static final double SMOOTHED_TOL = 1.0e-11;

    private static double relative(double a, double b) {
        return Math.abs(a - b) / Math.max(1.0, Math.abs(b));
    }

    private static void agreesWithTheJoint(String what, StateSpaceModel model, DMatrix y) {
        int n = y.numRows();
        int d = model.stateDimension();
        StackedJoint joint = new StackedJoint(model, n);
        KalmanFilter.Result filtered = KalmanFilter.filter(model, y);
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
        assertEquals(what, n, smoothed.length);

        for (int t = 0; t < n; ++t) {
            double[] wantMean = joint.conditionalMean(y, t, n);
            DMatrix wantCov = joint.conditionalCovariance(y, t, n);
            for (int i = 0; i < d; ++i) {
                assertEquals(what + ", smoothed mean at " + t, 0.0,
                        relative(smoothed.mean[t][i], wantMean[i]), SMOOTHED_TOL);
                for (int j = 0; j < d; ++j) {
                    assertEquals(what + ", smoothed covariance at " + t, 0.0,
                            relative(smoothed.covariance[t].get(i, j), wantCov.get(i, j)),
                            SMOOTHED_TOL);
                }
            }
            if (t < n - 1) {
                DMatrix wantLag = joint.conditionalCrossCovariance(y, t, t + 1, n);
                DMatrix gotLag = smoothed.lagOneCovariance(t);
                for (int i = 0; i < d; ++i) {
                    for (int j = 0; j < d; ++j) {
                        assertEquals(what + ", lag-one covariance at " + t, 0.0,
                                relative(gotLag.get(i, j), wantLag.get(i, j)), SMOOTHED_TOL);
                    }
                }
            }
        }
    }

    @Test
    public void theBackwardPassIsTheJointConditionedOnEverything() {
        agreesWithTheJoint("local level", Series.localLevel(0.5, 2.0, 10.0),
                Series.draw(Series.localLevel(0.5, 2.0, 10.0), 8, 11L));
        agreesWithTheJoint("local linear trend", Series.trend(), Series.draw(Series.trend(), 8, 22L));
        agreesWithTheJoint("coupled", Series.coupled(), Series.draw(Series.coupled(), 7, 33L));
        agreesWithTheJoint("seasonal, singular Q", Series.seasonal(),
                Series.draw(Series.seasonal(), 9, 44L));
        agreesWithTheJoint("static state", Series.staticState(),
                Series.draw(Series.staticState(), 7, 55L));
    }

    @Test
    public void gapsChangeNothingAboutThat() {
        DMatrix holed = Series.draw(Series.coupled(), 7, 66L);
        holed.set(2, 0, Double.NaN);
        holed.set(2, 1, Double.NaN);
        holed.set(4, 1, Double.NaN);
        holed.set(5, 0, Double.NaN);
        agreesWithTheJoint("one whole gap and two half gaps", Series.coupled(), holed);

        DMatrix scalar = Series.draw(Series.seasonal(), 9, 77L);
        scalar.set(0, 0, Double.NaN);
        scalar.set(3, 0, Double.NaN);
        scalar.set(8, 0, Double.NaN);
        agreesWithTheJoint("a scalar series with three gaps", Series.seasonal(), scalar);
    }

    @Test
    public void theBackwardPassFollowsAVaryingModel() {
        // the backward gain is C[t] = P[t] F(t+1)' P[t+1|t]^-1, and F(t+1)
        // rather than F(t) is the whole content of this test: the wrong index
        // is still a plausible matrix, and only the joint notices
        TimeVaryingModel irregular = Series.irregular(0.4, 1.0, 5.0, Series.gaps(9, 44L));
        agreesWithTheJoint("irregular sampling", irregular, Series.draw(irregular, 9, 55L));

        TimeVaryingModel sensor = Series.switchingSensor(8);
        agreesWithTheJoint("a sensor that changes", sensor, Series.draw(sensor, 8, 66L));

        TimeVaryingModel rotating = Series.rotating(7);
        agreesWithTheJoint("all four varying", rotating, Series.draw(rotating, 7, 77L));

        DMatrix holed = Series.draw(sensor, 8, 88L);
        holed.set(3, 0, Double.NaN);
        holed.set(6, 0, Double.NaN);
        agreesWithTheJoint("a changing sensor with two gaps", sensor, holed);
    }

    @Test
    public void aRepeatedScheduleSmoothsIdentically() {
        int n = 9;
        LinearGaussianModel constant = Series.trend();
        DMatrix y = Series.draw(constant, n, 12L);
        RtsSmoother.Result direct = RtsSmoother.smooth(constant, y);
        RtsSmoother.Result scheduled = RtsSmoother.smooth(Series.repeated(constant, n), y);
        for (int t = 0; t < n; ++t) {
            for (int i = 0; i < 2; ++i) {
                assertEquals("at " + t, direct.mean[t][i], scheduled.mean[t][i], 0.0);
                for (int j = 0; j < 2; ++j) {
                    assertEquals("at " + t, direct.covariance[t].get(i, j),
                            scheduled.covariance[t].get(i, j), 0.0);
                }
            }
        }
    }

    @Test
    public void theFilteredHistoryMayNotOutrunTheModel() {
        TimeVaryingModel sensor = Series.switchingSensor(8);
        KalmanFilter.Result longer = KalmanFilter.filter(Series.trend(), Series.draw(Series.trend(), 20, 13L));
        try {
            RtsSmoother.smooth(sensor, longer);
            fail("twenty time points against an eight step schedule");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("defined for 8"));
        }
    }

    @Test
    public void theLastIndexIsTheFilteredStateBitForBit() {
        LinearGaussianModel model = Series.trend();
        DMatrix y = Series.draw(model, 30, 88L);
        KalmanFilter.Result filtered = KalmanFilter.filter(model, y);
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
        int last = smoothed.length - 1;
        // nothing comes after it, so the recursion starts by copying it and
        // there is no arithmetic to disagree about
        for (int i = 0; i < 2; ++i) {
            assertEquals(filtered.filteredMean[last][i], smoothed.mean[last][i], 0.0);
            for (int j = 0; j < 2; ++j) {
                assertEquals(filtered.filteredCovariance[last].get(i, j),
                        smoothed.covariance[last].get(i, j), 0.0);
            }
        }
        assertNull("there is no state after the last one", smoothed.gain[last]);
        try {
            smoothed.lagOneCovariance(last);
            fail("there is no pair at the last index");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("last index"));
        }
    }

    @Test
    public void smoothingNeverWidensTheInterval() {
        StateSpaceModel[] models = { Series.localLevel(0.5, 2.0, 10.0), Series.trend(),
                Series.coupled(), Series.seasonal(), Series.switchingSensor(40), Series.rotating(40) };
        long seed = 101L;
        SymmetricJacobiEigen eigen = new SymmetricJacobiEigen();
        for (StateSpaceModel model : models) {
            int d = model.stateDimension();
            DMatrix y = Series.draw(model, 40, seed++);
            KalmanFilter.Result filtered = KalmanFilter.filter(model, y);
            RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
            for (int t = 0; t < smoothed.length; ++t) {
                DMatrix difference = filtered.filteredCovariance[t].minus(smoothed.covariance[t]);
                double scale = 0.0;
                for (int i = 0; i < d; ++i) {
                    scale = Math.max(scale, Math.abs(filtered.filteredCovariance[t].get(i, i)));
                    // no coordinate can get less certain by being told more
                    assertTrue("variance grew at " + t, smoothed.covariance[t].get(i, i)
                            <= filtered.filteredCovariance[t].get(i, i) * (1.0 + 1.0e-12));
                }
                SymmetricJacobiEigen.Result e = eigen.decompose(difference.getArrayUnsafe().clone(), d);
                for (int i = 0; i < d; ++i) {
                    // measured worst relative eigenvalue -1.5e-16, which is
                    // round-off in the difference and not a real direction
                    assertTrue("P_filtered - P_smoothed is not positive semidefinite at " + t,
                            e.lambda[i] >= -1.0e-12 * scale);
                }
            }
        }
    }

    @Test
    public void observingNothingLeavesThePropagatedPrior() {
        LinearGaussianModel model = Series.coupled();
        DMatrix blind = new DMatrix(5, 2);
        for (int t = 0; t < 5; ++t) {
            blind.set(t, 0, Double.NaN);
            blind.set(t, 1, Double.NaN);
        }
        KalmanFilter.Result filtered = KalmanFilter.filter(model, blind);
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
        for (int t = 0; t < 5; ++t) {
            for (int i = 0; i < 2; ++i) {
                // there is nothing to revise the past with, so the backward
                // pass has to leave it exactly where the forward pass put it
                assertEquals(filtered.filteredMean[t][i], smoothed.mean[t][i], 0.0);
            }
        }
    }

    @Test
    public void aStaticStateGivesTheSameAnswerAtEveryTime() {
        LinearGaussianModel model = Series.staticState();
        DMatrix y = Series.draw(model, 20, 123L);
        KalmanFilter.Result filtered = KalmanFilter.filter(model, y);
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
        // the state never moves, so every time point has the same posterior,
        // and it is the one the forward pass ends on
        for (int t = 0; t < 20; ++t) {
            for (int i = 0; i < 2; ++i) {
                assertEquals("at " + t, filtered.filteredMean[19][i], smoothed.mean[t][i], 0.0);
                for (int j = 0; j < 2; ++j) {
                    assertEquals("at " + t, filtered.filteredCovariance[19].get(i, j),
                            smoothed.covariance[t].get(i, j), 1.0e-12);
                }
            }
        }
    }

    @Test
    public void theBackwardPassActuallyRemovesUncertainty() {
        LinearGaussianModel model = Series.trend();
        DMatrix y = Series.draw(model, 100, 222L);
        KalmanFilter.Result filtered = KalmanFilter.filter(model, y);
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, filtered);
        // measured: 0.244 down to 0.130 at the start of the series and 0.132
        // down to 0.060 in the middle, where the whole future is available
        assertTrue(smoothed.covariance[0].get(0, 0) < 0.6 * filtered.filteredCovariance[0].get(0, 0));
        assertTrue(smoothed.covariance[50].get(0, 0) < 0.6 * filtered.filteredCovariance[50].get(0, 0));
        // and nothing at all at the end, where there is no future to use
        assertEquals(filtered.filteredCovariance[99].get(0, 0), smoothed.covariance[99].get(0, 0), 0.0);
    }

    @Test
    public void theConvenienceFormIsTheTwoStepForm() {
        LinearGaussianModel model = Series.seasonal();
        DMatrix y = Series.draw(model, 25, 333L);
        y.set(9, 0, Double.NaN);
        RtsSmoother.Result direct = RtsSmoother.smooth(model, y);
        RtsSmoother.Result staged = RtsSmoother.smooth(model, KalmanFilter.filter(model, y));
        for (int t = 0; t < 25; ++t) {
            for (int i = 0; i < 3; ++i) {
                assertEquals(staged.mean[t][i], direct.mean[t][i], 0.0);
                for (int j = 0; j < 3; ++j) {
                    assertEquals(staged.covariance[t].get(i, j), direct.covariance[t].get(i, j), 0.0);
                }
            }
        }
    }

    @Test
    public void theSmoothedStateComesBackAsTheDistributionItIs() {
        LinearGaussianModel model = Series.coupled();
        RtsSmoother.Result smoothed = RtsSmoother.smooth(model, Series.draw(model, 10, 444L));
        MultivariateNormal state = smoothed.smoothedState(4);
        assertEquals(2, state.dimension());
        double[] mean = new double[2];
        state.mean(mean);
        assertEquals(smoothed.mean[4][0], mean[0], 0.0);
        assertEquals(smoothed.mean[4][1], mean[1], 0.0);
        assertEquals(smoothed.covariance[4].get(0, 1), state.covariance().get(0, 1), 0.0);
    }

    @Test
    public void anUninvertiblePredictedCovarianceNamesTheStep() {
        // a transition that annihilates the state and no process noise, so the
        // predicted covariance is the zero matrix from the second step on
        LinearGaussianModel model = new LinearGaussianModel(new DMatrix(2, 2), new DMatrix(2, 2),
                Series.matrix(1, 2, 1.0, 0.0), Series.matrix(1, 1, 1.0), new double[] { 0.0, 0.0 },
                Series.matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
        try {
            RtsSmoother.smooth(model, Series.matrix(4, 1, 1.0, 2.0, 3.0, 4.0));
            fail("the backward gain is undefined here");
        } catch (ArithmeticException e) {
            // the last one, not the first: this pass runs from the end
            assertTrue(e.getMessage(), e.getMessage().contains("step 3"));
        }
    }

    @Test
    public void theRefusals() {
        LinearGaussianModel model = Series.coupled();
        KalmanFilter.Result filtered = KalmanFilter.filter(model, Series.draw(model, 5, 555L));
        try {
            RtsSmoother.smooth(null, filtered);
            fail("null model");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            RtsSmoother.smooth(model, (KalmanFilter.Result) null);
            fail("null history");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            RtsSmoother.smooth(model, (DMatrix) null);
            fail("null observations");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            // a history of the wrong model
            RtsSmoother.smooth(Series.seasonal(), filtered);
            fail("three state coordinates against a history of two");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("coordinates"));
        }
    }
}
