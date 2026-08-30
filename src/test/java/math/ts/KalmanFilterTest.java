package math.ts;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.linalg.SymmetricJacobiEigen;

/**
 * {@link KalmanFilter}, checked against oracles that share no line of code with
 * it.
 * <p>
 * The main one is {@link StackedJoint}: the whole series is one multivariate
 * normal whose mean and covariance the model determines without any recursion,
 * so the log likelihood is a single call to
 * {@link MultivariateNormal#logPdf(double[])} and the filtered state is what
 * conditioning that joint gives. The forward recursion and the closed form
 * agree to {@code 3.1e-15} relative on the likelihood, {@code 8.9e-15} on the
 * filtered mean and {@code 2.7e-14} on the covariances, measured over five
 * models including one with a singular {@code Q}.
 * <p>
 * Beside it, three closed forms: the local level filter <b>is</b> exponential
 * smoothing at its steady state and agrees with an EWMA to {@code 3.8e-14}
 * absolute; a state that never moves reproduces the batch normal-normal update
 * to {@code 2.8e-16}; and the steady state variance of an accurate sensor
 * matches its analytic value to {@code 2.1e-16} -- where the short covariance
 * update collapses it to exactly zero, which is what the Joseph form is here to
 * prevent.
 */
public final class KalmanFilterTest {

    /** Measured worst case 3.1e-15 relative, over the five models below. */
    private static final double LIKELIHOOD_TOL = 1.0e-12;
    /** Measured worst case 2.7e-14 relative, on the local linear trend. */
    private static final double STATE_TOL = 1.0e-11;

    private static double relative(double a, double b) {
        return Math.abs(a - b) / Math.max(1.0, Math.abs(b));
    }

    /** Every quantity of the forward pass against the stacked joint. */
    private static void agreesWithTheJoint(String what, StateSpaceModel model, DMatrix y) {
        int n = y.numRows();
        int d = model.stateDimension();
        StackedJoint joint = new StackedJoint(model, n);
        KalmanFilter.Result res = KalmanFilter.filter(model, y);

        assertEquals(what, 0.0, relative(res.logLikelihood, joint.logLikelihood(y)), LIKELIHOOD_TOL);
        assertEquals(what + ", the two static forms", res.logLikelihood,
                KalmanFilter.logLikelihood(model, y), 0.0);
        assertEquals(what, n, res.length);

        for (int t = 0; t < n; ++t) {
            double[] filteredMean = joint.conditionalMean(y, t, t + 1);
            DMatrix filteredCov = joint.conditionalCovariance(y, t, t + 1);
            double[] predictedMean = joint.conditionalMean(y, t, t);
            DMatrix predictedCov = joint.conditionalCovariance(y, t, t);
            for (int i = 0; i < d; ++i) {
                assertEquals(what + ", filtered mean at " + t, 0.0,
                        relative(res.filteredMean[t][i], filteredMean[i]), STATE_TOL);
                assertEquals(what + ", predicted mean at " + t, 0.0,
                        relative(res.predictedMean[t][i], predictedMean[i]), STATE_TOL);
                for (int j = 0; j < d; ++j) {
                    assertEquals(what + ", filtered covariance at " + t, 0.0,
                            relative(res.filteredCovariance[t].get(i, j), filteredCov.get(i, j)),
                            STATE_TOL);
                    assertEquals(what + ", predicted covariance at " + t, 0.0,
                            relative(res.predictedCovariance[t].get(i, j), predictedCov.get(i, j)),
                            STATE_TOL);
                }
            }
        }
    }

    @Test
    public void theForwardPassIsWhatTheStackedJointSays() {
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
    public void aMissingEntryIsStruckFromTheJointAndNothingElse() {
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
    public void observingNothingAtAllIsThePriorPropagated() {
        LinearGaussianModel model = Series.coupled();
        DMatrix blind = new DMatrix(5, 2);
        for (int t = 0; t < 5; ++t) {
            blind.set(t, 0, Double.NaN);
            blind.set(t, 1, Double.NaN);
        }
        KalmanFilter.Result res = KalmanFilter.filter(model, blind);
        // no evidence, so exactly no contribution
        assertEquals(0.0, res.logLikelihood, 0.0);

        double[] mu = new double[2];
        model.initialMean(mu);
        DMatrix f = model.transition();
        for (int t = 0; t < 5; ++t) {
            if (t > 0) {
                double a = f.get(0, 0) * mu[0] + f.get(0, 1) * mu[1];
                double b = f.get(1, 0) * mu[0] + f.get(1, 1) * mu[1];
                mu[0] = a;
                mu[1] = b;
            }
            // the deterministic propagation, and nothing has been subtracted
            // from it, so this is exact
            assertEquals(mu[0], res.filteredMean[t][0], 0.0);
            assertEquals(mu[1], res.filteredMean[t][1], 0.0);
            assertEquals(0, res.observedComponents[t]);
            assertEquals(0, res.innovation[t].length);
            assertNull(res.innovationCovariance[t]);
            assertTrue(Double.isNaN(res.squaredMahalanobis[t]));
        }
    }

    @Test
    public void aPermanentlyMissingComponentIsTheReducedModel() {
        DMatrix full = Series.draw(Series.coupled(), 12, 88L);
        DMatrix half = full.copy();
        DMatrix onlyFirst = new DMatrix(12, 1);
        for (int t = 0; t < 12; ++t) {
            half.set(t, 1, Double.NaN);
            onlyFirst.set(t, 0, full.get(t, 0));
        }
        LinearGaussianModel reduced = new LinearGaussianModel(Series.matrix(2, 2, 0.9, 0.2, -0.1, 0.8),
                Series.matrix(2, 2, 0.3, 0.1, 0.1, 0.2), Series.matrix(1, 2, 1.0, 0.5),
                Series.matrix(1, 1, 0.4), new double[] { 1.0, -2.0 },
                Series.matrix(2, 2, 2.0, 0.5, 0.5, 3.0));

        KalmanFilter.Result masked = KalmanFilter.filter(Series.coupled(), half);
        KalmanFilter.Result explicit = KalmanFilter.filter(reduced, onlyFirst);
        // striking the row out of H and R is the same arithmetic as never
        // having had it, so this holds bit for bit
        assertEquals(explicit.logLikelihood, masked.logLikelihood, 0.0);
        for (int t = 0; t < 12; ++t) {
            assertEquals(1, masked.observedComponents[t]);
            assertEquals(1, masked.innovation[t].length);
            for (int i = 0; i < 2; ++i) {
                assertEquals(explicit.filteredMean[t][i], masked.filteredMean[t][i], 0.0);
            }
        }
    }

    @Test
    public void theLocalLevelFilterIsExponentialSmoothing() {
        double q = 0.5;
        double r = 2.0;
        LinearGaussianModel model = Series.localLevel(q, r, 1.0e6);
        DMatrix y = Series.draw(Series.localLevel(q, r, 4.0), 400, 99L);
        KalmanFilter.Result res = KalmanFilter.filter(model, y);

        // the fixed point of the Riccati recursion P = P R / (P + R) + Q
        double pMinus = 0.5 * (q + Math.sqrt(q * q + 4.0 * q * r));
        double alpha = pMinus / (pMinus + r);
        assertEquals(0.0, relative(res.predictedCovariance[399].get(0, 0), pMinus), 1.0e-12);

        double ewma = res.filteredMean[0][0];
        for (int t = 1; t < 400; ++t) {
            ewma += alpha * (y.get(t, 0) - ewma);
            if (t >= 60) {
                // by then the gain has settled, and the two are the same rule
                // -- measured worst absolute gap 3.8e-14
                assertEquals("at " + t, ewma, res.filteredMean[t][0], 1.0e-10);
            }
        }
    }

    @Test
    public void theSteadyStateVarianceIsTheClosedFormOne() {
        double q = 1.0e-10;
        double r = 1.0e-10;
        LinearGaussianModel model = Series.localLevel(q, r, 1.0e8);
        KalmanFilter.Result res = KalmanFilter.filter(model, Series.draw(model, 2000, 111L));
        double pMinus = 0.5 * (q + Math.sqrt(q * q + 4.0 * q * r));
        double pFiltered = pMinus * r / (pMinus + r);
        for (int t = 200; t < 2000; ++t) {
            // 6.18e-11, which the short update (I - K H) P returns as exactly
            // zero at this scale. Measured worst relative gap 2.1e-16
            assertEquals("at " + t, 0.0,
                    Math.abs(res.filteredCovariance[t].get(0, 0) - pFiltered) / pFiltered, 1.0e-12);
        }
    }

    @Test
    public void aStaticStateIsTheBatchNormalNormalUpdate() {
        LinearGaussianModel model = Series.staticState();
        int n = 25;
        DMatrix y = Series.draw(model, n, 123L);
        KalmanFilter.Result res = KalmanFilter.filter(model, y);

        // F = I and Q = 0, so the state never moves and the posterior after n
        // observations is the ordinary conjugate update, written out here
        DMatrix h = model.observation();
        DMatrix rInv = model.observationNoise().inverse();
        DMatrix p0Inv = model.initialCovariance().inverse();
        DMatrix hrInv = h.transpose().mul(rInv);
        DMatrix precision = p0Inv.add(hrInv.mul(h).scale(n));
        double[] m0 = new double[2];
        model.initialMean(m0);
        double[] sumY = new double[2];
        for (int t = 0; t < n; ++t) {
            sumY[0] += y.get(t, 0);
            sumY[1] += y.get(t, 1);
        }
        double[] rhs = new double[2];
        for (int i = 0; i < 2; ++i) {
            rhs[i] = p0Inv.get(i, 0) * m0[0] + p0Inv.get(i, 1) * m0[1] + hrInv.get(i, 0) * sumY[0]
                    + hrInv.get(i, 1) * sumY[1];
        }
        double[] wantMean = new double[2];
        CholeskyDecomp.solve(CholeskyDecomp.cholesky(precision), rhs, wantMean);
        DMatrix wantCov = precision.inverse();

        for (int i = 0; i < 2; ++i) {
            // measured 2.8e-16 on the mean, 6.9e-18 on the covariance
            assertEquals(0.0, relative(res.filteredMean[n - 1][i], wantMean[i]), 1.0e-12);
            for (int j = 0; j < 2; ++j) {
                assertEquals(0.0,
                        relative(res.filteredCovariance[n - 1].get(i, j), wantCov.get(i, j)), 1.0e-12);
            }
        }
    }

    @Test
    public void theCovarianceStaysSymmetricAndFactorable() {
        // an accurate sensor, a vague prior and eighteen orders of magnitude
        // between the two: the short update leaves two of these 3000
        // covariances unfactorable, the Joseph form leaves none
        LinearGaussianModel hard = new LinearGaussianModel(Series.matrix(2, 2, 1.0, 1.0, 0.0, 1.0),
                Series.matrix(2, 2, 1.0e-16, 0.0, 0.0, 1.0e-18), Series.matrix(1, 2, 1.0, 0.0),
                Series.matrix(1, 1, 1.0e-12), new double[] { 0.0, 0.001 },
                Series.matrix(2, 2, 1.0e10, 0.0, 0.0, 1.0e10));
        KalmanFilter.Result res = KalmanFilter.filter(hard, Series.draw(hard, 3000, 222L));
        SymmetricJacobiEigen eigen = new SymmetricJacobiEigen();
        for (int t = 0; t < res.length; ++t) {
            DMatrix p = res.filteredCovariance[t];
            assertEquals("at " + t, p.get(0, 1), p.get(1, 0), 0.0);
            SymmetricJacobiEigen.Result e = eigen.decompose(p.getArrayUnsafe().clone(), 2);
            assertTrue("a negative eigenvalue at " + t, e.lambda[0] >= 0.0 && e.lambda[1] >= 0.0);
            CholeskyDecomp.cholesky(p);
        }
    }

    @Test
    public void theAnomalyScoreSeparatesASpike() {
        LinearGaussianModel model = Series.localLevel(0.1, 1.0, 4.0);
        DMatrix clean = Series.draw(model, 200, 333L);
        DMatrix spiked = clean.copy();
        spiked.set(150, 0, clean.get(150, 0) + 25.0);

        KalmanFilter.Result cleanRun = KalmanFilter.filter(model, clean);
        double worstClean = 0.0;
        for (int t = 0; t < 200; ++t) {
            worstClean = Math.max(worstClean, cleanRun.squaredMahalanobis[t]);
        }
        // measured 8.0 over 200 draws, which is where a chi-squared on one
        // degree of freedom belongs
        assertTrue("clean data should not look anomalous : " + worstClean, worstClean < 20.0);

        KalmanFilter.Result spikedRun = KalmanFilter.filter(model, spiked);
        // measured 485, sixty times the largest clean score
        assertTrue("the spike should stand out : " + spikedRun.squaredMahalanobis[150],
                spikedRun.squaredMahalanobis[150] > 20.0 * worstClean);

        // and it is the square of the standardized innovation, by construction
        double e = spikedRun.innovation[150][0];
        double s = spikedRun.innovationCovariance[150].get(0, 0);
        assertEquals(e * e / s, spikedRun.squaredMahalanobis[150], 1.0e-12);
    }

    @Test
    public void theSteppingObjectAndTheBatchAgreeBitForBit() {
        LinearGaussianModel model = Series.coupled();
        DMatrix y = Series.draw(model, 40, 444L);
        y.set(7, 0, Double.NaN);
        y.set(19, 0, Double.NaN);
        y.set(19, 1, Double.NaN);

        KalmanFilter.Result batch = KalmanFilter.filter(model, y);
        KalmanFilter kf = new KalmanFilter(model);
        double[] row = new double[2];
        double[] mean = new double[2];
        for (int t = 0; t < 40; ++t) {
            if (t > 0) {
                kf.predict();
            }
            row[0] = y.get(t, 0);
            row[1] = y.get(t, 1);
            kf.update(row);
            kf.mean(mean);
            for (int i = 0; i < 2; ++i) {
                assertEquals("at " + t, batch.filteredMean[t][i], mean[i], 0.0);
                for (int j = 0; j < 2; ++j) {
                    assertEquals("at " + t, batch.filteredCovariance[t].get(i, j),
                            kf.covariance().get(i, j), 0.0);
                }
            }
            assertEquals(batch.observedComponents[t], kf.observedComponents());
        }
        assertEquals(batch.logLikelihood, kf.logLikelihood(), 0.0);
        assertEquals(40, kf.steps());
    }

    @Test
    public void predictTwiceIsATwoStepAheadForecast() {
        LinearGaussianModel model = Series.trend();
        KalmanFilter kf = new KalmanFilter(model);
        DMatrix f = model.transition();
        DMatrix q = model.processNoise();

        double[] want = new double[2];
        model.initialMean(want);
        DMatrix wantCov = model.initialCovariance();
        for (int step = 1; step <= 3; ++step) {
            kf.predict();
            double a = f.get(0, 0) * want[0] + f.get(0, 1) * want[1];
            double b = f.get(1, 0) * want[0] + f.get(1, 1) * want[1];
            want[0] = a;
            want[1] = b;
            wantCov = f.mul(wantCov).mulBTrans(f).add(q);
            double[] got = new double[2];
            kf.mean(got);
            assertEquals(want[0], got[0], 1.0e-12);
            assertEquals(want[1], got[1], 1.0e-12);
            // the interval widens with every step and never narrows
            assertEquals(wantCov.get(0, 0), kf.covariance().get(0, 0), 1.0e-12);
        }
        // no observation was offered, so nothing was learned
        assertEquals(0.0, kf.logLikelihood(), 0.0);
        assertEquals(0, kf.steps());
    }

    @Test
    public void resetReturnsToThePrior() {
        LinearGaussianModel model = Series.trend();
        KalmanFilter kf = new KalmanFilter(model);
        DMatrix y = Series.draw(model, 20, 555L);
        double[] row = new double[1];
        for (int t = 0; t < 20; ++t) {
            if (t > 0) {
                kf.predict();
            }
            row[0] = y.get(t, 0);
            kf.update(row);
        }
        assertTrue(kf.logLikelihood() != 0.0);
        kf.reset();

        double[] prior = new double[2];
        model.initialMean(prior);
        double[] got = new double[2];
        kf.mean(got);
        assertEquals(prior[0], got[0], 0.0);
        assertEquals(prior[1], got[1], 0.0);
        assertEquals(model.initialCovariance().get(0, 0), kf.covariance().get(0, 0), 0.0);
        assertEquals(0.0, kf.logLikelihood(), 0.0);
        assertEquals(0, kf.steps());
        assertEquals(0, kf.observedComponents());
        assertEquals(0, kf.innovation().length);
        assertNull(kf.innovationCovariance());
        assertTrue(Double.isNaN(kf.squaredMahalanobisDistance()));
    }

    @Test
    public void correctingTwiceAgainstTheSamePredictionIsRefused() {
        KalmanFilter kf = new KalmanFilter(Series.localLevel(0.5, 2.0, 10.0));
        kf.update(new double[] { 1.0 });
        try {
            kf.update(new double[] { 1.0 });
            fail("the same evidence must not be counted twice");
        } catch (IllegalStateException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("predict"));
        }
        // and after a prediction it is allowed again
        kf.predict();
        kf.update(new double[] { 1.0 });
        assertEquals(2, kf.steps());
    }

    @Test
    public void theStatesComeBackAsTheDistributionsTheyAre() {
        LinearGaussianModel model = Series.coupled();
        DMatrix y = Series.draw(model, 6, 666L);
        KalmanFilter.Result res = KalmanFilter.filter(model, y);
        MultivariateNormal filtered = res.filteredState(3);
        assertEquals(2, filtered.dimension());
        double[] mean = new double[2];
        filtered.mean(mean);
        assertEquals(res.filteredMean[3][0], mean[0], 0.0);
        assertEquals(res.filteredCovariance[3].get(1, 1), filtered.covariance().get(1, 1), 0.0);

        MultivariateNormal predicted = res.predictedState(0);
        predicted.mean(mean);
        double[] prior = new double[2];
        model.initialMean(prior);
        // no transition runs before the first observation, so the predicted
        // state at index zero is the prior of the model itself
        assertEquals(prior[0], mean[0], 0.0);
        assertEquals(prior[1], mean[1], 0.0);

        KalmanFilter kf = new KalmanFilter(model);
        assertNotNull(kf.state());
        assertEquals(model.initialCovariance().get(0, 0), kf.state().covariance().get(0, 0), 0.0);
    }

    @Test
    public void theRefusals() {
        LinearGaussianModel model = Series.coupled();
        DMatrix y = Series.draw(model, 4, 777L);
        try {
            KalmanFilter.filter(null, y);
            fail("null model");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            KalmanFilter.filter(model, null);
            fail("null observations");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            KalmanFilter.logLikelihood(model, new DMatrix(4, 3));
            fail("the wrong number of columns");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("columns"));
        }
        try {
            new KalmanFilter(null);
            fail("null model");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        KalmanFilter kf = new KalmanFilter(model);
        try {
            kf.update(new double[] { 1.0 });
            fail("the wrong length");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("length"));
        }
        try {
            kf.update(new double[] { 1.0, Double.POSITIVE_INFINITY });
            fail("an infinite observation is neither a value nor a way of saying there is none");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("infinite"));
        }
        try {
            kf.mean(new double[3]);
            fail("the wrong length");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("length"));
        }
    }

    @Test
    public void aRepeatedScheduleIsTheConstantModelBitForBit() {
        // the cheapest check there is on the indexing, and the sharpest: an
        // off-by-one in F would still be a perfectly valid model, just not this
        // one, so only exact agreement rules it out
        LinearGaussianModel[] constants = { Series.localLevel(0.5, 2.0, 10.0), Series.trend(),
                Series.coupled(), Series.seasonal(), Series.staticState() };
        int n = 9;
        long seed = 11L;
        for (LinearGaussianModel constant : constants) {
            DMatrix y = Series.draw(constant, n, seed++);
            KalmanFilter.Result direct = KalmanFilter.filter(constant, y);
            KalmanFilter.Result scheduled = KalmanFilter.filter(Series.repeated(constant, n), y);
            assertEquals(direct.logLikelihood, scheduled.logLikelihood, 0.0);
            for (int t = 0; t < n; ++t) {
                for (int i = 0; i < constant.stateDimension(); ++i) {
                    assertEquals("at " + t, direct.filteredMean[t][i], scheduled.filteredMean[t][i], 0.0);
                    for (int j = 0; j < constant.stateDimension(); ++j) {
                        assertEquals("at " + t, direct.filteredCovariance[t].get(i, j),
                                scheduled.filteredCovariance[t].get(i, j), 0.0);
                    }
                }
            }
        }
    }

    @Test
    public void theFirstTransitionAndProcessNoiseAreNeverRead() {
        // F(t) leads into time t and nothing leads into the first time point,
        // so F[0] and Q[0] are there to keep all four schedules one length
        int n = 8;
        LinearGaussianModel trend = Series.trend();
        DMatrix y = Series.draw(trend, n, 33L);
        DMatrix[] f = TimeVaryingModel.repeat(trend.transition(), n);
        DMatrix[] q = TimeVaryingModel.repeat(trend.processNoise(), n);
        f[0] = Series.matrix(2, 2, 1.0e6, -3.0e5, 7.0, 42.0);
        q[0] = Series.matrix(2, 2, 9.0e9, 0.0, 0.0, 9.0e9);
        TimeVaryingModel poisoned = new TimeVaryingModel(f, q,
                TimeVaryingModel.repeat(trend.observation(), n),
                TimeVaryingModel.repeat(trend.observationNoise(), n), Series.meanOf(trend),
                trend.initialCovariance());

        KalmanFilter.Result sane = KalmanFilter.filter(Series.repeated(trend, n), y);
        KalmanFilter.Result loud = KalmanFilter.filter(poisoned, y);
        assertEquals(sane.logLikelihood, loud.logLikelihood, 0.0);
        for (int t = 0; t < n; ++t) {
            for (int i = 0; i < 2; ++i) {
                assertEquals("at " + t, sane.filteredMean[t][i], loud.filteredMean[t][i], 0.0);
            }
        }
    }

    @Test
    public void theVaryingModelsAreWhatTheStackedJointSays() {
        // unequal time steps: the state drifts in continuous time, so a step
        // spanning four units carries four times the process noise of one
        TimeVaryingModel irregular = Series.irregular(0.4, 1.0, 5.0, Series.gaps(9, 44L));
        agreesWithTheJoint("irregular sampling", irregular, Series.draw(irregular, 9, 55L));

        TimeVaryingModel sensor = Series.switchingSensor(8);
        agreesWithTheJoint("a sensor that changes", sensor, Series.draw(sensor, 8, 66L));

        TimeVaryingModel rotating = Series.rotating(7);
        agreesWithTheJoint("all four varying", rotating, Series.draw(rotating, 7, 77L));
    }

    @Test
    public void gapsInAVaryingModelChangeNothingAboutThat() {
        TimeVaryingModel sensor = Series.switchingSensor(8);
        DMatrix holed = Series.draw(sensor, 8, 88L);
        holed.set(3, 0, Double.NaN);
        holed.set(6, 0, Double.NaN);
        agreesWithTheJoint("a changing sensor with two gaps", sensor, holed);

        TimeVaryingModel rotating = Series.rotating(7);
        DMatrix partly = Series.draw(rotating, 7, 99L);
        partly.set(2, 0, Double.NaN);
        partly.set(2, 1, Double.NaN);
        partly.set(5, 1, Double.NaN);
        agreesWithTheJoint("all four varying, one gap and one half", rotating, partly);
    }

    @Test
    public void theOracleTellsTheScheduleFromAShiftOfIt() {
        // what makes the tests above worth running: the same transitions used
        // one step early are still a valid model, and the likelihood has to
        // notice. Measured -17.795 against -19.719, so nearly two nats apart
        TimeVaryingModel rotating = Series.rotating(7);
        DMatrix[] shifted = new DMatrix[7];
        DMatrix[] q = new DMatrix[7];
        DMatrix[] h = new DMatrix[7];
        DMatrix[] r = new DMatrix[7];
        for (int t = 0; t < 7; ++t) {
            shifted[t] = rotating.transition(Math.min(6, t + 1));
            q[t] = rotating.processNoise(t);
            h[t] = rotating.observation(t);
            r[t] = rotating.observationNoise(t);
        }
        TimeVaryingModel wrong = new TimeVaryingModel(shifted, q, h, r, Series.meanOf(rotating),
                rotating.initialCovariance());
        DMatrix y = Series.draw(rotating, 7, 77L);
        assertTrue(Math.abs(KalmanFilter.logLikelihood(rotating, y)
                - KalmanFilter.logLikelihood(wrong, y)) > 1.0);
    }

    @Test
    public void theSeriesMayNotOutrunTheModel() {
        TimeVaryingModel sensor = Series.switchingSensor(8);
        try {
            KalmanFilter.filter(sensor, new DMatrix(20, 1));
            fail("twenty time points through an eight step schedule");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("defined for 8"));
        }

        KalmanFilter kf = new KalmanFilter(sensor);
        double[] row = new double[1];
        for (int t = 0; t < 8; ++t) {
            if (t > 0) {
                kf.predict();
            }
            kf.update(row);
        }
        try {
            kf.predict();
            fail("a forecast the schedule has no matrices for");
        } catch (IllegalStateException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("step 8"));
        }
        // a constant model has no such limit, which is what a forecast needs
        KalmanFilter unbounded = new KalmanFilter(Series.trend());
        for (int t = 0; t < 100; ++t) {
            unbounded.predict();
        }
        assertTrue(unbounded.covariance().get(0, 0) > Series.trend().initialCovariance().get(0, 0));
    }

    @Test
    public void anImpossibleInnovationCovarianceNamesTheStep() {
        // no observation noise and a state that is already known exactly: S is
        // then the zero matrix and the step cannot be taken
        LinearGaussianModel model = new LinearGaussianModel(Series.matrix(1, 1, 1.0),
                new DMatrix(1, 1), Series.matrix(1, 1, 0.0), new DMatrix(1, 1), new double[] { 0.0 },
                Series.matrix(1, 1, 1.0));
        try {
            KalmanFilter.filter(model, Series.matrix(3, 1, 1.0, 2.0, 3.0));
            fail("a singular innovation covariance must be reported");
        } catch (ArithmeticException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("step 0"));
        }
    }
}
