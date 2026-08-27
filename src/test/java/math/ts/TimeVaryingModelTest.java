package math.ts;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.DMatrix;

/**
 * {@link TimeVaryingModel}, the schedule-backed {@link StateSpaceModel}.
 * <p>
 * The arithmetic of a varying model is checked where it belongs, against the
 * oracle in {@link KalmanFilterTest} and {@link RtsSmootherTest}. What is worth
 * testing here is what the constructor promises: that the four schedules have
 * to line up, that the dimensions may not drift with {@code t}, and that a
 * repeated instance is stored once rather than {@code n} times.
 */
public final class TimeVaryingModelTest {

    private static DMatrix matrix(int rows, int cols, double... values) {
        return Series.matrix(rows, cols, values);
    }

    private static DMatrix[] f(int n) {
        return TimeVaryingModel.repeat(matrix(2, 2, 1.0, 1.0, 0.0, 1.0), n);
    }

    private static DMatrix[] q(int n) {
        return TimeVaryingModel.repeat(matrix(2, 2, 0.04, 0.0, 0.0, 0.01), n);
    }

    private static DMatrix[] h(int n) {
        return TimeVaryingModel.repeat(matrix(1, 2, 1.0, 0.0), n);
    }

    private static DMatrix[] r(int n) {
        return TimeVaryingModel.repeat(matrix(1, 1, 0.25), n);
    }

    private static double[] m0() {
        return new double[] { 3.0, 0.5 };
    }

    private static DMatrix p0() {
        return matrix(2, 2, 10.0, 0.0, 0.0, 4.0);
    }

    @Test
    public void theScheduleIsWhatTheAccessorsReport() {
        int n = 6;
        DMatrix[] rs = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            rs[t] = matrix(1, 1, 0.1 * (t + 1));
        }
        TimeVaryingModel model = new TimeVaryingModel(f(n), q(n), h(n), rs, m0(), p0());
        assertEquals(2, model.stateDimension());
        assertEquals(1, model.observationDimension());
        assertEquals(n, model.length());
        for (int t = 0; t < n; ++t) {
            assertEquals(0.1 * (t + 1), model.observationNoise(t).get(0, 0), 1.0e-15);
        }
    }

    @Test
    public void aRepeatedInstanceIsStoredOnce() {
        // what keeps repeat(matrix, n) from costing n matrices, and what makes
        // "only H varies" pay for the varying part alone
        int n = 500;
        TimeVaryingModel model = new TimeVaryingModel(f(n), q(n), h(n), r(n), m0(), p0());
        for (int t = 1; t < n; ++t) {
            assertSame("at " + t, model.transition(0), model.transition(t));
            assertSame("at " + t, model.processNoise(0), model.processNoise(t));
            assertSame("at " + t, model.observation(0), model.observation(t));
            assertSame("at " + t, model.observationNoise(0), model.observationNoise(t));
        }
        // and it is still a copy of what went in, not the caller's matrix
        DMatrix[] given = f(n);
        TimeVaryingModel other = new TimeVaryingModel(given, q(n), h(n), r(n), m0(), p0());
        assertNotSame(given[0], other.transition(0));
        given[0].set(0, 1, 99.0);
        assertEquals(1.0, other.transition(3).get(0, 1), 0.0);
    }

    @Test
    public void entriesThatDifferAreStoredSeparately() {
        int n = 4;
        DMatrix[] fs = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            fs[t] = matrix(2, 2, 1.0, t, 0.0, 1.0);
        }
        TimeVaryingModel model = new TimeVaryingModel(fs, q(n), h(n), r(n), m0(), p0());
        for (int t = 0; t < n; ++t) {
            assertEquals(t, model.transition(t).get(0, 1), 0.0);
        }
    }

    @Test
    public void repeatRefusesWhatItCannotRepeat() {
        DMatrix a = matrix(1, 1, 1.0);
        assertEquals(3, TimeVaryingModel.repeat(a, 3).length);
        assertSame(a, TimeVaryingModel.repeat(a, 3)[2]);
        try {
            TimeVaryingModel.repeat(null, 3);
            fail("null");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            TimeVaryingModel.repeat(a, 0);
            fail("a schedule of no steps");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("positive"));
        }
    }

    @Test
    public void theFourSchedulesHaveToBeTheSameLength() {
        try {
            new TimeVaryingModel(f(5), q(4), h(5), r(5), m0(), p0());
            fail("Q is one short");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("not all equal"));
        }
        try {
            new TimeVaryingModel(new DMatrix[0], new DMatrix[0], new DMatrix[0], new DMatrix[0],
                    m0(), p0());
            fail("a series of no time points");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("at least one time point"));
        }
    }

    @Test
    public void noEntryMayBeNull() {
        DMatrix[] holed = f(4);
        holed[2] = null;
        try {
            new TimeVaryingModel(holed, q(4), h(4), r(4), m0(), p0());
            fail("a null entry");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("F[2]"));
        }
    }

    @Test
    public void theDimensionsMayNotDriftWithTime() {
        // a sensor that appears halfway through is a NaN in the data, not a
        // taller H -- the class refuses the second so nobody reaches for it
        DMatrix[] hs = h(4);
        hs[2] = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        try {
            new TimeVaryingModel(f(4), q(4), hs, r(4), m0(), p0());
            fail("H grew a row");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("H[2] has 2 rows, not 1"));
        }

        DMatrix[] fs = f(4);
        fs[3] = matrix(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1);
        try {
            new TimeVaryingModel(fs, q(4), h(4), r(4), m0(), p0());
            fail("F changed order");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("F[3] is of order 3, not 2"));
        }
    }

    @Test
    public void theCovariancePolicyIsTheSameOneAsForAConstantModel() {
        int n = 4;
        // a singular Q at one step is legitimate, as it is throughout
        DMatrix[] qs = q(n);
        qs[2] = matrix(2, 2, 0.0, 0.0, 0.0, 0.01);
        TimeVaryingModel model = new TimeVaryingModel(f(n), qs, h(n), r(n), m0(), p0());
        assertEquals(0.0, model.processNoise(2).get(0, 0), 0.0);

        // and an asymmetric one is refused with the index that named it
        DMatrix[] bad = q(n);
        bad[1] = matrix(2, 2, 1.0, 0.5, 0.4, 1.0);
        try {
            new TimeVaryingModel(f(n), bad, h(n), r(n), m0(), p0());
            fail("an asymmetric Q");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("Q[1] is not symmetric"));
        }

        // a singular P0 is refused here as it is there
        try {
            new TimeVaryingModel(f(n), q(n), h(n), r(n), m0(), matrix(2, 2, 1.0, 0.0, 0.0, 0.0));
            fail("a singular P0");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("P0"));
        }
    }

    @Test
    public void theRefusals() {
        try {
            new TimeVaryingModel(null, q(3), h(3), r(3), m0(), p0());
            fail("null schedule");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            new TimeVaryingModel(f(3), q(3), h(3), r(3), null, p0());
            fail("null m0");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("null"));
        }
        try {
            new TimeVaryingModel(f(3), q(3), h(3), r(3), new double[] { 1.0 }, p0());
            fail("m0 of the wrong length");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("length"));
        }
        DMatrix[] fs = f(3);
        fs[1] = matrix(2, 2, 1.0, Double.NaN, 0.0, 1.0);
        try {
            new TimeVaryingModel(fs, q(3), h(3), r(3), m0(), p0());
            fail("an entry that is not finite");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("finite"));
        }
    }

    @Test
    public void theInitialStateIsTheNormalItSays() {
        TimeVaryingModel model = new TimeVaryingModel(f(3), q(3), h(3), r(3), m0(), p0());
        assertEquals(2, model.initialState().dimension());
        assertEquals(10.0, model.initialState().covariance().get(0, 0), 0.0);
        double[] mean = new double[2];
        model.initialState().mean(mean);
        assertEquals(3.0, mean[0], 0.0);
        assertEquals(0.5, mean[1], 0.0);
    }
}
