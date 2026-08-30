package math.ts;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;

/**
 * {@link LinearGaussianModel}, the value type the two passes read.
 * <p>
 * What is worth a test here is not arithmetic but policy: which degeneracies
 * are legitimate and which are not. A singular {@code Q} is the normal case and
 * a singular {@code P0} is a broken model, and the class has to tell them
 * apart.
 */
public final class LinearGaussianModelTest {

    private static DMatrix matrix(int rows, int cols, double... values) {
        DMatrix a = new DMatrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                a.set(i, j, values[i * cols + j]);
            }
        }
        return a;
    }

    /** The local linear trend: a level that moves at a slope that drifts. */
    private static LinearGaussianModel trend() {
        DMatrix f = matrix(2, 2, 1.0, 1.0, 0.0, 1.0);
        DMatrix q = matrix(2, 2, 0.04, 0.0, 0.0, 0.01);
        DMatrix h = matrix(1, 2, 1.0, 0.0);
        DMatrix r = matrix(1, 1, 0.25);
        DMatrix p0 = matrix(2, 2, 10.0, 0.0, 0.0, 4.0);
        return new LinearGaussianModel(f, q, h, r, new double[] { 3.0, 0.5 }, p0);
    }

    @Test
    public void theAccessorsReportWhatWentIn() {
        LinearGaussianModel model = trend();
        assertEquals(2, model.stateDimension());
        assertEquals(1, model.observationDimension());
        assertEquals(1.0, model.transition().get(0, 1), 0.0);
        assertEquals(0.04, model.processNoise().get(0, 0), 0.0);
        assertEquals(1.0, model.observation().get(0, 0), 0.0);
        assertEquals(0.25, model.observationNoise().get(0, 0), 0.0);
        assertEquals(4.0, model.initialCovariance().get(1, 1), 0.0);
        double[] mean = new double[2];
        model.initialMean(mean);
        assertEquals(3.0, mean[0], 0.0);
        assertEquals(0.5, mean[1], 0.0);
    }

    @Test
    public void nothingTheCallerHoldsCanChangeTheModel() {
        DMatrix f = matrix(2, 2, 1.0, 1.0, 0.0, 1.0);
        DMatrix q = matrix(2, 2, 0.04, 0.0, 0.0, 0.01);
        DMatrix h = matrix(1, 2, 1.0, 0.0);
        DMatrix r = matrix(1, 1, 0.25);
        DMatrix p0 = matrix(2, 2, 10.0, 0.0, 0.0, 4.0);
        double[] m0 = new double[] { 3.0, 0.5 };
        LinearGaussianModel model = new LinearGaussianModel(f, q, h, r, m0, p0);

        // the arguments, afterwards
        f.set(0, 1, 99.0);
        q.set(0, 0, 99.0);
        h.set(0, 0, 99.0);
        r.set(0, 0, 99.0);
        p0.set(0, 0, 99.0);
        m0[0] = 99.0;

        assertEquals(1.0, model.transition().get(0, 1), 0.0);
        assertEquals(0.04, model.processNoise().get(0, 0), 0.0);
        assertEquals(1.0, model.observation().get(0, 0), 0.0);
        assertEquals(0.25, model.observationNoise().get(0, 0), 0.0);
        assertEquals(10.0, model.initialCovariance().get(0, 0), 0.0);
        double[] mean = new double[2];
        model.initialMean(mean);
        assertEquals(3.0, mean[0], 0.0);
    }

    @Test
    public void whatComesBackOutIsTheModelsOwn() {
        // the other half of the contract, and the one that changed when the
        // filter started reading the matrices at every step instead of once:
        // copying them out per step would allocate for nothing, so the model
        // hands out its own and asks not to be written to
        LinearGaussianModel model = trend();
        assertSame(model.transition(), model.transition());
        assertSame(model.transition(), model.transition(0));
        assertSame(model.transition(0), model.transition(7));
        assertSame(model.processNoise(), model.processNoise(3));
        assertSame(model.observation(), model.observation(3));
        assertSame(model.observationNoise(), model.observationNoise(3));
        assertSame(model.initialCovariance(), model.initialCovariance());
    }

    @Test
    public void aConstantModelDoesNotRunOut() {
        // the default of the interface: a bound check against it is one
        // comparison that never fails, which is what the filter relies on
        assertEquals(Integer.MAX_VALUE, trend().length());
    }

    @Test
    public void aSingularProcessNoiseIsTheNormalCase() {
        // the smooth trend model: the only shock is on the slope, so Q has
        // rank one, and the state space form of an ARMA model is worse still
        DMatrix q = matrix(2, 2, 0.0, 0.0, 0.0, 0.01);
        LinearGaussianModel model = new LinearGaussianModel(matrix(2, 2, 1.0, 1.0, 0.0, 1.0), q,
                matrix(1, 2, 1.0, 0.0), matrix(1, 1, 0.25), new double[] { 0.0, 0.0 },
                matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
        assertEquals(0.0, model.processNoise().get(0, 0), 0.0);
    }

    @Test
    public void aStateThatMovesWithoutNoiseIsAccepted() {
        DMatrix zero = new DMatrix(2, 2);
        LinearGaussianModel model = new LinearGaussianModel(matrix(2, 2, 1.0, 1.0, 0.0, 1.0), zero,
                matrix(1, 2, 1.0, 0.0), matrix(1, 1, 0.25), new double[] { 0.0, 0.0 },
                matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
        assertEquals(0.0, model.processNoise().get(1, 1), 0.0);
    }

    @Test
    public void aSingularInitialCovarianceIsRefused() {
        try {
            new LinearGaussianModel(matrix(2, 2, 1.0, 1.0, 0.0, 1.0), matrix(2, 2, 0.04, 0.0, 0.0, 0.01),
                    matrix(1, 2, 1.0, 0.0), matrix(1, 1, 0.25), new double[] { 0.0, 0.0 },
                    matrix(2, 2, 1.0, 0.0, 0.0, 0.0));
            fail("a singular P0 must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("P0"));
        }
    }

    @Test
    public void asymmetryWithinToleranceIsAveraged() {
        DMatrix q = matrix(2, 2, 1.0, 0.5, 0.5 + 1.0e-14, 1.0);
        LinearGaussianModel model = new LinearGaussianModel(matrix(2, 2, 1.0, 0.0, 0.0, 1.0), q,
                matrix(1, 2, 1.0, 0.0), matrix(1, 1, 1.0), new double[] { 0.0, 0.0 },
                matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
        DMatrix stored = model.processNoise();
        assertEquals(stored.get(0, 1), stored.get(1, 0), 0.0);
        assertEquals(0.5 + 0.5e-14, stored.get(0, 1), 1.0e-15);
    }

    @Test
    public void asymmetryBeyondToleranceIsRefused() {
        try {
            new LinearGaussianModel(matrix(2, 2, 1.0, 0.0, 0.0, 1.0), matrix(2, 2, 1.0, 0.5, 0.4, 1.0),
                    matrix(1, 2, 1.0, 0.0), matrix(1, 1, 1.0), new double[] { 0.0, 0.0 },
                    matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
            fail("an asymmetric Q must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("symmetric"));
        }
    }

    @Test
    public void theSymmetryToleranceIsRelativeAndHasToBe() {
        // a covariance at a scale of 1e6 whose two halves differ by 1e-5:
        // round-off for that scale, and four orders of magnitude beyond the
        // absolute 1e-10 that CholeskyDecomp itself applies
        DMatrix p0 = matrix(2, 2, 4.0e6, 1.0e6, 1.0e6 + 1.0e-5, 9.0e6);
        try {
            CholeskyDecomp.cholesky(p0);
            fail("the absolute check was expected to reject this");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("symmetric"));
        }
        LinearGaussianModel model = new LinearGaussianModel(matrix(2, 2, 1.0, 0.0, 0.0, 1.0),
                matrix(2, 2, 1.0, 0.0, 0.0, 1.0), matrix(1, 2, 1.0, 0.0), matrix(1, 1, 1.0),
                new double[] { 0.0, 0.0 }, p0);
        DMatrix stored = model.initialCovariance();
        assertEquals(stored.get(0, 1), stored.get(1, 0), 0.0);
    }

    @Test
    public void aNegativeVarianceIsRefused() {
        try {
            new LinearGaussianModel(matrix(2, 2, 1.0, 0.0, 0.0, 1.0), matrix(2, 2, -1.0, 0.0, 0.0, 1.0),
                    matrix(1, 2, 1.0, 0.0), matrix(1, 1, 1.0), new double[] { 0.0, 0.0 },
                    matrix(2, 2, 1.0, 0.0, 0.0, 1.0));
            fail("a negative variance must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("negative variance"));
        }
    }

    @Test
    public void theShapesHaveToAgree() {
        DMatrix f = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        DMatrix q = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        DMatrix h = matrix(1, 2, 1.0, 0.0);
        DMatrix r = matrix(1, 1, 1.0);
        DMatrix p0 = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        double[] m0 = new double[] { 0.0, 0.0 };

        refused(matrix(2, 3, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0), q, h, r, m0, p0, "square");
        refused(f, matrix(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1), h, r, m0, p0, "order");
        refused(f, q, matrix(1, 3, 1.0, 0.0, 0.0), r, m0, p0, "columns");
        refused(f, q, h, matrix(2, 2, 1.0, 0.0, 0.0, 1.0), m0, p0, "order");
        refused(f, q, h, r, new double[] { 0.0 }, p0, "length");
        refused(f, q, h, r, m0, matrix(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1), "order");
    }

    @Test
    public void nothingMayBeNullAndNothingMayBeInfinite() {
        DMatrix f = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        DMatrix q = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        DMatrix h = matrix(1, 2, 1.0, 0.0);
        DMatrix r = matrix(1, 1, 1.0);
        DMatrix p0 = matrix(2, 2, 1.0, 0.0, 0.0, 1.0);
        double[] m0 = new double[] { 0.0, 0.0 };

        refused(null, q, h, r, m0, p0, "null");
        refused(f, null, h, r, m0, p0, "null");
        refused(f, q, null, r, m0, p0, "null");
        refused(f, q, h, null, m0, p0, "null");
        refused(f, q, h, r, null, p0, "null");
        refused(f, q, h, r, m0, null, "null");

        refused(matrix(2, 2, Double.NaN, 0.0, 0.0, 1.0), q, h, r, m0, p0, "finite");
        refused(f, matrix(2, 2, Double.POSITIVE_INFINITY, 0.0, 0.0, 1.0), h, r, m0, p0, "finite");
        refused(f, q, matrix(1, 2, Double.NaN, 0.0), r, m0, p0, "finite");
        refused(f, q, h, matrix(1, 1, Double.NaN), m0, p0, "finite");
        refused(f, q, h, r, new double[] { Double.NaN, 0.0 }, p0, "finite");
        refused(f, q, h, r, m0, matrix(2, 2, 1.0, Double.NaN, Double.NaN, 1.0), "finite");
    }

    @Test
    public void theInitialStateIsTheNormalItSays() {
        LinearGaussianModel model = trend();
        MultivariateNormal prior = model.initialState();
        assertEquals(2, prior.dimension());
        double[] mean = new double[2];
        prior.mean(mean);
        assertEquals(3.0, mean[0], 0.0);
        assertEquals(0.5, mean[1], 0.0);
        assertEquals(10.0, prior.covariance().get(0, 0), 0.0);
        assertEquals(4.0, prior.covariance().get(1, 1), 0.0);
    }

    @Test
    public void theMeanIsWrittenIntoTheCallersArray() {
        LinearGaussianModel model = trend();
        try {
            model.initialMean(null);
            fail("null must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("null"));
        }
        try {
            model.initialMean(new double[3]);
            fail("the wrong length must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("length"));
        }
    }

    private static void refused(DMatrix f, DMatrix q, DMatrix h, DMatrix r, double[] m0, DMatrix p0,
            String word) {
        try {
            new LinearGaussianModel(f, q, h, r, m0, p0);
            fail("expected a refusal mentioning \"" + word + "\"");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains(word));
        }
    }
}
