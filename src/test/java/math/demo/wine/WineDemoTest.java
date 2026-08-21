package math.demo.wine;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import org.junit.Test;

import math.optim.LimitedMemoryBFGS;
import math.optim.Termination;

/**
 * The bounds in here were measured before they were asserted, over several
 * fold seeds where a seed could move them. Nothing is a golden value that a
 * rerun produced once.
 */
public class WineDemoTest {

    private static final double[] FOLD_SEEDS = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };

    @Test
    public void dataIsIntact() {
        assertEquals(178, Datasets.size());
        int[] counts = Datasets.classCounts();
        assertEquals(59, counts[0]);
        assertEquals(71, counts[1]);
        assertEquals(48, counts[2]);
        assertEquals(159975.296, Datasets.checksum(), 1.0e-3);

        double[][] x = Datasets.features();
        assertEquals(178, x.length);
        for (int i = 0; i < x.length; ++i) {
            assertEquals(WineDemo.P, x[i].length);
            for (int j = 0; j < WineDemo.P; ++j) {
                assertTrue("cell " + i + "," + j + " is not finite", !Double.isNaN(x[i][j])
                        && !Double.isInfinite(x[i][j]));
                assertTrue("cell " + i + "," + j + " is not positive", x[i][j] > 0.0);
            }
        }
        int[] labels = Datasets.labels();
        assertEquals(178, labels.length);
        for (int i = 0; i < labels.length; ++i) {
            assertTrue(labels[i] >= 0 && labels[i] < WineDemo.K);
        }
    }

    @Test
    public void accessorsHandOutCopies() {
        double[][] first = Datasets.features();
        assertNotSame(first, Datasets.features());
        double original = first[0][0];
        first[0][0] = -1.0;
        assertEquals(original, Datasets.features()[0][0], 0.0);

        int[] labels = Datasets.labels();
        labels[0] = 99;
        assertEquals(0, Datasets.labels()[0]);
    }

    @Test
    public void standardizingLeavesZeroMeanAndUnitVariance() {
        double[][] z = WineDemo.standardized();
        for (int j = 0; j < WineDemo.P; ++j) {
            double mean = 0.0;
            for (int i = 0; i < z.length; ++i) {
                mean += z[i][j];
            }
            mean /= z.length;
            double variance = 0.0;
            for (int i = 0; i < z.length; ++i) {
                variance += (z[i][j] - mean) * (z[i][j] - mean);
            }
            variance /= (z.length - 1);
            assertEquals("column " + j + " is not centred", 0.0, mean, 1.0e-12);
            assertEquals("column " + j + " is not scaled", 1.0, variance, 1.0e-12);
        }
    }

    /** On the raw table the leading component is one column in different units. */
    @Test
    public void rawComponentsAreProlineAlone() {
        WineDemo.Reduction raw = WineDemo.reduce(Datasets.features());
        assertTrue("PC1 explains " + raw.explainedRatio[0], raw.explainedRatio[0] > 0.99);
        assertEquals("proline", Datasets.featureName(raw.dominant));
        assertTrue("loading " + raw.dominantLoading, raw.dominantLoading > 0.999);
        assertTrue("next loading " + raw.nextLoading, raw.nextLoading < 0.05);
    }

    /** After standardizing no single column can carry the first component. */
    @Test
    public void standardizingSpreadsTheVariance() {
        WineDemo.Reduction standard = WineDemo.reduce(WineDemo.standardized());
        assertTrue("PC1 explains " + standard.explainedRatio[0],
                standard.explainedRatio[0] > 0.34 && standard.explainedRatio[0] < 0.39);
        double together = standard.explainedRatio[0] + standard.explainedRatio[1];
        assertTrue("PC1 and PC2 explain " + together, together > 0.53 && together < 0.58);
        double maxLoading = 0.0;
        for (int j = 0; j < WineDemo.P; ++j) {
            maxLoading = Math.max(maxLoading, Math.abs(standard.components[0][j]));
        }
        assertTrue("largest loading " + maxLoading, maxLoading < 0.5);

        double sum = 0.0;
        for (int k = 0; k < standard.explainedRatio.length; ++k) {
            sum += standard.explainedRatio[k];
        }
        assertEquals("the ratios must exhaust the variance", 1.0, sum, 1.0e-12);
    }

    /** Two exact routes and a randomized one, asked for the same two components. */
    @Test
    public void threeRoutesFindTheSameComponents() {
        WineDemo.Agreement agreement = WineDemo.pcaRoutes();
        assertTrue("components " + agreement.componentCovJacobi, agreement.componentCovJacobi < 1.0e-12);
        assertTrue("variance " + agreement.varianceCovJacobi, agreement.varianceCovJacobi < 1.0e-12);
        assertTrue("randomized components " + agreement.componentCovTruncated,
                agreement.componentCovTruncated < 1.0e-6);
        assertTrue("randomized variance " + agreement.varianceCovTruncated, agreement.varianceCovTruncated < 1.0e-12);
        assertTrue("eigenvalues " + agreement.eigenvalueAgreement, agreement.eigenvalueAgreement < 1.0e-12);
        assertTrue("the randomized route did not converge", agreement.truncatedConverged);
    }

    /** The analytic gradient against central differences, on a point nothing was fitted to. */
    @Test
    public void gradientMatchesCentralDifferences() {
        double[][] x = WineDemo.standardized();
        int[] y = Datasets.labels();
        WineDemo.MultinomialLogistic model = new WineDemo.MultinomialLogistic(x, y, WineDemo.L2);
        double[] w = new double[model.getNumParameters()];
        long lcg = 7L;
        for (int i = 0; i < w.length; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            w[i] = ((lcg >>> 11) * 0x1.0p-53) * 2.0 - 1.0;
        }
        model.setParameters(w);
        double[] gradient = new double[w.length];
        model.getValueGradient(gradient);

        double worst = 0.0;
        for (int i = 0; i < w.length; ++i) {
            double h = 1.0e-5;
            double[] point = w.clone();
            point[i] = w[i] + h;
            model.setParameters(point);
            double up = model.getValue();
            point[i] = w[i] - h;
            model.setParameters(point);
            double down = model.getValue();
            double difference = (up - down) / (2.0 * h);
            worst = Math.max(worst, Math.abs(difference - gradient[i]) / Math.max(1.0, Math.abs(gradient[i])));
        }
        assertTrue("worst relative deviation " + worst, worst < 1.0e-6);
    }

    /** The objective is a log-likelihood, so the optimizer must raise it. */
    @Test
    public void theOptimizerClimbs() {
        double[][] x = WineDemo.standardized();
        int[] y = Datasets.labels();
        WineDemo.MultinomialLogistic model = new WineDemo.MultinomialLogistic(x, y, WineDemo.L2);
        double atZero = model.getValue();
        new LimitedMemoryBFGS(model).optimize();
        double afterwards = model.getValue();
        assertTrue("value fell from " + atZero + " to " + afterwards, afterwards > atZero);
        assertTrue("a log-likelihood cannot be positive: " + afterwards, afterwards < 0.0);
    }

    /**
     * The point of section 4: without a penalty the maximum does not exist, the
     * coefficients grow with every tolerance that is tightened, and the
     * optimizer reports convergence throughout -- with
     * {@link Termination#GRADIENT_TOLERANCE}, the reason that sounds most like
     * an answer, because the gradient vanishes along the ray as well.
     */
    @Test
    public void unpenalizedSearchReportsConvergenceWithoutAMaximum() {
        WineDemo.Fit loose = WineDemo.fit(0.0);
        WineDemo.Fit tight = WineDemo.fitTightened(0.0);

        assertTrue("the loose search did not report convergence", loose.converged);
        assertTrue("the tightened search did not report convergence", tight.converged);
        assertEquals(Termination.GRADIENT_TOLERANCE, loose.termination);
        assertEquals(Termination.GRADIENT_TOLERANCE, tight.termination);
        assertTrue("a gradient exit must be stationary", loose.termination.isStationary());
        assertTrue("gradient norm " + loose.gradientNorm, loose.gradientNorm < 1.0e-3);
        assertTrue("gradient norm " + tight.gradientNorm, tight.gradientNorm < 1.0e-10);

        assertTrue("tightening moved ||w|| from " + loose.norm + " to " + tight.norm, tight.norm > 2.0 * loose.norm);
        assertTrue("the loose fit is not separating: " + loose.logLikelihoodPerWine,
                loose.logLikelihoodPerWine > -1.0e-4);
        assertTrue("the tight fit is not separating: " + tight.logLikelihoodPerWine,
                tight.logLikelihoodPerWine > -1.0e-10);
        assertEquals("a separating fit must classify its own data perfectly", 1.0, loose.accuracy, 0.0);
        assertEquals(1.0, tight.accuracy, 0.0);
    }

    /**
     * The control that makes the previous test mean something: on the penalized
     * problem, where the maximum exists, the same tightening leaves the
     * coefficients where they were.
     */
    @Test
    public void tighteningDoesNotMoveAFitThatExists() {
        WineDemo.Fit loose = WineDemo.fit(WineDemo.L2);
        WineDemo.Fit tight = WineDemo.fitTightened(WineDemo.L2);

        assertTrue(loose.converged);
        assertTrue(tight.converged);
        assertTrue("||w|| moved from " + loose.norm + " to " + tight.norm,
                Math.abs(tight.norm - loose.norm) / loose.norm < 0.01);
        assertTrue("the tightened run must come closer to stationarity: " + tight.gradientNorm + " against "
                + loose.gradientNorm, tight.gradientNorm < loose.gradientNorm);

        WineDemo.Fit unpenalizedLoose = WineDemo.fit(0.0);
        WineDemo.Fit unpenalizedTight = WineDemo.fitTightened(0.0);
        double penalizedDrift = Math.abs(tight.norm - loose.norm) / loose.norm;
        double unpenalizedDrift = Math.abs(unpenalizedTight.norm - unpenalizedLoose.norm) / unpenalizedLoose.norm;
        assertTrue("penalized drift " + penalizedDrift + " against unpenalized " + unpenalizedDrift,
                unpenalizedDrift > 100.0 * penalizedDrift);
    }

    /** With a penalty the same search settles at a point of finite size. */
    @Test
    public void thePenaltyBoundsTheCoefficients() {
        WineDemo.Fit unpenalized = WineDemo.fit(0.0);
        WineDemo.Fit penalized = WineDemo.fit(WineDemo.L2);
        assertTrue("the penalized search did not converge", penalized.converged);
        assertTrue("||w|| " + penalized.norm + " against " + unpenalized.norm,
                penalized.norm < unpenalized.norm / 4.0);
        assertTrue("accuracy " + penalized.accuracy, penalized.accuracy > 0.98);
        assertTrue("loglik/wine " + penalized.logLikelihoodPerWine, penalized.logLikelihoodPerWine < 0.0);
    }

    /** Three classifiers, several fold seeds, one floor that holds for all of them. */
    @Test
    public void crossValidationStaysAboveTheFloor() {
        for (int s = 0; s < FOLD_SEEDS.length; ++s) {
            long seed = (long) FOLD_SEEDS[s];
            WineDemo.CrossValidation cv = WineDemo.crossValidate(WineDemo.L2, seed);
            assertTrue("logistic " + cv.logistic + " at seed " + seed, cv.logistic > 0.95);
            assertTrue("lda " + cv.lda + " at seed " + seed, cv.lda > 0.95);
            assertTrue("qda " + cv.qda + " at seed " + seed, cv.qda > 0.95);

            assertEquals(178, cv.logisticHits.length);
            double sum = 0.0;
            for (int i = 0; i < cv.logisticHits.length; ++i) {
                assertTrue("an indicator must be 0 or 1", cv.logisticHits[i] == 0.0 || cv.logisticHits[i] == 1.0);
                sum += cv.logisticHits[i];
            }
            assertEquals("the reported accuracy must be the mean of the indicators it averages",
                    sum / cv.logisticHits.length, cv.logistic, 0.0);
        }
    }

    /** Every wine is tested exactly once, by a model that never saw it. */
    @Test
    public void everyWineIsHeldOutExactlyOnce() {
        int[] fold = WineDemo.foldOf(178, WineDemo.FOLDS, WineDemo.FOLD_SEED);
        int[] sizes = new int[WineDemo.FOLDS];
        for (int i = 0; i < fold.length; ++i) {
            assertTrue(fold[i] >= 0 && fold[i] < WineDemo.FOLDS);
            sizes[fold[i]]++;
        }
        for (int f = 0; f < WineDemo.FOLDS; ++f) {
            assertTrue("fold " + f + " holds " + sizes[f] + " wines", Math.abs(sizes[f] - 178 / WineDemo.FOLDS) <= 1);
        }
        int[] again = WineDemo.foldOf(178, WineDemo.FOLDS, WineDemo.FOLD_SEED);
        for (int i = 0; i < fold.length; ++i) {
            assertEquals("the same seed must give the same folds", fold[i], again[i]);
        }
        int[] other = WineDemo.foldOf(178, WineDemo.FOLDS, WineDemo.FOLD_SEED + 1L);
        boolean differs = false;
        for (int i = 0; i < fold.length; ++i) {
            differs |= fold[i] != other[i];
        }
        assertTrue("a different seed must give different folds", differs);
    }

    /**
     * The quadratic discriminant is admissible on this data -- the smallest
     * training class still has three times as many rows as columns -- so it is
     * reported rather than dropped.
     */
    @Test
    public void perClassCovariancesAreUsable() {
        WineDemo.CrossValidation cv = WineDemo.crossValidate(WineDemo.L2, WineDemo.FOLD_SEED);
        assertTrue("pooled condition " + cv.pooledCondition, cv.pooledCondition > 1.0 && cv.pooledCondition < 100.0);
        for (int k = 0; k < WineDemo.K; ++k) {
            assertTrue("class " + k + " condition " + cv.classCondition[k],
                    cv.classCondition[k] > 1.0 && cv.classCondition[k] < 1000.0);
        }
        assertTrue("the pooled covariance must be the better conditioned one",
                cv.pooledCondition < cv.classCondition[0]);
    }

    /**
     * Both L1 engines answer in exact zeros. The ROADMAP had OWL-QN reaching
     * one only by accident of an orthant crossing; the crossing is the
     * mechanism, and the counts below are the same whether they are taken
     * exactly or with a cutoff.
     */
    @Test
    public void bothL1EnginesProduceExactZeros() {
        WineDemo.Selection selection = WineDemo.selectFeatures(2);
        for (int t = 0; t < selection.l1Weights.length; ++t) {
            assertEquals("weight " + selection.l1Weights[t] + " has coefficients that are near zero but not zero",
                    selection.owlqnExactZeros[t], selection.owlqnNearZeros[t]);
            assertTrue("weight " + selection.l1Weights[t] + " zeroed nothing", selection.owlqnExactZeros[t] > 0);
            if (t > 0) {
                assertTrue("a heavier penalty must not keep more coefficients",
                        selection.owlqnExactZeros[t] >= selection.owlqnExactZeros[t - 1]);
                assertTrue("a heavier penalty must not enlarge the coefficients",
                        selection.owlqnNorm[t] < selection.owlqnNorm[t - 1]);
            }
        }
        assertTrue("the lasso zeroed nothing", selection.lassoExactZeros > 0);

        int droppedByBoth = 0;
        for (int j = 0; j < WineDemo.P; ++j) {
            boolean owlqn = true;
            boolean lasso = true;
            for (int k = 0; k < WineDemo.K; ++k) {
                owlqn &= selection.owlqnCoefficients[k][j] == 0.0;
                lasso &= selection.lassoCoefficients[k][j] == 0.0;
            }
            if (owlqn && lasso) {
                droppedByBoth++;
            }
        }
        assertTrue("the two engines agree on nothing", droppedByBoth > 0);
    }

    /**
     * The claim section 3 and section 8 both make: separable in all thirteen
     * measurements, not separable in the two the plot shows.
     */
    @Test
    public void theCultivarsAreSeparableInThirteenDimensionsAndNotInTwo() {
        WineDemo.Separability separable = WineDemo.separability();

        assertEquals("13 dimensions must separate the cultivars exactly", 0, separable.wrongInFullSpace);
        assertEquals("a separating fit is certain about every wine", 1.0, separable.smallestMargin, 1.0e-9);
        assertTrue("the plotted plane must not separate them: " + separable.wrongInThePlane,
                separable.wrongInThePlane > 0);
        assertTrue("the plane carries about half the variance: " + separable.plottedVariance,
                separable.plottedVariance > 0.5 && separable.plottedVariance < 0.6);
    }

    /**
     * OWL-QN stops on the value rule here, well short of the L1 optimality
     * condition, and tightening the tolerances that {@code math.optim} now
     * exposes moves the coefficients while leaving the selection almost alone.
     */
    @Test
    public void theL1SelectionIsSharperThanItsCoefficients() {
        WineDemo.Selection selection = WineDemo.selectFeatures(2);

        for (int t = 0; t < selection.l1Weights.length; ++t) {
            assertEquals("weight " + selection.l1Weights[t], Termination.VALUE_TOLERANCE,
                    selection.owlqnTermination[t]);
            assertTrue("pseudo gradient " + selection.owlqnPseudoGradient[t],
                    selection.owlqnPseudoGradient[t] > 1.0e-2);
        }
        assertTrue("tightening must come closer to optimality: " + selection.tightenedPseudoGradient,
                selection.tightenedPseudoGradient < 0.01 * selection.owlqnPseudoGradient[2]);
        assertTrue("tightening must cost iterations",
                selection.tightenedIterations > selection.owlqnIterations[2]);
        assertTrue("the support must be nearly stable: " + selection.supportChanges + " changes",
                selection.supportChanges <= 3);
        assertTrue("the coefficients must not be: " + selection.largestCoefficientChange,
                selection.largestCoefficientChange > 1.0e-3);
    }

    /** A seeded bootstrap has to reproduce exactly, or the demo cannot print it. */
    @Test
    public void theBootstrapReproduces() {
        WineDemo.CrossValidation cv = WineDemo.crossValidate(WineDemo.L2, WineDemo.FOLD_SEED);
        WineDemo.Interval first = WineDemo.interval(cv.logisticHits, 0.95);
        WineDemo.Interval second = WineDemo.interval(cv.logisticHits, 0.95);

        assertEquals(first.standardError, second.standardError, 0.0);
        assertEquals(first.percentile[0], second.percentile[0], 0.0);
        assertEquals(first.percentile[1], second.percentile[1], 0.0);
        assertEquals(first.bca[0], second.bca[0], 0.0);
        assertEquals(first.bca[1], second.bca[1], 0.0);

        assertTrue("the interval must contain the estimate",
                first.percentile[0] <= first.accuracy && first.accuracy <= first.percentile[1]);
        assertTrue("a share cannot leave [0, 1]", first.bca[0] >= 0.0 && first.bca[1] <= 1.0);
        assertTrue("the interval has no width", first.percentile[1] > first.percentile[0]);
        assertTrue("the bootstrap mean is far from the estimate: " + first.bootstrapMean,
                Math.abs(first.bootstrapMean - first.accuracy) < 0.01);
    }

    /** Everything the demo prints has to come out the same way twice. */
    @Test
    public void theDemoPrintsTheSameThingTwice() {
        String first = run();
        String second = run();
        assertEquals(first, second);
        assertTrue("the demo printed almost nothing", first.length() > 2000);
        assertFalse("a locale slipped into the output", first.contains("0,9"));
        assertTrue(first.contains("proline"));
        assertTrue("the demo must say what it established", first.contains("what this run established"));
        assertTrue("and it must not leave the separability implicit",
                first.contains("linearly separable"));
    }

    private static String run() {
        PrintStream out = System.out;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(buffer, true, "UTF-8"));
            WineDemo.main(new String[0]);
        } catch (UnsupportedEncodingException e) {
            throw new AssertionError(e);
        } finally {
            System.setOut(out);
        }
        try {
            return buffer.toString("UTF-8");
        } catch (UnsupportedEncodingException e) {
            throw new AssertionError(e);
        }
    }
}
