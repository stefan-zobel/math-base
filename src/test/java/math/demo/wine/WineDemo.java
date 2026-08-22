package math.demo.wine;

import java.util.Arrays;
import java.util.Locale;

import math.dl.MultiClassAccuracy;
import math.dl.Softmax;
import math.linalg.CholeskyDecomp;
import math.linalg.CovariancePCA;
import math.linalg.DMatrix;
import math.linalg.JacobiPCA;
import math.linalg.Lasso;
import math.linalg.Standardization;
import math.linalg.SymmetricJacobiEigen;
import math.linalg.TruncatedPCA;
import math.optim.LimitedMemoryBFGS;
import math.optim.Optimizable;
import math.optim.OrthantWiseLimitedMemoryBFGS;
import math.optim.Termination;
import math.probe.Bootstrap;
import math.probe.CovarianceAccumulator;
import math.probe.SampleStatistic;

/**
 * A worked example: 178 wines, thirteen chemical measurements and three
 * cultivars, classified with the pieces this library provides.
 * <p>
 * The question is the one most often asked of a numerical library -- can it do
 * machine learning -- and answering it needs five packages: {@code math.probe}
 * for the moments, {@code math.linalg} for the components and the lasso,
 * {@code math.optim} for the optimizer, {@code math.dl} for the softmax and
 * the accuracy, and {@code math.probe} again for an interval around the
 * result.
 * <p>
 * Three things in here are worth more than the classification, which is easy:
 * <ul>
 * <li><b>Standardizing is not a ritual.</b> On the raw table the first
 * principal component explains 99.8 percent of the variance and is proline
 * alone, because proline is measured in the hundreds and everything else
 * between 0.1 and 30. The demo runs the reduction both ways to show it.</li>
 * <li><b>{@link LimitedMemoryBFGS} maximizes.</b> {@link MultinomialLogistic}
 * is therefore written as a log-likelihood, not as a loss, and its javadoc
 * says so; getting that sign backwards is the single most likely way to
 * misuse {@link Optimizable.ByGradientValue}.</li>
 * <li><b>A converged flag is not a maximum, and neither is a small
 * gradient.</b> Wine is linearly separable after standardizing, so the
 * unpenalized maximum likelihood does not exist -- the coefficients grow
 * without bound. The optimizer reports convergence, and it reports
 * {@link Termination#GRADIENT_TOLERANCE} while doing so, because the gradient
 * vanishes along the ray the coefficients run off on. What gives it away is
 * running the same fit with tighter tolerances and watching the coefficients
 * move; section 4 does that, with the penalized fit beside it as the control
 * that does not move.</li>
 * </ul>
 * <p>
 * Everything is seeded: the fold permutation, the randomized PCA, the lasso
 * cross-validation and the bootstrap. Two runs of {@link #main(String[])}
 * produce identical output, on either source tree.
 */
public final class WineDemo {

    /** Formatting locale, so that the output does not depend on the machine. */
    private static final Locale L = Locale.ROOT;

    /** Number of measurements per wine. */
    static final int P = Datasets.FEATURES_PER_WINE;

    /** Number of cultivars. */
    static final int K = Datasets.CLASSES;

    /** Ridge penalty of the logistic fit, chosen in section 4. */
    static final double L2 = 1.0;

    /** Number of cross validation folds. */
    static final int FOLDS = 5;

    /** Seed of the fold permutation. */
    static final long FOLD_SEED = 20260821L;

    /** Seed of the randomized PCA. */
    static final long PCA_SEED = 42L;

    /** Seed of the lasso cross validation. */
    static final long LASSO_SEED = 12345L;

    /** Seed of the bootstrap. */
    static final long BOOT_SEED = 4242L;

    /** Bootstrap replications. */
    static final int BOOT_ITERATIONS = 20000;

    /** L1 weights section 5 walks through. */
    static final double[] L1_WEIGHTS = { 0.1, 1.0, 5.0, 20.0, 50.0 };

    /** The arithmetic mean, as a bootstrap statistic. */
    private static final SampleStatistic MEAN = new SampleStatistic() {
        @Override
        public double apply(double[] sample) {
            double sum = 0.0;
            for (int i = 0; i < sample.length; ++i) {
                sum += sample[i];
            }
            return sum / sample.length;
        }
    };

    // ---------------------------------------------------------------- 1. data

    /** What the table looks like before anything is fitted to it. */
    static final class Description {
        final double[] min;
        final double[] max;
        final double[] mean;
        final double[] sd;
        final int[] counts;
        /** Column with the smallest mean, and the one with the largest. */
        final int smallest;
        final int largest;

        Description(double[] min, double[] max, double[] mean, double[] sd, int[] counts, int smallest, int largest) {
            this.min = min;
            this.max = max;
            this.mean = mean;
            this.sd = sd;
            this.counts = counts;
            this.smallest = smallest;
            this.largest = largest;
        }

        /** Ratio of the largest mean to the smallest, the reason to standardize. */
        double scaleSpread() {
            return Math.abs(mean[largest]) / Math.abs(mean[smallest]);
        }
    }

    static Description describe() {
        double[][] x = Datasets.features();
        CovarianceAccumulator acc = new CovarianceAccumulator(P);
        acc.addAll(x);
        double[] mean = acc.mean();
        double[] var = acc.variance();
        double[] sd = new double[P];
        double[] min = new double[P];
        double[] max = new double[P];
        Arrays.fill(min, Double.MAX_VALUE);
        Arrays.fill(max, -Double.MAX_VALUE);
        for (int i = 0; i < x.length; ++i) {
            for (int j = 0; j < P; ++j) {
                min[j] = Math.min(min[j], x[i][j]);
                max[j] = Math.max(max[j], x[i][j]);
            }
        }
        int smallest = 0;
        int largest = 0;
        for (int j = 0; j < P; ++j) {
            sd[j] = Math.sqrt(var[j]);
            if (Math.abs(mean[j]) < Math.abs(mean[smallest])) {
                smallest = j;
            }
            if (Math.abs(mean[j]) > Math.abs(mean[largest])) {
                largest = j;
            }
        }
        return new Description(min, max, mean, sd, Datasets.classCounts(), smallest, largest);
    }

    /** The table with every column centered and scaled to unit variance. */
    static double[][] standardized() {
        double[][] x = Datasets.features();
        return Standardization.of(x).transform(x);
    }

    // ----------------------------------------------------------- 2. and 3. PCA

    /** A principal component analysis reduced to what this demo asks of it. */
    static final class Reduction {
        final double[] explainedRatio;
        final double[][] components;
        final double[][] scores;
        final int dominant;
        final double dominantLoading;
        final double nextLoading;

        Reduction(double[] explainedRatio, double[][] components, double[][] scores, int dominant,
                double dominantLoading, double nextLoading) {
            this.explainedRatio = explainedRatio;
            this.components = components;
            this.scores = scores;
            this.dominant = dominant;
            this.dominantLoading = dominantLoading;
            this.nextLoading = nextLoading;
        }
    }

    /** Principal components of {@code x}, keeping all of them so the ratios are complete. */
    static Reduction reduce(double[][] x) {
        JacobiPCA pca = new JacobiPCA();
        double[][] scores = pca.pca(x, P);
        double[] variance = pca.getExplainedVariance();
        double total = 0.0;
        for (int k = 0; k < variance.length; ++k) {
            total += variance[k];
        }
        double[] ratio = new double[variance.length];
        for (int k = 0; k < variance.length; ++k) {
            ratio[k] = variance[k] / total;
        }
        double[][] components = pca.getComponents();
        int dominant = 0;
        for (int j = 1; j < P; ++j) {
            if (Math.abs(components[0][j]) > Math.abs(components[0][dominant])) {
                dominant = j;
            }
        }
        double next = 0.0;
        for (int j = 0; j < P; ++j) {
            if (j != dominant) {
                next = Math.max(next, Math.abs(components[0][j]));
            }
        }
        return new Reduction(ratio, components, scores, dominant, Math.abs(components[0][dominant]), next);
    }

    /** Three routes to the same two components, and how far apart they land. */
    static final class Agreement {
        final double[] covarianceVariance;
        final double[] jacobiVariance;
        final double[] truncatedVariance;
        final double varianceCovJacobi;
        final double varianceCovTruncated;
        final double componentCovJacobi;
        final double componentCovTruncated;
        final double eigenvalueAgreement;
        final boolean truncatedConverged;
        final int truncatedIterations;

        Agreement(double[] covarianceVariance, double[] jacobiVariance, double[] truncatedVariance,
                double varianceCovJacobi, double varianceCovTruncated, double componentCovJacobi,
                double componentCovTruncated, double eigenvalueAgreement, boolean truncatedConverged,
                int truncatedIterations) {
            this.covarianceVariance = covarianceVariance;
            this.jacobiVariance = jacobiVariance;
            this.truncatedVariance = truncatedVariance;
            this.varianceCovJacobi = varianceCovJacobi;
            this.varianceCovTruncated = varianceCovTruncated;
            this.componentCovJacobi = componentCovJacobi;
            this.componentCovTruncated = componentCovTruncated;
            this.eigenvalueAgreement = eigenvalueAgreement;
            this.truncatedConverged = truncatedConverged;
            this.truncatedIterations = truncatedIterations;
        }
    }

    /**
     * The exact eigen route, the exact SVD route and the randomized route, all
     * asked for the same two components of the standardized table.
     */
    static Agreement pcaRoutes() {
        double[][] x = standardized();
        CovarianceAccumulator acc = new CovarianceAccumulator(P);
        acc.addAll(x);
        double[] covariance = acc.covariance();

        CovariancePCA fromCovariance = new CovariancePCA().fit(covariance, P, acc.mean(), 2);
        JacobiPCA fromMatrix = new JacobiPCA();
        fromMatrix.pca(x, 2);
        TruncatedPCA randomized = new TruncatedPCA(10, 1.0e-9, 200, PCA_SEED);
        randomized.pca(x, 2);

        double[] cv = fromCovariance.getExplainedVariance();
        double[] jv = fromMatrix.getExplainedVariance();
        double[] tv = new double[2];
        double[] tsv = randomized.getSingularValues();
        for (int k = 0; k < 2; ++k) {
            tv[k] = tsv[k] * tsv[k] / (x.length - 1);
        }

        SymmetricJacobiEigen.Result eigen = new SymmetricJacobiEigen().decompose(covariance, P);
        double[] sv = fromMatrix.getSingularValues();
        double eigenAgreement = 0.0;
        for (int k = 0; k < sv.length; ++k) {
            eigenAgreement = Math.max(eigenAgreement, relative(sv[k] * sv[k] / (x.length - 1), eigen.lambda[k]));
        }

        return new Agreement(cv, jv, tv, Math.max(relative(cv[0], jv[0]), relative(cv[1], jv[1])),
                Math.max(relative(cv[0], tv[0]), relative(cv[1], tv[1])),
                componentDistance(fromCovariance.getComponents(), fromMatrix.getComponents()),
                componentDistance(fromCovariance.getComponents(), randomized.getComponents()), eigenAgreement,
                randomized.converged(), randomized.getIterations());
    }

    /** Largest coefficient difference between two sets of components, each taken up to sign. */
    static double componentDistance(double[][] a, double[][] b) {
        double worst = 0.0;
        for (int k = 0; k < Math.min(a.length, b.length); ++k) {
            double dot = 0.0;
            for (int j = 0; j < a[k].length; ++j) {
                dot += a[k][j] * b[k][j];
            }
            double sign = dot < 0.0 ? -1.0 : 1.0;
            for (int j = 0; j < a[k].length; ++j) {
                worst = Math.max(worst, Math.abs(a[k][j] - sign * b[k][j]));
            }
        }
        return worst;
    }

    static double relative(double a, double b) {
        double scale = Math.max(Math.abs(a), Math.abs(b));
        return scale == 0.0 ? Math.abs(a - b) : Math.abs(a - b) / scale;
    }

    /** How separable the cultivars are, in the full space and in the plotted plane. */
    static final class Separability {
        final int wrongInFullSpace;
        final int wrongInThePlane;
        final double smallestMargin;
        final double plottedVariance;

        Separability(int wrongInFullSpace, int wrongInThePlane, double smallestMargin, double plottedVariance) {
            this.wrongInFullSpace = wrongInFullSpace;
            this.wrongInThePlane = wrongInThePlane;
            this.smallestMargin = smallestMargin;
            this.plottedVariance = plottedVariance;
        }
    }

    /**
     * The same classifier fitted twice, on all thirteen measurements and on the
     * two components the scatter above plots, both scored on the data they were
     * fitted to. The difference is what the picture costs.
     */
    static Separability separability() {
        double[][] x = standardized();
        int[] y = Datasets.labels();
        Reduction reduction = reduce(x);
        double[][] plane = new double[x.length][2];
        for (int i = 0; i < x.length; ++i) {
            plane[i][0] = reduction.scores[i][0];
            plane[i][1] = reduction.scores[i][1];
        }
        return new Separability(wrongInSample(x, y), wrongInSample(plane, y), margin(x, y),
                reduction.explainedRatio[0] + reduction.explainedRatio[1]);
    }

    /** Wines a separating fit gets wrong on the data it was fitted to. */
    private static int wrongInSample(double[][] x, int[] y) {
        MultinomialLogistic model = new MultinomialLogistic(x, y, 0.0);
        new LimitedMemoryBFGS(model, 10000, 1.0e-12, 1.0e-12, 5, 1.0e-14, 1.0e-14).optimize();
        int wrong = 0;
        double[] probability = new double[K];
        for (int i = 0; i < x.length; ++i) {
            model.predict(x[i], probability);
            int best = 0;
            for (int k = 1; k < K; ++k) {
                if (probability[k] > probability[best]) {
                    best = k;
                }
            }
            if (best != y[i]) {
                wrong++;
            }
        }
        return wrong;
    }

    /** Smallest gap between the probability of the right cultivar and the next one. */
    private static double margin(double[][] x, int[] y) {
        MultinomialLogistic model = new MultinomialLogistic(x, y, 0.0);
        new LimitedMemoryBFGS(model, 10000, 1.0e-12, 1.0e-12, 5, 1.0e-14, 1.0e-14).optimize();
        double smallest = Double.MAX_VALUE;
        double[] probability = new double[K];
        for (int i = 0; i < x.length; ++i) {
            model.predict(x[i], probability);
            double other = 0.0;
            for (int k = 0; k < K; ++k) {
                if (k != y[i]) {
                    other = Math.max(other, probability[k]);
                }
            }
            smallest = Math.min(smallest, probability[y[i]] - other);
        }
        return smallest;
    }

    // ------------------------------------------------------- 4. the classifier

    /**
     * Multinomial logistic regression over {@code K} classes, written against
     * {@link Optimizable.ByGradientValue}.
     * <p>
     * <b>{@link LimitedMemoryBFGS} and {@link OrthantWiseLimitedMemoryBFGS}
     * both maximize</b> -- the first search direction is the gradient itself.
     * {@link #getValue()} therefore returns the penalized <em>log-likelihood</em>
     * and {@link #getValueGradient(double[])} its gradient, both with the sign
     * that makes the optimizer climb. Handing either of them a cross-entropy
     * loss instead makes the search run away from the answer, and nothing in
     * the two APIs says which one they want.
     * <p>
     * The parameters are {@code K x (P + 1)}, class-major, the intercept last
     * in each class and left out of the penalty. That is one class more than
     * the model needs -- adding a constant to every class leaves the
     * likelihood unchanged -- and the ridge term is what picks one of those
     * equivalent solutions.
     */
    static final class MultinomialLogistic implements Optimizable.ByGradientValue {

        private final double[][] x;
        private final int[] y;
        private final double l2;
        /** Columns of {@code x}, so that the model can also be fitted to a reduction of it. */
        private final int p;
        /** Parameters per class: the features plus one intercept. */
        private final int d;
        private final double[] w;
        private final double[] score = new double[K];
        private final double[] probability = new double[K];

        MultinomialLogistic(double[][] x, int[] y, double l2) {
            this.x = x;
            this.y = y;
            this.l2 = l2;
            this.p = x[0].length;
            this.d = p + 1;
            this.w = new double[K * d];
        }

        /** Linear scores of one wine, before the softmax. */
        void score(double[] wine, double[] out) {
            for (int k = 0; k < K; ++k) {
                double s = w[k * d + p];
                for (int j = 0; j < p; ++j) {
                    s += w[k * d + j] * wine[j];
                }
                out[k] = s;
            }
        }

        /** Class probabilities of one wine, written into {@code out}. */
        double[] predict(double[] wine, double[] out) {
            score(wine, score);
            return Softmax.softmax(K, 0, score, 0, out);
        }

        @Override
        public double getValue() {
            double logLikelihood = 0.0;
            for (int i = 0; i < x.length; ++i) {
                score(x[i], score);
                // log sum exp, so that a confident fit does not underflow
                double max = score[0];
                for (int k = 1; k < K; ++k) {
                    max = Math.max(max, score[k]);
                }
                double sum = 0.0;
                for (int k = 0; k < K; ++k) {
                    sum += Math.exp(score[k] - max);
                }
                logLikelihood += score[y[i]] - max - Math.log(sum);
            }
            if (l2 > 0.0) {
                double squares = 0.0;
                for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < p; ++j) {
                        squares += w[k * d + j] * w[k * d + j];
                    }
                }
                logLikelihood -= 0.5 * l2 * squares;
            }
            return logLikelihood;
        }

        @Override
        public void getValueGradient(double[] buffer) {
            Arrays.fill(buffer, 0.0);
            for (int i = 0; i < x.length; ++i) {
                score(x[i], score);
                Softmax.softmax(K, 0, score, 0, probability);
                for (int k = 0; k < K; ++k) {
                    double residual = (k == y[i] ? 1.0 : 0.0) - probability[k];
                    for (int j = 0; j < p; ++j) {
                        buffer[k * d + j] += residual * x[i][j];
                    }
                    buffer[k * d + p] += residual;
                }
            }
            if (l2 > 0.0) {
                for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < p; ++j) {
                        buffer[k * d + j] -= l2 * w[k * d + j];
                    }
                }
            }
        }

        @Override
        public int getNumParameters() {
            return w.length;
        }

        @Override
        public void getParameters(double[] buffer) {
            System.arraycopy(w, 0, buffer, 0, w.length);
        }

        @Override
        public double getParameter(int index) {
            return w[index];
        }

        @Override
        public void setParameters(double[] parameters) {
            System.arraycopy(parameters, 0, w, 0, w.length);
        }

        @Override
        public void setParameter(int index, double value) {
            w[index] = value;
        }

        /** Coefficient of feature {@code j} for class {@code k}, intercept excluded. */
        double coefficient(int k, int j) {
            return w[k * d + j];
        }

        /** Euclidean norm of the penalized coefficients, the intercepts left out. */
        double norm() {
            double sum = 0.0;
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < p; ++j) {
                    sum += w[k * d + j] * w[k * d + j];
                }
            }
            return Math.sqrt(sum);
        }

        /** How many of the penalized coefficients, the intercepts aside, are exactly zero. */
        int exactZeros() {
            int zeros = 0;
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < p; ++j) {
                    if (coefficient(k, j) == 0.0) {
                        zeros++;
                    }
                }
            }
            return zeros;
        }

        /** How many of them are smaller than {@code eps} in absolute value. */
        int nearZeros(double eps) {
            int zeros = 0;
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < p; ++j) {
                    if (Math.abs(coefficient(k, j)) < eps) {
                        zeros++;
                    }
                }
            }
            return zeros;
        }
    }

    /** One trained logistic model and what became of the search that produced it. */
    static final class Fit {
        final MultinomialLogistic model;
        final boolean converged;
        final Termination termination;
        final double gradientNorm;
        final int iterations;
        final double norm;
        final double logLikelihoodPerWine;
        final double accuracy;

        Fit(MultinomialLogistic model, boolean converged, Termination termination, double gradientNorm, int iterations,
                double norm, double logLikelihoodPerWine, double accuracy) {
            this.model = model;
            this.converged = converged;
            this.termination = termination;
            this.gradientNorm = gradientNorm;
            this.iterations = iterations;
            this.norm = norm;
            this.logLikelihoodPerWine = logLikelihoodPerWine;
            this.accuracy = accuracy;
        }
    }

    /** Fits the model on all 178 wines with the default stopping rules. */
    static Fit fit(double l2) {
        return run(new LimitedMemoryBFGSFactory() {
            @Override
            public LimitedMemoryBFGS create(MultinomialLogistic model) {
                return new LimitedMemoryBFGS(model);
            }
        }, l2);
    }

    /**
     * The same fit with every tolerance tightened and the budget enlarged.
     * Where the maximum exists the coefficients stay where they were; where it
     * does not, they grow.
     */
    static Fit fitTightened(double l2) {
        return run(new LimitedMemoryBFGSFactory() {
            @Override
            public LimitedMemoryBFGS create(MultinomialLogistic model) {
                return new LimitedMemoryBFGS(model, 10000, 1.0e-12, 1.0e-12, 5, 1.0e-14, 1.0e-14);
            }
        }, l2);
    }

    /** How one of the two runs above builds its optimizer. */
    interface LimitedMemoryBFGSFactory {
        LimitedMemoryBFGS create(MultinomialLogistic model);
    }

    private static Fit run(LimitedMemoryBFGSFactory factory, double l2) {
        double[][] x = standardized();
        int[] y = Datasets.labels();
        MultinomialLogistic model = new MultinomialLogistic(x, y, l2);
        LimitedMemoryBFGS optimizer = factory.create(model);
        boolean converged = optimizer.optimize();
        return new Fit(model, converged, optimizer.getTermination(), optimizer.getGradientNorm(),
                optimizer.getIteration(), model.norm(), model.getValue() / x.length, accuracy(model, x, y));
    }

    /** Share of wines whose most probable class is the observed one. */
    static double accuracy(MultinomialLogistic model, double[][] x, int[] y) {
        MultiClassAccuracy measure = new MultiClassAccuracy();
        double[] probability = new double[K];
        double[] observed = new double[K];
        for (int i = 0; i < x.length; ++i) {
            model.predict(x[i], probability);
            Arrays.fill(observed, 0.0);
            observed[y[i]] = 1.0;
            measure.compare(K, probability, 0, observed);
        }
        return measure.getAverage();
    }

    // -------------------------------------------------------- 5. L1 selection

    /** Which measurements survive an L1 penalty, asked of two different engines. */
    static final class Selection {
        final double[] l1Weights;
        final int[] owlqnExactZeros;
        final int[] owlqnNearZeros;
        final double[] owlqnNorm;
        final double[] owlqnAccuracy;
        final Termination[] owlqnTermination;
        final double[] owlqnPseudoGradient;
        final int[] owlqnIterations;
        final double[][] owlqnCoefficients;
        final double[][] lassoCoefficients;
        final int lassoExactZeros;
        final double[] lassoLambda;
        /** The reported weight, refitted with the tolerances tightened. */
        final Termination tightenedTermination;
        final double tightenedPseudoGradient;
        final int tightenedIterations;
        final int tightenedExactZeros;
        final double largestCoefficientChange;
        final int supportChanges;

        Selection(double[] l1Weights, int[] owlqnExactZeros, int[] owlqnNearZeros, double[] owlqnNorm,
                double[] owlqnAccuracy, Termination[] owlqnTermination, double[] owlqnPseudoGradient,
                int[] owlqnIterations,
                double[][] owlqnCoefficients, double[][] lassoCoefficients, int lassoExactZeros, double[] lassoLambda,
                Termination tightenedTermination, double tightenedPseudoGradient, int tightenedIterations,
                int tightenedExactZeros, double largestCoefficientChange, int supportChanges) {
            this.l1Weights = l1Weights;
            this.owlqnExactZeros = owlqnExactZeros;
            this.owlqnNearZeros = owlqnNearZeros;
            this.owlqnNorm = owlqnNorm;
            this.owlqnAccuracy = owlqnAccuracy;
            this.owlqnTermination = owlqnTermination;
            this.owlqnPseudoGradient = owlqnPseudoGradient;
            this.owlqnIterations = owlqnIterations;
            this.owlqnCoefficients = owlqnCoefficients;
            this.lassoCoefficients = lassoCoefficients;
            this.lassoExactZeros = lassoExactZeros;
            this.lassoLambda = lassoLambda;
            this.tightenedTermination = tightenedTermination;
            this.tightenedPseudoGradient = tightenedPseudoGradient;
            this.tightenedIterations = tightenedIterations;
            this.tightenedExactZeros = tightenedExactZeros;
            this.largestCoefficientChange = largestCoefficientChange;
            this.supportChanges = supportChanges;
        }
    }

    /**
     * OWL-QN on the multinomial log-likelihood, and the lasso on a one against
     * the rest indicator per class. Both engines produce exact zeros: OWL-QN
     * because a step that would change a coefficient's sign is clamped to the
     * orthant boundary, the lasso because the soft threshold sets the
     * coordinate to zero outright.
     */
    static Selection selectFeatures(int reportAt) {
        double[][] x = standardized();
        int[] y = Datasets.labels();

        int[] exact = new int[L1_WEIGHTS.length];
        int[] near = new int[L1_WEIGHTS.length];
        double[] norm = new double[L1_WEIGHTS.length];
        double[] accuracy = new double[L1_WEIGHTS.length];
        Termination[] termination = new Termination[L1_WEIGHTS.length];
        int[] iterations = new int[L1_WEIGHTS.length];
        double[] pseudoGradient = new double[L1_WEIGHTS.length];
        double[][] owlqn = new double[K][P];
        for (int t = 0; t < L1_WEIGHTS.length; ++t) {
            MultinomialLogistic model = new MultinomialLogistic(x, y, 0.0);
            OrthantWiseLimitedMemoryBFGS optimizer = new OrthantWiseLimitedMemoryBFGS(model, L1_WEIGHTS[t]);
            optimizer.optimize();
            exact[t] = model.exactZeros();
            near[t] = model.nearZeros(1.0e-6);
            norm[t] = model.norm();
            accuracy[t] = accuracy(model, x, y);
            termination[t] = optimizer.getTermination();
            pseudoGradient[t] = optimizer.getGradientNorm();
            iterations[t] = optimizer.getIteration();
            if (t == reportAt) {
                for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < P; ++j) {
                        owlqn[k][j] = model.coefficient(k, j);
                    }
                }
            }
        }

        // the reported weight again, with the tolerances the class now exposes
        MultinomialLogistic sharpened = new MultinomialLogistic(x, y, 0.0);
        OrthantWiseLimitedMemoryBFGS sharp = new OrthantWiseLimitedMemoryBFGS(sharpened, L1_WEIGHTS[reportAt], 10000,
                1.0e-12, 1.0e-9, 4);
        sharp.optimize();
        double largestChange = 0.0;
        int supportChanges = 0;
        for (int k = 0; k < K; ++k) {
            for (int j = 0; j < P; ++j) {
                largestChange = Math.max(largestChange, Math.abs(owlqn[k][j] - sharpened.coefficient(k, j)));
                if ((owlqn[k][j] == 0.0) != (sharpened.coefficient(k, j) == 0.0)) {
                    supportChanges++;
                }
            }
        }

        DMatrix design = new DMatrix(x.length, P);
        for (int i = 0; i < x.length; ++i) {
            for (int j = 0; j < P; ++j) {
                design.setUnsafe(i, j, x[i][j]);
            }
        }
        double[][] lasso = new double[K][];
        double[] lambda = new double[K];
        int lassoZeros = 0;
        for (int k = 0; k < K; ++k) {
            DMatrix indicator = new DMatrix(x.length, 1);
            for (int i = 0; i < x.length; ++i) {
                indicator.setUnsafe(i, 0, y[i] == k ? 1.0 : 0.0);
            }
            Lasso.CvResult cv = Lasso.cv(design, indicator, 1.0, FOLDS, LASSO_SEED);
            Lasso.Result at1se = Lasso.estimate(design, indicator, cv.lambda1se);
            lasso[k] = at1se.beta.clone();
            lambda[k] = cv.lambda1se;
            lassoZeros += P - at1se.nonZeros;
        }
        return new Selection(L1_WEIGHTS.clone(), exact, near, norm, accuracy, termination, pseudoGradient, iterations,
                owlqn, lasso,
                lassoZeros, lambda, sharp.getTermination(), sharp.getGradientNorm(), sharp.getIteration(),
                sharpened.exactZeros(), largestChange, supportChanges);
    }

    // ------------------------------------------------------- 6. generative model

    /**
     * A Gaussian classifier fitted by moments rather than by an optimizer:
     * shared covariance is the linear discriminant, one covariance per class
     * the quadratic one.
     */
    static final class Gaussian {
        private final double[][] means;
        private final DMatrix[] factors;
        private final double[] logDeterminant;
        private final double[] logPrior;
        final double[] condition;

        Gaussian(double[][] means, DMatrix[] factors, double[] logDeterminant, double[] logPrior,
                double[] condition) {
            this.means = means;
            this.factors = factors;
            this.logDeterminant = logDeterminant;
            this.logPrior = logPrior;
            this.condition = condition;
        }

        /** Class probabilities of one wine, written into {@code out}. */
        double[] predict(double[] wine, double[] out) {
            double[] score = new double[K];
            double[] centred = new double[P];
            for (int k = 0; k < K; ++k) {
                for (int j = 0; j < P; ++j) {
                    centred[j] = wine[j] - means[k][j];
                }
                double[] z = forwardSolve(factors[k], centred);
                double quadratic = 0.0;
                for (int j = 0; j < P; ++j) {
                    quadratic += z[j] * z[j];
                }
                score[k] = -0.5 * quadratic - 0.5 * logDeterminant[k] + logPrior[k];
            }
            return Softmax.softmax(K, 0, score, 0, out);
        }
    }

    /** Fits the Gaussian classifier; {@code shared} selects LDA over QDA. */
    static Gaussian gaussian(double[][] x, int[] y, boolean shared) {
        double[][] means = new double[K][];
        double[][] classCovariance = new double[K][];
        int[] counts = new int[K];
        double[] pooled = new double[P * P];
        for (int k = 0; k < K; ++k) {
            CovarianceAccumulator acc = new CovarianceAccumulator(P);
            for (int i = 0; i < x.length; ++i) {
                if (y[i] == k) {
                    acc.add(x[i]);
                    counts[k]++;
                }
            }
            means[k] = acc.mean();
            classCovariance[k] = acc.covariance();
            for (int e = 0; e < pooled.length; ++e) {
                pooled[e] += (counts[k] - 1) * classCovariance[k][e];
            }
        }
        for (int e = 0; e < pooled.length; ++e) {
            pooled[e] /= (x.length - K);
        }

        DMatrix[] factors = new DMatrix[K];
        double[] logDeterminant = new double[K];
        double[] logPrior = new double[K];
        double[] condition = new double[shared ? 1 : K];
        DMatrix pooledFactor = shared ? CholeskyDecomp.cholesky(square(pooled)) : null;
        if (shared) {
            condition[0] = condition(pooled);
        }
        for (int k = 0; k < K; ++k) {
            logPrior[k] = Math.log(counts[k] / (double) x.length);
            if (shared) {
                factors[k] = pooledFactor;
            } else {
                factors[k] = CholeskyDecomp.cholesky(square(classCovariance[k]));
                condition[k] = condition(classCovariance[k]);
            }
            double sum = 0.0;
            for (int j = 0; j < P; ++j) {
                sum += Math.log(factors[k].getUnsafe(j, j));
            }
            logDeterminant[k] = 2.0 * sum;
        }
        return new Gaussian(means, factors, logDeterminant, logPrior, condition);
    }

    /** Ratio of the largest to the smallest eigenvalue of a symmetric matrix. */
    static double condition(double[] symmetric) {
        SymmetricJacobiEigen.Result result = new SymmetricJacobiEigen().decompose(symmetric, P);
        return result.lambda[0] / result.lambda[P - 1];
    }

    /** A flat column-major {@code P x P} matrix as a {@link DMatrix}. */
    static DMatrix square(double[] flat) {
        DMatrix matrix = new DMatrix(P, P);
        for (int col = 0; col < P; ++col) {
            for (int row = 0; row < P; ++row) {
                matrix.setUnsafe(row, col, flat[col * P + row]);
            }
        }
        return matrix;
    }

    /** Solves {@code L z = v} for a lower triangular {@code L}. */
    static double[] forwardSolve(DMatrix lower, double[] v) {
        double[] z = new double[v.length];
        for (int i = 0; i < v.length; ++i) {
            double sum = v[i];
            for (int j = 0; j < i; ++j) {
                sum -= lower.getUnsafe(i, j) * z[j];
            }
            z[i] = sum / lower.getUnsafe(i, i);
        }
        return z;
    }

    // ------------------------------------------------------- 7. cross validation

    /** Every wine classified exactly once, by a model that never saw it. */
    static final class CrossValidation {
        final double[] logisticHits;
        final double[] ldaHits;
        final double[] qdaHits;
        final double logistic;
        final double lda;
        final double qda;
        final double pooledCondition;
        final double[] classCondition;

        CrossValidation(double[] logisticHits, double[] ldaHits, double[] qdaHits, double logistic, double lda,
                double qda, double pooledCondition, double[] classCondition) {
            this.logisticHits = logisticHits;
            this.ldaHits = ldaHits;
            this.qdaHits = qdaHits;
            this.logistic = logistic;
            this.lda = lda;
            this.qda = qda;
            this.pooledCondition = pooledCondition;
            this.classCondition = classCondition;
        }
    }

    /** The fold every wine belongs to, from a seeded permutation. */
    static int[] foldOf(int n, int folds, long seed) {
        int[] permutation = new int[n];
        for (int i = 0; i < n; ++i) {
            permutation[i] = i;
        }
        long lcg = seed;
        for (int i = n - 1; i > 0; --i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            int j = (int) (((lcg >>> 11) * 0x1.0p-53) * (i + 1));
            int swap = permutation[i];
            permutation[i] = permutation[j];
            permutation[j] = swap;
        }
        int[] fold = new int[n];
        for (int i = 0; i < n; ++i) {
            fold[permutation[i]] = i % folds;
        }
        return fold;
    }

    static CrossValidation crossValidate(double l2, long seed) {
        double[][] x = standardized();
        int[] y = Datasets.labels();
        int[] fold = foldOf(x.length, FOLDS, seed);

        double[] logisticHits = new double[x.length];
        double[] ldaHits = new double[x.length];
        double[] qdaHits = new double[x.length];
        double pooledCondition = 0.0;
        double[] classCondition = new double[K];

        for (int f = 0; f < FOLDS; ++f) {
            int trainSize = 0;
            for (int i = 0; i < x.length; ++i) {
                if (fold[i] != f) {
                    trainSize++;
                }
            }
            double[][] xt = new double[trainSize][];
            int[] yt = new int[trainSize];
            int at = 0;
            for (int i = 0; i < x.length; ++i) {
                if (fold[i] != f) {
                    xt[at] = x[i];
                    yt[at] = y[i];
                    at++;
                }
            }

            MultinomialLogistic model = new MultinomialLogistic(xt, yt, l2);
            new LimitedMemoryBFGS(model).optimize();
            Gaussian linear = gaussian(xt, yt, true);
            Gaussian quadratic = gaussian(xt, yt, false);
            if (f == 0) {
                pooledCondition = linear.condition[0];
                classCondition = quadratic.condition.clone();
            }

            double[] probability = new double[K];
            for (int i = 0; i < x.length; ++i) {
                if (fold[i] == f) {
                    logisticHits[i] = hit(model.predict(x[i], probability), y[i]);
                    ldaHits[i] = hit(linear.predict(x[i], probability), y[i]);
                    qdaHits[i] = hit(quadratic.predict(x[i], probability), y[i]);
                }
            }
        }
        return new CrossValidation(logisticHits, ldaHits, qdaHits, MEAN.apply(logisticHits), MEAN.apply(ldaHits),
                MEAN.apply(qdaHits), pooledCondition, classCondition);
    }

    /** {@code 1.0} if the most probable class is the observed one, else {@code 0.0}. */
    static double hit(double[] probability, int observed) {
        MultiClassAccuracy measure = new MultiClassAccuracy();
        double[] one = new double[K];
        one[observed] = 1.0;
        measure.compare(K, probability, 0, one);
        return measure.getAverage();
    }

    // ------------------------------------------------------------ 8. interval

    /** An interval around an accuracy, from the per-wine indicators it averages. */
    static final class Interval {
        final double accuracy;
        final double bootstrapMean;
        final double standardError;
        final double[] percentile;
        final double[] bca;

        Interval(double accuracy, double bootstrapMean, double standardError, double[] percentile, double[] bca) {
            this.accuracy = accuracy;
            this.bootstrapMean = bootstrapMean;
            this.standardError = standardError;
            this.percentile = percentile;
            this.bca = bca;
        }
    }

    static Interval interval(double[] hits, double level) {
        Bootstrap bootstrap = new Bootstrap(hits, MEAN, BOOT_ITERATIONS, BOOT_SEED);
        return new Interval(MEAN.apply(hits), bootstrap.getMean(), bootstrap.getStdDev(),
                bootstrap.getConfidenceInterval(level), bootstrap.getConfidenceIntervalBCa(level));
    }

    // -------------------------------------------------------------- output

    /**
     * One row of the table in section 4. {@code decimals} is 1 for the
     * unpenalized runs, whose coefficient norm is round-off past that digit,
     * and 3 for the penalized ones, where it is a result.
     */
    private static void row(String label, Fit fit, int decimals) {
        System.out.println(String.format(L, "  %-26s %5d %6b  %-18s %10.1e %9." + decimals + "f", label,
                Integer.valueOf(fit.iterations), Boolean.valueOf(fit.converged), fit.termination,
                Double.valueOf(fit.gradientNorm), Double.valueOf(fit.norm)));
    }

    /** One cross validation result, as a share and as a count of wines. */
    private static void scored(String label, double accuracy, double[] hits) {
        System.out.println(String.format(L, "    %-34s %.4f   (%d of %d wrong)", label, Double.valueOf(accuracy),
                Integer.valueOf(missed(hits)), Integer.valueOf(hits.length)));
    }

    /** How many of the per-wine indicators are zero. */
    static int missed(double[] hits) {
        int wrong = 0;
        for (int i = 0; i < hits.length; ++i) {
            if (hits[i] == 0.0) {
                wrong++;
            }
        }
        return wrong;
    }

    /** An interval as one printable field. */
    private static String bracket(double[] interval) {
        return String.format(L, "[%.4f, %.4f]", Double.valueOf(interval[0]), Double.valueOf(interval[1]));
    }

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    /** Scatter of the first two components, one glyph per cultivar. */
    private static void scatter(double[][] scores, int[] y, int width, int height) {
        double xLo = Double.POSITIVE_INFINITY, xHi = Double.NEGATIVE_INFINITY;
        double yLo = Double.POSITIVE_INFINITY, yHi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < scores.length; ++i) {
            xLo = Math.min(xLo, scores[i][0]);
            xHi = Math.max(xHi, scores[i][0]);
            yLo = Math.min(yLo, scores[i][1]);
            yHi = Math.max(yHi, scores[i][1]);
        }
        char[] glyph = { 'a', 'b', 'c' };
        char[][] grid = new char[height][width];
        for (int r = 0; r < height; ++r) {
            Arrays.fill(grid[r], ' ');
        }
        for (int i = 0; i < scores.length; ++i) {
            int col = (int) ((scores[i][0] - xLo) / (xHi - xLo) * (width - 1));
            int row = (int) ((scores[i][1] - yLo) / (yHi - yLo) * (height - 1));
            char mark = glyph[y[i]];
            grid[height - 1 - row][col] = grid[height - 1 - row][col] == ' ' || grid[height - 1 - row][col] == mark
                    ? mark
                    : '+';
        }
        for (int r = 0; r < height; ++r) {
            System.out.println(String.format(L, "%8.2f |", Double.valueOf(yHi - (yHi - yLo) * r / (height - 1)))
                    + new String(grid[r]));
        }
        StringBuilder axis = new StringBuilder("         +");
        for (int c = 0; c < width; ++c) {
            axis.append('-');
        }
        System.out.println(axis);
        System.out.println(String.format(L, "          %-" + (width - 8) + ".2f%8.2f", Double.valueOf(xLo),
                Double.valueOf(xHi)));
        System.out.println("          a, b, c are the three cultivars, + is more than one wine in a cell");
    }

    public static void main(String[] args) {
        System.out.println("Wine recognition data, " + Datasets.size() + " wines, " + P + " measurements, " + K
                + " cultivars");
        System.out.println("Source: UCI Machine Learning Repository, retrieved 2026-08-21, CC BY 4.0");
        System.out.println("Seven steps; section 8 states what they established, with the numbers.");

        rule("1. the table  (math.probe)");
        Description d = describe();
        System.out.println(String.format(L, "  %d wines: %d of the first cultivar, %d of the second, %d of the third",
                Integer.valueOf(Datasets.size()), Integer.valueOf(d.counts[0]), Integer.valueOf(d.counts[1]),
                Integer.valueOf(d.counts[2])));
        System.out.println();
        System.out.println(String.format(L, "  %-22s %10s %10s %10s %10s", "measurement", "min", "max", "mean", "sd"));
        for (int j = 0; j < P; ++j) {
            System.out.println(String.format(L, "  %-22s %10.3f %10.3f %10.3f %10.3f", Datasets.featureName(j),
                    Double.valueOf(d.min[j]), Double.valueOf(d.max[j]), Double.valueOf(d.mean[j]),
                    Double.valueOf(d.sd[j])));
        }
        System.out.println();
        System.out.println(String.format(L,
                "  the columns are %.0f times apart in size, from %s", Double.valueOf(d.scaleSpread()),
                Datasets.featureName(d.smallest)));
        System.out.println(String.format(L, "  to %s -- which is the whole of section 2.",
                Datasets.featureName(d.largest)));

        rule("2. principal components of the raw table  (math.linalg)");
        Reduction raw = reduce(Datasets.features());
        System.out.println(String.format(L, "  PC1 explains %.4f of the variance, PC2 %.4f",
                Double.valueOf(raw.explainedRatio[0]), Double.valueOf(raw.explainedRatio[1])));
        System.out.println(String.format(L, "  and PC1 is %s: loading %.6f against %.6f for everything else",
                Datasets.featureName(raw.dominant), Double.valueOf(raw.dominantLoading),
                Double.valueOf(raw.nextLoading)));
        System.out.println("  That is not a finding about wine, it is a finding about units.");

        rule("3. standardize first, then three routes to the same components  (math.linalg)");
        double[][] x = standardized();
        Reduction standard = reduce(x);
        System.out.println(String.format(L, "  PC1 now explains %.4f, PC2 %.4f, together %.4f",
                Double.valueOf(standard.explainedRatio[0]), Double.valueOf(standard.explainedRatio[1]),
                Double.valueOf(standard.explainedRatio[0] + standard.explainedRatio[1])));
        Agreement agreement = pcaRoutes();
        System.out.println();
        System.out.println(String.format(L, "  %-34s %16s %16s", "route", "variance PC1", "variance PC2"));
        System.out.println(String.format(L, "  %-34s %16.12f %16.12f", "CovariancePCA on the covariance",
                Double.valueOf(agreement.covarianceVariance[0]), Double.valueOf(agreement.covarianceVariance[1])));
        System.out.println(String.format(L, "  %-34s %16.12f %16.12f", "JacobiPCA on the data matrix",
                Double.valueOf(agreement.jacobiVariance[0]), Double.valueOf(agreement.jacobiVariance[1])));
        System.out.println(String.format(L, "  %-34s %16.12f %16.12f", "TruncatedPCA, randomized",
                Double.valueOf(agreement.truncatedVariance[0]), Double.valueOf(agreement.truncatedVariance[1])));
        System.out.println();
        System.out.println(String.format(L, "  the two exact routes agree to %.1e in the components, %.1e in variance",
                Double.valueOf(agreement.componentCovJacobi), Double.valueOf(agreement.varianceCovJacobi)));
        System.out.println(String.format(L, "  the randomized one to %.1e, after %d subspace iterations",
                Double.valueOf(agreement.componentCovTruncated), Integer.valueOf(agreement.truncatedIterations)));
        System.out.println(String.format(L, "  squared singular values over n-1 match the eigenvalues to %.1e",
                Double.valueOf(agreement.eigenvalueAgreement)));
        System.out.println();
        scatter(standard.scores, Datasets.labels(), 68, 20);
        System.out.println();
        Separability separable = separability();
        System.out.println("  The cultivars overlap in that picture, and the picture is 2 of 13 axes,");
        System.out.println(String.format(L,
                "  %.1f percent of the variance. Fitted on those two components alone, a",
                Double.valueOf(100.0 * separable.plottedVariance)));
        System.out.println(String.format(L, "  classifier gets %d of %d wines wrong. Fitted on all thirteen it gets",
                Integer.valueOf(separable.wrongInThePlane), Integer.valueOf(Datasets.size())));
        System.out.println(String.format(L,
                "  %d of %d wrong, with a probability margin of %.3f: the cultivars are",
                Integer.valueOf(separable.wrongInFullSpace), Integer.valueOf(Datasets.size()),
                Double.valueOf(separable.smallestMargin)));
        System.out.println("  linearly separable, and the overlap is the projection's, not the data's.");

        rule("4. logistic regression against Optimizable  (math.optim + math.dl)");
        System.out.println("  Wine is linearly separable once standardized, so the unpenalized maximum");
        System.out.println("  likelihood does not exist: the coefficients run off along a ray and the");
        System.out.println("  likelihood creeps up to zero. Watch what the optimizer says about that.");
        System.out.println();
        Fit loose = fit(0.0);
        Fit tight = fitTightened(0.0);
        Fit penalized = fit(L2);
        Fit penalizedTight = fitTightened(L2);
        System.out.println(String.format(L, "  %-26s %5s %6s  %-18s %10s %9s", "run", "iter", "conv", "termination",
                "|gradient|", "||w||"));
        row("unpenalized, defaults", loose, 1);
        row("unpenalized, tightened", tight, 1);
        row("L2 = 1, defaults", penalized, 3);
        row("L2 = 1, tightened", penalizedTight, 3);
        System.out.println();
        System.out.println("  All four report convergence, and the two unpenalized runs report the reason");
        System.out.println("  that sounds most convincing: the gradient was small. It was -- and there is");
        System.out.println("  still no maximum. A supremum that is never attained is approached along a");
        System.out.println("  ray on which the gradient vanishes too, so no single run can tell.");
        System.out.println();
        System.out.println("  What tells is tightening the tolerances and looking at the parameters:");
        System.out.println(String.format(L, "  unpenalized they grow from %.1f to %.1f, penalized they stay at %.3f",
                Double.valueOf(loose.norm), Double.valueOf(tight.norm), Double.valueOf(penalized.norm)));
        System.out.println(String.format(L, "  against %.3f -- a fit that exists does not move when the rule does.",
                Double.valueOf(penalizedTight.norm)));
        System.out.println("  The unpenalized size is printed to one decimal on purpose: past that the");
        System.out.println("  digits differ between the scalar and the vectorized build, which is the same");
        System.out.println("  fact once more -- a point that does not exist cannot be located exactly.");
        System.out.println();
        System.out.println(String.format(L, "  the penalized fit: loglik/wine = %.6f, accuracy on its own data %.4f",
                Double.valueOf(penalized.logLikelihoodPerWine), Double.valueOf(penalized.accuracy)));
        System.out.println("  Now the objective is strictly concave and the answer is a point, not a ray.");

        rule("5. which measurements survive an L1 penalty  (math.optim + math.linalg)");
        Selection selection = selectFeatures(2);
        System.out.println(String.format(L, "  %-10s %12s %12s %9s %9s  %-16s %10s", "l1 weight", "exact zeros",
                "|w| < 1e-6", "||w||", "accuracy", "termination", "|pseudo g|"));
        for (int t = 0; t < selection.l1Weights.length; ++t) {
            System.out.println(String.format(L, "  %-10.1f %12d %12d %9.3f %9.4f  %-16s %10.1e",
                    Double.valueOf(selection.l1Weights[t]), Integer.valueOf(selection.owlqnExactZeros[t]),
                    Integer.valueOf(selection.owlqnNearZeros[t]), Double.valueOf(selection.owlqnNorm[t]),
                    Double.valueOf(selection.owlqnAccuracy[t]), selection.owlqnTermination[t],
                    Double.valueOf(selection.owlqnPseudoGradient[t])));
        }
        System.out.println(String.format(L,
                "  At every weight the two counts are the same: OWL-QN sets %d of %d coefficients to",
                Integer.valueOf(selection.owlqnExactZeros[2]), Integer.valueOf(K * P)));
        System.out.println("  exactly zero at weight 5, rather than merely close to it. A step that would");
        System.out.println("  change a coefficient's sign is clamped to the orthant boundary instead.");
        System.out.println();
        System.out.println(String.format(L, "  coefficients at l1 weight %.1f, and the lasso at its 1-SE penalty:",
                Double.valueOf(selection.l1Weights[2])));
        System.out.println(String.format(L, "  %-22s %24s  %24s", "measurement", "OWL-QN, per cultivar",
                "lasso, one against rest"));
        for (int j = 0; j < P; ++j) {
            StringBuilder line = new StringBuilder(String.format(L, "  %-22s", Datasets.featureName(j)));
            for (int k = 0; k < K; ++k) {
                line.append(String.format(L, " %7.3f", Double.valueOf(selection.owlqnCoefficients[k][j])));
            }
            line.append("  ");
            for (int k = 0; k < K; ++k) {
                line.append(String.format(L, " %7.3f", Double.valueOf(selection.lassoCoefficients[k][j])));
            }
            System.out.println(line);
        }
        System.out.println(String.format(L, "  exact zeros: OWL-QN %d of %d, lasso %d of %d",
                Integer.valueOf(selection.owlqnExactZeros[2]), Integer.valueOf(K * P),
                Integer.valueOf(selection.lassoExactZeros), Integer.valueOf(K * P)));
        System.out.println("  Different objectives -- a likelihood against three squared errors -- so the");
        System.out.println("  numbers differ, but both engines answer the question in exact zeros.");
        System.out.println();
        System.out.println("  Every one of those runs stopped on the value rule, with the pseudo gradient");
        System.out.println(String.format(L,
                "  still at %.2f -- not at the L1 optimality condition, which is what vanishes",
                Double.valueOf(selection.owlqnPseudoGradient[2])));
        System.out.println(String.format(L, "  there. Tightening the tolerances takes %d iterations instead of %d and",
                Integer.valueOf(selection.tightenedIterations), Integer.valueOf(selection.owlqnIterations[2])));
        System.out.println(String.format(L, "  reaches %.1e, and it moves the answer:",
                Double.valueOf(selection.tightenedPseudoGradient)));
        System.out.println(String.format(L,
                "  %d of %d support decisions change and the largest coefficient moves by %.3f.",
                Integer.valueOf(selection.supportChanges), Integer.valueOf(K * P),
                Double.valueOf(selection.largestCoefficientChange)));
        System.out.println("  The selection is worth reading; the third decimal of a coefficient is not.");

        rule("6. a second opinion, fitted by moments  (math.probe + math.linalg)");
        CrossValidation cv = crossValidate(L2, FOLD_SEED);
        System.out.println(String.format(L, "  pooled covariance of the first training split: condition %.2f",
                Double.valueOf(cv.pooledCondition)));
        for (int k = 0; k < K; ++k) {
            System.out.println(String.format(L, "  covariance of cultivar %d alone:                condition %.2f",
                    Integer.valueOf(k + 1), Double.valueOf(cv.classCondition[k])));
        }
        System.out.println("  All four factorize, so the quadratic discriminant is admissible here and");
        System.out.println("  is reported beside the linear one rather than dropped.");
        System.out.println();
        System.out.println(String.format(L, "  %d-fold cross validation, every wine classified once:",
                Integer.valueOf(FOLDS)));
        scored(String.format(L, "logistic regression, L2 = %.1f", Double.valueOf(L2)), cv.logistic, cv.logisticHits);
        scored("linear discriminant", cv.lda, cv.ldaHits);
        scored("quadratic discriminant", cv.qda, cv.qdaHits);
        System.out.println("  Two entirely different routes -- a numerical search against a closed form of");
        System.out.println("  means and covariances -- landing within a wine or two of each other.");

        rule("7. how sure is that number  (math.probe)");
        Interval logistic = interval(cv.logisticHits, 0.95);
        Interval lda = interval(cv.ldaHits, 0.95);
        System.out.println(String.format(L, "  %d bootstrap replications of the per-wine indicators, seed %d",
                Integer.valueOf(BOOT_ITERATIONS), Long.valueOf(BOOT_SEED)));
        System.out.println();
        System.out.println(String.format(L, "  %-22s %9s %10s %20s %20s", "classifier", "accuracy", "std error",
                "95% percentile", "95% BCa"));
        System.out.println(String.format(L, "  %-22s %9.4f %10.4f %20s %20s", "logistic",
                Double.valueOf(logistic.accuracy), Double.valueOf(logistic.standardError),
                bracket(logistic.percentile), bracket(logistic.bca)));
        System.out.println(String.format(L, "  %-22s %9.4f %10.4f %20s %20s", "linear discriminant",
                Double.valueOf(lda.accuracy), Double.valueOf(lda.standardError), bracket(lda.percentile),
                bracket(lda.bca)));
        System.out.println();
        System.out.println("  The intervals overlap almost entirely: on 178 wines the two classifiers are");
        System.out.println("  not distinguishable, and an accuracy quoted without one of these is a number");
        System.out.println("  with three digits of which one is real.");

        rule("8. what this run established");
        System.out.println(String.format(L,
                "  1. Standardizing is not optional here. On the raw table PC1 carries %.4f of",
                Double.valueOf(raw.explainedRatio[0])));
        System.out.println(String.format(L,
                "     the variance and is %s alone, loading %.4f against %.4f for every",
                Datasets.featureName(raw.dominant), Double.valueOf(raw.dominantLoading),
                Double.valueOf(raw.nextLoading)));
        System.out.println(String.format(L, "     other column. Standardized, PC1 carries %.4f.",
                Double.valueOf(standard.explainedRatio[0])));
        System.out.println(String.format(L,
                "  2. The three reduction routes are the same reduction: exact against exact to"));
        System.out.println(String.format(L, "     %.1e, randomized to %.1e. The components are the data's, not one",
                Double.valueOf(agreement.componentCovJacobi), Double.valueOf(agreement.componentCovTruncated)));
        System.out.println("     algorithm's.");
        System.out.println(String.format(L,
                "  3. The cultivars are linearly separable in 13 dimensions -- %d of %d wrong",
                Integer.valueOf(separable.wrongInFullSpace), Integer.valueOf(Datasets.size())));
        System.out.println(String.format(L,
                "     in sample, margin %.3f -- which is exactly why the unpenalized likelihood",
                Double.valueOf(separable.smallestMargin)));
        System.out.println("     has no maximum. In the two plotted components they are not separable.");
        System.out.println(String.format(L,
                "  4. Both unpenalized runs report convergence on the gradient rule, at %.1e",
                Double.valueOf(tight.gradientNorm)));
        System.out.println(String.format(L,
                "     and coefficients of %.1f. Tightening a tolerance is what exposes it: the",
                Double.valueOf(tight.norm)));
        System.out.println(String.format(L, "     unpenalized norm grows %.1f-fold, the penalized one moves %.1f%%.",
                Double.valueOf(tight.norm / loose.norm),
                Double.valueOf(100.0 * Math.abs(penalizedTight.norm - penalized.norm) / penalized.norm)));
        System.out.println(String.format(L,
                "  5. Both L1 engines answer in exact zeros -- OWL-QN %d of %d, the lasso %d --",
                Integer.valueOf(selection.owlqnExactZeros[2]), Integer.valueOf(K * P),
                Integer.valueOf(selection.lassoExactZeros)));
        System.out.println(String.format(L,
                "     but OWL-QN stops on the value rule at a pseudo gradient of %.2f, so its",
                Double.valueOf(selection.owlqnPseudoGradient[2])));
        System.out.println("     selection is trustworthy and its third decimal is not.");
        System.out.println(String.format(L,
                "  6. Out of sample %d of %d wines are misclassified by the logistic model and",
                Integer.valueOf(missed(cv.logisticHits)), Integer.valueOf(Datasets.size())));
        System.out.println(String.format(L,
                "     %d by the quadratic discriminant, and the bootstrap interval [%.4f, %.4f]",
                Integer.valueOf(missed(cv.qdaHits)), Double.valueOf(logistic.bca[0]),
                Double.valueOf(logistic.bca[1])));
        System.out.println("     is wide enough that the three classifiers cannot be told apart.");
        System.out.println();
        System.out.println("  None of that is a fact about wine. Items 3, 4 and 5 are facts about this");
        System.out.println("  library's optimizer, and they are the reason the demo exists.");
    }

    private WineDemo() {
        throw new AssertionError();
    }
}
