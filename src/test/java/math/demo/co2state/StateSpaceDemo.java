package math.demo.co2state;

import java.util.Locale;

import math.demo.maunaloa.Datasets;
import math.fun.DMultiFunction;
import math.fun.DMultiFunctionEval;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.optim.NelderMead;
import math.probe.ACF;
import math.ts.KalmanFilter;
import math.ts.LinearGaussianModel;
import math.ts.RtsSmoother;
import math.ts.TimeVaryingModel;

/**
 * A worked example: the same sixty-eight years of Mauna Loa CO2 that
 * {@code math.demo.maunaloa.MaunaLoaDemo} fits with a curve, fitted instead
 * with a state space model.
 * <p>
 * The two demos share a dataset and answer the same closing question, and the
 * reason for having both is a defect the first one diagnoses about itself.
 * Its step 5 finds that the residuals of a deterministic trend plus harmonics
 * are strongly autocorrelated -- {@code acf(1) = 0.906} -- and concludes that
 * the standard errors of its step 3 are therefore too small. That is not a
 * flaw in the fit; it is a statement that the model is wrong in a particular
 * way. A quadratic cannot bend when the atmosphere bends, so everything the
 * curve fails to follow is left in the residuals, and a residual that
 * remembers the previous residual is what a misspecified mean looks like.
 * <p>
 * A structural time series model puts that memory where it belongs. The level
 * and the slope are allowed to wander -- they are states, not coefficients --
 * so what is left over is only the measurement error, and if the model is
 * right the innovations are white. Measured here: {@code acf(1) = 0.030}
 * against the deterministic fit's {@code 0.906}, and 5 of 24 lags outside the
 * two-sigma band against 24 of 24.
 * <p>
 * <b>The page states that question and answers it first</b>, which it did not
 * always: until 2026-08-30 the paragraph above lived only here, the page opened
 * on a procedure -- how many harmonics, and what are the four variances -- and
 * the answer went past at position three in a dozen lines. Sections 4 to 7 also
 * say what they have in common, because a reader is owed the reason four
 * capabilities are in one demo: <b>a fitted curve can do none of them</b>. It
 * has nowhere to put a missing month, no scale on which a surprise is large or
 * small, no way of being told that one observation was better measured than
 * another, and nothing to say when the months are not evenly spaced.
 * <p>
 * Eight steps, across {@code math.ts}, {@code math.optim}, {@code math.probe}
 * and {@code math.linalg}:
 * <ol>
 * <li>test the innovations for whiteness against the deterministic residuals,
 * which is the question the demo exists to answer and so is answered first,</li>
 * <li>fit the four variances by maximum likelihood, let the likelihood choose
 * how many seasonal harmonics the series wants, and -- since three of the four
 * are three more than the page needs -- turn the one that matters from end to
 * end, from a level frozen into a coefficient to one that follows every
 * reading, and watch the innovations,</li>
 * <li>decompose the series into level, slope and season, and check the level
 * against the deseasonalized series NOAA publishes beside the data,</li>
 * <li>blank the two months NOAA interpolated and see what the smoother puts
 * back,</li>
 * <li>look for the 2022 eruption in the anomaly score,</li>
 * <li>replace the fitted observation variance with the uncertainty NOAA
 * publishes month by month, which is a {@code TimeVaryingModel},</li>
 * <li>throw away seven months in ten at random and recover the level from an
 * irregular schedule, where the transition and the process noise depend on the
 * gap,</li>
 * <li>forecast the 450 ppm crossing with an interval, and compare it with the
 * point estimate the deterministic demo gives.</li>
 * </ol>
 * <p>
 * <b>Three of the eight are worth the trip on their own.</b> Step 5 finds
 * nothing where something was expected and something where nothing was: the
 * eight months the observatory spent at Maunakea, 21 miles away, leave a mean
 * innovation of {@code +0.18} against a standard error of {@code 0.35}, which
 * is no offset at all, while the largest anomaly in the whole record falls in
 * 2016-04 -- the El Nino that produced the largest annual CO2 growth ever
 * measured. Step 6 buys 23 nats of log likelihood from NOAA's own published
 * uncertainties, at no extra fitted parameter. And step 7 discards 578 of the
 * 821 months, with gaps of up to 18 months, and still recovers the level to
 * 0.09 ppm.
 * <p>
 * <b>One number is deliberately missing from the page.</b> The demo does not
 * print how many iterations the fit took, because that is the one quantity on
 * it that is not the same on JDK 8 and JDK 25 -- measured, 439 against 432. The
 * filter is not the reason: {@link KalmanFilter#logLikelihood} at a fixed set
 * of variances returns the identical {@code double} on both, down to the last
 * bit. What differs is {@code Math.exp} and {@code Math.cos} used to build the
 * model, which are permitted an ulp of error and take different intrinsics on
 * the two runtimes; one function value out by an ulp flips a comparison in the
 * simplex and the two searches take different routes to the same place. They
 * arrive within {@code 6e-12} of each other, so every figure the demo prints is
 * identical -- but an iteration count is a property of the route and not of the
 * answer, and printing it would make a reproducible page unreproducible. It is
 * still on {@link Fit}, where a test can look at it.
 * <p>
 * Run {@link #main(String[])} to see it. Every step is a package-private method
 * returning a small value type, so {@code StateSpaceDemoTest} can assert what
 * this prints.
 *
 * @see math.demo.maunaloa.Datasets
 */
public final class StateSpaceDemo {

    private static final Locale L = Locale.ROOT;
    private static final double TWO_PI = 2.0 * Math.PI;

    /** Months in the seasonal period. */
    private static final int PERIOD = 12;

    /**
     * How many seasonal harmonics the fitted model uses. Step 1 is what picks
     * it; measured, the log likelihood runs -920.5, -263.3, -255.7, -252.4 for
     * one through four, and the lags outside the band run 24, 8, 5, 2 of 24. A
     * fourth harmonic still buys something and costs twice the fitting time, so
     * three is where this stops.
     */
    static final int HARMONICS = 3;

    /** What the deterministic demo uses: harmonics at 12 and 6 months. */
    static final int DETERMINISTIC_HARMONICS = 2;

    /** Number of terms in the deterministic design of step 3. */
    private static final int TERMS = 7;

    /** Where every fit starts, in log variance. */
    private static final double[] START = { Math.log(0.01), Math.log(1.0e-5), Math.log(1.0e-4),
            Math.log(0.1) };

    /** Seed for the thinning of step 7, so that the printed schedule is fixed. */
    static final long THINNING_SEED = 20260827L;

    /** Roughly what fraction of the months step 7 keeps. */
    static final double KEEP_FRACTION = 0.3;

    /** The first month measured at Maunakea after the 2022 eruption. */
    static final String MOVE_FIRST = "2022-12";

    /** The last month measured at Maunakea. */
    static final String MOVE_LAST = "2023-07";

    // -------------------------------------------------------------- the model

    /** Level, slope, and a cosine and a sine for each harmonic. */
    static int dimension(int harmonics) {
        return 2 + 2 * harmonics;
    }

    /**
     * The transition over a gap of {@code dt} months: the level takes up
     * {@code dt} times the slope, and each seasonal pair rotates by its own
     * angle times {@code dt}. At {@code dt == 1} this is the ordinary monthly
     * model.
     */
    static DMatrix transition(int harmonics, double dt) {
        int d = dimension(harmonics);
        DMatrix f = new DMatrix(d, d);
        f.set(0, 0, 1.0);
        f.set(0, 1, dt);
        f.set(1, 1, 1.0);
        for (int j = 0; j < harmonics; ++j) {
            double angle = TWO_PI * (j + 1) * dt / PERIOD;
            int at = 2 + 2 * j;
            f.set(at, at, Math.cos(angle));
            f.set(at, at + 1, Math.sin(angle));
            f.set(at + 1, at, -Math.sin(angle));
            f.set(at + 1, at + 1, Math.cos(angle));
        }
        return f;
    }

    /** What is observed: the level plus the cosine of every harmonic. */
    static DMatrix observation(int harmonics) {
        DMatrix h = new DMatrix(1, dimension(harmonics));
        h.set(0, 0, 1.0);
        for (int j = 0; j < harmonics; ++j) {
            h.set(0, 2 + 2 * j, 1.0);
        }
        return h;
    }

    /**
     * The process noise over a gap of {@code dt} months, scaled linearly in the
     * gap.
     * <p>
     * That is the simple scaling and not the exact one. A level that also
     * carries the slope's wander over the gap picks up a term in
     * {@code dt^3 / 3} and an off-diagonal in {@code dt^2 / 2}, which is what a
     * navigation filter would use. The linear form is here because it reduces
     * <em>exactly</em> to the monthly model at {@code dt == 1}, which is what
     * makes step 7 a comparison rather than two unrelated fits; how much the
     * difference costs is what step 7 measures.
     */
    static DMatrix processNoise(int harmonics, double[] logVariance, double dt) {
        int d = dimension(harmonics);
        DMatrix q = new DMatrix(d, d);
        q.set(0, 0, Math.exp(logVariance[0]) * dt);
        q.set(1, 1, Math.exp(logVariance[1]) * dt);
        double seasonal = Math.exp(logVariance[2]) * dt;
        for (int i = 2; i < d; ++i) {
            q.set(i, i, seasonal);
        }
        return q;
    }

    static DMatrix observationNoise(double[] logVariance) {
        DMatrix r = new DMatrix(1, 1);
        r.set(0, 0, Math.exp(logVariance[3]));
        return r;
    }

    /** The prior: the first observation for the level, zero for everything else. */
    static double[] priorMean(int harmonics) {
        double[] m0 = new double[dimension(harmonics)];
        m0[0] = Datasets.co2()[0];
        return m0;
    }

    /**
     * A deliberately vague prior. It has to be proper -- {@code P0} is what the
     * filter factors -- but 100 ppm of standard deviation on a level near 315
     * says nothing that 821 observations cannot overrule.
     */
    static DMatrix priorCovariance(int harmonics) {
        int d = dimension(harmonics);
        DMatrix p0 = new DMatrix(d, d);
        p0.set(0, 0, 1.0e4);
        p0.set(1, 1, 1.0);
        for (int i = 2; i < d; ++i) {
            p0.set(i, i, 10.0);
        }
        return p0;
    }

    static LinearGaussianModel model(int harmonics, double[] logVariance) {
        return new LinearGaussianModel(transition(harmonics, 1.0),
                processNoise(harmonics, logVariance, 1.0), observation(harmonics),
                observationNoise(logVariance), priorMean(harmonics), priorCovariance(harmonics));
    }

    static DMatrix column(double[] v) {
        DMatrix out = new DMatrix(v.length, 1);
        for (int i = 0; i < v.length; ++i) {
            out.set(i, 0, v[i]);
        }
        return out;
    }

    // -------------------------------------------- 1. fit it by maximum likelihood

    /** What the maximum likelihood fit of one model returned. */
    static final class Fit {

        /** How many seasonal harmonics this model carries. */
        public final int harmonics;
        /** The four variances, as logarithms, in the order the fit uses them. */
        public final double[] logVariance;
        /** The log likelihood at the optimum. */
        public final double logLikelihood;
        /** How many Nelder-Mead iterations it took. */
        public final int iterations;
        /** Whether the simplex converged rather than running out of budget. */
        public final boolean converged;

        Fit(int harmonics, double[] logVariance, double logLikelihood, int iterations,
                boolean converged) {
            this.harmonics = harmonics;
            this.logVariance = logVariance;
            this.logLikelihood = logLikelihood;
            this.iterations = iterations;
            this.converged = converged;
        }

        /** The standard deviation of one component, on the scale of the data. */
        public double sigma(int component) {
            return Math.exp(0.5 * logVariance[component]);
        }

        /**
         * How much the level is allowed to move relative to how badly it is
         * measured, {@code sigmaLevel^2 / sigmaObservation^2}.
         * <p>
         * The only thing about these four numbers that matters on its own: what
         * the filter does is decided by the ratio, not by the scale. A
         * deterministic curve is the corner where this is zero.
         */
        public double signalToNoise() {
            return Math.exp(logVariance[0] - logVariance[3]);
        }

        /** The model these variances describe. */
        public LinearGaussianModel model() {
            return StateSpaceDemo.model(harmonics, logVariance);
        }
    }

    /**
     * Maximizes {@link KalmanFilter#logLikelihood} over the four variances.
     * <p>
     * This is the prediction error decomposition being used for what it is for.
     * The search runs over the logarithms, because a variance is positive and
     * the likelihood is flat in it over orders of magnitude; the box is a
     * refusal rather than a penalty, since a variance of {@code e^-40} is not a
     * model anyone meant.
     */
    static Fit fit(final int harmonics) {
        final DMatrix y = column(Datasets.co2());
        DMultiFunction negativeLogLikelihood = new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                for (int i = 0; i < x.length; ++i) {
                    if (x[i] < -40.0 || x[i] > 20.0) {
                        return 1.0e30;
                    }
                }
                return -KalmanFilter.logLikelihood(model(harmonics, x), y);
            }
        };
        DMultiFunctionEval found = new NelderMead().minimize(negativeLogLikelihood, START.clone());
        return new Fit(harmonics, found.point, -found.value, found.iterations, found.converged);
    }

    /**
     * How much of one month's surprise the filter takes as a real move of the
     * level, once the gain has settled.
     * <p>
     * The Kalman gain against the level, {@code (P- H') / S} in its first
     * coordinate. This is what the four variances actually decide: at zero the
     * level never responds and the model is a rigid curve, at one it follows
     * every observation and there is no smoothing left.
     */
    static double steadyStateLevelGain(Fit fit) {
        KalmanFilter.Result filtered = KalmanFilter.filter(fit.model(), column(Datasets.co2()));
        int last = filtered.length - 1;
        DMatrix predicted = filtered.predictedCovariance[last];
        DMatrix h = observation(fit.harmonics);
        double ph = 0.0;
        for (int j = 0; j < dimension(fit.harmonics); ++j) {
            ph += predicted.get(0, j) * h.get(0, j);
        }
        return ph / filtered.innovationCovariance[last].get(0, 0);
    }

    /**
     * One row of the sweep in section 2: the same model with the level allowed
     * to move a different amount, and what that does to the innovations.
     */
    static final class Knob {

        /** How the row is labelled on the page. */
        public final String label;
        /** The share of a month's surprise the filter takes as a real move. */
        public final double gain;
        /** First autocorrelation of the standardized innovations. */
        public final double acf1;
        /** How many of the first 24 lags fall outside the two sigma band. */
        public final int outside;
        /** Whether this is the row the likelihood chose. */
        public final boolean fitted;

        Knob(String label, double gain, double acf1, int outside, boolean fitted) {
            this.label = label;
            this.gain = gain;
            this.acf1 = acf1;
            this.outside = outside;
            this.fitted = fitted;
        }
    }

    /**
     * The same model at six settings of the one knob that matters, from a rigid
     * curve to a level that follows every reading.
     * <p>
     * The frozen row is what a deterministic fit <i>is</i> in this vocabulary:
     * the level, the slope and the seasonal shape are all coefficients rather
     * than states, so the model has to explain the whole record with a straight
     * line plus fixed harmonics and everything it cannot follow is left in the
     * residuals. It is a worse curve than {@code MaunaLoaDemo}'s quadratic and
     * its residuals are correspondingly worse.
     */
    static Knob[] knobSweep(Fit fit) {
        double frozen = Math.log(1.0e-24);
        double[] v = fit.logVariance;
        return new Knob[] {
                knob("not at all (a rigid line)", fit.harmonics,
                        new double[] { frozen, frozen, frozen, v[3] }, false),
                knob("a tenth of the fitted", fit.harmonics, scaled(v, 0.1), false),
                knob("half", fit.harmonics, scaled(v, 0.5), false),
                knob("as fitted", fit.harmonics, v.clone(), true),
                knob("twice", fit.harmonics, scaled(v, 2.0), false),
                knob("ten times", fit.harmonics, scaled(v, 10.0), false) };
    }

    /** The fitted variances with the level's standard deviation scaled. */
    private static double[] scaled(double[] logVariance, double factor) {
        double[] out = logVariance.clone();
        out[0] = logVariance[0] + 2.0 * Math.log(factor);
        return out;
    }

    private static Knob knob(String label, int harmonics, double[] logVariance, boolean fitted) {
        Fit f = new Fit(harmonics, logVariance, 0.0, 0, true);
        double[] z = standardizedInnovations(f);
        double[] acf = ACF.acf(z, 24);
        double bound = 1.96 / Math.sqrt(z.length);
        return new Knob(label, steadyStateLevelGain(f), acf[1], countOutside(acf, bound), fitted);
    }

    // ------------------------------------------------------ 2. take it apart

    /** The series split into the parts the model says it is made of. */
    static final class Decomposition {

        /** The smoothed level, month by month. */
        public final double[] level;
        /** The smoothed slope, in ppm per month. */
        public final double[] slope;
        /** The smoothed seasonal component. */
        public final double[] season;
        /** Root mean square difference between the level and NOAA's own series. */
        public final double rmsAgainstNoaa;
        /** The largest such difference. */
        public final double worstAgainstNoaa;
        /** Where that largest difference falls. */
        public final int worstAt;
        /** The peak to trough range of the seasonal component in the last full year. */
        public final double seasonalAmplitude;

        Decomposition(double[] level, double[] slope, double[] season, double rmsAgainstNoaa,
                double worstAgainstNoaa, int worstAt, double seasonalAmplitude) {
            this.level = level;
            this.slope = slope;
            this.season = season;
            this.rmsAgainstNoaa = rmsAgainstNoaa;
            this.worstAgainstNoaa = worstAgainstNoaa;
            this.worstAt = worstAt;
            this.seasonalAmplitude = seasonalAmplitude;
        }
    }

    /**
     * Smooths the series and reads the states off.
     * <p>
     * The check is worth more than the decomposition: NOAA publishes its own
     * deseasonalized series beside the raw one, produced by a method that has
     * nothing to do with this, and the smoothed level has to agree with it.
     */
    static Decomposition decompose(Fit fit) {
        double[] co2 = Datasets.co2();
        int n = co2.length;
        LinearGaussianModel m = fit.model();
        RtsSmoother.Result smoothed = RtsSmoother.smooth(m, column(co2));
        DMatrix h = observation(fit.harmonics);

        double[] level = new double[n];
        double[] slope = new double[n];
        double[] season = new double[n];
        for (int t = 0; t < n; ++t) {
            level[t] = smoothed.mean[t][0];
            slope[t] = smoothed.mean[t][1];
            double s = 0.0;
            for (int j = 2; j < dimension(fit.harmonics); ++j) {
                s += h.get(0, j) * smoothed.mean[t][j];
            }
            season[t] = s;
        }

        double[] noaa = Datasets.deseasonalized();
        double sum = 0.0;
        double worst = 0.0;
        int worstAt = 0;
        for (int t = 0; t < n; ++t) {
            double gap = level[t] - noaa[t];
            sum += gap * gap;
            if (Math.abs(gap) > worst) {
                worst = Math.abs(gap);
                worstAt = t;
            }
        }

        double high = -Double.MAX_VALUE;
        double low = Double.MAX_VALUE;
        for (int t = n - PERIOD; t < n; ++t) {
            high = Math.max(high, season[t]);
            low = Math.min(low, season[t]);
        }
        return new Decomposition(level, slope, season, Math.sqrt(sum / n), worst, worstAt, high - low);
    }

    // ------------------------------------------------- 3. are they white now

    /** The autocorrelation of what each of the two models leaves behind. */
    static final class Whiteness {

        /** Two standard errors under the hypothesis that the series is white. */
        public final double bound;
        /** First autocorrelation of the standardized innovations. */
        public final double structuralAcf1;
        /** Twelfth autocorrelation of the standardized innovations. */
        public final double structuralAcf12;
        /** How many of the first 24 lags leave the band. */
        public final int structuralOutside;
        /** First autocorrelation of the deterministic fit's residuals. */
        public final double deterministicAcf1;
        /** How many of its first 24 lags leave the band. */
        public final int deterministicOutside;
        /** Mean of the standardized innovations, which should be near zero. */
        public final double innovationMean;
        /** Their standard deviation, which should be near one. */
        public final double innovationSd;

        Whiteness(double bound, double structuralAcf1, double structuralAcf12, int structuralOutside,
                double deterministicAcf1, int deterministicOutside, double innovationMean,
                double innovationSd) {
            this.bound = bound;
            this.structuralAcf1 = structuralAcf1;
            this.structuralAcf12 = structuralAcf12;
            this.structuralOutside = structuralOutside;
            this.deterministicAcf1 = deterministicAcf1;
            this.deterministicOutside = deterministicOutside;
            this.innovationMean = innovationMean;
            this.innovationSd = innovationSd;
        }
    }

    /** The standardized innovations of a fit, {@code e / sqrt(S)} at each step. */
    static double[] standardizedInnovations(Fit fit) {
        KalmanFilter.Result filtered = KalmanFilter.filter(fit.model(), column(Datasets.co2()));
        double[] z = new double[filtered.length];
        for (int t = 0; t < z.length; ++t) {
            if (filtered.observedComponents[t] == 0) {
                z[t] = Double.NaN;
            } else {
                z[t] = filtered.innovation[t][0]
                        / Math.sqrt(filtered.innovationCovariance[t].get(0, 0));
            }
        }
        return z;
    }

    /**
     * The residuals of the deterministic model: a quadratic in time plus
     * harmonics at 12 and 6 months, exactly the design
     * {@code MaunaLoaDemo} fits in its step 3.
     */
    static double[] deterministicResiduals() {
        double[] date = Datasets.decimalDate();
        int n = date.length;
        double mean = 0.0;
        for (int i = 0; i < n; ++i) {
            mean += date[i];
        }
        mean /= n;
        DMatrix x = new DMatrix(n, TERMS);
        for (int i = 0; i < n; ++i) {
            double t = date[i] - mean;
            x.set(i, 0, 1.0);
            x.set(i, 1, t);
            x.set(i, 2, t * t);
            x.set(i, 3, Math.sin(TWO_PI * t));
            x.set(i, 4, Math.cos(TWO_PI * t));
            x.set(i, 5, Math.sin(2.0 * TWO_PI * t));
            x.set(i, 6, Math.cos(2.0 * TWO_PI * t));
        }
        LSSummary summary = OLS.estimate(0.05, x, column(Datasets.co2()));
        return summary.getResiduals().toArray();
    }

    static Whiteness whiteness(Fit fit) {
        double[] z = standardizedInnovations(fit);
        int n = z.length;
        double bound = 1.96 / Math.sqrt(n);
        double[] structural = ACF.acf(z, 24);
        int structuralOutside = countOutside(structural, bound);

        double[] residuals = deterministicResiduals();
        double[] deterministic = ACF.acf(residuals, 24);
        int deterministicOutside = countOutside(deterministic, bound);

        double sum = 0.0;
        for (int t = 0; t < n; ++t) {
            sum += z[t];
        }
        double mean = sum / n;
        double variance = 0.0;
        for (int t = 0; t < n; ++t) {
            variance += (z[t] - mean) * (z[t] - mean);
        }
        return new Whiteness(bound, structural[1], structural[12], structuralOutside,
                deterministic[1], deterministicOutside, mean, Math.sqrt(variance / (n - 1)));
    }

    private static int countOutside(double[] acf, double bound) {
        int outside = 0;
        for (int k = 1; k <= 24; ++k) {
            if (Math.abs(acf[k]) > bound) {
                ++outside;
            }
        }
        return outside;
    }

    // -------------------------------- 4. the two months NOAA filled in itself

    /** What the smoother says about months that were never measured. */
    static final class Filled {

        /** The indices NOAA marks as interpolated. */
        public final int[] index;
        /** What NOAA put there. */
        public final double[] noaa;
        /** What the smoother puts there when the month is blanked out. */
        public final double[] smoothed;
        /** The standard deviation the smoother attaches to its own answer. */
        public final double[] sd;

        Filled(int[] index, double[] noaa, double[] smoothed, double[] sd) {
            this.index = index;
            this.noaa = noaa;
            this.smoothed = smoothed;
            this.sd = sd;
        }
    }

    /**
     * Blanks the two months NOAA interpolated and lets the smoother put them
     * back.
     * <p>
     * Nothing about the model changes -- a {@code NaN} is simply a step with no
     * correction -- so this is the whole of what handling missing data costs
     * here.
     */
    static Filled fillTheInterpolatedMonths(Fit fit) {
        double[] co2 = Datasets.co2();
        int[] index = Datasets.noaaInterpolated();
        double[] blanked = co2.clone();
        for (int i = 0; i < index.length; ++i) {
            blanked[index[i]] = Double.NaN;
        }
        LinearGaussianModel m = fit.model();
        RtsSmoother.Result smoothed = RtsSmoother.smooth(m, column(blanked));
        DMatrix h = observation(fit.harmonics);
        int d = dimension(fit.harmonics);

        double[] noaa = new double[index.length];
        double[] guess = new double[index.length];
        double[] sd = new double[index.length];
        for (int i = 0; i < index.length; ++i) {
            int at = index[i];
            noaa[i] = co2[at];
            double mean = 0.0;
            for (int j = 0; j < d; ++j) {
                mean += h.get(0, j) * smoothed.mean[at][j];
            }
            guess[i] = mean;
            double variance = 0.0;
            for (int a = 0; a < d; ++a) {
                for (int b = 0; b < d; ++b) {
                    variance += h.get(0, a) * smoothed.covariance[at].get(a, b) * h.get(0, b);
                }
            }
            sd[i] = Math.sqrt(variance);
        }
        return new Filled(index, noaa, guess, sd);
    }

    // ------------------------------------------------- 5. look for the move

    /** What the anomaly score found, and what it did not. */
    static final class Anomalies {

        /** The largest standardized innovation anywhere in the record. */
        public final double largest;
        /** Which month that was. */
        public final int largestAt;
        /** Mean standardized innovation over the months spent at Maunakea. */
        public final double meanDuringMove;
        /** The largest one over those months. */
        public final double largestDuringMove;
        /** Which month that was. */
        public final int largestDuringMoveAt;
        /** How many months those were. */
        public final int moveMonths;
        /** Mean standardized innovation over the rest of the record. */
        public final double meanElsewhere;
        /** How many months exceed three standard deviations. */
        public final int beyondThree;

        Anomalies(double largest, int largestAt, double meanDuringMove, double largestDuringMove,
                int largestDuringMoveAt, int moveMonths, double meanElsewhere, int beyondThree) {
            this.largest = largest;
            this.largestAt = largestAt;
            this.meanDuringMove = meanDuringMove;
            this.largestDuringMove = largestDuringMove;
            this.largestDuringMoveAt = largestDuringMoveAt;
            this.moveMonths = moveMonths;
            this.meanElsewhere = meanElsewhere;
            this.beyondThree = beyondThree;
        }

        /**
         * The standard error of {@link #meanDuringMove} under the model, which
         * is what says whether it is a step change or a coincidence: the
         * innovations are standard normal, so the mean of {@code moveMonths} of
         * them has standard deviation {@code 1 / sqrt(moveMonths)}.
         */
        public double moveStandardError() {
            return 1.0 / Math.sqrt(moveMonths);
        }
    }

    /** Whether a month was measured at Maunakea rather than at Mauna Loa. */
    static boolean measuredAtMaunakea(int index) {
        String label = Datasets.label(index);
        return label.compareTo(MOVE_FIRST) >= 0 && label.compareTo(MOVE_LAST) <= 0;
    }

    /**
     * Runs the anomaly score over the record.
     * <p>
     * The innovation comes with its own covariance, so {@code e / sqrt(S)} is
     * already scaled by how uncertain the prediction was. Under the model it is
     * standard normal, which makes it a score that can be compared against a
     * threshold rather than a residual that has to be eyeballed.
     */
    static Anomalies anomalies(Fit fit) {
        double[] z = standardizedInnovations(fit);
        int n = z.length;
        double largest = 0.0;
        int largestAt = 0;
        double duringSum = 0.0;
        double duringLargest = 0.0;
        int duringLargestAt = -1;
        int during = 0;
        double elsewhereSum = 0.0;
        int elsewhere = 0;
        int beyondThree = 0;
        for (int t = 0; t < n; ++t) {
            double magnitude = Math.abs(z[t]);
            if (magnitude > largest) {
                largest = magnitude;
                largestAt = t;
            }
            if (magnitude > 3.0) {
                ++beyondThree;
            }
            if (measuredAtMaunakea(t)) {
                duringSum += z[t];
                if (magnitude > duringLargest) {
                    duringLargest = magnitude;
                    duringLargestAt = t;
                }
                ++during;
            } else {
                elsewhereSum += z[t];
                ++elsewhere;
            }
        }
        return new Anomalies(largest, largestAt, duringSum / during, duringLargest, duringLargestAt,
                during, elsewhereSum / elsewhere, beyondThree);
    }

    // -------------------------- 6. the uncertainty NOAA publishes, month by month

    /** What replacing the fitted observation variance with the published one buys. */
    static final class Published {

        /** How many months carry a usable published uncertainty. */
        public final int months;
        /** The smallest, mean and largest published standard deviation. */
        public final double minSd;
        /** See {@link #minSd}. */
        public final double meanSd;
        /** See {@link #minSd}. */
        public final double maxSd;
        /** The single standard deviation the fit of step 1 chose instead. */
        public final double fittedSd;
        /** Log likelihood with one constant observation variance. */
        public final double logLikelihoodConstant;
        /** Log likelihood with the published variance where there is one. */
        public final double logLikelihoodPublished;
        /** The largest the smoothed level moves between the two. */
        public final double worstLevelShift;
        /** Where it moves most. */
        public final int worstAt;

        Published(int months, double minSd, double meanSd, double maxSd, double fittedSd,
                double logLikelihoodConstant, double logLikelihoodPublished, double worstLevelShift,
                int worstAt) {
            this.months = months;
            this.minSd = minSd;
            this.meanSd = meanSd;
            this.maxSd = maxSd;
            this.fittedSd = fittedSd;
            this.logLikelihoodConstant = logLikelihoodConstant;
            this.logLikelihoodPublished = logLikelihoodPublished;
            this.worstLevelShift = worstLevelShift;
            this.worstAt = worstAt;
        }

        /** What the published uncertainties are worth, in nats of log likelihood. */
        public double gain() {
            return logLikelihoodPublished - logLikelihoodConstant;
        }
    }

    /**
     * Uses NOAA's published monthly uncertainty as {@code R(t)}.
     * <p>
     * This is a {@code TimeVaryingModel} for the reason the class is for: the
     * measurement setup genuinely changes. 625 of the 821 months carry an
     * uncertainty and 196 do not -- the Keeling era at Scripps reported none --
     * so those keep the fitted constant. Nothing is refitted, so the gain is a
     * lower bound on what the published numbers are worth.
     */
    static Published publishedUncertainty(Fit fit) {
        double[] co2 = Datasets.co2();
        double[] uncertainty = Datasets.uncertainty();
        int n = co2.length;
        DMatrix y = column(co2);

        double minSd = Double.MAX_VALUE;
        double maxSd = 0.0;
        double sumSd = 0.0;
        int months = 0;
        double fittedVariance = Math.exp(fit.logVariance[3]);
        DMatrix[] r = new DMatrix[n];
        for (int t = 0; t < n; ++t) {
            DMatrix rt = new DMatrix(1, 1);
            if (Datasets.hasUncertainty(t)) {
                double sd = uncertainty[t];
                minSd = Math.min(minSd, sd);
                maxSd = Math.max(maxSd, sd);
                sumSd += sd;
                ++months;
                rt.set(0, 0, sd * sd);
            } else {
                rt.set(0, 0, fittedVariance);
            }
            r[t] = rt;
        }

        int h = fit.harmonics;
        TimeVaryingModel varying = new TimeVaryingModel(
                TimeVaryingModel.repeat(transition(h, 1.0), n),
                TimeVaryingModel.repeat(processNoise(h, fit.logVariance, 1.0), n),
                TimeVaryingModel.repeat(observation(h), n), r, priorMean(h), priorCovariance(h));

        LinearGaussianModel constant = fit.model();
        KalmanFilter.Result constantRun = KalmanFilter.filter(constant, y);
        KalmanFilter.Result varyingRun = KalmanFilter.filter(varying, y);
        RtsSmoother.Result constantSmoothed = RtsSmoother.smooth(constant, constantRun);
        RtsSmoother.Result varyingSmoothed = RtsSmoother.smooth(varying, varyingRun);

        double worst = 0.0;
        int worstAt = 0;
        for (int t = 0; t < n; ++t) {
            double shift = Math.abs(varyingSmoothed.mean[t][0] - constantSmoothed.mean[t][0]);
            if (shift > worst) {
                worst = shift;
                worstAt = t;
            }
        }
        return new Published(months, minSd, sumSd / months, maxSd, fit.sigma(3),
                constantRun.logLikelihood, varyingRun.logLikelihood, worst, worstAt);
    }

    // ------------------------------------------ 7. throw most of the record away

    /** What survives an irregular schedule. */
    static final class Irregular {

        /** How many months were kept. */
        public final int kept;
        /** Out of how many. */
        public final int total;
        /** The largest gap between two kept months. */
        public final int largestGap;
        /** Root mean square difference against the level from the whole record. */
        public final double rmsAgainstFull;
        /** The largest such difference. */
        public final double worstAgainstFull;
        /** Root mean square difference against NOAA's deseasonalized series. */
        public final double rmsAgainstNoaa;

        Irregular(int kept, int total, int largestGap, double rmsAgainstFull, double worstAgainstFull,
                double rmsAgainstNoaa) {
            this.kept = kept;
            this.total = total;
            this.largestGap = largestGap;
            this.rmsAgainstFull = rmsAgainstFull;
            this.worstAgainstFull = worstAgainstFull;
            this.rmsAgainstNoaa = rmsAgainstNoaa;
        }
    }

    /**
     * Keeps roughly three months in ten, at irregular spacing, and recovers the
     * level.
     * <p>
     * Each retained step spans a different number of months, so the transition
     * and the process noise depend on the gap -- this is the case that cannot
     * be written with constant matrices, and the one {@code TimeVaryingModel}
     * exists for. The gaps are not missing data: a missing month is a
     * {@code NaN} and needs none of this. The difference is that here the steps
     * themselves have different lengths.
     */
    static Irregular thinToAnIrregularSchedule(Fit fit) {
        double[] co2 = Datasets.co2();
        int n = co2.length;

        long lcg = THINNING_SEED;
        boolean[] keep = new boolean[n];
        int kept = 0;
        for (int t = 0; t < n; ++t) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double u = ((lcg >>> 11) + 0.5) * 0x1.0p-53;
            keep[t] = t == 0 || t == n - 1 || u < KEEP_FRACTION;
            if (keep[t]) {
                ++kept;
            }
        }
        int[] index = new int[kept];
        int at = 0;
        for (int t = 0; t < n; ++t) {
            if (keep[t]) {
                index[at++] = t;
            }
        }

        int h = fit.harmonics;
        DMatrix thinned = new DMatrix(kept, 1);
        DMatrix[] f = new DMatrix[kept];
        DMatrix[] q = new DMatrix[kept];
        int largestGap = 0;
        for (int i = 0; i < kept; ++i) {
            thinned.set(i, 0, co2[index[i]]);
            // the entry at zero is never read, so any valid matrix will do
            double dt = i == 0 ? 1.0 : index[i] - index[i - 1];
            largestGap = Math.max(largestGap, (int) dt);
            f[i] = transition(h, dt);
            q[i] = processNoise(h, fit.logVariance, dt);
        }
        TimeVaryingModel schedule = new TimeVaryingModel(f, q,
                TimeVaryingModel.repeat(observation(h), kept),
                TimeVaryingModel.repeat(observationNoise(fit.logVariance), kept), priorMean(h),
                priorCovariance(h));

        RtsSmoother.Result thin = RtsSmoother.smooth(schedule, thinned);
        RtsSmoother.Result full = RtsSmoother.smooth(fit.model(), column(co2));
        double[] noaa = Datasets.deseasonalized();

        double sumFull = 0.0;
        double worstFull = 0.0;
        double sumNoaa = 0.0;
        for (int i = 0; i < kept; ++i) {
            double againstFull = thin.mean[i][0] - full.mean[index[i]][0];
            sumFull += againstFull * againstFull;
            worstFull = Math.max(worstFull, Math.abs(againstFull));
            double againstNoaa = thin.mean[i][0] - noaa[index[i]];
            sumNoaa += againstNoaa * againstNoaa;
        }
        return new Irregular(kept, n, largestGap, Math.sqrt(sumFull / kept), worstFull,
                Math.sqrt(sumNoaa / kept));
    }

    // ---------------------------------------------------- 8. when is 450 ppm

    /** When the model expects a given level, and how sure it is. */
    static final class Crossing {

        /** The level asked about. */
        public final double target;
        /** Months past the end of the record for the earliest the band allows. */
        public final int earliest;
        /** Months past the end for the central estimate. */
        public final int central;
        /** Months past the end for the latest the band allows. */
        public final int latest;

        Crossing(double target, int earliest, int central, int latest) {
            this.target = target;
            this.earliest = earliest;
            this.central = central;
            this.latest = latest;
        }

        /** The calendar label of a horizon, as {@code yyyy-MM}. */
        public String label(int monthsOut) {
            int index = Datasets.size() - 1 + monthsOut;
            int year = Datasets.FIRST_YEAR + (Datasets.FIRST_MONTH - 1 + index) / 12;
            int month = (Datasets.FIRST_MONTH - 1 + index) % 12 + 1;
            return String.format(L, "%d-%02d", Integer.valueOf(year), Integer.valueOf(month));
        }

        /** How many years wide the interval is. */
        public double widthInYears() {
            return (latest - earliest) / 12.0;
        }
    }

    /**
     * Filters the whole record, then predicts forward until the level is
     * reached.
     * <p>
     * Forecasting is {@link KalmanFilter#predict()} called repeatedly with no
     * observation to correct against: the mean follows the transition and the
     * covariance grows, which is the widening interval. The deterministic demo
     * answers the same question by finding a root of its trend, and has no
     * interval to offer, because a curve fitted by least squares has no opinion
     * about how far the atmosphere may wander away from it.
     */
    static Crossing crossing(Fit fit, double target) {
        double[] co2 = Datasets.co2();
        int d = dimension(fit.harmonics);
        DMatrix h = observation(fit.harmonics);
        double observationVariance = Math.exp(fit.logVariance[3]);

        KalmanFilter kf = new KalmanFilter(fit.model());
        double[] row = new double[1];
        for (int t = 0; t < co2.length; ++t) {
            if (t > 0) {
                kf.predict();
            }
            row[0] = co2[t];
            kf.update(row);
        }

        int earliest = -1;
        int central = -1;
        int latest = -1;
        double[] mean = new double[d];
        for (int step = 1; step <= 1200 && latest < 0; ++step) {
            kf.predict();
            kf.mean(mean);
            double mu = 0.0;
            for (int j = 0; j < d; ++j) {
                mu += h.get(0, j) * mean[j];
            }
            DMatrix cov = kf.covariance();
            double variance = observationVariance;
            for (int a = 0; a < d; ++a) {
                for (int b = 0; b < d; ++b) {
                    variance += h.get(0, a) * cov.get(a, b) * h.get(0, b);
                }
            }
            double sd = Math.sqrt(variance);
            if (earliest < 0 && mu + 1.96 * sd >= target) {
                earliest = step;
            }
            if (central < 0 && mu >= target) {
                central = step;
            }
            if (latest < 0 && mu - 1.96 * sd >= target) {
                latest = step;
            }
        }
        return new Crossing(target, earliest, central, latest);
    }

    // ------------------------------------------------------------------ main

    private static void title(String text) {
        System.out.println();
        System.out.println("=== " + text);
    }

    /**
     * Runs the eight steps and prints them.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        System.out.println("Mauna Loa CO2 as a state space model");
        System.out.println(Datasets.size() + " monthly means, " + Datasets.label(0) + " to "
                + Datasets.label(Datasets.size() - 1));

        System.out.println();
        System.out.println("MaunaLoaDemo fits this record with a curve -- a quadratic and two");
        System.out.println("harmonics -- and then, in its own step 5, finds that the residuals of");
        System.out.println("that fit remember each other: acf(1) = 0.906. It draws the right");
        System.out.println("conclusion and prints it: its standard errors are therefore too small.");
        System.out.println();
        System.out.println("That is not a bad fit. It is a fit of the wrong shape. A quadratic");
        System.out.println("cannot bend when the atmosphere bends, so everything the curve fails to");
        System.out.println("follow is left in the residuals, and a residual that remembers the");
        System.out.println("previous one is what a misspecified mean looks like.");
        System.out.println();
        System.out.println("So: what happens if the shape is allowed to move? Here the level and");
        System.out.println("the slope are not coefficients fitted once but states that wander, and");
        System.out.println("if that is the right repair then what is left over is only measurement");
        System.out.println("error -- and measurement error has no memory. Section 1 is that test.");
        System.out.println("Everything after it is what the repair also buys.");

        // both fits are needed before any of it can be printed
        Fit two = fit(DETERMINISTIC_HARMONICS);
        Fit fit = fit(HARMONICS);

        title("1. the residuals the other demo could not get rid of");
        Whiteness white = whiteness(fit);
        System.out.printf(L, "  two standard errors for an autocorrelation here: %.4f%n", white.bound);
        System.out.printf(L, "  deterministic trend plus harmonics : acf(1) %7.4f, %2d of 24 lags outside%n",
                white.deterministicAcf1, Integer.valueOf(white.deterministicOutside));
        System.out.printf(L, "  state space innovations            : acf(1) %7.4f, %2d of 24 lags outside%n",
                white.structuralAcf1, Integer.valueOf(white.structuralOutside));
        System.out.printf(L, "  the innovations have mean %.4f and standard deviation %.4f,%n",
                white.innovationMean, white.innovationSd);
        System.out.println("  which is what the model claims they should be.");
        System.out.println();
        System.out.println("  So the repair works, and note what it is not: nothing here improves");
        System.out.println("  the other demo's fit. The correlation has not been removed, it has");
        System.out.printf(L, "  been moved into the model. What is left at lag 12 is %.4f, so the%n",
                white.structuralAcf12);
        System.out.println("  seasonal shape is still not quite captured, and section 2 is where");
        System.out.println("  the model says how much of each month it believes.");

        title("2. what the four variances are, and the one that matters");
        System.out.printf(L, "  %d harmonics (d = %d): log likelihood %12.4f%n", two.harmonics,
                dimension(two.harmonics), two.logLikelihood);
        System.out.printf(L, "  %d harmonics (d = %d): log likelihood %12.4f%n", fit.harmonics,
                dimension(fit.harmonics), fit.logLikelihood);
        System.out.printf(L, "  the third harmonic is worth %.2f nats, so the model keeps it.%n",
                fit.logLikelihood - two.logLikelihood);
        System.out.println("  The deterministic demo has harmonics at 12 and 6 months and stops");
        System.out.println("  there. The likelihood says the series wants one at 4 months too.");
        System.out.println();
        // the iteration count is deliberately not printed; see the class comment
        System.out.printf(L, "  Nelder-Mead over the four log variances, converged: %b%n",
                Boolean.valueOf(fit.converged));
        String[] names = { "level", "slope", "seasonal", "observation" };
        for (int i = 0; i < names.length; ++i) {
            System.out.printf(L, "    sigma %-12s %12.6f ppm%s%n", names[i], fit.sigma(i),
                    i == 1 ? " per month" : "");
        }
        System.out.println();
        System.out.println("  Four numbers is three more than the page needs. What the filter does");
        System.out.println("  is decided by one thing: how far the level may move from one month");
        System.out.printf(L, "  to the next, against how badly a month is measured. Here that is%n");
        System.out.printf(L, "  %.4f, and its consequence is the gain -- the share of each month's%n",
                fit.signalToNoise());
        System.out.printf(L, "  surprise the filter takes as a real move of the atmosphere: %.1f%%.%n",
                100.0 * steadyStateLevelGain(fit));
        System.out.println();
        System.out.println("  That one knob is the whole difference between a curve and a state.");
        System.out.println("  Turn it down and the level is not allowed to move at all, which is");
        System.out.println("  what a coefficient is: a state that may not move. Turn it up and the");
        System.out.println("  level follows every reading. Both ends are here, with the same");
        System.out.println("  machinery and only that number changed:");
        System.out.println();
        System.out.printf(L, "  %-26s %-9s %-10s %-12s%n", "the level may move", "gain", "acf(1)",
                "lags outside");
        Knob[] knobs = knobSweep(fit);
        for (int i = 0; i < knobs.length; ++i) {
            System.out.printf(L, "  %-26s %-9.4f %-10.4f %2d of 24%s%n", knobs[i].label,
                    knobs[i].gain, knobs[i].acf1, Integer.valueOf(knobs[i].outside),
                    knobs[i].fitted ? "   <- the likelihood picked this" : "");
        }
        System.out.println();
        System.out.println("  Read the middle column downwards. Frozen, the innovations remember");
        System.out.println("  each other almost perfectly and every lag is outside the band -- the");
        System.out.println("  model cannot follow the series so everything is left over. That is");
        System.out.println("  worse than MaunaLoaDemo's 0.906, and it should be: a straight line");
        System.out.println("  is a worse curve than a quadratic. Too free, and the correlation");
        System.out.println("  comes back with the other sign, because a level that chases the");
        System.out.println("  noise overshoots and the next surprise undoes the last one.");
        System.out.println();
        System.out.println("  The last column is the blunter instrument of the two and stops");
        System.out.println("  separating the rows at the free end: ten times the fitted movement");
        System.out.println("  leaves as few lags outside the band as the fit does, on an acf(1)");
        System.out.println("  seven times larger. Counting how many lags cross a threshold throws");
        System.out.println("  away how far they crossed it, which is what makes it blunt.");
        System.out.println();
        System.out.println("  White is the middle, and nothing above told the fit where the middle");
        System.out.println("  was. The likelihood chose it without ever being shown an");
        System.out.println("  autocorrelation, and the autocorrelation agrees. Two criteria that");
        System.out.println("  share no arithmetic land on the same number, which is the reason to");
        System.out.println("  believe section 1 rather than merely to read it.");


        title("3. level, slope and season");
        Decomposition parts = decompose(fit);
        int n = Datasets.size();
        System.out.printf(L, "  level at %s : %.3f ppm, rising %.4f ppm per month%n",
                Datasets.label(n - 1), parts.level[n - 1], parts.slope[n - 1]);
        System.out.printf(L, "  that is %.3f ppm per year, against %.3f ppm per year at %s%n",
                12.0 * parts.slope[n - 1], 12.0 * parts.slope[24], Datasets.label(24));
        System.out.printf(L, "  seasonal peak to trough over the last year: %.3f ppm%n",
                parts.seasonalAmplitude);
        System.out.println();
        System.out.println("  NOAA publishes its own deseasonalized series beside the data, made");
        System.out.println("  by a method with nothing in common with this one:");
        System.out.printf(L, "    root mean square difference %.4f ppm over %d months%n",
                parts.rmsAgainstNoaa, n);
        System.out.printf(L, "    worst %.4f ppm, at %s%n", parts.worstAgainstNoaa,
                Datasets.label(parts.worstAt));

        System.out.println();
        System.out.println("The four sections that follow have one thing in common, and it is why");
        System.out.println("they are here rather than in a list of features: a fitted curve can do");
        System.out.println("none of them. Not badly -- at all. A curve has nowhere to put a missing");
        System.out.println("month, no scale on which a surprise is large or small, no way of being");
        System.out.println("told that one observation was better measured than another, and nothing");
        System.out.println("to say when the months are not evenly spaced.");

        title("4. the two months NOAA interpolated");
        Filled filled = fillTheInterpolatedMonths(fit);
        System.out.println("  Blanked out and put back by the smoother, which needs no interpolation");
        System.out.println("  scheme at all -- a NaN is simply a step with no correction:");
        for (int i = 0; i < filled.index.length; ++i) {
            System.out.printf(L, "    %s   NOAA %.4f   smoother %.4f +- %.4f   difference %+.4f%n",
                    Datasets.label(filled.index[i]), filled.noaa[i], filled.smoothed[i], filled.sd[i],
                    filled.smoothed[i] - filled.noaa[i]);
        }

        title("5. the eruption, and what the anomaly score actually found");
        Anomalies found = anomalies(fit);
        System.out.println("  On 2022-11-29 the Mauna Loa observatory was cut off by lava and the");
        System.out.printf(L, "  next %d months were measured at Maunakea, 21 miles north.%n",
                Integer.valueOf(found.moveMonths));
        System.out.printf(L, "    mean standardized innovation over those months : %+.4f%n",
                found.meanDuringMove);
        System.out.printf(L, "    over the rest of the record                    : %+.4f%n",
                found.meanElsewhere);
        System.out.printf(L, "  A step change would show as a shifted mean, and the standard error%n");
        System.out.printf(L, "  of a mean over %d months is %.4f -- so %+.4f is half of one. The move%n",
                Integer.valueOf(found.moveMonths), found.moveStandardError(), found.meanDuringMove);
        System.out.println("  left no offset this model can see.");
        System.out.printf(L, "  One of the eight, %s, does reach %.4f. Over eight months that is%n",
                Datasets.label(found.largestDuringMoveAt), found.largestDuringMove);
        System.out.printf(L, "  a coincidence worth naming rather than a signal, and %d of the whole%n",
                Integer.valueOf(found.beyondThree));
        System.out.println("  record pass three, so the tails are heavier than normal throughout.");
        System.out.println();
        System.out.printf(L, "  The largest anomaly in sixty-eight years is %s, at %.4f standard%n",
                Datasets.label(found.largestAt), found.largest);
        System.out.println("  deviations. That is the 2015-16 El Nino, which produced the largest");
        System.out.println("  annual CO2 growth ever recorded -- so the score found the event");
        System.out.println("  nobody built it to look for, and not the one that was expected.");

        title("6. the uncertainty NOAA publishes, month by month");
        Published published = publishedUncertainty(fit);
        System.out.printf(L, "  %d of %d months carry one: %.3f to %.3f ppm, mean %.3f%n",
                Integer.valueOf(published.months), Integer.valueOf(n), published.minSd,
                published.maxSd, published.meanSd);
        System.out.printf(L, "  the fit of step 1 chose one constant %.4f ppm instead%n",
                published.fittedSd);
        System.out.printf(L, "    log likelihood, constant  : %12.4f%n", published.logLikelihoodConstant);
        System.out.printf(L, "    log likelihood, published : %12.4f%n", published.logLikelihoodPublished);
        System.out.printf(L, "  %.2f nats, and not one extra fitted parameter -- the numbers come%n",
                published.gain());
        System.out.println("  from outside. Nothing was refitted, so that is a lower bound.");
        System.out.printf(L, "  The level itself barely moves: at most %.4f ppm, at %s.%n",
                published.worstLevelShift, Datasets.label(published.worstAt));
        System.out.println("  Better measurements do not sharpen a level whose uncertainty comes");
        System.out.println("  from its own wandering.");

        title("7. seven months in ten thrown away");
        Irregular thin = thinToAnIrregularSchedule(fit);
        System.out.printf(L, "  kept %d of %d months at irregular spacing, largest gap %d months%n",
                Integer.valueOf(thin.kept), Integer.valueOf(thin.total),
                Integer.valueOf(thin.largestGap));
        System.out.println("  Every step now spans a different number of months, so F and Q depend");
        System.out.println("  on the gap. That is what constant matrices cannot express, and it is");
        System.out.println("  not the same thing as a missing month, which costs nothing.");
        System.out.printf(L, "    level against the full record  : rms %.4f ppm, worst %.4f%n",
                thin.rmsAgainstFull, thin.worstAgainstFull);
        System.out.printf(L, "    level against NOAA             : rms %.4f ppm%n", thin.rmsAgainstNoaa);
        System.out.printf(L, "  The whole record reaches %.4f against NOAA, so throwing away %d of%n",
                parts.rmsAgainstNoaa, Integer.valueOf(thin.total - thin.kept));
        System.out.printf(L, "  the %d months cost %.4f ppm.%n", Integer.valueOf(thin.total),
                thin.rmsAgainstNoaa - parts.rmsAgainstNoaa);

        title("8. when does it reach 450 ppm");
        Crossing crossing = crossing(fit, 450.0);
        System.out.printf(L, "  earliest the interval allows : %s%n", crossing.label(crossing.earliest));
        System.out.printf(L, "  the central estimate         : %s%n", crossing.label(crossing.central));
        System.out.printf(L, "  latest the interval allows   : %s%n", crossing.label(crossing.latest));
        System.out.printf(L, "  a window %.1f years wide%n", crossing.widthInYears());
        System.out.println();
        System.out.println("  MaunaLoaDemo answers 2034.42 by finding a root of its fitted trend.");
        System.out.println("  Two models with almost nothing in common -- a quadratic fitted by");
        System.out.println("  least squares, a random walk fitted by maximum likelihood -- land");
        System.out.println("  within a few months of each other, and only one of them says how");
        System.out.println("  wide the answer is.");
    }

    private StateSpaceDemo() {
        throw new AssertionError();
    }
}
