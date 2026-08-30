package math.demo.centralpark;

import java.util.Locale;

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

/**
 * One year of daily maxima and minima at one weather station, as a linear
 * Gaussian state space model.
 * <p>
 * {@link math.demo.co2state.StateSpaceDemo} is the other {@code math.ts} demo,
 * and it observes one number per month. This one observes two per day, which
 * is the whole point: {@code TMAX} and {@code TMIN} are not two series, they
 * are two noisy readings of one temperature, and everything interesting here
 * comes from what the model does with that.
 * <p>
 * The state is built out of parts that mean something on their own:
 *
 * <pre>
 *  0     level        a constant, carrying the mean
 *  1, 2  annual       a rotation at one cycle per year
 *  3, 4  semiannual   a rotation at two
 *  5     half range   a constant, half of TMAX minus TMIN
 *  6, 7  annual       the half range has a season of its own
 *  8, 9  semiannual
 * 10     anomaly      AR(1), the synoptic weather
 * </pre>
 *
 * Everything with a {@code +1} in both rows of {@code H} is the temperature
 * both readings share; everything with {@code +1} against {@code -1} is the
 * spread between them. The seasonal blocks have <b>no process noise at all</b>
 * -- one year is one cycle, and a seasonal component that is allowed to wander
 * cannot be told apart from the level over a single cycle -- so {@code Q} is
 * singular by construction, which {@link LinearGaussianModel} allows and says
 * so.
 * <p>
 * <b>What it finds.</b> The two columns sit side by side in the same file and
 * are worth wildly different amounts. The fit puts {@code 2.74} degrees of
 * error on the maximum and {@code 0.43} on the minimum, because the minimum is
 * set by the air mass before dawn -- which is the thing the state tracks --
 * and the maximum by a short afternoon peak that a cloud can move. Hide the
 * maximum and the weather is {@code 0.3} percent less certain; hide the
 * minimum and it is {@code 2.19} times less. The rest of the page is a
 * precipitation column that hides 35 days behind a flag and changes nothing, a
 * tenth of a degree that does not exist, and an archive that rewrites its own
 * last twelve months.
 * <p>
 * <b>Step 8 is where the model is chosen, and it is deliberately last.</b>
 * Four nested models are fitted by maximum likelihood:
 *
 * <pre>
 * M0  seasonal only, R diagonal          k = 2
 * M1  + the off-diagonal of R            k = 3
 * M2  + the AR(1) anomaly, R diagonal    k = 4
 * M3  + the off-diagonal again           k = 5
 * </pre>
 *
 * The off-diagonal of {@code R} is worth <b>153 nats at M1 and 1.7 at M3</b>,
 * and it changes sign on the way, from {@code +0.757} to exactly {@code -1}. At
 * M1 it is standing in for the common signal, because no state carries it yet.
 * Once the anomaly does, there is nothing left for it to explain, and its
 * optimum walks to the boundary of the parameter space. AIC still prefers M3,
 * by 1.4. <b>A fit that lands on a boundary is a statement about the
 * parameterization and not about the weather</b>, so the demo reports M3 and
 * works with M2.
 * <p>
 * Four further refinements were tried and refused by the data, each measured
 * rather than argued: a second AR(1) component (its standard deviation goes to
 * zero and the likelihood does not move), a level allowed to wander (to
 * {@code 1e-6}), an anomaly that reaches the two readings with different
 * strength (a factor of {@code 0.977}, worth {@code 0.06} nats), and three
 * years of data instead of one (the same story at five times the cost). A
 * model that stops growing is the more useful finding.
 *
 * @see Datasets
 * @see math.demo.co2state.StateSpaceDemo
 */
public final class WeatherDemo {

    private static final Locale L = Locale.ROOT;
    private static final double TWO_PI = 2.0 * Math.PI;

    /** Days in the seasonal period, the mean length of a Gregorian year. */
    static final double PERIOD = 365.25;

    /** The level, a constant that carries the mean of both readings. */
    static final int LEVEL = 0;
    /** First component of the annual rotation the two readings share. */
    static final int ANNUAL = 1;
    /** First component of the semiannual rotation the two readings share. */
    static final int SEMIANNUAL = 3;
    /** Half of TMAX minus TMIN, a constant. */
    static final int HALF_RANGE = 5;
    /** First component of the annual rotation of the half range. */
    static final int HALF_ANNUAL = 6;
    /** First component of the semiannual rotation of the half range. */
    static final int HALF_SEMIANNUAL = 8;
    /** The AR(1) anomaly, present only from rung 2 upwards. */
    static final int ANOMALY = 10;

    /** How many rungs the ladder has. */
    static final int RUNGS = 4;

    /**
     * The rung everything after step 1 works with: four parameters, a diagonal
     * {@code R}, and no parameter on a boundary.
     */
    static final int WORKING = 2;

    /** The rung whose {@code R} has a free off-diagonal on top of the anomaly. */
    static final int UNBOUNDED = 3;

    /**
     * How many terms the deterministic design of {@link #ordinaryLeastSquares}
     * carries: a constant and two harmonics, once shared and once differenced.
     */
    static final int TERMS = 10;

    /**
     * First day step 7 withholds the minimum for, {@code 0} based. Mid July,
     * chosen because the stretch contains the wettest day of the year, so the
     * spread between the two readings is doing something over it rather than
     * sitting on its seasonal value.
     */
    static final int BLANK_FROM = 190;

    /** How many days step 7 withholds. */
    static final int BLANK_LENGTH = 14;

    // -------------------------------------------------------------- the model

    /**
     * Size of the state at one rung: ten deterministic components, plus the
     * anomaly from rung 2 upwards.
     *
     * @param rung
     *            0 to 3
     * @return 10 or 11
     */
    static int dimension(int rung) {
        return rung >= WORKING ? ANOMALY + 1 : ANOMALY;
    }

    /**
     * How many parameters the maximum likelihood fit searches over at one rung.
     *
     * @param rung
     *            0 to 3
     * @return 2, 3, 4 or 5
     */
    static int parameterCount(int rung) {
        return rung == 0 ? 2 : rung == 1 ? 3 : rung == WORKING ? 4 : 5;
    }

    /**
     * Whether the rung lets {@code R} have an off-diagonal.
     *
     * @param rung
     *            0 to 3
     * @return {@code true} at rungs 1 and 3
     */
    static boolean correlatedNoise(int rung) {
        return rung == 1 || rung == UNBOUNDED;
    }

    /** Where the two entries of {@code R} start in the parameter vector. */
    private static int noiseAt(int rung) {
        return rung >= WORKING ? 2 : 0;
    }

    private static void rotate(DMatrix f, int at, double angle) {
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        f.set(at, at, c);
        f.set(at, at + 1, s);
        f.set(at + 1, at, -s);
        f.set(at + 1, at + 1, c);
    }

    /**
     * The AR(1) coefficient, squeezed into {@code (0, 1)} so that the search
     * cannot leave the stationary range and cannot propose a negative daily
     * persistence, which is not a weather model.
     *
     * @param u
     *            the unconstrained parameter
     * @return the coefficient, strictly between 0 and 1
     */
    static double persistence(double u) {
        return 1.0 / (1.0 + Math.exp(-u));
    }

    /**
     * The transition: everything holds its value except the two pairs of
     * rotations, which turn by one day, and the anomaly, which decays.
     *
     * @param rung
     *            0 to 3
     * @param p
     *            the parameters of that rung
     * @return a fresh transition matrix
     */
    static DMatrix transition(int rung, double[] p) {
        int d = dimension(rung);
        DMatrix f = new DMatrix(d, d);
        f.set(LEVEL, LEVEL, 1.0);
        rotate(f, ANNUAL, TWO_PI / PERIOD);
        rotate(f, SEMIANNUAL, 2.0 * TWO_PI / PERIOD);
        f.set(HALF_RANGE, HALF_RANGE, 1.0);
        rotate(f, HALF_ANNUAL, TWO_PI / PERIOD);
        rotate(f, HALF_SEMIANNUAL, 2.0 * TWO_PI / PERIOD);
        if (rung >= WORKING) {
            f.set(ANOMALY, ANOMALY, persistence(p[0]));
        }
        return f;
    }

    /**
     * The process noise, which is zero everywhere but on the anomaly. That is
     * what makes the seasonal part of this model a regression rather than a
     * component that wanders, and it is the reason one year is enough.
     *
     * @param rung
     *            0 to 3
     * @param p
     *            the parameters of that rung
     * @return a fresh process noise matrix, singular at every rung
     */
    static DMatrix processNoise(int rung, double[] p) {
        int d = dimension(rung);
        DMatrix q = new DMatrix(d, d);
        if (rung >= WORKING) {
            q.set(ANOMALY, ANOMALY, Math.exp(p[1]));
        }
        return q;
    }

    /**
     * What is observed. Row 0 is the daily maximum, row 1 the daily minimum;
     * they agree on the level, the two shared rotations and the anomaly, and
     * differ in the sign of the half range and its two rotations.
     *
     * @param rung
     *            0 to 3
     * @return a fresh two-row observation matrix
     */
    static DMatrix observation(int rung) {
        DMatrix h = new DMatrix(2, dimension(rung));
        int[] shared = { LEVEL, ANNUAL, SEMIANNUAL };
        for (int i = 0; i < shared.length; ++i) {
            h.set(0, shared[i], 1.0);
            h.set(1, shared[i], 1.0);
        }
        int[] opposed = { HALF_RANGE, HALF_ANNUAL, HALF_SEMIANNUAL };
        for (int i = 0; i < opposed.length; ++i) {
            h.set(0, opposed[i], 1.0);
            h.set(1, opposed[i], -1.0);
        }
        if (rung >= WORKING) {
            h.set(0, ANOMALY, 1.0);
            h.set(1, ANOMALY, 1.0);
        }
        return h;
    }

    /**
     * The observation noise, as {@code L L'} with {@code L} lower triangular
     * and a positive diagonal. That form is positive definite at every point
     * of the search space, so the fit needs no constraint on it and cannot
     * propose a covariance that is not one.
     *
     * @param rung
     *            0 to 3
     * @param p
     *            the parameters of that rung
     * @return a fresh two by two covariance
     */
    static DMatrix observationNoise(int rung, double[] p) {
        int at = noiseAt(rung);
        double l00 = Math.exp(0.5 * p[at]);
        double l11 = Math.exp(0.5 * p[at + 1]);
        double l10 = correlatedNoise(rung) ? p[at + 2] : 0.0;
        DMatrix r = new DMatrix(2, 2);
        r.set(0, 0, l00 * l00);
        r.set(0, 1, l00 * l10);
        r.set(1, 0, l00 * l10);
        r.set(1, 1, l10 * l10 + l11 * l11);
        return r;
    }

    /**
     * The prior mean: the first day's average for the level, half its range
     * for the half range, zero for everything else.
     *
     * @param rung
     *            0 to 3
     * @return a fresh prior mean
     */
    static double[] priorMean(int rung) {
        double[] m0 = new double[dimension(rung)];
        double first = Datasets.tmax()[0];
        double firstMin = Datasets.tmin()[0];
        m0[LEVEL] = 0.5 * (first + firstMin);
        m0[HALF_RANGE] = 0.5 * (first - firstMin);
        return m0;
    }

    /**
     * A deliberately vague prior. It has to be proper, because {@code P0} is
     * the one matrix the filter factors, but 100 degrees of standard deviation
     * on the level says nothing that 730 observations cannot overrule.
     *
     * @param rung
     *            0 to 3
     * @return a fresh prior covariance, positive definite
     */
    static DMatrix priorCovariance(int rung) {
        int d = dimension(rung);
        DMatrix p0 = new DMatrix(d, d);
        p0.set(LEVEL, LEVEL, 1.0e4);
        for (int i = ANNUAL; i < ANOMALY; ++i) {
            p0.set(i, i, 1.0e2);
        }
        if (rung >= WORKING) {
            p0.set(ANOMALY, ANOMALY, 1.0e2);
        }
        return p0;
    }

    /**
     * The model at one rung and one point of its parameter space.
     *
     * @param rung
     *            0 to 3
     * @param p
     *            the parameters of that rung
     * @return a fresh model
     */
    static LinearGaussianModel model(int rung, double[] p) {
        return new LinearGaussianModel(transition(rung, p), processNoise(rung, p), observation(rung),
                observationNoise(rung, p), priorMean(rung), priorCovariance(rung));
    }

    /**
     * The observations, one day per row: the maximum in column 0, the minimum
     * in column 1.
     *
     * @return a fresh 365 by 2 matrix
     */
    static DMatrix observations() {
        double[] tmax = Datasets.tmax();
        double[] tmin = Datasets.tmin();
        DMatrix y = new DMatrix(tmax.length, 2);
        for (int i = 0; i < tmax.length; ++i) {
            y.set(i, 0, tmax[i]);
            y.set(i, 1, tmin[i]);
        }
        return y;
    }

    // ------------------------------------------ 1. the ladder, by likelihood

    /** What the maximum likelihood fit of one rung returned. */
    static final class Fit {

        /** Which rung of the ladder this is, 0 to 3. */
        public final int rung;
        /** The parameters at the optimum, in the order that rung uses them. */
        public final double[] parameters;
        /** The log likelihood there. */
        public final double logLikelihood;
        /** Whether the simplex converged rather than running out of budget. */
        public final boolean converged;

        Fit(int rung, double[] parameters, double logLikelihood, boolean converged) {
            this.rung = rung;
            this.parameters = parameters;
            this.logLikelihood = logLikelihood;
            this.converged = converged;
        }

        /** The model these parameters describe. */
        public LinearGaussianModel model() {
            return WeatherDemo.model(rung, parameters);
        }

        /** Akaike's criterion, which charges two per parameter. */
        public double aic() {
            return 2.0 * parameterCount(rung) - 2.0 * logLikelihood;
        }

        /** Standard deviation of what the daily maximum gets wrong. */
        public double sigmaMaximum() {
            return Math.sqrt(observationNoise(rung, parameters).get(0, 0));
        }

        /** Standard deviation of what the daily minimum gets wrong. */
        public double sigmaMinimum() {
            return Math.sqrt(observationNoise(rung, parameters).get(1, 1));
        }

        /**
         * Correlation between what the two readings get wrong. Zero by
         * construction at the rungs that do not carry an off-diagonal.
         */
        public double noiseCorrelation() {
            DMatrix r = observationNoise(rung, parameters);
            return r.get(0, 1) / Math.sqrt(r.get(0, 0) * r.get(1, 1));
        }

        /**
         * The AR(1) coefficient of the anomaly, or {@code 0.0} at a rung that
         * has no anomaly.
         */
        public double phi() {
            return rung >= WORKING ? persistence(parameters[0]) : 0.0;
        }

        /** The standard deviation of one day's shock to the anomaly. */
        public double sigmaAnomaly() {
            return rung >= WORKING ? Math.exp(0.5 * parameters[1]) : 0.0;
        }

        /**
         * The standard deviation the anomaly settles at, which is what makes
         * {@link #sigmaAnomaly} comparable with the noise: a small daily shock
         * on a persistent component still adds up.
         */
        public double stationarySigmaAnomaly() {
            double p = phi();
            return sigmaAnomaly() / Math.sqrt(1.0 - p * p);
        }

        /**
         * How many days it takes the anomaly to decay to {@code 1/e} of a
         * shock, which is the readable form of {@link #phi}.
         */
        public double memoryInDays() {
            return -1.0 / Math.log(phi());
        }
    }

    private static double[] start(int rung, double[] previous) {
        switch (rung) {
        case 0:
            return new double[] { Math.log(18.0), Math.log(10.0) };
        case 1:
            return new double[] { previous[0], previous[1], 0.0 };
        case WORKING:
            return new double[] { 0.4, Math.log(6.0), Math.log(7.0), Math.log(0.2) };
        default:
            return new double[] { previous[0], previous[1], previous[2], previous[3], 0.0 };
        }
    }

    /**
     * Maximizes {@link KalmanFilter#logLikelihood} over the parameters of one
     * rung.
     * <p>
     * The search runs over logarithms of the variances, because a variance is
     * positive and the likelihood is flat in it over orders of magnitude, and
     * over the lower triangle of {@code R} rather than over {@code R} itself.
     * The box is a refusal rather than a penalty: a variance of {@code e^-40}
     * is not a model anyone meant.
     * <p>
     * The simplex is given <b>no restart</b>, unlike the default. That was
     * measured rather than assumed: five starting points spread over the
     * parameter space, with and without a restart, all reach the same optimum
     * to six decimal places, and dropping the restart cuts the fit of the four
     * rungs from 8.0 to 4.6 seconds.
     *
     * @param rung
     *            0 to 3
     * @param previous
     *            the optimum of the rung below, used as the starting point
     *            where the two are nested; may be {@code null} at rungs 0 and 2
     * @return what the search settled on
     */
    static Fit fit(final int rung, double[] previous) {
        return fit(rung, previous, observations());
    }

    /**
     * The same fit against a series the caller supplies, which is what lets
     * step 4 refit the model on deliberately coarsened readings.
     *
     * @param rung
     *            0 to 3
     * @param previous
     *            the optimum of the rung below, or {@code null}
     * @param y
     *            the observations to fit against. Not modified
     * @return what the search settled on
     */
    static Fit fit(final int rung, double[] previous, final DMatrix y) {
        DMultiFunction negativeLogLikelihood = new DMultiFunction() {
            @Override
            public double apply(double[] x) {
                for (int i = 0; i < x.length; ++i) {
                    if (x[i] < -40.0 || x[i] > 40.0) {
                        return 1.0e30;
                    }
                }
                return -KalmanFilter.logLikelihood(model(rung, x), y);
            }
        };
        DMultiFunctionEval found = new NelderMead(1.0e-10, 1.0e-8, 0, 0)
                .minimize(negativeLogLikelihood, start(rung, previous));
        return new Fit(rung, found.point, -found.value, found.converged);
    }

    /**
     * The four rungs, fitted in order so that each nested pair starts from the
     * optimum of the one below it.
     *
     * @return the four fits, indexed by rung
     */
    static Fit[] ladder() {
        Fit[] fits = new Fit[RUNGS];
        double[] previous = null;
        for (int rung = 0; rung < RUNGS; ++rung) {
            fits[rung] = fit(rung, previous);
            previous = fits[rung].parameters;
        }
        return fits;
    }

    // ------------------------------------------------------ 2. take it apart

    /** What the smoother makes of the year, component by component. */
    static final class Decomposition {

        /** The temperature both readings share, level plus season, per day. */
        public final double[] seasonal;
        /** The synoptic anomaly, per day. */
        public final double[] anomaly;
        /** Half the spread between the two readings, per day. */
        public final double[] halfRange;
        /** What the model says the maximum was, noise removed. */
        public final double[] maximum;
        /** What the model says the minimum was, noise removed. */
        public final double[] minimum;
        /** Peak to trough of {@link #seasonal} over the year. */
        public final double seasonalSwing;
        /** Peak to trough of {@link #halfRange} over the year. */
        public final double halfRangeSwing;
        /** The day {@link #halfRange} peaks, {@code 0} based. */
        public final int widestDay;
        /** The day {@link #seasonal} peaks, {@code 0} based. */
        public final int warmestDay;
        /** The largest anomaly of the year, in degrees. */
        public final double largestAnomaly;
        /** The day it fell on, {@code 0} based. */
        public final int largestAnomalyDay;

        Decomposition(double[] seasonal, double[] anomaly, double[] halfRange, double[] maximum,
                double[] minimum, double seasonalSwing, double halfRangeSwing, int widestDay,
                int warmestDay, double largestAnomaly, int largestAnomalyDay) {
            this.seasonal = seasonal;
            this.anomaly = anomaly;
            this.halfRange = halfRange;
            this.maximum = maximum;
            this.minimum = minimum;
            this.seasonalSwing = seasonalSwing;
            this.halfRangeSwing = halfRangeSwing;
            this.widestDay = widestDay;
            this.warmestDay = warmestDay;
            this.largestAnomaly = largestAnomaly;
            this.largestAnomalyDay = largestAnomalyDay;
        }
    }

    /**
     * Runs the smoother and splits the state into the three things it was
     * built out of.
     *
     * @param fit
     *            a fitted model; the anomaly is zero unless it carries one
     * @return the decomposition
     */
    static Decomposition decompose(Fit fit) {
        LinearGaussianModel m = fit.model();
        RtsSmoother.Result s = RtsSmoother.smooth(m, observations());
        int n = s.length;
        double[] seasonal = new double[n];
        double[] anomaly = new double[n];
        double[] halfRange = new double[n];
        double[] maximum = new double[n];
        double[] minimum = new double[n];
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;
        double loHalf = Double.MAX_VALUE;
        double hiHalf = -Double.MAX_VALUE;
        int warmest = 0;
        int widest = 0;
        double largest = 0.0;
        int largestAt = 0;
        for (int t = 0; t < n; ++t) {
            double[] x = s.mean[t];
            seasonal[t] = x[LEVEL] + x[ANNUAL] + x[SEMIANNUAL];
            halfRange[t] = x[HALF_RANGE] + x[HALF_ANNUAL] + x[HALF_SEMIANNUAL];
            anomaly[t] = fit.rung >= WORKING ? x[ANOMALY] : 0.0;
            maximum[t] = seasonal[t] + anomaly[t] + halfRange[t];
            minimum[t] = seasonal[t] + anomaly[t] - halfRange[t];
            if (seasonal[t] > hi) {
                hi = seasonal[t];
                warmest = t;
            }
            if (seasonal[t] < lo) {
                lo = seasonal[t];
            }
            if (halfRange[t] > hiHalf) {
                hiHalf = halfRange[t];
                widest = t;
            }
            loHalf = Math.min(loHalf, halfRange[t]);
            if (Math.abs(anomaly[t]) > Math.abs(largest)) {
                largest = anomaly[t];
                largestAt = t;
            }
        }
        return new Decomposition(seasonal, anomaly, halfRange, maximum, minimum, hi - lo,
                hiHalf - loHalf, widest, warmest, largest, largestAt);
    }

    // ------------------------------- the deterministic fit this is measured against

    /**
     * The same ten terms as rung 0, fitted by ordinary least squares over the
     * 730 readings stacked into one column.
     * <p>
     * It is the model of rung 0 written the other way round, and it is here as
     * an independent answer: a filter whose process noise is zero everywhere is
     * doing recursive least squares, so its final state has to be what
     * {@link OLS} returns from the same design.
     *
     * @return the least squares summary of the stacked design
     */
    static LSSummary ordinaryLeastSquares() {
        int n = Datasets.size();
        DMatrix x = new DMatrix(2 * n, TERMS);
        DMatrix y = new DMatrix(2 * n, 1);
        double[] tmax = Datasets.tmax();
        double[] tmin = Datasets.tmin();
        for (int t = 0; t < n; ++t) {
            double w = TWO_PI * t / PERIOD;
            double[] shared = { 1.0, Math.cos(w), Math.sin(w), Math.cos(2.0 * w), Math.sin(2.0 * w) };
            for (int j = 0; j < shared.length; ++j) {
                x.set(t, j, shared[j]);
                x.set(n + t, j, shared[j]);
                x.set(t, 5 + j, shared[j]);
                x.set(n + t, 5 + j, -shared[j]);
            }
            y.set(t, 0, tmax[t]);
            y.set(n + t, 0, tmin[t]);
        }
        return OLS.estimate(0.05, x, y);
    }

    /**
     * The residuals of {@link #ordinaryLeastSquares}, averaged over the two
     * readings: the series whose memory the state space model has to absorb.
     *
     * @return one residual per day
     */
    static double[] deterministicResiduals() {
        LSSummary summary = ordinaryLeastSquares();
        int n = Datasets.size();
        double[] common = new double[n];
        for (int t = 0; t < n; ++t) {
            common[t] = 0.5 * (summary.getResiduals().get(t) + summary.getResiduals().get(n + t));
        }
        return common;
    }

    // ------------------------------------------------- 3. are they white now

    /** How much memory is left in what the model could not predict. */
    static final class Whiteness {

        /** Two standard errors for an autocorrelation over this many days. */
        public final double bound;
        /** First autocorrelation of the least squares residuals, per reading. */
        public final double[] deterministicAcf1;
        /** How many of 20 lags leave the band, for the same residuals. */
        public final int[] deterministicOutside;
        /** First autocorrelation of the standardized innovations. */
        public final double[] innovationAcf1;
        /** How many of 20 lags leave the band, for the innovations. */
        public final int[] innovationOutside;
        /** Mean of the standardized innovations, which should be zero. */
        public final double[] innovationMean;
        /** Standard deviation of them, which should be one. */
        public final double[] innovationSd;

        Whiteness(double bound, double[] deterministicAcf1, int[] deterministicOutside,
                double[] innovationAcf1, int[] innovationOutside, double[] innovationMean,
                double[] innovationSd) {
            this.bound = bound;
            this.deterministicAcf1 = deterministicAcf1;
            this.deterministicOutside = deterministicOutside;
            this.innovationAcf1 = innovationAcf1;
            this.innovationOutside = innovationOutside;
            this.innovationMean = innovationMean;
            this.innovationSd = innovationSd;
        }
    }

    /**
     * The innovations, whitened by the Cholesky factor of their own covariance
     * so that the two readings become two independent standard normal series
     * if the model is right.
     * <p>
     * A univariate demo can divide by a standard deviation; with two readings
     * at once the innovations are correlated through {@code S}, and dividing
     * each by its own square root would leave that correlation in place. The
     * lower triangular factor is what removes it.
     *
     * @param filtered
     *            a filter run over the whole year
     * @return two rows of standardized innovations, maximum first
     */
    static double[][] standardizedInnovations(KalmanFilter.Result filtered) {
        int n = filtered.length;
        double[][] z = new double[2][n];
        for (int t = 0; t < n; ++t) {
            DMatrix s = filtered.innovationCovariance[t];
            double l00 = Math.sqrt(s.get(0, 0));
            double l10 = s.get(1, 0) / l00;
            double l11 = Math.sqrt(s.get(1, 1) - l10 * l10);
            z[0][t] = filtered.innovation[t][0] / l00;
            z[1][t] = (filtered.innovation[t][1] - l10 * z[0][t]) / l11;
        }
        return z;
    }

    /**
     * Compares what the state space model could not predict with what the
     * least squares fit of the same ten terms left over.
     *
     * @param fit
     *            the working model
     * @return the comparison
     */
    static Whiteness whiteness(Fit fit) {
        int n = Datasets.size();
        double bound = 2.0 / Math.sqrt(n);
        LSSummary summary = ordinaryLeastSquares();
        double[][] deterministic = new double[2][n];
        for (int t = 0; t < n; ++t) {
            deterministic[0][t] = summary.getResiduals().get(t);
            deterministic[1][t] = summary.getResiduals().get(n + t);
        }
        double[][] z = standardizedInnovations(KalmanFilter.filter(fit.model(), observations()));
        double[] detAcf = new double[2];
        int[] detOut = new int[2];
        double[] innAcf = new double[2];
        int[] innOut = new int[2];
        double[] mean = new double[2];
        double[] sd = new double[2];
        for (int c = 0; c < 2; ++c) {
            double[] a = ACF.acf(deterministic[c], 20);
            detAcf[c] = a[1];
            detOut[c] = countOutside(a, bound);
            double[] b = ACF.acf(z[c], 20);
            innAcf[c] = b[1];
            innOut[c] = countOutside(b, bound);
            mean[c] = mean(z[c]);
            sd[c] = standardDeviation(z[c]);
        }
        return new Whiteness(bound, detAcf, detOut, innAcf, innOut, mean, sd);
    }

    private static int countOutside(double[] acf, double bound) {
        int outside = 0;
        for (int lag = 1; lag < acf.length; ++lag) {
            if (Math.abs(acf[lag]) > bound) {
                ++outside;
            }
        }
        return outside;
    }

    // --------------------------------------------- 4. what the grid is worth

    /** What the resolution of the column costs the model. */
    static final class Grid {

        /** Standard deviation of rounding to the Fahrenheit grid. */
        public final double quantizationSd;
        /** Share of the fitted maximum noise variance that rounding explains. */
        public final double shareOfMaximum;
        /** The same share for the minimum. */
        public final double shareOfMinimum;
        /** Standard deviation of rounding to whole degrees Celsius. */
        public final double coarseSd;
        /**
         * Standard deviation the coarsening actually has, measured from the
         * differences rather than assumed to be {@code d / sqrt(12)}.
         */
        public final double measuredCoarseSd;
        /** What the minimum's noise should become on the coarser grid. */
        public final double predictedMinimum;
        /** The same prediction, from the coarsening that was measured. */
        public final double predictedFromMeasured;
        /** What it does become when the model is refitted there. */
        public final double refittedMinimum;
        /** And the maximum's, refitted. */
        public final double refittedMaximum;

        Grid(double quantizationSd, double shareOfMaximum, double shareOfMinimum, double coarseSd,
                double measuredCoarseSd, double predictedMinimum, double predictedFromMeasured,
                double refittedMinimum, double refittedMaximum) {
            this.predictedFromMeasured = predictedFromMeasured;
            this.quantizationSd = quantizationSd;
            this.shareOfMaximum = shareOfMaximum;
            this.shareOfMinimum = shareOfMinimum;
            this.coarseSd = coarseSd;
            this.measuredCoarseSd = measuredCoarseSd;
            this.predictedMinimum = predictedMinimum;
            this.refittedMinimum = refittedMinimum;
            this.refittedMaximum = refittedMaximum;
        }
    }

    /**
     * The observations rounded onto a coarser grid than the one they arrive
     * on, so that the model can be asked what a resolution costs it.
     *
     * @param step
     *            the spacing to round to, in degrees
     * @return a fresh 365 by 2 matrix
     */
    static DMatrix roundedObservations(double step) {
        DMatrix y = observations();
        for (int t = 0; t < y.numRows(); ++t) {
            for (int c = 0; c < 2; ++c) {
                y.set(t, c, step * Math.rint(y.get(t, c) / step));
            }
        }
        return y;
    }

    /**
     * Measures the resolution of the column against the noise the model
     * attributes to it, and checks the one quantitative prediction that
     * follows: a uniform rounding of width {@code d} adds {@code d^2 / 12} to
     * a variance, so refitting on a coarser grid has to move the fitted noise
     * by a known amount, or the model is not measuring what it claims to.
     *
     * @param fit
     *            the working model
     * @return the comparison
     */
    static Grid grid(Fit fit) {
        double fine = Datasets.FAHRENHEIT_STEP / 10.0 / Math.sqrt(12.0);
        double coarse = 1.0 / Math.sqrt(12.0);
        double varMax = fit.sigmaMaximum() * fit.sigmaMaximum();
        double varMin = fit.sigmaMinimum() * fit.sigmaMinimum();
        double predicted = Math.sqrt(varMin - fine * fine + coarse * coarse);
        DMatrix rounded = roundedObservations(1.0);
        DMatrix exact = observations();
        double[] shift = new double[2 * exact.numRows()];
        for (int t = 0; t < exact.numRows(); ++t) {
            shift[2 * t] = rounded.get(t, 0) - exact.get(t, 0);
            shift[2 * t + 1] = rounded.get(t, 1) - exact.get(t, 1);
        }
        double measured = standardDeviation(shift);
        Fit refitted = fit(WORKING, null, rounded);
        return new Grid(fine, fine * fine / varMax, fine * fine / varMin, coarse, measured, predicted,
                Math.sqrt(varMin - fine * fine + measured * measured), refitted.sigmaMinimum(),
                refitted.sigmaMaximum());
    }

    // ------------------------------------- 5. the column that changes nothing

    /** The three kinds of day the precipitation column really holds. */
    static final class Trace {

        /** How many days had measurable rain, a trace, and none. */
        public final int[] count;
        /** Mean spread between the readings, season removed, per kind. */
        public final double[] meanSpread;
        /** Standard error of each of those means. */
        public final double[] standardError;
        /** Trace against dry, in standard errors. */
        public final double traceAgainstDry;
        /** Trace against wet, in standard errors. */
        public final double traceAgainstWet;

        Trace(int[] count, double[] meanSpread, double[] standardError, double traceAgainstDry,
                double traceAgainstWet) {
            this.count = count;
            this.meanSpread = meanSpread;
            this.standardError = standardError;
            this.traceAgainstDry = traceAgainstDry;
            this.traceAgainstWet = traceAgainstWet;
        }
    }

    /**
     * Splits the year by what the precipitation column says, using the flag
     * as well as the number, and asks whether the spread between the two
     * readings can tell the three groups apart.
     * <p>
     * A wet day should have a narrow spread: cloud holds the maximum down and
     * the minimum up. The question is where a trace falls, and the answer
     * decides whether reading the number without the flag costs anything.
     *
     * @param parts
     *            a decomposition, whose half range removes the season
     * @return the three groups and the two comparisons
     */
    static Trace trace(Decomposition parts) {
        int n = Datasets.size();
        double[] tmax = Datasets.tmax();
        double[] tmin = Datasets.tmin();
        int[] count = new int[3];
        double[] sum = new double[3];
        double[] sumSquares = new double[3];
        for (int t = 0; t < n; ++t) {
            int kind = Datasets.isWet(t) ? 0 : Datasets.isTrace(t) ? 1 : 2;
            double residual = tmax[t] - tmin[t] - 2.0 * parts.halfRange[t];
            ++count[kind];
            sum[kind] += residual;
            sumSquares[kind] += residual * residual;
        }
        double[] mean = new double[3];
        double[] error = new double[3];
        for (int k = 0; k < 3; ++k) {
            mean[k] = sum[k] / count[k];
            double variance = (sumSquares[k] - count[k] * mean[k] * mean[k]) / (count[k] - 1);
            error[k] = Math.sqrt(variance / count[k]);
        }
        double againstDry = (mean[1] - mean[2]) / Math.hypot(error[1], error[2]);
        double againstWet = (mean[1] - mean[0]) / Math.hypot(error[1], error[0]);
        return new Trace(count, mean, error, againstDry, againstWet);
    }

    // ------------------------------------ 6. what the second reading is worth

    /** How sharply the weather is known from one reading and from two. */
    static final class Channels {

        /** Posterior standard deviation of the anomaly from both readings. */
        public final double fromBoth;
        /** From the daily maximum alone. */
        public final double fromMaximum;
        /** From the daily minimum alone. */
        public final double fromMinimum;

        Channels(double fromBoth, double fromMaximum, double fromMinimum) {
            this.fromBoth = fromBoth;
            this.fromMaximum = fromMaximum;
            this.fromMinimum = fromMinimum;
        }
    }

    /**
     * The observations with one column blanked out for the whole year, which
     * is what a {@code NaN} means to {@link KalmanFilter}.
     *
     * @param keep
     *            0 to keep the maximum, 1 to keep the minimum
     * @return a fresh 365 by 2 matrix with the other column all {@code NaN}
     */
    static DMatrix onlyOneReading(int keep) {
        DMatrix y = observations();
        for (int t = 0; t < y.numRows(); ++t) {
            y.set(t, 1 - keep, Double.NaN);
        }
        return y;
    }

    /**
     * How much of the weather each reading carries on its own.
     * <p>
     * The level and the half range both enter a single reading with the same
     * sign, so from one column they are not separately identified and their
     * posterior stays wide. The anomaly is the component that is identified
     * either way, and it is the one that means weather, so it is the honest
     * thing to compare.
     *
     * @param fit
     *            the working model
     * @return the three posterior widths
     */
    static Channels channels(Fit fit) {
        LinearGaussianModel m = fit.model();
        return new Channels(anomalyWidth(m, observations()), anomalyWidth(m, onlyOneReading(0)),
                anomalyWidth(m, onlyOneReading(1)));
    }

    private static double anomalyWidth(LinearGaussianModel m, DMatrix y) {
        RtsSmoother.Result s = RtsSmoother.smooth(m, y);
        double sum = 0.0;
        for (int t = 0; t < s.length; ++t) {
            sum += s.covariance[t].get(ANOMALY, ANOMALY);
        }
        return Math.sqrt(sum / s.length);
    }

    // ---------------------------------- 7. the state, and the reading that is not there

    /** What happens when one of the two readings goes missing for a while. */
    static final class Blanked {

        /** First day withheld, {@code 0} based. */
        public final int from;
        /** How many days were withheld. */
        public final int length;
        /** Root mean square error of the minimum the smoother implies. */
        public final double plainError;
        /** The same, after adding the term a correlated R would contribute. */
        public final double correctedError;
        /** The largest that term reaches over the gap, in degrees. */
        public final double largestTerm;
        /** The correlation in R that the term is proportional to. */
        public final double correlation;

        Blanked(int from, int length, double plainError, double correctedError, double largestTerm,
                double correlation) {
            this.from = from;
            this.length = length;
            this.plainError = plainError;
            this.correctedError = correctedError;
            this.largestTerm = largestTerm;
            this.correlation = correlation;
        }
    }

    /**
     * Withholds the daily minimum over one stretch and asks the model to put
     * it back.
     * <p>
     * The point is a distinction that is easy to miss: a Kalman filter given
     * only some components of an observation returns the <em>state</em>
     * correctly, because the state is what the observed components carry
     * information about. It does not return the conditional mean of the
     * component that was withheld. That one is
     * <code>H1 x + R10 / (H0 P H0' + R00) * innovation</code>, and the filter
     * cannot supply the second term, because the reduced update never forms a
     * residual for the row that is missing. With a diagonal {@code R} the term
     * is zero and nothing is lost; with a correlated one it is not.
     *
     * @param fit
     *            the model to reconstruct with
     * @param from
     *            first day to withhold
     * @param length
     *            how many days
     * @return what each answer cost
     */
    static Blanked blank(Fit fit, int from, int length) {
        LinearGaussianModel m = fit.model();
        DMatrix r = observationNoise(fit.rung, fit.parameters);
        DMatrix h = observation(fit.rung);
        DMatrix y = observations();
        double[] truth = Datasets.tmin();
        for (int t = from; t < from + length; ++t) {
            y.set(t, 1, Double.NaN);
        }
        KalmanFilter.Result filtered = KalmanFilter.filter(m, y);
        double plain = 0.0;
        double corrected = 0.0;
        double largest = 0.0;
        for (int t = from; t < from + length; ++t) {
            double implied = 0.0;
            double predicted = 0.0;
            for (int j = 0; j < dimension(fit.rung); ++j) {
                implied += h.get(1, j) * filtered.filteredMean[t][j];
                predicted += h.get(0, j) * filtered.predictedMean[t][j];
            }
            double variance = quadratic(h, filtered.predictedCovariance[t]) + r.get(0, 0);
            double term = r.get(1, 0) / variance * (y.get(t, 0) - predicted);
            largest = Math.max(largest, Math.abs(term));
            plain += (implied - truth[t]) * (implied - truth[t]);
            corrected += (implied + term - truth[t]) * (implied + term - truth[t]);
        }
        return new Blanked(from, length, Math.sqrt(plain / length), Math.sqrt(corrected / length),
                largest, fit.noiseCorrelation());
    }

    /** Row 0 of {@code H} times {@code P} times its transpose. */
    private static double quadratic(DMatrix h, DMatrix p) {
        double total = 0.0;
        for (int i = 0; i < p.numRows(); ++i) {
            for (int j = 0; j < p.numColumns(); ++j) {
                total += h.get(0, i) * p.get(i, j) * h.get(0, j);
            }
        }
        return total;
    }

    // --------------------------------------- 8. the day the source changed

    /** Whether the change of source left a mark the model can see. */
    static final class Break {

        /** Mean standardized innovation before the break, per reading. */
        public final double[] before;
        /** And after it. */
        public final double[] after;
        /** Standard error of the difference, per reading. */
        public final double[] standardError;
        /** The difference in standard errors, per reading. */
        public final double[] ratio;

        Break(double[] before, double[] after, double[] standardError, double[] ratio) {
            this.before = before;
            this.after = after;
            this.standardError = standardError;
            this.ratio = ratio;
        }
    }

    /**
     * Splits the standardized innovations at {@link Datasets#SOURCE_BREAK} and
     * asks whether the two stretches look like the same measurement process.
     * <p>
     * This is the one place where the demo knows in advance that something
     * changed, and knows exactly when. If a documented change of source leaves
     * no mark, that is worth as much as finding one.
     *
     * @param fit
     *            the working model
     * @return the comparison, per reading
     */
    static Break sourceBreak(Fit fit) {
        double[][] z = standardizedInnovations(KalmanFilter.filter(fit.model(), observations()));
        int split = Datasets.SOURCE_BREAK;
        int n = Datasets.size();
        double[] before = new double[2];
        double[] after = new double[2];
        double[] error = new double[2];
        double[] ratio = new double[2];
        for (int c = 0; c < 2; ++c) {
            double sumBefore = 0.0;
            double sumAfter = 0.0;
            for (int t = 0; t < split; ++t) {
                sumBefore += z[c][t];
            }
            for (int t = split; t < n; ++t) {
                sumAfter += z[c][t];
            }
            before[c] = sumBefore / split;
            after[c] = sumAfter / (n - split);
            error[c] = Math.sqrt(1.0 / split + 1.0 / (n - split));
            ratio[c] = (after[c] - before[c]) / error[c];
        }
        return new Break(before, after, error, ratio);
    }

    // ------------------------------------------------------------------ main

    private static void title(String text) {
        System.out.println();
        System.out.println("=== " + text);
    }

    /**
     * Runs the steps and prints them.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        System.out.println("A year of daily maxima and minima at one weather station");
        System.out.println(Datasets.size() + " days at " + Datasets.STATION + ", " + Datasets.NAME
                + ", " + Datasets.label(0) + " to " + Datasets.label(Datasets.size() - 1));
        System.out.println();
        System.out.println("The two columns are not two measurements. They are one temperature,");
        System.out.println("read twice a day, and every question below comes out of that.");

        Fit[] rungs = ladder();
        Fit working = rungs[WORKING];
        String[] readings = { "maximum", "minimum" };

        title("1. one temperature, read twice, and the two readings are not equal");
        System.out.println("  The model carries a temperature both readings see, a spread between");
        System.out.println("  them that the sun opens by day and the night closes, an annual cycle");
        System.out.println("  for each of those, and a wander for the weather. Four numbers are");
        System.out.println("  left for the data to fill in, and it fills them in like this:");
        System.out.println();
        System.out.printf(L, "    the maximum is wrong by         %7.4f degrees%n",
                Double.valueOf(working.sigmaMaximum()));
        System.out.printf(L, "    the minimum is wrong by         %7.4f degrees%n",
                Double.valueOf(working.sigmaMinimum()));
        System.out.printf(L, "    one day's push on the weather   %7.4f degrees%n",
                Double.valueOf(working.sigmaAnomaly()));
        System.out.printf(L, "    tomorrow keeps                  %7.4f of today%n",
                Double.valueOf(working.phi()));
        System.out.println();
        System.out.printf(L, "  The first two are the surprise: the night minimum is %.1f times%n",
                Double.valueOf(working.sigmaMaximum() / working.sigmaMinimum()));
        System.out.println("  sharper than the afternoon maximum. Not because of the thermometer.");
        System.out.println("  The minimum is set by the air mass before dawn, which is the thing");
        System.out.println("  the model is following; the maximum by a short afternoon peak that a");
        System.out.println("  passing cloud can move.");
        System.out.println();
        System.out.println("  What that is worth: hide one column and ask how well the model still");
        System.out.println("  knows what the weather did.");
        Channels pair = channels(working);
        System.out.println();
        System.out.printf(L, "    from both readings       %.4f degrees%n",
                Double.valueOf(pair.fromBoth));
        System.out.printf(L, "    from the maximum alone   %.4f%n", Double.valueOf(pair.fromMaximum));
        System.out.printf(L, "    from the minimum alone   %.4f%n", Double.valueOf(pair.fromMinimum));
        System.out.println();
        System.out.printf(L, "  Throw the maximum away and you lose %.1f percent. Throw the minimum%n",
                Double.valueOf(100.0 * (pair.fromMinimum / pair.fromBoth - 1.0)));
        System.out.printf(L, "  away and the weather is %.2f times less certain. Two columns that%n",
                Double.valueOf(pair.fromMaximum / pair.fromBoth));
        System.out.println("  look alike, side by side in the same file, worth wildly different");
        System.out.println("  amounts.");

        title("2. season, weather, and the spread between the readings");
        Decomposition parts = decompose(working);
        System.out.println("  The model takes the year apart into the three pieces it was built");
        System.out.println("  from, and each of them can be read on its own:");
        System.out.println();
        System.out.printf(L, "    the shared temperature  swings %.2f degrees, warmest about %s%n",
                Double.valueOf(parts.seasonalSwing), Datasets.label(parts.warmestDay));
        System.out.printf(L, "    the spread between them averages %.2f, plus %.2f for its season%n",
                Double.valueOf(2.0 * mean(parts.halfRange)),
                Double.valueOf(2.0 * parts.halfRangeSwing));
        System.out.printf(L, "    the weather             furthest out %+.2f on %s%n",
                Double.valueOf(parts.largestAnomaly), Datasets.label(parts.largestAnomalyDay));
        System.out.println();
        System.out.printf(L, "  The spread is widest around %s, when the sun is already%n",
                Datasets.label(parts.widestDay));
        System.out.println("  strong and the nights are still long, and the day the weather");
        System.out.println("  departs furthest is also the coldest night of the year.");
        System.out.println();
        System.out.printf(L, "  The weather settles at %.2f degrees either side of the season, and a%n",
                Double.valueOf(working.stationarySigmaAnomaly()));
        System.out.printf(L, "  push on it has faded to a third of itself after %.2f days. Both of%n",
                Double.valueOf(working.memoryInDays()));
        System.out.println("  those come out of the fit, not out of a choice.");

        title("3. thirty-five days the number column cannot describe");
        Trace kinds = trace(parts);
        System.out.println("  A zero in the precipitation column is two different days: it did not");
        System.out.println("  rain, or it rained too little to measure and a separate flag says so.");
        System.out.println("  Cloud holds the maximum down and the minimum up, so if a trace day is");
        System.out.println("  really a wet day, the spread between the readings should show it.");
        System.out.println();
        String[] kindNames = { "measurable rain", "a trace", "dry" };
        for (int k = 0; k < 3; ++k) {
            System.out.printf(L, "    %-16s %3d days, spread %+.3f degrees from its season (se %.3f)%n",
                    kindNames[k], Integer.valueOf(kinds.count[k]), Double.valueOf(kinds.meanSpread[k]),
                    Double.valueOf(kinds.standardError[k]));
        }
        System.out.println();
        System.out.printf(L, "  A trace sits %.2f standard errors from a dry day and %.2f from a wet%n",
                Double.valueOf(kinds.traceAgainstDry), Double.valueOf(kinds.traceAgainstWet));
        System.out.println("  one. So a trace is a dry day, and anyone who read the number column");
        System.out.println("  and ignored the flag put these 35 days where they belonged anyway.");
        System.out.println();
        System.out.println("  Worth saying plainly, because the other demo in this package found");
        System.out.println("  the opposite: there a column nobody had used moved a p value by a");
        System.out.println("  factor of 241. A hidden column that changes nothing and one that");
        System.out.println("  changes everything look exactly alike until you use them.");

        title("4. the tenth of a degree that is not there");
        Grid resolution = grid(working);
        System.out.println("  The columns are written in tenths of a degree Celsius, which looks");
        System.out.printf(L, "  like the resolution. It is not. The values are whole degrees%n");
        System.out.printf(L, "  Fahrenheit, %.4f degrees apart, so the real grid is five and a half%n",
                Double.valueOf(Datasets.FAHRENHEIT_STEP / 10.0));
        System.out.printf(L, "  times coarser, and rounding to it costs %.4f degrees on its own.%n",
                Double.valueOf(resolution.quantizationSd));
        System.out.println("  Against the error the model assigns to each reading:");
        System.out.println();
        System.out.printf(L, "    of the maximum's error, rounding is %5.2f percent%n",
                Double.valueOf(100.0 * resolution.shareOfMaximum));
        System.out.printf(L, "    of the minimum's error, rounding is %5.2f percent%n",
                Double.valueOf(100.0 * resolution.shareOfMinimum));
        System.out.println();
        System.out.println("  So for the minimum the way the number was written down is a seventh");
        System.out.println("  of everything that is wrong with it, and for the maximum it vanishes.");
        System.out.println();
        System.out.println("  There is a rule of thumb for this, d^2/12, and it can be tested:");
        System.out.println("  round the readings to whole degrees Celsius and fit again.");
        System.out.println();
        System.out.printf(L, "    the rule predicts the minimum's error goes to  %.4f%n",
                Double.valueOf(resolution.predictedMinimum));
        System.out.printf(L, "    measuring the rounding instead predicts        %.4f%n",
                Double.valueOf(resolution.predictedFromMeasured));
        System.out.printf(L, "    refitting says                                 %.4f%n",
                Double.valueOf(resolution.refittedMinimum));
        System.out.printf(L, "    and the maximum does not move, %.4f against %.4f%n",
                Double.valueOf(resolution.refittedMaximum), Double.valueOf(working.sigmaMaximum()));
        System.out.println();
        System.out.println("  Both predictions miss. The rule treats rounding as extra noise, and");
        System.out.println("  it is not noise: it is a fixed function of the value being rounded,");
        System.out.println("  so it is not independent of anything. It is a rule about continuous");
        System.out.println("  data, and these columns were already on a grid before we touched");
        System.out.println("  them.");

        title("5. the day the source changed");
        Whiteness white = whiteness(working);
        Break change = sourceBreak(working);
        System.out.printf(L, "  On %s this file stops coming from the airport's automatic%n",
                Datasets.label(Datasets.SOURCE_BREAK));
        System.out.printf(L, "  station and starts coming from a weather service form, %d days%n",
                Integer.valueOf(Datasets.SOURCE_BREAK));
        System.out.printf(L, "  against %d. The model was told nothing about it. Its one day ahead%n",
                Integer.valueOf(Datasets.size() - Datasets.SOURCE_BREAK));
        System.out.println("  errors, which should sit around zero if nothing changed:");
        System.out.println();
        for (int c = 0; c < 2; ++c) {
            System.out.printf(L, "    %-8s %+.4f before, %+.4f after, a shift of %+.2f standard errors%n",
                    readings[c], Double.valueOf(change.before[c]), Double.valueOf(change.after[c]),
                    Double.valueOf(change.ratio[c]));
        }
        System.out.println();
        System.out.println("  No step either way. The change of source is real but it is in the");
        System.out.println("  resolution, not in the level: the last 254 readings sit exactly on");
        System.out.println("  the Fahrenheit grid and the earlier ones do not quite.");
        System.out.println();
        System.out.println("  And in a year they will be gone. The archive replaces its own last");
        System.out.println("  twelve months as the authoritative source catches up, so downloading");
        System.out.println("  this file again returns different numbers for September. That is why");
        System.out.println("  the demo carries its own copy of the table rather than fetching one.");

        title("6. putting back a reading that is not there");
        System.out.printf(L, "  Hide the minimum for %d days from %s and let the model fill it in%n",
                Integer.valueOf(BLANK_LENGTH), Datasets.label(BLANK_FROM));
        System.out.println("  from the maximum and from the days on either side:");
        System.out.println();
        int[] lengths = { BLANK_LENGTH, 1 };
        for (int i = 0; i < lengths.length; ++i) {
            Blanked gap = blank(working, BLANK_FROM, lengths[i]);
            System.out.printf(L, "    %2d %-4s hidden: off by %.4f degrees%n", Integer.valueOf(lengths[i]),
                    lengths[i] == 1 ? "day" : "days", Double.valueOf(gap.plainError));
        }
        System.out.println();
        System.out.println("  One caveat, and it is the only piece of theory on this page. A filter");
        System.out.println("  handed part of a day's readings is exactly right about the thing it");
        System.out.println("  is tracking. It is not quite right about the reading it never saw:");
        System.out.println("  two readings taken the same day can share an error, and the filter");
        System.out.println("  has nothing to work with for the row that is missing, so it cannot");
        System.out.println("  use that. Here the fit says the two errors are unrelated, so nothing");
        System.out.printf(L, "  is lost at all. A model that said otherwise would lose %.4f degrees.%n",
                Double.valueOf(blank(rungs[UNBOUNDED], BLANK_FROM, BLANK_LENGTH).largestTerm));

        title("7. is anything left over");
        System.out.println("  Fit the season by ordinary least squares and the leftovers are still");
        System.out.println("  full of pattern: a warm day follows a warm day. A model that leaves");
        System.out.println("  that in its residuals is claiming more certainty than it has. What");
        System.out.println("  the state space model leaves over should be pattern free, and is:");
        System.out.println();
        System.out.printf(L, "    %-26s %8s %8s%n", "", "acf(1)", "outside");
        for (int c = 0; c < 2; ++c) {
            System.out.printf(L, "    least squares, %-11s %8.4f %8d%n", readings[c],
                    Double.valueOf(white.deterministicAcf1[c]),
                    Integer.valueOf(white.deterministicOutside[c]));
        }
        for (int c = 0; c < 2; ++c) {
            System.out.printf(L, "    the model,     %-11s %8.4f %8d%n", readings[c],
                    Double.valueOf(white.innovationAcf1[c]),
                    Integer.valueOf(white.innovationOutside[c]));
        }
        System.out.printf(L, "%n  How many of the first 20 lags leave the band that two standard%n");
        System.out.printf(L, "  errors, %.4f, would draw around zero.%n", Double.valueOf(white.bound));
        System.out.println();
        double[] acf = ACF.acf(deterministicResiduals(), 1);
        System.out.printf(L, "  The minimum comes out clean. The maximum does not quite, which is%n");
        System.out.println("  the asymmetry of step 1 again, from the other side.");
        System.out.println();
        System.out.printf(L, "  One number worth its own line: the leftovers say %.4f of today%n",
                Double.valueOf(acf[1]));
        System.out.printf(L, "  carries into tomorrow, and the model says %.4f. The model has to say%n",
                Double.valueOf(working.phi()));
        System.out.println("  more, because part of what the leftovers show is measurement error,");
        System.out.println("  and measurement error has no memory. Reading persistence straight off");
        System.out.println("  residuals always understates it. Separating those two is the reason a");
        System.out.println("  model like this exists.");

        title("8. how the model was chosen, and one uncomfortable thing");
        System.out.println("  Nothing above was decided by hand. Four models were fitted, each one");
        System.out.println("  the previous with something added, and the fit itself decided. Higher");
        System.out.println("  log L is better; AIC charges two for every parameter, so it is what a");
        System.out.println("  model has to beat to earn its extra number.");
        System.out.println();
        System.out.printf(L, "  %-32s %s %10s %9s %10s%n", "model", "k", "log L", "gain", "AIC");
        String[] names = { "  the season only", "  + linked errors", "  + weather that persists",
                "  + linked errors again" };
        for (int i = 0; i < RUNGS; ++i) {
            System.out.printf(L, "  %-32s %d %10.2f %9.2f %10.2f%n", names[i],
                    Integer.valueOf(parameterCount(i)), Double.valueOf(rungs[i].logLikelihood),
                    Double.valueOf(rungs[i].logLikelihood - rungs[i == 0 ? 0 : i - 1].logLikelihood),
                    Double.valueOf(rungs[i].aic()));
        }
        System.out.println();
        System.out.println("  Steps 1 to 7 all use the third model. The fourth is the uncomfortable");
        System.out.println("  one, and it is the same addition as the second: let the errors of the");
        System.out.printf(L, "  two readings be linked. Added to the first model that is worth %.2f;%n",
                Double.valueOf(rungs[1].logLikelihood - rungs[0].logLikelihood));
        System.out.printf(L, "  added to the third it is worth %.2f, and the link comes out%n",
                Double.valueOf(rungs[UNBOUNDED].logLikelihood - rungs[WORKING].logLikelihood));
        System.out.printf(L, "  %+.2f where before it was %+.2f.%n",
                Double.valueOf(rungs[UNBOUNDED].noiseCorrelation()),
                Double.valueOf(rungs[1].noiseCorrelation()));
        System.out.println();
        System.out.println("  The reason is that in the first model nothing else could explain the");
        System.out.println("  weather, so the link stood in for it. Once the third model tracks the");
        System.out.println("  weather properly, there is nothing left for the link to do, and it");
        System.out.println("  runs to -1, which is as far as a correlation can go. A number pinned");
        System.out.println("  against the end of its range is not an estimate, it is the fit saying");
        System.out.printf(L, "  it has run out of room. AIC prefers that model by %.2f anyway, and%n",
                Double.valueOf(rungs[WORKING].aic() - rungs[UNBOUNDED].aic()));
        System.out.println("  this demo declines it. That is the whole lesson of this step: what a");
        System.out.println("  parameter means depends on what else is in the model beside it.");
    }

    private static double mean(double[] v) {
        double sum = 0.0;
        for (int i = 0; i < v.length; ++i) {
            sum += v[i];
        }
        return sum / v.length;
    }

    private static double standardDeviation(double[] v) {
        double m = mean(v);
        double sum = 0.0;
        for (int i = 0; i < v.length; ++i) {
            sum += (v[i] - m) * (v[i] - m);
        }
        return Math.sqrt(sum / (v.length - 1));
    }

    private WeatherDemo() {
        throw new AssertionError();
    }
}
