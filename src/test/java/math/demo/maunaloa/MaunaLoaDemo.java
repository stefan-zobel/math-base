package math.demo.maunaloa;

import java.util.Locale;

import math.fit.AkimaInterpolator;
import math.fit.CubicSpline;
import math.fit.KrugerInterpolator;
import math.fit.SplineInterpolator;
import math.fft.ComplexArray;
import math.fft.Fourier;
import math.fun.DFunction;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.Wls;
import math.probe.ACF;
import math.probe.Bootstrap;
import math.probe.DoubleStatistics;
import math.probe.PartialACF;
import math.probe.SampleStatistic;
import math.solve.AdaptiveGaussKronrod;
import math.solve.RootFinder;

/**
 * A worked example: sixty-eight years of monthly CO2 measurements from Mauna
 * Loa, taken apart with the pieces this library provides.
 * <p>
 * The question is the obvious one -- how fast is CO2 rising, how sure can we
 * be, and when does the trend reach 450 ppm -- and answering it needs eight
 * packages, each doing the one thing it is for:
 * <ol>
 * <li>{@code math.probe} to describe the series,</li>
 * <li>{@code math.fit} to put back months that were removed,</li>
 * <li>{@code math.linalg} to fit a trend and a seasonal cycle,</li>
 * <li>{@code math.linalg} again to weight the fit by the measurement
 * uncertainty,</li>
 * <li>{@code math.probe} to show that the residuals are not independent, so
 * the standard errors of step 3 are too small,</li>
 * <li>{@code math.probe} once more for an interval that does not assume
 * they are,</li>
 * <li>{@code math.solve} and {@code math.fit} to integrate the trend two ways
 * and check them against each other,</li>
 * <li>{@code math.solve} and {@code math.fft} to find the 450 ppm crossing and
 * to confirm the seasonal cycle from the frequency side.</li>
 * </ol>
 * Run {@link #main(String[])} to see it. Every step is a package-private
 * method returning a small value type, so {@code MaunaLoaDemoTest} can assert
 * what this prints.
 *
 * @see Datasets
 */
public final class MaunaLoaDemo {

    private static final Locale L = Locale.ROOT;
    private static final double TWO_PI = 2.0 * Math.PI;

    /** Seed for the bootstrap, so that the printed interval is reproducible. */
    static final long BOOTSTRAP_SEED = 20260821L;

    /** Number of terms in the design matrix. */
    private static final int TERMS = 7;

    /** Number of recent years the growth rate is taken over. */
    static final int GROWTH_YEARS = 30;

    // ---------------------------------------------------------------- data

    /** Years elapsed since the first observation. Centred before it is squared. */
    static double[] centredTime() {
        double[] t = Datasets.decimalDate();
        double t0 = t[0];
        double[] tc = new double[t.length];
        for (int i = 0; i < t.length; ++i) {
            tc[i] = t[i] - t0;
        }
        return tc;
    }

    // ------------------------------------------------- 1. describe

    /** What the series looks like before anything is fitted to it. */
    static final class Description {
        final double first;
        final double last;
        final double min;
        final double max;
        final double meanAnnualIncrement;
        final int months;

        Description(double first, double last, double min, double max, double meanAnnualIncrement, int months) {
            this.first = first;
            this.last = last;
            this.min = min;
            this.max = max;
            this.meanAnnualIncrement = meanAnnualIncrement;
            this.months = months;
        }
    }

    static Description describe() {
        double[] y = Datasets.co2();
        DoubleStatistics stats = DoubleStatistics.newInstance();
        for (int i = 0; i < y.length; ++i) {
            stats.accept(y[i]);
        }
        double[] increments = annualIncrements();
        DoubleStatistics inc = DoubleStatistics.newInstance();
        for (int i = 0; i < increments.length; ++i) {
            inc.accept(increments[i]);
        }
        return new Description(y[0], y[y.length - 1], stats.getMin(), stats.getMax(), inc.getAverage(), y.length);
    }

    /** Mean CO2 of every calendar year that is complete in the record. */
    static double[] annualMeans() {
        double[] y = Datasets.co2();
        int firstFull = 12 - (Datasets.FIRST_MONTH - 1);
        int years = (y.length - firstFull) / 12;
        double[] means = new double[years];
        for (int k = 0; k < years; ++k) {
            double total = 0.0;
            for (int m = 0; m < 12; ++m) {
                total += y[firstFull + 12 * k + m];
            }
            means[k] = total / 12.0;
        }
        return means;
    }

    /** Year of the first complete calendar year in the record. */
    static int firstFullYear() {
        return Datasets.FIRST_YEAR + 1;
    }

    /** Year-over-year change of the annual means, in ppm per year. */
    static double[] annualIncrements() {
        double[] means = annualMeans();
        double[] increments = new double[means.length - 1];
        for (int i = 0; i < increments.length; ++i) {
            increments[i] = means[i + 1] - means[i];
        }
        return increments;
    }

    // ------------------------------------------------- 2. hold months out

    /** How well one interpolation scheme puts back a run of missing months. */
    static final class HoldOut {
        final String scheme;
        final int gapLength;
        final double rms;
        final double worst;
        final int cases;

        HoldOut(String scheme, int gapLength, double rms, double worst, int cases) {
            this.scheme = scheme;
            this.gapLength = gapLength;
            this.rms = rms;
            this.worst = worst;
            this.cases = cases;
        }
    }

    /** The interpolation schemes under test, by name. */
    static final String[] SCHEMES = { "natural spline", "Akima classic", "Akima modified", "Kruger" };

    private static CubicSpline interpolate(int scheme, double[] x, double[] y) {
        switch (scheme) {
        case 0:
            return SplineInterpolator.interpolate(x, y);
        case 1:
            return AkimaInterpolator.interpolate(x, y, AkimaInterpolator.Variant.CLASSIC);
        case 2:
            return AkimaInterpolator.interpolate(x, y, AkimaInterpolator.Variant.MODIFIED);
        default:
            return KrugerInterpolator.interpolate(x, y);
        }
    }

    /**
     * Removes a run of months, rebuilds it from what is left, and measures the
     * error. NOAA has filled its own two gaps already, so a gap has to be made
     * before it can be filled: that turns two anecdotes into a measurement over
     * every position the series admits.
     */
    static HoldOut holdOut(int scheme, int gapLength) {
        double[] t = Datasets.decimalDate();
        double[] y = Datasets.co2();
        int n = t.length;
        double sumSquared = 0.0;
        double worst = 0.0;
        int cases = 0;
        double[] xKept = new double[n - gapLength];
        double[] yKept = new double[n - gapLength];
        for (int start = 1; start + gapLength <= n - 1; ++start) {
            int k = 0;
            for (int i = 0; i < n; ++i) {
                if (i < start || i >= start + gapLength) {
                    xKept[k] = t[i];
                    yKept[k] = y[i];
                    ++k;
                }
            }
            CubicSpline spline = interpolate(scheme, xKept, yKept);
            for (int i = start; i < start + gapLength; ++i) {
                double error = spline.value(t[i]) - y[i];
                sumSquared += error * error;
                worst = Math.max(worst, Math.abs(error));
                ++cases;
            }
        }
        return new HoldOut(SCHEMES[scheme], gapLength, Math.sqrt(sumSquared / cases), worst, cases);
    }

    /** Our reconstruction of the two months NOAA interpolated, against theirs. */
    static double[] versusNoaa(int scheme) {
        double[] t = Datasets.decimalDate();
        double[] y = Datasets.co2();
        int[] months = Datasets.noaaInterpolated();
        double[] errors = new double[months.length];
        int n = t.length;
        double[] xKept = new double[n - 1];
        double[] yKept = new double[n - 1];
        for (int m = 0; m < months.length; ++m) {
            int drop = months[m];
            int k = 0;
            for (int i = 0; i < n; ++i) {
                if (i != drop) {
                    xKept[k] = t[i];
                    yKept[k] = y[i];
                    ++k;
                }
            }
            errors[m] = interpolate(scheme, xKept, yKept).value(t[drop]) - y[drop];
        }
        return errors;
    }

    // ------------------------------------------------- 3. trend and season

    /** A fitted trend plus an annual and a semi-annual harmonic. */
    static final class Fit {
        final double[] beta;
        final double[] standardErrors;
        final double[] tValues;
        final double rSquared;
        final double[] residuals;
        final double annualAmplitude;
        final double semiAnnualAmplitude;
        final double seasonalPeakToPeak;
        final int observations;

        Fit(double[] beta, double[] standardErrors, double[] tValues, double rSquared, double[] residuals,
                int observations) {
            this.beta = beta;
            this.standardErrors = standardErrors;
            this.tValues = tValues;
            this.rSquared = rSquared;
            this.residuals = residuals;
            this.observations = observations;
            this.annualAmplitude = Math.hypot(beta[3], beta[4]);
            this.semiAnnualAmplitude = Math.hypot(beta[5], beta[6]);
            double lo = Double.POSITIVE_INFINITY;
            double hi = Double.NEGATIVE_INFINITY;
            for (int k = 0; k < 2000; ++k) {
                double s = seasonal(k / 2000.0);
                lo = Math.min(lo, s);
                hi = Math.max(hi, s);
            }
            this.seasonalPeakToPeak = hi - lo;
        }

        /** The smooth part: quadratic in years since 1958. */
        double trend(double tc) {
            return beta[0] + beta[1] * tc + beta[2] * tc * tc;
        }

        /** The cyclic part, which averages to zero over a year. */
        double seasonal(double tc) {
            return beta[3] * Math.sin(TWO_PI * tc) + beta[4] * Math.cos(TWO_PI * tc)
                    + beta[5] * Math.sin(2.0 * TWO_PI * tc) + beta[6] * Math.cos(2.0 * TWO_PI * tc);
        }

        /** Rate of rise of the trend at that time, in ppm per year. */
        double growth(double tc) {
            return beta[1] + 2.0 * beta[2] * tc;
        }
    }

    /**
     * The design matrix. Time is measured from the first observation and not
     * from year zero: an uncentred abscissa around 2000 makes the squared
     * column of order 4e6 and ruins the conditioning of a fit that is
     * otherwise trivial.
     */
    static DMatrix design(double[] tc) {
        DMatrix X = new DMatrix(tc.length, TERMS);
        for (int i = 0; i < tc.length; ++i) {
            double t = tc[i];
            X.set(i, 0, 1.0);
            X.set(i, 1, t);
            X.set(i, 2, t * t);
            X.set(i, 3, Math.sin(TWO_PI * t));
            X.set(i, 4, Math.cos(TWO_PI * t));
            X.set(i, 5, Math.sin(2.0 * TWO_PI * t));
            X.set(i, 6, Math.cos(2.0 * TWO_PI * t));
        }
        return X;
    }

    private static DMatrix column(double[] v) {
        DMatrix y = new DMatrix(v.length, 1);
        for (int i = 0; i < v.length; ++i) {
            y.set(i, 0, v[i]);
        }
        return y;
    }

    private static Fit toFit(LSSummary summary, int observations) {
        return new Fit(summary.getBeta().toArray(), summary.getCoefficientStandardErrors().toArray(),
                summary.getTValues().toArray(), summary.getRSquared(), summary.getResiduals().toArray(), observations);
    }

    static Fit fitOrdinary() {
        double[] tc = centredTime();
        return toFit(OLS.estimate(0.05, design(tc), column(Datasets.co2())), tc.length);
    }

    /** Indices of the months that carry a usable measurement uncertainty. */
    static int[] weightedMonths() {
        int count = 0;
        for (int i = 0; i < Datasets.size(); ++i) {
            if (Datasets.hasUncertainty(i)) {
                ++count;
            }
        }
        int[] idx = new int[count];
        int k = 0;
        for (int i = 0; i < Datasets.size(); ++i) {
            if (Datasets.hasUncertainty(i)) {
                idx[k++] = i;
            }
        }
        return idx;
    }

    /** The same fit on the months that have an uncertainty, weighted or not. */
    static Fit fitOnUncertainMonths(boolean weighted) {
        int[] idx = weightedMonths();
        double[] tc = centredTime();
        double[] y = Datasets.co2();
        double[] unc = Datasets.uncertainty();
        double[] tSub = new double[idx.length];
        double[] ySub = new double[idx.length];
        double[] w = new double[idx.length];
        for (int k = 0; k < idx.length; ++k) {
            tSub[k] = tc[idx[k]];
            ySub[k] = y[idx[k]];
            w[k] = 1.0 / (unc[idx[k]] * unc[idx[k]]);
        }
        DMatrix X = design(tSub);
        DMatrix yv = column(ySub);
        LSSummary summary = weighted ? Wls.estimate(0.05, X, yv, w) : OLS.estimate(0.05, X, yv);
        return toFit(summary, idx.length);
    }

    /** Root mean square gap between our seasonal adjustment and NOAA's. */
    static double versusNoaaDeseasonalized(Fit fit) {
        double[] tc = centredTime();
        double[] y = Datasets.co2();
        double[] noaa = Datasets.deseasonalized();
        double sumSquared = 0.0;
        for (int i = 0; i < y.length; ++i) {
            double ours = y[i] - fit.seasonal(tc[i]);
            double d = ours - noaa[i];
            sumSquared += d * d;
        }
        return Math.sqrt(sumSquared / y.length);
    }

    // ------------------------------------------------- 5. residual structure

    /** Autocorrelation of the residuals, and of the annual increments. */
    static final class Correlation {
        final double[] acf;
        final double[] pacf;
        final double band;
        final double incrementLag1;
        final double incrementBand;
        final double allIncrementsLag1;

        Correlation(double[] acf, double[] pacf, double band, double incrementLag1, double incrementBand,
                double allIncrementsLag1) {
            this.acf = acf;
            this.pacf = pacf;
            this.band = band;
            this.incrementLag1 = incrementLag1;
            this.incrementBand = incrementBand;
            this.allIncrementsLag1 = allIncrementsLag1;
        }
    }

    static Correlation correlate(Fit fit) {
        int lags = 36;
        double[] acf = ACF.acf(fit.residuals, lags);
        double[] pacf = PartialACF.partialAutocorrelation(fit.residuals, lags);
        double band = PartialACF.getConfidenceInterval(fit.residuals.length)[1];
        // the window that is actually resampled in step 6, not the whole
        // record: over all 66 increments the lag 1 value mostly measures the
        // long run rise of the increments themselves
        double[] recent = recentIncrements(GROWTH_YEARS);
        return new Correlation(acf, pacf, band, ACF.acf(recent, 1)[1],
                PartialACF.getConfidenceInterval(GROWTH_YEARS)[1], ACF.acf(annualIncrements(), 1)[1]);
    }

    // ------------------------------------------------- 6. growth rate

    /** An interval for the recent growth rate that assumes no model. */
    static final class Growth {
        final double observed;
        final double[] percentile;
        final double[] bca;
        final int years;
        final double classicalHalfWidth;

        Growth(double observed, double[] percentile, double[] bca, int years, double classicalHalfWidth) {
            this.observed = observed;
            this.percentile = percentile;
            this.bca = bca;
            this.years = years;
            this.classicalHalfWidth = classicalHalfWidth;
        }
    }

    private static final SampleStatistic MEAN = new SampleStatistic() {
        @Override
        public double apply(double[] sample) {
            double total = 0.0;
            for (int i = 0; i < sample.length; ++i) {
                total += sample[i];
            }
            return total / sample.length;
        }
    };

    /** The most recent annual increments, over which the trend is nearly flat. */
    static double[] recentIncrements(int years) {
        double[] all = annualIncrements();
        double[] recent = new double[years];
        System.arraycopy(all, all.length - years, recent, 0, years);
        return recent;
    }

    static Growth growth(int years) {
        double[] recent = recentIncrements(years);
        Bootstrap bs = new Bootstrap(recent, MEAN, 20000, BOOTSTRAP_SEED);
        double mean = MEAN.apply(recent);
        double sumSquared = 0.0;
        for (int i = 0; i < recent.length; ++i) {
            double d = recent[i] - mean;
            sumSquared += d * d;
        }
        double standardError = Math.sqrt(sumSquared / (recent.length - 1)) / Math.sqrt(recent.length);
        return new Growth(mean, bs.getConfidenceInterval(0.95), bs.getConfidenceIntervalBCa(0.95), years,
                1.96 * standardError);
    }

    // ------------------------------------------------- 7. cumulative excess

    /** The area between the seasonally adjusted record and its 1958 level. */
    static final class Excess {
        final double closedForm;
        final double quadrature;
        final double relativeDifference;

        Excess(double closedForm, double quadrature) {
            this.closedForm = closedForm;
            this.quadrature = quadrature;
            this.relativeDifference = Math.abs(closedForm - quadrature) / Math.abs(closedForm);
        }
    }

    /** A spline through the seasonally adjusted observations. */
    static CubicSpline deseasonalizedSpline(Fit fit) {
        double[] tc = centredTime();
        double[] y = Datasets.co2();
        double[] adjusted = new double[y.length];
        for (int i = 0; i < y.length; ++i) {
            adjusted[i] = y[i] - fit.seasonal(tc[i]);
        }
        return SplineInterpolator.interpolate(tc, adjusted);
    }

    /**
     * Integrates that spline twice over: once by the closed form the spline
     * carries, once by summing Gauss-Kronrod over its segments. K15 is exact
     * for a cubic on one interval, so the two must agree to round-off, and
     * neither knows anything about the other.
     */
    static Excess excess(Fit fit) {
        final CubicSpline spline = deseasonalizedSpline(fit);
        final double baseline = spline.value(spline.lowerBound());
        double[] knots = spline.knots();
        double span = spline.upperBound() - spline.lowerBound();
        double closed = spline.integrate(spline.lowerBound(), spline.upperBound()) - baseline * span;
        DFunction shifted = new DFunction() {
            @Override
            public double apply(double x) {
                return spline.value(x) - baseline;
            }
        };
        double quadrature = 0.0;
        for (int j = 0; j + 1 < knots.length; ++j) {
            quadrature += AdaptiveGaussKronrod.integrate1D(AdaptiveGaussKronrod.G7_K15.POINTS_15, shifted, knots[j],
                    knots[j + 1]).value;
        }
        return new Excess(closed, quadrature);
    }

    // ------------------------------------------------- 8. crossing and spectrum

    /** Where the fitted trend reaches a given level. */
    static final class Crossing {
        final double level;
        final double year;
        final double residual;

        Crossing(double level, double year, double residual) {
            this.level = level;
            this.year = year;
            this.residual = residual;
        }
    }

    static Crossing crossing(final Fit fit, final double level) {
        double t0 = Datasets.decimalDate()[0];
        DFunction f = new DFunction() {
            @Override
            public double apply(double tc) {
                return fit.trend(tc) - level;
            }
        };
        double root = RootFinder.brentDekker(0.0, 200.0, f, 1.0e-12);
        return new Crossing(level, t0 + root, f.apply(root));
    }

    /** The strongest cycles left once the trend is taken out. */
    static final class Spectrum {
        final double[] periods;
        final double[] power;
        final int months;

        Spectrum(double[] periods, double[] power, int months) {
            this.periods = periods;
            this.power = power;
            this.months = months;
        }
    }

    /**
     * A periodogram of the detrended series. The length is trimmed to a whole
     * number of years so that a twelve month cycle lands exactly on a bin
     * instead of straddling two; 816 is not a power of two, so this goes
     * through the Bluestein path.
     */
    static Spectrum spectrum(Fit fit) {
        double[] tc = centredTime();
        double[] y = Datasets.co2();
        int months = (y.length / 12) * 12;
        double[] detrended = new double[months];
        for (int i = 0; i < months; ++i) {
            detrended[i] = y[i] - fit.trend(tc[i]);
        }
        double mean = 0.0;
        for (int i = 0; i < months; ++i) {
            mean += detrended[i];
        }
        mean /= months;
        for (int i = 0; i < months; ++i) {
            detrended[i] -= mean;
        }
        ComplexArray spectrumOf = Fourier.forwardDFT(detrended);
        double[] squared = spectrumOf.absSquaredScaled();
        int half = months / 2;
        double[] periods = new double[half];
        double[] power = new double[half];
        for (int k = 1; k < half; ++k) {
            periods[k] = (double) months / k;
            power[k] = squared[k];
        }
        return new Spectrum(periods, power, months);
    }

    /** Index of the n-th strongest line in a spectrum, skipping the DC bin. */
    static int strongest(Spectrum s, int rank) {
        boolean[] taken = new boolean[s.power.length];
        int best = -1;
        for (int r = 0; r <= rank; ++r) {
            best = -1;
            for (int k = 1; k < s.power.length; ++k) {
                if (!taken[k] && (best < 0 || s.power[k] > s.power[best])) {
                    best = k;
                }
            }
            taken[best] = true;
        }
        return best;
    }

    // ---------------------------------------------------------------- output

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    /** A crude line plot, because a column of numbers hides the shape. */
    private static void plot(double[] x, double[] y, int width, int height) {
        double xLo = x[0];
        double xHi = x[x.length - 1];
        double yLo = Double.POSITIVE_INFINITY;
        double yHi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; ++i) {
            yLo = Math.min(yLo, y[i]);
            yHi = Math.max(yHi, y[i]);
        }
        int[] row = new int[width];
        double[] sum = new double[width];
        int[] count = new int[width];
        for (int i = 0; i < y.length; ++i) {
            int col = (int) ((x[i] - xLo) / (xHi - xLo) * (width - 1));
            sum[col] += y[i];
            count[col]++;
        }
        for (int c = 0; c < width; ++c) {
            double v = count[c] > 0 ? sum[c] / count[c] : Double.NaN;
            row[c] = count[c] > 0 ? (int) ((v - yLo) / (yHi - yLo) * (height - 1)) : -1;
        }
        for (int r = height - 1; r >= 0; --r) {
            StringBuilder line = new StringBuilder();
            line.append(String.format(L, "%8.1f |", Double.valueOf(yLo + (yHi - yLo) * r / (height - 1))));
            for (int c = 0; c < width; ++c) {
                line.append(row[c] == r ? '*' : ' ');
            }
            System.out.println(line);
        }
        StringBuilder axis = new StringBuilder("         +");
        for (int c = 0; c < width; ++c) {
            axis.append('-');
        }
        System.out.println(axis);
        System.out.println(String.format(L, "          %-" + (width - 6) + ".0f%6.0f", Double.valueOf(xLo),
                Double.valueOf(xHi)));
    }

    private static void correlogram(String title, double[] r, double band, int lags) {
        System.out.println("  " + title + "   (band +/-" + String.format(L, "%.3f", Double.valueOf(band)) + ")");
        for (int k = 1; k <= lags; ++k) {
            StringBuilder bar = new StringBuilder();
            int len = (int) Math.round(Math.abs(r[k]) * 40.0);
            for (int i = 0; i < len; ++i) {
                bar.append('#');
            }
            System.out.println(String.format(L, "   lag %2d %s%-40s %s%.3f", Integer.valueOf(k), r[k] < 0 ? "-" : "+",
                    bar.toString(), Math.abs(r[k]) > band ? "*" : " ", Double.valueOf(r[k])));
        }
    }

    public static void main(String[] args) {
        System.out.println("Mauna Loa CO2, " + Datasets.size() + " monthly means from " + Datasets.label(0) + " to "
                + Datasets.label(Datasets.size() - 1));
        System.out.println("Source: NOAA Global Monitoring Laboratory, retrieved 2026-08-21");

        rule("1. the record");
        Description d = describe();
        System.out.println(String.format(L, "  %d months, from %.2f to %.2f ppm (min %.2f, max %.2f)",
                Integer.valueOf(d.months), Double.valueOf(d.first), Double.valueOf(d.last), Double.valueOf(d.min),
                Double.valueOf(d.max)));
        System.out.println(String.format(L, "  mean annual increment over the whole record: %.3f ppm/yr",
                Double.valueOf(d.meanAnnualIncrement)));
        System.out.println();
        plot(Datasets.decimalDate(), Datasets.co2(), 64, 16);

        rule("2. take months out, put them back  (math.fit)");
        System.out.println("  NOAA has already filled its own two gaps, so a gap has to be made before");
        System.out.println("  it can be filled. Every admissible position is used; the error is in ppm.");
        System.out.println();
        int[] gaps = { 1, 3, 6, 12 };
        System.out.println(String.format(L, "  %-16s %10s %10s %10s %10s", "scheme", "gap 1", "gap 3", "gap 6",
                "gap 12"));
        for (int s = 0; s < SCHEMES.length; ++s) {
            StringBuilder line = new StringBuilder(String.format(L, "  %-16s", SCHEMES[s]));
            for (int g = 0; g < gaps.length; ++g) {
                line.append(String.format(L, " %10.4f", Double.valueOf(holdOut(s, gaps[g]).rms)));
            }
            System.out.println(line + "   RMS");
        }
        System.out.println();
        int[] noaaMonths = Datasets.noaaInterpolated();
        System.out.println("  the two months NOAA interpolated itself, ours minus theirs, in ppm:");
        for (int s = 0; s < SCHEMES.length; ++s) {
            double[] e = versusNoaa(s);
            StringBuilder line = new StringBuilder(String.format(L, "  %-16s", SCHEMES[s]));
            for (int m = 0; m < e.length; ++m) {
                line.append(String.format(L, " %s %+8.4f", Datasets.label(noaaMonths[m]), Double.valueOf(e[m])));
            }
            System.out.println(line);
        }

        rule("3. a trend and a seasonal cycle  (math.linalg)");
        Fit fit = fitOrdinary();
        String[] names = { "intercept", "t", "t^2", "sin 2pi t", "cos 2pi t", "sin 4pi t", "cos 4pi t" };
        System.out.println(String.format(L, "  %-12s %12s %12s %10s", "term", "estimate", "std error", "t"));
        for (int j = 0; j < TERMS; ++j) {
            System.out.println(String.format(L, "  %-12s %12.6f %12.6f %10.1f", names[j], Double.valueOf(fit.beta[j]),
                    Double.valueOf(fit.standardErrors[j]), Double.valueOf(fit.tValues[j])));
        }
        System.out.println(String.format(L, "  R^2 = %.6f on %d observations", Double.valueOf(fit.rSquared),
                Integer.valueOf(fit.observations)));
        System.out.println(String.format(L, "  seasonal cycle: annual amplitude %.3f ppm, semi-annual %.3f ppm,",
                Double.valueOf(fit.annualAmplitude), Double.valueOf(fit.semiAnnualAmplitude)));
        System.out.println(String.format(L, "                  peak to peak %.3f ppm", Double.valueOf(
                fit.seasonalPeakToPeak)));
        double[] tc = centredTime();
        System.out.println(String.format(L, "  growth of the trend: %.3f ppm/yr in 1958, %.3f ppm/yr now",
                Double.valueOf(fit.growth(0.0)), Double.valueOf(fit.growth(tc[tc.length - 1]))));
        System.out.println(String.format(L,
                "  our seasonal adjustment against NOAA's own: %.4f ppm RMS over %d months",
                Double.valueOf(versusNoaaDeseasonalized(fit)), Integer.valueOf(Datasets.size())));

        rule("4. weighting by the measurement uncertainty  (math.linalg.Wls)");
        int[] idx = weightedMonths();
        System.out.println("  Only " + idx.length + " of " + Datasets.size() + " months carry an uncertainty; the");
        System.out.println("  Scripps era before 1974-05 has none. Both fits use the same " + idx.length + " months,");
        System.out.println("  so the only difference is the weighting.");
        Fit unweighted = fitOnUncertainMonths(false);
        Fit weighted = fitOnUncertainMonths(true);
        System.out.println();
        System.out.println(String.format(L, "  %-12s %14s %14s %12s", "term", "OLS", "WLS", "shift/se"));
        for (int j = 0; j < TERMS; ++j) {
            double shift = (weighted.beta[j] - unweighted.beta[j]) / unweighted.standardErrors[j];
            System.out.println(String.format(L, "  %-12s %14.6f %14.6f %12.2f", names[j],
                    Double.valueOf(unweighted.beta[j]), Double.valueOf(weighted.beta[j]), Double.valueOf(shift)));
        }

        rule("5. the residuals are not independent  (math.probe)");
        Correlation c = correlate(fit);
        System.out.println("  If they were, the standard errors printed in step 3 would mean what");
        System.out.println("  they say. Lag 1 is " + String.format(L, "%.3f", Double.valueOf(c.acf[1]))
                + ", so they do not.");
        System.out.println();
        correlogram("autocorrelation of the residuals", c.acf, c.band, 24);
        System.out.println();
        System.out.println("  partial autocorrelation, first six lags:");
        for (int k = 1; k <= 6; ++k) {
            System.out.println(String.format(L, "   lag %2d %8.3f %s", Integer.valueOf(k), Double.valueOf(c.pacf[k]),
                    Math.abs(c.pacf[k]) > c.band ? "significant" : ""));
        }

        rule("6. an interval that assumes no model  (math.probe.Bootstrap)");
        Growth g = growth(GROWTH_YEARS);
        System.out.println("  Step 5 rules out resampling the monthly residuals, so the statistic is");
        System.out.println("  the mean of the last " + GROWTH_YEARS + " annual increments instead. Their lag 1");
        System.out.println(String.format(L, "  autocorrelation is %.3f against a band of %.3f, so treating them as",
                Double.valueOf(c.incrementLag1), Double.valueOf(c.incrementBand)));
        System.out.println(String.format(L,
                "  independent holds up. Over all %d increments it would be %.3f, but that",
                Integer.valueOf(annualIncrements().length), Double.valueOf(c.allIncrementsLag1)));
        System.out.println("  number is the long run rise of the increments, not dependence.");
        System.out.println();
        System.out.println(String.format(L, "  observed mean increment  %.4f ppm/yr", Double.valueOf(g.observed)));
        System.out.println(String.format(L, "  percentile 95%%           [%.4f, %.4f]", Double.valueOf(g.percentile[0]),
                Double.valueOf(g.percentile[1])));
        System.out.println(String.format(L, "  BCa 95%%                  [%.4f, %.4f]", Double.valueOf(g.bca[0]),
                Double.valueOf(g.bca[1])));
        System.out.println(String.format(L, "  classical 1.96 se        [%.4f, %.4f]",
                Double.valueOf(g.observed - g.classicalHalfWidth), Double.valueOf(g.observed + g.classicalHalfWidth)));
        System.out.println("  (seed " + BOOTSTRAP_SEED + ", so these lines are the same on every run and machine)");
        System.out.println();
        System.out.println("  The three intervals agree, and that agreement is the point: on a");
        System.out.println("  statistic whose dependence has been checked, the classical formula and");
        System.out.println("  the bootstrap say the same thing. The standard errors in step 3 are the");
        System.out.println("  ones with nothing behind them.");

        rule("7. cumulative excess, integrated two ways  (math.fit and math.solve)");
        Excess e = excess(fit);
        System.out.println("  Area between the seasonally adjusted record and its 1958 level.");
        System.out.println(String.format(L, "  closed form from the spline   %18.8f ppm*yr",
                Double.valueOf(e.closedForm)));
        System.out.println(String.format(L, "  Gauss-Kronrod over the same   %18.8f ppm*yr",
                Double.valueOf(e.quadrature)));
        System.out.println(String.format(L, "  relative difference           %18.2e", Double.valueOf(
                e.relativeDifference)));

        rule("8. when does the trend reach 450 ppm, and is the cycle really annual");
        Crossing x = crossing(fit, 450.0);
        double t0 = Datasets.decimalDate()[0];
        System.out.println(String.format(L, "  Brent-Dekker root: %.2f, where the trend is %.9f ppm",
                Double.valueOf(x.year), Double.valueOf(fit.trend(x.year - t0))));
        System.out.println("  (the residual there is at round-off level, so it is not printed as a");
        System.out.println("  number: it differs in its last digits between the scalar and the");
        System.out.println("  vectorized build, which is round-off and not a result)");
        System.out.println("  (a quadratic extrapolated forward, which is arithmetic and not a forecast)");
        System.out.println();
        Spectrum sp = spectrum(fit);
        System.out.println("  Periodogram of the detrended series over " + sp.months + " months:");
        for (int rank = 0; rank < 3; ++rank) {
            int k = strongest(sp, rank);
            System.out.println(String.format(L, "   %d. strongest line at %6.2f months", Integer.valueOf(rank + 1),
                    Double.valueOf(sp.periods[k])));
        }
        System.out.println("  The harmonics of step 3 were put in by hand at 12 and 6 months; the");
        System.out.println("  transform did not know that and found them anyway.");
        System.out.println();
    }

    private MaunaLoaDemo() {
        throw new AssertionError();
    }
}
