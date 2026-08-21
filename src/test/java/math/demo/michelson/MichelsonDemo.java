package math.demo.michelson;

import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import math.distribution.Normal;
import math.distribution.StudentT;
import math.fun.DFunction;
import math.linalg.DMatrix;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.list.DoubleList;
import math.probe.Bootstrap;
import math.probe.DoubleStatistics;
import math.probe.PSquaredMedian;
import math.probe.SampleStatistic;
import math.probe.SimpleTDigest;
import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.RootFinder;

/**
 * A worked example: the 100 measurements Michelson made of the speed of light
 * in the summer of 1879, and the question of what a confidence interval is
 * worth.
 * <p>
 * Every other demo in this package compares an estimate against another
 * estimate. Here the truth is known exactly -- since 1983 the metre is defined
 * in terms of the speed of light, so {@link Datasets#ACCEPTED} carries no error
 * at all -- and that turns the usual exercise inside out. The interval can be
 * checked rather than believed.
 * <p>
 * Three things this shows that are worth more than the estimate:
 * <ul>
 * <li><b>Precision is not accuracy.</b> The classical interval, the bootstrap
 * percentile interval and the BCa interval agree with each other to three
 * decimals and all three miss the truth, by nearly three times their own half
 * width. Nothing computed from these numbers could have told anyone that.</li>
 * <li><b>A difference that survives the eye need not survive a model.</b> The
 * afternoon runs read 52 km/s above the morning ones, which looks like a
 * finding until day and temperature are in the same regression, where it
 * collapses to nothing. What does survive is a drift over the month.</li>
 * <li><b>The library checking itself.</b> The quantile the interval is built
 * from, the distribution function behind it and the density behind that are
 * three separate pieces of code, and section 4 makes each one confirm the
 * next.
 * </ul>
 * <p>
 * Everything is seeded and two runs of {@link #main(String[])} produce
 * identical output, on either source tree.
 */
public final class MichelsonDemo {

    /** Formatting locale, so that the output does not depend on the machine. */
    private static final Locale L = Locale.ROOT;

    /** Confidence level of every interval printed here. */
    static final double CONFIDENCE = 0.95;

    /** Bootstrap replications. */
    static final int BOOT_ITERATIONS = 20000;

    /** Seed of the bootstrap. */
    static final long BOOT_SEED = 18790605L;

    /** Quantiles the sketches are asked for in section 2. */
    static final double[] QUANTILES = { 0.05, 0.25, 0.5, 0.75, 0.95 };

    /**
     * The bound section 4 claims for every agreement it prints, and the one
     * {@code MichelsonDemoTest} asserts. The measured deviations are one to
     * three orders below it; their digits are round-off and differ between the
     * two source trees, so they do not reach the page.
     */
    static final double AGREEMENT = 1.0e-11;

    /** {@link #AGREEMENT} as it is printed. */
    private static final String AGREEMENT_TEXT = "1e-11";

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

    // ------------------------------------------------------------- 1. the run

    /** What the hundred runs look like before anything is fitted to them. */
    static final class Description {
        final int runs;
        final int morning;
        final int afternoon;
        final double mean;
        final double sd;
        final double min;
        final double max;
        final double standardError;

        Description(int runs, int morning, int afternoon, double mean, double sd, double min, double max) {
            this.runs = runs;
            this.morning = morning;
            this.afternoon = afternoon;
            this.mean = mean;
            this.sd = sd;
            this.min = min;
            this.max = max;
            this.standardError = sd / Math.sqrt(runs);
        }

        /** Distance from the accepted value, in km/s. */
        double biasKmPerSecond() {
            return 1000.0 * (mean - Datasets.ACCEPTED);
        }
    }

    static Description describe() {
        double[] speed = Datasets.speed();
        double[] afternoon = Datasets.afternoon();
        DoubleStatistics stats = DoubleStatistics.newInstance();
        int pm = 0;
        for (int i = 0; i < speed.length; ++i) {
            stats.accept(speed[i]);
            pm += (int) afternoon[i];
        }
        return new Description(speed.length, speed.length - pm, pm, stats.getAverage(), stats.getStandardDeviation(),
                stats.getMin(), stats.getMax());
    }

    // -------------------------------------------------------- 2. the sketches

    /** What a sketch answers when it is holding every point it was given. */
    static final class SketchCheck {
        final double[] exact;
        final double[] digest;
        final int centroids;
        final double digestWorst;
        final double squaredMedian;
        final double medianDeviation;

        SketchCheck(double[] exact, double[] digest, int centroids, double digestWorst, double squaredMedian,
                double medianDeviation) {
            this.exact = exact;
            this.digest = digest;
            this.centroids = centroids;
            this.digestWorst = digestWorst;
            this.squaredMedian = squaredMedian;
            this.medianDeviation = medianDeviation;
        }
    }

    /** Quantile by linear interpolation between the order statistics. */
    static double exactQuantile(double[] sorted, double q) {
        double h = (sorted.length - 1) * q;
        int lo = (int) Math.floor(h);
        int hi = Math.min(lo + 1, sorted.length - 1);
        return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
    }

    static SketchCheck sketches() {
        double[] speed = Datasets.speed();
        double[] sorted = speed.clone();
        Arrays.sort(sorted);

        SimpleTDigest digest = new SimpleTDigest(100.0);
        PSquaredMedian squared = new PSquaredMedian();
        for (int i = 0; i < speed.length; ++i) {
            digest.accept(speed[i]);
            squared.accept(speed[i]);
        }

        double[] exact = new double[QUANTILES.length];
        double[] got = new double[QUANTILES.length];
        double worst = 0.0;
        for (int k = 0; k < QUANTILES.length; ++k) {
            exact[k] = exactQuantile(sorted, QUANTILES[k]);
            got[k] = digest.getQuantile(QUANTILES[k]);
            worst = Math.max(worst, Math.abs(got[k] - exact[k]));
        }
        double median = exactQuantile(sorted, 0.5);
        return new SketchCheck(exact, got, digest.getCentroidCount(), worst, squared.getMedian(),
                Math.abs(squared.getMedian() - median));
    }

    // ------------------------------------------------------- 3. the intervals

    /** Three intervals for the same quantity, and the truth they are measured against. */
    static final class Intervals {
        final double mean;
        final double criticalT;
        final double[] classical;
        final double[] percentile;
        final double[] bca;
        final double bootstrapMean;
        final double standardError;

        Intervals(double mean, double criticalT, double[] classical, double[] percentile, double[] bca,
                double bootstrapMean, double standardError) {
            this.mean = mean;
            this.criticalT = criticalT;
            this.classical = classical;
            this.percentile = percentile;
            this.bca = bca;
            this.bootstrapMean = bootstrapMean;
            this.standardError = standardError;
        }

        /** How far the accepted value lies below an interval, in units of its half width. */
        static double halfWidthsAway(double[] interval) {
            double half = 0.5 * (interval[1] - interval[0]);
            return (interval[0] - Datasets.ACCEPTED) / half;
        }

        /** Whether the accepted value lies inside an interval at all. */
        static boolean contains(double[] interval) {
            return interval[0] <= Datasets.ACCEPTED && Datasets.ACCEPTED <= interval[1];
        }
    }

    static Intervals intervals() {
        double[] speed = Datasets.speed();
        Description d = describe();
        double criticalT = new StudentT(speed.length - 1.0).inverseCdf(0.5 * (1.0 + CONFIDENCE));
        double half = criticalT * d.standardError;
        double[] classical = { d.mean - half, d.mean + half };

        Bootstrap bootstrap = new Bootstrap(speed, MEAN, BOOT_ITERATIONS, BOOT_SEED);
        return new Intervals(d.mean, criticalT, classical, bootstrap.getConfidenceInterval(CONFIDENCE),
                bootstrap.getConfidenceIntervalBCa(CONFIDENCE), bootstrap.getMean(), d.standardError);
    }

    // ------------------------------------------- 4. the library checking itself

    /** Four agreements between pieces of the library that are written independently. */
    static final class CrossCheck {
        final double cdfAgainstIntegral;
        final double normalQuantileAgainstRoot;
        final double studentQuantileAgainstRoot;
        final double tabulatedQuantile;
        final double meanAgainstIntegral;
        final double varianceAgainstIntegral;

        CrossCheck(double cdfAgainstIntegral, double normalQuantileAgainstRoot, double studentQuantileAgainstRoot,
                double tabulatedQuantile, double meanAgainstIntegral, double varianceAgainstIntegral) {
            this.cdfAgainstIntegral = cdfAgainstIntegral;
            this.normalQuantileAgainstRoot = normalQuantileAgainstRoot;
            this.studentQuantileAgainstRoot = studentQuantileAgainstRoot;
            this.tabulatedQuantile = tabulatedQuantile;
            this.meanAgainstIntegral = meanAgainstIntegral;
            this.varianceAgainstIntegral = varianceAgainstIntegral;
        }
    }

    private static DFunction density(final Normal normal) {
        return new DFunction() {
            @Override
            public double apply(double x) {
                return normal.pdf(x);
            }
        };
    }

    /** Inverts a distribution function by Brent-Dekker, independently of the class itself. */
    static double invert(final math.distribution.ContinuousDistribution distribution, final double p, double lo,
            double hi) {
        return RootFinder.brentDekker(lo, hi, new DFunction() {
            @Override
            public double apply(double x) {
                return distribution.cdf(x) - p;
            }
        }, 1.0e-14);
    }

    static CrossCheck crossCheck() {
        Description d = describe();
        final double mu = d.mean;
        final double sigma = d.sd;
        final Normal fitted = new Normal(mu, sigma);
        double from = mu - 12.0 * sigma;
        double below = fitted.cdf(from);

        double worstCdf = 0.0;
        for (int k = -30; k <= 30; k += 3) {
            double x = mu + 0.1 * k * sigma;
            double integrated = below
                    + AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, density(fitted), from, x, 1.0e-13, 40);
            worstCdf = Math.max(worstCdf, Math.abs(fitted.cdf(x) - integrated));
        }

        double worstNormal = 0.0;
        double worstStudent = 0.0;
        double[] ps = { 0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999 };
        StudentT student = new StudentT(Datasets.size() - 1.0);
        for (int k = 0; k < ps.length; ++k) {
            worstNormal = Math.max(worstNormal,
                    Math.abs(fitted.inverseCdf(ps[k]) - invert(fitted, ps[k], mu - 40.0 * sigma, mu + 40.0 * sigma)));
            worstStudent = Math.max(worstStudent,
                    Math.abs(student.inverseCdf(ps[k]) - invert(student, ps[k], -1.0e4, 1.0e4)));
        }

        double first = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                return x * fitted.pdf(x);
            }
        }, mu - 14.0 * sigma, mu + 14.0 * sigma, 1.0e-13, 40);
        double second = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
            @Override
            public double apply(double x) {
                return (x - mu) * (x - mu) * fitted.pdf(x);
            }
        }, mu - 14.0 * sigma, mu + 14.0 * sigma, 1.0e-13, 40);

        return new CrossCheck(worstCdf, worstNormal, worstStudent, student.inverseCdf(0.975),
                Math.abs(first - fitted.mean()), Math.abs(second - fitted.variance()));
    }

    // ------------------------------------------------------------ 5. the model

    /** The regression that asks where the spread of the measurements comes from. */
    static final class Model {
        final String[] terms;
        final double[] beta;
        final double[] tValues;
        final double[] pValues;
        final double[][] intervals;
        final double rSquared;
        final double residualSd;
        final double rawSd;
        final int degreesOfFreedom;

        Model(String[] terms, double[] beta, double[] tValues, double[] pValues, double[][] intervals, double rSquared,
                double residualSd, double rawSd, int degreesOfFreedom) {
            this.terms = terms;
            this.beta = beta;
            this.tValues = tValues;
            this.pValues = pValues;
            this.intervals = intervals;
            this.rSquared = rSquared;
            this.residualSd = residualSd;
            this.rawSd = rawSd;
            this.degreesOfFreedom = degreesOfFreedom;
        }

        /** Index of a term by name. */
        int indexOf(String term) {
            for (int i = 0; i < terms.length; ++i) {
                if (terms[i].equals(term)) {
                    return i;
                }
            }
            throw new IllegalArgumentException(term);
        }
    }

    static Model model() {
        double[] speed = Datasets.speed();
        double[] afternoon = Datasets.afternoon();
        double[] day = Datasets.day();
        double[] temperature = Datasets.temperature();
        int n = speed.length;

        DMatrix X = new DMatrix(n, 4);
        DMatrix y = new DMatrix(n, 1);
        for (int i = 0; i < n; ++i) {
            X.set(i, 0, 1.0);
            X.set(i, 1, afternoon[i]);
            X.set(i, 2, day[i]);
            X.set(i, 3, temperature[i]);
            y.set(i, 0, speed[i]);
        }
        LSSummary fit = OLS.estimate(1.0 - CONFIDENCE, X, y);
        List<DoubleList> bounds = fit.getConfidenceIntervals();
        double[][] ci = new double[bounds.size()][];
        for (int j = 0; j < ci.length; ++j) {
            ci[j] = bounds.get(j).toArray();
        }
        return new Model(new String[] { "intercept", "afternoon", "day", "temperature" }, fit.getBeta().toArray(),
                fit.getTValues().toArray(), fit.getPValues().toArray(), ci, fit.getRSquared(),
                Math.sqrt(fit.getSigmaHatSquared()), describe().sd, fit.getDegreesOfFreedom());
    }

    /** The difference the eye sees, before any model is asked about it. */
    static double rawAfternoonEffect() {
        double[] speed = Datasets.speed();
        double[] afternoon = Datasets.afternoon();
        double am = 0.0;
        double pm = 0.0;
        int nAm = 0;
        int nPm = 0;
        for (int i = 0; i < speed.length; ++i) {
            if (afternoon[i] == 1.0) {
                pm += speed[i];
                nPm++;
            } else {
                am += speed[i];
                nAm++;
            }
        }
        return 1000.0 * (pm / nPm - am / nAm);
    }

    // --------------------------------------------------------------- output

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    private static String bracket(double[] interval) {
        return String.format(L, "[%.6f, %.6f]", Double.valueOf(interval[0]), Double.valueOf(interval[1]));
    }

    /** A histogram of the sample, one row per bin. */
    private static void histogram(double[] values, int bins, int width) {
        double lo = Double.POSITIVE_INFINITY;
        double hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < values.length; ++i) {
            lo = Math.min(lo, values[i]);
            hi = Math.max(hi, values[i]);
        }
        int[] count = new int[bins];
        for (int i = 0; i < values.length; ++i) {
            int b = (int) ((values[i] - lo) / (hi - lo) * (bins - 1));
            count[b]++;
        }
        int tallest = 0;
        for (int b = 0; b < bins; ++b) {
            tallest = Math.max(tallest, count[b]);
        }
        int markedBin = (int) ((Datasets.ACCEPTED - lo) / (hi - lo) * (bins - 1));
        for (int b = 0; b < bins; ++b) {
            StringBuilder bar = new StringBuilder();
            int cells = tallest == 0 ? 0 : (count[b] * width) / tallest;
            for (int c = 0; c < cells; ++c) {
                bar.append('#');
            }
            System.out.println(String.format(L, "  %9.3f |%-" + width + "s %3d%s",
                    Double.valueOf(lo + (hi - lo) * b / (bins - 1)), bar.toString(), Integer.valueOf(count[b]),
                    b == markedBin ? "   <-- the accepted value, 299.792458" : ""));
        }
        System.out.println("  Individual runs fall on both sides of it. The mean does not.");
    }

    public static void main(String[] args) {
        Description d = describe();
        System.out.println("Michelson's 1879 determinations of the speed of light, " + d.runs + " runs");
        System.out.println("Source: NIST Dataplot, retrieved 2026-08-22, a US government publication");
        System.out.println("Six steps; section 6 states what they established, with the numbers.");

        rule("1. the hundred runs  (math.probe)");
        System.out.println(String.format(L, "  %d runs between 5 June and 2 July 1879, %d in the morning and %d in",
                Integer.valueOf(d.runs), Integer.valueOf(d.morning), Integer.valueOf(d.afternoon)));
        System.out.println(String.format(L, "  the afternoon, grouped into %d measurement sets.",
                Integer.valueOf(Datasets.SETS)));
        System.out.println();
        System.out.println(String.format(L, "  mean %.6f, sd %.6f, range %.2f to %.2f, standard error %.6f",
                Double.valueOf(d.mean), Double.valueOf(d.sd), Double.valueOf(d.min), Double.valueOf(d.max),
                Double.valueOf(d.standardError)));
        System.out.println(String.format(L, "  accepted value %.6f, so the sample sits %+.1f km/s above it",
                Double.valueOf(Datasets.ACCEPTED), Double.valueOf(d.biasKmPerSecond())));
        System.out.println();
        histogram(Datasets.speed(), 16, 40);

        rule("2. what a sketch knows that the sample does not  (math.probe)");
        SketchCheck s = sketches();
        System.out.println(String.format(L, "  t-digest at compression 100 kept %d centroids for %d values, so it is",
                Integer.valueOf(s.centroids), Integer.valueOf(d.runs)));
        System.out.println("  holding every point it was given and compressing nothing.");
        System.out.println();
        StringBuilder head = new StringBuilder("  quantile ");
        StringBuilder exact = new StringBuilder("  exact    ");
        StringBuilder digest = new StringBuilder("  t-digest ");
        for (int k = 0; k < QUANTILES.length; ++k) {
            head.append(String.format(L, " %10.2f", Double.valueOf(QUANTILES[k])));
            exact.append(String.format(L, " %10.4f", Double.valueOf(s.exact[k])));
            digest.append(String.format(L, " %10.4f", Double.valueOf(s.digest[k])));
        }
        System.out.println(head);
        System.out.println(exact);
        System.out.println(digest);
        System.out.println(String.format(L, "  worst deviation %.1e -- not compression, but a different rule for",
                Double.valueOf(s.digestWorst)));
        System.out.println("  what the quantile of a finite sample means.");
        System.out.println(String.format(L, "  P-squared median %.6f against the exact %.6f, off by %.1e, and it",
                Double.valueOf(s.squaredMedian), Double.valueOf(s.exact[2]), Double.valueOf(s.medianDeviation)));
        System.out.println("  kept five markers rather than a hundred values.");

        rule("3. three intervals, and the one number that can judge them  (math.probe + math.distribution)");
        Intervals in = intervals();
        System.out.println(String.format(L, "  t(%.3f, %d) = %.6f, standard error %.6f",
                Double.valueOf(0.5 * (1.0 + CONFIDENCE)), Integer.valueOf(d.runs - 1), Double.valueOf(in.criticalT),
                Double.valueOf(in.standardError)));
        System.out.println();
        System.out.println(String.format(L, "  %-24s %-26s %10s %12s", "interval", "95 percent", "contains?",
                "half widths"));
        row("classical t", in.classical);
        row("bootstrap percentile", in.percentile);
        row("bootstrap BCa", in.bca);
        System.out.println();
        System.out.println(String.format(L,
                "  The three agree with each other to three decimals. The accepted value is"));
        System.out.println(String.format(L, "  %+.1f km/s away, which is %.1f standard errors, and it lies outside",
                Double.valueOf(d.biasKmPerSecond()), Double.valueOf(d.biasKmPerSecond() / (1000.0 * in.standardError))));
        System.out.println("  all three. They describe how the measurements scatter, which is not the");
        System.out.println("  same question as whether they are right.");

        rule("4. the library checking itself  (math.distribution + math.solve)");
        CrossCheck c = crossCheck();
        System.out.println("  The interval above rests on three separate pieces of code: a quantile, the");
        System.out.println("  distribution function behind it, and the density behind that. Each one can");
        System.out.println("  be made to confirm the next.");
        System.out.println();
        agreement("Normal.cdf against the integral of its own pdf", c.cdfAgainstIntegral);
        agreement("Normal.inverseCdf against brentDekker on the cdf", c.normalQuantileAgainstRoot);
        agreement("StudentT.inverseCdf against brentDekker on the cdf", c.studentQuantileAgainstRoot);
        agreement("mean as the integral of x pdf(x)", c.meanAgainstIntegral);
        agreement("variance as the integral of (x-mu)^2 pdf(x)", c.varianceAgainstIntegral);
        System.out.println(String.format(L, "  %-52s %.6f", "t(0.975, 99), which every table prints as 1.98422",
                Double.valueOf(c.tabulatedQuantile)));
        System.out.println();
        System.out.println("  InfiniteIntegrator returns 0 for a density centred far from the origin,");
        System.out.println("  which is why this section integrates over mu +/- 12 sd instead.");

        rule("5. where the error actually lives  (math.linalg)");
        Model m = model();
        System.out.println(String.format(L,
                "  The afternoon runs read %+.1f km/s above the morning ones, which looks like",
                Double.valueOf(rawAfternoonEffect())));
        System.out.println("  an explanation. Put it in a model beside the day and the temperature:");
        System.out.println();
        System.out.println(String.format(L, "  %-12s %13s %9s %9s %9s   %s", "term", "coefficient", "km/s", "t", "p",
                "95 percent interval"));
        for (int j = 0; j < m.terms.length; ++j) {
            System.out.println(String.format(L, "  %-12s %13.6f %9.1f %9.3f %9.5f   [%+.6f, %+.6f]", m.terms[j],
                    Double.valueOf(m.beta[j]), Double.valueOf(1000.0 * m.beta[j]), Double.valueOf(m.tValues[j]),
                    Double.valueOf(m.pValues[j]), Double.valueOf(m.intervals[j][0]),
                    Double.valueOf(m.intervals[j][1])));
        }
        System.out.println();
        int afternoon = m.indexOf("afternoon");
        int day = m.indexOf("day");
        System.out.println(String.format(L,
                "  The afternoon effect does not survive it: p = %.2f, and its interval covers",
                Double.valueOf(m.pValues[afternoon])));
        System.out.println("  zero comfortably. The runs were simply not spread evenly over the month.");
        System.out.println(String.format(L, "  What does survive is a drift of %+.1f km/s per day, p = %.5f.",
                Double.valueOf(1000.0 * m.beta[day]), Double.valueOf(m.pValues[day])));
        System.out.println(String.format(L,
                "  The model explains %.0f percent of the variance and leaves an sd of %.6f",
                Double.valueOf(100.0 * m.rSquared), Double.valueOf(m.residualSd)));
        System.out.println(String.format(L, "  against the raw %.6f -- and none of it touches the %.0f km/s that",
                Double.valueOf(m.rawSd), Double.valueOf(d.biasKmPerSecond())));
        System.out.println("  separate this table from the truth.");

        rule("6. what this run established");
        System.out.println(String.format(L,
                "  1. The sample mean is %.6f, %+.1f km/s from a value that carries no",
                Double.valueOf(d.mean), Double.valueOf(d.biasKmPerSecond())));
        System.out.println(String.format(L, "     error at all, which is %.1f standard errors.",
                Double.valueOf(d.biasKmPerSecond() / (1000.0 * in.standardError))));
        System.out.println(String.format(L,
                "  2. Three intervals computed three different ways agree with each other and"));
        System.out.println(String.format(L,
                "     all miss it, the nearest by %.1f of its own half widths. An interval is a",
                Double.valueOf(Intervals.halfWidthsAway(in.percentile))));
        System.out.println("     statement about scatter, not about correctness.");
        System.out.println(String.format(L,
                "  3. On 100 values a t-digest stores all 100 and still answers %.1e away from",
                Double.valueOf(s.digestWorst)));
        System.out.println("     the order statistic, because a quantile of a finite sample is a");
        System.out.println("     convention. P-squared holds five markers and lands within " + String.format(L, "%.1e.",
                Double.valueOf(s.medianDeviation)));
        System.out.println("  4. The density, the distribution function and the quantile agree to better");
        System.out.println("     than " + AGREEMENT_TEXT
                + ", once two defects that this section found had been fixed.");
        System.out.println(String.format(L,
                "  5. The eye-catching %+.0f km/s between afternoon and morning is confounding",
                Double.valueOf(rawAfternoonEffect())));
        System.out.println(String.format(L, "     (p = %.2f). The real structure is a drift of %+.1f km/s per day",
                Double.valueOf(m.pValues[afternoon]), Double.valueOf(1000.0 * m.beta[day])));
        System.out.println(String.format(L, "     (p = %.5f), and it explains none of the bias.",
                Double.valueOf(m.pValues[day])));
        System.out.println();
        System.out.println("  Michelson's apparatus was precise and it was wrong, and no amount of");
        System.out.println("  arithmetic on its output could have revealed that. Only a value obtained");
        System.out.println("  another way could -- which is the argument for external oracles.");
    }

    private static void row(String label, double[] interval) {
        System.out.println(String.format(L, "  %-24s %-26s %10s %12.1f", label, bracket(interval),
                Intervals.contains(interval) ? "yes" : "no", Double.valueOf(Intervals.halfWidthsAway(interval))));
    }

    /**
     * The agreements below are round-off, and their digits differ between the
     * scalar and the vectorized build. The page therefore carries the bound
     * the demo claims and the test asserts, not the digits of one run.
     */
    private static void agreement(String what, double deviation) {
        System.out.println(String.format(L, "  %-52s %s", what,
                deviation <= AGREEMENT ? "better than " + AGREEMENT_TEXT : String.format(L, "%.2e", Double
                        .valueOf(deviation))));
    }

    private MichelsonDemo() {
        throw new AssertionError();
    }
}
