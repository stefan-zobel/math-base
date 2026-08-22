package math.demo.quakes;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import math.distribution.ContinuousDistribution;
import math.distribution.Exponential;
import math.distribution.LogNormal;
import math.distribution.Weibull;
import math.fun.DFunction;
import math.optim.LimitedMemoryBFGS;
import math.optim.Optimizable;
import math.probe.Bootstrap;
import math.probe.CountMinSketch;
import math.probe.PSquaredMedian;
import math.probe.SampleStatistic;
import math.probe.SimpleTDigest;
import math.solve.RootFinder;

/**
 * A worked example on a year of earthquakes: 2115 events of magnitude 5.0 and
 * above, and the question of how often a very large one happens.
 * <p>
 * The Gutenberg-Richter law says the number of events above a magnitude falls
 * off by a constant factor per unit of magnitude, and the exponent {@code b} of
 * that decay is close to 1 nearly everywhere on Earth. Estimating it from a
 * catalog takes one line of arithmetic. Getting it right takes rather more,
 * and this demo is about the difference:
 * <ul>
 * <li><b>The magnitude column is not one measurement.</b> The catalog mixes
 * magnitude scales, and the mixture yields {@code b = 1.20} -- close enough to
 * the textbook value to pass unexamined, and wrong.</li>
 * <li><b>A catalog is not a sample from a continuous distribution.</b>
 * Magnitudes are reported to one decimal and depths are pinned to 10 km
 * whenever they cannot be resolved. Both facts change what a quantile means,
 * and one of the two sketches in {@code math.probe} is visibly hurt by
 * it.</li>
 * <li><b>An extrapolation cannot know where its model stops.</b> The fitted law
 * puts the hundred year event at magnitude 10.2. No fault on Earth is long
 * enough for that, and nothing in the arithmetic says so.</li>
 * </ul>
 * <p>
 * Everything is seeded and two runs of {@link #main(String[])} produce
 * identical output, on either source tree.
 */
public final class QuakeDemo {

    /** Formatting locale, so that the output does not depend on the machine. */
    private static final Locale L = Locale.ROOT;

    /** Natural logarithm of ten, which turns a decay rate into a b value. */
    static final double LN10 = Math.log(10.0);

    /** Half the width of a magnitude bin: magnitudes are reported to one decimal. */
    static final double HALF_BIN = 0.05;

    /** Completeness threshold of the headline fit, chosen in section 3. */
    static final double MC = 5.5;

    /** Thresholds section 3 walks through. */
    static final double[] THRESHOLDS = { 5.0, 5.2, 5.4, 5.5, 5.6, 5.8, 6.0, 6.2 };

    /** Bootstrap replications. */
    static final int BOOT_ITERATIONS = 20000;

    /** Seed of the bootstrap, the year of the catalog. */
    static final long BOOT_SEED = 2025L;

    /**
     * Seed of the frequency sketch. Without one it draws fresh hash functions
     * per instance, and two runs of this demo would not agree.
     */
    static final long SKETCH_SEED = 20260822L;

    /** Confidence level of every interval printed here. */
    static final double CONFIDENCE = 0.95;

    /** Years the return level extrapolates to. */
    static final double RETURN_YEARS = 100.0;

    /** The strongest earthquake ever recorded, Chile, 22 May 1960. */
    static final double STRONGEST_EVER = 9.5;

    /**
     * The bound this demo claims for agreements between two routes to the same
     * number. The measured deviations are orders below it; their digits are
     * round-off and differ between the two source trees, so they do not reach
     * the page.
     */
    static final double AGREEMENT = 1.0e-10;

    // ------------------------------------------------------------ the b value

    /** Magnitudes at or above a threshold, optionally restricted to one scale. */
    static double[] above(double mc, String onlyScale) {
        double[] magnitude = Datasets.magnitude();
        int count = 0;
        for (int i = 0; i < magnitude.length; ++i) {
            if (magnitude[i] >= mc - 1.0e-9 && (onlyScale == null || onlyScale.equals(Datasets.scale(i)))) {
                count++;
            }
        }
        double[] out = new double[count];
        int at = 0;
        for (int i = 0; i < magnitude.length; ++i) {
            if (magnitude[i] >= mc - 1.0e-9 && (onlyScale == null || onlyScale.equals(Datasets.scale(i)))) {
                out[at++] = magnitude[i];
            }
        }
        return out;
    }

    /**
     * The maximum likelihood estimate of {@code b} in closed form. Magnitudes
     * above the threshold are exponential, so the mean excess determines the
     * decay; the half bin corrects for magnitudes being bin centers rather than
     * exact values.
     */
    static double akiB(double[] magnitudes, double mc) {
        double mean = 0.0;
        for (int i = 0; i < magnitudes.length; ++i) {
            mean += magnitudes[i];
        }
        mean /= magnitudes.length;
        return 1.0 / (LN10 * (mean - mc + HALF_BIN));
    }

    /**
     * The same likelihood written for the optimizer, so that the closed form has
     * something independent to be checked against.
     * <p>
     * {@link LimitedMemoryBFGS} maximizes, so {@link #getValue()} returns the
     * log-likelihood of the decay rate {@code beta} and not a loss. The b value
     * is {@code beta / ln 10}.
     */
    static final class TruncatedExponential implements Optimizable.ByGradientValue {

        private final double[] excess;
        private final double total;
        private double beta = 1.0;

        TruncatedExponential(double[] magnitudes, double mc) {
            excess = new double[magnitudes.length];
            double sum = 0.0;
            for (int i = 0; i < magnitudes.length; ++i) {
                excess[i] = magnitudes[i] - mc + HALF_BIN;
                sum += excess[i];
            }
            total = sum;
        }

        @Override
        public double getValue() {
            return excess.length * Math.log(beta) - beta * total;
        }

        @Override
        public void getValueGradient(double[] buffer) {
            buffer[0] = excess.length / beta - total;
        }

        @Override
        public int getNumParameters() {
            return 1;
        }

        @Override
        public void getParameters(double[] buffer) {
            buffer[0] = beta;
        }

        @Override
        public double getParameter(int index) {
            return beta;
        }

        @Override
        public void setParameters(double[] params) {
            beta = params[0];
        }

        @Override
        public void setParameter(int index, double value) {
            beta = value;
        }

        double bValue() {
            return beta / LN10;
        }
    }

    /** The headline fit, by both routes, with what the search reported about itself. */
    static final class Fit {
        final int events;
        final double closedForm;
        final double optimized;
        final double agreement;
        final int iterations;
        final String termination;
        final double gradientNorm;
        final double[] percentile;
        final double[] bca;

        Fit(int events, double closedForm, double optimized, double agreement, int iterations, String termination,
                double gradientNorm, double[] percentile, double[] bca) {
            this.events = events;
            this.closedForm = closedForm;
            this.optimized = optimized;
            this.agreement = agreement;
            this.iterations = iterations;
            this.termination = termination;
            this.gradientNorm = gradientNorm;
            this.percentile = percentile;
            this.bca = bca;
        }
    }

    static Fit fit() {
        final double[] magnitudes = above(MC, Datasets.MOMENT_MAGNITUDE);
        double closed = akiB(magnitudes, MC);

        TruncatedExponential model = new TruncatedExponential(magnitudes, MC);
        LimitedMemoryBFGS optimizer = new LimitedMemoryBFGS(model, 10000, 1.0e-12, 1.0e-12, 5, 1.0e-14, 1.0e-14);
        optimizer.optimize();

        Bootstrap bootstrap = new Bootstrap(magnitudes, new SampleStatistic() {
            @Override
            public double apply(double[] sample) {
                return akiB(sample, MC);
            }
        }, BOOT_ITERATIONS, BOOT_SEED);

        return new Fit(magnitudes.length, closed, model.bValue(), Math.abs(model.bValue() - closed),
                optimizer.getIteration(), optimizer.getTermination().toString(), optimizer.getGradientNorm(),
                bootstrap.getConfidenceInterval(CONFIDENCE), bootstrap.getConfidenceIntervalBCa(CONFIDENCE));
    }

    // ------------------------------------------------- the three tail models

    /** Three candidate laws for the tail, scored the way quantized data require. */
    static final class Models {
        final String[] names;
        final double[] binned;
        final double[] continuous;
        final double[] aic;
        final int best;

        Models(String[] names, double[] binned, double[] continuous, double[] aic, int best) {
            this.names = names;
            this.binned = binned;
            this.continuous = continuous;
            this.aic = aic;
            this.best = best;
        }
    }

    static Models models() {
        double[] magnitudes = above(MC, Datasets.MOMENT_MAGNITUDE);
        double[] excess = new double[magnitudes.length];
        for (int i = 0; i < magnitudes.length; ++i) {
            excess[i] = magnitudes[i] - MC + HALF_BIN;
        }

        double mean = 0.0;
        double logMean = 0.0;
        for (int i = 0; i < excess.length; ++i) {
            mean += excess[i];
            logMean += Math.log(excess[i]);
        }
        mean /= excess.length;
        logMean /= excess.length;
        double logVariance = 0.0;
        for (int i = 0; i < excess.length; ++i) {
            logVariance += (Math.log(excess[i]) - logMean) * (Math.log(excess[i]) - logMean);
        }
        logVariance /= excess.length;

        double shape = Math.PI / Math.sqrt(6.0 * logVariance);
        double scale = Math.exp(logMean + 0.5772156649015329 / shape);

        ContinuousDistribution[] candidates = { new Exponential(1.0 / mean), new Weibull(scale, shape),
                new LogNormal(logMean, Math.sqrt(logVariance)) };
        String[] names = { "exponential", "weibull", "lognormal" };
        int[] parameters = { 1, 2, 2 };

        double[] binned = new double[candidates.length];
        double[] continuous = new double[candidates.length];
        double[] aic = new double[candidates.length];
        int best = 0;
        for (int k = 0; k < candidates.length; ++k) {
            for (int i = 0; i < excess.length; ++i) {
                double probability = candidates[k].cdf(excess[i] + HALF_BIN) - candidates[k].cdf(excess[i] - HALF_BIN);
                binned[k] += Math.log(Math.max(probability, Double.MIN_VALUE));
                continuous[k] += Math.log(Math.max(candidates[k].pdf(excess[i]), Double.MIN_VALUE));
            }
            aic[k] = 2.0 * parameters[k] - 2.0 * binned[k];
            if (aic[k] < aic[best]) {
                best = k;
            }
        }
        return new Models(names, binned, continuous, aic, best);
    }

    // -------------------------------------------------------- the sketches

    /** Quantile by linear interpolation between the order statistics. */
    static double exactQuantile(double[] sorted, double q) {
        double h = (sorted.length - 1) * q;
        int lo = (int) Math.floor(h);
        int hi = Math.min(lo + 1, sorted.length - 1);
        return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
    }

    /** What the two quantile sketches answer on the depths, against the exact values. */
    static final class SketchCheck {
        final double[] quantiles;
        final double[] exact;
        final double[] digest;
        final int centroids;
        final double digestWorst;
        final double squaredMedian;
        final double squaredMiss;
        final int distinctDepths;
        final int unresolved;

        SketchCheck(double[] quantiles, double[] exact, double[] digest, int centroids, double digestWorst,
                double squaredMedian, double squaredMiss, int distinctDepths, int unresolved) {
            this.quantiles = quantiles;
            this.exact = exact;
            this.digest = digest;
            this.centroids = centroids;
            this.digestWorst = digestWorst;
            this.squaredMedian = squaredMedian;
            this.squaredMiss = squaredMiss;
            this.distinctDepths = distinctDepths;
            this.unresolved = unresolved;
        }
    }

    static SketchCheck sketches(double compression) {
        double[] depth = Datasets.depth();
        double[] sorted = depth.clone();
        Arrays.sort(sorted);

        SimpleTDigest digest = new SimpleTDigest(compression);
        PSquaredMedian squared = new PSquaredMedian();
        int unresolved = 0;
        for (int i = 0; i < depth.length; ++i) {
            digest.accept(depth[i]);
            squared.accept(depth[i]);
            if (depth[i] == Datasets.UNRESOLVED_DEPTH) {
                unresolved++;
            }
        }

        double[] quantiles = { 0.05, 0.25, 0.5, 0.75, 0.95, 0.99 };
        double[] exact = new double[quantiles.length];
        double[] got = new double[quantiles.length];
        double worst = 0.0;
        for (int k = 0; k < quantiles.length; ++k) {
            exact[k] = exactQuantile(sorted, quantiles[k]);
            got[k] = digest.getQuantile(quantiles[k]);
            worst = Math.max(worst, Math.abs(got[k] - exact[k]));
        }
        int distinct = 1;
        for (int i = 1; i < sorted.length; ++i) {
            if (sorted[i] != sorted[i - 1]) {
                distinct++;
            }
        }
        double median = exactQuantile(sorted, 0.5);
        return new SketchCheck(quantiles, exact, got, digest.getCentroidCount(), worst, squared.getMedian(),
                Math.abs(squared.getMedian() - median), distinct, unresolved);
    }

    /** What the frequency sketch answers about the regions, against the exact counts. */
    static final class RegionSketch {
        final int depth;
        final int width;
        final int bytes;
        final long worstOverestimate;
        final double meanOverestimate;
        final int undercounts;
        final List<String> topTen;
        final long[] estimated;
        final int[] exact;
        final int keys;

        RegionSketch(int depth, int width, int bytes, long worstOverestimate, double meanOverestimate,
                int undercounts, List<String> topTen, long[] estimated, int[] exact, int keys) {
            this.depth = depth;
            this.width = width;
            this.bytes = bytes;
            this.worstOverestimate = worstOverestimate;
            this.meanOverestimate = meanOverestimate;
            this.undercounts = undercounts;
            this.topTen = topTen;
            this.estimated = estimated;
            this.exact = exact;
            this.keys = keys;
        }
    }

    static Map<String, Integer> exactRegionCounts() {
        Map<String, Integer> counts = new HashMap<String, Integer>();
        for (int i = 0; i < Datasets.size(); ++i) {
            String region = Datasets.region(i);
            Integer had = counts.get(region);
            counts.put(region, Integer.valueOf(had == null ? 1 : had.intValue() + 1));
        }
        return counts;
    }

    static RegionSketch regionSketch(int depth, int width) {
        CountMinSketch<String> sketch = new CountMinSketch<String>(depth, width, 10, SKETCH_SEED);
        for (int i = 0; i < Datasets.size(); ++i) {
            sketch.add(Datasets.region(i));
        }
        Map<String, Integer> exact = exactRegionCounts();

        long worst = 0L;
        long total = 0L;
        int under = 0;
        for (Map.Entry<String, Integer> entry : exact.entrySet()) {
            long over = sketch.estimateCount(entry.getKey()) - entry.getValue().intValue();
            if (over < 0L) {
                under++;
            }
            worst = Math.max(worst, over);
            total += over;
        }

        List<String> top = sketch.getTopK();
        long[] estimated = new long[top.size()];
        int[] trueCounts = new int[top.size()];
        for (int i = 0; i < top.size(); ++i) {
            estimated[i] = sketch.estimateCount(top.get(i));
            trueCounts[i] = exact.get(top.get(i)).intValue();
        }
        return new RegionSketch(depth, width, depth * width * 8, worst, total / (double) exact.size(), under, top,
                estimated, trueCounts, exact.size());
    }

    // ------------------------------------------------------ the return level

    /** The magnitude the fitted law puts once in a hundred years, by two routes. */
    static final class ReturnLevel {
        final double rate;
        final double byClosedForm;
        final double byRoot;
        final double agreement;
        final double lower;
        final double upper;
        final double largestObserved;

        ReturnLevel(double rate, double byClosedForm, double byRoot, double agreement, double lower, double upper,
                double largestObserved) {
            this.rate = rate;
            this.byClosedForm = byClosedForm;
            this.byRoot = byRoot;
            this.agreement = agreement;
            this.lower = lower;
            this.upper = upper;
            this.largestObserved = largestObserved;
        }
    }

    static ReturnLevel returnLevel(Fit fit) {
        double rate = fit.events; // events per year at or above MC
        final double target = 1.0 / (rate * RETURN_YEARS);
        double byClosedForm = magnitudeFor(target, fit.closedForm);

        final Exponential tail = new Exponential(fit.closedForm * LN10);
        double excess = RootFinder.brentDekker(0.0, 20.0, new DFunction() {
            @Override
            public double apply(double x) {
                return (1.0 - tail.cdf(x)) - target;
            }
        }, 1.0e-14);
        double byRoot = MC - HALF_BIN + excess;

        double largest = 0.0;
        double[] magnitude = Datasets.magnitude();
        for (int i = 0; i < magnitude.length; ++i) {
            largest = Math.max(largest, magnitude[i]);
        }
        return new ReturnLevel(rate, byClosedForm, byRoot, Math.abs(byClosedForm - byRoot),
                magnitudeFor(target, fit.percentile[1]), magnitudeFor(target, fit.percentile[0]), largest);
    }

    /** The magnitude whose exceedance probability is {@code target} under a b value. */
    static double magnitudeFor(double target, double b) {
        return MC - HALF_BIN + Math.log(1.0 / target) / (b * LN10);
    }

    // ----------------------------------------------------------------- output

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    /** A deviation as the bound the demo claims, since its digits are round-off. */
    private static String agreement(double deviation) {
        return deviation <= AGREEMENT ? "better than 1e-10" : String.format(L, "%.2e", Double.valueOf(deviation));
    }

    public static void main(String[] args) {
        System.out.println("A year of earthquakes: " + Datasets.size() + " events of magnitude 5.0 and above, 2025");
        System.out.println("Source: the USGS earthquake catalog, retrieved 2026-08-22, public domain");
        System.out.println("Six steps; section 6 states what they established, with the numbers.");

        rule("1. the catalog  (math.probe)");
        double[] magnitude = Datasets.magnitude();
        double[] depth = Datasets.depth();
        SketchCheck s = sketches(100.0);
        System.out.println(String.format(L, "  %d events, magnitude %.1f to %.1f, depth %.1f to %.1f km",
                Integer.valueOf(Datasets.size()), Double.valueOf(min(magnitude)), Double.valueOf(max(magnitude)),
                Double.valueOf(min(depth)), Double.valueOf(max(depth))));
        System.out.println(String.format(L, "  in %d named regions, of which the five largest hold %.0f percent",
                Integer.valueOf(Datasets.regionNames().length), Double.valueOf(100.0 * topFiveShare())));
        System.out.println();
        System.out.println("  Neither column is continuous, and neither is rounded by accident:");
        System.out.println(String.format(L,
                "  magnitudes are reported to one decimal, so %d events take %d distinct values;",
                Integer.valueOf(Datasets.size()), Integer.valueOf(distinct(magnitude))));
        System.out.println(String.format(L,
                "  and %d of them carry a depth of exactly %.0f km, which is what the USGS",
                Integer.valueOf(s.unresolved), Double.valueOf(Datasets.UNRESOLVED_DEPTH)));
        System.out.println("  assigns when the depth cannot be resolved. Section 4 pays for both.");

        rule("2. the column that is not one measurement  (math.probe)");
        System.out.println(String.format(L, "  %-6s %8s %9s %9s %9s", "scale", "events", "mean", "largest", "b"));
        String[] scales = Datasets.scaleNames();
        for (int k = 0; k < scales.length; ++k) {
            double[] m = Datasets.magnitudesOn(scales[k]);
            if (m.length < 20) {
                continue;
            }
            System.out.println(String.format(L, "  %-6s %8d %9.4f %9.1f %9.4f", scales[k], Integer.valueOf(m.length),
                    Double.valueOf(mean(m)), Double.valueOf(max(m)), Double.valueOf(akiB(m, 5.0))));
        }
        double[] mixed = above(5.0, null);
        System.out.println(String.format(L, "  %-6s %8d %9.4f %9.1f %9.4f", "all", Integer.valueOf(mixed.length),
                Double.valueOf(mean(mixed)), Double.valueOf(max(mixed)), Double.valueOf(akiB(mixed, 5.0))));
        System.out.println();
        System.out.println("  The body wave magnitude saturates -- it cannot express a large earthquake");
        System.out.println("  and stops at 6.3, where the moment magnitude reaches 8.8. The USGS prefers");
        System.out.println("  one for small events and the other for large ones, so the column mixes two");
        System.out.println(String.format(L,
                "  scales by construction. Read as one quantity it gives b = %.4f, which is",
                Double.valueOf(akiB(mixed, 5.0))));
        System.out.println("  close enough to the textbook 1.0 that nobody would look twice.");

        rule("3. where the fit becomes honest  (math.optim + math.distribution)");
        System.out.println("  Using the moment magnitude alone costs most of the catalog and raises");
        System.out.println("  the question of where it is complete: below some threshold the small");
        System.out.println("  events were given the other scale and are missing here.");
        System.out.println();
        System.out.println(String.format(L, "  %6s %8s %9s %11s", "Mc", "events", "b", "std error"));
        for (int k = 0; k < THRESHOLDS.length; ++k) {
            double[] m = above(THRESHOLDS[k], Datasets.MOMENT_MAGNITUDE);
            double b = akiB(m, THRESHOLDS[k]);
            System.out.println(String.format(L, "  %6.1f %8d %9.4f %11.4f%s", Double.valueOf(THRESHOLDS[k]),
                    Integer.valueOf(m.length), Double.valueOf(b), Double.valueOf(b / Math.sqrt(m.length)),
                    THRESHOLDS[k] == MC ? "   <-- the fit below" : ""));
        }
        System.out.println();
        System.out.println("  The bias falls away as the threshold rises and the standard error grows to");
        System.out.println("  meet it. There is no reading off the answer from one fit; the shape of the");
        System.out.println(String.format(L, "  column is the argument, and %.1f is where it settles.", Double.valueOf(MC)));
        System.out.println();
        Fit f = fit();
        System.out.println(String.format(L, "  b = %.4f from %d events, bootstrap percentile [%.4f, %.4f],",
                Double.valueOf(f.closedForm), Integer.valueOf(f.events), Double.valueOf(f.percentile[0]),
                Double.valueOf(f.percentile[1])));
        System.out.println(String.format(L, "  BCa [%.4f, %.4f]", Double.valueOf(f.bca[0]), Double.valueOf(f.bca[1])));
        System.out.println(String.format(L,
                "  the same number by LimitedMemoryBFGS on the likelihood: agrees %s,", agreement(f.agreement)));
        System.out.println(String.format(L, "  %d iterations, %s, gradient norm %.0e", Integer.valueOf(f.iterations),
                f.termination, Double.valueOf(f.gradientNorm)));
        System.out.println();
        Models m = models();
        System.out.println("  Is the tail exponential at all? Scored on the bins the data actually have:");
        System.out.println(String.format(L, "  %-14s %14s %14s %10s", "model", "binned", "continuous", "AIC"));
        for (int k = 0; k < m.names.length; ++k) {
            System.out.println(String.format(L, "  %-14s %14.2f %14.2f %10.1f%s", m.names[k],
                    Double.valueOf(m.binned[k]), Double.valueOf(m.continuous[k]), Double.valueOf(m.aic[k]),
                    k == m.best ? "   <-- best" : ""));
        }
        System.out.println("  The two columns disagree about the winner. The binned one is the correct");
        System.out.println("  comparison here, because a magnitude of 5.7 means the interval from 5.65 to");
        System.out.println("  5.75 and not the point; scored as points, a lognormal wins on the pile-up");
        System.out.println("  near the threshold that the rounding itself produced.");

        rule("4. what a sketch costs and what it buys  (math.probe)");
        System.out.println(String.format(L, "  t-digest at compression 100 on the %d depths, %d centroids:",
                Integer.valueOf(Datasets.size()), Integer.valueOf(s.centroids)));
        StringBuilder head = new StringBuilder("  quantile ");
        StringBuilder exact = new StringBuilder("  exact    ");
        StringBuilder digest = new StringBuilder("  t-digest ");
        for (int k = 0; k < s.quantiles.length; ++k) {
            head.append(String.format(L, " %9.2f", Double.valueOf(s.quantiles[k])));
            exact.append(String.format(L, " %9.3f", Double.valueOf(s.exact[k])));
            digest.append(String.format(L, " %9.3f", Double.valueOf(s.digest[k])));
        }
        System.out.println(head);
        System.out.println(exact);
        System.out.println(digest);
        System.out.println(String.format(L, "  worst deviation %.3f km on a range of %.0f km.",
                Double.valueOf(s.digestWorst), Double.valueOf(max(depth) - min(depth))));
        System.out.println();
        System.out.println(String.format(L, "  P-squared holds five markers and answers %.4f for the median,",
                Double.valueOf(s.squaredMedian)));
        System.out.println(String.format(L, "  where the exact value is %.4f -- a miss of %.2f km. It is not a",
                Double.valueOf(s.exact[2]), Double.valueOf(s.squaredMiss)));
        System.out.println(String.format(L,
                "  defect: %d of the %d depths are exactly %.0f, the median falls inside that",
                Integer.valueOf(s.unresolved), Integer.valueOf(Datasets.size()),
                Double.valueOf(Datasets.UNRESOLVED_DEPTH)));
        System.out.println("  block of ties, and a parabola fitted through five markers has nothing to");
        System.out.println("  interpolate there. The t-digest, which stores the values, gets it exactly.");
        System.out.println();
        RegionSketch r = regionSketch(4, 128);
        System.out.println(String.format(L,
                "  Count-Min sketch on the %d regions, depth %d width %d, %d bytes:", Integer.valueOf(r.keys),
                Integer.valueOf(r.depth), Integer.valueOf(r.width), Integer.valueOf(r.bytes)));
        for (int i = 0; i < r.topTen.size() && i < 5; ++i) {
            System.out.println(String.format(L, "    %-32s sketch %5d   exact %5d", r.topTen.get(i),
                    Long.valueOf(r.estimated[i]), Integer.valueOf(r.exact[i])));
        }
        System.out.println(String.format(L,
                "  worst overestimate %d, mean %.2f, undercounts %d -- a Count-Min sketch",
                Long.valueOf(r.worstOverestimate), Double.valueOf(r.meanOverestimate),
                Integer.valueOf(r.undercounts)));
        System.out.println("  can only overcount, and that is the whole of its guarantee.");
        System.out.println(String.format(L,
                "  But %d keys is far too few for it to pay: an exact map of %d entries costs",
                Integer.valueOf(r.keys), Integer.valueOf(r.keys)));
        System.out.println(String.format(L,
                "  about what these %d bytes do and answers exactly. The sketch earns its place",
                Integer.valueOf(r.bytes)));
        System.out.println("  when the key count runs past its width, which is not this catalog.");

        rule("5. the event nobody has seen  (math.solve)");
        ReturnLevel level = returnLevel(f);
        System.out.println(String.format(L, "  %.0f events a year above magnitude %.1f, and a decay of b = %.4f.",
                Double.valueOf(level.rate), Double.valueOf(MC), Double.valueOf(f.closedForm)));
        System.out.println(String.format(L, "  The magnitude exceeded once in %.0f years is then",
                Double.valueOf(RETURN_YEARS)));
        System.out.println(String.format(L, "    %.4f  from the closed form", Double.valueOf(level.byClosedForm)));
        System.out.println(String.format(L, "    %.4f  by inverting the fitted tail with brentDekker",
                Double.valueOf(level.byRoot)));
        System.out.println(String.format(L, "  which agree %s, and the interval on b puts it between %.2f and %.2f.",
                agreement(level.agreement), Double.valueOf(level.lower), Double.valueOf(level.upper)));
        System.out.println();
        System.out.println(String.format(L,
                "  The largest event in this catalog is %.1f. The largest ever recorded is %.1f,",
                Double.valueOf(level.largestObserved), Double.valueOf(STRONGEST_EVER)));
        System.out.println("  in Chile in 1960, and no fault on Earth is long enough to produce the");
        System.out.println("  number above. The arithmetic is sound; the model has simply been asked");
        System.out.println("  about a magnitude where it no longer applies, and it cannot say so.");

        rule("6. what this run established");
        System.out.println(String.format(L,
                "  1. Read as one quantity the magnitude column gives b = %.4f. It mixes two",
                Double.valueOf(akiB(mixed, 5.0))));
        System.out.println(String.format(L,
                "     scales, one of which saturates at %.1f, and the plausible number is the one",
                Double.valueOf(max(Datasets.magnitudesOn("mb")))));
        System.out.println("     nobody should use.");
        System.out.println(String.format(L,
                "  2. On the moment magnitude alone, b settles at %.4f once the threshold clears",
                Double.valueOf(f.closedForm)));
        System.out.println(String.format(L,
                "     incompleteness, with a bootstrap interval of [%.2f, %.2f]. Two independent",
                Double.valueOf(f.percentile[0]), Double.valueOf(f.percentile[1])));
        System.out.println(String.format(L, "     routes to that number agree %s.", agreement(f.agreement)));
        System.out.println("  3. The tail is exponential, as the law says -- but only when the models are");
        System.out.println("     scored on the bins the data have. Scored as if the magnitudes were");
        System.out.println("     exact, a lognormal wins on an artefact of the rounding.");
        System.out.println(String.format(L,
                "  4. A t-digest reproduces the depth quantiles to %.2f km; P-squared misses the",
                Double.valueOf(s.digestWorst)));
        System.out.println(String.format(L,
                "     median by %.2f because %d depths are pinned to one value. A Count-Min sketch",
                Double.valueOf(s.squaredMiss), Integer.valueOf(s.unresolved)));
        System.out.println(String.format(L, "     never undercounts, and on %d keys it saves nothing.",
                Integer.valueOf(r.keys)));
        System.out.println(String.format(L,
                "  5. The fitted law puts the hundred year event at %.2f, which is larger than any",
                Double.valueOf(level.byClosedForm)));
        System.out.println("     earthquake in recorded history and physically impossible.");
        System.out.println();
        System.out.println("  Every number in section 3 is defensible and the one in section 5 is not,");
        System.out.println("  and nothing in the arithmetic distinguishes them. What does is knowing");
        System.out.println("  what the data are: two scales, two roundings, and a law with an upper");
        System.out.println("  limit that the catalog is too short to show.");
    }

    // ------------------------------------------------------------- helpers

    static double topFiveShare() {
        Map<String, Integer> counts = exactRegionCounts();
        int[] all = new int[counts.size()];
        int at = 0;
        for (Map.Entry<String, Integer> entry : counts.entrySet()) {
            all[at++] = entry.getValue().intValue();
        }
        Arrays.sort(all);
        int top = 0;
        for (int i = all.length - 1; i >= 0 && i >= all.length - 5; --i) {
            top += all[i];
        }
        return top / (double) Datasets.size();
    }

    static int distinct(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        int n = 1;
        for (int i = 1; i < sorted.length; ++i) {
            if (sorted[i] != sorted[i - 1]) {
                n++;
            }
        }
        return n;
    }

    static double mean(double[] values) {
        double sum = 0.0;
        for (int i = 0; i < values.length; ++i) {
            sum += values[i];
        }
        return sum / values.length;
    }

    static double min(double[] values) {
        double lo = Double.POSITIVE_INFINITY;
        for (int i = 0; i < values.length; ++i) {
            lo = Math.min(lo, values[i]);
        }
        return lo;
    }

    static double max(double[] values) {
        double hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < values.length; ++i) {
            hi = Math.max(hi, values[i]);
        }
        return hi;
    }

    private QuakeDemo() {
        throw new AssertionError();
    }
}
