package math.demo.hubble;

import java.util.Arrays;
import java.util.Locale;

import math.linalg.DMatrix;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.SvdLeastSquares;
import math.probe.Bootstrap;
import math.probe.SampleStatistic;

/**
 * A worked example on the smallest table in this package: the 24 nebulae with
 * which Hubble established that velocity grows with distance.
 * <p>
 * Fitting a line to 24 points is the first thing anyone learns to do, and this
 * table is where that turns out not to be one question but four. Distance and
 * velocity are both measured with error, and once that is true there is no
 * single least squares problem to solve:
 * <ul>
 * <li>minimizing the vertical residuals gives one slope,</li>
 * <li>minimizing the horizontal ones gives another, larger by a factor of
 * 1.6,</li>
 * <li>minimizing the perpendicular distance -- total least squares -- gives a
 * third, which lands on the second here <b>and on neither if the columns are
 * scaled first</b>, because it is the one estimator of the four that is not
 * scale invariant,</li>
 * <li>and the geometric mean of the first two, which is scale invariant, gives
 * a fourth.</li>
 * </ul>
 * <p>
 * All four are defensible, all four come with confidence intervals, and their
 * disagreement is larger than any of those intervals. That disagreement is a
 * measurement of something the intervals cannot see, and it is the point of the
 * demo. The value known today is about {@link Datasets#MODERN}, which is a
 * sixth of the smallest of them: the error is in the distance scale of 1929 and
 * is invisible to every statistic computed from the table.
 * <p>
 * Everything is seeded and two runs of {@link #main(String[])} produce
 * identical output, on either source tree.
 */
public final class HubbleDemo {

    /** Formatting locale, so that the output does not depend on the machine. */
    private static final Locale L = Locale.ROOT;

    /** Confidence level of every interval printed here. */
    static final double CONFIDENCE = 0.95;

    /** Bootstrap replications. */
    static final int BOOT_ITERATIONS = 20000;

    /** Seed of the bootstrap, the year of the paper. */
    static final long BOOT_SEED = 1929L;

    /** Names of the four estimators, in the order {@link #slopes(int[])} returns them. */
    static final String[] ESTIMATORS = { "ordinary, v on r", "reverse, r on v", "total least squares",
            "geometric mean" };

    // ---------------------------------------------------------------- the fits

    /** Centered sums of squares and products over a selection of rows. */
    static double[] moments(double[] x, double[] y, int[] rows) {
        double mx = 0.0;
        double my = 0.0;
        for (int k = 0; k < rows.length; ++k) {
            mx += x[rows[k]];
            my += y[rows[k]];
        }
        mx /= rows.length;
        my /= rows.length;
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        for (int k = 0; k < rows.length; ++k) {
            double dx = x[rows[k]] - mx;
            double dy = y[rows[k]] - my;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
        return new double[] { mx, my, sxx, syy, sxy };
    }

    /** All 24 rows, in order. */
    static int[] everything() {
        int[] rows = new int[Datasets.size()];
        for (int i = 0; i < rows.length; ++i) {
            rows[i] = i;
        }
        return rows;
    }

    /**
     * The four slopes, in the order of {@link #ESTIMATORS}: ordinary least
     * squares, the reverse regression, total least squares in closed form, and
     * the geometric mean of the first two.
     */
    static double[] slopes(int[] rows) {
        double[] m = moments(Datasets.distance(), Datasets.velocity(), rows);
        double sxx = m[2];
        double syy = m[3];
        double sxy = m[4];
        double ordinary = sxy / sxx;
        double reverse = syy / sxy;
        double b = (syy - sxx) / (2.0 * sxy);
        double total = b + Math.sqrt(b * b + 1.0);
        double geometric = Math.signum(sxy) * Math.sqrt(Math.abs(ordinary * reverse));
        return new double[] { ordinary, reverse, total, geometric };
    }

    /**
     * Total least squares the way a caller would compute it: the right singular
     * vector belonging to the smallest singular value of the centered two-column
     * matrix is normal to the fitted line, so the slope is minus the ratio of
     * its components.
     */
    static double totalLeastSquaresBySvd(double[] x, double[] y) {
        int n = x.length;
        double mx = 0.0;
        double my = 0.0;
        for (int i = 0; i < n; ++i) {
            mx += x[i];
            my += y[i];
        }
        mx /= n;
        my /= n;
        double[] a = new double[2 * n]; // column major, n x 2
        for (int i = 0; i < n; ++i) {
            a[i] = x[i] - mx;
            a[n + i] = y[i] - my;
        }
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(a, n, 2);
        int smallest = svd.sigma[0] <= svd.sigma[1] ? 0 : 1;
        return -svd.V[smallest * 2] / svd.V[smallest * 2 + 1];
    }

    /** Centered and scaled to unit variance. */
    static double[] standardized(double[] v) {
        double mean = 0.0;
        for (int i = 0; i < v.length; ++i) {
            mean += v[i];
        }
        mean /= v.length;
        double squares = 0.0;
        for (int i = 0; i < v.length; ++i) {
            squares += (v[i] - mean) * (v[i] - mean);
        }
        double sd = Math.sqrt(squares / (v.length - 1));
        double[] z = new double[v.length];
        for (int i = 0; i < v.length; ++i) {
            z[i] = (v[i] - mean) / sd;
        }
        return z;
    }

    static double standardDeviation(double[] v) {
        double mean = 0.0;
        for (int i = 0; i < v.length; ++i) {
            mean += v[i];
        }
        mean /= v.length;
        double squares = 0.0;
        for (int i = 0; i < v.length; ++i) {
            squares += (v[i] - mean) * (v[i] - mean);
        }
        return Math.sqrt(squares / (v.length - 1));
    }

    // ------------------------------------------------------- the ordinary fit

    /** The ordinary fit with the inference that comes with it. */
    static final class Ordinary {
        final double intercept;
        final double slope;
        final double slopeStandardError;
        final double tValue;
        final double pValue;
        final double[] slopeInterval;
        final double rSquared;
        final double slopeBySvd;

        Ordinary(double intercept, double slope, double slopeStandardError, double tValue, double pValue,
                double[] slopeInterval, double rSquared, double slopeBySvd) {
            this.intercept = intercept;
            this.slope = slope;
            this.slopeStandardError = slopeStandardError;
            this.tValue = tValue;
            this.pValue = pValue;
            this.slopeInterval = slopeInterval;
            this.rSquared = rSquared;
            this.slopeBySvd = slopeBySvd;
        }
    }

    static Ordinary ordinary() {
        double[] r = Datasets.distance();
        double[] v = Datasets.velocity();
        int n = r.length;

        DMatrix X = new DMatrix(n, 2);
        DMatrix y = new DMatrix(n, 1);
        double[] flat = new double[2 * n];
        for (int i = 0; i < n; ++i) {
            X.set(i, 0, 1.0);
            X.set(i, 1, r[i]);
            y.set(i, 0, v[i]);
            flat[i] = 1.0;
            flat[n + i] = r[i];
        }
        LSSummary fit = OLS.estimate(1.0 - CONFIDENCE, X, y);

        // the same fit from the decomposition, which is the route OLS takes
        // internally and which a caller can now take directly
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(flat, n, 2);
        double[] beta = SvdLeastSquares.solve(svd, y.getArrayUnsafe(), 0.0);

        return new Ordinary(fit.getBeta().toArray()[0], fit.getBeta().toArray()[1],
                fit.getCoefficientStandardErrors().toArray()[1], fit.getTValues().toArray()[1],
                fit.getPValues().toArray()[1], fit.getConfidenceIntervals().get(1).toArray(), fit.getRSquared(),
                beta[1]);
    }

    // -------------------------------------------------------- the intervals

    /** A percentile and a BCa interval for one estimator. */
    static final class Interval {
        final String name;
        final double point;
        final double[] percentile;
        final double[] bca;

        Interval(String name, double point, double[] percentile, double[] bca) {
            this.name = name;
            this.point = point;
            this.percentile = percentile;
            this.bca = bca;
        }

        boolean covers(double value) {
            return percentile[0] <= value && value <= percentile[1];
        }
    }

    /**
     * A case bootstrap: the galaxies are resampled, not the residuals, because
     * both columns carry error and there is no direction in which the residual
     * is the error. {@link Bootstrap} resamples a vector of doubles, so it is
     * handed the row numbers and the statistic reads them back as rows.
     */
    static Interval interval(final int which) {
        double[] indices = new double[Datasets.size()];
        for (int i = 0; i < indices.length; ++i) {
            indices[i] = i;
        }
        Bootstrap bootstrap = new Bootstrap(indices, new SampleStatistic() {
            @Override
            public double apply(double[] sample) {
                int[] rows = new int[sample.length];
                for (int k = 0; k < sample.length; ++k) {
                    rows[k] = (int) sample[k];
                }
                return slopes(rows)[which];
            }
        }, BOOT_ITERATIONS, BOOT_SEED);
        return new Interval(ESTIMATORS[which], slopes(everything())[which],
                bootstrap.getConfidenceInterval(CONFIDENCE), bootstrap.getConfidenceIntervalBCa(CONFIDENCE));
    }

    // ----------------------------------------------------------------- output

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    private static String bracket(double[] interval) {
        return String.format(L, "[%7.1f, %7.1f]", Double.valueOf(interval[0]), Double.valueOf(interval[1]));
    }

    /** The 24 objects, with the ordinary and the reverse line drawn through them. */
    private static void scatter(int width, int height) {
        double[] r = Datasets.distance();
        double[] v = Datasets.velocity();
        double[] s = slopes(everything());
        double[] m = moments(r, v, everything());

        double xLo = 0.0;
        double xHi = 2.1;
        double yLo = -300.0;
        double yHi = 1150.0;

        char[][] grid = new char[height][width];
        for (int row = 0; row < height; ++row) {
            Arrays.fill(grid[row], ' ');
        }
        for (int col = 0; col < width; ++col) {
            double x = xLo + (xHi - xLo) * col / (width - 1);
            mark(grid, col, row(m[1] + s[0] * (x - m[0]), yLo, yHi, height), 'o');
            mark(grid, col, row(m[1] + s[1] * (x - m[0]), yLo, yHi, height), '+');
        }
        for (int i = 0; i < r.length; ++i) {
            int col = (int) ((r[i] - xLo) / (xHi - xLo) * (width - 1));
            mark(grid, col, row(v[i], yLo, yHi, height), '*');
        }
        for (int row = 0; row < height; ++row) {
            System.out.println(String.format(L, "  %7.0f |%s", Double.valueOf(yHi - (yHi - yLo) * row / (height - 1)),
                    new String(grid[row])));
        }
        StringBuilder axis = new StringBuilder("          +");
        for (int col = 0; col < width; ++col) {
            axis.append('-');
        }
        System.out.println(axis);
        System.out.println(String.format(L, "           %-" + (width - 4) + ".1f%4.1f", Double.valueOf(xLo),
                Double.valueOf(xHi)));
        System.out.println("           distance in millions of parsecs; velocity in km/s");
        System.out.println("           * a nebula, o the ordinary fit, + the reverse fit");
    }

    private static int row(double value, double lo, double hi, int height) {
        return (int) Math.round((hi - value) / (hi - lo) * (height - 1));
    }

    private static void mark(char[][] grid, int col, int row, char glyph) {
        if (row < 0 || row >= grid.length || col < 0 || col >= grid[0].length) {
            return;
        }
        grid[row][col] = grid[row][col] == ' ' || grid[row][col] == glyph ? glyph : '*';
    }

    public static void main(String[] args) {
        double[] r = Datasets.distance();
        double[] v = Datasets.velocity();
        double[] m = moments(r, v, everything());
        double[] s = slopes(everything());

        System.out.println("Hubble 1929, table 1: " + Datasets.size() + " nebulae with a distance and a velocity");
        System.out.println("Source: Proceedings of the National Academy of Sciences 15(3):168-173");
        System.out.println("Five steps; section 5 states what they established, with the numbers.");

        rule("1. twenty-four points  (math.probe)");
        System.out.println(String.format(L, "  distances %.3f to %.3f million parsecs, velocities %+.0f to %+.0f km/s",
                Double.valueOf(min(r)), Double.valueOf(max(r)), Double.valueOf(min(v)), Double.valueOf(max(v))));
        System.out.println(String.format(L, "  correlation %.4f, and four of the objects share a distance of 2.0",
                Double.valueOf(m[4] / Math.sqrt(m[2] * m[3]))));
        System.out.println("  because their distance is the mean of the Virgo cluster rather than their own.");
        System.out.println();
        scatter(66, 18);

        rule("2. the ordinary fit, and what it says about itself  (math.linalg)");
        Ordinary o = ordinary();
        System.out.println(String.format(L, "  v = %+.2f %+.3f r", Double.valueOf(o.intercept),
                Double.valueOf(o.slope)));
        System.out.println(String.format(L, "  slope %.3f km/s/Mpc, standard error %.3f, t = %.3f, p = %.6f",
                Double.valueOf(o.slope), Double.valueOf(o.slopeStandardError), Double.valueOf(o.tValue),
                Double.valueOf(o.pValue)));
        System.out.println(String.format(L, "  95 percent interval for the slope %s, R^2 = %.4f",
                bracket(o.slopeInterval), Double.valueOf(o.rSquared)));
        System.out.println(String.format(L,
                "  the same slope through SvdLeastSquares.solve on the decomposition: %.3f",
                Double.valueOf(o.slopeBySvd)));
        System.out.println();
        System.out.println("  A p of four in a million and an interval that excludes zero: by every");
        System.out.println("  ordinary reading this settles the question. Section 3 asks it again.");

        rule("3. four ways to fit the same line  (math.linalg)");
        double sdR = standardDeviation(r);
        double sdV = standardDeviation(v);
        double totalRaw = totalLeastSquaresBySvd(r, v);
        double totalScaled = totalLeastSquaresBySvd(standardized(r), standardized(v)) * sdV / sdR;
        System.out.println(String.format(L, "  %-30s %12s", "estimator", "km/s/Mpc"));
        System.out.println(String.format(L, "  %-30s %12.3f", "ordinary, v on r", Double.valueOf(s[0])));
        System.out.println(String.format(L, "  %-30s %12.3f", "reverse, r on v", Double.valueOf(s[1])));
        System.out.println(String.format(L, "  %-30s %12.3f", "total least squares", Double.valueOf(totalRaw)));
        System.out.println(String.format(L, "  %-30s %12.3f", "the same, columns scaled first",
                Double.valueOf(totalScaled)));
        System.out.println(String.format(L, "  %-30s %12.3f", "geometric mean of the first two",
                Double.valueOf(s[3])));
        System.out.println();
        System.out.println(String.format(L, "  The spread of the velocities is %.1f times the spread of the distances,",
                Double.valueOf(sdV / sdR)));
        System.out.println("  and that number decides what total least squares returns. Perpendicular");
        System.out.println("  distance is not a quantity until the axes have units, so on the raw");
        System.out.println(String.format(L, "  columns it collapses onto the reverse regression -- %.3f against %.3f --",
                Double.valueOf(totalRaw), Double.valueOf(s[1])));
        System.out.println("  and on standardized columns it becomes the geometric mean exactly. One");
        System.out.println("  estimator, two answers, chosen by the units and by nothing else.");

        rule("4. intervals that cannot see what is wrong  (math.probe)");
        Interval[] intervals = new Interval[4];
        for (int k = 0; k < intervals.length; ++k) {
            intervals[k] = interval(k);
        }
        System.out.println(String.format(L, "  %d resamples of the 24 galaxies, seed %d", Integer.valueOf(BOOT_ITERATIONS),
                Long.valueOf(BOOT_SEED)));
        System.out.println();
        System.out.println(String.format(L, "  %-22s %8s %19s %19s %7s %6s", "estimator", "point", "95% percentile",
                "95% BCa", "465?", "70?"));
        int coveringHubble = 0;
        int coveringModern = 0;
        for (int k = 0; k < intervals.length; ++k) {
            boolean hubble = intervals[k].covers(Datasets.HUBBLE_K);
            boolean modern = intervals[k].covers(Datasets.MODERN);
            coveringHubble += hubble ? 1 : 0;
            coveringModern += modern ? 1 : 0;
            System.out.println(String.format(L, "  %-22s %8.1f %19s %19s %7s %6s", intervals[k].name,
                    Double.valueOf(intervals[k].point), bracket(intervals[k].percentile), bracket(intervals[k].bca),
                    hubble ? "yes" : "no", modern ? "yes" : "no"));
        }
        System.out.println();
        System.out.println(String.format(L, "  Hubble's own value is %.0f +/- %.0f, from a fit that removes the motion",
                Double.valueOf(Datasets.HUBBLE_K), Double.valueOf(Datasets.HUBBLE_K_ERROR)));
        System.out.println("  of the sun at the same time; table 1 has no coordinates, so that fit cannot");
        System.out.println(String.format(L, "  be repeated here. %d of the four intervals cover it and %d do not,",
                Integer.valueOf(coveringHubble), Integer.valueOf(intervals.length - coveringHubble)));
        System.out.println("  which is already the whole difficulty: the table does not decide.");
        System.out.println(String.format(L,
                "  The value known today is about %.0f, and %s interval above covers it. The",
                Double.valueOf(Datasets.MODERN), coveringModern == 0 ? "no" : "some"));
        System.out.println(String.format(L, "  nearest bound of any of them is %.0f, a factor of %.1f away.",
                Double.valueOf(nearest(intervals)), Double.valueOf(nearest(intervals) / Datasets.MODERN)));

        rule("5. what this run established");
        System.out.println(String.format(L,
                "  1. The ordinary fit is %.1f km/s/Mpc with a standard error of %.1f, a p of",
                Double.valueOf(o.slope), Double.valueOf(o.slopeStandardError)));
        System.out.println(String.format(L, "     %.6f and an interval of [%.0f, %.0f]. Everything a regression",
                Double.valueOf(o.pValue), Double.valueOf(o.slopeInterval[0]), Double.valueOf(o.slopeInterval[1])));
        System.out.println("     reports about itself here is healthy.");
        System.out.println(String.format(L,
                "  2. Fitting the same 24 points the other three ways gives %.0f, %.0f and %.0f.",
                Double.valueOf(s[1]), Double.valueOf(totalRaw), Double.valueOf(s[3])));
        System.out.println(String.format(L,
                "     The spread between defensible estimators is %.0f km/s/Mpc, which is %.1f",
                Double.valueOf(s[1] - s[0]), Double.valueOf((s[1] - s[0]) / o.slopeStandardError)));
        System.out.println("     times the standard error the ordinary fit quotes.");
        System.out.println("  3. Total least squares is not scale invariant, and on this table that is");
        System.out.println(String.format(L, "     not a subtlety: %.0f on the raw columns, %.0f on scaled ones.",
                Double.valueOf(totalRaw), Double.valueOf(totalScaled)));
        System.out.println(String.format(L,
                "  4. Every interval in section 4 excludes %.0f, the value known today; the nearest",
                Double.valueOf(Datasets.MODERN)));
        System.out.println(String.format(L, "     bound any of them reaches is %.0f, a factor of %.1f away. The error is in",
                Double.valueOf(nearest(intervals)), Double.valueOf(nearest(intervals) / Datasets.MODERN)));
        System.out.println("     the distance scale of 1929, and no arithmetic on this table can reach it.");
        System.out.println();
        System.out.println("  The disagreement between the estimators is the only quantity on this page");
        System.out.println("  that points at the truth, and it is the one no confidence interval");
        System.out.println("  contains. Where a model is uncertain, fitting it several ways measures");
        System.out.println("  more than fitting it once and reporting the standard error.");
    }

    /** The lowest bound any of the intervals reaches, percentile or BCa. */
    static double nearest(Interval[] intervals) {
        double best = Double.MAX_VALUE;
        for (int k = 0; k < intervals.length; ++k) {
            best = Math.min(best, Math.min(intervals[k].percentile[0], intervals[k].bca[0]));
        }
        return best;
    }

    private static double min(double[] v) {
        double lo = Double.POSITIVE_INFINITY;
        for (int i = 0; i < v.length; ++i) {
            lo = Math.min(lo, v[i]);
        }
        return lo;
    }

    private static double max(double[] v) {
        double hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; ++i) {
            hi = Math.max(hi, v[i]);
        }
        return hi;
    }

    private HubbleDemo() {
        throw new AssertionError();
    }
}
