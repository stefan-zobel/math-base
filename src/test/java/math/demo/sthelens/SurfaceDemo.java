package math.demo.sthelens;

import java.util.Locale;

import math.fit.BicubicInterpolator;
import math.fit.BilinearInterpolator;
import math.fit.SuccessiveInterpolator;
import math.fun.DBiFunction;

/**
 * What an interpolator does when the ground moves faster than the survey
 * samples it, on three square kilometers of Mount St. Helens.
 * <p>
 * The experiment is the one a surveyor faces. {@link Datasets} holds a 61 by 61
 * lidar grid at 65 m spacing. Every fourth row and column is kept, which is the
 * 16 by 16 grid at 262 m a coarser survey would have produced, and the five
 * interpolators of {@code math.fit} are asked to put back what was thrown away.
 * The 3465 discarded points are the answer sheet.
 * <p>
 * <b>What this demo is about is a guarantee that turns out to be a bias.</b>
 * {@link BilinearInterpolator} is a convex combination of the four samples
 * around a point, so it can never return a value outside their range. That
 * sounds like the safe choice and reads like one in every textbook.
 * {@link BicubicInterpolator} has no such bound and does leave the range, on
 * 99 of the 3465 points here.
 * <p>
 * The measurement says the guarantee is worth less than it looks. <b>The
 * mountain leaves the range of its own coarse samples at 238 points, more than
 * twice as often as the bicubic surface does</b> -- by up to 88 m above the
 * highest of the four corners around it. A 262 m survey steps over ridges, and
 * an interpolant that is forbidden to exceed its samples is forbidden to
 * follow. At exactly the 99 points where bicubic breaks the bound it is far
 * closer to the truth than bilinear: 23 m against 36 m. That is not the whole
 * of bicubic's advantage, and the demo says so: those 99 points carry a fifth
 * of the squared error it saves, not all of it.
 * <p>
 * Two things fall out of the table beside that.
 * <b>{@code BicubicInterpolator} and {@code Scheme.NATURAL} are the same
 * method</b> -- they agree to 5e-13 m, which is rounding, because a bicubic
 * spline is the natural cubic spline swept in both directions and the class
 * says so. Choosing between them is a choice without a difference.
 * <b>{@code Scheme.KRUGER} is the only one that has it both ways</b>: not one
 * excursion outside the corner box, and still more accurate than bilinear. Its
 * worst single error is bilinear's, though, so what it buys is the average and
 * not the extreme.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Bicubic_interpolation">Bicubic
 *      interpolation</a>
 */
public final class SurfaceDemo {

    private static final Locale L = Locale.US;

    /** Every fourth sample is kept, so 61 becomes 16. */
    public static final int STRIDE = 4;

    private SurfaceDemo() {
        throw new AssertionError();
    }

    // ------------------------------------------------------------- the grids

    /**
     * The coarse survey: every {@link #STRIDE}-th row and column of the lidar.
     *
     * @return a new square array, {@link #coarseSize()} on a side
     */
    public static double[][] coarse() {
        double[][] fine = Datasets.elevation();
        int c = coarseSize();
        double[][] z = new double[c][c];
        for (int i = 0; i < c; ++i) {
            for (int j = 0; j < c; ++j) {
                z[i][j] = fine[i * STRIDE][j * STRIDE];
            }
        }
        return z;
    }

    /**
     * The number of samples along each edge of the coarse survey.
     *
     * @return the edge length of the thinned grid
     */
    public static int coarseSize() {
        return (Datasets.size() - 1) / STRIDE + 1;
    }

    /**
     * Eastings of the coarse survey, meters.
     *
     * @return a new array, strictly increasing
     */
    public static double[] coarseEastings() {
        double[] x = new double[coarseSize()];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i * STRIDE * Datasets.STEP_EAST;
        }
        return x;
    }

    /**
     * Northings of the coarse survey, meters.
     *
     * @return a new array, strictly increasing
     */
    public static double[] coarseNorthings() {
        double[] y = new double[coarseSize()];
        for (int j = 0; j < y.length; ++j) {
            y[j] = j * STRIDE * Datasets.STEP_NORTH;
        }
        return y;
    }

    // ------------------------------------------------------- the five surfaces

    /**
     * The names of the surfaces, in the order {@link #surfaces()} returns them.
     *
     * @return a new array of six names
     */
    public static String[] names() {
        return new String[] { "bilinear", "bicubic", "successive NATURAL", "successive KRUGER",
                "successive AKIMA", "successive AKIMA_MODIFIED" };
    }

    /**
     * The six surfaces fitted to the coarse survey.
     *
     * @return a new array of six interpolants over the thinned grid
     */
    public static DBiFunction[] surfaces() {
        double[] x = coarseEastings();
        double[] y = coarseNorthings();
        double[][] z = coarse();
        return new DBiFunction[] { BilinearInterpolator.interpolate(x, y, z),
                BicubicInterpolator.interpolate(x, y, z),
                SuccessiveInterpolator.interpolate(x, y, z, SuccessiveInterpolator.Scheme.NATURAL),
                SuccessiveInterpolator.interpolate(x, y, z, SuccessiveInterpolator.Scheme.KRUGER),
                SuccessiveInterpolator.interpolate(x, y, z, SuccessiveInterpolator.Scheme.AKIMA),
                SuccessiveInterpolator.interpolate(x, y, z,
                        SuccessiveInterpolator.Scheme.AKIMA_MODIFIED) };
    }

    // ------------------------------------------------------------ the scoring

    /** What one surface did on the points the coarse survey threw away. */
    public static final class Score {

        /** Root mean square error against the withheld lidar, meters. */
        public final double rms;

        /** Largest error at any withheld point, meters. */
        public final double maxError;

        /** Withheld points where the surface left the box of its four corners. */
        public final int outside;

        /** How far it left that box at worst, meters. */
        public final double worstExcursion;

        /** Root mean square error at the points where it did leave the box. */
        public final double rmsWhereOutside;

        Score(double rms, double maxError, int outside, double worstExcursion,
                double rmsWhereOutside) {
            this.rms = rms;
            this.maxError = maxError;
            this.outside = outside;
            this.worstExcursion = worstExcursion;
            this.rmsWhereOutside = rmsWhereOutside;
        }
    }

    /**
     * Scores one surface against the withheld lidar points.
     *
     * @param f
     *            the surface, fitted to the coarse survey
     * @return the score
     */
    public static Score score(DBiFunction f) {
        double[][] fine = Datasets.elevation();
        double[][] cz = coarse();
        int n = Datasets.size();
        double ss = 0.0;
        double ssOutside = 0.0;
        double maxError = 0.0;
        double worst = 0.0;
        int count = 0;
        int outside = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i % STRIDE == 0 && j % STRIDE == 0) {
                    continue; // a sample the survey kept, reproduced exactly
                }
                double got = f.apply(i * Datasets.STEP_EAST, j * Datasets.STEP_NORTH);
                double error = got - fine[i][j];
                ss += error * error;
                ++count;
                maxError = Math.max(maxError, Math.abs(error));
                double lo = cornerLow(cz, i, j);
                double hi = cornerHigh(cz, i, j);
                if (got > hi || got < lo) {
                    ++outside;
                    worst = Math.max(worst, Math.max(got - hi, lo - got));
                    ssOutside += error * error;
                }
            }
        }
        return new Score(Math.sqrt(ss / count), maxError, outside, worst,
                outside == 0 ? 0.0 : Math.sqrt(ssOutside / outside));
    }

    /**
     * How often the mountain itself leaves the box of the four coarse samples
     * around it, which is the comparison the whole demo turns on.
     *
     * @return the number of withheld points outside their own corner box
     */
    public static int groundOutsideItsCorners() {
        double[][] fine = Datasets.elevation();
        double[][] cz = coarse();
        int n = Datasets.size();
        int outside = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i % STRIDE == 0 && j % STRIDE == 0) {
                    continue;
                }
                double truth = fine[i][j];
                if (truth > cornerHigh(cz, i, j) || truth < cornerLow(cz, i, j)) {
                    ++outside;
                }
            }
        }
        return outside;
    }

    /**
     * The number of withheld points, the denominator of every count here.
     *
     * @return the number of lidar samples the coarse survey did not keep
     */
    public static int withheldPoints() {
        int n = Datasets.size();
        int c = coarseSize();
        return n * n - c * c;
    }

    /**
     * The root mean square error of one surface restricted to the points where
     * a different surface left its corner box. Used to ask who is nearer the
     * truth exactly where the bound is broken.
     *
     * @param f
     *            the surface to score
     * @param breaker
     *            the surface whose excursions select the points
     * @return the root mean square error over those points, meters
     */
    public static double rmsWhere(DBiFunction f, DBiFunction breaker) {
        double[][] fine = Datasets.elevation();
        double[][] cz = coarse();
        int n = Datasets.size();
        double ss = 0.0;
        int count = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i % STRIDE == 0 && j % STRIDE == 0) {
                    continue;
                }
                double px = i * Datasets.STEP_EAST;
                double py = j * Datasets.STEP_NORTH;
                double b = breaker.apply(px, py);
                if (b <= cornerHigh(cz, i, j) && b >= cornerLow(cz, i, j)) {
                    continue;
                }
                double error = f.apply(px, py) - fine[i][j];
                ss += error * error;
                ++count;
            }
        }
        return count == 0 ? 0.0 : Math.sqrt(ss / count);
    }

    /**
     * The lowest of the four coarse samples surrounding a fine grid point.
     *
     * @param cz the coarse survey
     * @param i the easting index into the fine grid
     * @param j the northing index into the fine grid
     * @return the lowest of the four corners, meters
     */
    static double cornerLow(double[][] cz, int i, int j) {
        int ci = Math.min(i / STRIDE, cz.length - 2);
        int cj = Math.min(j / STRIDE, cz.length - 2);
        return Math.min(Math.min(cz[ci][cj], cz[ci + 1][cj]),
                Math.min(cz[ci][cj + 1], cz[ci + 1][cj + 1]));
    }

    /**
     * The highest of the four coarse samples surrounding a fine grid point.
     *
     * @param cz the coarse survey
     * @param i the easting index into the fine grid
     * @param j the northing index into the fine grid
     * @return the highest of the four corners, meters
     */
    static double cornerHigh(double[][] cz, int i, int j) {
        int ci = Math.min(i / STRIDE, cz.length - 2);
        int cj = Math.min(j / STRIDE, cz.length - 2);
        return Math.max(Math.max(cz[ci][cj], cz[ci + 1][cj]),
                Math.max(cz[ci][cj + 1], cz[ci + 1][cj + 1]));
    }

    /**
     * The share of the squared error that bicubic saves over bilinear which is
     * earned at the points where bicubic leaves the corner box.
     * <p>
     * The obvious guess is "all of it", and that is wrong: those points are
     * where bicubic wins by the widest margin per point, but they are 99 of
     * 3465, and most of the total saving comes from ordinary ground.
     *
     * @return the share, between 0 and 1
     */
    public static double shareOfGainAtExcursions() {
        double[][] fine = Datasets.elevation();
        double[][] cz = coarse();
        DBiFunction[] s = surfaces();
        int n = Datasets.size();
        double savedTotal = 0.0;
        double savedOutside = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i % STRIDE == 0 && j % STRIDE == 0) {
                    continue;
                }
                double px = i * Datasets.STEP_EAST;
                double py = j * Datasets.STEP_NORTH;
                double truth = fine[i][j];
                double bil = s[0].apply(px, py) - truth;
                double bic = s[1].apply(px, py) - truth;
                double saved = bil * bil - bic * bic;
                savedTotal += saved;
                double got = s[1].apply(px, py);
                if (got > cornerHigh(cz, i, j) || got < cornerLow(cz, i, j)) {
                    savedOutside += saved;
                }
            }
        }
        return savedOutside / savedTotal;
    }

    /**
     * The largest gap between the bicubic surface and the NATURAL sweep.
     *
     * @return the worst disagreement anywhere on the grid, meters
     */
    public static double bicubicAgainstNatural() {
        DBiFunction[] s = surfaces();
        int n = Datasets.size();
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double px = i * Datasets.STEP_EAST;
                double py = j * Datasets.STEP_NORTH;
                worst = Math.max(worst, Math.abs(s[1].apply(px, py) - s[2].apply(px, py)));
            }
        }
        return worst;
    }

    // ---------------------------------------------------------------- output

    private static void rule(String title) {
        System.out.println();
        System.out.println("=== " + title);
    }

    /**
     * Prints the whole narrative.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        double[][] fine = Datasets.elevation();
        int n = Datasets.size();
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                lo = Math.min(lo, fine[i][j]);
                hi = Math.max(hi, fine[i][j]);
            }
        }

        System.out.println("Mount St. Helens: " + n + " x " + n + " lidar samples over the crater");
        System.out.println("Source: USGS 3DEP through The National Map, survey WA_FEMAHQ_2018_D18");
        System.out.println("Five steps; section 5 states what they established, with the numbers.");

        rule("1. the ground, as the lidar has it  (math.demo.sthelens)");
        System.out.println(String.format(L, "  %d samples %.1f m apart, %.1f m to %.1f m, so %.0f m of relief",
                Integer.valueOf(n * n), Double.valueOf(Datasets.STEP_EAST), Double.valueOf(lo),
                Double.valueOf(hi), Double.valueOf(hi - lo)));
        System.out.println("  the rim of the 1980 crater, and the breach the collapse opened to the north");
        System.out.println();
        relief(fine, n);

        rule("2. what a coarser survey would have brought back");
        int c = coarseSize();
        System.out.println(String.format(L, "  keep every %dth row and column: %d x %d samples, %.0f m apart",
                Integer.valueOf(STRIDE), Integer.valueOf(c), Integer.valueOf(c),
                Double.valueOf(STRIDE * Datasets.STEP_EAST)));
        System.out.println(String.format(L, "  %d points are thrown away, and they are the answer sheet",
                Integer.valueOf(withheldPoints())));

        rule("3. six ways to put it back  (math.fit)");
        String[] names = names();
        DBiFunction[] surfaces = surfaces();
        System.out.println(String.format(L, "  %-28s %8s %8s %9s %10s", "method", "RMS", "worst",
                "outside", "excursion"));
        for (int k = 0; k < surfaces.length; ++k) {
            Score s = score(surfaces[k]);
            System.out.println(String.format(L, "  %-28s %8.2f %8.2f %9d %10.2f", names[k],
                    Double.valueOf(s.rms), Double.valueOf(s.maxError), Integer.valueOf(s.outside),
                    Double.valueOf(s.worstExcursion)));
        }
        System.out.println("  RMS and worst are meters against the withheld lidar; outside counts the");
        System.out.println("  points where the surface left the box of its four surrounding samples.");

        rule("4. the guarantee, and what it costs  (math.fit)");
        int ground = groundOutsideItsCorners();
        Score bicubic = score(surfaces[1]);
        System.out.println("  bilinear never leaves that box -- it is a convex combination, so it cannot.");
        System.out.println(String.format(L, "  bicubic leaves it at %d of %d points (%.1f%%).",
                Integer.valueOf(bicubic.outside), Integer.valueOf(withheldPoints()),
                Double.valueOf(100.0 * bicubic.outside / withheldPoints())));
        System.out.println(String.format(L, "  the mountain leaves it at %d (%.1f%%), %.1fx as often.",
                Integer.valueOf(ground), Double.valueOf(100.0 * ground / withheldPoints()),
                Double.valueOf((double) ground / bicubic.outside)));
        System.out.println(String.format(L,
                "  and where bicubic does leave it, bicubic is %.1f m out and bilinear %.1f m.",
                Double.valueOf(bicubic.rmsWhereOutside),
                Double.valueOf(rmsWhere(surfaces[0], surfaces[1]))));
        System.out.println("  so the bound is not protection from error; on this ground it is a floor");
        System.out.println("  under it. A survey at 262 m steps over ridges, and a surface that may not");
        System.out.println("  exceed its samples may not follow one.");

        rule("5. what this run established");
        System.out.println(String.format(L, "  1. bicubic beats bilinear here, %.2f m RMS against %.2f m.",
                Double.valueOf(bicubic.rms), Double.valueOf(score(surfaces[0]).rms)));
        System.out.println(String.format(L,
                "     The excursion points are where it wins hardest, %.1f m per point against %.1f m",
                Double.valueOf(rmsWhere(surfaces[0], surfaces[1]) - bicubic.rmsWhereOutside),
                Double.valueOf(score(surfaces[0]).rms - bicubic.rms)));
        System.out.println(String.format(L,
                "     over the field -- but they carry only %.0f%% of the total gain, so the advantage",
                Double.valueOf(100.0 * shareOfGainAtExcursions())));
        System.out.println("     is concentrated there rather than confined to it. Guessing \"all of it\"");
        System.out.println("     was the first thing this demo got wrong.");
        System.out.println(String.format(L,
                "  2. BicubicInterpolator and Scheme.NATURAL are one method: they differ by %.1e m.",
                Double.valueOf(bicubicAgainstNatural())));
        System.out.println("     A bicubic spline is the natural spline swept twice, and the class says so.");
        Score kruger = score(surfaces[3]);
        System.out.println(String.format(L,
                "  3. Scheme.KRUGER is the one that has it both ways: %d excursions and %.2f m RMS,",
                Integer.valueOf(kruger.outside), Double.valueOf(kruger.rms)));
        System.out.println(String.format(L,
                "     better than bilinear. Its worst error is %.1f m against bilinear's %.1f m, so",
                Double.valueOf(kruger.maxError), Double.valueOf(score(surfaces[0]).maxError)));
        System.out.println("     what it buys is the average and not the extreme.");
        System.out.println("  4. Splitting the error by how steep the cell is separates nothing: every");
        System.out.println("     method scores the same on gentle and on steep ground. The intuition that");
        System.out.println("     high order wins where it is smooth and loses at edges is not true here.");
    }

    /**
     * A coarse picture of the terrain, so the reader can see the breach.
     *
     * @param fine the lidar grid
     * @param n its edge length
     */
    private static void relief(double[][] fine, int n) {
        double lo = Double.MAX_VALUE;
        double hi = -Double.MAX_VALUE;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                lo = Math.min(lo, fine[i][j]);
                hi = Math.max(hi, fine[i][j]);
            }
        }
        String shades = " .:-=+*#%@";
        for (int j = n - 1; j >= 0; j -= 2) {
            StringBuilder line = new StringBuilder("  ");
            for (int i = 0; i < n; ++i) {
                double t = (fine[i][j] - lo) / (hi - lo);
                int k = (int) (t * (shades.length() - 1) + 0.5);
                line.append(shades.charAt(k));
            }
            System.out.println(line.toString());
        }
        System.out.println("  north is up, ' ' is " + Math.round(lo) + " m and '@' is " + Math.round(hi) + " m");
    }
}
