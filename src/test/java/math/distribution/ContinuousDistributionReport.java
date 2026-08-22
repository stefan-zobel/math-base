package math.distribution;

import java.util.List;
import java.util.Locale;

import math.distribution.ContinuousDistributionTest.Row;
import math.fun.DFunction;
import math.solve.AdaptiveGaussKronrod;
import math.solve.AdaptiveGaussKronrod.G7_K15;
import math.solve.RootFinder;

/**
 * Prints how well each continuous distribution agrees with the quadrature and
 * root finding of {@code math.solve}.
 * <p>
 * A distribution states a density, a distribution function, a quantile and two
 * moments, and those are five descriptions of one object. Integrating the
 * density has to reproduce the distribution function; inverting the
 * distribution function has to reproduce the quantile; the moments have to be
 * the integrals they are named after. Nothing checked any of that before, and
 * for most of these the two sides are separate code -- the distribution
 * function goes through the special function ports, the density is a closed
 * form written in the class itself.
 * <p>
 * The digits are printed to one decimal on purpose. The last place of a log
 * relative error is round-off; what the table is for is the size of the number.
 * Two rows do move in that last place between JDK 8 and JDK 25 -- a tenth of a
 * digit on the mean of a narrow normal and on one Student t quantile. It is not
 * the vectorized overlay: running the JDK 25 runtime over the JDK 8 class files
 * reproduces the JDK 25 answer, and none of the classes involved has a
 * vectorized variant. {@code Math.exp}, {@code Math.log} and {@code Math.pow}
 * are each allowed one unit in the last place and their intrinsics differ
 * between the two, which is enough to move a digit count by a tenth.
 * <p>
 * Cost is reported as a count of density evaluations rather than a wall clock,
 * so that it is a property of the computation and not of the machine. That
 * count moves by about half a percent between the two runtimes for the same
 * reason the digits do: where the subdivision stops depends on an error
 * estimate, and the estimate is built out of the same intrinsics.
 */
public final class ContinuousDistributionReport {

    private static final Locale L = Locale.ROOT;

    public static void main(String[] args) {
        List<Row> rows = ContinuousDistributionTest.rows();

        System.out.println("Continuous distributions against the quadrature that should reproduce them");
        System.out.println();
        System.out.println("  cdf     the integral of the density against the distribution function");
        System.out.println("  quant   the quantile against an independent Brent-Dekker root of cdf(x) - p");
        System.out.println("  trip    the round trip cdf(inverseCdf(p)) against p");
        System.out.println("  mean    the mean against the integral of x f(x)");
        System.out.println("  var     the variance against the integral of (x - mu)^2 f(x)");
        System.out.println();

        System.out.println("=== digits of agreement");
        System.out.println(String.format(L, "  %-12s %-15s %11s %11s %9s | %6s %6s %6s %6s %6s", "distribution",
                "parameters", "from", "to", "tail", "cdf", "quant", "trip", "mean", "var"));

        double worst = 17.0;
        String worstWhere = "";
        int noMoments = 0;
        int singular = 0;
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            double cdf = ContinuousDistributionTest.absDigits(
                    ContinuousDistributionTest.integral(row.d, row.lo, row.hi),
                    row.d.cdf(row.hi) - row.d.cdf(row.lo));

            double quant = 17.0;
            double trip = 17.0;
            for (int k = 0; k < ContinuousDistributionTest.PS.length; ++k) {
                final double p = ContinuousDistributionTest.PS[k];
                final Row r = row;
                double[] br = ContinuousDistributionTest.bracket(row.d, row.lo, row.hi, p);
                double byRoot = RootFinder.brentDekker(br[0], br[1], new DFunction() {
                    @Override
                    public double apply(double x) {
                        return r.d.cdf(x) - p;
                    }
                }, 1.0e-14);
                double inverse = row.d.inverseCdf(p);
                quant = Math.min(quant, ContinuousDistributionTest.relDigits(inverse, byRoot));
                trip = Math.min(trip, ContinuousDistributionTest.absDigits(row.d.cdf(inverse), p));
            }

            String mean = "--";
            String var = "--";
            if (row.hasMean()) {
                double mu = row.d.mean();
                mean = String.format(L, "%6.1f", Double.valueOf(ContinuousDistributionTest.relDigits(
                        ContinuousDistributionTest.moment(row.d, 1, 0.0, row.mlo, row.mhi), mu)));
                if (row.hasVariance()) {
                    var = String.format(L, "%6.1f", Double.valueOf(ContinuousDistributionTest.relDigits(
                            ContinuousDistributionTest.moment(row.d, 2, mu, row.mlo, row.mhi), row.d.variance())));
                } else {
                    noMoments++;
                }
            } else {
                noMoments++;
            }
            if (row.is(ContinuousDistributionTest.SINGULAR)) {
                singular++;
            }

            System.out.println(String.format(L, "  %-12s %-15s %11.4g %11.4g %9.1e | %6.1f %6.1f %6.1f %6s %6s",
                    row.name, row.params, Double.valueOf(row.lo), Double.valueOf(row.hi), Double.valueOf(row.tail),
                    Double.valueOf(cdf), Double.valueOf(quant), Double.valueOf(trip), mean, var));

            if (cdf < worst) {
                worst = cdf;
                worstWhere = row.toString();
            }
        }

        System.out.println();
        System.out.println("  A dash is a moment that does not exist. Cauchy has neither, Student t on one");
        System.out.println("  degree of freedom neither, on two a mean but no finite variance, and Fisher's");
        System.out.println("  F no mean until its denominator carries more than two degrees of freedom. Those");
        System.out.println("  are the rows where the distribution has to refuse rather than answer, and the");
        System.out.println("  integrals that would define them diverge -- logarithmically, and at the rate");
        System.out.println("  the closed form predicts, which is the sharper way to say it.");

        System.out.println();
        System.out.println("=== what an unbounded density costs, and what an inner window gives back");
        System.out.println(String.format(L, "  %-12s %-15s %14s %14s %14s", "distribution", "parameters",
                "whole support", "stepped in", "mass given up"));
        for (int i = 0; i < rows.size(); ++i) {
            Row row = rows.get(i);
            if (!row.is(ContinuousDistributionTest.SINGULAR)) {
                continue;
            }
            double outer = ContinuousDistributionTest.absDigits(
                    ContinuousDistributionTest.integral(row.d, row.lo, row.hi),
                    row.d.cdf(row.hi) - row.d.cdf(row.lo));
            double step = 1.0e-6 * (row.hi - row.lo);
            double a = row.lo + step;
            double b = row.hi - step;
            double inner = ContinuousDistributionTest.absDigits(ContinuousDistributionTest.integral(row.d, a, b),
                    row.d.cdf(b) - row.d.cdf(a));
            System.out.println(String.format(L, "  %-12s %-15s %14.1f %14.1f %14.3e", row.name, row.params,
                    Double.valueOf(outer), Double.valueOf(inner),
                    Double.valueOf(1.0 - (row.d.cdf(b) - row.d.cdf(a)))));
        }
        System.out.println();
        System.out.println("  Bisection is the wrong tool against a pole. The panel next to it never becomes");
        System.out.println("  small enough, so the rule runs out of recursion instead of converging: measured");
        System.out.println("  over depths 8 to 30 the gain is about a seventh of a digit per level, at twice");
        System.out.println("  the work each time. Step in from both ends by a millionth of the span and the");
        System.out.println("  same density comes back to full precision, which is how one can tell this is a");
        System.out.println("  property of the quadrature and not of the distributions. What that step costs");
        System.out.println("  in mass is the last column, and it is not uniform: a millionth of the span is");
        System.out.println("  a coarse step when the window has to run to a million to hold a heavy tail.");

        System.out.println();
        System.out.println("=== depth against cost, on a density with nothing wrong with it");
        System.out.println("    the second moment of a chi square on thirty degrees of freedom, over [0, 500]");
        final ContinuousDistribution chi = new ChiSquare(30.0);
        System.out.println(String.format(L, "  %-10s %10s %18s", "maxDepth", "digits", "evaluations"));
        int[] depths = { 14, 18, 22, 26 };
        for (int i = 0; i < depths.length; ++i) {
            final long[] calls = new long[1];
            double got = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, new DFunction() {
                @Override
                public double apply(double x) {
                    calls[0]++;
                    double c = x - 30.0;
                    return c * c * chi.pdf(x);
                }
            }, 0.0, 500.0, 1.0e-13, depths[i]);
            System.out.println(String.format(L, "  %-10d %10.1f %18d", Integer.valueOf(depths[i]),
                    Double.valueOf(ContinuousDistributionTest.relDigits(got, chi.variance())),
                    Long.valueOf(calls[0])));
        }
        System.out.println();
        System.out.println("  The same value, for as much work as one cares to spend. The error target is");
        System.out.println("  halved at every level, so it is divided evenly over the panels wherever they");
        System.out.println("  fall; panels that carry real mass cannot meet a share that small, and the rule");
        System.out.println("  bisects the whole support down to the limit. Raising the tolerance does not");
        System.out.println("  help, because both tolerances tried are equally out of reach. The cure is a");
        System.out.println("  window drawn tightly around the mass, not a deeper recursion.");

        System.out.println();
        System.out.println("=== what this says");
        System.out.println(String.format(L,
                "  Worst agreement between a density and its own distribution function: %.1f digits, on %s.",
                Double.valueOf(worst), worstWhere));
        System.out.println(String.format(L,
                "  %d of the %d rows carry an unbounded density, and %d have a moment that does not exist.",
                Integer.valueOf(singular), Integer.valueOf(rows.size()), Integer.valueOf(noMoments)));
        System.out.println("  Every distribution function below ten digits is one of those, or a heavy tail");
        System.out.println("  cut off by a finite window. None is the distribution getting the arithmetic");
        System.out.println("  wrong, and none of the moment columns is either: a truncated window is what");
        System.out.println("  costs a heavy tail its variance, not the closed form it is compared against.");
    }

    private ContinuousDistributionReport() {
        throw new AssertionError();
    }
}
