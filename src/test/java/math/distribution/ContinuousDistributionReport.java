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
        System.out.println("=== tolerance and depth against cost, on a density with nothing wrong with it");
        System.out.println("    the second moment of a chi square on thirty degrees of freedom, over [0, 500],");
        System.out.println("    whose exact value is 60. Digits first, evaluations under them.");
        final ContinuousDistribution chi = new ChiSquare(30.0);
        final long[] calls = new long[1];
        DFunction integrand = new DFunction() {
            @Override
            public double apply(double x) {
                calls[0]++;
                double c = x - 30.0;
                return c * c * chi.pdf(x);
            }
        };
        double[] tolerances = { 1.0e-10, 1.0e-12, 1.0e-13 };
        int[] depths = { 14, 18, 22, 26 };
        System.out.print(String.format(L, "  %-12s", "tolerance"));
        for (int j = 0; j < depths.length; ++j) {
            System.out.print(String.format(L, " %18s", "depth " + depths[j]));
        }
        System.out.println();
        for (int i = 0; i < tolerances.length; ++i) {
            StringBuilder digits = new StringBuilder(String.format(L, "  %-12.0e", Double.valueOf(tolerances[i])));
            StringBuilder cost = new StringBuilder(String.format(L, "  %-12s", ""));
            for (int j = 0; j < depths.length; ++j) {
                calls[0] = 0;
                double got = AdaptiveGaussKronrod.integrate1DAdaptive(G7_K15.POINTS_15, integrand, 0.0, 500.0,
                        tolerances[i], depths[j]);
                digits.append(String.format(L, " %18.1f",
                        Double.valueOf(ContinuousDistributionTest.relDigits(got, chi.variance()))));
                cost.append(String.format(L, " %18d", Long.valueOf(calls[0])));
            }
            System.out.println(digits);
            System.out.println(cost);
        }
        System.out.println();
        System.out.println("  The accuracy is the same everywhere in that grid. The cost is not, and the");
        System.out.println("  tolerance decides it: between 1e-12 and 1e-13 the work jumps by a factor of");
        System.out.println("  twenty-four, and only past that edge does maxDepth start to matter at all.");
        System.out.println();
        System.out.println("  A panel is accepted when |K - G| falls under its share of the budget, and the");
        System.out.println("  budget is halved at every level. Once a panel is small enough for the integrand");
        System.out.println("  to look like a polynomial on it -- around the seventh level here -- |K - G| is");
        System.out.println("  no longer approximation error but round-off, and round-off is proportional to");
        System.out.println("  what the panel holds. Halving the panel halves what it holds and halves its");
        System.out.println("  share of the budget in the same breath, so the ratio between them is frozen:");
        System.out.println("  measured on the panel over the peak it sits between 2 and 20 from the seventh");
        System.out.println("  level to the twenty-fourth and never trends downwards. If that ratio is above");
        System.out.println("  one, no amount of depth will ever bring it under, and the rule bisects to the");
        System.out.println("  limit; if it is below one, the panel is accepted immediately and depth costs");
        System.out.println("  nothing. Which side it falls on is fixed by the tolerance against the size of");
        System.out.println("  the integral, and not by maxDepth at all.");
        System.out.println();
        System.out.println("  Here 1e-13 on an integral of 60 asks for about seven units in the last place,");
        System.out.println("  and the panels over the peak carry some fifteen times an average panel's share");
        System.out.println("  of the mass, so their round-off lands a few times above what they are allowed.");
        System.out.println("  Ask for 1e-12 and the same 14.7 digits arrive in 750 evaluations. The number to");
        System.out.println("  take away is that a tolerance below roughly 1e-14 times the value of the");
        System.out.println("  integral cannot be met by double arithmetic, and asking for it anyway buys");
        System.out.println("  nothing while costing everything.");

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
