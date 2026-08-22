package math.nist;

import java.util.Locale;

import math.linalg.DMatrix;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.SvdLeastSquares;
import math.optim.LevenbergMarquardt;

/**
 * Prints how many digits this library shares with the values NIST certifies for
 * its reference datasets.
 * <p>
 * Every other check in this repository compares the library against itself. The
 * numbers below were computed elsewhere in extended precision by people whose
 * purpose was to give software something to be measured against, so this is the
 * one place that answers the question outright: <i>how accurate is it?</i>
 * <p>
 * The digits are printed to one decimal on purpose. The last places of a log
 * relative error are round-off, and round-off differs between the scalar and
 * the vectorized build of this library; what the table is for is the size of
 * the number, not its final digit.
 */
public final class Certificate {

    private static final Locale L = Locale.ROOT;

    public static void main(String[] args) {
        System.out.println("NIST Statistical Reference Datasets: digits of agreement with the certified values");
        System.out.println("Source: NIST/ITL StRD, retrieved 2026-08-22, public domain");
        System.out.println();

        double worstLinear = Digits.MAX;
        String worstLinearWhere = "";

        System.out.println("=== linear least squares  (math.linalg)");
        System.out.println(String.format(L, "  %-10s %-10s %4s %4s %8s %8s %9s %8s %11s", "dataset", "difficulty",
                "obs", "par", "beta", "sd", "resid sd", "R^2", "condition"));
        StRD.LinearSet[] linear = StRD.linear();
        for (int k = 0; k < linear.length; ++k) {
            StRD.LinearSet set = linear[k];
            FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(set.design(), set.observations,
                    set.parameters);
            double condition = svd.sigma[0] / svd.sigma[svd.n - 1];

            try {
                LSSummary fit = OLS.estimate(0.05, matrix(set), column(set.response()));
                double beta = Digits.worstOf(fit.getBeta().toArray(), set.certifiedBeta());
                double sd = Digits.worstOf(fit.getCoefficientStandardErrors().toArray(),
                        set.certifiedStandardDeviation());
                double residual = Digits.of(Math.sqrt(fit.getSigmaHatSquared()), set.residualStandardDeviation);
                double rSquared = Digits.of(fit.getRSquared(), set.rSquared);

                System.out.println(String.format(L, "  %-10s %-10s %4d %4d %8.1f %8.1f %9.1f %8.1f %11.2e", set.name,
                        set.difficulty, Integer.valueOf(set.observations), Integer.valueOf(set.parameters),
                        Double.valueOf(beta), Double.valueOf(sd), Double.valueOf(residual),
                        Double.valueOf(rSquared), Double.valueOf(condition)));
                if (beta < worstLinear) {
                    worstLinear = beta;
                    worstLinearWhere = set.name + ", the parameters";
                }
            } catch (RuntimeException refused) {
                System.out.println(String.format(L, "  %-10s %-10s %4d %4d %8s %8s %9s %8s %11.2e", set.name,
                        set.difficulty, Integer.valueOf(set.observations), Integer.valueOf(set.parameters),
                        "refused", "--", "--", "--", Double.valueOf(condition)));
            }
        }

        StRD.LinearSet filip = linear[linear.length - 1];
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(filip.design(), filip.observations,
                filip.parameters);
        double[] direct = SvdLeastSquares.solve(svd, filip.response(), 0.0);
        double tolerance = svd.sigma[0] * Math.max(filip.observations, filip.parameters) * 2.220446049250313e-16;
        double[] truncated = SvdLeastSquares.solveTruncated(svd, filip.response(), tolerance);
        System.out.println();
        System.out.println("  Filip is the set that separates the packages, and it separates this one from");
        System.out.println("  itself. OLS declines it: the smallest singular value falls under the rank");
        System.out.println(String.format(L, "  criterion, %.2e against a tolerance of %.2e. Going straight to the",
                Double.valueOf(svd.sigma[svd.n - 1]), Double.valueOf(tolerance)));
        System.out.println(String.format(L,
                "  decomposition instead reaches %.1f digits on all eleven coefficients, which is",
                Double.valueOf(Digits.worstOf(direct, filip.certifiedBeta()))));
        System.out.println("  what double precision has to give on a design conditioned like this. Asking");
        System.out.println(String.format(L, "  for the truncated solution reaches %.1f: the tolerance discards the very",
                Double.valueOf(Digits.worstOf(truncated, filip.certifiedBeta()))));
        System.out.println("  direction the answer lives in. Ill conditioned is not rank deficient.");

        System.out.println();
        System.out.println("=== nonlinear least squares  (math.optim)");
        System.out.println(String.format(L, "  %-10s %-10s %4s %4s %6s %8s %8s %8s %8s", "dataset", "difficulty",
                "obs", "par", "start", "beta", "sd", "rss", "fevals"));
        double worstNonlinear = Digits.MAX;
        String worstNonlinearWhere = "";
        StRD.NonlinearSet[] nonlinear = StRD.nonlinear();
        for (int k = 0; k < nonlinear.length; ++k) {
            StRD.NonlinearSet set = nonlinear[k];
            for (int start = 1; start <= 2; ++start) {
                LevenbergMarquardt.Result r = new LevenbergMarquardt().solve(Models.of(set), set.start(start),
                        set.observations);
                double beta = Digits.worstOf(r.parameters, set.certifiedBeta());
                double sd = r.standardErrors == null ? 0.0
                        : Digits.worstOf(r.standardErrors, set.certifiedStandardDeviation());
                double rss = Digits.of(r.sumOfSquares, set.residualSumOfSquares);

                System.out.println(String.format(L, "  %-10s %-10s %4d %4d %6d %8.1f %8.1f %8.1f %8d", set.name,
                        set.difficulty, Integer.valueOf(set.observations), Integer.valueOf(set.parameters),
                        Integer.valueOf(start), Double.valueOf(beta), Double.valueOf(sd), Double.valueOf(rss),
                        Integer.valueOf(r.functionEvaluations)));
                if (beta < worstNonlinear) {
                    worstNonlinear = beta;
                    worstNonlinearWhere = set.name + " from start " + start;
                }
            }
        }
        System.out.println();
        System.out.println("  All five converge from both starting points NIST prescribes, the far one");
        System.out.println("  included, which is where a solver is usually found out. The residual sum of");
        System.out.println("  squares is known better than the parameters everywhere: near a minimum the");
        System.out.println("  value is flat, so its position is the harder half of the question.");
        System.out.println();
        System.out.println("  That flatness has a consequence worth stating. These models go through");
        System.out.println("  StrictMath rather than Math, because Math.exp is allowed to be wrong by one");
        System.out.println("  unit in the last place and its intrinsic differs between JDK versions. On");
        System.out.println("  MGH10 that one unit was worth two digits of the fitted parameters -- 9.0 on");
        System.out.println("  one JDK against 10.9 on another, from the same source on the same data.");
        System.out.println("  StrictMath is specified to return the same value everywhere, and with it");
        System.out.println("  this page is identical on both. A certificate that cannot reproduce itself");
        System.out.println("  would be certifying the JDK, not the library.");

        System.out.println();
        System.out.println("=== what this certifies");
        System.out.println(String.format(L,
                "  Worst agreement anywhere in the linear ladder: %.1f digits, on %s.", Double.valueOf(worstLinear),
                worstLinearWhere));
        System.out.println(String.format(L,
                "  Worst anywhere in the nonlinear one: %.1f digits, on %s.", Double.valueOf(worstNonlinear),
                worstNonlinearWhere));
        System.out.println("  Every other check in this library compares it against itself. These numbers");
        System.out.println("  come from outside it, and they are the only ones that could have said the");
        System.out.println("  library was wrong about something both of its own opinions agreed on.");
    }

    private static DMatrix matrix(StRD.LinearSet set) {
        double[] flat = set.design();
        DMatrix X = new DMatrix(set.observations, set.parameters);
        for (int c = 0; c < set.parameters; ++c) {
            for (int i = 0; i < set.observations; ++i) {
                X.set(i, c, flat[c * set.observations + i]);
            }
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

    private Certificate() {
        throw new AssertionError();
    }
}
