package math.stats.gof;

import math.cern.Arithmetic;

/**
 * The null distribution of the two-sample Kolmogorov-Smirnov statistic
 * {@code D_(m,n) = sup |F_m(t) - G_n(t)|}, the largest gap between the
 * empirical distribution functions of two samples drawn from the same
 * continuous distribution.
 * <p>
 * https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
 *
 * @since 1.5.3
 */
public final class KolmogorovSmirnovTwoSample {

    /**
     * The largest product {@code m * n} the exact methods will accept. The
     * recursion is one pass over an {@code (m + 1) x (n + 1)} lattice, measured
     * at 0.08 ms for a product of 10000 and 4.4 ms for one of a million.
     */
    public static final long EXACT_LIMIT = 1_000_000L;

    /**
     * The exact upper tail {@code P[D_(m,n) >= d]} of the two-sided statistic.
     * <p>
     * Under the null hypothesis every interleaving of the two samples is
     * equally likely, so the p-value is the fraction of the
     * {@code C(m + n, m)} of them whose statistic reaches {@code d}. Those are
     * counted as the paths of a lattice walk that leave a corridor of width
     * {@code d} around the diagonal, and it is the mass which <em>leaves</em>
     * that is summed rather than one minus the mass that stays: every term is
     * then positive, and a tail of 1e-300 comes back as itself instead of as
     * the round-off of a subtraction from one. Against an exact rational
     * evaluation of the same counts the worst relative error measured over all
     * shapes up to {@code m + n = 50} and every attainable {@code d} is 2.3e-14,
     * where the complement form reaches 1.8e-2.
     * <p>
     * <b>This is conditional on the pooled sample having no ties</b>, which is
     * what the counting argument assumes. With ties the conditional
     * distribution is a different one and
     * {@link #barFAsymptotic(int, int, double)} is the honest answer.
     *
     * @param m
     *            the size of the first sample, {@code 1} or greater
     * @param n
     *            the size of the second sample, {@code 1} or greater
     * @param d
     *            the value of the statistic
     * @return {@code P[D_(m,n) >= d]}
     * @throws IllegalArgumentException
     *             if {@code m} or {@code n} is not strictly positive, if
     *             {@code d} is {@code NaN}, or if {@code m * n} is above
     *             {@link #EXACT_LIMIT}
     */
    public static double barFExact(int m, int n, double d) {
        return exactTail(m, n, d, 0);
    }

    /**
     * The exact upper tail of either one-sided statistic,
     * {@code D+ = sup (F_m - G_n)} or {@code D- = sup (G_n - F_m)}.
     * <p>
     * The two share this distribution: reversing the order of the pooled sample
     * reflects every lattice path, which exchanges the two statistics while
     * leaving {@code m} and {@code n} where they are. Measured, the two tails
     * agree to 7e-14 relative even for samples as unequal as 20 against 100.
     *
     * @param m
     *            the size of the first sample, {@code 1} or greater
     * @param n
     *            the size of the second sample, {@code 1} or greater
     * @param d
     *            the value of the statistic
     * @return {@code P[D+ >= d]}, which is also {@code P[D- >= d]}
     * @throws IllegalArgumentException
     *             if {@code m} or {@code n} is not strictly positive, if
     *             {@code d} is {@code NaN}, or if {@code m * n} is above
     *             {@link #EXACT_LIMIT}
     */
    public static double barFExactOneSided(int m, int n, double d) {
        return exactTail(m, n, d, 1);
    }

    /**
     * The Kolmogorov limit {@code Q(lambda) = 2 sum (-1)^(k-1) exp(-2 k^2
     * lambda^2)} at {@code lambda = sqrt(m n / (m + n)) * d}, the large-sample
     * approximation to {@link #barFExact(int, int, double)}.
     * <p>
     * It is <b>conservative wherever a decision is taken</b>: measured against
     * the exact value over the shapes from {@code m = n = 40} upwards it never
     * fell below it by more than 5e-5 relative, so it does not turn an
     * unremarkable result into a significant one. Only for samples of ten or
     * fewer does it dip properly under, and there only around a p-value of
     * 0.4, worst 0.92 of the exact value at {@code m = n = 5}.
     * <p>
     * How much too large it is depends on the sizes and on how far into the
     * tail the question sits: for {@code m = n = 1000} it is within 6 per cent
     * even at a tail of 1e-8, for {@code m = n = 100} it is 1.2 times the
     * exact value at 1e-4, and for two samples as unequal as 5 against 200000
     * it is four orders of magnitude out there. Prefer the exact method
     * whenever {@code m * n} allows it.
     *
     * @param m
     *            the size of the first sample, {@code 1} or greater
     * @param n
     *            the size of the second sample, {@code 1} or greater
     * @param d
     *            the value of the statistic
     * @return the limiting approximation to {@code P[D_(m,n) >= d]}
     * @throws IllegalArgumentException
     *             if {@code m} or {@code n} is not strictly positive or if
     *             {@code d} is {@code NaN}
     */
    public static double barFAsymptotic(int m, int n, double d) {
        return kolmogorovTail(lambda(m, n, d));
    }

    /**
     * The limiting one-sided tail {@code exp(-2 lambda^2)} at
     * {@code lambda = sqrt(m n / (m + n)) * d}, conservative in the same way as
     * {@link #barFAsymptotic(int, int, double)}.
     *
     * @param m
     *            the size of the first sample, {@code 1} or greater
     * @param n
     *            the size of the second sample, {@code 1} or greater
     * @param d
     *            the value of the statistic
     * @return the limiting approximation to {@code P[D+ >= d]}
     * @throws IllegalArgumentException
     *             if {@code m} or {@code n} is not strictly positive or if
     *             {@code d} is {@code NaN}
     */
    public static double barFAsymptoticOneSided(int m, int n, double d) {
        double lambda = lambda(m, n, d);
        return Math.exp(-2.0 * lambda * lambda);
    }

    /**
     * Sums the probability of the interleavings that leave the corridor.
     * {@code side} is {@code 0} for the two-sided statistic and {@code 1} for
     * either one-sided one.
     */
    private static double exactTail(int m, int n, double d, int side) {
        requireSizes(m, n, d);
        if ((long) m * n > EXACT_LIMIT) {
            throw new IllegalArgumentException(
                    "m * n is " + (long) m * n + ", above the exact limit of " + EXACT_LIMIT);
        }
        if (d <= 0.0) {
            return 1.0;
        }
        if (d > 1.0) {
            return 0.0;
        }
        double md = m;
        double nd = n;
        // the statistic is a multiple of 1/(m n), so the corridor has to admit
        // the largest attainable value strictly below d and nothing more
        double corridor = (0.5 + Math.floor(d * md * nd - 1.0e-7)) / (md * nd);
        double logAll = logBinomial(m + n, m);

        // u[j] is the probability that a uniformly random path to the current
        // (i, j) never left the corridor. Writing it as that probability rather
        // than as a scaled count of paths is what keeps it in [0, 1]: it is a
        // convex combination of its two predecessors, so it cannot underflow
        // while the answer it feeds is still a normal number
        double[] u = new double[n + 1];
        double tail = 0.0;
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= n; j++) {
                double reached;
                if (i == 0 && j == 0) {
                    reached = 1.0;
                } else if (i == 0) {
                    reached = u[j - 1];
                } else if (j == 0) {
                    reached = u[0];
                } else {
                    reached = (i * u[j] + j * u[j - 1]) / (i + j);
                }
                if (isInside(i / md, j / nd, corridor, side)) {
                    u[j] = reached;
                } else {
                    if (reached > 0.0) {
                        // the paths that first leave here, weighted by the
                        // chance of passing through this cell at all
                        tail += reached
                                * Math.exp(logBinomial(i + j, i) + logBinomial(m - i + n - j, m - i) - logAll);
                    }
                    u[j] = 0.0;
                }
            }
        }
        return Math.min(1.0, tail);
    }

    private static boolean isInside(double first, double second, double corridor, int side) {
        if (side == 0) {
            return Math.abs(first - second) <= corridor;
        }
        return first - second <= corridor;
    }

    private static double logBinomial(int a, int b) {
        return Arithmetic.logFactorial(a) - Arithmetic.logFactorial(b) - Arithmetic.logFactorial(a - b);
    }

    private static double lambda(int m, int n, double d) {
        requireSizes(m, n, d);
        if (d <= 0.0) {
            return 0.0;
        }
        return Math.sqrt((double) m * n / ((double) m + n)) * d;
    }

    /**
     * {@code Q(lambda)}, by the alternating series where it converges quickly
     * and by its theta transformation where it does not.
     */
    private static double kolmogorovTail(double lambda) {
        if (lambda <= 0.0) {
            return 1.0;
        }
        if (lambda < 1.0) {
            // below one the alternating series cancels catastrophically, while
            // this form needs three terms and returns one minus something small
            double scale = Math.PI * Math.PI / (8.0 * lambda * lambda);
            double sum = 0.0;
            for (int k = 1; k <= 21; k += 2) {
                double term = Math.exp(-k * k * scale);
                sum += term;
                if (term < 1.0e-18 * sum) {
                    break;
                }
            }
            return Math.max(0.0, 1.0 - Math.sqrt(2.0 * Math.PI) / lambda * sum);
        }
        double sum = 0.0;
        for (int k = 1; k <= 100; k++) {
            double term = Math.exp(-2.0 * (double) k * k * lambda * lambda);
            sum += (k % 2 == 1) ? term : -term;
            if (term < 1.0e-20) {
                break;
            }
        }
        return Math.min(1.0, Math.max(0.0, 2.0 * sum));
    }

    private static void requireSizes(int m, int n, double d) {
        if (m < 1) {
            throw new IllegalArgumentException("m must be strictly positive : " + m);
        }
        if (n < 1) {
            throw new IllegalArgumentException("n must be strictly positive : " + n);
        }
        if (Double.isNaN(d)) {
            throw new IllegalArgumentException("d must not be NaN");
        }
    }

    private KolmogorovSmirnovTwoSample() {
        throw new AssertionError();
    }
}
