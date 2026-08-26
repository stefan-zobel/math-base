package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

import math.cern.FastGamma;

final class BinomialSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    /**
     * The value of {@code n min(p, 1-p)} above which the transformed rejection
     * is used instead of inversion.
     * <p>
     * Both are exact, so this decides speed and nothing else: inversion costs
     * about {@code n min(p, 1-p)} steps and rejection a constant few.
     * <p>
     * The literature puts this at ten. Measured here it belongs near fifty -- at
     * {@code n = 100}, {@code p = 0.4} inversion takes 66 nanoseconds a draw
     * against rejection's 77, and at {@code n = 150} it is 89 against 68 -- for
     * the reason the Poisson sampler gives: the rare path here costs two
     * log-gammas rather than a table lookup.
     */
    static final double REJECTION_THRESHOLD = 50.0;

    final int n;
    final double p;
    /** Whether the count is taken as failures and turned round at the end. */
    final boolean reflect;
    /** The derived constants, built once and shared with every split, or {@code null} below the threshold. */
    final Rejection rejection;
    final PseudoRandom prng;

    BinomialSpliterator(PseudoRandom prng, long index, long fence, int n, double p) {
        super(index, fence);
        checkParameters(n, p);
        this.n = n;
        this.p = p;
        this.reflect = p > 0.5;
        double q = reflect ? 1.0 - p : p;
        this.rejection = (n > 0 && q > 0.0 && n * q >= REJECTION_THRESHOLD) ? new Rejection(n, q) : null;
        this.prng = prng;
    }

    private BinomialSpliterator(PseudoRandom prng, long index, long fence, int n, double p, boolean reflect,
            Rejection rejection) {
        super(index, fence);
        this.n = n;
        this.p = p;
        this.reflect = reflect;
        this.rejection = rejection;
        this.prng = prng;
    }

    static void checkParameters(int n, double p) {
        if (n < 0) {
            throw new IllegalArgumentException("n < 0 : " + n);
        }
        if (Double.isNaN(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in [0, 1] : " + p);
        }
    }

    @Override
    public Spliterator.OfInt trySplit() {
        long s = splitPoint();
        if (s < 0L) {
            return null;
        }
        PseudoRandom half = detach(prng);
        if (half == null) {
            // the source cannot hand out an independent generator
            return null;
        }
        long idx = index;
        index = s;
        // the constants are immutable and depend only on n and p, so the two
        // halves share them rather than deriving them again
        return new BinomialSpliterator(half, idx, s, n, p, reflect, rejection);
    }

    private int next(PseudoRandom pr) {
        if (n == 0 || p <= 0.0) {
            return 0;
        }
        if (p >= 1.0) {
            return n;
        }
        double q = reflect ? 1.0 - p : p;
        int k = rejection == null ? byInversion(pr, n, q) : rejection.sample(pr);
        return reflect ? n - k : k;
    }

    @Override
    public boolean tryAdvance(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(next(prng));
            index = idx + 1;
            return true;
        } else {
            return false;
        }
    }

    @Override
    public void forEachRemaining(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            index = fence_;
            PseudoRandom pr = prng;
            do {
                consumer.accept(next(pr));
            } while (++idx < fence_);
        }
    }

    /**
     * One draw, deriving whatever constants it needs.
     * <p>
     * Only {@code p <= 0.5} is implemented; the other half of the range is
     * reached by reflection, counting the failures instead of the successes.
     * That is exact -- the two are the same experiment read the other way round
     * -- and it halves both the code and the range any approximation must hold
     * over. The streams above do not use this, since they keep the constants.
     */
    static int sample(PseudoRandom prng, int n, double p) {
        if (n == 0 || p <= 0.0) {
            return 0;
        }
        if (p >= 1.0) {
            return n;
        }
        boolean reflect = p > 0.5;
        double q = reflect ? 1.0 - p : p;
        int k = n * q < REJECTION_THRESHOLD ? byInversion(prng, n, q) : new Rejection(n, q).sample(prng);
        return reflect ? n - k : k;
    }

    /**
     * Inversion by the mass function's own recurrence,
     * {@code f(k+1) = f(k) (n-k) p / ((k+1) (1-p))}.
     * <p>
     * Exact, and it costs about {@code n p} steps, which is why it is only used
     * below the crossover. It starts from {@code (1-p)^n}, which stays a normal
     * number as long as {@code n p} does: below the threshold that quantity is
     * at least {@code exp(-10)}.
     */
    static int byInversion(PseudoRandom prng, int n, double p) {
        double q = 1.0 - p;
        double ratio = p / q;
        double numerator = (n + 1) * ratio;
        double term = Math.pow(q, n);
        double u = prng.nextDouble();
        int k = 0;
        while (u > term) {
            u -= term;
            k++;
            if (k >= n) {
                return n;
            }
            term *= (numerator / k - ratio);
        }
        return k;
    }

    /** One draw by transformed rejection, deriving the constants first. */
    static int byTransformedRejection(PseudoRandom prng, int n, double p) {
        return new Rejection(n, p).sample(prng);
    }

    /**
     * Hoermann's BTRS: transformed rejection with the same hyperbolic hat the
     * Poisson sampler uses, accepting in constant expected time.
     * <p>
     * The constants depend only on {@code n} and {@code p}, so they are derived
     * once and kept. That matters more here than for the Poisson: the
     * acceptance needs the mass function at the mode, and computing its two
     * log-gammas per draw made the whole thing four times slower than the
     * inversion it was supposed to beat.
     * <p>
     * As there, the two cheap tests come first and the slow path compares
     * against the logarithm of the mass function itself. The normalization
     * {@code log(n!)} is common to the candidate and the mode and cancels, so
     * two log-gammas per candidate suffice and no Stirling correction table is
     * needed.
     * <p>
     * <b>Expects {@code p <= 0.5} and {@code n p} of about ten or more</b>,
     * which is what the constants were derived for.
     */
    static final class Rejection {

        private final int n;
        private final double b;
        private final double a;
        private final double c;
        private final double vr;
        private final double alpha;
        private final double logRatio;
        private final double logPmfAtMode;

        Rejection(int n, double p) {
            double q = 1.0 - p;
            double spq = Math.sqrt(n * p * q);
            this.n = n;
            this.b = 1.15 + 2.53 * spq;
            this.a = -0.0873 + 0.0248 * b + 0.01 * p;
            this.c = n * p + 0.5;
            this.vr = 0.92 - 4.2 / b;
            this.alpha = (2.83 + 5.1 / b) * spq;
            this.logRatio = Math.log(p / q);
            this.logPmfAtMode = logPmf(Math.floor((n + 1.0) * p), n, logRatio);
        }

        int sample(PseudoRandom prng) {
            while (true) {
                double u = prng.nextDouble() - 0.5;
                double v = prng.nextDouble();
                double us = 0.5 - Math.abs(u);
                double kd = Math.floor((2.0 * a / us + b) * u + c);

                if (kd < 0.0 || kd > n) {
                    continue;
                }
                if (us >= 0.07 && v <= vr) {
                    return (int) kd;
                }
                if (Math.log(v * alpha / (a / (us * us) + b)) <= logPmf(kd, n, logRatio) - logPmfAtMode) {
                    return (int) kd;
                }
            }
        }
    }

    /**
     * The logarithm of the binomial mass function at {@code k}, up to the
     * constant {@code log(n!) + n log(1-p)} that cancels in every ratio taken
     * of it here.
     */
    private static double logPmf(double k, int n, double logRatio) {
        return k * logRatio - FastGamma.logGamma(k + 1.0) - FastGamma.logGamma(n - k + 1.0);
    }
}
