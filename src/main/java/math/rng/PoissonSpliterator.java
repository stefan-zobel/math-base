package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

import math.cern.FastGamma;

final class PoissonSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    /**
     * The mean above which the transformed rejection is used instead of
     * multiplying uniforms.
     * <p>
     * Both are exact, so this decides speed and nothing else: multiplication
     * costs about {@code lambda + 1} uniforms and rejection a constant few.
     * <p>
     * The literature puts this at ten. Measured here it belongs near thirty --
     * at {@code lambda = 20} multiplication takes 41 nanoseconds a draw against
     * rejection's 60, and at 40 it is 71 against 49 -- because the acceptance
     * below compares against a log-gamma rather than a Stirling expansion, which
     * buys a shorter file at the price of a slower rare path. Inversion is cheap
     * enough to carry the extra range.
     */
    static final double REJECTION_THRESHOLD = 30.0;

    /**
     * The largest mean that is offered. A count is returned as an {@code int},
     * and a mean of a billion sits about seventy thousand standard deviations
     * below {@link Integer#MAX_VALUE}, so nothing this side of it can overflow.
     */
    static final double MAX_MEAN = 1.0e9;

    final double lambda;
    /** The derived constants, built once and shared with every split, or {@code null} below the threshold. */
    final Rejection rejection;
    final PseudoRandom prng;

    PoissonSpliterator(PseudoRandom prng, long index, long fence, double lambda) {
        super(index, fence);
        checkMean(lambda);
        this.lambda = lambda;
        this.rejection = lambda < REJECTION_THRESHOLD ? null : new Rejection(lambda);
        this.prng = prng;
    }

    private PoissonSpliterator(PseudoRandom prng, long index, long fence, double lambda, Rejection rejection) {
        super(index, fence);
        this.lambda = lambda;
        this.rejection = rejection;
        this.prng = prng;
    }

    static void checkMean(double lambda) {
        if (Double.isNaN(lambda) || lambda < 0.0) {
            throw new IllegalArgumentException("lambda < 0.0 : " + lambda);
        }
        if (lambda > MAX_MEAN) {
            throw new IllegalArgumentException("lambda > " + MAX_MEAN + " does not fit an int : " + lambda);
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
        // the constants are immutable and depend only on lambda, so the two
        // halves share them rather than deriving them again
        return new PoissonSpliterator(half, idx, s, lambda, rejection);
    }

    private int next(PseudoRandom pr) {
        if (lambda <= 0.0) {
            return 0;
        }
        return rejection == null ? byMultiplication(pr, lambda) : rejection.sample(pr);
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
     * One draw, deriving whatever constants it needs. The streams above do not
     * use this -- they keep the constants -- so it is here for a caller that
     * wants a single count, and for the tests.
     */
    static int sample(PseudoRandom prng, double lambda) {
        if (lambda <= 0.0) {
            return 0;
        }
        return lambda < REJECTION_THRESHOLD ? byMultiplication(prng, lambda)
                : new Rejection(lambda).sample(prng);
    }

    /**
     * Knuth's method: multiply uniforms until the product falls below
     * {@code exp(-lambda)}, and count how many it took.
     * <p>
     * Exact, and it costs about {@code lambda + 1} uniforms, which is why it is
     * only used for a small mean. It also needs {@code exp(-lambda)} to be a
     * normal number, which it stops being around {@code lambda = 745} -- far
     * above where the other branch takes over.
     */
    static int byMultiplication(PseudoRandom prng, double lambda) {
        double limit = Math.exp(-lambda);
        double product = prng.nextDouble();
        int k = 0;
        while (product > limit) {
            product *= prng.nextDouble();
            k++;
        }
        return k;
    }

    /** One draw by transformed rejection, deriving the constants first. */
    static int byTransformedRejection(PseudoRandom prng, double lambda) {
        return new Rejection(lambda).sample(prng);
    }

    /**
     * Hoermann's PTRS: transformed rejection with a hat that is a scaled
     * hyperbola, accepting in constant expected time.
     * <p>
     * The constants depend only on the mean, so they are derived once and kept.
     * Measured, deriving them per draw instead cost more than the draw: a
     * square root, a logarithm and four divisions against two uniforms.
     * <p>
     * The two cheap tests come first -- a wide central region accepted outright
     * and a narrow one rejected outright -- so the logarithm and the log-gamma
     * are reached by only a few percent of the draws. The acceptance compares
     * against the logarithm of the mass function itself rather than a Stirling
     * expansion of it, which costs one {@link FastGamma#logGamma} on that slow
     * path and keeps a table of correction terms out of this file.
     * <p>
     * <b>The constants were derived for {@code lambda >= 10}</b> and are not
     * merely inaccurate below it: at {@code lambda = 1} the acceptance bound
     * {@code vr} comes out negative, and the sampler then draws from something
     * that is not a Poisson distribution at all.
     */
    static final class Rejection {

        private final double lambda;
        private final double b;
        private final double a;
        private final double inverseAlpha;
        private final double vr;
        private final double logLambda;

        Rejection(double lambda) {
            this.lambda = lambda;
            this.b = 0.931 + 2.53 * Math.sqrt(lambda);
            this.a = -0.059 + 0.02483 * b;
            this.inverseAlpha = 1.1239 + 1.1328 / (b - 3.4);
            this.vr = 0.9277 - 3.6224 / (b - 2.0);
            this.logLambda = Math.log(lambda);
        }

        int sample(PseudoRandom prng) {
            while (true) {
                double u = prng.nextDouble() - 0.5;
                double v = prng.nextDouble();
                double us = 0.5 - Math.abs(u);
                double kd = Math.floor((2.0 * a / us + b) * u + lambda + 0.43);

                if (us >= 0.07 && v <= vr) {
                    return (int) kd;
                }
                if (kd < 0.0 || (us < 0.013 && v > us)) {
                    continue;
                }
                if (Math.log(v * inverseAlpha / (a / (us * us) + b)) <= kd * logLambda - lambda
                        - FastGamma.logGamma(kd + 1.0)) {
                    return (int) kd;
                }
            }
        }
    }
}
