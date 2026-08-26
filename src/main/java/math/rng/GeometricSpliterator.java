package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

final class GeometricSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    /**
     * The smallest success probability that is offered.
     * <p>
     * The count is the number of failures before the first success and comes
     * back as an {@code int}, and the distribution has an exponential tail:
     * the chance of a draw beyond {@link Integer#MAX_VALUE} is
     * {@code (1-p)^(2^31)}, which at this bound is about {@code 1e-94}. Below
     * it the tail stops being negligible -- at {@code p = 1e-11} nearly every
     * draw would overflow.
     */
    static final double MIN_PROBABILITY = 1.0e-7;

    final double p;
    /** {@code 1 / log(1-p)}, the only derived quantity, kept rather than taken per draw. */
    final double inverseLogQ;
    final PseudoRandom prng;

    GeometricSpliterator(PseudoRandom prng, long index, long fence, double p) {
        super(index, fence);
        checkProbability(p);
        this.p = p;
        this.inverseLogQ = p >= 1.0 ? 0.0 : 1.0 / Math.log1p(-p);
        this.prng = prng;
    }

    static void checkProbability(double p) {
        if (Double.isNaN(p) || p > 1.0) {
            throw new IllegalArgumentException("p must be in [" + MIN_PROBABILITY + ", 1] : " + p);
        }
        if (p < MIN_PROBABILITY) {
            throw new IllegalArgumentException(
                    "p < " + MIN_PROBABILITY + " does not fit an int : " + p);
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
        return new GeometricSpliterator(half, idx, s, p);
    }

    @Override
    public boolean tryAdvance(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, inverseLogQ));
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
            double factor = inverseLogQ;
            do {
                consumer.accept(sample(pr, factor));
            } while (++idx < fence_);
        }
    }

    /**
     * Inversion, which for this distribution is exact and closed:
     * {@code P[K >= k]} is {@code (1-p)^k}, so {@code floor(log(U) / log(1-p))}
     * has the right law however small the uniform comes out.
     * <p>
     * The logarithm of {@code 1-p} is taken with {@code log1p}, which keeps its
     * digits for a small {@code p} where {@code log(1 - p)} would have lost
     * them to the subtraction.
     */
    static int sample(PseudoRandom prng, double inverseLogQ) {
        if (inverseLogQ == 0.0) {
            // p == 1: the first trial always succeeds, so there are no failures
            return 0;
        }
        double u;
        do {
            u = prng.nextDouble();
        } while (u == 0.0);
        double k = Math.floor(Math.log(u) * inverseLogQ);
        return k >= Integer.MAX_VALUE ? Integer.MAX_VALUE : (int) k;
    }

    /** One draw for the given probability, deriving the factor first. */
    static int sampleFor(PseudoRandom prng, double p) {
        return sample(prng, p >= 1.0 ? 0.0 : 1.0 / Math.log1p(-p));
    }
}
