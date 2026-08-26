package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

final class NegativeBinomialSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    final int r;
    final double p;
    /** The scale of the mixing gamma, {@code (1-p)/p}. */
    final double scale;
    final PseudoRandom prng;

    NegativeBinomialSpliterator(PseudoRandom prng, long index, long fence, int r, double p) {
        super(index, fence);
        checkParameters(r, p);
        this.r = r;
        this.p = p;
        this.scale = p >= 1.0 ? 0.0 : (1.0 - p) / p;
        this.prng = prng;
    }

    static void checkParameters(int r, double p) {
        if (r < 1) {
            throw new IllegalArgumentException("r < 1 : " + r);
        }
        if (Double.isNaN(p) || p <= 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in (0, 1] : " + p);
        }
        double mean = p >= 1.0 ? 0.0 : r * (1.0 - p) / p;
        if (mean > PoissonSpliterator.MAX_MEAN) {
            throw new IllegalArgumentException("the mean r (1-p) / p = " + mean + " > "
                    + PoissonSpliterator.MAX_MEAN + " does not fit an int");
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
        return new NegativeBinomialSpliterator(half, idx, s, r, p);
    }

    @Override
    public boolean tryAdvance(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, r, scale));
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
            int successes = r;
            double scale_ = scale;
            do {
                consumer.accept(sample(pr, successes, scale_));
            } while (++idx < fence_);
        }
    }

    /**
     * The gamma-Poisson mixture: a Poisson count whose own mean is drawn from
     * {@code Gamma(r, (1-p)/p)} has exactly the negative binomial law, the
     * number of failures before the {@code r}-th success.
     * <p>
     * That identity is why this file is short. Both halves already exist and
     * are tested in their own right, so nothing here is a third algorithm --
     * and the test checks the mixture against the definition it stands for,
     * counting failures against Bernoulli draws, rather than against itself.
     * <p>
     * The Poisson mean changes with every draw, so its constants cannot be kept
     * the way the other samplers keep theirs.
     */
    static int sample(PseudoRandom prng, int r, double scale) {
        if (scale == 0.0) {
            // p == 1: every trial succeeds, so there are no failures
            return 0;
        }
        double lambda = GammaSpliterator.sample(prng, r, scale);
        if (lambda >= PoissonSpliterator.MAX_MEAN) {
            // the mixing draw wandered past what an int can hold, which the
            // constructor made improbable rather than impossible
            lambda = PoissonSpliterator.MAX_MEAN;
        }
        return PoissonSpliterator.sample(prng, lambda);
    }

    /** One draw for the given parameters. */
    static int sampleFor(PseudoRandom prng, int r, double p) {
        return sample(prng, r, p >= 1.0 ? 0.0 : (1.0 - p) / p);
    }
}
