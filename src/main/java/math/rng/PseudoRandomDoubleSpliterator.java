package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.DoubleConsumer;

final class PseudoRandomDoubleSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double min;
    final double max;
    final PseudoRandom prng;

    PseudoRandomDoubleSpliterator(PseudoRandom prng, long index, long fence, double min, double max) {
        super(index, fence);
        this.min = min;
        this.max = max;
        this.prng = prng;
    }

    @Override
    public Spliterator.OfDouble trySplit() {
        long idx = index;
        long s = (idx + fence) >>> 1;
        if (s <= idx) {
            return null;
        }
        index = s;
        return new PseudoRandomDoubleSpliterator(prng, idx, s, min, max);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(prng.nextDouble(min, max));
            index = idx + 1;
            return true;
        } else {
            return false;
        }
    }

    @Override
    public void forEachRemaining(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            index = fence_;
            PseudoRandom r = prng;
            double min_ = min;
            double max_ = max;
            do {
                consumer.accept(r.nextDouble(min_, max_));
            } while (++idx < fence_);
        }
    }
}
