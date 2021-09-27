package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.LongConsumer;

final class PseudoRandomLongSpliterator extends PseudoRandomSpliterator implements Spliterator.OfLong {

    final long min;
    final long max;
    final PseudoRandom prng;

    PseudoRandomLongSpliterator(PseudoRandom prng, long index, long fence, long min, long max) {
        super(index, fence);
        this.min = min;
        this.max = max;
        this.prng = prng;
    }

    @Override
    public Spliterator.OfLong trySplit() {
        long idx = index;
        long s = (idx + fence) >>> 1;
        if (s <= idx) {
            return null;
        }
        index = s;
        return new PseudoRandomLongSpliterator(prng, idx, s, min, max);
    }

    @Override
    public boolean tryAdvance(LongConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(prng.nextLong(min, max));
            index = idx + 1;
            return true;
        } else {
            return false;
        }
    }

    @Override
    public void forEachRemaining(LongConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            index = fence_;
            PseudoRandom r = prng;
            long min_ = min;
            long max_ = max;
            do {
                consumer.accept(r.nextLong(min_, max_));
            } while (++idx < fence_);
        }
    }
}
