package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

final class PseudoRandomIntSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    final int min;
    final int max;
    final PseudoRandom prng;

    PseudoRandomIntSpliterator(PseudoRandom prng, long index, long fence, int min, int max) {
        super(index, fence);
        this.min = min;
        this.max = max;
        this.prng = prng;
    }

    @Override
    public Spliterator.OfInt trySplit() {
        long idx = index;
        long s = (idx + fence) >>> 1;
        if (s <= idx) {
            return null;
        }
        index = s;
        return new PseudoRandomIntSpliterator(prng, idx, s, min, max);
    }

    @Override
    public boolean tryAdvance(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(prng.nextInt(min, max));
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
            PseudoRandom r = prng;
            int min_ = min;
            int max_ = max;
            do {
                consumer.accept(r.nextInt(min_, max_));
            } while (++idx < fence_);
        }
    }
}
