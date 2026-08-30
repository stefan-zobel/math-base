package math.rng;

import java.util.Objects;
import java.util.Spliterator;
import java.util.function.IntConsumer;

final class CategoricalSpliterator extends PseudoRandomSpliterator implements Spliterator.OfInt {

    // built once for the whole stream and shared with every half of every
    // split, which is what an immutable table is for: the other spliterators
    // in this package carry primitive parameters and have nothing to share
    final AliasTable table;
    final PseudoRandom prng;

    CategoricalSpliterator(PseudoRandom prng, long index, long fence, AliasTable table) {
        super(index, fence);
        this.table = table;
        this.prng = prng;
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
        return new CategoricalSpliterator(half, idx, s, table);
    }

    @Override
    public boolean tryAdvance(IntConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(table.draw(prng));
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
            AliasTable table_ = table;
            do {
                consumer.accept(table_.draw(pr));
            } while (++idx < fence_);
        }
    }
}
