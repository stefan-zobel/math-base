package math.rng;

import java.util.Spliterator;

abstract class PseudoRandomSpliterator {

    long index;
    final long fence;

    PseudoRandomSpliterator(long index, long fence) {
        this.index = index;
        this.fence = fence;
    }

    public long estimateSize() {
        return fence - index;
    }

    public int characteristics() {
        return (Spliterator.SIZED | Spliterator.SUBSIZED | Spliterator.NONNULL | Spliterator.IMMUTABLE);
    }
}
