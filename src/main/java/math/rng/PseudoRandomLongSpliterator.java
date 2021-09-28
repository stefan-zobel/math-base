/*
 * Copyright 2021 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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
            PseudoRandom pr = prng;
            long min_ = min;
            long max_ = max;
            do {
                consumer.accept(pr.nextLong(min_, max_));
            } while (++idx < fence_);
        }
    }
}
