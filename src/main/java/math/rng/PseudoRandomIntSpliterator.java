/*
 * Copyright 2021, 2026 Stefan Zobel
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
        return new PseudoRandomIntSpliterator(half, idx, s, min, max);
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
            PseudoRandom pr = prng;
            int min_ = min;
            int max_ = max;
            do {
                consumer.accept(pr.nextInt(min_, max_));
            } while (++idx < fence_);
        }
    }
}
