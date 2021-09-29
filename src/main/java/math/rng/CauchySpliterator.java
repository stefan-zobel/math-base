/*
 * Copyright 2013, 2021 Stefan Zobel
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
import java.util.function.DoubleConsumer;

final class CauchySpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double location;
    final double scale;
    final PseudoRandom prng;

    CauchySpliterator(PseudoRandom prng, long index, long fence, double location, double scale) {
        super(index, fence);
        if (scale <= 0.0) {
            throw new IllegalArgumentException("scale <= 0.0 (" + scale + ")");
        }
        this.location = location;
        this.scale = scale;
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
        return new CauchySpliterator(prng, idx, s, location, scale);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, location, scale));
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
            PseudoRandom pr = prng;
            double loc_ = location;
            double sc_ = scale;
            do {
                consumer.accept(sample(pr, loc_, sc_));
            } while (++idx < fence_);
        }
    }

    private static double sample(PseudoRandom prng, double loc, double scale) {
        return loc + (scale * Math.tan(Math.PI * (prng.nextDouble() - 0.5)));
    }
}
