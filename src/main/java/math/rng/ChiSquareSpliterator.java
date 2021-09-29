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

final class ChiSquareSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    private final double degreesOfFreedom;
    final PseudoRandom prng;

    ChiSquareSpliterator(PseudoRandom prng, long index, long fence, double degreesOfFreedom) {
        super(index, fence);
        if (degreesOfFreedom <= 0.0) {
            throw new IllegalArgumentException("degreesOfFreedom <= 0.0 (" + degreesOfFreedom + ")");
        }
        this.degreesOfFreedom = degreesOfFreedom;
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
        return new ChiSquareSpliterator(prng, idx, s, degreesOfFreedom);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, degreesOfFreedom));
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
            double df = degreesOfFreedom;
            do {
                consumer.accept(sample(pr, df));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng, double degreesOfFreedom) {
        return GammaSpliterator.sample(prng, degreesOfFreedom / 2.0, 2.0);
    }
}
