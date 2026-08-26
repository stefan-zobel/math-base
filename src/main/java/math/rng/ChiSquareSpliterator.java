/*
 * Copyright 2013, 2026 Stefan Zobel
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

    final double degreesOfFreedom;
    final PseudoRandom prng;

    /** The parameter check, shared by the stream and the single draw. */
    static void checkDegreesOfFreedom(double degreesOfFreedom) {
        if (!(degreesOfFreedom > 0.0)) {
            throw new IllegalArgumentException(
                    "degreesOfFreedom must be greater than zero (" + degreesOfFreedom + ")");
        }
    }

    ChiSquareSpliterator(PseudoRandom prng, long index, long fence, double degreesOfFreedom) {
        super(index, fence);
        checkDegreesOfFreedom(degreesOfFreedom);
        this.degreesOfFreedom = degreesOfFreedom;
        this.prng = prng;
    }

    @Override
    public Spliterator.OfDouble trySplit() {
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
        return new ChiSquareSpliterator(half, idx, s, degreesOfFreedom);
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
