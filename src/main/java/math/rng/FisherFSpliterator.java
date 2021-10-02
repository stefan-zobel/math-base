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

final class FisherFSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final int d1;
    final int d2;
    final PseudoRandom prng_U;
    final PseudoRandom prng_V;

    FisherFSpliterator(PseudoRandom prng, long index, long fence, int numeratorDF, int denominatorDF) {
        super(index, fence);
        if (numeratorDF < 1) {
            throw new IllegalArgumentException("numeratorDF < 1 (" + numeratorDF + ")");
        }
        if (denominatorDF < 1) {
            throw new IllegalArgumentException("denominatorDF < 1 (" + denominatorDF + ")");
        }
        this.d1 = numeratorDF;
        this.d2 = denominatorDF;
        this.prng_U = prng;
        if (prng instanceof SplittablePseudoRandom) {
            this.prng_V = ((SplittablePseudoRandom) prng).split();
        } else {
            // this not only requires that 'prng' is an AbstractRng64 but
            // also that is has a public constructor taking a seed as a long
            this.prng_V = ((AbstractRng64) prng).newInstance();
        }
    }

    private FisherFSpliterator(PseudoRandom prng_u, PseudoRandom prng_v, long index, long fence, int numeratorDF,
            int denominatorDF) {
        super(index, fence);
        this.d1 = numeratorDF;
        this.d2 = denominatorDF;
        this.prng_U = prng_u;
        this.prng_V = prng_v;
    }

    @Override
    public Spliterator.OfDouble trySplit() {
        long idx = index;
        long s = (idx + fence) >>> 1;
        if (s <= idx) {
            return null;
        }
        index = s;
        return new FisherFSpliterator(prng_U, prng_V, idx, s, d1, d2);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng_U, prng_V, d1, d2));
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
            PseudoRandom pru = prng_U;
            PseudoRandom prv = prng_V;
            double numeratorDF = d1;
            double denominatorDF = d2;
            do {
                consumer.accept(sample(pru, prv, numeratorDF, denominatorDF));
            } while (++idx < fence_);
        }
    }

    private static double sample(PseudoRandom prng_U, PseudoRandom prng_V, double d1, double d2) {
        double y = BetaSpliterator.sample(prng_U, prng_V, d1 / 2.0, d2 / 2.0);
        return (y * d2) / (d1 - d1 * y);
    }
}
