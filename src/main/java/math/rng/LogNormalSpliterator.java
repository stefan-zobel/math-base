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

final class LogNormalSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double mu;
    final double sigma;
    final PseudoRandom prng;

    LogNormalSpliterator(PseudoRandom prng, long index, long fence, double mu, double sigma) {
        super(index, fence);
        if (sigma <= 0.0) {
            throw new IllegalArgumentException("Standard deviation must be positive (" + sigma + ")");
        }
        this.mu = mu;
        this.sigma = sigma;
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
        return new LogNormalSpliterator(prng, idx, s, mu, sigma);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, mu, sigma));
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
            double mu_ = mu;
            double sigma_ = sigma;
            do {
                consumer.accept(sample(pr, mu_, sigma_));
            } while (++idx < fence_);
        }
    }

    private static double sample(PseudoRandom prng, double mu, double sigma) {
        double stdNormal = prng.nextGaussian();
        return Math.exp(mu + sigma * stdNormal);
    }
}
