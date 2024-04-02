/*
 * Copyright 2024 Stefan Zobel
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

final class InverseGammaSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double alpha;
    final double beta;
    final double inverse_scale_theta;
    final PseudoRandom prng;

    InverseGammaSpliterator(PseudoRandom prng, long index, long fence, double alpha, double beta) {
        super(index, fence);
        if (alpha <= 0.0) {
            throw new IllegalArgumentException("alpha <= 0.0 (" + alpha + ")");
        }
        if (beta <= 0.0) {
            throw new IllegalArgumentException("beta <= 0.0 (" + beta + ")");
        }
        this.alpha = alpha;
        this.beta = beta;
        this.inverse_scale_theta = 1.0 / beta;
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
        return new InverseGammaSpliterator(prng, idx, s, alpha, beta);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, alpha, inverse_scale_theta));
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
            double a = alpha;
            double theta_inv = inverse_scale_theta;
            do {
                consumer.accept(sample(pr, a, theta_inv));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng, double alpha, double inverse_scale_theta) {
        double gamma = 0.0;
        do {
            gamma = GammaSpliterator.sample(prng, alpha, inverse_scale_theta);
        } while (gamma == 0.0);
        return 1.0 / gamma;
    }
}
