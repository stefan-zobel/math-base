/*
 * Copyright 2024, 2026 Stefan Zobel
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

    /** The parameter check, shared by the stream and the single draw. */
    static void checkParameters(double alpha, double beta) {
        if (!(alpha > 0.0)) {
            throw new IllegalArgumentException("alpha must be greater than zero (" + alpha + ")");
        }
        if (!(beta > 0.0)) {
            throw new IllegalArgumentException("beta must be greater than zero (" + beta + ")");
        }
    }

    /**
     * One draw for the scale as the caller states it. The core below wants the
     * reciprocal, which the stream derives once in its constructor.
     */
    static double sampleFor(PseudoRandom prng, double alpha, double beta) {
        checkParameters(alpha, beta);
        return sample(prng, alpha, 1.0 / beta);
    }

    InverseGammaSpliterator(PseudoRandom prng, long index, long fence, double alpha, double beta) {
        super(index, fence);
        checkParameters(alpha, beta);
        this.alpha = alpha;
        this.beta = beta;
        this.inverse_scale_theta = 1.0 / beta;
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
        return new InverseGammaSpliterator(half, idx, s, alpha, beta);
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
