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

final class WeibullSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double scale_lambda;
    final double shape_k;
    final PseudoRandom prng;

    /** The parameter check, shared by the stream and the single draw. */
    static void checkParameters(double scale, double shape) {
        if (!(scale > 0.0)) {
            throw new IllegalArgumentException("scale must be greater than zero (" + scale + ")");
        }
        if (!(shape > 0.0)) {
            throw new IllegalArgumentException("shape must be greater than zero (" + shape + ")");
        }
    }

    WeibullSpliterator(PseudoRandom prng, long index, long fence, double scale, double shape) {
        super(index, fence);
        checkParameters(scale, shape);
        this.scale_lambda = scale;
        this.shape_k = shape;
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
        return new WeibullSpliterator(half, idx, s, scale_lambda, shape_k);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, scale_lambda, shape_k));
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
            double scale = scale_lambda;
            double shape = shape_k;
            do {
                consumer.accept(sample(pr, scale, shape));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng, double scale, double shape) {
        double u = prng.nextDouble();
        if (u == 0.0 || u == 1.0) {
            return 0.0;
        }
        return scale * Math.pow(-Math.log1p(-u), 1.0 / shape);
    }
}
