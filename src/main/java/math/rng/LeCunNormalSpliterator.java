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

import math.MathConsts;
import math.cern.ProbabilityFuncs;

// https://github.com/google/jax/blob/8815b236b656f494171131301d1d81e84cf4c67c/jax/_src/nn/initializers.py#L465
// https://github.com/google/jax/blob/8815b236b656f494171131301d1d81e84cf4c67c/jax/_src/nn/initializers.py#L265
// https://github.com/google/jax/blob/8815b236b656f494171131301d1d81e84cf4c67c/jax/_src/random.py#L818
// https://github.com/google/jax/blob/8815b236b656f494171131301d1d81e84cf4c67c/jax/_src/random.py#L861
final class LeCunNormalSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    // ProbabilityFuncs.errorFunction(-2.0 / MathConsts.SQRT_TWO)
    private static final double A = -0.9544997361036416;
    private static final double B = -A;

    final double sigma;
    final double stdDev;
    final PseudoRandom prng;

    LeCunNormalSpliterator(PseudoRandom prng, long index, long fence, double sigma) {
        super(index, fence);
        if (sigma <= 0.0) {
            throw new IllegalArgumentException("Standard deviation must be positive (" + sigma + ")");
        }
        this.sigma = sigma;
        // constant is stddev of standard normal truncated to (-2.0, 2.0)
        this.stdDev = sigma / 0.87962566103423978;
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
        return new LeCunNormalSpliterator(prng, idx, s, sigma);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(leCun());
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
            do {
                consumer.accept(leCun());
            } while (++idx < fence_);
        }
    }

    private double leCun() {
        double u = prng.nextDouble(A, B);
        double out = MathConsts.SQRT_TWO * ProbabilityFuncs.errorFunctionInverse(u);
        // Clamp the value to the open interval (-2.0, 2.0) to make sure that
        // rounding doesn't push us outside of the range
        if (out <= -2.0) {
            // Math.nextUp(-2.0)
            out = -1.9999999999999998;
        }
        if (out >= 2.0) {
            // Math.nextDown(2.0)
            out = 1.9999999999999998;
        }
        return stdDev * out;
    }
}
