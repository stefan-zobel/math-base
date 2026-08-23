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

import math.MathConsts;
import math.cern.ProbabilityFuncs;

// https://jax.readthedocs.io/en/latest/_autosummary/jax.random.truncated_normal.html
final class TruncatedNormalSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double lower;
    final double upper;
    final double a;
    final double b;
    final PseudoRandom prng;

    TruncatedNormalSpliterator(PseudoRandom prng, long index, long fence, double lower, double upper) {
        super(index, fence);
        this.lower = lower;
        this.upper = upper;
        this.a = ProbabilityFuncs.errorFunction(lower / MathConsts.SQRT_TWO);
        if (upper == -lower) {
            this.b = -a;
        } else {
            this.b = ProbabilityFuncs.errorFunction(upper / MathConsts.SQRT_TWO);
        }
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
        return new TruncatedNormalSpliterator(half, idx, s, lower, upper);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(nextTruncatedNormal());
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
                consumer.accept(nextTruncatedNormal());
            } while (++idx < fence_);
        }
    }

    private double nextTruncatedNormal() {
        double u = prng.nextDouble(a, b);
        double out = MathConsts.SQRT_TWO * ProbabilityFuncs.errorFunctionInverse(u);
        // Clamp the value to the open interval (lower, upper) to make sure that
        // rounding doesn't push us outside of the range
        if (out == lower) {
            out = Math.nextUp(lower);
        }
        if (out == upper) {
            out = Math.nextDown(upper);
        }
        return out;
    }
}
