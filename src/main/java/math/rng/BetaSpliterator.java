/*
 * Copyright 2014, 2021 Stefan Zobel
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

final class BetaSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double alpha;
    final double beta;
    final PseudoRandom prng_U;
    final PseudoRandom prng_V;

    BetaSpliterator(PseudoRandom prng, long index, long fence, double alpha, double beta) {
        super(index, fence);
        if (alpha <= 0.0) {
            throw new IllegalArgumentException("alpha <= 0.0 (" + alpha + ")");
        }
        if (beta <= 0.0) {
            throw new IllegalArgumentException("beta <= 0.0 (" + beta + ")");
        }
        this.alpha = alpha;
        this.beta = beta;
        this.prng_U = prng;
        if (prng instanceof SplittablePseudoRandom) {
            this.prng_V = ((SplittablePseudoRandom) prng).split();
        } else {
            // this not only requires that 'prng' is an AbstractRng64 but
            // also that is has a public constructor taking a seed as a long
            this.prng_V = ((AbstractRng64) prng).newInstance();
        }
    }

    private BetaSpliterator(PseudoRandom prng_u, PseudoRandom prng_v, long index, long fence, double alpha,
            double beta) {
        super(index, fence);
        this.alpha = alpha;
        this.beta = beta;
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
        return new BetaSpliterator(prng_U, prng_V, idx, s, alpha, beta);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng_U, prng_V, alpha, beta));
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
            double alpha_ = alpha;
            double beta_ = beta;
            do {
                consumer.accept(sample(pru, prv, alpha_, beta_));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng_U, PseudoRandom prng_V, double alpha, double beta) {
        // This may not be the most efficient solution,
        // but it doesn't get any simpler. The problem is
        // alpha and beta must not be too small, especially
        // a beta < 1 paired with a very large alpha is numerically
        // inaccurate. But this seems to be true for all algorithms
        // (commons.math appears to be even more inaccurate than this
        // simple implementation - not to mention that it is much slower)
        //
        // An alpha and/or beta of 0.125 (1/8) should be ok, values below are
        // not. If you need to have a beta in the range 1/8 <= beta < 1 then
        // alpha must not be too large. A ratio of beta : alpha of 1 : 1800
        // should be ok (e.g. alpha = 225 for beta = 0.125). Don't go above
        // that. If only alpha is small (but not less than code 1/8) there seems
        // to exist no practically relevant limit for the magnitude of beta
        // (other than the lower bound of 1/8).
        double u = GammaSpliterator.sample(prng_U, alpha, 1.0);
        return u / (u + GammaSpliterator.sample(prng_V, beta, 1.0));
    }
}
