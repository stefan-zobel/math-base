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
        // this distribution draws from two independent generators
        PseudoRandom second = detach(prng);
        if (second == null) {
            throw new IllegalArgumentException("cannot derive the second generator a beta distribution needs from "
                    + prng.getAlgorithm() + ": it is neither splittable nor re-creatable from its current state");
        }
        this.prng_V = second;
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
        long s = splitPoint();
        if (s < 0L) {
            return null;
        }
        PseudoRandom halfU = detach(prng_U);
        PseudoRandom halfV = detach(prng_V);
        if (halfU == null || halfV == null) {
            // the source cannot hand out independent generators
            return null;
        }
        long idx = index;
        index = s;
        return new BetaSpliterator(halfU, halfV, idx, s, alpha, beta);
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

    /**
     * Draws from {@code Beta(alpha, beta)} through the gamma ratio
     * {@code X / (X + Y)}, computed as {@code 1 / (1 + exp(logY - logX))}.
     * <p>
     * The identity is exact, but forming {@code X} and {@code Y} first is not:
     * for a shape below one a {@code Gamma(k, 1)} draw spans hundreds of
     * decades and falls off the bottom of the {@code double} range. One
     * underflowed draw drags the ratio to an exact zero where the true value is
     * representable -- at {@code alpha = beta = 0.005} that doubled the rate of
     * exact zeros, 0.0236 measured against the 0.0121 the distribution has --
     * and two underflowed draws made it {@code 0/0}, which is where the 22 per
     * cent of {@code NaN} at {@code alpha = beta = 0.001} came from. Taking the
     * difference of the logarithms instead keeps the lower tail down to
     * {@code 4.9e-324}, exactly as far as a {@code double} reaches.
     * <p>
     * What no arrangement of the arithmetic can buy back is the upper end. A
     * {@code double} resolves {@code 4.9e-324} above zero but only
     * {@code 1.1e-16} below one, so for small shapes a large share of the
     * result is exactly {@code 1.0} -- 34.6 per cent at
     * {@code alpha = beta = 0.01}. That is the correctly rounded answer and not
     * a defect: the distribution really does put that much mass within half an
     * ulp of one. It does mean that a caller who needs the upper tail should
     * swap the two shapes and read the lower one, because {@code 1 - X} has
     * nothing left to say.
     *
     * @param prng_U
     *            the generator for the first gamma variate
     * @param prng_V
     *            an independent generator for the second
     * @param alpha
     *            the first shape, greater than zero
     * @param beta
     *            the second shape, greater than zero
     * @return a {@code Beta(alpha, beta)} variate
     */
    static double sample(PseudoRandom prng_U, PseudoRandom prng_V, double alpha, double beta) {
        double logX = GammaSpliterator.logSample(prng_U, alpha);
        double logY = GammaSpliterator.logSample(prng_V, beta);
        double d = logY - logX;
        // stay on whichever side of one half the result falls: for d >= 0 it is
        // exp(-d) all the way down to the smallest subnormal
        if (d >= 0.0) {
            double e = Math.exp(-d);
            return e / (1.0 + e);
        }
        return 1.0 / (1.0 + Math.exp(d));
    }
}
