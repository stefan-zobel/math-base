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

final class StudentTSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double df;
    final PseudoRandom prng;

    StudentTSpliterator(PseudoRandom prng, long index, long fence, double df) {
        super(index, fence);
        if (df <= 0.0) {
            throw new IllegalArgumentException("df <= 0.0 : " + df);
        }
        this.df = df;
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
        return new StudentTSpliterator(prng, idx, s, df);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, df));
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
            double degrees = df;
            do {
                consumer.accept(sample(pr, degrees));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng, double df) {
        /*
         * Marsaglia's formulation of the Box/Muller polar method for generating
         * Normal variates is adapted to the Student-t distribution. The two
         * generated variates are not independent and the expected number of
         * uniforms per variate is 2.5464.
         * 
         * Reference:
         * 
         * R.W. Bailey (1994): Polar generation of random variates with the
         * t-distribution, Mathematics of Computation 62, 779-781.
         */
        double u1, u2, q;
        do {
            u1 = 2.0 * prng.nextDouble() - 1.0; // between -1 and 1
            u2 = 2.0 * prng.nextDouble() - 1.0; // between -1 and 1
            q = u1 * u1 + u2 * u2;
        } while (q > 1.0);
        return u1 * Math.sqrt(df * (Math.exp(-2.0 / df * Math.log(q)) - 1.0) / q);
    }
}
