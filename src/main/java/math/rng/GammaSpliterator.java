/*
 * Copyright 2021 Stefan Zobel
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

final class GammaSpliterator extends PseudoRandomSpliterator implements Spliterator.OfDouble {

    final double shape_k;
    final double scale_theta;
    final PseudoRandom prng;

    GammaSpliterator(PseudoRandom prng, long index, long fence, double k, double theta) {
        super(index, fence);
        if (k <= 0.0) {
            throw new IllegalArgumentException("k <= 0.0 (" + k + ")");
        }
        if (theta <= 0.0) {
            throw new IllegalArgumentException("theta <= 0.0 (" + theta + ")");
        }
        this.shape_k = k;
        this.scale_theta = theta;
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
        return new GammaSpliterator(prng, idx, s, shape_k, scale_theta);
    }

    @Override
    public boolean tryAdvance(DoubleConsumer consumer) {
        Objects.requireNonNull(consumer);
        long idx = index;
        long fence_ = fence;
        if (idx < fence_) {
            consumer.accept(sample(prng, shape_k, scale_theta));
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
            double k = shape_k;
            double theta = scale_theta;
            do {
                consumer.accept(sample(pr, k, theta));
            } while (++idx < fence_);
        }
    }

    static double sample(PseudoRandom prng, double shape_k, double scale_theta) {
        if (shape_k < 1.0) {
            // [1]: p. 228, Algorithm GS
            final double bGS = 1.0 + shape_k / Math.E;

            while (true) {
                // Step 1:
                double u = prng.nextDouble();
                double p = bGS * u;

                if (p <= 1.0) {
                    // Step 2:

                    double x = Math.pow(p, 1.0 / shape_k);
                    double u2 = prng.nextDouble();

                    if (u2 > Math.exp(-x)) {
                        // reject
                        continue;
                    } else {
                        return scale_theta * x;
                    }
                } else {
                    // Step 3:

                    double x = -1.0 * Math.log((bGS - p) / shape_k);
                    double u2 = prng.nextDouble();

                    if (u2 > Math.pow(x, shape_k - 1.0)) {
                        // reject
                        continue;
                    } else {
                        return scale_theta * x;
                    }
                }
            }
        }

        // shape >= 1
        final double d = shape_k - 0.333333333333333333;
        final double c = 1.0 / (3.0 * Math.sqrt(d));

        while (true) {
            double x = prng.nextGaussian();
            double cx = 1.0 + c * x;
            double v = cx * cx * cx;

            if (v <= 0.0) {
                continue;
            }

            double x2 = x * x;
            double u = prng.nextDouble();

            // squeeze
            if (u < 1.0 - 0.0331 * x2 * x2) {
                return scale_theta * d * v;
            }

            if (Math.log(u) < 0.5 * x2 + d * (1.0 - v + Math.log(v))) {
                return scale_theta * d * v;
            }
        }
    }
}
