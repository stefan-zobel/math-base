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
        return new GammaSpliterator(half, idx, s, shape_k, scale_theta);
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

    /**
     * Draws a {@code Gamma(shape_k, 1)} variate and returns its natural
     * logarithm without ever forming the variate itself.
     * <p>
     * For a shape below one, step 2 of Algorithm GS returns
     * {@code x = p^(1/shape_k)}, and that spans hundreds of decades: at
     * {@code shape_k = 0.001} more than nine of ten draws fall below the
     * smallest positive {@code double} and come back from {@link #sample} as an
     * exact zero. A caller that only needs a ratio of two such variates can
     * take the difference of their logarithms instead, and {@code log(p) /
     * shape_k} is that logarithm exactly.
     * <p>
     * The accept and reject decisions are the ones {@link #sample} makes, from
     * the same draws in the same order, so the two agree on which variate comes
     * next. The one exception is a {@code p} of zero, which {@link #sample}
     * accepts as a zero variate and this method has to reject, because its
     * logarithm is not a number the caller can subtract.
     *
     * @param prng
     *            the generator to draw from
     * @param shape_k
     *            the shape, greater than zero
     * @return the natural logarithm of a {@code Gamma(shape_k, 1)} variate
     */
    static double logSample(PseudoRandom prng, double shape_k) {
        if (shape_k < 1.0) {
            // [1]: p. 228, Algorithm GS
            final double bGS = 1.0 + shape_k / Math.E;

            while (true) {
                // Step 1:
                double u = prng.nextDouble();
                double p = bGS * u;

                if (p == 0.0) {
                    // log(0) is not a logarithm a difference can use
                    continue;
                }

                if (p <= 1.0) {
                    // Step 2: x is needed for the acceptance test only, and
                    // that test is unharmed when it underflows -- exp(-0) is
                    // one and accepts
                    double x = Math.pow(p, 1.0 / shape_k);
                    double u2 = prng.nextDouble();

                    if (u2 > Math.exp(-x)) {
                        // reject
                        continue;
                    } else {
                        return Math.log(p) / shape_k;
                    }
                } else {
                    // Step 3: x is of order one here, its logarithm is plain
                    double x = -1.0 * Math.log((bGS - p) / shape_k);
                    double u2 = prng.nextDouble();

                    if (u2 > Math.pow(x, shape_k - 1.0)) {
                        // reject
                        continue;
                    } else {
                        return Math.log(x);
                    }
                }
            }
        }

        // shape >= 1: the Marsaglia-Tsang branch cannot underflow far enough
        // to matter, so the variate can be formed and then taken apart
        return Math.log(sample(prng, shape_k, 1.0));
    }
}
