/*
 * Copyright 2017 Stefan Zobel
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
package math.distribution;

/**
 * The uniform distribution.
 * <p>
 * <b>See</b> <a
 * href="https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)">Wikipedia
 * Uniform distribution</a>.
 */
public class Uniform implements ContinuousDistribution {

    private final double a;
    private final double b;

    public Uniform(double a, double b) {
        if (b <= a) {
            throw new IllegalArgumentException("b <= a");
        }
        this.a = a;
        this.b = b;
    }

    public Uniform() {
        this(0.0, 1.0);
    }

    @Override
    public double pdf(double x) {
        if (x < a || x > b) {
            return 0.0;
        }
        return 1.0 / (b - a);
    }

    @Override
    public double cdf(double x) {
        if (x <= a) {
            return 0.0;
        }
        if (x >= b) {
            return 1.0;
        }
        return (x - a) / (b - a);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return a;
        }
        if (probability >= 1.0) {
            return b;
        }
        return a + probability * (b - a);
    }

    @Override
    public double mean() {
        return (a + b) / 2.0;
    }

    @Override
    public double variance() {
        return ((b - a) * (b - a)) / 12.0;
    }
}
