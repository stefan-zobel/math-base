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

/**
 * {@code 4} daisy-chained {@link Stc64} generators that are called in
 * round-robin fashion, i.e., first {@code gen_0}, then {@code gen_1}, then
 * {@code gen_2}, then {@code gen_3}, then again {@code gen_0} and so on. This
 * should still be a little faster than a {@link MersenneTwister64} but is a lot
 * slower (about {@code 2.5x}) than working with a single {@code Stc64}
 * generator. The generators following {@code gen_0} are constructed by calling
 * {@link SplittablePseudoRandom#split()} on the respective predecessor
 * generator.
 */
public class Interleaved4Stc64 extends AbstractRng64 implements SplittablePseudoRandom {

    private static final int SIZE = 4;
    private static final Interleaved4Stc64 defaultRng = new Interleaved4Stc64();
    private int pos = 0;
    private final Stc64[] prng = new Stc64[SIZE];

    public Interleaved4Stc64() {
        this(new Stc64());
    }

    public Interleaved4Stc64(long seed) {
        this(new Stc64(seed));
    }

    public Interleaved4Stc64(long[] seed) {
        this(new Stc64(seed));
    }

    private Interleaved4Stc64(Stc64 gen0) {
        this(gen0, gen0.split(), gen0.split(), gen0.split());
    }

    protected Interleaved4Stc64(Stc64 gen0, Stc64 gen1, Stc64 gen2, Stc64 gen3) {
        prng[0] = gen0;
        prng[1] = gen1;
        prng[2] = gen2;
        prng[3] = gen3;
        long[] seed0 = gen0.getSeed();
        long[] seed = new long[4 * seed0.length];
        System.arraycopy(seed0, 0, seed, 0, seed0.length);
        System.arraycopy(gen1.getSeed(), 0, seed, seed0.length, seed0.length);
        System.arraycopy(gen2.getSeed(), 0, seed, 2 * seed0.length, seed0.length);
        System.arraycopy(gen3.getSeed(), 0, seed, 3 * seed0.length, seed0.length);
        saveSeed(seed);
    }

    @Override
    public final long nextLong() {
        int idx = pos;
        long rnd = prng[idx].nextLong();
        if (idx == SIZE - 1) {
            pos = 0;
        } else {
            pos = idx + 1;
        }
        return rnd;
    }

    @Override
    public Interleaved4Stc64 split() {
        return new Interleaved4Stc64(prng[0].split(), prng[1].split(), prng[2].split(), prng[3].split());
    }

    public static SplittablePseudoRandom getDefault() {
        return defaultRng;
    }
}
