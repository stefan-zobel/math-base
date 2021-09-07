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
 * Abstract base class for {@link XorShiftRot256PlusPlus} and
 * {@link XorShiftRot256StarStar}.
 */
abstract class Xoshiro256 extends AbstractRng64 {

    protected long x0;
    protected long x1;
    protected long x2;
    protected long x3;

    public Xoshiro256() {
        x0 = SplitMix64Seed.seed();
        x1 = SplitMix64Seed.seed();
        x2 = SplitMix64Seed.seed();
        x3 = SplitMix64Seed.seed();
        escape();
    }

    public Xoshiro256(long seed) {
        this(new XorShift64Star(seed));
    }

    public Xoshiro256(long[] seed) {
        this(new XorShift64Star(seed));
    }

    private Xoshiro256(XorShift64Star seeder) {
        x0 = seeder.nextLong();
        x1 = seeder.nextLong();
        x2 = seeder.nextLong();
        x3 = seeder.nextLong();
        escape();
    }

    private void escape() {
        saveSeed(new long[] { x0, x1, x2, x3 });
        long l = 0L;
        for (int i = 0; i < 20; ++i) {
            l = nextLong();
        }
        if (l == 0L) {
            unused = (byte) l;
        }
    }
}
