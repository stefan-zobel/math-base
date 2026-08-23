/*
 * Copyright 2021, 2026 Stefan Zobel
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

import java.util.Spliterator;

abstract class PseudoRandomSpliterator {

    long index;
    final long fence;

    PseudoRandomSpliterator(long index, long fence) {
        this.index = index;
        this.fence = fence;
    }

    public long estimateSize() {
        return fence - index;
    }

    public int characteristics() {
        return (Spliterator.SIZED | Spliterator.SUBSIZED | Spliterator.NONNULL | Spliterator.IMMUTABLE);
    }

    /**
     * Returns the index at which this spliterator would be halved, or
     * {@code -1} if it is too small to be halved. Leaves {@code index}
     * untouched: a caller that goes on to split must assign the returned value
     * itself, and only once it has secured a generator for the new half.
     *
     * @return the middle index, or {@code -1} if this spliterator is too small
     */
    final long splitPoint() {
        long mid = (index + fence) >>> 1;
        return (mid <= index) ? -1L : mid;
    }

    /**
     * Returns a generator that shares no mutable state with {@code prng}, for
     * one half of a split, or {@code null} if none can be had. The two halves
     * of a split are handed to different threads and none of the generators in
     * this package is thread safe, so a spliterator that cannot obtain an
     * independent generator must refuse to split.
     *
     * @param prng
     *            the generator this spliterator draws from
     * @return an independent generator, or {@code null} if {@code prng} can
     *         neither be split nor re-created from its current state
     */
    static PseudoRandom detach(PseudoRandom prng) {
        if (prng instanceof SplittablePseudoRandom) {
            return ((SplittablePseudoRandom) prng).split();
        }
        if (prng instanceof AbstractRng64) {
            try {
                // derives a fresh generator from the current state; needs a
                // constructor taking a long, which not every implementation has
                return ((AbstractRng64) prng).newInstance();
            } catch (RuntimeException e) {
                return null;
            }
        }
        return null;
    }
}
