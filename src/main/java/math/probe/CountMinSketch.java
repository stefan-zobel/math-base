/*
 * Copyright 2026 Stefan Zobel
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
package math.probe;

import math.rng.BitMix;
import math.rng.Seed;
import math.rng.SplitMix64Seed;

/**
 * Count-Min Sketch implementation. Estimates frequencies of discrete elements
 * in a stream using minimal memory.
 *
 * @since 1.4.1
 */
public final class CountMinSketch<T> {

    private static final long NULL_HASH = Seed.get1Constant()[0];

    private final int rows; // depth
    private final int cols; // width
    private long totalCount = 0L;
    private final long table[];
    private final long[] hashSeeds;

    /**
     * @param depth
     *            Number of hash functions (rows)
     * @param width
     *            Size of the counter arrays (columns)
     */
    public CountMinSketch(int depth, int width) {
        this.rows = depth;
        this.cols = width;
        this.table = new long[checkArrayLength(depth, width)];
        this.hashSeeds = new long[depth];
        for (int i = 0; i < depth; ++i) {
            hashSeeds[i] = SplitMix64Seed.seed();
        }
    }

    /**
     * Increments the frequency of a T element.
     */
    public void add(T item) {
        totalCount++;
        for (int row = 0; row < rows; ++row) {
            int col = hash(item, row);
            table[idx(row, col)]++;
        }
    }

    /**
     * Estimates the frequency of an item. Guaranteed to be >= actual frequency.
     */
    public long estimateCount(T item) {
        long min = Long.MAX_VALUE;
        for (int row = 0; row < rows; ++row) {
            int col = hash(item, row);
            min = Math.min(min, table[idx(row, col)]);
        }
        return min;
    }

    public long getTotalCount() {
        return totalCount;
    }

    // A simple version of a universal hash
    private int hash(T item, int row) {
        long itemHash = (item == null) ? NULL_HASH : item.hashCode();
        long h = ((itemHash == 0L) ? NULL_HASH : itemHash) ^ hashSeeds[row];
        int col = (int) (BitMix.rrmxmx(h) % cols);
        return Math.abs(col);
    }

    private final int idx(int row, int col) {
        return col * rows + row;
    }

    private static int checkArrayLength(int rows, int cols) {
        long length = (long) checkRows(rows) * (long) checkCols(cols);
        if (length > (long) Integer.MAX_VALUE) {
            throw new IllegalArgumentException(
                    "rows x cols (= " + length + ") exceeds the maximal possible length (= 2147483647) of an array");
        }
        return (int) length;
    }

    private static int checkRows(int rows) {
        if (rows <= 0) {
            throw new IllegalArgumentException("number of rows must be strictly positive : " + rows);
        }
        return rows;
    }

    private static int checkCols(int cols) {
        if (cols <= 0) {
            throw new IllegalArgumentException("number of columns must be strictly positive : " + cols);
        }
        return cols;
    }
}
