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
package math.gemm;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

final class GemmParallelSupport {

    static final int MAX_PARALLEL_TASKS = 4;

    static int taskCount(int columns, int blockSize) {
        if (columns <= 0) {
            return 1;
        }
        int available = Math.max(1, Math.min(MAX_PARALLEL_TASKS, Runtime.getRuntime().availableProcessors()));
        int blocks = ceilDiv(columns, blockSize);
        return Math.max(1, Math.min(available, blocks));
    }

    static int taskCountUncapped(int columns, int blockSize) {
        if (columns <= 0) {
            return 1;
        }
        int available = Math.max(1, Runtime.getRuntime().availableProcessors());
        int blocks = ceilDiv(columns, blockSize);
        return Math.max(1, Math.min(available, blocks));
    }

    static int blockStart(int taskIndex, int taskCount, int columns, int blockSize) {
        int blocks = ceilDiv(columns, blockSize);
        int blockFrom = (taskIndex * blocks) / taskCount;
        return blockFrom * blockSize;
    }

    static int blockEnd(int taskIndex, int taskCount, int columns, int blockSize) {
        int blocks = ceilDiv(columns, blockSize);
        int blockTo = ((taskIndex + 1) * blocks) / taskCount;
        return Math.min(columns, blockTo * blockSize);
    }

    static void awaitAll(List<? extends Future<?>> futures) {
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException(e);
            } catch (ExecutionException e) {
                Throwable cause = e.getCause();
                if (cause instanceof RuntimeException) {
                    throw (RuntimeException) cause;
                }
                if (cause instanceof Error) {
                    throw (Error) cause;
                }
                throw new RuntimeException(cause);
            }
        }
    }

    private static int ceilDiv(int value, int divisor) {
        return (value + divisor - 1) / divisor;
    }

    private GemmParallelSupport() {
        throw new AssertionError();
    }
}
