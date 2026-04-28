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
package math.linalg;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import math.gemm.DgemmMRxNR;

/**
 * Shared thread-pool executor for matrix multiplication in {@link DMatrix} and
 * {@link FMatrix}.
 *
 * <h3>Thread-count rationale</h3>
 * The two gemm back-ends impose different parallelism bounds:
 * <ul>
 *   <li><b>Baseline</b> ({@code DgemmBaseline} / {@code SgemmBaseline}):
 *       internally caps the number of column-range tasks at
 *       {@code MAX_PARALLEL_TASKS = 4} via {@code GemmParallelSupport.taskCount()}.
 *       The work threshold is {@code m*n*k >= 1_000_000}.</li>
 *   <li><b>Lehn / BLIS-style MRxNR</b> ({@code DgemmMRxNR} / {@code SgemmMRxNR}):
 *       uses all logical processors ({@code GemmParallelSupport.taskCountUncapped()}).
 *       The work threshold is {@code m*n*k >= 25_000_000}.</li>
 * </ul>
 *
 * <h3>Conditional pool sizing</h3>
 * {@code GemmSwitch} only routes to the Lehn/MRxNR path when <em>both</em>
 * conditions hold at runtime:
 * <ol>
 *   <li>The {@code jdk.incubator.vector} module is present
 *       ({@code --add-modules jdk.incubator.vector}), <em>and</em></li>
 *   <li>the MRxNR kernel reports itself as actually vectorized
 *       (i.e. we are running on Java 25 with the Vector API enabled).</li>
 * </ol>
 * When neither condition is met, {@code GemmSwitch} always selects Baseline,
 * which caps its own task count at 4. In that case a pool larger than 4
 * threads would be wasted. The pool size is therefore chosen as:
 * <ul>
 *   <li>{@code availableProcessors()} &mdash; when the vectorized MRxNR path
 *       is reachable (Java 25 + Vector API enabled)</li>
 *   <li>{@code min(4, availableProcessors())} &mdash; otherwise (only Baseline
 *       can ever run)</li>
 * </ul>
 *
 * <h3>Thread lifetime</h3>
 * Worker threads are created on demand (lazily) and terminate automatically
 * after {@value #KEEP_ALIVE_SECONDS} seconds of idleness
 * ({@code allowCoreThreadTimeOut(true)}). They will be recreated transparently
 * on the next matrix multiplication that requires parallelism.
 * All worker threads are daemon threads so they do not prevent JVM shutdown.
 */
final class GemmExecutor {

    private static final long KEEP_ALIVE_SECONDS = 60L;

    private static final ExecutorService EXECUTOR = createExecutor();

    /**
     * Returns the shared executor. Never {@code null} and never shut down by
     * this class.
     */
    static ExecutorService get() {
        return EXECUTOR;
    }

    private static ExecutorService createExecutor() {
        int nThreads = computeThreadCount();
        ThreadPoolExecutor pool = new ThreadPoolExecutor(
                nThreads,                       // corePoolSize
                nThreads,                       // maximumPoolSize (fixed)
                KEEP_ALIVE_SECONDS, TimeUnit.SECONDS,
                new LinkedBlockingQueue<>(),    // unbounded queue (same as newFixedThreadPool)
                new DaemonThreadFactory("gemm-worker"));
        pool.allowCoreThreadTimeOut(true);      // idle threads die after keepAliveTime
        return pool;
    }

    /**
     * Determines the pool size based on whether the vectorized Lehn/MRxNR path
     * is reachable at runtime.
     * <p>
     * The MRxNR path uses {@code GemmParallelSupport.taskCountUncapped()} and
     * can therefore submit tasks up to {@code availableProcessors()}. It is
     * only reachable when the {@code jdk.incubator.vector} incubator module is
     * present <em>and</em> the MRxNR kernel is actually vectorized. If that
     * is not the case, {@code GemmSwitch} always routes to Baseline, which
     * internally caps its task count at {@code MAX_PARALLEL_TASKS = 4}.
     */
    private static int computeThreadCount() {
        int cpus = Math.max(1, Runtime.getRuntime().availableProcessors());
        boolean mrxnrVectorized = DgemmMRxNR.isVectorApiPresent() && DgemmMRxNR.isVectorized();
        if (mrxnrVectorized) {
            // Lehn/MRxNR path is reachable: pool must cover all logical cores.
            return cpus;
        }
        // Only Baseline can run: cap at 4 to avoid an oversized idle pool.
        return Math.min(4, cpus);
    }

    private static final class DaemonThreadFactory implements ThreadFactory {
        private final String namePrefix;
        private final AtomicInteger count = new AtomicInteger(0);

        DaemonThreadFactory(String namePrefix) {
            this.namePrefix = namePrefix;
        }

        @Override
        public Thread newThread(Runnable r) {
            Thread t = new Thread(r, namePrefix + "-" + count.incrementAndGet());
            t.setDaemon(true);
            t.setPriority(Thread.NORM_PRIORITY);
            return t;
        }
    }

    private GemmExecutor() {
        throw new AssertionError();
    }
}
