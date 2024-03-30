/*
 * Copyright 2018 Stefan Zobel
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

import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

// simplest thing that can possibly work
final class ThreadPool {

    static Future<Integer> submit(Callable<Integer> computation) {
        return exec.submit(computation);
    }

    private static final ThreadPoolExecutor exec = new ThreadPoolExecutor(DgemmTasks.availableCores(),
            2 * DgemmTasks.availableCores(), 2, TimeUnit.MINUTES, new SynchronousQueue<Runnable>(),
            new ThreadFactory() {
                @Override
                public Thread newThread(Runnable r) {
                    Thread t = new Thread(r);
                    t.setName(ThreadPool.class.getName() + "-" + cnt.getAndIncrement());
                    t.setDaemon(true);
                    return t;
                }
            });

    /* package */ static final AtomicInteger cnt = new AtomicInteger(1);
}
