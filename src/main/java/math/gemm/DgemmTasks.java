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

import java.util.Arrays;

final class DgemmTasks {

    static final class TaskConfig {
        int nb_lo = 0;
        int nb_hi = 1;
        boolean nb_hi_is_last = true;
        int kb_lo = 0;
        int kb_hi = 1;
        boolean kb_hi_is_last = true;
        int mb_lo = 0;
        int mb_hi = 1;
        boolean mb_hi_is_last = true;

        String loop = "";
        int tasks;

        @Override
        public String toString() {
            return "{nb_lo=" + nb_lo + ", nb_hi=" + nb_hi + " (last=" + nb_hi_is_last + "), kb_lo=" + kb_lo + ", kb_hi="
                    + kb_hi + " (last=" + kb_hi_is_last + "), mb_lo=" + mb_lo + ", mb_hi=" + mb_hi + " (last="
                    + mb_hi_is_last + ")}";
        }
    }

    private static final class LoopBounds {
        int lo = 0;
        int hi = 1;
        boolean hi_is_last = true;

        /* package */ LoopBounds() {}

        @Override
        public String toString() {
            return "for [" + lo + ".." + hi + ") (hi_is_last=" + hi_is_last + ")";
        }
    }

    static TaskConfig[] split(int nb, int kb, int mb) {
        int cores = availableCores();
        if (cores <= 1) {
            throw new AssertionError("#cores : " + cores);
        }

        LoopBounds[] loopsMb = divideLoopBounds(mb, cores);
        LoopBounds[] loopsNb = null;
        int threads = loopsMb.length;

        if (loopsMb.length <= 1) {
            // can't divide mb
            loopsNb = divideLoopBounds(nb, cores);
            threads = loopsNb.length;
            if (loopsNb.length <= 1) {
                throw new AssertionError("At least one loop must be divisible!");
            }
        }
        TaskConfig[] threadCfgs = new TaskConfig[threads];
        // initialization
        for (int i = 0; i < threadCfgs.length; ++i) {
            TaskConfig cfg = new TaskConfig();
            cfg.nb_hi = nb;
            cfg.kb_hi = kb;
            cfg.mb_hi = mb;
            threadCfgs[i] = cfg;
        }
        // correct task configs for the particular loop to be parallelized
        if (loopsMb != null && loopsMb.length > 1 && loopsMb.length == threadCfgs.length) {
            for (int i = 0; i < threadCfgs.length; ++i) {
                TaskConfig cfg = threadCfgs[i];
                cfg.mb_lo = loopsMb[i].lo;
                cfg.mb_hi = loopsMb[i].hi;
                cfg.mb_hi_is_last = loopsMb[i].hi_is_last;
                cfg.loop = "mb";
                cfg.tasks = loopsMb.length;
            }
        } else if (loopsNb.length > 1 && loopsNb.length == threadCfgs.length) {
            for (int i = 0; i < threadCfgs.length; ++i) {
                TaskConfig cfg = threadCfgs[i];
                cfg.nb_lo = loopsNb[i].lo;
                cfg.nb_hi = loopsNb[i].hi;
                cfg.nb_hi_is_last = loopsNb[i].hi_is_last;
                cfg.loop = "nb";
                cfg.tasks = loopsNb.length;
            }
        } else {
            throw new IllegalStateException(
                    "loopsNb: " + Arrays.toString(loopsNb) + ", loopsMb: " + Arrays.toString(loopsMb));
        }
        return threadCfgs;
    }

    private static LoopBounds[] divideLoopBounds(int count, int cores) {
        int div = count / cores;
        int rest = count % cores;
        if (div == 0 && rest == 1) {
            // can't divide
            int threads = 1;
            LoopBounds[] bounds = new LoopBounds[threads];
            LoopBounds bound = new LoopBounds();
            bounds[0] = bound;
            return bounds;
        } else if (div == 0 && rest > 1) {
            // can divide in less tasks than we have cores
            int threads = rest;
            LoopBounds[] bounds = new LoopBounds[threads];
            for (int i = 0; i < bounds.length; ++i) {
                LoopBounds bound = new LoopBounds();
                bound.lo = i;
                bound.hi = i + 1;
                bound.hi_is_last = false;
                if (i == bounds.length - 1) {
                    bound.hi_is_last = true;
                }
                bounds[i] = bound;
            }
            return bounds;
        } else {
            // div >= 1
            int threads = cores;
            LoopBounds[] bounds = new LoopBounds[threads];
            for (int i = 0; i < bounds.length; ++i) {
                LoopBounds bound = new LoopBounds();
                bound.lo = i * div;
                bound.hi = (i + 1) * div;
                bound.hi_is_last = false;
                if (i == bounds.length - 1) {
                    bound.hi_is_last = true;
                }
                bounds[i] = bound;
            }
            if (rest > 0) {
                // distribution of the remainder
                for (int i = 0; i < bounds.length; ++i) {
                    if (i == 0) {
                        bounds[i].hi += 1;
                    } else {
                        int delta = bounds[i].hi - bounds[i].lo;
                        bounds[i].lo = bounds[i - 1].hi;
                        if (rest > 0) {
                            bounds[i].hi = bounds[i].lo + delta + 1;
                        } else {
                            bounds[i].hi = bounds[i].lo + delta;
                        }
                    }
                    --rest;
                }
            }
            return bounds;
        }
    }

    static int availableCores() {
        int cores = Runtime.getRuntime().availableProcessors();
        if (cores <= 2) {
            return cores;
        }
        if (cores % 2 != 0) {
            return cores;
        }
        // use all cores (for now)
        return (cores / 2) + (cores / 2);
    }
}
