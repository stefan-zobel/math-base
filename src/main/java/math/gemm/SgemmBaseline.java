package math.gemm;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/**
 * Single-precision dense matrix multiplication implementation for primitive
 * column-major {@code float[]} storage.
 * <p>
 * Supports all four transpose combinations of {@code A} and {@code B} and is
 * intended as a straightforward pure-Java baseline implementation for
 * correctness checks, benchmarking, and further experimentation.
 */
public final class SgemmBaseline {

    private static final int UNROLL = 8;
    private static final int NC = 2048;
    private static final long PARALLEL_WORK_THRESHOLD = 22_500_000L;

    /**
     * Performs the single-precision matrix-matrix operation
     * {@code C := alpha * op(A) * op(B) + beta * C} using the baseline
     * implementation.
     * 
     * @param notATransposed
     *            {@code true} if {@code A} is used as-is, {@code false} if
     *            {@code A} is treated as transposed
     * @param notBTransposed
     *            {@code true} if {@code B} is used as-is, {@code false} if
     *            {@code B} is treated as transposed
     * @param m
     *            the number of rows of op(A) and of the matrix {@code C}
     * @param n
     *            the number of columns of op(B) and the number of columns of
     *            the matrix {@code C}
     * @param k
     *            the number of columns of op(A) and the number of rows of
     *            op(B)
     * @param alpha
     *            scalar factor applied to {@code op(A) * op(B)}
     * @param a
     *            array containing matrix {@code A} in column-major storage
     * @param _a_offset
     *            offset into {@code a} where matrix {@code A} starts
     * @param lda
     *            leading dimension of matrix {@code A}
     * @param b
     *            array containing matrix {@code B} in column-major storage
     * @param _b_offset
     *            offset into {@code b} where matrix {@code B} starts
     * @param ldb
     *            leading dimension of matrix {@code B}
     * @param beta
     *            scalar factor applied to the existing contents of {@code C}
     * @param c
     *            array containing matrix {@code C} in column-major storage
     * @param _c_offset
     *            offset into {@code c} where matrix {@code C} starts
     * @param ldc
     *            leading dimension of matrix {@code C}
     */
    public static void sgemm(boolean notATransposed, boolean notBTransposed, int m, int n, int k, float alpha, float[] a, int _a_offset, int lda,
            float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {

        sgemm(notATransposed, notBTransposed, m, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc, null);
    }

    /**
     * Variant of
     * {@link #sgemm(boolean, boolean, int, int, int, float, float[], int, int, float[], int, int, float, float[], int, int)}
     * that optionally parallelizes the computation over column ranges of
     * {@code C} using the supplied {@link ExecutorService}.
     * <p>
     * If {@code executor} is {@code null}, shut down, or the internal workload
     * heuristic decides against parallel execution, the computation falls back
     * to the sequential path.
     * 
     * @param notATransposed
     *            {@code true} if {@code A} is used as-is, {@code false} if
     *            {@code A} is treated as transposed
     * @param notBTransposed
     *            {@code true} if {@code B} is used as-is, {@code false} if
     *            {@code B} is treated as transposed
     * @param m
     *            the number of rows of op(A) and of the matrix {@code C}
     * @param n
     *            the number of columns of op(B) and the number of columns of
     *            the matrix {@code C}
     * @param k
     *            the number of columns of op(A) and the number of rows of
     *            op(B)
     * @param alpha
     *            scalar factor applied to {@code op(A) * op(B)}
     * @param a
     *            array containing matrix {@code A} in column-major storage
     * @param _a_offset
     *            offset into {@code a} where matrix {@code A} starts
     * @param lda
     *            leading dimension of matrix {@code A}
     * @param b
     *            array containing matrix {@code B} in column-major storage
     * @param _b_offset
     *            offset into {@code b} where matrix {@code B} starts
     * @param ldb
     *            leading dimension of matrix {@code B}
     * @param beta
     *            scalar factor applied to the existing contents of {@code C}
     * @param c
     *            array containing matrix {@code C} in column-major storage
     * @param _c_offset
     *            offset into {@code c} where matrix {@code C} starts
     * @param ldc
     *            leading dimension of matrix {@code C}
     * @param executor
     *            executor used for optional parallel execution over the
     *            {@code N} dimension; may be {@code null}
     */
    public static void sgemm(boolean notATransposed, boolean notBTransposed, int m, int n, int k, float alpha, float[] a, int _a_offset, int lda,
            float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc, ExecutorService executor) {

        // Quick return if alpha == 0 or if op(A)*op(B) is empty
        if (alpha == 0.0f || k == 0) {
            if (beta == 0.0f) {
                int v = 1;
                for (int o = n; o > 0; o--) {
                    int cOffset = (v - 1) * ldc + _c_offset;
                    for (int i = 0; i < m; i++) {
                        c[i + cOffset] = 0.0f;
                    }
                    v++;
                }
            } else {
                int v = 1;
                for (int o = n; o > 0; o--) {
                    int cOffset = (v - 1) * ldc + _c_offset;
                    for (int i = 0; i < m; i++) {
                        c[i + cOffset] = beta * c[i + cOffset];
                    }
                    v++;
                }
            }
            return;
        }

        // Start the operations

        long work = (long) m * (long) n * (long) k;
        if (shouldParallelize(executor, n, work)) {
            parallelizeOverN(executor, notATransposed, notBTransposed, m, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c,
                    _c_offset, ldc);
            return;
        }

        sgemmRange(notATransposed, notBTransposed, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
    }

    private static boolean shouldParallelize(ExecutorService executor, int n, long work) {
        return executor != null && !executor.isShutdown() && work >= PARALLEL_WORK_THRESHOLD
                && GemmParallelSupport.taskCount(n, NC) >= 3;
    }

    private static void parallelizeOverN(ExecutorService executor, boolean notATransposed, boolean notBTransposed, int m, int n, int k,
            float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset,
            int ldc) {

        int taskCount = GemmParallelSupport.taskCount(n, NC);
        if (taskCount <= 1) {
            sgemmRange(notATransposed, notBTransposed, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset,
                    ldc);
            return;
        }

        ArrayList<Future<?>> futures = new ArrayList<>(taskCount);
        for (int taskIndex = 0; taskIndex < taskCount; taskIndex++) {
            int jFrom = GemmParallelSupport.blockStart(taskIndex, taskCount, n, NC);
            int jTo = GemmParallelSupport.blockEnd(taskIndex, taskCount, n, NC);
            futures.add(executor.submit(new ColumnRangeTask(notATransposed, notBTransposed, m, jFrom, jTo, k, alpha, a, _a_offset, lda,
                    b, _b_offset, ldb, beta, c, _c_offset, ldc)));
        }
        GemmParallelSupport.awaitAll(futures);
    }

    private static void sgemmRange(boolean notATransposed, boolean notBTransposed, int m, int jFrom, int jTo, int k, float alpha,
            float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {

        if (notBTransposed) {
            if (notATransposed) {
                // Form C := alpha*A*B + beta*C
                aTimesBRange(m, jFrom, jTo, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
            } else {
                // Form C := alpha*A**T*B + beta*C
                aTransTimesBRange(m, jFrom, jTo, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
            }
        } else if (notATransposed) {
            // Form C := alpha*A*B**T + beta*C
            aTimesBTransRange(m, jFrom, jTo, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
        } else {
            // Form C := alpha*A**T*B**T + beta*C
            aTransTimesBTransRange(m, jFrom, jTo, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
        }
    }

    // Form C := alpha*A*B + beta*C
    private static void aTimesB(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        sgemmRange(true, true, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
    }

    private static void aTimesBRange(int m, int jFrom, int jTo, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        final int KC = 256;
        final int MC = 128;

        for (int j = jFrom; j < jTo; j += NC) {
            int jb = Math.min(jTo - j, NC);
            for (int p = 0; p < k; p += KC) {
                int pb = Math.min(k - p, KC);
                // On the first K-tile we apply `beta`, on subsequent K-tiles we must strictly accumulate
                float currentBeta = (p == 0) ? beta : 1.0f;
                for (int i = 0; i < m; i += MC) {
                    int ib = Math.min(m - i, MC);
                    int aTileOffset = _a_offset + i + p * lda;
                    int bTileOffset = _b_offset + p + j * ldb;
                    int cTileOffset = _c_offset + i + j * ldc;
                    kernel_aTimesB(ib, jb, pb, alpha, a, aTileOffset, lda, b, bTileOffset, ldb, currentBeta, c, cTileOffset, ldc);
                }
            }
        }
    }

    // Inner block kernel
    private static void kernel_aTimesB(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        // Form C := alpha*A*B + beta*C
        final int mLimit = m - UNROLL + 1;
        int u = 1;
        for (int o = n; o > 0; o--) {
            int cOffset = (u - 1) * ldc + _c_offset;
            if (beta == 0.0f) {
                for (int i = 0; i < m; i++) {
                    c[i + cOffset] = 0.0f;
                }
            } else if (beta != 1.0f) {
                for (int i = 0; i < m; i++) {
                    c[i + cOffset] = beta * c[i + cOffset];
                }
            }
            int bColOffset = (u - 1) * ldb + _b_offset;
            int w = 1;
            int p = k;
            if (alpha == 1.0f) {
                for (; p >= 4; p -= 4) {
                    int aOffset0 = (w - 1) * lda + _a_offset;
                    int aOffset1 = w * lda + _a_offset;
                    int aOffset2 = (w + 1) * lda + _a_offset;
                    int aOffset3 = (w + 2) * lda + _a_offset;
                    float tmp0 = b[(w - 1) + bColOffset];
                    float tmp1 = b[w + bColOffset];
                    float tmp2 = b[(w + 1) + bColOffset];
                    float tmp3 = b[(w + 2) + bColOffset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp0 * a[i      + aOffset0] + tmp1 * a[i      + aOffset1] + tmp2 * a[i      + aOffset2] + tmp3 * a[i      + aOffset3];
                        c[i +  1 + cOffset] += tmp0 * a[i +  1 + aOffset0] + tmp1 * a[i +  1 + aOffset1] + tmp2 * a[i +  1 + aOffset2] + tmp3 * a[i +  1 + aOffset3];
                        c[i +  2 + cOffset] += tmp0 * a[i +  2 + aOffset0] + tmp1 * a[i +  2 + aOffset1] + tmp2 * a[i +  2 + aOffset2] + tmp3 * a[i +  2 + aOffset3];
                        c[i +  3 + cOffset] += tmp0 * a[i +  3 + aOffset0] + tmp1 * a[i +  3 + aOffset1] + tmp2 * a[i +  3 + aOffset2] + tmp3 * a[i +  3 + aOffset3];
                        c[i +  4 + cOffset] += tmp0 * a[i +  4 + aOffset0] + tmp1 * a[i +  4 + aOffset1] + tmp2 * a[i +  4 + aOffset2] + tmp3 * a[i +  4 + aOffset3];
                        c[i +  5 + cOffset] += tmp0 * a[i +  5 + aOffset0] + tmp1 * a[i +  5 + aOffset1] + tmp2 * a[i +  5 + aOffset2] + tmp3 * a[i +  5 + aOffset3];
                        c[i +  6 + cOffset] += tmp0 * a[i +  6 + aOffset0] + tmp1 * a[i +  6 + aOffset1] + tmp2 * a[i +  6 + aOffset2] + tmp3 * a[i +  6 + aOffset3];
                        c[i +  7 + cOffset] += tmp0 * a[i +  7 + aOffset0] + tmp1 * a[i +  7 + aOffset1] + tmp2 * a[i +  7 + aOffset2] + tmp3 * a[i +  7 + aOffset3];
                    }
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp0 * a[i + aOffset0] + tmp1 * a[i + aOffset1] + tmp2 * a[i + aOffset2] + tmp3 * a[i + aOffset3];
                    }
                    w += 4;
                }
                for (; p > 0; p--) {
                    int aOffset = (w - 1) * lda + _a_offset;
                    float tmp = b[(w - 1) + bColOffset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp * a[i      + aOffset];
                        c[i +  1 + cOffset] += tmp * a[i +  1 + aOffset];
                        c[i +  2 + cOffset] += tmp * a[i +  2 + aOffset];
                        c[i +  3 + cOffset] += tmp * a[i +  3 + aOffset];
                        c[i +  4 + cOffset] += tmp * a[i +  4 + aOffset];
                        c[i +  5 + cOffset] += tmp * a[i +  5 + aOffset];
                        c[i +  6 + cOffset] += tmp * a[i +  6 + aOffset];
                        c[i +  7 + cOffset] += tmp * a[i +  7 + aOffset];
                    }
                    // handle remainder
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp * a[i + aOffset];
                    }
                    w++;
                }
            } else {
                for (; p >= 4; p -= 4) {
                    int aOffset0 = (w - 1) * lda + _a_offset;
                    int aOffset1 = w * lda + _a_offset;
                    int aOffset2 = (w + 1) * lda + _a_offset;
                    int aOffset3 = (w + 2) * lda + _a_offset;
                    float tmp0 = alpha * b[(w - 1) + bColOffset];
                    float tmp1 = alpha * b[w + bColOffset];
                    float tmp2 = alpha * b[(w + 1) + bColOffset];
                    float tmp3 = alpha * b[(w + 2) + bColOffset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp0 * a[i      + aOffset0] + tmp1 * a[i      + aOffset1] + tmp2 * a[i      + aOffset2] + tmp3 * a[i      + aOffset3];
                        c[i +  1 + cOffset] += tmp0 * a[i +  1 + aOffset0] + tmp1 * a[i +  1 + aOffset1] + tmp2 * a[i +  1 + aOffset2] + tmp3 * a[i +  1 + aOffset3];
                        c[i +  2 + cOffset] += tmp0 * a[i +  2 + aOffset0] + tmp1 * a[i +  2 + aOffset1] + tmp2 * a[i +  2 + aOffset2] + tmp3 * a[i +  2 + aOffset3];
                        c[i +  3 + cOffset] += tmp0 * a[i +  3 + aOffset0] + tmp1 * a[i +  3 + aOffset1] + tmp2 * a[i +  3 + aOffset2] + tmp3 * a[i +  3 + aOffset3];
                        c[i +  4 + cOffset] += tmp0 * a[i +  4 + aOffset0] + tmp1 * a[i +  4 + aOffset1] + tmp2 * a[i +  4 + aOffset2] + tmp3 * a[i +  4 + aOffset3];
                        c[i +  5 + cOffset] += tmp0 * a[i +  5 + aOffset0] + tmp1 * a[i +  5 + aOffset1] + tmp2 * a[i +  5 + aOffset2] + tmp3 * a[i +  5 + aOffset3];
                        c[i +  6 + cOffset] += tmp0 * a[i +  6 + aOffset0] + tmp1 * a[i +  6 + aOffset1] + tmp2 * a[i +  6 + aOffset2] + tmp3 * a[i +  6 + aOffset3];
                        c[i +  7 + cOffset] += tmp0 * a[i +  7 + aOffset0] + tmp1 * a[i +  7 + aOffset1] + tmp2 * a[i +  7 + aOffset2] + tmp3 * a[i +  7 + aOffset3];
                    }
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp0 * a[i + aOffset0] + tmp1 * a[i + aOffset1] + tmp2 * a[i + aOffset2] + tmp3 * a[i + aOffset3];
                    }
                    w += 4;
                }
                for (; p > 0; p--) {
                    int aOffset = (w - 1) * lda + _a_offset;
                    float tmp = alpha * b[(w - 1) + bColOffset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp * a[i      + aOffset];
                        c[i +  1 + cOffset] += tmp * a[i +  1 + aOffset];
                        c[i +  2 + cOffset] += tmp * a[i +  2 + aOffset];
                        c[i +  3 + cOffset] += tmp * a[i +  3 + aOffset];
                        c[i +  4 + cOffset] += tmp * a[i +  4 + aOffset];
                        c[i +  5 + cOffset] += tmp * a[i +  5 + aOffset];
                        c[i +  6 + cOffset] += tmp * a[i +  6 + aOffset];
                        c[i +  7 + cOffset] += tmp * a[i +  7 + aOffset];
                    }
                    // handle remainder
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp * a[i + aOffset];
                    }
                    w++;
                }
            }
            u++;
        }
    }

    // Form C := alpha*A**T*B + beta*C
    private static void aTransTimesB(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        sgemmRange(false, true, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
    }

    private static void aTransTimesBRange(int m, int jFrom, int jTo, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        // Form C := alpha*A**T*B + beta*C
        final int kLimit = k - UNROLL + 1;
        for (int u = jFrom + 1; u <= jTo; u++) {
            int bOffset = (u - 1) * ldb + _b_offset;
            int cOffset = (u - 1) * ldc + _c_offset;
            int w = 0;
            // process 4 rows of A**T at a time, sharing each b[i] across all 4
            for (; w <= m - 4; w += 4) {
                int aOffset0 =  w      * lda + _a_offset;
                int aOffset1 = (w + 1) * lda + _a_offset;
                int aOffset2 = (w + 2) * lda + _a_offset;
                int aOffset3 = (w + 3) * lda + _a_offset;
                float s0a = 0.0f, s0b = 0.0f;
                float s1a = 0.0f, s1b = 0.0f;
                float s2a = 0.0f, s2b = 0.0f;
                float s3a = 0.0f, s3b = 0.0f;
                int i = 0;
                for (; i < kLimit; i += UNROLL) {
                    float b0  = b[i      + bOffset];
                    float b1  = b[i +  1 + bOffset];
                    float b2  = b[i +  2 + bOffset];
                    float b3  = b[i +  3 + bOffset];
                    float b4  = b[i +  4 + bOffset];
                    float b5  = b[i +  5 + bOffset];
                    float b6  = b[i +  6 + bOffset];
                    float b7  = b[i +  7 + bOffset];
                    s0a += a[i      + aOffset0]*b0  + a[i +  1 + aOffset0]*b1  + a[i +  2 + aOffset0]*b2 + a[i +  3 + aOffset0]*b3;
                    s0b += a[i +  4 + aOffset0]*b4  + a[i +  5 + aOffset0]*b5  + a[i +  6 + aOffset0]*b6 + a[i +  7 + aOffset0]*b7;
                    s1a += a[i      + aOffset1]*b0  + a[i +  1 + aOffset1]*b1  + a[i +  2 + aOffset1]*b2 + a[i +  3 + aOffset1]*b3;
                    s1b += a[i +  4 + aOffset1]*b4  + a[i +  5 + aOffset1]*b5  + a[i +  6 + aOffset1]*b6 + a[i +  7 + aOffset1]*b7;
                    s2a += a[i      + aOffset2]*b0  + a[i +  1 + aOffset2]*b1  + a[i +  2 + aOffset2]*b2 + a[i +  3 + aOffset2]*b3;
                    s2b += a[i +  4 + aOffset2]*b4  + a[i +  5 + aOffset2]*b5  + a[i +  6 + aOffset2]*b6 + a[i +  7 + aOffset2]*b7;
                    s3a += a[i      + aOffset3]*b0  + a[i +  1 + aOffset3]*b1  + a[i +  2 + aOffset3]*b2 + a[i +  3 + aOffset3]*b3;
                    s3b += a[i +  4 + aOffset3]*b4  + a[i +  5 + aOffset3]*b5  + a[i +  6 + aOffset3]*b6 + a[i +  7 + aOffset3]*b7;
                }
                // handle remainder
                for (; i < k; i++) {
                    float bv = b[i + bOffset];
                    s0a += a[i + aOffset0] * bv;
                    s1a += a[i + aOffset1] * bv;
                    s2a += a[i + aOffset2] * bv;
                    s3a += a[i + aOffset3] * bv;
                }
                float tmp0 = s0a + s0b;
                float tmp1 = s1a + s1b;
                float tmp2 = s2a + s2b;
                float tmp3 = s3a + s3b;
                if (beta == 0.0f) {
                    c[w     + cOffset] = alpha * tmp0;
                    c[w + 1 + cOffset] = alpha * tmp1;
                    c[w + 2 + cOffset] = alpha * tmp2;
                    c[w + 3 + cOffset] = alpha * tmp3;
                } else {
                    c[w     + cOffset] = alpha * tmp0 + beta * c[w     + cOffset];
                    c[w + 1 + cOffset] = alpha * tmp1 + beta * c[w + 1 + cOffset];
                    c[w + 2 + cOffset] = alpha * tmp2 + beta * c[w + 2 + cOffset];
                    c[w + 3 + cOffset] = alpha * tmp3 + beta * c[w + 3 + cOffset];
                }
            }
            // handle remainder rows
            for (; w < m; w++) {
                int aOffset = w * lda + _a_offset;
                float tmp = 0.0f;
                int i = 0;
                for (; i < kLimit; i += UNROLL) {
                    tmp += a[i      + aOffset] * b[i      + bOffset];
                    tmp += a[i +  1 + aOffset] * b[i +  1 + bOffset];
                    tmp += a[i +  2 + aOffset] * b[i +  2 + bOffset];
                    tmp += a[i +  3 + aOffset] * b[i +  3 + bOffset];
                    tmp += a[i +  4 + aOffset] * b[i +  4 + bOffset];
                    tmp += a[i +  5 + aOffset] * b[i +  5 + bOffset];
                    tmp += a[i +  6 + aOffset] * b[i +  6 + bOffset];
                    tmp += a[i +  7 + aOffset] * b[i +  7 + bOffset];
                }
                for (; i < k; i++) {
                    tmp += a[i + aOffset] * b[i + bOffset];
                }
                if (beta == 0.0f) {
                    c[w + cOffset] = alpha * tmp;
                } else {
                    c[w + cOffset] = alpha * tmp + beta * c[w + cOffset];
                }
            }
        }
    }

    // Form C := alpha*A*B**T + beta*C
    private static void aTimesBTrans(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        sgemmRange(true, false, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
    }

    private static void aTimesBTransRange(int m, int jFrom, int jTo, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        final int KC = 256;
        final int MC = 128;

        for (int j = jFrom; j < jTo; j += NC) {
            int jb = Math.min(jTo - j, NC);
            for (int p = 0; p < k; p += KC) {
                int pb = Math.min(k - p, KC);
                // On the first K-tile we apply `beta`, on subsequent K-tiles we must strictly accumulate
                float currentBeta = (p == 0) ? beta : 1.0f;
                for (int i = 0; i < m; i += MC) {
                    int ib = Math.min(m - i, MC);
                    int aTileOffset = _a_offset + i + p * lda;
                    int bTileOffset = _b_offset + j + p * ldb;
                    int cTileOffset = _c_offset + i + j * ldc;
                    kernel_aTimesBTrans(ib, jb, pb, alpha, a, aTileOffset, lda, b, bTileOffset, ldb, currentBeta, c, cTileOffset, ldc);
                }
            }
        }
    }

    // Inner block kernel
    private static void kernel_aTimesBTrans(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        final int mLimit = m - UNROLL + 1;
        int u = 1;
        for (int o = n; o > 0; o--) {
            int cOffset = (u - 1) * ldc + _c_offset;
            if (beta == 0.0f) {
                for (int i = 0; i < m; i++) {
                    c[i + cOffset] = 0.0f;
                }
            } else if (beta != 1.0f) {
                for (int i = 0; i < m; i++) {
                    c[i + cOffset] = beta * c[i + cOffset];
                }
            }
            int w = 1;
            int p = k;
            if (alpha == 1.0f) {
                for (; p >= 4; p -= 4) {
                    int aOffset0 = (w - 1) * lda + _a_offset;
                    int aOffset1 = w * lda + _a_offset;
                    int aOffset2 = (w + 1) * lda + _a_offset;
                    int aOffset3 = (w + 2) * lda + _a_offset;
                    float tmp0 = b[(u - 1) + (w - 1) * ldb + _b_offset];
                    float tmp1 = b[(u - 1) + w * ldb + _b_offset];
                    float tmp2 = b[(u - 1) + (w + 1) * ldb + _b_offset];
                    float tmp3 = b[(u - 1) + (w + 2) * ldb + _b_offset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp0 * a[i      + aOffset0] + tmp1 * a[i      + aOffset1] + tmp2 * a[i      + aOffset2] + tmp3 * a[i      + aOffset3];
                        c[i +  1 + cOffset] += tmp0 * a[i +  1 + aOffset0] + tmp1 * a[i +  1 + aOffset1] + tmp2 * a[i +  1 + aOffset2] + tmp3 * a[i +  1 + aOffset3];
                        c[i +  2 + cOffset] += tmp0 * a[i +  2 + aOffset0] + tmp1 * a[i +  2 + aOffset1] + tmp2 * a[i +  2 + aOffset2] + tmp3 * a[i +  2 + aOffset3];
                        c[i +  3 + cOffset] += tmp0 * a[i +  3 + aOffset0] + tmp1 * a[i +  3 + aOffset1] + tmp2 * a[i +  3 + aOffset2] + tmp3 * a[i +  3 + aOffset3];
                        c[i +  4 + cOffset] += tmp0 * a[i +  4 + aOffset0] + tmp1 * a[i +  4 + aOffset1] + tmp2 * a[i +  4 + aOffset2] + tmp3 * a[i +  4 + aOffset3];
                        c[i +  5 + cOffset] += tmp0 * a[i +  5 + aOffset0] + tmp1 * a[i +  5 + aOffset1] + tmp2 * a[i +  5 + aOffset2] + tmp3 * a[i +  5 + aOffset3];
                        c[i +  6 + cOffset] += tmp0 * a[i +  6 + aOffset0] + tmp1 * a[i +  6 + aOffset1] + tmp2 * a[i +  6 + aOffset2] + tmp3 * a[i +  6 + aOffset3];
                        c[i +  7 + cOffset] += tmp0 * a[i +  7 + aOffset0] + tmp1 * a[i +  7 + aOffset1] + tmp2 * a[i +  7 + aOffset2] + tmp3 * a[i +  7 + aOffset3];
                    }
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp0 * a[i + aOffset0] + tmp1 * a[i + aOffset1] + tmp2 * a[i + aOffset2] + tmp3 * a[i + aOffset3];
                    }
                    w += 4;
                }
                for (; p > 0; p--) {
                    int aOffset = (w - 1) * lda + _a_offset;
                    float tmp = b[(u - 1) + (w - 1) * ldb + _b_offset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp * a[i      + aOffset];
                        c[i +  1 + cOffset] += tmp * a[i +  1 + aOffset];
                        c[i +  2 + cOffset] += tmp * a[i +  2 + aOffset];
                        c[i +  3 + cOffset] += tmp * a[i +  3 + aOffset];
                        c[i +  4 + cOffset] += tmp * a[i +  4 + aOffset];
                        c[i +  5 + cOffset] += tmp * a[i +  5 + aOffset];
                        c[i +  6 + cOffset] += tmp * a[i +  6 + aOffset];
                        c[i +  7 + cOffset] += tmp * a[i +  7 + aOffset];
                    }
                    // handle remainder
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp * a[i + aOffset];
                    }
                    w++;
                }
            } else {
                for (; p >= 4; p -= 4) {
                    int aOffset0 = (w - 1) * lda + _a_offset;
                    int aOffset1 = w * lda + _a_offset;
                    int aOffset2 = (w + 1) * lda + _a_offset;
                    int aOffset3 = (w + 2) * lda + _a_offset;
                    float tmp0 = alpha * b[(u - 1) + (w - 1) * ldb + _b_offset];
                    float tmp1 = alpha * b[(u - 1) + w * ldb + _b_offset];
                    float tmp2 = alpha * b[(u - 1) + (w + 1) * ldb + _b_offset];
                    float tmp3 = alpha * b[(u - 1) + (w + 2) * ldb + _b_offset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp0 * a[i      + aOffset0] + tmp1 * a[i      + aOffset1] + tmp2 * a[i      + aOffset2] + tmp3 * a[i      + aOffset3];
                        c[i +  1 + cOffset] += tmp0 * a[i +  1 + aOffset0] + tmp1 * a[i +  1 + aOffset1] + tmp2 * a[i +  1 + aOffset2] + tmp3 * a[i +  1 + aOffset3];
                        c[i +  2 + cOffset] += tmp0 * a[i +  2 + aOffset0] + tmp1 * a[i +  2 + aOffset1] + tmp2 * a[i +  2 + aOffset2] + tmp3 * a[i +  2 + aOffset3];
                        c[i +  3 + cOffset] += tmp0 * a[i +  3 + aOffset0] + tmp1 * a[i +  3 + aOffset1] + tmp2 * a[i +  3 + aOffset2] + tmp3 * a[i +  3 + aOffset3];
                        c[i +  4 + cOffset] += tmp0 * a[i +  4 + aOffset0] + tmp1 * a[i +  4 + aOffset1] + tmp2 * a[i +  4 + aOffset2] + tmp3 * a[i +  4 + aOffset3];
                        c[i +  5 + cOffset] += tmp0 * a[i +  5 + aOffset0] + tmp1 * a[i +  5 + aOffset1] + tmp2 * a[i +  5 + aOffset2] + tmp3 * a[i +  5 + aOffset3];
                        c[i +  6 + cOffset] += tmp0 * a[i +  6 + aOffset0] + tmp1 * a[i +  6 + aOffset1] + tmp2 * a[i +  6 + aOffset2] + tmp3 * a[i +  6 + aOffset3];
                        c[i +  7 + cOffset] += tmp0 * a[i +  7 + aOffset0] + tmp1 * a[i +  7 + aOffset1] + tmp2 * a[i +  7 + aOffset2] + tmp3 * a[i +  7 + aOffset3];
                    }
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp0 * a[i + aOffset0] + tmp1 * a[i + aOffset1] + tmp2 * a[i + aOffset2] + tmp3 * a[i + aOffset3];
                    }
                    w += 4;
                }
                for (; p > 0; p--) {
                    int aOffset = (w - 1) * lda + _a_offset;
                    float tmp = alpha * b[(u - 1) + (w - 1) * ldb + _b_offset];
                    int i = 0;
                    for (; i < mLimit; i += UNROLL) {
                        c[i      + cOffset] += tmp * a[i      + aOffset];
                        c[i +  1 + cOffset] += tmp * a[i +  1 + aOffset];
                        c[i +  2 + cOffset] += tmp * a[i +  2 + aOffset];
                        c[i +  3 + cOffset] += tmp * a[i +  3 + aOffset];
                        c[i +  4 + cOffset] += tmp * a[i +  4 + aOffset];
                        c[i +  5 + cOffset] += tmp * a[i +  5 + aOffset];
                        c[i +  6 + cOffset] += tmp * a[i +  6 + aOffset];
                        c[i +  7 + cOffset] += tmp * a[i +  7 + aOffset];
                    }
                    // handle remainder
                    for (; i < m; i++) {
                        c[i + cOffset] += tmp * a[i + aOffset];
                    }
                    w++;
                }
            }
            u++;
        }
    }

    // Form C := alpha*A**T*B**T + beta*C
    private static void aTransTimesBTrans(int m, int n, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        sgemmRange(false, false, m, 0, n, k, alpha, a, _a_offset, lda, b, _b_offset, ldb, beta, c, _c_offset, ldc);
    }

    private static void aTransTimesBTransRange(int m, int jFrom, int jTo, int k, float alpha, float[] a, int _a_offset, int lda, float[] b, int _b_offset, int ldb, float beta, float[] c, int _c_offset, int ldc) {
        // Form C := alpha*A**T*B**T + beta*C
        final int kLimit = k - UNROLL + 1;
        for (int u = jFrom + 1; u <= jTo; u++) {
            int cOffset = (u - 1) * ldc + _c_offset;
            int bOffset = (u - 1) + _b_offset;
            int w = 0;
            // process 4 rows of A**T at a time, sharing each b[bOffset + i*ldb] across all 4
            for (; w <= m - 4; w += 4) {
                int aOffset0 =  w      * lda + _a_offset;
                int aOffset1 = (w + 1) * lda + _a_offset;
                int aOffset2 = (w + 2) * lda + _a_offset;
                int aOffset3 = (w + 3) * lda + _a_offset;
                float s0a = 0.0f, s0b = 0.0f;
                float s1a = 0.0f, s1b = 0.0f;
                float s2a = 0.0f, s2b = 0.0f;
                float s3a = 0.0f, s3b = 0.0f;
                int i = 0;
                for (; i < kLimit; i += UNROLL) {
                    float b0  = b[bOffset + (i     ) * ldb];
                    float b1  = b[bOffset + (i +  1) * ldb];
                    float b2  = b[bOffset + (i +  2) * ldb];
                    float b3  = b[bOffset + (i +  3) * ldb];
                    float b4  = b[bOffset + (i +  4) * ldb];
                    float b5  = b[bOffset + (i +  5) * ldb];
                    float b6  = b[bOffset + (i +  6) * ldb];
                    float b7  = b[bOffset + (i +  7) * ldb];
                    s0a += a[i      + aOffset0]*b0  + a[i +  1 + aOffset0]*b1  + a[i +  2 + aOffset0]*b2 + a[i +  3 + aOffset0]*b3;
                    s0b += a[i +  4 + aOffset0]*b4  + a[i +  5 + aOffset0]*b5  + a[i +  6 + aOffset0]*b6 + a[i +  7 + aOffset0]*b7;
                    s1a += a[i      + aOffset1]*b0  + a[i +  1 + aOffset1]*b1  + a[i +  2 + aOffset1]*b2 + a[i +  3 + aOffset1]*b3;
                    s1b += a[i +  4 + aOffset1]*b4  + a[i +  5 + aOffset1]*b5  + a[i +  6 + aOffset1]*b6 + a[i +  7 + aOffset1]*b7;
                    s2a += a[i      + aOffset2]*b0  + a[i +  1 + aOffset2]*b1  + a[i +  2 + aOffset2]*b2 + a[i +  3 + aOffset2]*b3;
                    s2b += a[i +  4 + aOffset2]*b4  + a[i +  5 + aOffset2]*b5  + a[i +  6 + aOffset2]*b6 + a[i +  7 + aOffset2]*b7;
                    s3a += a[i      + aOffset3]*b0  + a[i +  1 + aOffset3]*b1  + a[i +  2 + aOffset3]*b2 + a[i +  3 + aOffset3]*b3;
                    s3b += a[i +  4 + aOffset3]*b4  + a[i +  5 + aOffset3]*b5  + a[i +  6 + aOffset3]*b6 + a[i +  7 + aOffset3]*b7;
                }
                // handle remainder
                for (; i < k; i++) {
                    float bv = b[bOffset + i * ldb];
                    s0a += a[i + aOffset0] * bv;
                    s1a += a[i + aOffset1] * bv;
                    s2a += a[i + aOffset2] * bv;
                    s3a += a[i + aOffset3] * bv;
                }
                float tmp0 = s0a + s0b;
                float tmp1 = s1a + s1b;
                float tmp2 = s2a + s2b;
                float tmp3 = s3a + s3b;
                if (beta == 0.0f) {
                    c[w     + cOffset] = alpha * tmp0;
                    c[w + 1 + cOffset] = alpha * tmp1;
                    c[w + 2 + cOffset] = alpha * tmp2;
                    c[w + 3 + cOffset] = alpha * tmp3;
                } else {
                    c[w     + cOffset] = alpha * tmp0 + beta * c[w     + cOffset];
                    c[w + 1 + cOffset] = alpha * tmp1 + beta * c[w + 1 + cOffset];
                    c[w + 2 + cOffset] = alpha * tmp2 + beta * c[w + 2 + cOffset];
                    c[w + 3 + cOffset] = alpha * tmp3 + beta * c[w + 3 + cOffset];
                }
            }
            // handle remainder rows
            for (; w < m; w++) {
                int aOffset = w * lda + _a_offset;
                float tmp = 0.0f;
                int i = 0;
                for (; i < kLimit; i += UNROLL) {
                    tmp += a[i      + aOffset] * b[bOffset + (i     ) * ldb];
                    tmp += a[i +  1 + aOffset] * b[bOffset + (i +  1) * ldb];
                    tmp += a[i +  2 + aOffset] * b[bOffset + (i +  2) * ldb];
                    tmp += a[i +  3 + aOffset] * b[bOffset + (i +  3) * ldb];
                    tmp += a[i +  4 + aOffset] * b[bOffset + (i +  4) * ldb];
                    tmp += a[i +  5 + aOffset] * b[bOffset + (i +  5) * ldb];
                    tmp += a[i +  6 + aOffset] * b[bOffset + (i +  6) * ldb];
                    tmp += a[i +  7 + aOffset] * b[bOffset + (i +  7) * ldb];
                }
                for (; i < k; i++) {
                    tmp += a[i + aOffset] * b[bOffset + i * ldb];
                }
                if (beta == 0.0f) {
                    c[w + cOffset] = alpha * tmp;
                } else {
                    c[w + cOffset] = alpha * tmp + beta * c[w + cOffset];
                }
            }
        }
    }

    private static final class ColumnRangeTask implements Runnable {

        private final boolean notATransposed;
        private final boolean notBTransposed;
        private final int m;
        private final int jFrom;
        private final int jTo;
        private final int k;
        private final float alpha;
        private final float[] a;
        private final int aOffset;
        private final int lda;
        private final float[] b;
        private final int bOffset;
        private final int ldb;
        private final float beta;
        private final float[] c;
        private final int cOffset;
        private final int ldc;

        private ColumnRangeTask(boolean notATransposed, boolean notBTransposed, int m, int jFrom, int jTo, int k, float alpha, float[] a,
                int aOffset, int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc) {
            this.notATransposed = notATransposed;
            this.notBTransposed = notBTransposed;
            this.m = m;
            this.jFrom = jFrom;
            this.jTo = jTo;
            this.k = k;
            this.alpha = alpha;
            this.a = a;
            this.aOffset = aOffset;
            this.lda = lda;
            this.b = b;
            this.bOffset = bOffset;
            this.ldb = ldb;
            this.beta = beta;
            this.c = c;
            this.cOffset = cOffset;
            this.ldc = ldc;
        }

        @Override
        public void run() {
            sgemmRange(notATransposed, notBTransposed, m, jFrom, jTo, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset,
                    ldc);
        }
    }
}
