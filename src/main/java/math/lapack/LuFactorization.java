package math.lapack;

import java.util.Objects;

/**
 * An LU factorization with partial pivoting that survives its right hand side:
 * {@link #factor(int, double[], int[]) factor} once, then
 * {@link #solve(int, int, double[], int[], double[]) solve} as often as there
 * are vectors to solve for.
 * <p>
 * {@link Dgesv} fuses the two, so a caller whose right hand sides arrive one
 * after another pays for the factorization again each time. Where the matrix
 * outlives the vector -- a Rosenbrock step solving six systems with one matrix,
 * a Newton iteration that holds its Jacobian, an implicit method of any kind --
 * that is the cubic part of the work repeated for nothing, and this class is
 * the same two LAPACK routines {@code Dgesv} calls, separated.
 * <p>
 * The matrix is <em>flat and column-major</em>, as everywhere else in this
 * library: entry {@code (i, j)} sits at {@code a[j * lda + i]}.
 * {@link #factor(int, double[], int[]) factor} overwrites it with {@code L} and
 * {@code U} of {@code P A = L U}, the unit diagonal of {@code L} not stored,
 * and writes the row interchanges into {@code ipiv}. Both are needed by
 * {@link #solve(int, int, double[], int[], double[]) solve} and neither may be
 * disturbed in between.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/LU_decomposition">Wikipedia LU
 * decomposition</a>.
 *
 * @since 1.5.3
 */
public final class LuFactorization {

    /**
     * Factors an {@code n} by {@code n} matrix stored contiguously, the common
     * case of {@link #factor(int, double[], int, int, int[], int)} with no
     * offsets and a leading dimension of {@code n}.
     *
     * @param n
     *            the order of the matrix, zero or greater
     * @param a
     *            a <code>double[]</code> of length at least {@code n * n}
     *            holding the matrix in column-major order, overwritten with its
     *            factors (modified)
     * @param ipiv
     *            an <code>int[]</code> of length at least {@code n} receiving
     *            the pivot indices (modified)
     * @return {@code true} if the factorization succeeded, {@code false} if
     *         {@code U} came out exactly singular, in which case nothing can be
     *         solved with it
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, {@code n} is negative or an
     *             array is too short
     */
    public static boolean factor(int n, double[] a, int[] ipiv) {
        return factor(n, a, 0, Math.max(1, n), ipiv, 0);
    }

    /**
     * Factors an {@code n} by {@code n} matrix as {@code P A = L U} by
     * Gaussian elimination with partial pivoting, in place.
     *
     * @param n
     *            the order of the matrix, zero or greater
     * @param a
     *            a <code>double[]</code> holding the matrix in column-major
     *            order, overwritten with its factors (modified)
     * @param aOffset
     *            offset of the first entry of the matrix in {@code a}
     * @param lda
     *            the leading dimension of the matrix, at least
     *            {@code max(1, n)}
     * @param ipiv
     *            an <code>int[]</code> receiving the pivot indices, one per
     *            row, in the convention {@link Dgesv} uses (modified)
     * @param ipivOffset
     *            offset of the first pivot index in {@code ipiv}
     * @return {@code true} if the factorization succeeded, {@code false} if
     *         {@code U} came out exactly singular, in which case nothing can be
     *         solved with it
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, {@code n}, an offset or
     *             {@code lda} is out of range, or an array is too short
     */
    public static boolean factor(int n, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset) {
        checkMatrix(n, a, aOffset, lda, "a", "lda");
        checkPivots(n, ipiv, ipivOffset);
        if (n == 0) {
            return true;
        }
        intW info = new intW(0);
        Dgetrf.dgetrf(n, n, a, aOffset, lda, ipiv, ipivOffset, info);
        if (info.val < 0) {
            throw new IllegalArgumentException("illegal argument at position " + (-info.val));
        }
        return info.val == 0;
    }

    /**
     * Solves {@code A X = B} for {@code X} with a factorization already
     * computed, the common case of
     * {@link #solve(int, int, double[], int, int, int[], int, double[], int, int)}
     * with no offsets and leading dimensions of {@code n}.
     *
     * @param n
     *            the order of the matrix, zero or greater
     * @param nrhs
     *            the number of right hand sides, that is columns of {@code b},
     *            zero or greater
     * @param lu
     *            a <code>double[]</code> holding the factors left by
     *            {@link #factor(int, double[], int[])} (not modified)
     * @param ipiv
     *            the <code>int[]</code> pivot indices left by the same call
     *            (not modified)
     * @param b
     *            a <code>double[]</code> of length at least {@code n * nrhs}
     *            holding the right hand sides in column-major order,
     *            overwritten with the solutions (modified)
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, a size is negative or an
     *             array is too short
     */
    public static void solve(int n, int nrhs, double[] lu, int[] ipiv, double[] b) {
        solve(n, nrhs, lu, 0, Math.max(1, n), ipiv, 0, b, 0, Math.max(1, n));
    }

    /**
     * Solves {@code A X = B} for {@code X} with a factorization already
     * computed, by forward and back substitution and the row interchanges that
     * produced it. The factors are read and not written, so the same
     * factorization serves any number of calls.
     *
     * @param n
     *            the order of the matrix, zero or greater
     * @param nrhs
     *            the number of right hand sides, that is columns of {@code b},
     *            zero or greater
     * @param lu
     *            a <code>double[]</code> holding the factors left by
     *            {@link #factor(int, double[], int, int, int[], int)} (not
     *            modified)
     * @param luOffset
     *            offset of the first entry of the factors in {@code lu}
     * @param ldlu
     *            the leading dimension of the factors, at least
     *            {@code max(1, n)}
     * @param ipiv
     *            the <code>int[]</code> pivot indices left by the same call
     *            (not modified)
     * @param ipivOffset
     *            offset of the first pivot index in {@code ipiv}
     * @param b
     *            a <code>double[]</code> holding the right hand sides in
     *            column-major order, overwritten with the solutions (modified)
     * @param bOffset
     *            offset of the first entry of the right hand sides in {@code b}
     * @param ldb
     *            the leading dimension of the right hand sides, at least
     *            {@code max(1, n)}
     * @throws IllegalArgumentException
     *             if an argument is {@code null}, a size, an offset or a
     *             leading dimension is out of range, or an array is too short
     */
    public static void solve(int n, int nrhs, double[] lu, int luOffset, int ldlu, int[] ipiv, int ipivOffset,
            double[] b, int bOffset, int ldb) {
        checkMatrix(n, lu, luOffset, ldlu, "lu", "ldlu");
        checkPivots(n, ipiv, ipivOffset);
        if (nrhs < 0) {
            throw new IllegalArgumentException("nrhs must not be negative, but is " + nrhs);
        }
        Objects.requireNonNull(b, "b");
        if (ldb < Math.max(1, n)) {
            throw new IllegalArgumentException("ldb must be at least " + Math.max(1, n) + ", but is " + ldb);
        }
        if (bOffset < 0) {
            throw new IllegalArgumentException("bOffset must not be negative, but is " + bOffset);
        }
        if (nrhs > 0 && b.length - bOffset < (long) ldb * (nrhs - 1) + n) {
            throw new IllegalArgumentException("b is too short for " + nrhs + " right hand sides of length " + n);
        }
        if (n == 0 || nrhs == 0) {
            return;
        }
        intW info = new intW(0);
        Dgetrs.dgetrs("No transpose", n, nrhs, lu, luOffset, ldlu, ipiv, ipivOffset, b, bOffset, ldb, info);
        if (info.val < 0) {
            throw new IllegalArgumentException("illegal argument at position " + (-info.val));
        }
    }

    private static void checkMatrix(int n, double[] a, int offset, int lda, String name, String ldName) {
        Objects.requireNonNull(a, name);
        if (n < 0) {
            throw new IllegalArgumentException("n must not be negative, but is " + n);
        }
        if (lda < Math.max(1, n)) {
            throw new IllegalArgumentException(ldName + " must be at least " + Math.max(1, n) + ", but is " + lda);
        }
        if (offset < 0) {
            throw new IllegalArgumentException(name + "Offset must not be negative, but is " + offset);
        }
        if (n > 0 && a.length - offset < (long) lda * (n - 1) + n) {
            throw new IllegalArgumentException(name + " is too short for a matrix of order " + n);
        }
    }

    private static void checkPivots(int n, int[] ipiv, int offset) {
        Objects.requireNonNull(ipiv, "ipiv");
        if (offset < 0) {
            throw new IllegalArgumentException("ipivOffset must not be negative, but is " + offset);
        }
        if (ipiv.length - offset < n) {
            throw new IllegalArgumentException("ipiv must hold at least " + n + " indices");
        }
    }

    private LuFactorization() {
        throw new AssertionError();
    }
}
