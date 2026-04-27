package math.gemm;

/*package*/ class GemmSwitch {

    private static final boolean IS_VECTOR_API_PRESENT = DgemmMRxNR.isVectorApiPresent();
    private static final long TINY_WORK_THRESHOLD = 1_000_000L;

    // Single runtime dispatch rule:
    // true  -> use MRxNR API
    // false -> use Baseline API
    static boolean useMrxNrAtRuntime(
            boolean isFloat,
            boolean notATransposed,
            boolean notBTransposed,
            int m, int n, int k) {

        // Java 8+ path (non-vectorized MRxNR classes): Baseline wins consistently.
        // Also the Java 25 path when "--add-modules jdk.incubator.vector" is not used
        boolean vectorized = isFloat ? (IS_VECTOR_API_PRESENT && SgemmMRxNR.isVectorized())
                : (IS_VECTOR_API_PRESENT && DgemmMRxNR.isVectorized());
        if (!vectorized) {
            return false;
        }

        long work = (long) m * (long) n * (long) k;

        // Java 25 tiny/small region: Baseline is faster in all transpose cases.
        if (work <= TINY_WORK_THRESHOLD) {
            return false;
        }

        // Java 25 SGEMM: MRxNR wins from medium sizes onward (NN/TT/NT/TN).
        if (isFloat) {
            return true;
        }

        // Java 25 DGEMM: TT/NT/TN are MRxNR-favored once outside tiny region.
        boolean nn = notATransposed && notBTransposed;
        if (!nn) {
            return true;
        }

        // Conservative DGEMM-NN exception from measured crossover region.
        if (m <= 128 && n <= 256 && k <= 128) {
            return false;
        }

        return true;
    }
}
