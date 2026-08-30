package math.stats;

/**
 * Which departure from the null hypothesis a test is looking for.
 * <p>
 * The choice belongs to the question, not to the data: it has to be made
 * before the sample is seen, and it changes what the p-value means. A
 * one-sided p-value of 0.03 and a two-sided one of 0.03 are different
 * statements, which is why the answer travels back inside
 * {@link TestResult#alternative} rather than only into the call.
 *
 * @since 1.5.3
 */
public enum Alternative {

    /** The parameter differs from the null value, in either direction. */
    TWO_SIDED,

    /** The parameter is smaller than the null value. */
    LESS,

    /** The parameter is larger than the null value. */
    GREATER
}
