package math.nist;

/**
 * The log relative error, which is how NIST scores a package against its
 * reference datasets: the number of significant digits a computed value shares
 * with the certified one.
 * <p>
 * {@code -log10(|got - certified| / |certified|)}, so 7.5 means the two agree
 * to seven and a half digits. It is capped at 15, the number of digits NIST
 * publishes, and floored at 0, since a negative count says only "wrong" and the
 * distance between kinds of wrong is not interesting.
 * <p>
 * Where the certified value is exactly zero the relative form has no meaning --
 * Wampler1 and Wampler2 fit their data exactly, so every residual quantity is
 * zero -- and the absolute error is reported in its place.
 */
public final class Digits {

    /** The most any certified value here is published to. */
    public static final double MAX = 15.0;

    /**
     * Digits of agreement between a computed and a certified value.
     *
     * @param got
     *            what the library produced
     * @param certified
     *            what NIST certifies
     * @return the log relative error, between 0 and {@link #MAX}
     */
    public static double of(double got, double certified) {
        if (got == certified) {
            return MAX;
        }
        if (Double.isNaN(got) || Double.isInfinite(got)) {
            return 0.0;
        }
        if (certified == 0.0) {
            return clamp(-Math.log10(Math.abs(got)));
        }
        return clamp(-Math.log10(Math.abs(got - certified) / Math.abs(certified)));
    }

    /**
     * The worst agreement over a vector, which is the number that describes a
     * fit: a set of coefficients is only as good as its weakest one.
     *
     * @param got
     *            what the library produced
     * @param certified
     *            what NIST certifies, of the same length
     * @return the smallest log relative error over the entries
     */
    public static double worstOf(double[] got, double[] certified) {
        double worst = MAX;
        for (int i = 0; i < certified.length; ++i) {
            worst = Math.min(worst, of(got[i], certified[i]));
        }
        return worst;
    }

    private static double clamp(double digits) {
        return Math.max(0.0, Math.min(MAX, digits));
    }

    private Digits() {
        throw new AssertionError();
    }
}
