package math.demo.hubble;

/**
 * The data behind {@link HubbleDemo}, embedded as literals so that the demo
 * runs without a network, without a resource path, and against a table that
 * cannot shift under it.
 * <p>
 * <b>Table 1 of Hubble (1929).</b> "A relation between distance and radial
 * velocity among extra-galactic nebulae", <i>Proceedings of the National
 * Academy of Sciences</i> 15(3):168-173, communicated 17 January 1929. A United
 * States publication from before 1930 and therefore in the public domain. The
 * 24 rows here were read out of the paper itself and agree with an independent
 * transcription hosted by NASA in every value.
 * <p>
 * The columns are the paper's own, with its footnotes:
 * <ul>
 * <li>{@link #distance()}, printed as {@code r}, in units of 10^6 parsecs. The
 * first two are Shapley's values.</li>
 * <li>{@link #velocity()}, printed as {@code v}, measured radial velocities in
 * km/s. Those of N.G.C. 6822, 221, 224 and 5457 were recent determinations by
 * Humason.</li>
 * <li>{@link #brightestStarMagnitude()}, printed as {@code ms}, the photographic
 * magnitude of the brightest stars involved, from which most of the distances
 * were derived. It is {@code NaN} for the ten objects where the paper prints
 * two dots: the six nearest, whose distances rest on individual stars, and the
 * four in the Virgo cluster, whose distances rest on the cluster mean.</li>
 * <li>{@link #visualMagnitude()}, printed as {@code mt}, Holetschek's visual
 * magnitude as corrected by Hopmann; for the first three objects an estimate by
 * Hubble.</li>
 * </ul>
 * <p>
 * The paper's fifth column, the absolute magnitude {@code Mt}, is left out
 * because it is not independent: it follows from {@code mt} and {@code r}.
 * <p>
 * <b>What the table is not.</b> The slope Hubble himself reports is
 * {@link #HUBBLE_K} and its companion {@code 513 +/- 60}, both from a solution
 * that removes the motion of the sun at the same time. Table 1 carries no
 * coordinates, so that solution cannot be reproduced from these numbers, and
 * {@link HubbleDemo} does not pretend to: it fits the same table the ways a
 * reader would fit it today. {@link #MODERN} is what the constant is now
 * believed to be, and the distance between it and everything on this page is
 * the demo's subject.
 * <p>
 * Checksums: the distances sum to 21.873, the velocities to 8955.
 */
public final class Datasets {

    /**
     * The slope Hubble reports in the paper, in km/s per 10^6 parsecs, from the
     * solution that treats the nebulae individually. The companion solution,
     * which groups them, gives {@code 513 +/- 60}.
     */
    public static final double HUBBLE_K = 465.0;

    /** The uncertainty Hubble quotes with {@link #HUBBLE_K}. */
    public static final double HUBBLE_K_ERROR = 50.0;

    /**
     * The Hubble constant as it is understood today, in km/s per megaparsec.
     * Measurements disagree at the level of a few units; the figure here is
     * only ever used to show that everything derived from Table 1 is an order
     * of magnitude away from it.
     */
    public static final double MODERN = 70.0;

    /** Object names as the paper prints them. */
    private static final String[] NAME = { "S. Mag.", "L. Mag.", "N.G.C. 6822", "N.G.C. 598", "N.G.C. 221",
            "N.G.C. 224", "N.G.C. 5457", "N.G.C. 4736", "N.G.C. 5194", "N.G.C. 4449", "N.G.C. 4214", "N.G.C. 3031",
            "N.G.C. 3627", "N.G.C. 4826", "N.G.C. 5236", "N.G.C. 1068", "N.G.C. 5055", "N.G.C. 7331", "N.G.C. 4258",
            "N.G.C. 4151", "N.G.C. 4382", "N.G.C. 4472", "N.G.C. 4486", "N.G.C. 4649" };

    /** Distance r, in units of 10^6 parsecs. */
    private static final double[] DISTANCE = { 0.032, 0.034, 0.214, 0.263, 0.275, 0.275, 0.45, 0.5, 0.5, 0.63, 0.8,
            0.9, 0.9, 0.9, 0.9, 1.0, 1.1, 1.1, 1.4, 1.7, 2.0, 2.0, 2.0, 2.0 };

    /** Radial velocity v, in km/s. */
    private static final double[] VELOCITY = { 170.0, 290.0, -130.0, -70.0, -185.0, -220.0, 200.0, 290.0, 270.0,
            200.0, 300.0, -30.0, 650.0, 150.0, 500.0, 920.0, 450.0, 500.0, 500.0, 960.0, 500.0, 850.0, 800.0, 1090.0 };

    /** Photographic magnitude ms of the brightest stars, NaN where the paper prints two dots. */
    private static final double[] BRIGHTEST_STAR = { Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            Double.NaN, 17.0, 17.3, 17.3, 17.8, 18.3, 18.5, 18.5, 18.5, 18.5, 18.7, 19.0, 19.0, 19.5, 20.0,
            Double.NaN, Double.NaN, Double.NaN, Double.NaN };

    /** Visual magnitude mt. */
    private static final double[] VISUAL = { 1.5, 0.5, 9.0, 7.0, 8.8, 5.0, 9.9, 8.4, 7.4, 9.5, 11.3, 8.3, 9.1, 9.0,
            10.4, 9.1, 9.6, 10.4, 8.7, 12.0, 10.0, 8.8, 9.7, 9.5 };

    /**
     * Number of objects in the table.
     *
     * @return 24
     */
    public static int size() {
        return DISTANCE.length;
    }

    /**
     * Distances in units of 10^6 parsecs, in the order the paper prints them.
     *
     * @return a fresh copy of the 24 distances
     */
    public static double[] distance() {
        return DISTANCE.clone();
    }

    /**
     * Radial velocities in km/s, in the order the paper prints them.
     *
     * @return a fresh copy of the 24 velocities
     */
    public static double[] velocity() {
        return VELOCITY.clone();
    }

    /**
     * Photographic magnitude of the brightest stars involved, {@code NaN} for
     * the ten objects the paper leaves blank.
     *
     * @return a fresh copy of the 24 magnitudes
     */
    public static double[] brightestStarMagnitude() {
        return BRIGHTEST_STAR.clone();
    }

    /**
     * Visual magnitudes.
     *
     * @return a fresh copy of the 24 magnitudes
     */
    public static double[] visualMagnitude() {
        return VISUAL.clone();
    }

    /**
     * Name of one object as the paper prints it.
     *
     * @param index
     *            row index, {@code 0} to {@code 23}
     * @return the name
     */
    public static String name(int index) {
        return NAME[index];
    }

    /**
     * Sum of {@link #distance()}, the checksum quoted in the class comment.
     *
     * @return 21.873
     */
    public static double distanceChecksum() {
        return sum(DISTANCE);
    }

    /**
     * Sum of {@link #velocity()}, the checksum quoted in the class comment.
     *
     * @return 8955
     */
    public static double velocityChecksum() {
        return sum(VELOCITY);
    }

    private static double sum(double[] values) {
        double total = 0.0;
        for (int i = 0; i < values.length; ++i) {
            total += values[i];
        }
        return total;
    }

    private Datasets() {
        throw new AssertionError();
    }
}
