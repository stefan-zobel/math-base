package math.demo.michelson;

/**
 * The data behind {@link MichelsonDemo}, embedded as literals so that the demo
 * runs without a network, without a resource path, and against a table that
 * cannot shift under it.
 * <p>
 * <b>Michelson's 1879 determination of the speed of light.</b> Source:
 * <a href="https://www.itl.nist.gov/div898/software/dataplot/data/MICHELSO.DAT">NIST
 * Dataplot</a>, retrieved 2026-08-22, which prints Michelson's runs as
 * tabulated by Dorsey (1944), <i>Transactions of the American Philosophical
 * Society</i> 34, part 1, table 22. A publication of the United States
 * government, and measurements are facts rather than authorship.
 * <p>
 * 100 determinations made between 5 June and 2 July 1879, in five columns: the
 * value in millions of metres per second, the air temperature in degrees
 * Fahrenheit, the elapsed day counted from 1 on 5 June, whether the run was
 * made in the afternoon, and which of the 24 measurement sets it belongs to.
 * 37 runs were made in the morning and 63 in the afternoon; the temperature
 * ranges from 58 to 90 degrees.
 * <p>
 * NIST heads the file "speed of light in air", but the numbers are the values
 * reduced to vacuum, and the table itself says so: its mean is 299852.4 km/s,
 * which is 59.9 km/s <em>above</em> the vacuum value, whereas light in air at
 * sea level travels about 88 km/s <em>below</em> it. Comparing these values
 * against {@link #ACCEPTED} is therefore the comparison Michelson intended.
 * <p>
 * {@link #ACCEPTED} is exact rather than measured: since 1983 the metre has
 * been defined in terms of the speed of light, so the constant is fixed by
 * definition at 299792458 m/s. That is what makes this table worth a demo --
 * the truth is known without error, so every interval computed from these
 * numbers can be checked against it rather than against another estimate.
 * <p>
 * The sum of {@link #speed()} is 29985.24, which is the cheapest way to tell
 * whether a refresh of the file changed the table or only reformatted it.
 */
public final class Datasets {

    /**
     * The speed of light in vacuum in millions of metres per second, fixed by
     * the 1983 definition of the metre at exactly 299792458 m/s.
     */
    public static final double ACCEPTED = 299.792458;

    /** Number of measurement sets the 100 runs are grouped into. */
    public static final int SETS = 24;

    /** The determination of each run, in millions of metres per second. */
    private static final double[] SPEED = {
            299.85, 299.74, 299.90, 300.07, 299.93, 299.85, 299.95, 299.98, 299.98, 299.88,
            300.00, 299.98, 299.93, 299.65, 299.76, 299.81, 300.00, 300.00, 299.96, 299.96,
            299.96, 299.94, 299.96, 299.94, 299.88, 299.80, 299.85, 299.88, 299.90, 299.84,
            299.83, 299.79, 299.81, 299.88, 299.88, 299.83, 299.80, 299.79, 299.76, 299.80,
            299.88, 299.88, 299.88, 299.86, 299.72, 299.72, 299.62, 299.86, 299.97, 299.95,
            299.88, 299.91, 299.85, 299.87, 299.84, 299.84, 299.85, 299.84, 299.84, 299.84,
            299.89, 299.81, 299.81, 299.82, 299.80, 299.77, 299.76, 299.74, 299.75, 299.76,
            299.91, 299.92, 299.89, 299.86, 299.88, 299.72, 299.84, 299.85, 299.85, 299.78,
            299.89, 299.84, 299.78, 299.81, 299.76, 299.81, 299.79, 299.81, 299.82, 299.85,
            299.87, 299.87, 299.81, 299.74, 299.81, 299.94, 299.95, 299.80, 299.81, 299.87
    };

    /** Air temperature during each run, in degrees Fahrenheit. */
    private static final int[] TEMPERATURE = {
            76, 72, 72, 72, 72, 72, 83, 83, 83, 83, 83, 90, 90, 71, 71, 71, 72, 72, 72, 79,
            79, 79, 79, 79, 79, 79, 64, 64, 65, 66, 67, 84, 85, 84, 84, 84, 62, 63, 64, 77,
            77, 77, 77, 77, 58, 58, 59, 75, 75, 75, 60, 61, 62, 63, 78, 79, 80, 79, 79, 79,
            61, 62, 63, 64, 65, 80, 81, 82, 82, 81, 89, 89, 90, 90, 90, 72, 73, 74, 75, 76,
            86, 86, 73, 74, 75, 75, 76, 76, 85, 86, 86, 86, 83, 84, 86, 86, 86, 86, 86, 85
    };

    /** Day of each run, counted from 1 on 5 June 1879. */
    private static final int[] DAY = {
            1, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 6, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 13, 13, 13, 13,
            13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
            17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20,
            22, 22, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28
    };

    /** 1 for a run made in the afternoon, 0 for one made in the morning. */
    private static final int[] AFTERNOON = {
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1,
            1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
            1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };

    /** The measurement set each run belongs to. */
    private static final int[] SET = {
            1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7,
            7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 11,
            11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15,
            16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19,
            20, 20, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24
    };

    /**
     * Number of runs in the table.
     *
     * @return 100
     */
    public static int size() {
        return SPEED.length;
    }

    /**
     * Michelson's determinations, in millions of metres per second.
     *
     * @return a fresh copy of the 100 values, in the order they were made
     */
    public static double[] speed() {
        return SPEED.clone();
    }

    /**
     * Air temperature during each run, in degrees Fahrenheit.
     *
     * @return a fresh copy of the 100 temperatures
     */
    public static double[] temperature() {
        return toDouble(TEMPERATURE);
    }

    /**
     * Day of each run, counted from 1 on 5 June 1879.
     *
     * @return a fresh copy of the 100 day numbers
     */
    public static double[] day() {
        return toDouble(DAY);
    }

    /**
     * Whether each run was made in the afternoon, as an indicator that can go
     * straight into a design matrix: {@code 1.0} for afternoon, {@code 0.0}
     * for morning.
     *
     * @return a fresh copy of the 100 indicators
     */
    public static double[] afternoon() {
        return toDouble(AFTERNOON);
    }

    /**
     * The measurement set each run belongs to, {@code 1} to {@link #SETS}.
     *
     * @return a fresh copy of the 100 set numbers
     */
    public static int[] set() {
        return SET.clone();
    }

    /**
     * Sum of {@link #speed()}, the checksum quoted in the class comment.
     * Recomputed rather than stored, so that it verifies the literals above
     * instead of repeating them.
     *
     * @return the sum of the 100 determinations
     */
    public static double checksum() {
        double sum = 0.0;
        for (int i = 0; i < SPEED.length; ++i) {
            sum += SPEED[i];
        }
        return sum;
    }

    private static double[] toDouble(int[] values) {
        double[] copy = new double[values.length];
        for (int i = 0; i < values.length; ++i) {
            copy[i] = values[i];
        }
        return copy;
    }

    private Datasets() {
        throw new AssertionError();
    }
}
