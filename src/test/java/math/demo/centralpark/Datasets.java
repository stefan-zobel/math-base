package math.demo.centralpark;

/**
 * The data behind {@link WeatherDemo}, embedded as literals so that the demo
 * runs without a network, without a resource path, and against a table that
 * cannot shift under it.
 * <p>
 * <b>One year of daily weather at one station.</b> Source: NOAA NCEI
 * Global Historical Climatology Network - Daily, the per-station file
 * {@code https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/USW00094728.csv},
 * retrieved 2026-08-27. The dataset carries a DOI, {@code 10.7289/V5D21VHZ},
 * and the readme asks that it be cited along with Menne, M. J., I. Durre,
 * R. S. Vose, B. E. Gleason and T. G. Houston, 2012: <i>An overview of the
 * Global Historical Climatology Network-Daily Database</i>, Journal of
 * Atmospheric and Oceanic Technology 29, 897-910,
 * {@code doi:10.1175/JTECH-D-11-00103.1}.
 * <p>
 * The station is {@code USW00094728}, NY CITY CENTRAL PARK, at
 * {@code 40.77898} north and {@code 73.96925} west, {@code 42.7} m above sea
 * level, and it is a United States station on purpose. GHCN-Daily also carries
 * stations contributed by other countries under WMO Resolution 40, which
 * restricts redistribution; a US station has none of that and is in the public
 * domain outright.
 * <p>
 * The calendar year 2025, all 365 days, with no missing value in any of the
 * three columns kept here and no quality flag set anywhere in the year. That
 * is unusual enough to be worth saying: {@link WeatherDemo} has to create its
 * own gaps, because the record has none.
 * <p>
 * <b>The values are stored in the units NOAA publishes</b> -- tenths of a
 * degree Celsius for temperature, tenths of a millimetre for precipitation --
 * and as {@code int}, not as {@code double}. The resolution of these columns
 * is something the demo measures, so the literals have to be exactly what the
 * file says rather than something a conversion has already rounded.
 * <p>
 * <b>The tenth of a degree is fiction.</b> The 730 temperature values take
 * only 93 distinct values, and 73 of the 92 gaps between consecutive distinct
 * values are five or six tenths: that is the whole-Fahrenheit grid, since one
 * degree Fahrenheit is {@code 5/9} of a degree Celsius, or 5.56 tenths. 84 of
 * the 93 sit on that grid exactly, and from 27 August to the end of the year,
 * where the source is the CF6 form, all 254 values do. The nine that do not
 * are every one of them below freezing, out of the ASOS stretch, where the
 * conversion rounds a negative differently. Every one of the 121 positive
 * precipitation values is a whole hundredth of an inch. Anything that reads
 * these columns as though they resolved a tenth of a degree is working with a
 * grid five and a half times coarser than it thinks.
 * <p>
 * <b>The upstream file rewrites its own last twelve months.</b> The year has
 * two source flags: {@code W}, the WBAN/ASOS summary from NOAA's Integrated
 * Surface Data, up to 26 August, and then {@code 1}, the CF6 daily climate
 * summary of the US National Weather Service, which is published in whole
 * degrees Fahrenheit. The changeover sits exactly 365 days before the file was
 * last modified, and the four preceding years in the same download are all
 * {@code W} throughout -- their CF6 stretches have already been replaced. A
 * download taken three months from now returns different numbers for September
 * 2025. See {@link #SOURCE_BREAK}.
 * <p>
 * <b>A zero in the precipitation column is not always a dry day.</b> 35 of the
 * 244 zeros carry the measurement flag {@code T}, a trace: precipitation fell
 * but stayed below the smallest amount the gauge records. So there are 121
 * days with measurable rain, 35 with a trace, and 209 that were dry, and only
 * {@link #isTrace(int)} can tell the last two apart.
 * <p>
 * Checksums, all recomputed from the literals rather than stored: the maxima
 * sum to 6229.8 degrees, the minima to 3397.3, and the precipitation to
 * 1006.1 mm. Those three sums are the cheapest way to tell whether a refresh
 * of the file changed the table or only reformatted it.
 */
public final class Datasets {

    /** The GHCN-Daily identifier of the station these values come from. */
    public static final String STATION = "USW00094728";

    /** The name GHCN-Daily gives the station. */
    public static final String NAME = "NY CITY CENTRAL PARK, NY US";

    /** Latitude of the station in degrees north. */
    public static final double LATITUDE = 40.77898;

    /** Longitude of the station in degrees east, so negative here. */
    public static final double LONGITUDE = -73.96925;

    /** Elevation of the station in metres above sea level. */
    public static final double ELEVATION = 42.7;

    /** The calendar year these 365 days belong to. Not a leap year. */
    public static final int YEAR = 2025;

    /**
     * The first index whose source flag is the CF6 form rather than the ASOS
     * summary, {@code 0} based: 238, which is 27 August. Everything from here
     * to the end of the year is the stretch the upstream file will replace,
     * and it is the stretch that lies exactly on the Fahrenheit grid.
     */
    public static final int SOURCE_BREAK = 238;

    /**
     * One degree Fahrenheit in tenths of a degree Celsius, {@code 50/9}. The
     * spacing of the grid the temperature columns are really on.
     */
    public static final double FAHRENHEIT_STEP = 50.0 / 9.0;

    /** Daily maximum temperature in tenths of a degree Celsius. */
    private static final int[] TMAX = {
             106,   61,   39,    6,    6,  -10,    6,  -10,    6,   39,
              22,   56,   56,   17,    0,  -10,   56,   78,   50,  -32,
             -71,  -66,  -21,    6,    6,   61,   44,   56,  111,   33,
              89,   89,   -5,   94,   83,   22,   33,   61,   17,   22,
              28,   17,   22,   78,   22,   17,   83,   22,  -27,  -10,
             -16,   17,   44,   78,   94,  139,  133,  122,  106,  178,
               6,   39,  111,  122,  128,  100,   72,  117,  189,  183,
             128,  100,  139,   94,  161,  150,  161,  172,  111,  106,
             178,   89,  100,  144,  117,  117,  172,  272,   89,  228,
             144,   94,  228,  228,  150,  133,  100,   72,   94,  106,
              78,   44,  144,  189,  206,  100,  178,  200,  294,  228,
             139,  278,  244,  250,  250,  233,  150,  244,  261,  228,
             217,  289,  294,  206,  161,  178,  233,  256,  156,  222,
             250,  222,  206,  167,  206,  222,  278,  228,  206,  194,
             150,  106,  172,  178,  189,  228,  233,  183,  194,  239,
             194,  189,  217,  261,  283,  306,  289,  244,  244,  178,
             228,  272,  306,  256,  200,  178,  206,  194,  289,  311,
             278,  300,  311,  356,  372,  356,  294,  222,  283,  317,
             322,  317,  289,  311,  283,  289,  306,  306,  339,  322,
             289,  289,  283,  283,  294,  300,  306,  322,  278,  272,
             311,  294,  272,  267,  306,  350,  294,  289,  344,  361,
             350,  317,  228,  267,  289,  317,  300,  261,  272,  278,
             294,  322,  317,  328,  317,  317,  317,  289,  328,  239,
             250,  211,  222,  283,  278,  267,  300,  261,  239,  256,
             256,  244,  261,  256,  261,  261,  250,  278,  300,  206,
             228,  228,  228,  261,  256,  250,  272,  261,  233,  211,
             283,  289,  228,  217,  228,  283,  267,  250,  278,  261,
             283,  256,  272,  206,  189,  217,  278,  289,  267,  267,
             222,  150,  167,  178,  172,  144,  161,  200,  161,  172,
             183,  211,  200,  189,  167,  156,  150,  144,  150,  128,
             122,  133,  172,  144,  144,  150,  150,  150,  189,  139,
             161,  183,  172,  150,   50,  100,  106,  106,  100,  133,
              67,   72,   83,   83,  111,  111,   94,  111,  144,  150,
              94,   61,   56,   83,   61,   44,   50,   50,    0,   56,
              61,   33,   17,   83,   50,   17,   50,    6,  -17,    0,
              83,  100,  144,   33,   78,   39,   44,   78,   94,  -17,
              -6,   39,   94,    6,    6
    };

    /** Daily minimum temperature in tenths of a degree Celsius. */
    private static final int[] TMIN = {
              39,    6,   -5,  -21,  -21,  -55,  -71,  -49,  -55,  -27,
             -21,  -10,    6,  -38,  -49,  -49,  -32,   17,  -43,  -82,
            -116, -121,  -82,  -49,  -66,   -5,   -5,    0,   -5,  -38,
              22,  -60,  -82,   -5,  -10,  -32,  -10,  -10,  -32,  -10,
             -16,  -21,   -5,   17,  -16,  -21,    6,  -32,  -71,  -82,
             -77,  -60,  -60,    6,   17,   61,   61,   67,   50,  -27,
             -66,  -49,   11,   56,   22,   11,   22,    0,   72,   67,
              56,   39,   44,   50,   83,   28,   17,   61,   56,   33,
              44,   -5,   33,   83,   28,   11,   78,   61,   56,   72,
              50,   28,   67,  139,   61,   61,   39,   17,   -5,   44,
              44,   22,   39,   89,   89,   67,   56,   89,  122,  128,
             100,  117,  122,  133,  156,  117,   89,  106,  133,  150,
             111,  139,  172,  156,  139,  139,  128,  156,  106,  100,
             144,  139,  144,  133,  150,  161,  172,  161,  106,   94,
              94,   89,   83,   89,  111,  111,  139,  111,  128,  167,
             117,  100,  117,  133,  167,  211,  217,  189,  172,  167,
             167,  183,  211,  194,  150,  150,  156,  167,  183,  211,
             194,  222,  217,  267,  272,  272,  189,  167,  172,  250,
             217,  222,  217,  206,  194,  200,  222,  244,  228,  228,
             211,  228,  233,  217,  228,  222,  250,  239,  228,  206,
             239,  217,  194,  206,  211,  244,  233,  228,  244,  256,
             239,  183,  172,  172,  178,  194,  217,  178,  178,  172,
             178,  206,  222,  217,  217,  228,  233,  233,  228,  189,
             183,  144,  150,  144,  194,  194,  211,  172,  161,  161,
             183,  139,  150,  178,  167,  161,  178,  189,  183,  150,
             133,  144,  156,  161,  172,  178,  178,  189,  189,  172,
             178,  200,  161,  144,  139,  183,  211,  211,  200,  189,
             178,  194,  178,  133,  117,  122,  144,  178,  178,  189,
             133,   89,   78,  144,  128,  117,  117,  111,   78,   83,
              89,  128,  122,  100,  117,  100,   78,   72,   89,   61,
              61,   89,  100,  106,   89,   61,  100,   83,   89,   56,
              50,  122,  122,   17,    6,   33,   56,   44,   44,   50,
              28,   28,   39,   33,   67,   44,   22,   61,   78,   94,
              22,    6,    0,   28,   11,   11,   -6,  -44,  -67,  -11,
              11,  -56,  -72,    0,  -22,  -39,  -28,  -72,  -72,  -44,
             -11,   22,    0,  -11,  -17,  -22,   17,   17,  -17,  -67,
             -39,  -56,   -6,  -28,  -28
    };

    /**
     * Daily precipitation in tenths of a millimetre. A zero here is either a
     * dry day or a trace; {@link #isTrace(int)} is what separates them.
     */
    private static final int[] PRCP = {
               0,    0,    0,    0,    0,   20,    0,    0,    0,    0,
               5,    0,    0,    0,    0,    0,    0,   23,   66,    0,
               0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
              41,    0,   28,    0,    0,    0,   86,    0,  155,   15,
               0,   30,   28,   38,    0,   79,  196,    0,    0,    0,
               3,    0,    0,    0,    0,    0,    0,    3,    0,    0,
               0,    0,    0,  434,    0,    0,    0,    0,    0,    0,
               0,    0,    0,    0,    3,  234,    0,    0,  269,    8,
               0,    0,  122,    0,    0,    0,    8,    0,   10,  315,
              20,    0,   13,   23,   48,  193,   33,    0,    0,    3,
             104,  160,    0,   18,    3,    0,    0,    0,    0,    0,
               0,    0,    0,    0,    0,  208,    0,    0,    0,    0,
               0,    0,  241,   38,  188,   61,    5,    0,  175,    0,
               0,    0,   41,  310,    0,   28,    0,    0,    0,    0,
              23,  229,   66,    0,    0,    0,    0,  124,   25,    0,
             117,    0,    0,    0,    0,    0,    5,  178,    0,    3,
             102,    0,    0,    0,   58,    8,   38,    5,    8,  165,
               0,    0,   25,    0,    0,    0,    0,   18,   13,    0,
               0,  109,    8,    0,    0,    0,    0,    0,   23,   13,
              46,    0,    0,    0,  671,    0,    0,    0,    0,    0,
               0,    0,    0,    0,    0,    3,    0,    0,    0,    0,
              71,   81,    0,    0,    0,    0,    0,   15,    0,    0,
               0,    0,    0,    0,  180,    0,    0,    0,  168,    0,
               0,  130,   69,    0,    0,    0,    0,    0,    0,    0,
               0,    0,    0,    0,    0,    0,  102,   43,  257,  160,
               0,    0,   20,    0,    0,    0,    0,    0,    0,    0,
               3,    0,    0,    0,    0,    0,   20,   69,    0,    3,
              25,    0,    0,    0,    0,    0,    0,    0,    0,    0,
              69,    0,    0,   23,  168,  229,    0,    0,    0,    0,
               0,    0,    3,    8,   69,    0,    0,    0,    0,    0,
               0,    5,  465,    0,    0,    0,   10,    0,    0,    0,
               0,   41,   48,   56,    0,    0,    3,    0,   18,   15,
               0,    0,   43,    0,    0,   51,    0,    0,  188,   30,
               0,    0,    0,   28,    0,  198,    0,    0,    0,    5,
               0,    0,    0,   18,    0,    0,    8,  130,    0,    0,
               0,    5,  295,    0,    0,    0,   56,    0,    0,   69,
              56,    0,   20,    0,    0
    };

    /**
     * The 35 days whose precipitation carries the measurement flag {@code T},
     * {@code 0} based. Every one of them has a zero in {@link #PRCP}.
     */
    private static final int[] TRACE = {
               0,  15,  19,  31,  56,  80,  84,  87,  91, 110, 114, 127,
             129, 136, 143, 149, 158, 176, 180, 183, 187, 193, 196, 197,
             199, 200, 207, 212, 218, 259, 308, 310, 314, 361, 364
    };

    /** Days before the first of each month in a non-leap year. */
    private static final int[] MONTH_START = {
            0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
    };

    /**
     * How many days the year holds.
     *
     * @return 365
     */
    public static int size() {
        return TMAX.length;
    }

    /**
     * Daily maximum temperature in degrees Celsius.
     *
     * @return a fresh copy of the 365 maxima
     */
    public static double[] tmax() {
        return scaled(TMAX);
    }

    /**
     * Daily minimum temperature in degrees Celsius.
     *
     * @return a fresh copy of the 365 minima
     */
    public static double[] tmin() {
        return scaled(TMIN);
    }

    /**
     * Daily precipitation in millimetres.
     *
     * @return a fresh copy of the 365 totals
     */
    public static double[] precipitation() {
        return scaled(PRCP);
    }

    /**
     * Daily maximum temperature in the tenths of a degree the file publishes.
     * This is what anything asking about the resolution of the column has to
     * work on, because the conversion to degrees is where the grid stops being
     * visible.
     *
     * @return a fresh copy of the 365 maxima, in tenths
     */
    public static int[] tmaxTenths() {
        return TMAX.clone();
    }

    /**
     * Daily minimum temperature in the tenths of a degree the file publishes.
     *
     * @return a fresh copy of the 365 minima, in tenths
     */
    public static int[] tminTenths() {
        return TMIN.clone();
    }

    /**
     * Daily precipitation in the tenths of a millimetre the file publishes.
     *
     * @return a fresh copy of the 365 totals, in tenths
     */
    public static int[] precipitationTenths() {
        return PRCP.clone();
    }

    /**
     * The days flagged as a trace of precipitation.
     *
     * @return a fresh copy of the 35 indices, ascending
     */
    public static int[] trace() {
        return TRACE.clone();
    }

    /**
     * Whether the given day carries a trace flag, so that its zero in
     * {@link #precipitation()} means "too little to record" rather than "none".
     *
     * @param index
     *            the day, {@code 0} based
     * @return {@code true} if precipitation fell but was not measurable
     * @throws IndexOutOfBoundsException
     *             if {@code index} is outside the year
     */
    public static boolean isTrace(int index) {
        checkIndex(index);
        for (int i = 0; i < TRACE.length; ++i) {
            if (TRACE[i] == index) {
                return true;
            }
        }
        return false;
    }

    /**
     * Whether the given day had measurable precipitation, which is neither a
     * dry day nor a trace.
     *
     * @param index
     *            the day, {@code 0} based
     * @return {@code true} if the gauge recorded more than zero
     * @throws IndexOutOfBoundsException
     *             if {@code index} is outside the year
     */
    public static boolean isWet(int index) {
        checkIndex(index);
        return PRCP[index] > 0;
    }

    /**
     * Calendar month of a day, {@code 1} based.
     *
     * @param index
     *            the day, {@code 0} based
     * @return the month, 1 to 12
     * @throws IndexOutOfBoundsException
     *             if {@code index} is outside the year
     */
    public static int month(int index) {
        checkIndex(index);
        int m = 0;
        while (index >= MONTH_START[m + 1]) {
            ++m;
        }
        return m + 1;
    }

    /**
     * Day of the month, {@code 1} based.
     *
     * @param index
     *            the day, {@code 0} based
     * @return the day of the month
     * @throws IndexOutOfBoundsException
     *             if {@code index} is outside the year
     */
    public static int dayOfMonth(int index) {
        return index - MONTH_START[month(index) - 1] + 1;
    }

    /**
     * The date of a day as {@code YYYY-MM-DD}.
     *
     * @param index
     *            the day, {@code 0} based
     * @return the date, printable
     * @throws IndexOutOfBoundsException
     *             if {@code index} is outside the year
     */
    public static String label(int index) {
        int m = month(index);
        int d = dayOfMonth(index);
        return YEAR + "-" + (m < 10 ? "0" : "") + m + "-" + (d < 10 ? "0" : "") + d;
    }

    /**
     * Sum of {@link #tmax()}, one of the three checksums quoted in the class
     * comment. Recomputed rather than stored, so that it verifies the literals
     * instead of repeating them.
     *
     * @return the sum of the 365 maxima, in degrees Celsius
     */
    public static double checksumTmax() {
        return sum(TMAX);
    }

    /**
     * Sum of {@link #tmin()}, the second checksum.
     *
     * @return the sum of the 365 minima, in degrees Celsius
     */
    public static double checksumTmin() {
        return sum(TMIN);
    }

    /**
     * Sum of {@link #precipitation()}, the third checksum.
     *
     * @return the sum of the 365 totals, in millimetres
     */
    public static double checksumPrecipitation() {
        return sum(PRCP);
    }

    private static double sum(int[] tenths) {
        long total = 0L;
        for (int i = 0; i < tenths.length; ++i) {
            total += tenths[i];
        }
        return total / 10.0;
    }

    private static double[] scaled(int[] tenths) {
        double[] copy = new double[tenths.length];
        for (int i = 0; i < tenths.length; ++i) {
            copy[i] = tenths[i] / 10.0;
        }
        return copy;
    }

    private static void checkIndex(int index) {
        if (index < 0 || index >= TMAX.length) {
            throw new IndexOutOfBoundsException("not a day of " + YEAR + ": " + index);
        }
    }

    private Datasets() {
        throw new AssertionError();
    }
}
