package math.demo.measles;

import java.util.Locale;

import math.fun.DFunction;
import math.fun.DVectorField;
import math.ode.ButcherTableau;
import math.ode.Event;
import math.ode.ExplicitRungeKutta;
import math.ode.OdeEvent;
import math.ode.OdeIntegrator;
import math.ode.StepController;
import math.probe.ACF;
import math.solve.RootFinder;

/**
 * Thirty-five years of weekly measles reports before the vaccine, and one
 * question a solver cannot answer for you: <b>is the number it just handed you a
 * property of your equation, or of the step sizes it happened to choose?</b>
 * <p>
 * The demo asks it at a place the record picks rather than at one it invents.
 * The standard model of measles has exactly one free parameter -- how far the
 * contact rate swings over a year, because school terms do -- and the record
 * measures it: a heavy measles year against the light year beside it is a ratio
 * the model reproduces monotonically, so inverting it gives the swing. Doing
 * that year by year shows the record is not in one state. Massachusetts in the
 * 1930s calibrates to about {@code 0.06} and in the early 1950s to about
 * {@code 0.17}, and {@code 0.17} turns out to be a short step from where the
 * model stops having an answer at all.
 * <p>
 * <b>Section 3 says what is being read and shows that it exists</b>, because
 * every number after it is taken off a long trajectory and a reader is owed
 * that. Run forward, the model forgets where it started: two starts a factor of
 * three apart in the infectious count agree to twelve digits after two hundred
 * years, losing their difference at about a factor of a million per century.
 * What is left is a property of the equation and not of anybody's guess at 1928,
 * and it is where the demo starts looking. A <i>peak</i> is then the moment
 * {@code I' = 0} on the way down, found as an event on the interpolant rather
 * than sampled off a grid; the <i>trough</i> is the smallest the infectious
 * count gets; and the <i>pattern</i> is what the sequence of peak sizes does --
 * sizes and not spacings, since a two year cycle stays locked to the school year
 * and its mean interval reads {@code 1.0000} exactly as an annual one does.
 * <p>
 * Past that point the cycle breaks up. A trough is then a property of the
 * particular trajectory and not of the attractor, so tightening the tolerance a
 * decade does not converge on an answer -- there is none to converge on -- and
 * the runs carry the model to under one infective, which a population cannot do.
 * The record says the same from the other side, in the shape Bartlett named:
 * across four states the epidemic cycle weakens monotonically as the population
 * falls, because below some size measles does not cycle, it arrives, burns and
 * goes out.
 * <p>
 * <b>Where the edge is cannot be read off the parameter</b>, which is the
 * argument for running the check rather than for knowing when to. Between about
 * {@code 0.19} and {@code 0.25} orderly windows sit inside the disorder, so a
 * swing that answers can lie between two that do not.
 * <p>
 * Section 6 is the number none of that needed and a public health department
 * needs first, and it is where the eight days of latency finally matter for
 * their own sake: the recurrence period of the model's spiral into equilibrium
 * gives {@code R0}, and {@code 1 - 1/R0} is the fraction of a population that
 * has to be immune. One epidemic wave cannot give it -- measured on Kermack and
 * McKendrick's 1905-06 Bombay plague counts, where an excellent fit leaves
 * {@code R0} anywhere between {@code 1.07} and {@code 3.20} -- and thirty-five
 * years of recurrence can.
 * <p>
 * Two runs of {@link #main(String[])} produce identical output, on either source
 * tree and on either runtime. That is a constraint on <b>what</b> may be printed
 * and not only on how. Transcendentals in the models go through
 * {@code StrictMath}, but no such care reaches inside {@link StepController},
 * which picks a step size with {@code Math.pow} -- allowed one unit in the last
 * place, and taking a different intrinsic on JDK 8 and JDK 25. On a cycle that
 * is harmless. Past the edge it is not, so <b>a quantity whose two tolerances
 * disagree is not printed at all</b>: a cell that carries a number is a cell
 * whose number was checked.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)">
 * Wikipedia compartmental models</a>.
 */
public final class MeaslesDemo {

    /** Formatting locale, so that the output does not depend on the machine. */
    private static final Locale L = Locale.ROOT;

    /** Rate out of the latent class, per year: eight days. */
    static final double SIGMA = 365.0 / 8.0;

    /** Rate out of the infectious class, per year: five days. */
    static final double GAMMA = 365.0 / 5.0;

    /** Birth and death rate per year, a fifty year expectation of life. */
    static final double MU = 1.0 / 50.0;

    /** Rate out of the single infected class of the SIR, per year: thirteen days. */
    static final double GAMMA_SIR = 365.0 / 13.0;

    /** Population of Massachusetts in the middle of the record, to within a percent. */
    static final double MASSACHUSETTS_PEOPLE = 4.7e6;

    /** Years discarded before an attractor is read, so that the transient is gone. */
    static final double SETTLE_YEARS = 200.0;

    /** Years of attractor examined after that. */
    static final double WATCH_YEARS = 40.0;

    /** Relative tolerance of every integration here. */
    static final double RTOL = 1.0e-12;

    /** Absolute tolerance, well under one person in a population of millions. */
    static final double ATOL = 1.0e-18;

    /** The forcing amplitudes stepped through in sections 3 and 4. */
    static final double[] AMPLITUDES = { 0.02, 0.05, 0.08, 0.15, 0.20, 0.25 };

    /** Basic reproduction number used wherever a model has to be posed rather than fitted. */
    static final double R0_POSED = 17.0;

    private MeaslesDemo() {
        throw new AssertionError();
    }

    /**
     * Runs the whole page.
     *
     * @param args
     *            ignored
     */
    public static void main(String[] args) {
        rule("MEASLES: A RECORD THAT WALKED TO THE EDGE OF ITS OWN MODEL");
        theRecord();
        whatTheRecordSaysAboutItself();
        whatTheModelSettlesInto();
        theOneKnob();
        theEdgeBesideIt();
        theSmallerStates();
        theNumberAndTheThreshold();
    }

    // ---------------------------------------------------------------- 1

    /** Annual totals, and the vaccine, with no model in sight. */
    static void theRecord() {
        rule("1. The record");
        double[] annual = annualTotals(Datasets.massachusetts());
        System.out.println("Massachusetts, reported measles cases per 100,000 per year.");
        System.out.println("The vaccine was licensed in " + Datasets.VACCINE_YEAR + ".");
        System.out.println();
        for (int i = 0; i < annual.length; ++i) {
            int year = Datasets.FIRST_YEAR + i;
            System.out.printf(L, "  %d %7.0f %s%s%n", year, annual[i], bar(annual[i] / 25.0),
                    year == Datasets.VACCINE_YEAR ? "   <- vaccine" : "");
        }
        double before = mean(annual, 1958 - Datasets.FIRST_YEAR, 1962 - Datasets.FIRST_YEAR);
        double after = mean(annual, 1966 - Datasets.FIRST_YEAR, 1969 - Datasets.FIRST_YEAR);
        System.out.println();
        System.out.printf(L, "  1958-1962 averages %.0f per 100,000 per year, 1966-1969 averages%n",
                before);
        System.out.printf(L, "  %.0f, which is a factor of %.0f.%n", after, before / after);
        System.out.println();
        System.out.println("  No model is needed to see that. What follows is about the thirty");
        System.out.println("  five years before it, which are not one thing.");
    }

    // ---------------------------------------------------------------- 2

    /** The record's own measure of itself: how uneven consecutive years are. */
    static void whatTheRecordSaysAboutItself() {
        rule("2. The record is not in one state");
        System.out.println("Measles alternated: a heavy year, a light year. How heavy against how");
        System.out.println("light is a number the record carries and no model is needed for it.");
        System.out.println("Every year that is larger than both its neighbors, against the");
        System.out.println("smaller of them:");
        System.out.println();
        int[] years = peakYears();
        double[] ratios = peakRatios();
        for (int i = 0; i < years.length; ++i) {
            System.out.printf(L, "  %d   %5.1f  %s%n", years[i], ratios[i], bar(ratios[i] * 2.0));
        }
        System.out.println();
        System.out.printf(L, "  It runs from %.1f to %.1f and it is not noise around a mean: the%n",
                min(ratios), max(ratios));
        System.out.println("  1930s and 1940s sit low, the early 1950s are three times anything");
        System.out.println("  before them, and the record settles back after. Whatever the state");
        System.out.println("  of Massachusetts was doing, it changed, and section 3 asks what");
        System.out.println("  changed by asking what a model would need for each of them.");
    }

    // ---------------------------------------------------------------- 3

    /**
     * What the model settles into, shown rather than asserted, and how every
     * number below it is read off a trajectory.
     */
    static void whatTheModelSettlesInto() {
        rule("3. What the model settles into, and how that is measured");
        System.out.println("Everything below reads numbers off a model, so it is worth being");
        System.out.println("plain about which numbers those are.");
        System.out.println();
        System.out.println("The model is four compartments -- susceptible, exposed, infectious,");
        System.out.println("recovered -- with births replacing the susceptibles and a contact rate");
        System.out.println("that rises and falls once a year, because school terms do. Started");
        System.out.println("anywhere sensible and run forward, it forgets where it started. What");
        System.out.println("is left after that is the same whatever the start was, and is a");
        System.out.println("property of the equation rather than of anybody's guess at 1928.");
        System.out.println();
        System.out.println("That is a claim and it is cheap to check. Two runs, the second with");
        System.out.println("three times the infectious people and ten percent more susceptibles,");
        System.out.println("compared on the size of their next epidemic:");
        System.out.println();
        double amp = forcingFor(loudRatio());
        DVectorField f = seirField(amp);
        double[] startA = seirStart();
        double[] startB = { startA[0] * 1.10, startA[1] * 0.5, startA[2] * 3.0 };
        System.out.printf(L, "  %-14s %-16s %-16s %-14s%n", "after", "peak, start A",
                "peak, start B", "relative gap");
        for (double years : new double[] { 0.0, 10.0, 25.0, 50.0, 100.0, SETTLE_YEARS }) {
            double a = nextPeak(f, advance(f, startA, years), years);
            double b = nextPeak(f, advance(f, startB, years), years);
            System.out.printf(L, "  %-14s %-16s %-16s %-14s%n",
                    years == 0.0 ? "no time" : String.format(L, "%.0f years", years), people(a),
                    people(b), gapBetween(a, b));
        }
        System.out.println();
        System.out.printf(L,
                "  The two forget their difference at roughly a factor of a million per%n");
        System.out.printf(L,
                "  century, and after %.0f years to more digits than are worth printing --%n",
                SETTLE_YEARS);
        System.out.println("  past nine the next one is decided by the step sequence rather than");
        System.out.println("  by the equation, which is section 5's subject arriving early. So");
        System.out.printf(L, "  %.0f years is where this demo starts looking, and everything%n",
                SETTLE_YEARS);
        System.out.println("  below is read from the " + (int) WATCH_YEARS + " years after it.");
        System.out.println();
        System.out.println("Three quantities are read there, and each is a measurement on the");
        System.out.println("trajectory rather than a sample of it:");
        System.out.println();
        System.out.println("  a PEAK is where the infectious count stops rising -- the moment");
        System.out.println("  I' = 0 on the way down. It is found as an event, by bracketing the");
        System.out.println("  sign change and solving on the interpolant of the step that brackets");
        System.out.println("  it, so its time is as precise as the trajectory and owes nothing to");
        System.out.println("  where the step sizes happened to fall. An epidemic does not peak on");
        System.out.println("  a Monday and a weekly grid cannot see when it did.");
        System.out.println();
        System.out.println("  the TROUGH is the smallest the infectious count gets in those");
        System.out.printf(L, "  %d years. It is the number section 5 asks about, and the number%n",
                (int) WATCH_YEARS);
        System.out.println("  that decides whether the model is still describing a population.");
        System.out.println();
        System.out.println("  the PATTERN is what the sequence of peak sizes does. If every peak");
        System.out.println("  matches the one a year before it the model is running an annual");
        System.out.println("  cycle; if it matches the one two years before but not one year");
        System.out.println("  before, a two year cycle, and so on. Sizes and not spacings: a two");
        System.out.println("  year cycle stays locked to the school year, so its peaks are still");
        System.out.println("  twelve months apart and the mean interval reads 1.0000 either way.");
        System.out.println("  The statistic that looks like it measures the period cannot.");
    }

    /**
     * How far apart two epidemics are, relatively, floored where the answer
     * stops being reproducible.
     * <p>
     * The floor is this demo applying its own rule to its own table. Once two
     * runs agree to nine digits the tenth is decided by the step sequence, the
     * step sequence by a {@code Math.pow}, and that by which runtime is asked --
     * so the exact figure would differ between JDK 8 and JDK 25 and mean
     * nothing. What the row is for survives the floor: the gap keeps falling.
     */
    static String gapBetween(double a, double b) {
        double gap = Math.abs(a - b) / Math.max(a, b);
        return gap < 1.0e-09 ? "under 1e-09" : String.format(L, "%.0e", gap);
    }

    /** Runs the field forward without looking, to let the start be forgotten. */
    static double[] advance(DVectorField f, double[] start, double years) {
        if (years <= 0.0) {
            return start.clone();
        }
        return new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DOP853, f, 3),
                new StepController(RTOL, ATOL)).solve(0.0, start, years).finalState();
    }

    /** The size of the next epidemic peak after {@code t0}, as an event. */
    static double nextPeak(final DVectorField f, double[] y0, double t0) {
        OdeEvent peak = new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                double[] dydt = new double[3];
                f.valueAt(t, y, dydt);
                return dydt[2];
            }
        };
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, 3),
                new StepController(RTOL, ATOL),
                new Event[] { new Event(peak, Event.Direction.DECREASING, false, 1.0e-11, 1000) })
                        .solve(t0, y0, t0 + 4.0);
        double largest = 0.0;
        for (int i = 0; i < r.eventTimes.length; ++i) {
            largest = Math.max(largest, r.eventStates[i][2]);
        }
        return largest;
    }

    // ---------------------------------------------------------------- 4

    /** The record's ratio, turned into the one parameter that produces it. */
    static void theOneKnob() {
        rule("4. One knob, and the record sets it");
        System.out.println("Every rate in that model is a measured constant -- eight days latent,");
        System.out.println("five infectious, fifty years of life -- except one: the size of the");
        System.out.println("annual swing in contact, which nobody observes directly.");
        System.out.println();
        System.out.println("It is the only thing free, and the record measures it, because the");
        System.out.println("heavy-to-light ratio of section 2 grows with it. Settled and read as");
        System.out.println("section 3 describes:");
        System.out.println();
        System.out.printf(L, "  %-12s %-16s %-16s%n", "swing", "ratio, with the", "ratio, without");
        System.out.printf(L, "  %-12s %-16s %-16s%n", "", "latent class", "it");
        for (double amp = 0.06; amp <= 0.181; amp += 0.03) {
            System.out.printf(L, "  %-12.2f %-16.1f %-16.1f%n", amp, seirRatio(amp),
                    sirRatio(amp));
        }
        System.out.println();
        System.out.println("  The two columns are close and get closer, which is what lets the");
        System.out.println("  calibration below stand without the eight days of latency in it.");
        System.out.println("  Section 4 is where they stop being interchangeable.");
        System.out.println();
        System.out.println("Inverting it with Brent-Dekker, for the quietest stretch of the record");
        System.out.println("and for its loudest, through both models so that the disagreement is");
        System.out.println("visible rather than argued:");
        System.out.println();
        System.out.printf(L, "  %-22s %-10s %-14s %-14s%n", "stretch", "ratio", "swing, latent",
                "swing, without");
        System.out.printf(L, "  %-22s %-10.1f %-14.3f %-14.3f%n", "1930 to 1948", quietRatio(),
                forcingFor(quietRatio()), sirForcingFor(quietRatio()));
        System.out.printf(L, "  %-22s %-10.1f %-14.3f %-14.3f%n", "1952 and 1955", loudRatio(),
                forcingFor(loudRatio()), sirForcingFor(loudRatio()));
        System.out.println();
        System.out.printf(L,
                "  So the seasonal swing in contact roughly %s between the 1930s and the%n",
                "trebled");
        System.out.println("  early 1950s. That is a claim about schools and cities and not about");
        System.out.println("  arithmetic, and this demo does not defend it. What it does next is");
        System.out.println("  point out where on the model's own map those two numbers sit.");
    }

    // ---------------------------------------------------------------- 5

    /** Where the calibrated record sits, and what is immediately past it. */
    static void theEdgeBesideIt() {
        rule("5. Is that number the model's, or the step size's?");
        double loud = forcingFor(loudRatio());
        System.out.printf(L, "The loud stretch calibrates to a swing of %.3f. Before anything is%n",
                loud);
        System.out.println("read off a model there, one question has to be settled, and no solver");
        System.out.println("can settle it for you: is the number it just handed you a property of");
        System.out.println("your equation, or of the step sizes it happened to choose?");
        System.out.println();
        System.out.println("It costs one extra run to ask. Tighten the tolerance by a decade and");
        System.out.println("look at the same quantity again -- here the deepest trough of the");
        System.out.println("cycle, in people. If the two agree, the number is the model's.");
        System.out.println();
        System.out.printf(L, "  %-10s %-24s %-24s%n", "swing", "with the latent class",
                "without it");
        System.out.printf(L, "  %-10s %-11s %-12s %-11s %-12s%n", "", "trough", "same?", "trough",
                "same?");
        for (double amp : new double[] { 0.12, 0.167, 0.20, 0.24, 0.28 }) {
            double withLoose = troughAtTolerance(amp, RTOL, ATOL);
            double withTight = troughAtTolerance(amp, 1.0e-13, 1.0e-20);
            double withoutLoose = sirTroughAtTolerance(amp, RTOL, ATOL);
            double withoutTight = sirTroughAtTolerance(amp, 1.0e-13, 1.0e-20);
            boolean withAgrees = same(withLoose, withTight);
            boolean withoutAgrees = same(withoutLoose, withoutTight);
            System.out.printf(L, "  %-10.3f %-11s %-12s %-11s %-12s%s%n", amp,
                    withAgrees ? people(withLoose) : "--", withAgrees ? "yes" : "no",
                    withoutAgrees ? people(withoutLoose) : "--", withoutAgrees ? "yes" : "no",
                    Math.abs(amp - loud) < 0.005 ? "  <- the record" : "");
        }
        System.out.println();
        System.out.println("  The record sits where both models still answer. A little further on");
        System.out.println("  the one with the latent class stops answering while the one without");
        System.out.println("  it still does, and that gap is the whole of what the eight days buy");
        System.out.println("  here -- section 3 showed the calibration does not need them, and");
        System.out.println("  this is where they matter instead.");
        System.out.println();
        System.out.println("  What happens past the edge is that the cycle breaks up. A trough is");
        System.out.println("  then a property of the particular trajectory rather than of the");
        System.out.println("  attractor, so a tighter tolerance does not converge on an answer --");
        System.out.println("  there is no single answer to converge on -- and the runs also carry");
        System.out.println("  the model down to under one infective, which a population cannot do.");
        System.out.println();
        System.out.println("  Two honest warnings about that table. Where the edge is cannot be");
        System.out.println("  read off the parameter: between about 0.19 and 0.25 the model has");
        System.out.println("  orderly windows sitting inside the disorder, so a swing that answers");
        System.out.println("  can sit between two that do not, and stepping the parameter finely");
        System.out.println("  does not produce a boundary anyone could quote. That is an argument");
        System.out.println("  for running the check, not for knowing when to.");
        System.out.println();
        System.out.println("  And the pattern this demo names elsewhere is tested only against");
        System.out.println("  cycles of one, two and four years, so what it calls irregular means");
        System.out.println("  no more than that. The tolerance test above does not rest on it.");
    }

    /**
     * Whether two troughs a decade of tolerance apart are the same answer. A
     * factor of two is generous on purpose: what separates the two regimes here
     * is orders of magnitude, so nothing rests on where exactly the line is.
     */
    static boolean same(double loose, double tight) {
        return Math.max(loose, tight) / Math.min(loose, tight) < 2.0;
    }

    /** Whether the trough survives a decade of tolerance, for one model. */
    static boolean agrees(double amplitude, boolean withLatentClass) {
        double loose = withLatentClass ? troughAtTolerance(amplitude, RTOL, ATOL)
                : sirTroughAtTolerance(amplitude, RTOL, ATOL);
        double tight = withLatentClass ? troughAtTolerance(amplitude, 1.0e-13, 1.0e-20)
                : sirTroughAtTolerance(amplitude, 1.0e-13, 1.0e-20);
        return same(loose, tight);
    }

    /** The same trough for the model without a latent class. */
    static double sirTroughAtTolerance(double amplitude, double rtol, double atol) {
        DVectorField f = sirField(amplitude);
        double[] settled = new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DOP853, f, 2),
                new StepController(rtol, atol)).solve(0.0, sirStart(), SETTLE_YEARS).finalState();
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, 2),
                new StepController(rtol, atol)).solve(SETTLE_YEARS, settled,
                        SETTLE_YEARS + WATCH_YEARS);
        double trough = Double.MAX_VALUE;
        for (int i = 0; i < r.length; ++i) {
            trough = Math.min(trough, r.y[i][1]);
        }
        return trough;
    }

    // ---------------------------------------------------------------- 6

    /** The same question asked of the record, in three more states. */
    static void theSmallerStates() {
        rule("6. And the record says it too, in the states that are too small");
        System.out.println("Past the edge the model's epidemics run down to a trough of under one");
        System.out.println("infective, which a population cannot do: it carries nought or one, and");
        System.out.println("if it carries nought the disease is gone until somebody brings it");
        System.out.println("back. So the model's boundary should show up in the record as a");
        System.out.println("population below which measles stops cycling.");
        System.out.println();
        System.out.println("Four states in decreasing order of population, autocorrelation of the");
        System.out.println("log weekly incidence at one, two and three years:");
        System.out.println();
        double[][] series = Datasets.byDecreasingPopulation();
        String[] names = Datasets.byDecreasingPopulationNames();
        System.out.printf(L, "  %-16s %-10s %-10s %-10s %s%n", "state", "1 year", "2 years",
                "3 years", "no report");
        for (int i = 0; i < series.length; ++i) {
            double[] a = acfOf(series[i]);
            System.out.printf(L, "  %-16s %-10.3f %-10.3f %-10.3f %d weeks%n", names[i],
                    a[Datasets.WEEKS_PER_YEAR], a[2 * Datasets.WEEKS_PER_YEAR],
                    a[3 * Datasets.WEEKS_PER_YEAR], Datasets.weeksWithNoReport(series[i]));
        }
        System.out.println();
        System.out.println("  Monotone at all three lags, and the direction is the interesting");
        System.out.println("  part. It is not that the small state's record is noisier -- it is");
        System.out.println("  that the cycle is not there. Montana's incidence is the most violent");
        System.out.println("  of the four, because a hundred thousand is a large share of it, and");
        System.out.println("  its recurrence the least regular. Those are one fact, and they are");
        System.out.println("  Bartlett's critical community size: below some number of people");
        System.out.println("  measles does not cycle, it arrives, burns and goes out.");
        System.out.println();
        System.out.println("  The model said the same thing from the other side and neither knew");
        System.out.println("  about the other.");
    }

    // ---------------------------------------------------------------- 7

    /** The number the calibration did not need, and what it was worth. */
    static void theNumberAndTheThreshold() {
        rule("7. The number that decides the campaign");
        System.out.println("None of the above needed the basic reproduction number, and a public");
        System.out.println("health department needs nothing else: the fraction of a population");
        System.out.println("that has to be immune is 1 - 1/R0 and there is no way round it.");
        System.out.println();
        System.out.println("A single outbreak cannot give it. That was measured on Kermack and");
        System.out.println("McKendrick's 1905-06 Bombay plague counts while this demo was planned:");
        System.out.println("their model fits their curve to under one percent of 9043 deaths and");
        System.out.println("still leaves R0 anywhere between 1.07 and 3.20, because the shape of");
        System.out.println("one wave cannot separate how fast people infect each other from how");
        System.out.println("long they stay infectious.");
        System.out.println();
        System.out.println("Thirty-five years of recurrence can. An unforced model spirals into");
        System.out.println("its equilibrium and the period of that spiral is set by R0, so the");
        System.out.println("record's own recurrence is the second measurement the single wave was");
        System.out.println("missing. From the autocorrelation of section 5, Massachusetts recurs");
        System.out.printf(L, "at about two years, and the model's spiral says:%n");
        System.out.println();
        System.out.printf(L, "  %-10s %-16s%n", "R0", "spiral period");
        for (double r0 : new double[] { 8.0, 12.0, 17.0, 20.0 }) {
            System.out.printf(L, "  %-10.1f %-16.2f years%n", r0, ringingPeriod(r0));
        }
        System.out.println();
        System.out.printf(L, "  %-18s %-12s %-14s%n", "recurrence", "R0", "1 - 1/R0");
        for (double period : new double[] { 2.0, 2.5, 3.0 }) {
            double r0 = reproductionNumberFor(period);
            System.out.printf(L, "  %-18.1f %-12.1f %.0f%%%n", period, r0,
                    100.0 * (1.0 - 1.0 / r0));
        }
        System.out.println();
        System.out.printf(L,
                "  Read carefully: the spiral period is defined only to about %.0e years,%n",
                PERIOD_NOISE);
        System.out.println("  because it is an event time rather than an integral and the late");
        System.out.println("  peaks of a decaying oscillation are ill conditioned. Stepping R0 by");
        System.out.println("  a millionth moves it as much as stepping it by a thousandth does, so");
        System.out.println("  that is a floor and not a resolution, Brent-Dekker is asked for a");
        System.out.println("  hundredth rather than the 1e-10 it would report reaching, and R0 is");
        System.out.println("  worth a decimal place. The estimator also doubles whatever error the");
        System.out.println("  recurrence carries.");
        System.out.println();
        double r0 = reproductionNumberFor(2.0);
        System.out.printf(L, "  Taking the two year reading: R0 near %.0f, and %.0f%% of a%n", r0,
                100.0 * (1.0 - 1.0 / r0));
        System.out.println("  population immune. That is what the campaign had to reach, and here");
        System.out.println("  is what the record did:");
        System.out.println();
        double[] annual = annualTotals(Datasets.massachusetts());
        for (int year = 1960; year <= Datasets.MASSACHUSETTS_LAST_YEAR; ++year) {
            int i = year - Datasets.FIRST_YEAR;
            System.out.printf(L, "  %d %7.0f %s%s%n", year, annual[i], bar(annual[i] / 25.0),
                    year == Datasets.VACCINE_YEAR ? "   <- vaccine licensed" : "");
        }
        System.out.println();
        System.out.println("  Four years from licensure to a record with no epidemics in it, on a");
        System.out.println("  threshold that coverage did not reach in four years. Both are true.");
        System.out.println("  The threshold is a statement about the equilibrium and not about how");
        System.out.println("  long the equilibrium takes to arrive -- and section 5 is why they");
        System.out.println("  are compatible: what ends transmission is not incidence reaching");
        System.out.println("  zero but the trough falling below one person, and the vaccine did");
        System.out.println("  that long before it reached the threshold. The threshold is what");
        System.out.println("  keeps it away afterwards, and a record of four good years does not");
        System.out.println("  test it.");
    }

    // ------------------------------------------------------- the calibration

    /** Peak height ratio of the model with a latent class, at one forcing. */
    static double seirRatio(double amplitude) {
        Attractor at = attractorOf(seirField(amplitude), 3, seirStart(), amplitude);
        return at.period == 0 ? Double.NaN : at.largestPeak / at.smallestPeak;
    }

    /** Peak height ratio of the model without one. */
    static double sirRatio(double amplitude) {
        Attractor at = attractorOf(sirField(amplitude), 2, sirStart(), amplitude);
        return at.period == 0 ? Double.NaN : at.largestPeak / at.smallestPeak;
    }

    /**
     * The annual swing in contact that produces an observed heavy-to-light
     * ratio. Bracketed inside the window where the attractor is a cycle, because
     * outside it the ratio is not a function of the forcing at all.
     */
    static double forcingFor(final double observedRatio) {
        DFunction gap = new DFunction() {
            @Override
            public double apply(double amplitude) {
                return seirRatio(amplitude) - observedRatio;
            }
        };
        return RootFinder.brentDekker(0.055, 0.185, gap, 1.0e-04);
    }

    /** The same inversion through the model without a latent class. */
    static double sirForcingFor(final double observedRatio) {
        DFunction gap = new DFunction() {
            @Override
            public double apply(double amplitude) {
                return sirRatio(amplitude) - observedRatio;
            }
        };
        return RootFinder.brentDekker(0.055, 0.185, gap, 1.0e-04);
    }

    /** Years whose total exceeds both neighbors, before the vaccine. */
    static int[] peakYears() {
        double[] annual = annualTotals(Datasets.massachusetts());
        int[] out = new int[annual.length];
        int n = 0;
        for (int i = 1; i + 1 < annual.length; ++i) {
            int year = Datasets.FIRST_YEAR + i;
            if (year >= Datasets.VACCINE_YEAR) {
                break;
            }
            if (annual[i] > annual[i - 1] && annual[i] > annual[i + 1]) {
                out[n++] = year;
            }
        }
        int[] trimmed = new int[n];
        System.arraycopy(out, 0, trimmed, 0, n);
        return trimmed;
    }

    /** Each of those years against the smaller of its neighbors. */
    static double[] peakRatios() {
        double[] annual = annualTotals(Datasets.massachusetts());
        int[] years = peakYears();
        double[] out = new double[years.length];
        for (int k = 0; k < years.length; ++k) {
            int i = years[k] - Datasets.FIRST_YEAR;
            out[k] = annual[i] / Math.min(annual[i - 1], annual[i + 1]);
        }
        return out;
    }

    /** The mean ratio over the quiet stretch, 1930 to 1948. */
    static double quietRatio() {
        return meanRatioBetween(1930, 1948);
    }

    /** The mean ratio over the two loudest years the record has. */
    static double loudRatio() {
        return meanRatioBetween(1952, 1955);
    }

    private static double meanRatioBetween(int from, int to) {
        int[] years = peakYears();
        double[] ratios = peakRatios();
        double sum = 0.0;
        int n = 0;
        for (int i = 0; i < years.length; ++i) {
            if (years[i] >= from && years[i] <= to) {
                sum += ratios[i];
                ++n;
            }
        }
        return sum / n;
    }

    static double min(double[] v) {
        double m = Double.MAX_VALUE;
        for (int i = 0; i < v.length; ++i) {
            m = Math.min(m, v[i]);
        }
        return m;
    }

    static double max(double[] v) {
        double m = -Double.MAX_VALUE;
        for (int i = 0; i < v.length; ++i) {
            m = Math.max(m, v[i]);
        }
        return m;
    }

    /** The autocorrelation of one series, gaps filled and logarithms taken. */
    static double[] acfOf(double[] series) {
        return ACF.acf(logIncidence(fillSingleGaps(series)),
                3 * Datasets.WEEKS_PER_YEAR + 1);
    }

    /**
     * A fraction of the population as a count, to three significant figures.
     * <p>
     * Three and not more, because one unit in the last place of the
     * {@code Math.pow} inside the step controller is enough to move the last
     * digit of a count near a rounding boundary, and this page has to be the
     * same page on JDK 8 and on JDK 25.
     */
    static String people(double fraction) {
        double n = fraction * MASSACHUSETTS_PEOPLE;
        if (n < 1.0) {
            return "under one";
        }
        long exact = Math.round(n);
        long scale = 1L;
        while (exact / scale >= 1000L) {
            scale *= 10L;
        }
        return Long.toString(Math.round((double) exact / scale) * scale);
    }

    /**
     * What to print for a quantity taken off an irregular run, which is nothing.
     * A number read off a chaotic trajectory is a property of that trajectory
     * and not of the model, and it survives neither a change of tolerance nor a
     * change of runtime.
     */
    static String cell(Attractor at, double value) {
        return at.period == 0 ? "--" : people(value);
    }

    /** The trough of one attractor at a stated tolerance, to show what it is. */
    static double troughAtTolerance(double amplitude, double rtol, double atol) {
        DVectorField f = seirField(amplitude);
        double[] settled = new OdeIntegrator(new ExplicitRungeKutta(ButcherTableau.DOP853, f, 3),
                new StepController(rtol, atol)).solve(0.0, seirStart(), SETTLE_YEARS).finalState();
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, 3),
                new StepController(rtol, atol)).solve(SETTLE_YEARS, settled,
                        SETTLE_YEARS + WATCH_YEARS);
        double trough = Double.MAX_VALUE;
        for (int i = 0; i < r.length; ++i) {
            trough = Math.min(trough, r.y[i][2]);
        }
        return trough;
    }

    // ---------------------------------------------------------------- models

    /** A seasonally forced SIR over (S, I). */
    static DVectorField sirField(final double amplitude) {
        final double beta0 = R0_POSED * (GAMMA_SIR + MU);
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double beta = beta0 * (1.0 + amplitude * StrictMath.cos(2.0 * Math.PI * t));
                double infection = beta * y[0] * y[1];
                dydt[0] = MU - infection - MU * y[0];
                dydt[1] = infection - (GAMMA_SIR + MU) * y[1];
            }
        };
    }

    /** A seasonally forced SEIR over (S, E, I). */
    static DVectorField seirField(final double amplitude) {
        final double beta0 = betaFor(R0_POSED);
        return new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double beta = beta0 * (1.0 + amplitude * StrictMath.cos(2.0 * Math.PI * t));
                double infection = beta * y[0] * y[2];
                dydt[0] = MU - infection - MU * y[0];
                dydt[1] = infection - (SIGMA + MU) * y[1];
                dydt[2] = SIGMA * y[1] - (GAMMA + MU) * y[2];
            }
        };
    }

    /** The contact rate that gives a basic reproduction number of {@code r0}. */
    static double betaFor(double r0) {
        return r0 * (GAMMA + MU) * (SIGMA + MU) / SIGMA;
    }

    static double[] sirStart() {
        double sStar = 1.0 / R0_POSED;
        double iStar = MU * (R0_POSED - 1.0) / (R0_POSED * (GAMMA_SIR + MU));
        return new double[] { sStar * 1.02, iStar * 0.9 };
    }

    static double[] seirStart() {
        double sStar = 1.0 / R0_POSED;
        double iStar = MU * (R0_POSED - 1.0) / (R0_POSED * (GAMMA + MU));
        double eStar = (GAMMA + MU) * iStar / SIGMA;
        return new double[] { sStar * 1.02, eStar * 0.9, iStar };
    }

    // ---------------------------------------------------------------- attractor

    /** What one settled run turned out to be doing. */
    static final class Attractor {

        final int peaks;
        final double meanGap;
        final double smallestPeak;
        final double largestPeak;
        final double trough;

        /**
         * How many years the pattern of peaks takes to repeat, or zero if it
         * does not repeat inside the window. <b>This is read off the peak
         * heights and not off their spacing</b>, which is the whole trap here: a
         * two year cycle locks its peaks to the forcing, so they stay one year
         * apart and only their size alternates. A mean interval of exactly
         * 1.0000 is as true of the two year cycle as of the one year cycle and
         * says nothing about which is running.
         */
        final int period;

        Attractor(int peaks, double meanGap, double smallestPeak, double largestPeak, double trough,
                int period) {
            this.peaks = peaks;
            this.meanGap = meanGap;
            this.smallestPeak = smallestPeak;
            this.largestPeak = largestPeak;
            this.trough = trough;
            this.period = period;
        }

        String describe() {
            switch (period) {
            case 1:
                return "annual";
            case 2:
                return "biennial";
            case 0:
                return "irregular";
            default:
                return period + " years";
            }
        }
    }

    /**
     * The number of years the sequence of peak heights takes to repeat, tried
     * against one, two and four, or zero if none of them holds.
     */
    static int periodOf(double[] heights) {
        for (int p : new int[] { 1, 2, 4 }) {
            if (heights.length < 2 * p + 1) {
                continue;
            }
            boolean holds = true;
            for (int i = p; i < heights.length; ++i) {
                double a = heights[i];
                double b = heights[i - p];
                if (Math.abs(a - b) > 1.0e-3 * Math.max(a, b)) {
                    holds = false;
                    break;
                }
            }
            if (holds) {
                return p;
            }
        }
        return 0;
    }

    /**
     * Settles the field onto its attractor and reports what the attractor is.
     * The peaks are events, so their times come from the interpolant and not
     * from wherever the step sizes happened to land.
     */
    static Attractor attractorOf(final DVectorField f, final int dimension, double[] start,
            double amplitude) {
        final int last = dimension - 1;
        double[] settled = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, dimension),
                new StepController(RTOL, ATOL)).solve(0.0, start, SETTLE_YEARS).finalState();

        OdeEvent peak = new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                double[] dydt = new double[dimension];
                f.valueAt(t, y, dydt);
                return dydt[last];
            }
        };
        Event event = new Event(peak, Event.Direction.DECREASING, false, 1.0e-11, 100000);
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, dimension),
                new StepController(RTOL, ATOL), new Event[] { event }).solve(SETTLE_YEARS, settled,
                        SETTLE_YEARS + WATCH_YEARS);

        double trough = Double.MAX_VALUE;
        for (int i = 0; i < r.length; ++i) {
            trough = Math.min(trough, r.y[i][last]);
        }
        int n = r.eventTimes.length;
        double[] heights = new double[n];
        double smallest = Double.MAX_VALUE;
        double largest = 0.0;
        for (int i = 0; i < n; ++i) {
            heights[i] = r.eventStates[i][last];
            smallest = Math.min(smallest, heights[i]);
            largest = Math.max(largest, heights[i]);
        }
        double meanGap = n < 2 ? Double.NaN : (r.eventTimes[n - 1] - r.eventTimes[0]) / (n - 1);
        return new Attractor(n, meanGap, smallest, largest, trough, periodOf(heights));
    }

    /** The period of the damped oscillation about the endemic equilibrium. */
    static double ringingPeriod(double r0) {
        final double beta = betaFor(r0);
        DVectorField f = new DVectorField() {
            @Override
            public void valueAt(double t, double[] y, double[] dydt) {
                double infection = beta * y[0] * y[2];
                dydt[0] = MU - infection - MU * y[0];
                dydt[1] = infection - (SIGMA + MU) * y[1];
                dydt[2] = SIGMA * y[1] - (GAMMA + MU) * y[2];
            }
        };
        double sStar = 1.0 / r0;
        double iStar = MU * (r0 - 1.0) / (r0 * (GAMMA + MU));
        double eStar = (GAMMA + MU) * iStar / SIGMA;
        double[] y0 = { sStar, eStar * 1.05, iStar * 1.05 };

        OdeEvent peak = new OdeEvent() {
            @Override
            public double valueAt(double t, double[] y) {
                return SIGMA * y[1] - (GAMMA + MU) * y[2];
            }
        };
        Event event = new Event(peak, Event.Direction.DECREASING, false, 1.0e-12, 1000);
        OdeIntegrator.Result r = new OdeIntegrator(
                new ExplicitRungeKutta(ButcherTableau.DOP853, f, 3), new StepController(RTOL, ATOL),
                new Event[] { event }).solve(0.0, y0, 120.0);
        int n = r.eventTimes.length;
        if (n < 3) {
            return Double.NaN;
        }
        return (r.eventTimes[n - 1] - r.eventTimes[1]) / (n - 2);
    }

    /**
     * The noise floor of {@link #ringingPeriod(double)}, in years, measured
     * rather than assumed: stepping {@code R0} by {@code 1e-06} moves the period
     * by {@code 6.5e-04}, and stepping it by {@code 1e-03} moves it by the same
     * amount. That is a floor and not a derivative.
     * <p>
     * It is there because the period is an <b>event time</b> and not an
     * integral. A quantity integrated over a trajectory is smoothed by the
     * integration -- the same question asked of a least squares fit on the
     * Bombay plague data found the fitted parameters tracking the solver
     * tolerance down to {@code 1e-09} with nothing rough about them. A peak
     * location on a decaying oscillation has no such smoothing: the last peaks
     * are small, their times are ill conditioned, and a different step sequence
     * puts them somewhere slightly else.
     */
    static final double PERIOD_NOISE = 6.5e-04;

    /**
     * Inverts {@link #ringingPeriod(double)}.
     * <p>
     * The bracket is asked for no more than the function has. Near the answer
     * the period moves about {@code 0.07} years per unit of {@code R0}, so
     * {@link #PERIOD_NOISE} is about a hundredth of {@code R0} and that is the
     * tolerance. Asking Brent-Dekker for {@code 1e-10} converges on the noise
     * and returns a root whose residual is four orders larger than the
     * tolerance it was given -- which is not the solver's fault and is what a
     * convergence flag cannot tell you.
     */
    static double reproductionNumberFor(final double period) {
        DFunction gap = new DFunction() {
            @Override
            public double apply(double r0) {
                return ringingPeriod(r0) - period;
            }
        };
        return RootFinder.brentDekker(3.0, 22.0, gap, 1.0e-02);
    }

    // ---------------------------------------------------------------- data

    /** The part of a series that lies before the vaccine year. */
    static double[] preVaccine(double[] weekly) {
        int weeks = (Datasets.VACCINE_YEAR - Datasets.FIRST_YEAR) * Datasets.WEEKS_PER_YEAR;
        double[] out = new double[Math.min(weeks, weekly.length)];
        System.arraycopy(weekly, 0, out, 0, out.length);
        return out;
    }

    /**
     * Fills a week with no report by interpolating its neighbors. Almost every
     * gap in this record is one week wide, so this is a small thing, but it is a
     * choice and the demo makes it visibly rather than by dropping the weeks and
     * pretending the series is still evenly spaced -- which is what an
     * autocorrelation would then be computing.
     */
    static double[] fillSingleGaps(double[] series) {
        double[] out = series.clone();
        for (int i = 0; i < out.length; ++i) {
            if (!Double.isNaN(out[i])) {
                continue;
            }
            int before = i - 1;
            int after = i + 1;
            while (after < out.length && Double.isNaN(series[after])) {
                ++after;
            }
            if (before >= 0 && after < out.length) {
                double span = after - before;
                out[i] = out[before] + (series[after] - out[before]) * (i - before) / span;
            } else {
                out[i] = 0.0;
            }
        }
        return out;
    }

    /** The logarithm of the incidence, floored so that a reported zero survives it. */
    static double[] logIncidence(double[] series) {
        double[] out = new double[series.length];
        for (int i = 0; i < out.length; ++i) {
            out[i] = StrictMath.log(Math.max(series[i], 0.01));
        }
        return out;
    }

    /** Sums each calendar year, skipping the weeks that carry no report. */
    static double[] annualTotals(double[] weekly) {
        int years = weekly.length / Datasets.WEEKS_PER_YEAR;
        double[] out = new double[years];
        for (int y = 0; y < years; ++y) {
            double sum = 0.0;
            for (int w = 0; w < Datasets.WEEKS_PER_YEAR; ++w) {
                double v = weekly[y * Datasets.WEEKS_PER_YEAR + w];
                if (!Double.isNaN(v)) {
                    sum += v;
                }
            }
            out[y] = sum;
        }
        return out;
    }

    static double mean(double[] values, int from, int to) {
        double sum = 0.0;
        for (int i = from; i <= to; ++i) {
            sum += values[i];
        }
        return sum / (to - from + 1);
    }

    // ---------------------------------------------------------------- page

    static void rule(String title) {
        System.out.println();
        System.out.println("=========================================================================");
        System.out.println(title);
        System.out.println("=========================================================================");
        System.out.println();
    }

    static String bar(double value) {
        int n = (int) Math.round(Math.abs(value));
        if (n > 45) {
            n = 45;
        }
        StringBuilder sb = new StringBuilder();
        char c = value < 0.0 ? '-' : '#';
        for (int i = 0; i < n; ++i) {
            sb.append(c);
        }
        return sb.toString();
    }
}
