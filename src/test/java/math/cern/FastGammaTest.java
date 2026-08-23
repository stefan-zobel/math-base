package math.cern;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * {@link FastGamma#logGamma(double)} used to answer every argument at or below
 * zero with the literal {@code -999}, which reads as a perfectly ordinary
 * logarithm -- {@code Math.exp(-999)} is {@code 0.0}, so a caller got a density
 * of zero rather than an error. Above zero it shifted the argument up to eleven
 * and subtracted the logarithm of the shift, which for {@code x = 1} means
 * taking two values of about {@code 15.1} apart to get zero: the absolute error
 * was about five ulps everywhere below thirteen, and since {@code log(gamma)}
 * has zeros at {@code x = 1} and {@code x = 2} the relative error there reached
 * {@code 1e-10}. And {@code logGamma(POSITIVE_INFINITY)} was {@code NaN},
 * because {@code (inf - 0.5) * log(inf) - inf} is.
 * <p>
 * The oracles here are all independent of any gamma implementation: closed
 * forms at half-integers and integers, the functional equation
 * {@code gamma(x + 1) = x gamma(x)}, and for {@link FastGamma#logGammaRatio}
 * the asymptotic {@code 0.5 ln z - 1/(8z)}, whose own error is
 * {@code O(1/z^3)}.
 */
public class FastGammaTest {

    /** log(gamma(x)) at points where it has a closed form */
    private static final double[][] EXACT = { { 0.5, 0.5 * Math.log(Math.PI) }, { 1.0, 0.0 },
            { 2.0, 0.0 }, { 3.0, Math.log(2.0) }, { 4.0, Math.log(6.0) }, { 5.0, Math.log(24.0) },
            { 6.0, Math.log(120.0) }, { 7.0, Math.log(720.0) }, { 11.0, Math.log(3628800.0) },
            { 12.0, Math.log(39916800.0) }, { 13.0, Math.log(479001600.0) } };

    /**
     * The two zeros must come back as exactly zero -- before the change they
     * were {@code 1.78e-15}. The other closed forms are only asked for one
     * ulp, because {@code Math.log(720.0)} is itself a rounded value and
     * demanding equality would demand that the two round the same way; they
     * were off by up to {@code 3.6e-15}, which is five ulps.
     */
    @Test
    public void theClosedFormsAreExact() {
        assertEquals("logGamma(1.0)", 0.0, FastGamma.logGamma(1.0), 0.0);
        assertEquals("logGamma(2.0)", 0.0, FastGamma.logGamma(2.0), 0.0);
        for (int i = 0; i < EXACT.length; ++i) {
            double want = EXACT[i][1];
            assertEquals("logGamma(" + EXACT[i][0] + ")", want, FastGamma.logGamma(EXACT[i][0]),
                    Math.ulp(want));
        }
    }

    /**
     * {@code logGamma(x + 1) - logGamma(x) - log(x)} is zero for every
     * {@code x}, whichever way the two are computed. This walks the positive
     * line including the handover at thirteen.
     */
    @Test
    public void theFunctionalEquationHoldsAboveZero() {
        double worst = 0.0;
        double at = 0.0;
        for (int i = 1; i <= 40000; ++i) {
            double x = 30.0 * i / 40000.0;
            double d = Math.abs(FastGamma.logGamma(x + 1.0) - FastGamma.logGamma(x) - Math.log(x));
            if (d > worst) {
                worst = d;
                at = x;
            }
        }
        assertTrue("worst residual " + worst + " at x = " + at, worst < 1.0e-13);
    }

    /** The same across zero, where the reflection formula takes over. */
    @Test
    public void theFunctionalEquationHoldsBelowZero() {
        double worst = 0.0;
        double at = 0.0;
        for (int i = 0; i <= 40000; ++i) {
            double x = -6.0 + 5.9 * i / 40000.0;
            if (Math.abs(x - Math.rint(x)) < 5.0e-3) {
                continue; // too close to a pole for the difference to mean anything
            }
            double d = Math.abs(FastGamma.logGamma(x + 1.0) - FastGamma.logGamma(x) - Math.log(Math.abs(x)));
            if (d > worst) {
                worst = d;
                at = x;
            }
        }
        assertTrue("worst residual " + worst + " at x = " + at, worst < 1.0e-13);
    }

    /**
     * The negative half plane carries the logarithm of the magnitude, because
     * gamma alternates sign there. Checked against the downward recurrence
     * from a positive argument.
     */
    @Test
    public void theNegativeHalfPlaneIsTheLogOfTheMagnitude() {
        double[] xs = { -0.5, -1.5, -2.5, -3.5, -4.7, -0.1, -5.25, -40.25, -1.0e-8 };
        for (int i = 0; i < xs.length; ++i) {
            double x = xs[i];
            double n = Math.ceil(-x);
            double v = Math.exp(FastGamma.logGamma(x + n));
            for (double k = 0.0; k < n; ++k) {
                v /= (x + k);
            }
            assertEquals("logGamma(" + x + ")", Math.log(Math.abs(v)), FastGamma.logGamma(x),
                    1.0e-12 * Math.max(1.0, Math.abs(Math.log(Math.abs(v)))));
        }
    }

    /**
     * A magnitude that underflows on the way is the case the reflection has to
     * be written for: {@code q * sin(pi q)} is below the smallest normal
     * double from about {@code q = 1e-162}, while each factor on its own is
     * not. For a tiny {@code q}, {@code gamma(-q)} is about {@code -1/q}.
     */
    @Test
    public void theReflectionSurvivesATinyArgument() {
        for (double q : new double[] { 1.0e-100, 1.0e-200, 1.0e-300, 1.0e-320 }) {
            assertEquals("logGamma(" + (-q) + ")", -Math.log(q), FastGamma.logGamma(-q), 1.0e-12);
        }
    }

    /** The poles, and the two ends of the line. */
    @Test
    public void theEdgesOfTheDomainAreAnswered() {
        for (double pole : new double[] { 0.0, -1.0, -2.0, -3.0, -100.0 }) {
            assertEquals("pole at " + pole, Double.POSITIVE_INFINITY, FastGamma.logGamma(pole), 0.0);
            assertTrue("gamma has no value at " + pole, Double.isNaN(FastGamma.gamma(pole)));
        }
        assertEquals("logGamma(+inf)", Double.POSITIVE_INFINITY,
                FastGamma.logGamma(Double.POSITIVE_INFINITY), 0.0);
        assertTrue("logGamma(-inf) has no limit to converge to",
                Double.isNaN(FastGamma.logGamma(Double.NEGATIVE_INFINITY)));
        assertTrue("logGamma(NaN)", Double.isNaN(FastGamma.logGamma(Double.NaN)));
        assertTrue("gamma(NaN)", Double.isNaN(FastGamma.gamma(Double.NaN)));
    }

    /** {@code gamma} has to put the sign back on below zero. */
    @Test
    public void gammaAlternatesSignBelowZero() {
        double[][] signed = { { -0.5, -1.0 }, { -1.5, 1.0 }, { -2.5, -1.0 }, { -3.5, 1.0 }, { -4.5, -1.0 },
                { -0.1, -1.0 }, { -1.9, 1.0 }, { -4.7, -1.0 } };
        for (int i = 0; i < signed.length; ++i) {
            double g = FastGamma.gamma(signed[i][0]);
            assertEquals("sign of gamma(" + signed[i][0] + ") = " + g, signed[i][1], Math.signum(g), 0.0);
        }
        // and the functional equation holds for the signed value too
        double worst = 0.0;
        for (int i = 0; i <= 20000; ++i) {
            double x = -6.0 + 5.9 * i / 20000.0;
            if (Math.abs(x - Math.rint(x)) < 5.0e-3) {
                continue;
            }
            double lhs = FastGamma.gamma(x + 1.0);
            double rhs = x * FastGamma.gamma(x);
            worst = Math.max(worst, Math.abs(lhs - rhs) / Math.abs(rhs));
        }
        assertTrue("gamma(x+1) = x gamma(x), worst relative " + worst, worst < 1.0e-13);
    }

    /**
     * {@link FastGamma#logGammaRatio(double, double)} against the asymptotic
     * {@code 0.5 ln z - 1/(8z)} for {@code delta = 1/2}, whose own error is
     * {@code O(1/z^3)} and therefore below the last bit from {@code z = 1e8}
     * up. Subtracting two logarithms gives {@code 20.0} at {@code z = 5e14}
     * and {@code 0.0} from {@code 5e15}, where the answers are 16.92 and
     * 18.07.
     */
    @Test
    public void theRatioSurvivesALargeArgument() {
        for (double z : new double[] { 1.0e8, 1.0e10, 5.0e13, 5.0e14, 5.0e15, 5.0e16, 1.0e18 }) {
            double asymptotic = 0.5 * Math.log(z) - 1.0 / (8.0 * z);
            assertEquals("logGammaRatio(" + z + ", 0.5)", asymptotic, FastGamma.logGammaRatio(z, 0.5),
                    1.0e-14 * asymptotic);
        }
    }

    /**
     * And below the Stirling threshold, where it falls back to the plain
     * difference, against the closed forms.
     */
    @Test
    public void theRatioIsRightForSmallArguments() {
        for (int i = 0; i < EXACT.length; ++i) {
            for (int j = 0; j < EXACT.length; ++j) {
                double x = EXACT[i][0];
                double y = EXACT[j][0];
                assertEquals("logGammaRatio(" + x + ", " + (y - x) + ")", EXACT[j][1] - EXACT[i][1],
                        FastGamma.logGammaRatio(x, y - x), 1.0e-14);
            }
        }
        assertEquals("a zero offset is a zero ratio", 0.0, FastGamma.logGammaRatio(7.0, 0.0), 0.0);
    }
}
