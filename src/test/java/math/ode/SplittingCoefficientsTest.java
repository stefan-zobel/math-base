package math.ode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

/**
 * The coefficients of the symplectic methods, checked against the one step map
 * they induce on {@code q'' = -q}.
 * <p>
 * That map is the whole test apparatus and it is exact: for the harmonic
 * oscillator a drift is
 * <code>[[1, a h], [0, 1]]</code> and a kick is <code>[[1, 0], [-b h, 1]]</code>,
 * so one step of any method here is a {@code 2 x 2} matrix whose entries are
 * polynomials in {@code h}, and the flow it approximates is
 * <code>[[cos h, sin h], [-sin h, cos h]]</code>. Three things fall out of the
 * comparison without a single trajectory being integrated:
 * <ul>
 * <li><b>Symplecticity</b>, as the determinant. Measured over forty step sizes
 * spanning six decades, it leaves one by at most {@code 3.3e-16} -- rounding,
 * and nothing else. This is structural: a drift and a kick each have
 * determinant one whatever their coefficient is, so the property survives even
 * a mistyped coefficient. Only the <em>order</em> does not.</li>
 * <li><b>The order</b>, as the power of {@code h} at which the two matrices
 * part company: {@code 3.00} for Verlet, {@code 5.00} for the three fourth
 * order methods and {@code 7.00} for the sixth order one.</li>
 * <li><b>The error constant</b>, as the coefficient of that power, and this is
 * the number that decides which method is worth using. Verlet's is {@code 1/6}
 * exactly. The three fourth order methods come out at {@code 1.0e-01}
 * (Yoshida's triple jump), {@code 4.5e-03} (Suzuki's five-fold composition) and
 * {@code 1.6e-05} (Blanes-Moan) -- <b>four orders of magnitude at the same
 * order of accuracy</b>, and the reason order alone is not an argument for
 * anything.</li>
 * </ul>
 * <p>
 * Turned into what a method costs to reach a given accuracy, that ranks them
 * without ambiguity among the fourth order three: Blanes-Moan reaches any
 * accuracy for about a fifth of the triple jump's cost and a third of Suzuki's,
 * whatever the accuracy is. The sixth order method cannot be put on the same
 * scale, since a higher order gains as the accuracy tightens -- it is already
 * the cheapest at {@code 1e-6} and pulls away from there.
 * <p>
 * <b>The Blanes-Moan coefficients are the only copied numbers in the package.</b>
 * Everything else is computed from cube roots. The order check above is what
 * makes copying them safe: a single wrong digit does not perturb the answer, it
 * collapses the order, and the test stops the build.
 */
public final class SplittingCoefficientsTest {

    /** Measured worst case 3.3e-16 for the determinant and the sums. */
    private static final double EXACT = 1.0e-14;

    private static SplittingCoefficients[] all() {
        return new SplittingCoefficients[] { SplittingCoefficients.VERLET, SplittingCoefficients.YOSHIDA_4,
                SplittingCoefficients.SUZUKI_4, SplittingCoefficients.BLANES_MOAN_4,
                SplittingCoefficients.BLANES_MOAN_6 };
    }

    /**
     * A step size small enough for the error constant to have converged and
     * large enough that the difference of the two matrices is still above the
     * rounding. The window closes as the order rises: at order six the error at
     * {@code h = 0.05} is already {@code 1.6e-16}.
     */
    private static double settled(SplittingCoefficients c) {
        return (c.order() >= 6) ? 0.1 : 0.0125;
    }

    @Test
    public void testTheDriftsAndTheKicksEachSumToOne() {
        SplittingCoefficients[] methods = all();
        for (int q = 0; q < methods.length; ++q) {
            assertEquals(methods[q].name(), 1.0, sum(methods[q].a()), EXACT);
            assertEquals(methods[q].name(), 1.0, sum(methods[q].b()), EXACT);
        }
    }

    @Test
    public void testEveryMethodReadsTheSameForwardsAndBackwards() {
        SplittingCoefficients[] methods = all();
        for (int q = 0; q < methods.length; ++q) {
            assertTrue(methods[q].name() + " must be palindromic", methods[q].isPalindromic());
            assertEquals("a palindromic method has an even order", 0, methods[q].order() % 2);
        }
    }

    /**
     * A method that begins with a drift costs one evaluation per kick; one that
     * begins with a kick costs one less, because the kick that opens a step and
     * the one that closed the step before fall at the same place and the same
     * time.
     */
    @Test
    public void testWhatAStepCostsAndWhichMethodsGetOneForFree() {
        assertEquals(1, SplittingCoefficients.VERLET.evaluations());
        assertEquals(3, SplittingCoefficients.YOSHIDA_4.evaluations());
        assertEquals(5, SplittingCoefficients.SUZUKI_4.evaluations());
        assertEquals("seven kicks, six evaluations", 6, SplittingCoefficients.BLANES_MOAN_4.evaluations());
        assertEquals("twelve kicks, eleven evaluations", 11,
                SplittingCoefficients.BLANES_MOAN_6.evaluations());

        SplittingCoefficients[] methods = all();
        for (int q = 0; q < methods.length; ++q) {
            int counted = 0;
            double[] b = methods[q].b();
            for (int i = 0; i < b.length; ++i) {
                if (b[i] != 0.0) {
                    ++counted;
                }
            }
            int shared = methods[q].isVelocityFirst() ? 1 : 0;
            assertEquals(methods[q].name(), counted - shared, methods[q].evaluations());
            assertEquals("a method begins with a kick exactly when its first drift is zero",
                    methods[q].isVelocityFirst(), methods[q].a()[0] == 0.0);
        }
        assertFalse(SplittingCoefficients.VERLET.isVelocityFirst());
        assertFalse(SplittingCoefficients.SUZUKI_4.isVelocityFirst());
        assertTrue(SplittingCoefficients.BLANES_MOAN_4.isVelocityFirst());
        assertTrue(SplittingCoefficients.BLANES_MOAN_6.isVelocityFirst());
    }

    /**
     * A symmetric composition of a symmetric second order method is of order
     * four exactly when the cubes of its fractions cancel. Both closed forms
     * here satisfy that identically -- Yoshida because
     * <code>w0 = -2^(1/3) w1</code> and the cube of {@code 2^(1/3)} is two,
     * Suzuki because <code>1 - 4w = -4^(1/3) w</code> -- so what is measured
     * below is only that the arithmetic in the class agrees with the algebra.
     */
    @Test
    public void testTheFourthOrderCompositionsCancelTheCubesOfTheirFractions() {
        assertEquals("Yoshida", 0.0, sumOfCubes(SplittingCoefficients.YOSHIDA_4.b()), EXACT);
        assertEquals("Suzuki", 0.0, sumOfCubes(SplittingCoefficients.SUZUKI_4.b()), EXACT);
        assertEquals("Verlet is a single step and cannot cancel anything", 1.0,
                sumOfCubes(SplittingCoefficients.VERLET.b()), EXACT);
    }

    /**
     * The property the whole round exists for, and it holds whatever the
     * coefficients are: a drift and a kick have determinant one apiece, so
     * every product of them does.
     */
    @Test
    public void testTheOneStepMapPreservesArea() {
        SplittingCoefficients[] methods = all();
        for (int q = 0; q < methods.length; ++q) {
            double worst = 0.0;
            for (int i = 0; i < 40; ++i) {
                double h = 2.0 / Math.pow(2.0, i / 2.0);
                double[] m = map(methods[q], h);
                worst = Math.max(worst, Math.abs(m[0] * m[3] - m[1] * m[2] - 1.0));
            }
            assertEquals(methods[q].name() + " must preserve area at every step size", 0.0, worst, EXACT);
        }
    }

    /**
     * The window in which the order can be measured closes from both ends and
     * closes sooner the better the method is: a small error constant reaches
     * the rounding at a larger step. Rather than pick a range per method, halve
     * until the difference of the two matrices would be lost in the rounding
     * and use what is left, dropping the first two as pre-asymptotic.
     */
    @Test
    public void testTheOrderIsWhatEachMethodClaims() {
        SplittingCoefficients[] methods = all();
        for (int q = 0; q < methods.length; ++q) {
            double[] errors = usableMapErrors(methods[q]);
            assertTrue(methods[q].name() + " leaves only " + errors.length + " usable step sizes",
                    errors.length >= 4);
            for (int i = 2; i < errors.length; ++i) {
                double slope = Math.log(errors[i - 1] / errors[i]) / Math.log(2.0);
                assertEquals(methods[q].name() + ", halving " + i, methods[q].order() + 1.0, slope, 0.05);
            }
        }
    }

    /**
     * The two transcribed methods, and the only place in this package where a
     * digit was copied from a table rather than computed. A single wrong digit
     * shows up as the order collapsing, so the check above is the check on the
     * transcription and this one only records what came out.
     */
    @Test
    public void testTheBlanesMoanCoefficientsAreWhatTheirTableSays() {
        SplittingCoefficients four = SplittingCoefficients.BLANES_MOAN_4;
        assertEquals(4, four.order());
        assertEquals("six drifts and seven kicks", 7, four.stages());
        assertEquals(0.0, four.a()[0], 0.0);
        assertEquals("the middle kick is what the sum leaves", 0.119524194013151, four.b()[3], 1.0e-15);
        assertEquals("and the middle drift likewise", -0.350171622895351, four.a()[3], 1.0e-15);
        assertEquals(1.55e-05, errorConstant(four), 1.0e-7);

        SplittingCoefficients six = SplittingCoefficients.BLANES_MOAN_6;
        assertEquals(6, six.order());
        assertEquals("eleven drifts and twelve kicks", 12, six.stages());
        assertEquals(0.0, six.a()[0], 0.0);
        assertEquals("here the middle kick is repeated", six.b()[5], six.b()[6], 0.0);
        assertEquals(0.236669924786931, six.b()[5], 1.0e-15);
        assertEquals(0.204777054291470, six.a()[6], 1.0e-15);
        assertEquals(2.22e-07, errorConstant(six), 1.0e-8);
    }

    /**
     * Verlet's leading error term is <code>h^3 / 6</code>, which is what the
     * expansion of the map says it should be, and a pleasant thing to be able
     * to check against a whole number.
     */
    @Test
    public void testTheErrorConstantOfVerletIsExactlyOneSixth() {
        assertEquals(1.0 / 6.0, errorConstant(SplittingCoefficients.VERLET), 1.0e-4);
    }

    /**
     * The answer to "is there anything better than Yoshida of order four", and
     * it is not close. The three fourth order methods have error constants of
     * {@code 1.0e-01}, {@code 4.5e-03} and {@code 1.6e-05}; after paying for
     * their extra evaluations, Suzuki reaches a given accuracy for about
     * three quarters of Yoshida's cost and Blanes-Moan for about a fifth. The
     * ratios do not depend on the accuracy, because the order is the same.
     */
    @Test
    public void testTheThreeFourthOrderMethodsRankedByWhatTheyCost() {
        double yoshida = errorConstant(SplittingCoefficients.YOSHIDA_4);
        double suzuki = errorConstant(SplittingCoefficients.SUZUKI_4);
        double blanesMoan = errorConstant(SplittingCoefficients.BLANES_MOAN_4);
        assertEquals(0.1042, yoshida, 5.0e-4);
        assertEquals(0.0045, suzuki, 5.0e-4);
        assertEquals(1.55e-05, blanesMoan, 1.0e-7);
        assertTrue("four orders of magnitude at the same order of accuracy",
                yoshida / blanesMoan > 5000.0);

        double[] accuracies = { 1.0e-6, 1.0e-9, 1.0e-12 };
        for (int i = 0; i < accuracies.length; ++i) {
            double byYoshida = cost(SplittingCoefficients.YOSHIDA_4, accuracies[i]);
            double bySuzuki = cost(SplittingCoefficients.SUZUKI_4, accuracies[i]);
            double byBlanesMoan = cost(SplittingCoefficients.BLANES_MOAN_4, accuracies[i]);
            assertTrue("Suzuki beats the triple jump at " + accuracies[i], bySuzuki < byYoshida);
            assertTrue("and Blanes-Moan beats them both", byBlanesMoan < bySuzuki);
            assertEquals("by about a factor of four and a half, whatever the accuracy", 4.53,
                    byYoshida / byBlanesMoan, 0.1);
        }
    }

    /**
     * The sixth order method cannot be ranked against the others by one number,
     * because a higher order gains as the accuracy demanded tightens. What can
     * be said is that it is already cheaper at a millionth and pulls further
     * ahead from there.
     */
    @Test
    public void testTheSixthOrderMethodGainsAsTheAccuracyTightens() {
        double loose = cost(SplittingCoefficients.BLANES_MOAN_4, 1.0e-6)
                / cost(SplittingCoefficients.BLANES_MOAN_6, 1.0e-6);
        double tight = cost(SplittingCoefficients.BLANES_MOAN_4, 1.0e-12)
                / cost(SplittingCoefficients.BLANES_MOAN_6, 1.0e-12);
        assertTrue("already ahead at a millionth: " + loose, loose > 1.0);
        assertTrue("and further ahead at a millionth of a millionth: " + tight, tight > 3.0 * loose);

        // and it is the cheapest of the five at every accuracy tried
        SplittingCoefficients[] methods = all();
        double[] accuracies = { 1.0e-6, 1.0e-9, 1.0e-12 };
        for (int i = 0; i < accuracies.length; ++i) {
            double best = cost(SplittingCoefficients.BLANES_MOAN_6, accuracies[i]);
            for (int q = 0; q < methods.length; ++q) {
                if (methods[q] != SplittingCoefficients.BLANES_MOAN_6) {
                    assertTrue(methods[q].name() + " at " + accuracies[i],
                            cost(methods[q], accuracies[i]) > best);
                }
            }
        }
    }

    /**
     * Why the fourth order methods differ so much: every method above order two
     * runs part of its step backwards, and the one-norm is how far. Verlet
     * travels one unit of time to advance one; Yoshida travels {@code 4.40},
     * Suzuki {@code 2.32}, and Blanes-Moan only {@code 1.16}. The order of the
     * error constants is the order of these.
     */
    @Test
    public void testTheOneNormIsHowFarTheMethodTravelsToAdvanceOneStep() {
        assertEquals(1.0, SplittingCoefficients.VERLET.oneNorm(), EXACT);
        assertEquals(4.4048, SplittingCoefficients.YOSHIDA_4.oneNorm(), 1.0e-4);
        assertEquals(2.3159, SplittingCoefficients.SUZUKI_4.oneNorm(), 1.0e-4);
        assertEquals(1.1562, SplittingCoefficients.BLANES_MOAN_4.oneNorm(), 1.0e-4);
        assertEquals(1.2061, SplittingCoefficients.BLANES_MOAN_6.oneNorm(), 1.0e-4);

        assertTrue("Yoshida's middle step runs backwards by more than one whole step",
                SplittingCoefficients.YOSHIDA_4.b()[1] < -1.7);
        assertTrue("Suzuki's runs backwards by two thirds of one",
                SplittingCoefficients.SUZUKI_4.b()[2] < -0.65
                        && SplittingCoefficients.SUZUKI_4.b()[2] > -0.66);

        // the ranking of the one-norms is the ranking of the error constants
        assertTrue(SplittingCoefficients.YOSHIDA_4.oneNorm() > SplittingCoefficients.SUZUKI_4.oneNorm());
        assertTrue(SplittingCoefficients.SUZUKI_4.oneNorm() > SplittingCoefficients.BLANES_MOAN_4.oneNorm());
        assertTrue(errorConstant(SplittingCoefficients.YOSHIDA_4) > errorConstant(
                SplittingCoefficients.SUZUKI_4));
        assertTrue(errorConstant(SplittingCoefficients.SUZUKI_4) > errorConstant(
                SplittingCoefficients.BLANES_MOAN_4));
    }

    @Test
    public void testAVelocityFirstMethodIsBuiltFromHalfItsCoefficients() {
        // as many kicks as drifts: the middle drift stands alone
        SplittingCoefficients even = SplittingCoefficients.velocityFirst("even", 2,
                new double[] { 0.25, 0.5 }, new double[] { 0.25, 0.25 });
        assertEquals(4, even.stages());
        assertTrue(even.isVelocityFirst());
        assertTrue(even.isPalindromic());

        // one kick more than drifts: the middle kick stands alone
        SplittingCoefficients odd = SplittingCoefficients.velocityFirst("odd", 2,
                new double[] { 0.5 }, new double[] { 0.25, 0.5 });
        assertEquals(3, odd.stages());
        assertTrue(odd.isVelocityFirst());
        assertTrue(odd.isPalindromic());

        try {
            SplittingCoefficients.velocityFirst("wrong", 2, new double[] { 0.5 },
                    new double[] { 0.2, 0.3, 0.3, 0.2 });
            fail("expected a refusal about the two lengths");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("as many kicks as drifts"));
        }
        try {
            SplittingCoefficients.velocityFirst("wrong", 2, null, new double[] { 0.5 });
            fail("expected a refusal about the halves");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("at least one drift"));
        }
    }

    @Test
    public void testACompositionMergesTheDriftsBetweenItsSteps() {
        // three Verlet steps of a third each is Verlet at a third of the step,
        // three times over: four drifts rather than six, and three kicks
        SplittingCoefficients thrice = SplittingCoefficients.composition("thirds", 2,
                new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 });
        assertEquals(4, thrice.stages());
        assertEquals(3, thrice.evaluations());
        assertEquals(1.0 / 6.0, thrice.a()[0], EXACT);
        assertEquals(1.0 / 3.0, thrice.a()[1], EXACT);
        assertEquals(1.0 / 6.0, thrice.a()[3], EXACT);
        assertTrue(thrice.isPalindromic());
    }

    @Test
    public void testTheAccessorsHandOutCopies() {
        SplittingCoefficients c = SplittingCoefficients.YOSHIDA_4;
        double[] a = c.a();
        double kept = a[0];
        a[0] = 42.0;
        assertEquals(kept, c.a()[0], 0.0);
        double[] b = c.b();
        double keptB = b[0];
        b[0] = 42.0;
        assertEquals(keptB, c.b()[0], 0.0);
    }

    @Test
    public void testTheCoefficientsAreChecked() {
        double[] good = { 0.5, 0.5 };
        double[] kicks = { 1.0, 0.0 };
        refuse(null, 2, good, kicks, "name");
        refuse("x", 0, good, kicks, "order");
        refuse("x", 2, good, new double[] { 1.0 }, "must agree");
        refuse("x", 2, new double[] { 0.5, 0.4 }, kicks, "drift coefficients sum");
        refuse("x", 2, good, new double[] { 0.9, 0.0 }, "kick coefficients sum");
        refuse("x", 2, new double[] { Double.NaN, 0.5 }, kicks, "not finite");
        try {
            SplittingCoefficients.composition("x", 2, new double[] { 0.5, 0.4 });
            fail("expected a refusal about the sum");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("sum"));
        }
        try {
            SplittingCoefficients.composition("x", 2, new double[0]);
            fail("expected a refusal about the fractions");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("fraction"));
        }
    }

    @Test
    public void testWhatTheCoefficientsSayAboutThemselves() {
        SplittingCoefficients v = SplittingCoefficients.VERLET;
        assertEquals(2, v.order());
        assertEquals(2, v.stages());
        assertEquals("Stoermer-Verlet", v.name());
        assertTrue(v.toString().contains("order 2"));
        assertEquals(v, new SplittingCoefficients("under another name", 2, v.a(), v.b()));
        assertEquals(v.hashCode(), new SplittingCoefficients("x", 2, v.a(), v.b()).hashCode());
        assertFalse(v.equals(SplittingCoefficients.YOSHIDA_4));
    }

    // ------------------------------------------------------------------

    /**
     * One step of the method applied to {@code q'' = -q}, as the exact
     * {@code 2 x 2} matrix {@code {m00, m01, m10, m11}}.
     */
    private static double[] map(SplittingCoefficients c, double h) {
        double[] a = c.a();
        double[] b = c.b();
        double m00 = 1.0;
        double m01 = 0.0;
        double m10 = 0.0;
        double m11 = 1.0;
        for (int i = 0; i < a.length; ++i) {
            double drift = a[i] * h;
            m00 += drift * m10;
            m01 += drift * m11;
            double kick = -b[i] * h;
            m10 += kick * m00;
            m11 += kick * m01;
        }
        return new double[] { m00, m01, m10, m11 };
    }

    /** The worst entry of the map against the exact flow at one step size. */
    private static double mapError(SplittingCoefficients c, double h) {
        double[] m = map(c, h);
        double[] exact = { Math.cos(h), Math.sin(h), -Math.sin(h), Math.cos(h) };
        double worst = 0.0;
        for (int j = 0; j < 4; ++j) {
            worst = Math.max(worst, Math.abs(m[j] - exact[j]));
        }
        return worst;
    }

    /** Halved step sizes from 0.8, stopping where the rounding takes over. */
    private static double[] usableMapErrors(SplittingCoefficients c) {
        double[] buffer = new double[24];
        int count = 0;
        double h = 0.8;
        while (count < buffer.length) {
            double error = mapError(c, h);
            if (error < 1.0e-14) {
                break;
            }
            buffer[count++] = error;
            h *= 0.5;
        }
        return Arrays.copyOf(buffer, count);
    }

    /** The coefficient of the leading error term, at a settled step size. */
    private static double errorConstant(SplittingCoefficients c) {
        double h = settled(c);
        double[] m = map(c, h);
        double[] exact = { Math.cos(h), Math.sin(h), -Math.sin(h), Math.cos(h) };
        double worst = 0.0;
        for (int j = 0; j < 4; ++j) {
            worst = Math.max(worst, Math.abs(m[j] - exact[j]));
        }
        return worst / Math.pow(h, c.order() + 1);
    }

    /**
     * What a method costs to reach a given error, up to a factor common to all
     * of them: it has to take steps of {@code (eps / C)^(1/p)} and each costs
     * its evaluations.
     * <p>
     * <b>Only comparable across methods of the same order at a fixed
     * accuracy</b>, because the accuracy enters through the order. Two methods
     * of the same order have a cost ratio that does not depend on it; a higher
     * order method gains as the accuracy demanded tightens.
     */
    private static double cost(SplittingCoefficients c, double accuracy) {
        return c.evaluations() * Math.pow(errorConstant(c) / accuracy, 1.0 / c.order());
    }

    private static double sum(double[] x) {
        double total = 0.0;
        for (int i = 0; i < x.length; ++i) {
            total += x[i];
        }
        return total;
    }

    private static double sumOfCubes(double[] x) {
        double total = 0.0;
        for (int i = 0; i < x.length; ++i) {
            total += x[i] * x[i] * x[i];
        }
        return total;
    }

    private static void refuse(String name, int order, double[] a, double[] b, String hint) {
        try {
            new SplittingCoefficients(name, order, a, b);
            fail("expected a refusal mentioning " + hint);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains(hint));
        }
    }
}
