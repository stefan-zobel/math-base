package math.nist;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import math.linalg.DMatrix;
import math.linalg.FlatParallelJacobiSVD;
import math.linalg.LSSummary;
import math.linalg.OLS;
import math.linalg.SvdLeastSquares;

/**
 * The linear half of the certification: nine reference sets from Norris, which
 * anything can fit, to Filip, a degree ten polynomial whose design has a
 * condition number of about {@code 2e15}.
 * <p>
 * The floors asserted here were measured, not chosen. They say what this
 * library achieved when the certification was written, and they exist so that a
 * later change cannot quietly lose digits. They are not claims that the library
 * cannot do better.
 */
public class LinearLeastSquaresTest {

    /** Digits the fit has to reach: parameters, standard deviations, residual sd, R^2. */
    private static final double[][] FLOORS = {
            // Norris, Pontius, Longley
            { 12.0, 13.0, 13.0, 14.0 }, { 11.0, 12.0, 12.0, 14.0 }, { 10.0, 11.0, 11.0, 13.0 },
            // Wampler1 to Wampler5, the same design under growing noise
            { 9.0, 9.0, 9.0, 14.0 }, { 11.0, 13.0, 13.0, 14.0 }, { 9.0, 13.0, 13.0, 14.0 },
            { 8.0, 13.0, 14.0, 14.0 }, { 6.0, 13.0, 14.0, 14.0 },
            // Filip, which the default tolerance refuses; see its own test below
            { 7.0, 7.0, 8.0, 6.0 } };

    @Test
    public void theCollectionIsIntact() {
        StRD.LinearSet[] sets = StRD.linear();
        assertEquals(9, sets.length);
        assertEquals(FLOORS.length, sets.length);
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            assertEquals(set.name + " rows", set.observations, set.response().length);
            assertEquals(set.name + " design", set.observations * set.parameters, set.design().length);
            assertEquals(set.name + " certified parameters", set.parameters, set.certifiedBeta().length);
            assertEquals(set.name + " certified deviations", set.parameters,
                    set.certifiedStandardDeviation().length);
            assertTrue(set.name + " difficulty", set.difficulty.length() > 0);
            assertTrue(set.name + " R^2", set.rSquared > 0.0 && set.rSquared <= 1.0);
        }
        assertEquals("Norris", sets[0].name);
        assertEquals("Filip", sets[8].name);
        assertTrue("Filip is a polynomial fit", sets[8].isPolynomial());
        assertTrue("Longley is not", !sets[2].isPolynomial());
    }

    /**
     * All nine, against the digits they reached when this was written. Eight go
     * through the default entry point; Filip is the one design the default
     * tolerance declines, so it goes through the overload that lets the caller
     * set the tolerance, at {@code 0.0} -- accept anything that is not exactly
     * singular.
     */
    @Test
    public void theLadderKeepsItsDigits() {
        StRD.LinearSet[] sets = StRD.linear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            LSSummary fit = set.name.equals("Filip")
                    ? OLS.estimate(0.05, matrix(set), column(set.response()), 0.0)
                    : OLS.estimate(0.05, matrix(set), column(set.response()));

            double beta = Digits.worstOf(fit.getBeta().toArray(), set.certifiedBeta());
            double sd = Digits.worstOf(fit.getCoefficientStandardErrors().toArray(),
                    set.certifiedStandardDeviation());
            double residual = Digits.of(Math.sqrt(fit.getSigmaHatSquared()), set.residualStandardDeviation);
            double rSquared = Digits.of(fit.getRSquared(), set.rSquared);

            assertTrue(set.name + ": parameters agree to " + beta + " digits", beta >= FLOORS[k][0]);
            assertTrue(set.name + ": deviations agree to " + sd + " digits", sd >= FLOORS[k][1]);
            assertTrue(set.name + ": residual sd agrees to " + residual + " digits", residual >= FLOORS[k][2]);
            assertTrue(set.name + ": R^2 agrees to " + rSquared + " digits", rSquared >= FLOORS[k][3]);
        }
    }

    /**
     * Filip is the set that separates a package that factorizes from one that
     * forms normal equations, and this library does neither by default: `OLS`
     * refuses it as rank deficient. The refusal stays -- the smallest singular
     * value is below the level at which rounding alone could have produced it
     * -- but a caller who knows the problem has a certified answer now has two
     * ways through: `SvdLeastSquares` for the coefficients alone, and the
     * overload of `OLS.estimate` that takes the tolerance for the whole
     * summary. See {@link #theLadderKeepsItsDigits()} for the second.
     */
    @Test
    public void filipIsRefusedByTheFrontDoorAndSolvedByTheBack() {
        StRD.LinearSet filip = StRD.linear()[8];
        assertEquals("Filip", filip.name);

        try {
            OLS.estimate(0.05, matrix(filip), column(filip.response()));
            fail("OLS used to refuse Filip as rank deficient; if it no longer does, this test is the place to say so");
        } catch (RuntimeException refused) {
            assertTrue("the refusal should name the reason: " + refused.getMessage(),
                    refused.getMessage().contains("rank deficient"));
            assertTrue("and the conditioning it found: " + refused.getMessage(),
                    refused.getMessage().contains("cond(X)"));
            assertTrue("and where to go: " + refused.getMessage(),
                    refused.getMessage().contains("rankTolerance"));
        }

        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(filip.design(), filip.observations,
                filip.parameters);
        double condition = svd.sigma[0] / svd.sigma[svd.n - 1];
        assertTrue("the design is ill conditioned: " + condition, condition > 1.0e14);

        double[] beta = SvdLeastSquares.solve(svd, filip.response(), 0.0);
        double digits = Digits.worstOf(beta, filip.certifiedBeta());
        assertTrue("the singular value route reaches " + digits + " digits on all eleven coefficients",
                digits >= 7.0);
    }

    /**
     * Truncating at the customary rank tolerance throws away the direction that
     * carries Filip's answer: its smallest singular value is far above the
     * rounding level, so the design is ill conditioned rather than rank
     * deficient, and the truncated pseudo-inverse is the wrong tool for it.
     */
    @Test
    public void truncationDiscardsWhatFilipNeeds() {
        StRD.LinearSet filip = StRD.linear()[8];
        FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(filip.design(), filip.observations,
                filip.parameters);

        // the library's own tolerance rather than a second copy of the rule
        double tolerance = svd.sigma[0] * SvdLeastSquares.defaultRankTolerance(svd);
        assertTrue("the tolerance sits above the smallest singular value",
                tolerance > svd.sigma[filip.parameters - 1]);

        double[] truncated = SvdLeastSquares.solveTruncated(svd, filip.response(), tolerance);
        double[] full = SvdLeastSquares.solve(svd, filip.response(), 0.0);
        assertTrue("truncating loses the answer", Digits.worstOf(truncated, filip.certifiedBeta()) < 1.0);
        assertTrue("keeping everything finds it", Digits.worstOf(full, filip.certifiedBeta()) > 7.0);
    }

    /**
     * The overload is not a second policy: handed the tolerance the default
     * uses, it refuses exactly what the default refuses, because the two are
     * one code path.
     */
    @Test
    public void theOverloadAtTheDefaultToleranceRefusesFilipToo() {
        StRD.LinearSet filip = StRD.linear()[8];
        DMatrix X = matrix(filip);
        try {
            OLS.estimate(0.05, X, column(filip.response()), OLS.defaultRankTolerance(X));
            fail("the overload accepted at a tolerance the default rejects at");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("rank deficient"));
        }
    }

    /**
     * And on a design the default accepts, the two agree bit for bit -- the
     * assertion that adding the tolerance moved nothing.
     */
    @Test
    public void theTwoOverloadsAgreeBitForBitWhereBothAccept() {
        StRD.LinearSet[] sets = StRD.linear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            if (set.name.equals("Filip")) {
                continue;
            }
            DMatrix X = matrix(set);
            LSSummary a = OLS.estimate(0.05, X, column(set.response()));
            LSSummary b = OLS.estimate(0.05, X, column(set.response()), OLS.defaultRankTolerance(X));
            double[] betaA = a.getBeta().toArray();
            double[] betaB = b.getBeta().toArray();
            double[] seA = a.getCoefficientStandardErrors().toArray();
            double[] seB = b.getCoefficientStandardErrors().toArray();
            for (int j = 0; j < betaA.length; ++j) {
                assertEquals(set.name + ": beta " + j, betaA[j], betaB[j], 0.0);
                assertEquals(set.name + ": standard error " + j, seA[j], seB[j], 0.0);
            }
            assertEquals(set.name + ": residual variance", a.getSigmaHatSquared(), b.getSigmaHatSquared(), 0.0);
            assertEquals(set.name + ": R^2", a.getRSquared(), b.getRSquared(), 0.0);
            assertEquals(set.name + ": condition number", a.getConditionNumber(), b.getConditionNumber(), 0.0);
        }
    }

    /**
     * The summary carries the number that says how much of an ill conditioned
     * fit to believe, and it is the one the decomposition gives.
     */
    @Test
    public void theSummaryReportsTheConditioning() {
        StRD.LinearSet[] sets = StRD.linear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            LSSummary fit = set.name.equals("Filip")
                    ? OLS.estimate(0.05, matrix(set), column(set.response()), 0.0)
                    : OLS.estimate(0.05, matrix(set), column(set.response()));
            FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(set.design(), set.observations,
                    set.parameters);
            assertEquals(set.name + ": condition number", svd.sigma[0] / svd.sigma[set.parameters - 1],
                    fit.getConditionNumber(), 0.0);
        }
        // the two ends of the ladder, so that a plausible-looking wrong number
        // could not pass the equality above
        double norris = OLS.estimate(0.05, matrix(sets[0]), column(sets[0].response())).getConditionNumber();
        double filip = OLS.estimate(0.05, matrix(sets[8]), column(sets[8].response()), 0.0).getConditionNumber();
        assertTrue("Norris is well conditioned: " + norris, norris < 1.0e4);
        assertTrue("Filip is at the edge of what double arithmetic can do: " + filip, filip > 1.0e15);
    }

    /** {@code rankTolerance} is a fraction of the largest singular value, so it lives in {@code [0, 1)}. */
    @Test
    public void theRankToleranceIsValidated() {
        StRD.LinearSet norris = StRD.linear()[0];
        DMatrix X = matrix(norris);
        DMatrix y = column(norris.response());
        double[] rejected = { -1.0e-16, -1.0, 1.0, 2.0, Double.NaN, Double.POSITIVE_INFINITY };
        for (int i = 0; i < rejected.length; ++i) {
            try {
                OLS.estimate(0.05, X, y, rejected[i]);
                fail("accepted a rank tolerance of " + rejected[i]);
            } catch (IllegalArgumentException expected) {
                assertTrue(expected.getMessage(), expected.getMessage().contains("[0, 1)"));
            }
        }
    }

    /**
     * A tolerance of zero is the loosest setting there is and still refuses a
     * design that is exactly singular, because the ordinary least squares
     * filter would silently return a truncated fit for it rather than fail.
     */
    @Test
    public void anExactlySingularDesignIsRefusedEvenAtZero() {
        int n = 20;
        DMatrix X = new DMatrix(n, 3);
        DMatrix y = new DMatrix(n, 1);
        long lcg = 12345L;
        for (int i = 0; i < n; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double a = (lcg >>> 11) * 0x1.0p-53;
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            double b = (lcg >>> 11) * 0x1.0p-53;
            X.set(i, 0, 1.0);
            X.set(i, 1, a);
            // the third column is the second one over again
            X.set(i, 2, a);
            y.set(i, 0, b);
        }
        try {
            OLS.estimate(0.05, X, y, 0.0);
            fail("an exactly singular design was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("exactly singular"));
        }
    }

    /** The two routes into the same solver have to agree exactly, since one calls the other. */
    @Test
    public void bothRoutesGiveTheSameFit() {
        StRD.LinearSet[] sets = StRD.linear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            if (set.name.equals("Filip")) {
                continue;
            }
            double[] viaOls = OLS.estimate(0.05, matrix(set), column(set.response())).getBeta().toArray();
            FlatParallelJacobiSVD.Result svd = new FlatParallelJacobiSVD().decompose(set.design(), set.observations,
                    set.parameters);
            double[] viaSvd = SvdLeastSquares.solve(svd, set.response(), 0.0);
            for (int j = 0; j < viaOls.length; ++j) {
                assertEquals(set.name + " coefficient " + j, viaOls[j], viaSvd[j], 0.0);
            }
        }
    }

    /** Wampler1 to Wampler5 are one design under growing noise, and the digits fall with it. */
    @Test
    public void theWamplerLadderDegradesWithTheNoise() {
        StRD.LinearSet[] sets = StRD.linear();
        double[] digits = new double[5];
        for (int k = 0; k < 5; ++k) {
            StRD.LinearSet set = sets[3 + k];
            assertTrue(set.name.startsWith("Wampler"));
            digits[k] = Digits.worstOf(
                    OLS.estimate(0.05, matrix(set), column(set.response())).getBeta().toArray(),
                    set.certifiedBeta());
        }
        assertTrue("the noisiest of them is the least accurate: " + digits[4], digits[4] < digits[1]);
        assertTrue("and it still reaches six digits: " + digits[4], digits[4] > 6.0);
    }

    private static DMatrix matrix(StRD.LinearSet set) {
        double[] flat = set.design();
        DMatrix X = new DMatrix(set.observations, set.parameters);
        for (int c = 0; c < set.parameters; ++c) {
            for (int i = 0; i < set.observations; ++i) {
                X.set(i, c, flat[c * set.observations + i]);
            }
        }
        return X;
    }

    private static DMatrix column(double[] v) {
        DMatrix y = new DMatrix(v.length, 1);
        for (int i = 0; i < v.length; ++i) {
            y.set(i, 0, v[i]);
        }
        return y;
    }
}
