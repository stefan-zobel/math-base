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
            // Filip, which OLS refuses; see its own test below
            { 0.0, 0.0, 0.0, 0.0 } };

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

    /** Every set except Filip, against the digits it reached when this was written. */
    @Test
    public void theLadderKeepsItsDigits() {
        StRD.LinearSet[] sets = StRD.linear();
        for (int k = 0; k < sets.length; ++k) {
            StRD.LinearSet set = sets[k];
            if (set.name.equals("Filip")) {
                continue;
            }
            LSSummary fit = OLS.estimate(0.05, matrix(set), column(set.response()));

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
     * refuses it as rank deficient. The refusal is defensible -- the smallest
     * singular value is a hundred million times the largest one over the
     * reciprocal of the machine epsilon -- but a caller who knows the problem
     * has a certified answer needs a way through, and `SvdLeastSquares` is it.
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

        double tolerance = svd.sigma[0] * Math.max(filip.observations, filip.parameters) * 2.220446049250313e-16;
        assertTrue("the tolerance sits above the smallest singular value",
                tolerance > svd.sigma[filip.parameters - 1]);

        double[] truncated = SvdLeastSquares.solveTruncated(svd, filip.response(), tolerance);
        double[] full = SvdLeastSquares.solve(svd, filip.response(), 0.0);
        assertTrue("truncating loses the answer", Digits.worstOf(truncated, filip.certifiedBeta()) < 1.0);
        assertTrue("keeping everything finds it", Digits.worstOf(full, filip.certifiedBeta()) > 7.0);
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
