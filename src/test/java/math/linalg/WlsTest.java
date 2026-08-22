package math.linalg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import math.distribution.StudentT;
import math.list.DoubleList;

/**
 * Weighted least squares. The estimator has no closed form to compare against
 * that is not itself least squares, so the argument rests on invariants and on
 * two solvers that share no code with it: ordinary least squares on the data
 * set with rows repeated, and the active set method of
 * {@link BoundedLeastSquares} on the scaled problem.
 */
public class WlsTest {

    private static final double ALPHA = 0.05;

    /** Deterministic uniform noise in {@code [-0.5, 0.5]}. */
    private static final class Lcg {
        private long state;

        Lcg(long seed) {
            state = seed;
        }

        double next() {
            state = state * 6364136223846793005L + 1442695040888963407L;
            return ((state >>> 11) * 0x1.0p-53) - 0.5;
        }

        /** Roughly normal, by the sum of twelve uniforms. */
        double gauss() {
            double s = 0.0;
            for (int i = 0; i < 12; i++) {
                s += next();
            }
            return s;
        }
    }

    private static DMatrix design(int n, int p, Lcg rng) {
        DMatrix X = new DMatrix(n, p);
        for (int i = 0; i < n; i++) {
            X.set(i, 0, 1.0);
            for (int j = 1; j < p; j++) {
                X.set(i, j, rng.next());
            }
        }
        return X;
    }

    private static DMatrix response(DMatrix X, double[] beta, double noise, Lcg rng) {
        DMatrix y = new DMatrix(X.numRows(), 1);
        for (int i = 0; i < X.numRows(); i++) {
            double v = 0.0;
            for (int j = 0; j < X.numColumns(); j++) {
                v += beta[j] * X.get(i, j);
            }
            y.set(i, 0, v + noise * rng.next());
        }
        return y;
    }

    /** Weights spread log-uniformly over {@code [1, spread]}. */
    private static double[] weights(int n, double spread, Lcg rng) {
        double[] w = new double[n];
        for (int i = 0; i < n; i++) {
            w[i] = Math.exp(Math.log(spread) * (rng.next() + 0.5));
        }
        return w;
    }

    private static DMatrix scaleRows(DMatrix M, double[] w) {
        DMatrix out = new DMatrix(M.numRows(), M.numColumns());
        for (int i = 0; i < M.numRows(); i++) {
            double s = Math.sqrt(w[i]);
            for (int j = 0; j < M.numColumns(); j++) {
                out.set(i, j, s * M.get(i, j));
            }
        }
        return out;
    }

    private static long bits(double d) {
        return Double.doubleToRawLongBits(d);
    }

    private static void assertSameBits(String what, DoubleList a, DoubleList b) {
        assertEquals(what + " size", a.size(), b.size());
        for (int i = 0; i < a.size(); i++) {
            assertEquals(what + "[" + i + "] : " + a.get(i) + " vs " + b.get(i), bits(a.get(i)), bits(b.get(i)));
        }
    }

    @Test
    public void testUnitWeightsReproduceOrdinaryLeastSquaresBitForBit() {
        // Wls and OLS share one body; at w == 1 every scaling is skipped rather
        // than multiplied by one, so the two paths must not differ by a bit
        int[][] shapes = { { 200, 4 }, { 150, 3 }, { 256, 4 }, { 64, 8 }, { 17, 3 }, { 4, 3 }, { 33, 5 } };
        for (int s = 0; s < shapes.length; s++) {
            int n = shapes[s][0];
            int p = shapes[s][1];
            Lcg rng = new Lcg(2000L + s);
            double[] beta = new double[p];
            for (int j = 0; j < p; j++) {
                beta[j] = 1.0 + 0.5 * j;
            }
            DMatrix X = design(n, p, rng);
            DMatrix y = response(X, beta, 0.05, rng);
            double[] w = new double[n];
            Arrays.fill(w, 1.0);

            LSSummary o = OLS.estimate(ALPHA, X, y);
            LSSummary g = Wls.estimate(ALPHA, X, y, w);

            String what = "n=" + n + " p=" + p;
            assertSameBits(what + " beta", o.getBeta(), g.getBeta());
            assertSameBits(what + " yHat", o.getYHat(), g.getYHat());
            assertSameBits(what + " residuals", o.getResiduals(), g.getResiduals());
            assertSameBits(what + " standardErrors", o.getCoefficientStandardErrors(),
                    g.getCoefficientStandardErrors());
            assertSameBits(what + " tValues", o.getTValues(), g.getTValues());
            assertSameBits(what + " pValues", o.getPValues(), g.getPValues());
            assertEquals(what + " yBar", bits(o.getYBar()), bits(g.getYBar()));
            assertEquals(what + " rSquared", bits(o.getRSquared()), bits(g.getRSquared()));
            assertEquals(what + " sigmaHatSquared", bits(o.getSigmaHatSquared()), bits(g.getSigmaHatSquared()));
            assertEquals(what + " degreesOfFreedom", o.getDegreesOfFreedom(), g.getDegreesOfFreedom());
            DMatrix va = o.getVarianceCovarianceMatrix();
            DMatrix vb = g.getVarianceCovarianceMatrix();
            for (int i = 0; i < va.numRows(); i++) {
                for (int j = 0; j < va.numColumns(); j++) {
                    assertEquals(what + " varCov", bits(va.get(i, j)), bits(vb.get(i, j)));
                }
            }
            List<DoubleList> ca = o.getConfidenceIntervals();
            List<DoubleList> cb = g.getConfidenceIntervals();
            for (int i = 0; i < ca.size(); i++) {
                assertEquals(what + " ci lower", bits(ca.get(i).get(0)), bits(cb.get(i).get(0)));
                assertEquals(what + " ci upper", bits(ca.get(i).get(1)), bits(cb.get(i).get(1)));
            }
        }
    }

    @Test
    public void testExactDataIsRecoveredUnderAnyPositiveWeights() {
        // the residual is zero, so it does not matter what it is multiplied by
        Lcg rng = new Lcg(55L);
        double[] beta = { 2.0, -1.5, 0.75, 3.25 };
        DMatrix X = design(80, 4, rng);
        DMatrix y = response(X, beta, 0.0, rng);
        for (int k = 0; k < 5; k++) {
            double[] w = weights(80, Math.pow(10.0, 3 * k), rng);
            LSSummary s = Wls.estimate(ALPHA, X, y, w);
            for (int j = 0; j < beta.length; j++) {
                assertEquals("spread 1e" + (3 * k) + ", beta[" + j + "]", beta[j], s.getBeta().get(j), 1.0e-12);
            }
            assertEquals("R^2 on a perfect fit", 1.0, s.getRSquared(), 1.0e-12);
        }
    }

    @Test
    public void testTheWeightedNormalEquationsHoldAtTheAnswer() {
        // X' W (y - X beta) = 0 is what being the weighted least squares
        // solution means, and it is checkable without solving anything
        int[][] shapes = { { 120, 4 }, { 400, 9 }, { 40, 3 } };
        double[] spreads = { 1.0, 1.0e4, 1.0e8, 1.0e12 };
        for (int sh = 0; sh < shapes.length; sh++) {
            for (int k = 0; k < spreads.length; k++) {
                int n = shapes[sh][0];
                int p = shapes[sh][1];
                Lcg rng = new Lcg(31L + 7 * sh + k);
                double[] beta = new double[p];
                for (int j = 0; j < p; j++) {
                    beta[j] = 1.0 - 0.3 * j;
                }
                DMatrix X = design(n, p, rng);
                DMatrix y = response(X, beta, 0.1, rng);
                double[] w = weights(n, spreads[k], rng);
                LSSummary fit = Wls.estimate(ALPHA, X, y, w);
                for (int j = 0; j < p; j++) {
                    double g = 0.0;
                    double scale = 0.0;
                    for (int i = 0; i < n; i++) {
                        double t = X.get(i, j) * w[i] * fit.getResiduals().get(i);
                        g += t;
                        scale += Math.abs(t);
                    }
                    assertEquals("n=" + n + " spread " + spreads[k] + ", gradient " + j, 0.0, g / scale, 1.0e-11);
                }
            }
        }
    }

    @Test
    public void testIntegerWeightsAgreeWithTheDataSetThatRepeatsRows() {
        // a weight of 3 is three copies of the observation, which is an
        // ordinary least squares problem and shares no code with the weighting
        int[] sizes = { 30, 80 };
        for (int s = 0; s < sizes.length; s++) {
            int n = sizes[s];
            int p = 4;
            Lcg rng = new Lcg(4242L + s);
            DMatrix X = design(n, p, rng);
            DMatrix y = response(X, new double[] { 1.0, 2.0, -3.0, 0.5 }, 0.2, rng);
            int[] counts = new int[n];
            int total = 0;
            double[] w = new double[n];
            for (int i = 0; i < n; i++) {
                counts[i] = 1 + (int) (5.0 * (rng.next() + 0.5));
                w[i] = counts[i];
                total += counts[i];
            }
            DMatrix Xr = new DMatrix(total, p);
            DMatrix yr = new DMatrix(total, 1);
            int row = 0;
            for (int i = 0; i < n; i++) {
                for (int c = 0; c < counts[i]; c++) {
                    for (int j = 0; j < p; j++) {
                        Xr.set(row, j, X.get(i, j));
                    }
                    yr.set(row, 0, y.get(i, 0));
                    row++;
                }
            }

            LSSummary weighted = Wls.estimate(ALPHA, X, y, w);
            LSSummary repeated = OLS.estimate(ALPHA, Xr, yr);

            for (int j = 0; j < p; j++) {
                assertEquals("n=" + n + ", beta[" + j + "]", repeated.getBeta().get(j), weighted.getBeta().get(j),
                        1.0e-11);
            }
            assertEquals("n=" + n + ", yBar", repeated.getYBar(), weighted.getYBar(), 1.0e-11);
            assertEquals("n=" + n + ", rSquared", repeated.getRSquared(), weighted.getRSquared(), 1.0e-11);
            // the residual variance is the one thing that differs, and it
            // differs by exactly the ratio of the degrees of freedom
            double ratio = weighted.getSigmaHatSquared() / repeated.getSigmaHatSquared();
            double dfRatio = (double) repeated.getDegreesOfFreedom() / weighted.getDegreesOfFreedom();
            assertEquals("n=" + n + ", sigmaHat^2 ratio", dfRatio, ratio, 1.0e-12 * dfRatio);
            assertEquals("n=" + n + ", df", n - p, weighted.getDegreesOfFreedom());
            assertEquals("n=" + n + ", repeated df", total - p, repeated.getDegreesOfFreedom());
        }
    }

    @Test
    public void testAnIndependentSolverReachesTheSameCoefficients() {
        // BoundedLeastSquares with the bounds opened all the way is an
        // unconstrained least squares solve by an active set method
        int[][] shapes = { { 120, 4 }, { 60, 8 } };
        double[] spreads = { 1.0, 1.0e6 };
        for (int sh = 0; sh < shapes.length; sh++) {
            for (int k = 0; k < spreads.length; k++) {
                int n = shapes[sh][0];
                int p = shapes[sh][1];
                Lcg rng = new Lcg(808L + 3 * sh + k);
                double[] beta = new double[p];
                for (int j = 0; j < p; j++) {
                    beta[j] = 0.5 + j;
                }
                DMatrix X = design(n, p, rng);
                DMatrix y = response(X, beta, 0.1, rng);
                double[] w = weights(n, spreads[k], rng);
                LSSummary fit = Wls.estimate(ALPHA, X, y, w);

                double[] lo = new double[p];
                double[] hi = new double[p];
                Arrays.fill(lo, Double.NEGATIVE_INFINITY);
                Arrays.fill(hi, Double.POSITIVE_INFINITY);
                BoundedLeastSquares.Result r = BoundedLeastSquares.bounded(scaleRows(X, w), scaleRows(y, w), lo,
                        hi);

                assertTrue("the active set solver did not converge", r.converged);
                for (int j = 0; j < p; j++) {
                    assertEquals("n=" + n + " spread " + spreads[k] + ", beta[" + j + "]", r.solution[j],
                            fit.getBeta().get(j), 1.0e-11);
                }
            }
        }
    }

    @Test
    public void testRescalingEveryWeightMovesOnlyTheResidualVariance() {
        // w and c*w are the same problem; only sigmaHat^2 carries the factor
        Lcg rng = new Lcg(6161L);
        int n = 150;
        int p = 4;
        DMatrix X = design(n, p, rng);
        DMatrix y = response(X, new double[] { 2.0, -1.5, 0.75, 3.25 }, 0.1, rng);
        double[] w = weights(n, 1.0e4, rng);
        LSSummary base = Wls.estimate(ALPHA, X, y, w);

        double[] factors = { 2.0, 1024.0, 1.0e6, 3.0, 0.1, 1.0e-8 };
        for (int f = 0; f < factors.length; f++) {
            double c = factors[f];
            double[] wc = new double[n];
            for (int i = 0; i < n; i++) {
                wc[i] = c * w[i];
            }
            LSSummary s = Wls.estimate(ALPHA, X, y, wc);
            boolean powerOfTwo = (c == 1024.0);
            for (int j = 0; j < p; j++) {
                if (powerOfTwo) {
                    // a power of two rescales the singular values exactly
                    assertEquals("c = " + c + ", beta[" + j + "]", bits(base.getBeta().get(j)),
                            bits(s.getBeta().get(j)));
                } else {
                    assertEquals("c = " + c + ", beta[" + j + "]", base.getBeta().get(j), s.getBeta().get(j),
                            1.0e-13 * Math.abs(base.getBeta().get(j)));
                }
                assertEquals("c = " + c + ", standard error " + j, base.getCoefficientStandardErrors().get(j),
                        s.getCoefficientStandardErrors().get(j),
                        1.0e-13 * base.getCoefficientStandardErrors().get(j));
                assertEquals("c = " + c + ", p value " + j, base.getPValues().get(j), s.getPValues().get(j),
                        1.0e-10 * base.getPValues().get(j));
            }
            assertEquals("c = " + c + ", sigmaHat^2", c, s.getSigmaHatSquared() / base.getSigmaHatSquared(),
                    1.0e-12 * c);
            assertEquals("c = " + c + ", rSquared", base.getRSquared(), s.getRSquared(), 1.0e-12);
        }
    }

    @Test
    public void testEqualWeightsThatAreNotOneChangeNothingButTheVariance() {
        Lcg rng = new Lcg(999L);
        int n = 60;
        int p = 3;
        DMatrix X = design(n, p, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, -1.0 }, 0.3, rng);
        double[] w = new double[n];
        Arrays.fill(w, 7.5);

        LSSummary flat = Wls.estimate(ALPHA, X, y, w);
        LSSummary plain = OLS.estimate(ALPHA, X, y);

        for (int j = 0; j < p; j++) {
            assertEquals("beta[" + j + "]", plain.getBeta().get(j), flat.getBeta().get(j),
                    1.0e-12 * Math.abs(plain.getBeta().get(j)));
        }
        assertEquals("rSquared", plain.getRSquared(), flat.getRSquared(), 1.0e-12);
        assertEquals("sigmaHat^2", 7.5, flat.getSigmaHatSquared() / plain.getSigmaHatSquared(), 1.0e-12);
    }

    @Test
    public void testAnOverwhelmingWeightPullsTheFitOntoThatObservation() {
        Lcg rng = new Lcg(1234L);
        int n = 60;
        int p = 3;
        DMatrix X = design(n, p, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, -1.0 }, 0.3, rng);
        double[] w = new double[n];
        Arrays.fill(w, 1.0);
        w[13] = 1.0e14;

        LSSummary s = Wls.estimate(ALPHA, X, y, w);

        double pinned = Math.abs(s.getResiduals().get(13));
        double elsewhere = 0.0;
        for (int i = 0; i < n; i++) {
            if (i != 13) {
                elsewhere = Math.max(elsewhere, Math.abs(s.getResiduals().get(i)));
            }
        }
        assertTrue("the heavy observation was not fitted: " + pinned + " against " + elsewhere,
                pinned < elsewhere * 1.0e-10);
    }

    @Test
    public void testWeightingRecoversTheTruthBetterOnHeteroskedasticData() {
        // w = 1/sigma^2 is the efficient choice and ordinary least squares is
        // not; deterministic, so this is a fixed comparison and not a gamble
        int n = 200;
        int p = 4;
        double[] beta = { 2.0, -1.5, 0.75, 3.25 };
        int seeds = 60;
        int wins = 0;
        double sumOls = 0.0;
        double sumWls = 0.0;
        for (int s = 0; s < seeds; s++) {
            Lcg rng = new Lcg(20260821L + s);
            DMatrix X = design(n, p, rng);
            double[] w = new double[n];
            DMatrix y = new DMatrix(n, 1);
            for (int i = 0; i < n; i++) {
                double sigma = Math.pow(10.0, 1.5 * rng.next());
                w[i] = 1.0 / (sigma * sigma);
                double v = 0.0;
                for (int j = 0; j < p; j++) {
                    v += beta[j] * X.get(i, j);
                }
                y.set(i, 0, v + sigma * rng.gauss());
            }
            LSSummary o = OLS.estimate(ALPHA, X, y);
            LSSummary g = Wls.estimate(ALPHA, X, y, w);
            double eo = 0.0;
            double eg = 0.0;
            for (int j = 0; j < p; j++) {
                double a = o.getBeta().get(j) - beta[j];
                double b = g.getBeta().get(j) - beta[j];
                eo += a * a;
                eg += b * b;
            }
            sumOls += eo;
            sumWls += eg;
            if (eg < eo) {
                wins++;
            }
        }
        assertTrue("weighting won only " + wins + " of " + seeds, wins >= 55);
        assertTrue("aggregate error ratio was only " + (sumOls / sumWls), sumOls / sumWls > 10.0);
    }

    @Test
    public void testTheNormalEquationsAreNotUsableWhereThisStillIs() {
        // weights concentrated on a few rows are the case that drives the
        // condition number of the scaled design towards its bound, and they
        // are also the ordinary case: a handful of precise measurements
        Lcg rng = new Lcg(4711L);
        int n = 120;
        int p = 4;
        DMatrix X = design(n, p, rng);
        DMatrix y = response(X, new double[] { 2.0, -1.5, 0.75, 3.25 }, 0.05, rng);
        double[] w = new double[n];
        Arrays.fill(w, 1.0);
        for (int k = 0; k < p - 1; k++) {
            w[3 * k + 5] = 1.0e12;
        }

        LSSummary fit = Wls.estimate(ALPHA, X, y, w);
        for (int j = 0; j < p; j++) {
            double b = fit.getBeta().get(j);
            assertTrue("beta[" + j + "] is not finite : " + b, !Double.isNaN(b) && !Double.isInfinite(b));
        }
        for (int j = 0; j < p; j++) {
            double g = 0.0;
            double scale = 0.0;
            for (int i = 0; i < n; i++) {
                double t = X.get(i, j) * w[i] * fit.getResiduals().get(i);
                g += t;
                scale += Math.abs(t);
            }
            assertEquals("gradient " + j, 0.0, g / scale, 1.0e-2);
        }

        // the same fit through (X'WX)^-1 X'W y, which is what this class does
        // not do: it either refuses outright or is far off
        DMatrix Xw = scaleRows(X, w);
        DMatrix yw = scaleRows(y, w);
        boolean unusable;
        try {
            DMatrix ne = Xw.transpose().mul(Xw).inverse().mul(Xw.transpose()).mul(yw);
            double worst = 0.0;
            for (int j = 0; j < p; j++) {
                worst = Math.max(worst, Math.abs(ne.get(j, 0) - fit.getBeta().get(j)));
            }
            unusable = worst > 1.0e-3;
        } catch (RuntimeException expected) {
            unusable = true;
        }
        assertTrue("the normal equations coped after all, which makes this test pointless", unusable);
    }

    @Test
    public void testTheSummaryIsInternallyConsistent() {
        Lcg rng = new Lcg(31337L);
        int n = 90;
        int p = 4;
        DMatrix X = design(n, p, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, -1.0, 0.5 }, 0.2, rng);
        double[] w = weights(n, 1.0e5, rng);

        LSSummary s = Wls.estimate(ALPHA, X, y, w);

        assertNotNull("a weighted fit has weights", s.getWeights());
        assertNull("an unweighted fit has none", OLS.estimate(ALPHA, X, y).getWeights());
        assertEquals(n, s.getWeights().size());

        double sumW = 0.0;
        double sumWy = 0.0;
        double rss = 0.0;
        for (int i = 0; i < n; i++) {
            assertEquals("weights[" + i + "]", bits(w[i]), bits(s.getWeights().get(i)));
            // the residuals are the raw ones, not the weighted ones
            assertEquals("residual " + i, y.get(i, 0) - s.getYHat().get(i), s.getResiduals().get(i), 0.0);
            sumW += w[i];
            sumWy += w[i] * y.get(i, 0);
            rss += w[i] * s.getResiduals().get(i) * s.getResiduals().get(i);
        }
        assertEquals("degrees of freedom", n - p, s.getDegreesOfFreedom());
        assertEquals("sigmaHat^2 is the weighted residual variance", rss / (n - p), s.getSigmaHatSquared(),
                1.0e-11 * s.getSigmaHatSquared());

        double yBar = sumWy / sumW;
        assertEquals("yBar is the weighted mean", yBar, s.getYBar(), 1.0e-11 * Math.abs(yBar));
        double sqe = 0.0;
        double sqt = 0.0;
        for (int i = 0; i < n; i++) {
            double a = s.getYHat().get(i) - yBar;
            double b = y.get(i, 0) - yBar;
            sqe += w[i] * a * a;
            sqt += w[i] * b * b;
        }
        assertEquals("rSquared uses the weighted mean", sqe / sqt, s.getRSquared(), 1.0e-11);

        double t = new StudentT(s.getDegreesOfFreedom()).inverseCdf(1.0 - ALPHA / 2.0);
        for (int j = 0; j < p; j++) {
            double lo = s.getConfidenceIntervals().get(j).get(0);
            double hi = s.getConfidenceIntervals().get(j).get(1);
            double se = s.getCoefficientStandardErrors().get(j);
            assertEquals("interval " + j + " is centred on beta", s.getBeta().get(j), 0.5 * (lo + hi), 1.0e-9);
            assertEquals("interval " + j + " has half width t*se", t * se, 0.5 * (hi - lo), 1.0e-9 * t * se);
            assertEquals("t value " + j, s.getBeta().get(j) / se, s.getTValues().get(j), 0.0);
            assertEquals("standard error " + j, Math.sqrt(s.getVarianceCovarianceMatrix().get(j, j)), se, 0.0);
        }
    }

    @Test
    public void testClearTemporariesAlsoDropsTheWeights() {
        Lcg rng = new Lcg(17L);
        DMatrix X = design(40, 3, rng);
        DMatrix y = response(X, new double[] { 1.0, 1.0, 1.0 }, 0.1, rng);
        double[] w = weights(40, 100.0, rng);
        LSSummary s = Wls.estimate(ALPHA, X, y, w);
        assertNotNull(s.getWeights());
        s.clearTemporaries();
        assertNull("the weights are a per-observation temporary like the residuals", s.getWeights());
        assertNull(s.getResiduals());
        assertNotNull("the coefficients are not", s.getBeta());
    }

    @Test
    public void testMalformedWeightsAreRejected() {
        Lcg rng = new Lcg(5L);
        DMatrix X = design(20, 3, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, 3.0 }, 0.1, rng);
        double[] good = new double[20];
        Arrays.fill(good, 1.0);

        expectIae(X, y, null, "null weights");
        expectIae(X, y, new double[19], "too few weights");
        expectIae(X, y, new double[21], "too many weights");
        expectIae(X, y, at(good, 4, 0.0), "a zero weight");
        expectIae(X, y, at(good, 4, -1.0), "a negative weight");
        expectIae(X, y, at(good, 4, -0.0), "a negative zero weight");
        expectIae(X, y, at(good, 4, Double.NaN), "a NaN weight");
        expectIae(X, y, at(good, 4, Double.POSITIVE_INFINITY), "an infinite weight");
        try {
            Wls.estimate(ALPHA, null, y, good);
            fail("a null design was accepted");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }

    /**
     * The same code path as the four-argument form, and the tolerance applies
     * to the scaled design {@code sqrt(W) X} rather than to {@code X}, which is
     * the matrix the fit goes through and the one the weights make worse.
     */
    @Test
    public void testTheRankToleranceOverloadChangesNothingAtTheDefault() {
        Lcg rng = new Lcg(733L);
        DMatrix X = design(30, 3, rng);
        DMatrix y = response(X, new double[] { 1.0, -2.0, 0.5 }, 0.1, rng);
        double[] w = new double[30];
        for (int i = 0; i < w.length; i++) {
            w[i] = 1.0 + i;
        }

        LSSummary a = Wls.estimate(ALPHA, X, y, w);
        LSSummary b = Wls.estimate(ALPHA, X, y, w, OLS.defaultRankTolerance(X));
        for (int j = 0; j < X.numColumns(); j++) {
            assertEquals("beta " + j, a.getBeta().get(j), b.getBeta().get(j), 0.0);
            assertEquals("standard error " + j, a.getCoefficientStandardErrors().get(j),
                    b.getCoefficientStandardErrors().get(j), 0.0);
        }
        assertEquals(a.getRSquared(), b.getRSquared(), 0.0);
        assertEquals(a.getConditionNumber(), b.getConditionNumber(), 0.0);

        // and it is the conditioning of the scaled design that is reported:
        // weights this spread cannot leave it where the unweighted one was
        double unweighted = OLS.estimate(ALPHA, X, y).getConditionNumber();
        assertTrue("the weights left the conditioning untouched, which they cannot have",
                Math.abs(a.getConditionNumber() - unweighted) > 1.0e-6 * unweighted);

        assertTrue("the scaled design is not conditioned as this test assumes", a.getConditionNumber() > 2.0);
        try {
            Wls.estimate(ALPHA, X, y, w, 0.5);
            fail("a tolerance above the smallest relative singular value has to refuse");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("rank deficient"));
        }
    }

    @Test
    public void testTheInheritedChecksStillFire() {
        Lcg rng = new Lcg(6L);
        DMatrix X = design(20, 3, rng);
        DMatrix y = response(X, new double[] { 1.0, 2.0, 3.0 }, 0.1, rng);
        double[] w = new double[20];
        Arrays.fill(w, 1.0);
        double[] w2 = new double[2];
        Arrays.fill(w2, 1.0);

        expectIaeAlpha(X, y, w, 0.0, "alpha of 0");
        expectIaeAlpha(X, y, w, 1.0, "alpha of 1");
        expectIaeAlpha(X, y, w, -0.5, "a negative alpha");
        // fewer rows than columns
        DMatrix wide = design(2, 3, rng);
        DMatrix ys = response(wide, new double[] { 1.0, 2.0, 3.0 }, 0.1, rng);
        try {
            Wls.estimate(ALPHA, wide, ys, w2);
            fail("a design with fewer rows than columns was accepted");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
        // mismatched y
        try {
            Wls.estimate(ALPHA, X, response(design(19, 3, rng), new double[] { 1.0, 1.0, 1.0 }, 0.0, rng), w);
            fail("a mismatched regressand was accepted");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }

    private static double[] at(double[] w, int i, double v) {
        double[] copy = w.clone();
        copy[i] = v;
        return copy;
    }

    private static void expectIae(DMatrix X, DMatrix y, double[] w, String what) {
        try {
            Wls.estimate(ALPHA, X, y, w);
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            assertTrue(what + " threw without a message", expected.getMessage() != null);
        }
    }

    private static void expectIaeAlpha(DMatrix X, DMatrix y, double[] w, double alpha, String what) {
        try {
            Wls.estimate(alpha, X, y, w);
            fail(what + " was accepted");
        } catch (IllegalArgumentException expected) {
            // as intended
        }
    }
}
