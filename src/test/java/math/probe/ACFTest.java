/*
 * Copyright 2025 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.probe;

import org.junit.Assert;
import org.junit.Test;

import math.rng.DefaultRng;
import math.rng.PseudoRandom;

public final class ACFTest {

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
    }

    /**
     * The estimator as it is defined, which is what both code paths have to
     * reproduce: {@code sum (x_t - m)(x_{t+k} - m) / sum (x_t - m)^2}.
     */
    private static double[] byDefinition(double[] x, int lags) {
        int n = x.length;
        double m = 0.0;
        for (int i = 0; i < n; ++i) {
            m += x[i];
        }
        m /= n;
        double gamma0 = 0.0;
        for (int i = 0; i < n; ++i) {
            gamma0 += (x[i] - m) * (x[i] - m);
        }
        double[] r = new double[lags + 1];
        for (int k = 0; k <= lags; ++k) {
            double s = 0.0;
            for (int t = 0; t + k < n; ++t) {
                s += (x[t] - m) * (x[t + k] - m);
            }
            r[k] = s / gamma0;
        }
        return r;
    }

    private static double[] whiteNoise(int n, long seed) {
        Lcg rnd = new Lcg(seed);
        double[] x = new double[n];
        for (int i = 0; i < n; ++i) {
            x[i] = rnd.next();
        }
        return x;
    }

    private static double[] ar1(int n, double phi, long seed) {
        Lcg rnd = new Lcg(seed);
        double[] x = new double[n];
        double prev = 0.0;
        for (int i = 0; i < n; ++i) {
            prev = phi * prev + rnd.next();
            x[i] = prev;
        }
        return x;
    }

    private static double[] seasonal(int n) {
        double[] x = new double[n];
        for (int i = 0; i < n; ++i) {
            x[i] = Math.sin(2.0 * Math.PI * i / 12.0) + 0.01 * i;
        }
        return x;
    }

    /** A series of order 1e6 whose variation is of order 1. */
    private static double[] largeOffset(int n, long seed) {
        double[] x = whiteNoise(n, seed);
        for (int i = 0; i < n; ++i) {
            x[i] += 1.0e6;
        }
        return x;
    }

    private static double[][] series(int n) {
        return new double[][] { whiteNoise(n, 11L), ar1(n, 0.8, 12L), ar1(n, -0.6, 13L), seasonal(n),
                largeOffset(n, 14L) };
    }

    @Test
    public void testBasics() {
        try {
            ACF.acf(null);
            Assert.fail("null data should have been rejected");
        } catch (IllegalArgumentException expected) {
            // that is the contract
        }
        try {
            ACF.acf(null, 4);
            Assert.fail("null data should have been rejected");
        } catch (IllegalArgumentException expected) {
            // that is the contract
        }
        Assert.assertArrayEquals(new double[] {}, ACF.acf(new double[] {}), 1e-9);
        Assert.assertArrayEquals(new double[] { 1.0 }, ACF.acf(new double[] { 99.0 }), 1e-9);
        Assert.assertArrayEquals(new double[] { 1.0 }, ACF.acf(new double[] { 99.0 }, 1), 1e-9);
        Assert.assertArrayEquals(new double[] { 1.0 }, ACF.acf(new double[] { 99.0 }, 2), 1e-9);
        Assert.assertArrayEquals(new double[] { 1.0 }, ACF.acf(new double[] { 99.0 }, 3), 1e-9);
        Assert.assertArrayEquals(new double[] { 1.0 }, ACF.acf(new double[] { 1, 2, 3 }, 0), 1e-9);

        double[] data = new double[] { 5, -5, 5, -5, 5, -5, 5, -5, 5, -5 };
        Assert.assertArrayEquals(new double[] { 1.0, -0.9, 0.8, -0.7, 0.6, -0.5, 0.4, -0.3, 0.2, -0.1 }, ACF.acf(data),
                1e-3);
    }

    @Test
    public void testIncreasingSeq() {
        double[] data = new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        double[] expected = new double[] { 1.0, 0.70000000000000000, 0.41212121212121210, 0.14848484848484850,
                -0.078787878787878780, -0.25757575757575757, -0.37575757575757573, -0.42121212121212120,
                -0.38181818181818183, -0.24545454545454545 };
        Assert.assertArrayEquals(expected, ACF.acf(data), 1e-9);
    }

    @Test
    public void testUniformRandom() {
        // The seed is fixed. Unseeded this drew a fresh sample every run, so
        // the test was not reproducible -- it was not flaky, though, and the
        // sweep is what says so: over 5000 seeds the worst of the eleven lags
        // has a median of 0.034 and a maximum of 0.076 against the bound of
        // 0.1, and not one seed failed. This one sits at 0.350 of the bound.
        PseudoRandom rng = DefaultRng.newPseudoRandom(20260823005L);
        double[] uniform = rng.doubles(3000, -1000.0, 1000.0).toArray();
        double[] acf = ACF.acf(uniform);

        for (int i = 1; i <= 11; ++i) {
            Assert.assertTrue(Math.abs(acf[i]) < 0.1);
        }
    }

    @Test
    public void testBoxJenkins() {
        double[] boxJenkins = new double[] { 47, 64, 23, 71, 38, 64, 55, 41, 59, 48, 71, 35, 57, 40, 58, 44, 80, 55, 37,
                74, 51, 57, 50, 60, 45, 57, 50, 45, 25, 59, 50, 71, 56, 74, 50, 58, 45, 54, 36, 54, 48, 55, 45, 57, 50,
                62, 44, 64, 43, 52, 38, 59, 55, 41, 53, 49, 34, 35, 54, 45, 68, 38, 50, 60, 39, 59, 40, 57, 54, 23 };

        double[] expected = new double[] { 1.0, -0.38987831859560580, 0.30439408199642010, -0.16555471742498490,
                0.070719321101724810, -0.097039288117799480, -0.047057692464103170, 0.035373111886121224,
                -0.043458198567985350, -0.0047961622193493660, 0.014393137274734425, 0.10991720025367761,
                -0.068778491544323790, 0.14803448880650790, 0.035768581142729540, -0.0066778059016117090 };

        double[] acf = ACF.acf(boxJenkins, 15);

        Assert.assertArrayEquals(expected, acf, 1e-9);
    }

    @Test
    public void testAgreesWithTheDefinition() {
        // the class picks between a double sum and a transform by the number of
        // lags asked for; both have to reproduce the estimator as defined
        int[] lengths = { 16, 50, 128, 200, 512, 1000 };
        for (int l = 0; l < lengths.length; ++l) {
            int n = lengths[l];
            double[][] sets = series(n);
            for (int s = 0; s < sets.length; ++s) {
                int lags = n - 1;
                double[] expected = byDefinition(sets[s], lags);
                double[] actual = ACF.acf(sets[s]);
                Assert.assertEquals("n=" + n + ", set " + s, expected.length, actual.length);
                for (int k = 0; k <= lags; ++k) {
                    Assert.assertEquals("n=" + n + ", set " + s + ", lag " + k, expected[k], actual[k], 1e-12);
                }
            }
        }
    }

    @Test
    public void testTheTwoPathsAgreeWhereTheyMeet() {
        // one series, asked for few lags and for many, so both routes run on
        // it; the shared lags must not depend on which one did the work
        int[] lengths = { 200, 512, 1000, 2048 };
        for (int l = 0; l < lengths.length; ++l) {
            int n = lengths[l];
            double[][] sets = series(n);
            for (int s = 0; s < sets.length; ++s) {
                double[] few = ACF.acf(sets[s], 20);
                double[] many = ACF.acf(sets[s]);
                for (int k = 0; k < few.length; ++k) {
                    Assert.assertEquals("n=" + n + ", set " + s + ", lag " + k, many[k], few[k], 1e-12);
                }
            }
        }
    }

    @Test
    public void testLagZeroIsExactlyOne() {
        for (int n = 2; n <= 200; ++n) {
            double[][] sets = series(n);
            for (int s = 0; s < sets.length; ++s) {
                Assert.assertEquals("n=" + n + ", set " + s, 1.0, ACF.acf(sets[s], 3)[0], 0.0);
                Assert.assertEquals("n=" + n + ", set " + s, 1.0, ACF.acf(sets[s])[0], 0.0);
            }
        }
    }

    @Test
    public void testTooManyLagsAreClampedAndTheCountIsRight() {
        double[] x = whiteNoise(20, 41L);
        Assert.assertEquals(20, ACF.acf(x).length);
        Assert.assertEquals(20, ACF.acf(x, 100).length);
        Assert.assertEquals(20, ACF.acf(x, 19).length);
        Assert.assertEquals(20, ACF.acf(x, -5).length);
        Assert.assertEquals(6, ACF.acf(x, 5).length);
        Assert.assertEquals(1, ACF.acf(x, 0).length);
    }

    @Test
    public void testTheTransformPathSurvivesSmallAmplitudeData() {
        // the two paths have to agree however small the variation is. They did
        // not: the transform used to clean its output against an absolute
        // threshold, so from an amplitude of 1e-8 the spectrum was flushed to
        // zero and this returned 0.0, then NaN. The threshold in math.fft is
        // relative now, and this is the witness for it from outside that package
        int n = 512;
        double[] base = whiteNoise(n, 987L);
        for (int e = 0; e <= 20; e += 4) {
            double c = Math.pow(10.0, -e);
            double[] x = new double[n];
            for (int i = 0; i < n; ++i) {
                x[i] = c * base[i];
            }
            double[] allLags = ACF.acf(x);
            double[] fewLags = ACF.acf(x, 5);
            for (int k = 0; k <= 5; ++k) {
                Assert.assertFalse("lag " + k + " at amplitude 1e-" + e + " is NaN", Double.isNaN(allLags[k]));
                Assert.assertEquals("lag " + k + " at amplitude 1e-" + e, fewLags[k], allLags[k], 1e-12);
            }
        }
    }

    @Test
    public void testAConstantSeriesHasNoAutocorrelation() {
        // gamma_0 is zero, so every lag is 0/0; documented rather than hidden
        double[] constant = new double[] { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 };
        double[] r = ACF.acf(constant, 3);
        Assert.assertEquals(4, r.length);
        for (int k = 0; k < r.length; ++k) {
            Assert.assertTrue("lag " + k + " was " + r[k], Double.isNaN(r[k]));
        }
    }
}
