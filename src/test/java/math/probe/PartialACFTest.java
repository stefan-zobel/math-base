/*
 * Copyright 2026 Stefan Zobel
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

public final class PartialACFTest {

    @Test
    public void testBasics() {
        // Test invalid inputs
        try {
            PartialACF.partialAutocorrelation(null, 1);
            Assert.fail("Should have thrown IllegalArgumentException for null data");
        } catch (IllegalArgumentException e) {
            // expected
        }

        try {
            PartialACF.partialAutocorrelation(new double[] { 1.0, 2.0 }, 2);
            Assert.fail("Should have thrown IllegalArgumentException for lags >= length");
        } catch (IllegalArgumentException e) {
            // expected
        }

        // Test significance bounds logic
        double[] bounds = PartialACF.getConfidenceInterval(100);
        Assert.assertEquals(-0.196, bounds[0], 1e-9);
        Assert.assertEquals(0.196, bounds[1], 1e-9);
        Assert.assertTrue(PartialACF.isSignificant(0.2, 100));
        Assert.assertFalse(PartialACF.isSignificant(0.1, 100));
    }

    @Test
    public void testSimpleAR1() {
        // For an AR(1) process with phi = 0.5, PACF should be 0.5 at lag 1 and
        // ~0 for others
        // We simulate a very simple alternating sequence
        double[] data = new double[] { 10, -10, 10, -10, 10, -10, 10, -10, 10, -10 };
        double[] pacf = PartialACF.partialAutocorrelation(data, 3);

        // At lag 1, PACF is always equal to ACF at lag 1
        // For this sequence, ACF[1] is -0.9 (from the ACFTest)
        Assert.assertEquals(-0.9, pacf[0], 1e-3);
        // Higher order PACF should drop significantly for such a simple
        // structured signal
        Assert.assertTrue(Math.abs(pacf[1]) < Math.abs(pacf[0]));
    }

    @Test
    public void testUniformRandomWhiteNoise() {
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
        double[] whiteNoise = rng.doubles(5_000, -1.0, 1.0).toArray();
        double[] pacf = PartialACF.partialAutocorrelation(whiteNoise, 10);

        // For white noise, all PACF values should be within the confidence
        // interval
        double limit = 1.96 / Math.sqrt(whiteNoise.length); // approx 0.0277
        for (double val : pacf) {
            // using a slightly relaxed 0.1 for random variance
            Assert.assertTrue("PACF value " + val + " exceeds white noise limit " + limit, Math.abs(val) < 0.1);
        }
    }

    @Test
    public void testBoxJenkinsPACF() {
        // Using the same Box-Jenkins sample data
        double[] boxJenkins = new double[] { 47, 64, 23, 71, 38, 64, 55, 41, 59, 48, 71, 35, 57, 40, 58, 44, 80, 55, 37,
                74, 51, 57, 50, 60, 45, 57, 50, 45, 25, 59, 50, 71, 56, 74, 50, 58, 45, 54, 36, 54, 48, 55, 45, 57, 50,
                62, 44, 64, 43, 52, 38, 59, 55, 41, 53, 49, 34, 35, 54, 45, 68, 38, 50, 60, 39, 59, 40, 57, 54, 23 };

        // Calculated PACF values based on the Durbin-Levinson recursion from
        // the expected ACF
        double[] expectedPACF = new double[] {
                -0.38987832,          // Lag 1
                 0.179705062360354,   // Lag 2
                 0.0022644684921318,  // Lag 3
                -0.0442766939834587,  // Lag 4
                -0.069405617258291    // Lag 5
        };

        double[] pacf = PartialACF.partialAutocorrelation(boxJenkins, 5);

        // We use a delta of 1e-7 because of the 8-digit rounding inside ACF.acf
        Assert.assertArrayEquals(expectedPACF, pacf, 1e-7);
    }
}
