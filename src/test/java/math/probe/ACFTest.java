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

    @Test
    public void testBasics() {
        Assert.assertNull(ACF.acf(null));
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
        double[] expected = new double[] { 1.0, 0.7, 0.41212121, 0.14848485, -0.07878788, -0.25757576, -0.37575758,
                -0.42121212, -0.38181818, -0.24545455 };
        Assert.assertArrayEquals(expected, ACF.acf(data), 1e-9);
    }

    @Test
    public void testUniformRandom() {
        PseudoRandom rng = DefaultRng.getGlobalPseudoRandom();
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

        double[] expected = new double[] { 1.0, -0.38987832, 0.30439408, -0.16555472, 0.07071932, -0.09703929,
                -0.04705769, 0.03537311, -0.0434582, -0.00479616, 0.01439314, 0.1099172, -0.06877849, 0.14803449,
                0.03576858, -0.00667781 };

        double[] acf = ACF.acf(boxJenkins, 15);

        Assert.assertArrayEquals(expected, acf, 1e-9);
    }
}
