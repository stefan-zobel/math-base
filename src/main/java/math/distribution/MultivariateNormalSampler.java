/*
 * Copyright 2024 Stefan Zobel
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
package math.distribution;

import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;
import math.rng.Lcg64Xor1024Mix;
import math.rng.PseudoRandom;
import math.rng.SplitMix64Seed;

/**
 * A random sampler for arbitrary multivariate normal distributions.
 * <p>
 * Each sampler owns its generator, so a sampler must not be shared between
 * threads. Give each thread its own.
 */
public class MultivariateNormalSampler {

    private static final double EPS = 1.0e-7;

    private final PseudoRandom rng;
    private final DMatrix mean;
    private final DMatrix choleskyL;

    /**
     * Creates a sampler for a multivariate normal distribution identified by
     * {@code mean} and {@code covarianceMatrix}, seeded from the package seed
     * source, so its samples cannot be reproduced. Use
     * {@link #MultivariateNormalSampler(DMatrix, DMatrix, long)} when they have
     * to be.
     *
     * @param mean
     *            column vector containing the means
     * @param covarianceMatrix
     *            positive semidefinite covariance matrix
     */
    public MultivariateNormalSampler(DMatrix mean, DMatrix covarianceMatrix) {
        this(mean, covarianceMatrix, SplitMix64Seed.seed());
    }

    /**
     * Creates a sampler for a multivariate normal distribution identified by
     * {@code mean} and {@code covarianceMatrix}, seeded from {@code seed}. The
     * same seed, mean and covariance matrix reproduce the same samples exactly.
     *
     * @param mean
     *            column vector containing the means
     * @param covarianceMatrix
     *            positive semidefinite covariance matrix
     * @param seed
     *            the seed the samples are drawn from
     * @since 1.5.2
     */
    public MultivariateNormalSampler(DMatrix mean, DMatrix covarianceMatrix, long seed) {
        if (mean.numRows() != covarianceMatrix.numRows()) {
            throw new IllegalArgumentException("inconsistent matrix dimensions");
        }
        if (mean.numColumns() != 1) {
            throw new IllegalArgumentException("mean is not a column vector");
        }
        // add small perturbation to ensure numerical stability
        this.choleskyL = CholeskyDecomp
                .cholesky(covarianceMatrix.add(DMatrix.identity(covarianceMatrix.numRows()).scaleInplace(EPS)));
        this.mean = mean.copy();
        this.rng = new Lcg64Xor1024Mix(seed);
    }

    /**
     * Returns a {@code (d x numSamples)} matrix where {@code d} is the
     * dimension of the multivariate normal and {@code numSamples} is the
     * requested amount of random samples.
     * 
     * @param numSamples
     *            number of required samples
     * @return {@code (d x numSamples)} matrix containing the samples
     */
    public DMatrix sample(int numSamples) {
        DMatrix stdNormal = new DMatrix(choleskyL.numColumns(), numSamples);
        for (int c = 0; c < stdNormal.numColumns(); ++c) {
            for (int r = 0; r < stdNormal.numRows(); ++r) {
                stdNormal.setUnsafe(r, c, rng.nextGaussian());
            }
        }
        return choleskyL.mul(stdNormal).addBroadcastedVectorInplace(mean);
    }
}
