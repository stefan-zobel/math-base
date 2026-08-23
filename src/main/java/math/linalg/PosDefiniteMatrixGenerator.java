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
package math.linalg;

import math.MathConsts;
import math.rng.PseudoRandom;
import math.rng.SplitMix64Seed;
import math.rng.Stc64;

/**
 * Random generator for positive definite matrices for testing purposes.
 */
public final class PosDefiniteMatrixGenerator {

    /**
     * Generate a random positive definite matrix of dimension
     * {@code dim x dim}, seeded from the package seed source, so the matrix
     * cannot be reproduced. Use {@link #generate(int, long)} when it has to be.
     *
     * @param dim
     *            dimension of the desired matrix
     * @return a random positive definite matrix
     */
    public static DMatrix generate(int dim) {
        return generate(dim, SplitMix64Seed.seed());
    }

    /**
     * Generate a random positive definite matrix of dimension
     * {@code dim x dim} from the given seed. The same seed and dimension
     * reproduce the same matrix exactly.
     *
     * @param dim
     *            dimension of the desired matrix
     * @param seed
     *            the seed the matrix is drawn from
     * @return a random positive definite matrix
     * @since 1.5.2
     */
    public static DMatrix generate(int dim, long seed) {
        // a generator per call rather than one shared between them: two
        // threads calling this at the same time would otherwise race on the
        // same state, and the gaussian cache in AbstractRng64 answers NaN when
        // it does
        PseudoRandom rng = new Stc64(seed);
        DMatrix Q = new DMatrix(dim, dim);
        for (int c = 0; c < Q.numColumns(); ++c) {
            for (int r = 0; r < Q.numRows(); ++r) {
                Q.setUnsafe(r, c, rng.nextGaussian());
            }
        }
        double d = rng.nextDouble(0.0, MathConsts.TWO_PI);
        DMatrix D = DMatrix.diag(dim, d);
        for (int i = 0; i < dim; ++i) {
            double v = D.get(i, i);
            v += rng.nextGaussian();
            D.set(i, i, v);
        }
        D = D.absInplace();
        return Q.transpose().mul(D).mul(Q);
    }

    private PosDefiniteMatrixGenerator() {
        throw new AssertionError();
    }
}
