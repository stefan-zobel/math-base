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

import math.rng.PseudoRandom;
import math.rng.Stc64;

/**
 * Random generator for positive definite matrices for testing purposes.
 */
public final class PosDefiniteMatrixGenerator {

    private static final PseudoRandom rng = Stc64.getDefault();

    /**
     * Generate a random positive definite matrix of dimension
     * {@code dim x dim}.
     * 
     * @param dim
     *            dimension of the desired matrix
     * @return a random positive definite matrix
     */
    public static DMatrix generate(int dim) {
        DMatrix Q = new DMatrix(dim, dim);
        for (int c = 0; c < Q.numColumns(); ++c) {
            for (int r = 0; r < Q.numRows(); ++r) {
                Q.setUnsafe(r, c, rng.nextGaussian());
            }
        }
        double d = rng.nextDouble(0.0, 10.0);
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
