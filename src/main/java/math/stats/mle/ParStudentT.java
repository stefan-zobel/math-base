/*
 * Copyright 2018 Stefan Zobel
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
package math.stats.mle;

import math.cern.Arithmetic;
import math.distribution.StudentT;
import math.stats.Validity;

/**
 * MLE for the parameters of a {@link StudentT} distribution: the location
 * &mu;, the scale &sigma; and the degrees of freedom &nu;, which are
 * double-valued here.
 *
 * @since 1.5.3
 */
public final class ParStudentT implements Validity {

    /** location (&mu;) */
    public final double location;
    /** scale (&sigma;) */
    public final double scale;
    /** degrees of freedom (&nu;) */
    public final double df;
    /** whether the fit settled */
    public final boolean converged;

    ParStudentT(double location, double scale, double df, boolean converged) {
        this.location = location;
        this.scale = scale;
        this.df = df;
        this.converged = converged;
    }

    /**
     * The estimate is valid if the fit settled, the location is a number, the
     * scale is strictly positive, and &nu; is strictly positive <em>and
     * finite</em>.
     * <p>
     * {@link #df} of {@link Double#POSITIVE_INFINITY} is a result, not a
     * failure: it says the sample does not distinguish the t from its normal
     * limit, so the likelihood has no finite maximizer and {@link #location}
     * and {@link #scale} are the normal estimates. It is reported as invalid
     * because it is not a parameter any {@link StudentT} can be built from,
     * and a caller who wants to know which of the two happened has to look.
     */
    @Override
    public boolean isValid() {
        return converged && !Arithmetic.isBadNum(location) && scale > 0.0
                && !Arithmetic.isBadNum(scale) && df > 0.0 && df < Double.POSITIVE_INFINITY;
    }
}
