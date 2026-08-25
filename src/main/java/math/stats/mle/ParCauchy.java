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
import math.distribution.Cauchy;
import math.stats.Validity;

/**
 * MLE for the parameters of the {@link Cauchy} distribution: the location
 * {@code x0} and the scale &gamma;.
 *
 * @since 1.5.3
 */
public final class ParCauchy implements Validity {

    /** location ({@code x0}) */
    public final double location;
    /** scale (&gamma;) */
    public final double scale;
    /** whether the fit settled */
    public final boolean converged;

    ParCauchy(double location, double scale, boolean converged) {
        this.location = location;
        this.scale = scale;
        this.converged = converged;
    }

    /**
     * The estimate is valid if the fit settled, the location is a number and
     * the scale is strictly positive.
     * <p>
     * {@link #converged} is worth reading rather than assuming. The iteration
     * can crawl for very small samples -- measured at {@code n = 4}, one in
     * twelve samples needs more than the iteration budget, and where it does,
     * the point it stopped at is far from the maximum.
     */
    @Override
    public boolean isValid() {
        return converged && !Arithmetic.isBadNum(location) && scale > 0.0 && !Arithmetic.isBadNum(scale);
    }
}
