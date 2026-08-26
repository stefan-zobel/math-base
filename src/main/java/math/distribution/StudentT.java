/*
 * Copyright 2013, 2026 Stefan Zobel
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

import math.cern.FastGamma;
import math.cern.ProbabilityFuncs;

/**
 * StudentT distribution (a.k.a "T distribution").
 * <p>
 * <tt>p(x) = const  *  (1 + x^2/&nu;) ^ -(&nu;+1)/2</tt> where
 * <tt>const = &Gamma;((&nu;+1)/2) / (&radic;(&Pi;*&nu;) * &Gamma;(&nu;/2))</tt>
 * and <tt>&Gamma;(a)</tt> being the Gamma function and <tt>&nu;</tt> being the
 * degrees of freedom.
 * </p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Student%27s_t-distribution">Wikipedia
 * Student's t-distribution</a>.
 */
public class StudentT implements ContinuousDistribution {

    private final double df;
    private final double pdfConst;
    private final double logPdfConst;

    public StudentT(double df) {
        if (df <= 0.0) {
            throw new IllegalArgumentException("df <= 0.0 : " + df);
        }
        // the normalizing constant is log(gamma((df+1)/2) / gamma(df/2)), and
        // subtracting the two logarithms throws it away for a large df: they
        // are of order df*log(df) while their difference is of order log(df),
        // so from df = 1e13 onward the answer arrives quantized -- 20.0 where
        // 16.92 is right at 1e15 -- and from df = 1e16 it is zero, because
        // (df + 1) / 2 and df / 2 are then the same double
        this.pdfConst = Math.exp(FastGamma.logGammaRatio(df / 2.0, 0.5)) / Math.sqrt(Math.PI * df);
        // the logarithm of the line above, reached without the round trip
        // through exp and back, and without forming Math.PI * df, which
        // overflows for degrees of freedom above about 5.7e307
        this.logPdfConst = FastGamma.logGammaRatio(df / 2.0, 0.5) - 0.5 * (Math.log(Math.PI) + Math.log(df));
        this.df = df;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The power is taken through the logarithm because {@code 1 + x*x/df}
     * rounds to one long before the density does anything: at
     * {@code df = 1e14} the increment is 45 ulps of one and raising it to
     * {@code -(df+1)/2} multiplies that rounding by half the degrees of
     * freedom. {@code Math.log1p} keeps the increment at its own size.
     */
    @Override
    public double pdf(double x) {
        return pdfConst * Math.exp(-(df + 1.0) * 0.5 * Math.log1p(x * x / df));
    }

    /**
     * {@inheritDoc}
     * <p>
     * The expression {@link #pdf(double)} evaluates, with the normalizing
     * constant taken as a logarithm rather than as its exponential. The
     * density of a {@code t(3)} underflows to zero from about
     * {@code x = 1.05e81}, where this is still answering.
     */
    @Override
    public double logPdf(double x) {
        return logPdfConst - (df + 1.0) * 0.5 * Math.log1p(x * x / df);
    }

    @Override
    public double cdf(double x) {
        return ProbabilityFuncs.studentT(df, x);
    }

    @Override
    public double inverseCdf(double probability) {
        if (probability <= 0.0) {
            return supportLowerBound();
        }
        if (probability >= 1.0) {
            return supportUpperBound();
        }
        // Start at 0, not at mean(): the mean does not exist for df <= 1 and is
        // reported as NaN there, which poisoned every iteration below and made
        // this method return NaN for a quantile that is perfectly well defined
        // -- with one degree of freedom this distribution is the Cauchy, whose
        // quantile is tan(pi * (p - 1/2)). Zero is the median of every
        // symmetric t and therefore a valid starting point for any df.
        return findRoot(probability, 0.0, -Double.MAX_VALUE, Double.MAX_VALUE);
    }

    @Override
    public double mean() {
        if (df <= 1.0) {
            return Double.NaN;
        }
        return 0.0;
    }

    @Override
    public double variance() {
        if (df > 2.0) {
            return df / ((double) df - 2.0);
        }
        if (df == 2.0) {
            return Double.POSITIVE_INFINITY;
        }
        return Double.NaN;
    }

    public double getDegreesOfFreedom() {
        return df;
    }
}
