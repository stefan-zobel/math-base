/*
 * Copyright 2013 Stefan Zobel
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
/**
 * Continuous and discrete distributions.
 * <p>
 * Every one of them answers in logarithms as well as in probabilities:
 * {@link math.distribution.ContinuousDistribution#logPdf(double)} and
 * {@link math.distribution.DiscreteDistribution#logPmf(int)} go on answering
 * where the density or the mass has underflowed to zero, which is where a
 * likelihood over more than a few hundred observations lives. Both are
 * {@code default} methods taking the logarithm of the density, and every
 * implementation here overrides that with a closed form.
 * <p>
 * {@link math.distribution.Dirichlet} stands outside both interfaces, and
 * deliberately. Over a simplex there is no distribution function in closed
 * form and no quantile at all, a mean is a vector and a variance a matrix,
 * so of what those interfaces declare only the density survives the move to
 * several dimensions. It is a plain class: a density, a marginal that is a
 * {@link math.distribution.Beta} and brings the whole one dimensional
 * apparatus back with it, and the conjugate update of a multinomial.
 */
package math.distribution;
