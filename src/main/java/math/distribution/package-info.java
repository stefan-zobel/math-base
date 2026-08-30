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
 * {@link math.distribution.Dirichlet} and {@link math.distribution.Multinomial}
 * stand outside both interfaces, and deliberately. An outcome of either is a
 * vector, so there is no distribution function in closed form and no quantile
 * at all, and a mean is a vector and a variance a matrix; of what those
 * interfaces declare only the density or the mass survives the move to several
 * dimensions. Both are plain classes with a marginal that brings the whole one
 * dimensional apparatus back one component at a time -- a
 * {@link math.distribution.Beta} for the Dirichlet, a
 * {@link math.distribution.Binomial} for the multinomial.
 * <p>
 * The two are the halves of one conjugate pair. {@code Multinomial.logPmf} is
 * the likelihood of a count vector for known proportions;
 * {@code Dirichlet.posterior} is what observing those counts does to a belief
 * about the proportions, and {@code Dirichlet.logMarginalLikelihood} is the
 * same likelihood with the proportions integrated out.
 * <p>
 * {@link math.distribution.Categorical} is the multinomial of a single draw,
 * and it does <em>not</em> stand outside: an outcome of one draw is one
 * category, so it is an ordinary
 * {@link math.distribution.DiscreteDistribution} with a distribution function,
 * a quantile and both moments. It is the observation a
 * {@code Dirichlet.posterior} is updated with when observations arrive one at
 * a time rather than already counted, and it is the one law in this package
 * whose masses are given rather than computed from a formula -- its shape is
 * an argument.
 */
package math.distribution;
