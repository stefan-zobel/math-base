/**
 * Bayesian inference without a sampler.
 * <p>
 * This is the Bayesian counterpart of {@link math.stats.mle}, which fits the
 * same parameters by maximum likelihood. The difference is what comes back: a
 * point estimate there, a whole distribution here.
 * <p>
 * There is deliberately no Markov chain anywhere in this package, and there
 * are three things in its place, in increasing order of how much they assume:
 * <ul>
 * <li>{@link math.stats.bayes.ScalarPosterior} and
 * {@link math.stats.bayes.BivariatePosterior} normalize one and two parameters
 * by quadrature and assume nothing about the shape of the posterior. An
 * integral is better behaved than a chain in the one way that matters: it can
 * report how far off it is, where a sampler can only report that it has not
 * obviously failed to converge. Their {@code errorEstimate()} is that
 * report.</li>
 * <li>{@link math.stats.bayes.LaplaceApproximation} fits a Gaussian at the mode
 * for any number of parameters, which is what three and above can afford --
 * measured, three parameters cost about 520 evaluations of the caller's
 * function against four to seven million for one quadrature integral. What it
 * assumes is that the posterior is roughly Gaussian, and <b>how wrong that is
 * has been measured</b> against the two classes above rather than asserted; its
 * class comment carries the figures.</li>
 * <li>{@link math.stats.bayes.BayesianLinearRegression} assumes the model
 * outright, and in exchange every part of the answer is closed form: no
 * integration, no approximation, and an evidence that picks the shrinkage.</li>
 * </ul>
 * <p>
 * The three are not alternatives so much as a ladder, and they check each
 * other: the quadrature classes are what licensed the approximation, and the
 * approximation is exact for the regression, which is how the regression is
 * checked a third time.
 * <p>
 * The posterior a caller supplies is a <b>log</b> density, which is what makes
 * the whole thing work at realistic sample sizes:
 * {@link math.distribution.ContinuousDistribution#logPdf(double)} and
 * {@link math.distribution.DiscreteDistribution#logPmf(int)} go on answering
 * long after the density itself has underflowed to zero, and a posterior over
 * a few thousand observations lives entirely in that region.
 * <p>
 * The conjugate families need none of this -- a Beta prior against a binomial
 * likelihood is a Beta again, by arithmetic -- and the classes in
 * {@link math.distribution} that carry those updates should be preferred where
 * they apply. What is here is for the posterior that has no closed form, and
 * the conjugate cases serve as the check on it.
 */
package math.stats.bayes;
