/**
 * Bayesian inference without a sampler.
 * <p>
 * This is the Bayesian counterpart of {@link math.stats.mle}, which fits the
 * same parameters by maximum likelihood. The difference is what comes back: a
 * point estimate there, a whole distribution here.
 * <p>
 * There is deliberately no Markov chain anywhere in this package. For one
 * parameter -- and, with the machinery in {@link math.solve}, for two and
 * three -- the normalizing constant is an ordinary integral, and an integral
 * is better behaved than a chain in the one way that matters: it can report
 * how far off it is, where a sampler can only report that it has not obviously
 * failed to converge. {@link math.stats.bayes.ScalarPosterior#errorEstimate()}
 * is that report.
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
