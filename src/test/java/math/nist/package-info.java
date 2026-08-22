/**
 * Certification of this library against the NIST Statistical Reference
 * Datasets.
 * <p>
 * Every other test in this repository compares the library against itself: an
 * invariant, a second algorithm, a closed form. Both sides of such a check can
 * be wrong together. The values used here were computed elsewhere, in extended
 * precision, by people whose purpose was to give numerical software something
 * to be measured against, so this package is the one place that answers the
 * question outright -- <i>how many correct digits does this library
 * produce?</i>
 * <p>
 * {@link math.nist.Certificate} prints the answer as a table.
 * {@link math.nist.StRD} holds the datasets and their certified values,
 * {@link math.nist.Models} the five nonlinear models, and
 * {@link math.nist.Digits} the measure NIST itself scores packages with: the
 * log relative error {@code -log10(|got - certified| / |certified|)}, capped
 * at the fifteen digits NIST publishes.
 *
 * <p><b>What a digit count claims.</b> Not that the library computes
 * correctly, but how much of the answer double precision holds for a given
 * problem, and whether the code reaches it. For a least squares problem the
 * condition number of the design sets the ceiling: roughly
 * {@code log10(cond * eps)} digits are gone before any code is written. The
 * useful reading of the table is therefore never "6.5 is worse than 12.7", it
 * is "does the number reached match what this problem allows".
 *
 * <p><b>Wampler is where that can be checked</b>, because its five sets share
 * one design ({@code cond = 6.4e6}) and one exact solution, all six
 * coefficients equal to 1, and differ only in the noise added to the response.
 * With {@code cond * eps = 1.4e-9} the ceiling sits near nine digits, and the
 * worse of the two noise-free sets lands at 9.3: the solver takes what is
 * there. From there the residual standard deviation grows by a factor of a
 * hundred per step -- 2360, then 236015, then 23601450 -- and the digits fall
 * 9.3, 8.6, 6.5. That is the shape perturbation theory predicts for least
 * squares: with a vanishing residual the error is governed by the condition
 * number, with a growing one a term in the <i>square</i> of it enters.
 * Wampler5's 6.5 digits are a property of those twenty-one rows of data, not a
 * weakness of the solver. A package claiming twelve digits there would be the
 * suspicious one.
 *
 * <p><b>MGH09 is the same lesson, stated by the data.</b> Its certified
 * standard deviations, taken relative to the parameters they belong to, are
 * 5.9%, 102.6%, 65.7% and 66.2%: eleven observations for four parameters,
 * three of which the data barely determine and one of which is not pinned down
 * even in sign. Four digits of agreement on a quantity the experiment locates
 * only to within its own magnitude is already more precision than the question
 * contains. The residual sum of squares meanwhile agrees to 8.7 digits, and it
 * is the better known quantity on every nonlinear set: near a minimum the
 * objective is flat, so its value is determined long before its position is.
 *
 * <p><b>Filip is the one finding, and it is not about arithmetic.</b> At
 * {@code cond = 1.8e15} the standard bound guarantees no correct digits at
 * all, and {@link math.linalg.SvdLeastSquares} returns 7.5 of them on all
 * eleven coefficients. What fails is a default: the rank criterion in
 * {@link math.linalg.OLS} cannot tell ill conditioned from rank deficient, so
 * it declines a problem that has a certified answer. That is a question about
 * an API, not about a computation. Across the whole ladder this certification
 * found no defect in the arithmetic.
 *
 * <p><b>Reproducibility is part of what is certified.</b> The nonlinear models
 * go through {@code StrictMath} rather than {@code Math}, because the latter
 * is allowed one unit in the last place and its intrinsic differs between JDK
 * versions; on MGH10 that single bit was worth two digits of the fitted
 * parameters. A certificate that cannot print the same page twice is
 * certifying the runtime rather than the library.
 *
 * <p><b>What this does certify.</b> Two paths are covered, to the last
 * digit: the linear one through {@code OLS}, {@code LeastSquaresFit},
 * {@code FlatParallelJacobiSVD}, {@code SvdLeastSquares} and
 * {@code LSSummary}, and the nonlinear one through {@code LevenbergMarquardt},
 * {@code BoundedLevenbergMarquardt} and the four MINPACK classes beneath them.
 * <p>
 * Source: NIST/ITL Statistical Reference Datasets,
 * <a href="https://www.itl.nist.gov/div898/strd/">itl.nist.gov/div898/strd</a>,
 * retrieved 2026-08-22. Published by the United States government, in the
 * public domain, and released expressly for the testing of numerical software.
 */
package math.nist;
