/**
 * Worked examples: six published datasets, each taken apart with the pieces
 * this library provides.
 * <p>
 * Every subpackage holds the same three files -- {@code Datasets} with the
 * numbers as literals, a {@code main} that answers one question about them in
 * numbered steps, and a JUnit test that guards what the {@code main} prints.
 * The data are measurements rather than simulations, so an answer can be wrong
 * in ways generated data never is, and each demo needs several packages at
 * once for a question no single class answers. That is where the seams show,
 * and finding them is what these are for.
 * <ul>
 * <li>{@link math.demo.maunaloa.MaunaLoaDemo} -- sixty-eight years of monthly
 * CO2, in eight steps across {@code probe}, {@code fit}, {@code linalg},
 * {@code solve} and {@code fft}: trend and seasonal cycle, residuals that turn
 * out not to be independent, and when the curve reaches 450 ppm.</li>
 * <li>{@link math.demo.wine.WineDemo} -- 178 wines and thirteen chemical
 * measurements classified by cultivar, through {@code linalg} (principal
 * components and lasso), {@code optim} (L-BFGS and OWL-QN) and {@code dl}.
 * Shows what happens when the columns are not standardized, and that the
 * optimizer maximizes.</li>
 * <li>{@link math.demo.michelson.MichelsonDemo} -- the 100 measurements
 * Michelson made of the speed of light in 1879, against a value that has been
 * exact by definition since 1983. The classical, percentile and BCa intervals
 * agree with each other to three decimals and all three miss the truth.</li>
 * <li>{@link math.demo.hubble.HubbleDemo} -- Hubble's 24 nebulae. With error
 * in both columns, fitting a line is four questions rather than one, and the
 * four answers disagree by more than any of their confidence intervals.</li>
 * <li>{@link math.demo.quakes.QuakeDemo} -- a year of 2115 earthquakes and the
 * Gutenberg-Richter exponent: a magnitude column that mixes scales, quantized
 * data that changes what a quantile means, and an extrapolation with no way of
 * knowing where its model stops.</li>
 * <li>{@link math.demo.sthelens.SurfaceDemo} -- a lidar grid over the crater of
 * Mount St. Helens, thinned to a coarser survey and put back six ways through
 * {@code fit}. Bilinear interpolation cannot return a value outside the four
 * samples around it, which sounds like safety; the mountain leaves that range
 * twice as often as the bicubic surface does, and where the bicubic surface
 * leaves it, it is the nearer of the two.</li>
 * </ul>
 */
package math.demo;
