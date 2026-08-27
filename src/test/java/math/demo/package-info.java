/**
 * Worked examples: seven published datasets, each taken apart with the pieces
 * this library provides.
 * <p>
 * Every subpackage holds the same three files -- {@code Datasets} with the
 * numbers as literals, a {@code main} that answers one question about them in
 * numbered steps, and a JUnit test that guards what the {@code main} prints.
 * {@link math.demo.co2state} is the one exception and holds no data of its own:
 * it revisits the Mauna Loa record with a different model, and a second copy of
 * 821 literals could only ever drift away from the first.
 * <p>
 * The data are measurements rather than simulations, so an answer can be wrong
 * in ways generated data never is, and each demo needs several packages at
 * once for a question no single class answers. That is where the seams show,
 * and finding them is what these are for.
 * <ul>
 * <li>{@link math.demo.maunaloa.MaunaLoaDemo} -- sixty-eight years of monthly
 * CO2, in eight steps across {@code probe}, {@code fit}, {@code linalg},
 * {@code solve} and {@code fft}: trend and seasonal cycle, residuals that turn
 * out not to be independent, and when the curve reaches 450 ppm.</li>
 * <li>{@link math.demo.co2state.StateSpaceDemo} -- the same CO2 record as a
 * linear Gaussian state space model, through {@code math.ts} and
 * {@code math.optim}. It exists because of a defect the demo above diagnoses
 * about itself: its residuals are autocorrelated at {@code 0.906}, so its
 * standard errors mean less than they say. Letting the level and the slope
 * wander puts that memory in the model instead, and the innovations come out
 * at {@code 0.030}. The move to Maunakea after the 2022 eruption leaves no
 * trace; the largest anomaly in sixty-eight years is the 2016 El Nino; and
 * seven months in ten can be thrown away for a hundredth of a ppm.</li>
 * <li>{@link math.demo.centralpark.WeatherDemo} -- a year of daily maxima and
 * minima at one weather station, and the first demo whose observation has two
 * components rather than one. The two readings are not two series: they are one
 * temperature read twice, and the model has to decide how much of each is
 * weather. It decides that the night minimum is nearly all of the information
 * and the afternoon maximum nearly all of the noise -- discard the maximum and
 * the weather is 0.3 percent less certain, discard the minimum and it is 2.19
 * times less. Also a tenth of a degree that does not exist, 35 days of rain
 * hidden behind a flag that turn out not to matter, and an archive that
 * rewrites its own last twelve months. Its closing step is the one for
 * modellers: the off-diagonal of {@code R} is worth 153 nats in one model and
 * 1.7 in the next, and lands on the boundary of its own parameter space, which
 * AIC prefers and a modeller should not.</li>
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
