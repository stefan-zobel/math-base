/*
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
 * Routines for function fitting and data interpolation.
 * <p>
 * All three interpolators here produce a {@link math.fit.CubicSpline}, which
 * evaluates the value, the first and second derivative and the definite
 * integral in closed form. Choose {@link math.fit.SplineInterpolator} for
 * smooth data, where its two continuous derivatives make it the smoothest and,
 * on coarse grids, the most accurate; {@link math.fit.KrugerInterpolator}
 * where the shape of the data has to be preserved, since it is the only one of
 * the three that cannot leave the range of the values; and
 * {@link math.fit.AkimaInterpolator} for data with a few outliers, since it is
 * the only one whose slopes are local, so a single bad value disturbs six
 * segments rather than every segment there is.
 */
package math.fit;
