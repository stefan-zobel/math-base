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
package math.optim;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import math.fun.DiffDFunction;
import math.fun.DiffDMultiFunction;
import math.fun.DMultiFunctionEval;

/**
 * Conjugate-gradient implementation based on the code in "Numerical Recipes in
 * C" (see p. 423 and others).
 * <p>
 * As of now, it requires a differentiable function
 * {@link DiffDMultiFunction} as input.
 * <p>
 * The basic way to use the CGOptimizer is with the simple {@code minimize}
 * method:
 * <p>
 * <code>DiffDMultiFunction dmf = new SomeDiffDMultiFunction();</code>
 * <br>
 * <code>double[] initial = getInitialGuess();</code> <br>
 * <code>DMultiFunctionEval minimum = CGOptimizer.minimize(dmf, initial);</code>
 * 
 * @author <a HREF="mailto:klein@cs.stanford.edu">Dan Klein</a>
 */
public final class CGOptimizer {

    private static final DecimalFormat NF = new DecimalFormat("0.000E0", DecimalFormatSymbols.getInstance(Locale.ROOT));

    private static final boolean SIMPLE_GD = false;
    private static final boolean CHECK_SIMPLE_GD_CONVERGENCE = true;

    // constants
    private static final int ITMAX = 10001;
    private static final double EPS = 1.0e-30;

    // The one-dimensional minimizer this used to carry as private helpers. Its
    // settings are the ones that were hard-coded in the extracted dbrent, so
    // the line search behaves as it always did.
    private static final BrentMinimizer LINE_MINIMIZER = new BrentMinimizer(1.0e-4, 100, 200);
    // where a line search starts looking, in units of the direction vector
    private static final double INITIAL_STEP = 0.01;

    private static final int RESET_FREQ = 10;
    // default function tolerance
    private static final double FUNC_DEFAULT_TOL = 1e-10;

    static final class Minimand implements DiffDMultiFunction {
        private final DiffDMultiFunction f;

        Minimand(DiffDMultiFunction f) {
            this.f = f;
        }

        @Override
        public double apply(double[] x) {
            return -f.apply(x);
        }

        @Override
        public void derivativeAt(double[] x, double[] grad) {
            f.derivativeAt(x, grad);
            for (int i = 0; i < grad.length; ++i) {
                grad[i] = -grad[i];
            }
        }
    } // Minimand

    /** The multivariate function restricted to a line, as seen by the minimizer. */
    static final class OneDimDiffFunction implements DiffDFunction {
        private final DiffDMultiFunction function;
        private final double[] initial;
        private final double[] direction;
        private final double[] currVector;
        private final double[] currGradient;

        OneDimDiffFunction(DiffDMultiFunction function, double[] initial,
                double[] direction) {
            this.function = function;
            this.initial = initial.clone();
            this.direction = direction.clone();
            this.currVector = new double[initial.length];
            this.currGradient = new double[initial.length];
        }

        double[] vectorOf(double x) {
            for (int i = 0; i < initial.length; i++) {
                currVector[i] = initial[i] + (x * direction[i]);
            }
            return currVector;
        }

        @Override
        public double apply(double x) {
            return function.apply(vectorOf(x));
        }

        @Override
        public double derivativeAt(double x) {
            function.derivativeAt(vectorOf(x), currGradient);
            double d = 0.0;
            for (int i = 0; i < currGradient.length; i++) {
                d += currGradient[i] * direction[i];
            }
            return d;
        }
    } // OneDimDiffFunction

    private CGOptimizer() {
    }

    public static DMultiFunctionEval maximize(
            DiffDMultiFunction function, double[] initial) {
        return maximize(function, FUNC_DEFAULT_TOL, initial);
    }

    public static DMultiFunctionEval maximize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial) {
        return maximize(function, functionTolerance, initial, ITMAX);
    }

    public static DMultiFunctionEval maximize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial, int maxIterations) {
        return maximize(function, functionTolerance, initial, maxIterations,
                true);
    }

    public static DMultiFunctionEval maximize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial, int maxIterations, boolean silent) {
        return minimize(new Minimand(function), functionTolerance, initial,
                maxIterations, silent, true);
    }

    public static DMultiFunctionEval minimize(
            DiffDMultiFunction function, double[] initial) {
        return minimize(function, FUNC_DEFAULT_TOL, initial);
    }

    public static DMultiFunctionEval minimize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial) {
        return minimize(function, functionTolerance, initial, ITMAX);
    }

    public static DMultiFunctionEval minimize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial, int maxIterations) {
        return minimize(function, functionTolerance, initial, maxIterations,
                true);
    }

    public static DMultiFunctionEval minimize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial, int maxIterations, boolean silent) {
        return minimize(function, functionTolerance, initial, maxIterations,
                silent, false);
    }

    private static DMultiFunctionEval minimize(
            DiffDMultiFunction function, double functionTolerance,
            double[] initial, int maxIterations, boolean silent,
            boolean isMaximization) {

        int dimension = initial.length;
        double sign = isMaximization ? -1.0 : 1.0;

        // evaluate function
        double fp = function.apply(initial);
        double[] xi = new double[dimension];
        function.derivativeAt(initial, xi);

        // make some vectors
        double[] g = new double[dimension];
        double[] h = new double[dimension];
        double[] p = new double[dimension];
        for (int j = 0; j < dimension; j++) {
            g[j] = -xi[j];
            xi[j] = g[j];
            h[j] = g[j];
            p[j] = initial[j];
        }

        // iterations
        boolean simpleGDStep = false;
        int iter = 1;
        while (iter < maxIterations) {
            if (!silent) {
                System.err.print("Iter " + iter + ' ');
            }
            // do a line min along descent direction
            double[] p2 = lineMinimize(function, p, xi);
            if (p2 == null) {
                // The function falls without bound along this direction, so
                // there is nothing to minimize onto. Stop here rather than
                // iterate the budget away on infinities.
                return new DMultiFunctionEval(p, sign * function.apply(p), iter,
                        false);
            }
            double fp2 = function.apply(p2);

            if (!silent) {
                System.err.printf(Locale.ROOT, " %s (delta: %s)\n", NF.format(fp2),
                        NF.format(fp - fp2));
            }

            // check convergence
            if (2.0 * Math.abs(fp2 - fp) <= functionTolerance
                    * (Math.abs(fp2) + Math.abs(fp) + EPS)) {
                // convergence
                if (!CHECK_SIMPLE_GD_CONVERGENCE || simpleGDStep || SIMPLE_GD) {
                    return new DMultiFunctionEval(p, sign * fp2, iter);
                }
                simpleGDStep = true;
            } else {
                simpleGDStep = false;
            }
            // shift variables
            for (int j = 0; j < dimension; j++) {
                xi[j] = p2[j] - p[j];
                p[j] = p2[j];
            }
            fp = fp2;
            // find the new gradient
            function.derivativeAt(p, xi);

            if (!simpleGDStep && !SIMPLE_GD && (iter % RESET_FREQ != 0)) {
                // do the magic -- part i)
                // (calculate some dot products we'll need)
                double dgg = 0.0;
                double gg = 0.0;
                for (int j = 0; j < dimension; j++) {
                    // g dot g
                    gg += g[j] * g[j];
                    // grad dot grad
                    // FR method is:
                    // dgg += x[j]*x[j];
                    // PR method is:
                    dgg += (xi[j] + g[j]) * xi[j];
                }

                // check for miraculous convergence
                if (gg == 0.0) {
                    return new DMultiFunctionEval(p, sign
                            * function.apply(p), iter);
                }

                // do the magic -- part ii)
                // (update the sequence in a way that tries to preserve
                // conjugacy)
                double gam = dgg / gg;
                for (int j = 0; j < dimension; j++) {
                    g[j] = -xi[j];
                    h[j] = g[j] + gam * h[j];
                    xi[j] = h[j];
                }
            } else {
                // miraculous simpleGD convergence
                double xixi = 0.0;
                for (int j = 0; j < dimension; j++) {
                    xixi += xi[j] * xi[j];
                }
                // reset cgd
                for (int j = 0; j < dimension; j++) {
                    g[j] = -xi[j];
                    xi[j] = g[j];
                    h[j] = g[j];
                }
                if (xixi == 0.0) {
                    return new DMultiFunctionEval(p, sign
                            * function.apply(p), iter);
                }
            }
            ++iter;
        } // while

        // too many iterations; the returned converged flag says so
        return new DMultiFunctionEval(p, sign * function.apply(p),
                iter, false);
    }

    /**
     * Minimizes {@code function} along the line through {@code initial} in the
     * given {@code direction}. Returns {@code null} if the function has no
     * minimum along it, which the caller must treat as the end of the run.
     */
    private static double[] lineMinimize(DiffDMultiFunction function,
            double[] initial, double[] direction) {
        OneDimDiffFunction oneDim = new OneDimDiffFunction(function, initial,
                direction);
        BrentMinimizer.Bracket bracket = LINE_MINIMIZER.bracket(oneDim, 0.0,
                INITIAL_STEP);
        if (!bracket.bracketed) {
            return null;
        }
        BrentMinimizer.Result min = LINE_MINIMIZER.minimize(oneDim, bracket);
        double xmin = min.x;
        // A search that ran out of iterations has not shown that it improved on
        // where the line started, so do not step unless it did. On the
        // converged path the improvement is established and this costs nothing.
        if (!min.converged && !(min.value < oneDim.apply(0.0))) {
            xmin = 0.0;
        }
        return oneDim.vectorOf(xmin);
    }
}
