package math.ode;

import math.fun.DVectorField;

/**
 * An explicit Runge-Kutta step, for whichever {@link ButcherTableau} it is
 * constructed with.
 * <p>
 * One class rather than one per method: the tableau is the only thing that
 * differs, so a further method is a set of coefficients and the code that
 * executes them stays the single place where a step can be got wrong. The cost
 * is a loop over the stage coefficients where straight-line code would do, and
 * it is paid against an evaluation of the field, which for any equation worth
 * solving numerically dominates it.
 * <p>
 * <b>Two evaluations are saved where they can be.</b> A rejected step is
 * retried from the same point, and the first stage there is what it already
 * was, so it is not computed twice. A method that is first same as last --
 * {@link ButcherTableau#isFsal()}, which Dormand-Prince is -- ends a step with
 * the field evaluated at exactly the state and time the next step begins at,
 * and that stage carries over. Both are recognized by comparing the point,
 * never by a protocol the driver has to follow, so a driver that steps in an
 * order this class did not anticipate loses the saving and not the answer.
 * <p>
 * <b>Interpolation is always available and the endpoints are exact.</b> With a
 * continuous extension in the tableau it is a polynomial over the stages, and
 * for Dormand-Prince those are the stages the step computed anyway, so it costs
 * nothing further. Without one it is the cubic through the two states and the
 * two derivatives. A method whose interpolant needs stages of its own --
 * {@link ButcherTableau#denseStages()} above {@link ButcherTableau#stages()},
 * which DOP853 is by three -- pays for them, and so does the cubic for the
 * derivative at the end of the step. <b>Either way that cost is deferred to the
 * first interior point asked for inside a step and paid once</b>, so a caller
 * who never interpolates never pays it, and one who asks for a hundred points
 * in a step pays it once.
 * <p>
 * Either way the interpolant is <b>one order below the step it spans</b>, which
 * is what the continuous extension of Dormand-Prince buys over the cubic: it is
 * one order below a fifth order step rather than one order below a fourth order
 * one. Measured on the harmonic oscillator over five halvings, the worst error
 * inside a step falls by {@code 2^5.0} for Dormand-Prince, by {@code 2^4.0} for
 * the classical method and by {@code 2^8.0} for DOP853 -- the last being what
 * its three extra stages are bought with, and the reason they are worth
 * evaluating when they are wanted.
 * <p>
 * Instances are stateful and cannot be shared between threads.
 * <p>
 * <b>See</b>
 * <a href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">Wikipedia
 * Runge-Kutta methods</a>.
 *
 * @since 1.5.3
 */
public final class ExplicitRungeKutta implements OdeStepper {

    private final ButcherTableau tableau;
    private final DVectorField field;
    private final int n;

    private final double[] c;
    private final double[][] a;
    private final double[] b;
    private final double[] bStar;
    private final double[] bStarSecondary;
    private final double[][] denseB;
    private final int propagating;
    private final int allStages;
    private final boolean fsal;

    private final double[][] k;
    private final double[] stage;
    private final double[] secondError;
    private final double[] endDerivative;
    private final double[] stiffArgument;
    private final int stiffEarly;
    private final int stiffLate;

    private final double[] startY;
    private final double[] endY;
    private final double[] lastYOut;

    private double startT;
    private double endT;
    private double lastH;
    private boolean startValid;
    private boolean endValid;
    private boolean stepped;
    private boolean endDerivativeValid;
    private boolean denseStagesValid;
    private boolean secondErrorValid;
    private long evaluations;

    /**
     * A stepper for the given method and equation.
     *
     * @param tableau
     *            the coefficients of the method
     * @param field
     *            the right hand side of the equation
     * @param dimension
     *            the number of components of the state
     * @throws IllegalArgumentException
     *             if an argument is {@code null} or the dimension is not
     *             positive
     */
    public ExplicitRungeKutta(ButcherTableau tableau, DVectorField field, int dimension) {
        if (tableau == null) {
            throw new IllegalArgumentException("tableau must not be null");
        }
        if (field == null) {
            throw new IllegalArgumentException("field must not be null");
        }
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be at least 1, not " + dimension);
        }
        this.tableau = tableau;
        this.field = field;
        this.n = dimension;
        this.c = tableau.c();
        this.a = tableau.a();
        this.b = tableau.b();
        this.bStar = tableau.bStar();
        this.bStarSecondary = tableau.bStarSecondary();
        this.denseB = tableau.dense();
        this.propagating = tableau.stages();
        this.allStages = tableau.denseStages();
        this.fsal = tableau.isFsal();
        this.k = new double[allStages][dimension];
        this.stage = new double[dimension];
        this.secondError = (bStarSecondary == null) ? null : new double[dimension];
        this.endDerivative = new double[dimension];
        this.startY = new double[dimension];
        this.endY = new double[dimension];
        this.lastYOut = new double[dimension];

        // the two last stages that sit at the end of the step: their states
        // differ and their times do not, which is a difference quotient of the
        // right hand side and so an estimate of the dominant eigenvalue
        int late = -1;
        int early = -1;
        for (int i = propagating - 1; i >= 0 && early < 0; --i) {
            if (Math.abs(c[i] - 1.0) <= 1.0e-12) {
                if (late < 0) {
                    late = i;
                } else {
                    early = i;
                }
            }
        }
        this.stiffLate = late;
        this.stiffEarly = early;
        this.stiffArgument = (early < 0) ? null : new double[dimension];
    }

    /**
     * The coefficients this stepper executes.
     *
     * @return the tableau given to the constructor
     */
    public ButcherTableau tableau() {
        return tableau;
    }

    @Override
    public int dimension() {
        return n;
    }

    @Override
    public int order() {
        return tableau.order();
    }

    @Override
    public boolean hasErrorEstimate() {
        return bStar != null;
    }

    /**
     * Whether the method carries a continuous extension of its own, over stages
     * it has computed anyway. {@link #interpolate(double, double[])} works
     * either way; a {@code false} here means it is the cubic fallback, which
     * costs one evaluation and is of order three.
     *
     * @return {@code true} if the tableau has dense output coefficients
     */
    @Override
    public boolean hasDenseOutput() {
        return denseB != null;
    }

    @Override
    public void derivative(double t, double[] y, double[] dydt) {
        if (y == null || y.length != n) {
            throw new IllegalArgumentException("y must be of length " + n);
        }
        if (dydt == null || dydt.length != n) {
            throw new IllegalArgumentException("dydt must be of length " + n);
        }
        field.valueAt(t, y, dydt);
        ++evaluations;
    }

    @Override
    public void step(double t, double[] y, double h, double[] yOut, double[] errOut) {
        if (y == null || y.length != n) {
            throw new IllegalArgumentException("y must be of length " + n);
        }
        if (yOut == null || yOut.length != n) {
            throw new IllegalArgumentException("yOut must be of length " + n);
        }
        if (yOut == y) {
            throw new IllegalArgumentException("yOut must not be the array passed as y");
        }
        if (errOut != null && errOut.length != n) {
            throw new IllegalArgumentException("errOut must be of length " + n);
        }
        if (h == 0.0) {
            throw new IllegalArgumentException("h must not be zero");
        }

        if (startValid && t == startT && sameContent(y, startY)) {
            // k[0] is still the field at this very point, from the step that
            // was rejected
        } else if (fsal && endValid && t == endT && sameContent(y, endY)) {
            System.arraycopy(k[propagating - 1], 0, k[0], 0, n);
        } else {
            field.valueAt(t, y, k[0]);
            ++evaluations;
        }
        startT = t;
        System.arraycopy(y, 0, startY, 0, n);
        startValid = true;

        // only what advances the solution: any stage beyond that exists for the
        // interpolant alone and is evaluated when one is actually asked for
        denseStagesValid = allStages == propagating;
        for (int i = 1; i < propagating; ++i) {
            System.arraycopy(y, 0, stage, 0, n);
            double[] ai = a[i];
            for (int j = 0; j < i; ++j) {
                double aij = ai[j];
                if (aij != 0.0) {
                    double w = h * aij;
                    double[] kj = k[j];
                    for (int m = 0; m < n; ++m) {
                        stage[m] += w * kj[m];
                    }
                }
            }
            if (i == stiffEarly) {
                System.arraycopy(stage, 0, stiffArgument, 0, n);
            }
            field.valueAt(t + c[i] * h, stage, k[i]);
            ++evaluations;
        }

        System.arraycopy(y, 0, yOut, 0, n);
        for (int i = 0; i < propagating; ++i) {
            double bi = b[i];
            if (bi != 0.0) {
                double w = h * bi;
                double[] ki = k[i];
                for (int m = 0; m < n; ++m) {
                    yOut[m] += w * ki[m];
                }
            }
        }

        secondErrorValid = false;
        if (errOut != null && bStar != null) {
            weigh(errOut, bStar, h);
            if (bStarSecondary != null) {
                weigh(secondError, bStarSecondary, h);
                secondErrorValid = true;
            }
        }

        lastH = h;
        System.arraycopy(yOut, 0, lastYOut, 0, n);
        stepped = true;
        endDerivativeValid = false;
        if (fsal) {
            endT = t + h;
            System.arraycopy(yOut, 0, endY, 0, n);
            endValid = true;
        }
    }

    @Override
    public void interpolate(double theta, double[] out) {
        if (!stepped) {
            throw new IllegalStateException("no step has been taken");
        }
        if (out == null || out.length != n) {
            throw new IllegalArgumentException("out must be of length " + n);
        }
        if (theta == 0.0) {
            System.arraycopy(startY, 0, out, 0, n);
            return;
        }
        if (theta == 1.0) {
            System.arraycopy(lastYOut, 0, out, 0, n);
            return;
        }
        if (denseB != null) {
            ensureDenseStages();
            System.arraycopy(startY, 0, out, 0, n);
            for (int i = 0; i < allStages; ++i) {
                double[] p = denseB[i];
                double w = 0.0;
                for (int j = p.length - 1; j >= 0; --j) {
                    w = w * theta + p[j];
                }
                w *= theta * lastH;
                if (w != 0.0) {
                    double[] ki = k[i];
                    for (int m = 0; m < n; ++m) {
                        out[m] += w * ki[m];
                    }
                }
            }
            return;
        }

        if (!endDerivativeValid) {
            if (fsal) {
                System.arraycopy(k[propagating - 1], 0, endDerivative, 0, n);
            } else {
                field.valueAt(startT + lastH, lastYOut, endDerivative);
                ++evaluations;
            }
            endDerivativeValid = true;
        }
        double u2 = theta * theta;
        double u3 = u2 * theta;
        double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
        double h10 = u3 - 2.0 * u2 + theta;
        double h01 = -2.0 * u3 + 3.0 * u2;
        double h11 = u3 - u2;
        double[] k0 = k[0];
        for (int m = 0; m < n; ++m) {
            out[m] = h00 * startY[m] + h01 * lastYOut[m] + lastH * (h10 * k0[m] + h11 * endDerivative[m]);
        }
    }

    /**
     * The difference quotient between the last two stages that sit at the end
     * of the step, times the length of the step.
     * <p>
     * Both are evaluations of the field at the same time and at states that
     * differ only by the error of the lower half of the embedded pair, so their
     * ratio is a directional derivative in the direction the error happens to
     * point -- which is the direction the fastest mode grows in, since that is
     * what makes the error. Dormand-Prince has such a pair and the classical
     * method does not, so this answers {@link Double#NaN} there.
     *
     * @return the estimate, or {@link Double#NaN} if the tableau has no two
     *         stages at the end of the step, or no step has been taken
     */
    @Override
    public double stiffnessMeasure() {
        if (stiffEarly < 0 || !stepped) {
            return Double.NaN;
        }
        double numerator = 0.0;
        double denominator = 0.0;
        double[] late = k[stiffLate];
        double[] early = k[stiffEarly];
        for (int m = 0; m < n; ++m) {
            double dk = late[m] - early[m];
            double dy = lastYOut[m] - stiffArgument[m];
            numerator += dk * dk;
            denominator += dy * dy;
        }
        if (!(denominator > 0.0)) {
            return Double.NaN;
        }
        return Math.abs(lastH) * Math.sqrt(numerator / denominator);
    }

    @Override
    public double[] secondaryError() {
        return secondErrorValid ? secondError : null;
    }

    @Override
    public double stiffnessThreshold() {
        return tableau.stiffnessThreshold();
    }

    @Override
    public long evaluations() {
        return evaluations;
    }

    @Override
    public void reset() {
        startValid = false;
        endValid = false;
        stepped = false;
        endDerivativeValid = false;
        denseStagesValid = false;
        secondErrorValid = false;
    }

    /**
     * The method and the equation this stepper advances.
     */
    @Override
    public String toString() {
        return "ExplicitRungeKutta[" + tableau + ", dimension " + n + "]";
    }

    /**
     * The difference between the advancing solution and an embedded one, over
     * the stages that propagate, which is that pair's estimate of the error of
     * the step.
     */
    private void weigh(double[] out, double[] embedded, double h) {
        for (int m = 0; m < n; ++m) {
            out[m] = 0.0;
        }
        for (int i = 0; i < propagating; ++i) {
            double di = b[i] - embedded[i];
            if (di != 0.0) {
                double w = h * di;
                double[] ki = k[i];
                for (int m = 0; m < n; ++m) {
                    out[m] += w * ki[m];
                }
            }
        }
    }

    /**
     * Evaluates the stages that only the continuous extension needs, once per
     * step and only if an interior value is asked for.
     */
    private void ensureDenseStages() {
        if (denseStagesValid) {
            return;
        }
        for (int i = propagating; i < allStages; ++i) {
            System.arraycopy(startY, 0, stage, 0, n);
            double[] ai = a[i];
            for (int j = 0; j < i; ++j) {
                double aij = ai[j];
                if (aij != 0.0) {
                    double w = lastH * aij;
                    double[] kj = k[j];
                    for (int m = 0; m < n; ++m) {
                        stage[m] += w * kj[m];
                    }
                }
            }
            field.valueAt(startT + c[i] * lastH, stage, k[i]);
            ++evaluations;
        }
        denseStagesValid = true;
    }

    private static boolean sameContent(double[] x, double[] y) {
        for (int i = 0; i < x.length; ++i) {
            if (x[i] != y[i]) {
                return false;
            }
        }
        return true;
    }
}
