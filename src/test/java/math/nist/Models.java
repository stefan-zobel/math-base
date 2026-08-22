package math.nist;

import math.fun.DiffDVectorFunction;

/**
 * The five nonlinear models of {@link StRD#nonlinear()}, each as the residuals
 * and the analytic Jacobian a Levenberg-Marquardt solver wants.
 * <p>
 * The formulas are the ones the NIST files print, and the derivatives are taken
 * by hand rather than by difference, because a certification of the solver
 * should not also be a certification of a difference quotient. The Jacobian is
 * column-major, {@code jacobian[j * m + i] = d r_i / d b_j}, which is the
 * layout {@code math.optim} uses throughout.
 * <p>
 * <b>The exponentials go through {@code StrictMath}, not {@code Math}.</b>
 * {@code Math.exp} is allowed to be wrong by one unit in the last place and its
 * intrinsic differs between JDK versions; {@code StrictMath.exp} is specified
 * to return the same value everywhere. On MGH10 that one unit is worth two
 * digits of the fitted parameters -- 9.0 against 10.9 against the certified
 * values, on JDK 8 and JDK 25 respectively -- because the objective is flat
 * where the minimum sits, so the last bit of the model moves the answer far
 * more than it moves the residual. With {@code StrictMath} the fit is
 * bit-identical on both. A certification that cannot reproduce itself is not
 * one.
 */
public final class Models {

    /**
     * The residuals and Jacobian of one reference set.
     *
     * @param set
     *            the set, whose {@code name} selects the formula
     * @return the model, ready for {@code LevenbergMarquardt.solve}
     * @throws IllegalArgumentException
     *             if the set is not one of the five
     */
    public static DiffDVectorFunction of(final StRD.NonlinearSet set) {
        final double[] x = set.x();
        final double[] y = set.y();
        final int m = x.length;
        final String name = set.name;
        if (!known(name)) {
            throw new IllegalArgumentException("no model for " + name);
        }
        return new DiffDVectorFunction() {
            @Override
            public void valueAt(double[] b, double[] residuals) {
                for (int i = 0; i < m; ++i) {
                    residuals[i] = value(name, b, x[i]) - y[i];
                }
            }

            @Override
            public void jacobianAt(double[] b, double[] jacobian) {
                double[] row = new double[b.length];
                for (int i = 0; i < m; ++i) {
                    gradient(name, b, x[i], row);
                    for (int j = 0; j < row.length; ++j) {
                        jacobian[j * m + i] = row[j];
                    }
                }
            }
        };
    }

    private static boolean known(String name) {
        return name.equals("Misra1a") || name.equals("Chwirut1") || name.equals("Thurber") || name.equals("MGH09")
                || name.equals("MGH10");
    }

    /** The model value at one predictor. */
    static double value(String name, double[] b, double x) {
        if (name.equals("Misra1a")) {
            // y = b1 (1 - exp(-b2 x))
            return b[0] * (1.0 - StrictMath.exp(-b[1] * x));
        }
        if (name.equals("Chwirut1")) {
            // y = exp(-b1 x) / (b2 + b3 x)
            return StrictMath.exp(-b[0] * x) / (b[1] + b[2] * x);
        }
        if (name.equals("Thurber")) {
            // y = (b1 + b2 x + b3 x^2 + b4 x^3) / (1 + b5 x + b6 x^2 + b7 x^3)
            return numerator(b, x) / denominator(b, x);
        }
        if (name.equals("MGH09")) {
            // y = b1 (x^2 + b2 x) / (x^2 + b3 x + b4)
            return b[0] * (x * x + x * b[1]) / (x * x + x * b[2] + b[3]);
        }
        // MGH10: y = b1 exp(b2 / (x + b3))
        return b[0] * StrictMath.exp(b[1] / (x + b[2]));
    }

    /** The gradient of the model value with respect to the parameters. */
    static void gradient(String name, double[] b, double x, double[] out) {
        if (name.equals("Misra1a")) {
            double e = StrictMath.exp(-b[1] * x);
            out[0] = 1.0 - e;
            out[1] = b[0] * x * e;
            return;
        }
        if (name.equals("Chwirut1")) {
            double e = StrictMath.exp(-b[0] * x);
            double d = b[1] + b[2] * x;
            out[0] = -x * e / d;
            out[1] = -e / (d * d);
            out[2] = -x * e / (d * d);
            return;
        }
        if (name.equals("Thurber")) {
            double n = numerator(b, x);
            double d = denominator(b, x);
            double dd = d * d;
            out[0] = 1.0 / d;
            out[1] = x / d;
            out[2] = x * x / d;
            out[3] = x * x * x / d;
            out[4] = -n * x / dd;
            out[5] = -n * x * x / dd;
            out[6] = -n * x * x * x / dd;
            return;
        }
        if (name.equals("MGH09")) {
            double n = x * x + x * b[1];
            double d = x * x + x * b[2] + b[3];
            out[0] = n / d;
            out[1] = b[0] * x / d;
            out[2] = -b[0] * n * x / (d * d);
            out[3] = -b[0] * n / (d * d);
            return;
        }
        double e = StrictMath.exp(b[1] / (x + b[2]));
        double shifted = x + b[2];
        out[0] = e;
        out[1] = b[0] * e / shifted;
        out[2] = -b[0] * e * b[1] / (shifted * shifted);
    }

    private static double numerator(double[] b, double x) {
        return b[0] + x * (b[1] + x * (b[2] + x * b[3]));
    }

    private static double denominator(double[] b, double x) {
        return 1.0 + x * (b[4] + x * (b[5] + x * b[6]));
    }

    private Models() {
        throw new AssertionError();
    }
}
