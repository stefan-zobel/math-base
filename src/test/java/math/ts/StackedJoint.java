package math.ts;

import math.distribution.MultivariateNormal;
import math.linalg.CholeskyDecomp;
import math.linalg.DMatrix;

/**
 * The oracle for {@link KalmanFilter} and {@link RtsSmoother}: the whole series
 * written out as one multivariate normal.
 * <p>
 * A linear Gaussian state space model makes the stacked vector
 * <code>(x[0], ..., x[n-1], y[0], ..., y[n-1])</code> jointly normal, and its
 * mean and covariance follow from the model without any recursion over the
 * data:
 *
 * <pre>
 * mu_x[t]        = F(t) mu_x[t-1]
 * Cov(x[s],x[t]) = Cov(x[s],x[t-1]) F(t)'  for t &gt; s
 * Cov(x[t],x[t]) = F(t) Cov(x[t-1],x[t-1]) F(t)' + Q(t)
 * Cov(y[s],y[t]) = H(s) Cov(x[s],x[t]) H(t)' + [s == t] R(t)
 * </pre>
 *
 * So the log likelihood of the series is <b>one</b> call to
 * {@link MultivariateNormal#logPdf(double[])} on a vector of length
 * {@code n * m}, and the filtered and smoothed states are what conditioning
 * that joint on the observations gives. Neither route shares a line of code
 * with the classes under test, which is the point: the forward recursion is
 * checked against a formula that does not recurse, and the backward pass
 * against ordinary Gaussian conditioning.
 * <p>
 * It costs a factorization of an {@code n * m} matrix per query, so the series
 * has to stay short. Missing entries are handled by striking their rows and
 * columns out, which is why the same oracle covers the {@code NaN} cases.
 */
final class StackedJoint {

    private final int n;
    private final int d;
    private final int m;
    private final double[][] muX;
    private final DMatrix[][] covX;
    private final DMatrix[][] crossXY;
    private final double[] muY;
    private final DMatrix sigma;

    StackedJoint(StateSpaceModel model, int n) {
        this.n = n;
        this.d = model.stateDimension();
        this.m = model.observationDimension();

        muX = new double[n][d];
        model.initialMean(muX[0]);
        for (int t = 1; t < n; ++t) {
            DMatrix f = model.transition(t);
            for (int i = 0; i < d; ++i) {
                double sum = 0.0;
                for (int j = 0; j < d; ++j) {
                    sum += f.get(i, j) * muX[t - 1][j];
                }
                muX[t][i] = sum;
            }
        }

        covX = new DMatrix[n][n];
        covX[0][0] = model.initialCovariance();
        for (int t = 1; t < n; ++t) {
            DMatrix f = model.transition(t);
            covX[t][t] = f.mul(covX[t - 1][t - 1]).mulBTrans(f).add(model.processNoise(t));
        }
        for (int s = 0; s < n; ++s) {
            for (int t = s + 1; t < n; ++t) {
                covX[s][t] = covX[s][t - 1].mulBTrans(model.transition(t));
                covX[t][s] = covX[s][t].transpose();
            }
        }

        crossXY = new DMatrix[n][n];
        for (int s = 0; s < n; ++s) {
            for (int t = 0; t < n; ++t) {
                crossXY[s][t] = covX[s][t].mulBTrans(model.observation(t));
            }
        }

        muY = new double[n * m];
        for (int t = 0; t < n; ++t) {
            DMatrix h = model.observation(t);
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int j = 0; j < d; ++j) {
                    sum += h.get(i, j) * muX[t][j];
                }
                muY[t * m + i] = sum;
            }
        }

        sigma = new DMatrix(n * m, n * m);
        for (int s = 0; s < n; ++s) {
            for (int t = 0; t < n; ++t) {
                DMatrix block = model.observation(s).mul(crossXY[s][t]);
                DMatrix r = model.observationNoise(t);
                for (int i = 0; i < m; ++i) {
                    for (int j = 0; j < m; ++j) {
                        double v = block.get(i, j) + (s == t ? r.get(i, j) : 0.0);
                        sigma.set(s * m + i, t * m + j, v);
                    }
                }
            }
        }
    }

    /** The flat positions of the entries of {@code y} that were observed. */
    private int[] observed(DMatrix y, int upto) {
        int k = 0;
        for (int t = 0; t < upto; ++t) {
            for (int i = 0; i < m; ++i) {
                if (!Double.isNaN(y.get(t, i))) {
                    ++k;
                }
            }
        }
        int[] out = new int[k];
        int at = 0;
        for (int t = 0; t < upto; ++t) {
            for (int i = 0; i < m; ++i) {
                if (!Double.isNaN(y.get(t, i))) {
                    out[at++] = t * m + i;
                }
            }
        }
        return out;
    }

    private DMatrix subCovariance(int[] idx) {
        DMatrix out = new DMatrix(idx.length, idx.length);
        for (int a = 0; a < idx.length; ++a) {
            for (int b = 0; b < idx.length; ++b) {
                out.set(a, b, sigma.get(idx[a], idx[b]));
            }
        }
        return out;
    }

    private double[] residual(DMatrix y, int[] idx) {
        double[] out = new double[idx.length];
        for (int a = 0; a < idx.length; ++a) {
            out[a] = y.get(idx[a] / m, idx[a] % m) - muY[idx[a]];
        }
        return out;
    }

    /** The log density of everything that was actually observed. */
    double logLikelihood(DMatrix y) {
        int[] idx = observed(y, n);
        if (idx.length == 0) {
            return 0.0;
        }
        double[] mean = new double[idx.length];
        double[] value = new double[idx.length];
        for (int a = 0; a < idx.length; ++a) {
            mean[a] = muY[idx[a]];
            value[a] = y.get(idx[a] / m, idx[a] % m);
        }
        return new MultivariateNormal(mean, subCovariance(idx)).logPdf(value);
    }

    /**
     * <code>E[x[t] | the observed entries among the first upto time points]</code>,
     * by ordinary Gaussian conditioning.
     */
    double[] conditionalMean(DMatrix y, int t, int upto) {
        int[] idx = observed(y, upto);
        double[] out = muX[t].clone();
        if (idx.length == 0) {
            return out;
        }
        DMatrix sub = subCovariance(idx);
        double[] alpha = new double[idx.length];
        CholeskyDecomp.solve(CholeskyDecomp.cholesky(sub), residual(y, idx), alpha);
        double[][] cross = crossRow(t, idx);
        for (int j = 0; j < d; ++j) {
            double sum = 0.0;
            for (int a = 0; a < idx.length; ++a) {
                sum += cross[j][a] * alpha[a];
            }
            out[j] += sum;
        }
        return out;
    }

    /**
     * <code>Cov(x[t] | the observed entries among the first upto time points)</code>.
     */
    DMatrix conditionalCovariance(DMatrix y, int t, int upto) {
        int[] idx = observed(y, upto);
        DMatrix out = covX[t][t].copy();
        if (idx.length == 0) {
            return out;
        }
        DMatrix factor = CholeskyDecomp.cholesky(subCovariance(idx));
        double[][] cross = crossRow(t, idx);
        double[][] solved = new double[d][];
        for (int j = 0; j < d; ++j) {
            solved[j] = new double[idx.length];
            CholeskyDecomp.solve(factor, cross[j], solved[j]);
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int a = 0; a < idx.length; ++a) {
                    sum += cross[i][a] * solved[j][a];
                }
                out.set(i, j, out.get(i, j) - sum);
            }
        }
        return out;
    }

    /**
     * <code>Cov(x[s], x[t] | the observed entries among the first upto time
     * points)</code>. Not symmetric unless {@code s == t}, which is the point
     * of having it: it pins which way round the lag-one covariance goes.
     */
    DMatrix conditionalCrossCovariance(DMatrix y, int s, int t, int upto) {
        int[] idx = observed(y, upto);
        DMatrix out = covX[s][t].copy();
        if (idx.length == 0) {
            return out;
        }
        DMatrix factor = CholeskyDecomp.cholesky(subCovariance(idx));
        double[][] crossS = crossRow(s, idx);
        double[][] crossT = crossRow(t, idx);
        double[][] solved = new double[d][];
        for (int j = 0; j < d; ++j) {
            solved[j] = new double[idx.length];
            CholeskyDecomp.solve(factor, crossT[j], solved[j]);
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int a = 0; a < idx.length; ++a) {
                    sum += crossS[i][a] * solved[j][a];
                }
                out.set(i, j, out.get(i, j) - sum);
            }
        }
        return out;
    }

    /** <code>Cov(x[t], y)</code> over the selected entries, as {@code d} rows. */
    private double[][] crossRow(int t, int[] idx) {
        double[][] out = new double[d][idx.length];
        for (int a = 0; a < idx.length; ++a) {
            int s = idx[a] / m;
            int i = idx[a] % m;
            for (int j = 0; j < d; ++j) {
                out[j][a] = crossXY[t][s].get(j, i);
            }
        }
        return out;
    }
}
