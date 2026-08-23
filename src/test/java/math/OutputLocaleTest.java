package math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.function.Supplier;

import org.junit.Test;

import math.fun.DFunction;
import math.linalg.DMatrix;
import math.linalg.FMatrix;
import math.linalg.Lasso;
import math.linalg.Ridge;
import math.list.DoubleArrayList;
import math.optim.BrentMinimizer;
import math.probe.Bootstrap;
import math.probe.DoubleStatistics;
import math.probe.SampleStatistic;
import math.solve.AdaptiveGaussKronrod;
import math.solve.ClenshawCurtis;
import math.solve.DoubleExponential;

/**
 * Everything this library prints has to read the same on every machine. A
 * {@code String.format} without a {@link Locale} takes the platform default.
 * <p>
 * The check is not a pattern match on the text but the invariant itself: the
 * same call under {@link Locale#GERMANY} and under {@link Locale#ROOT} has to
 * produce the identical string.
 */
public class OutputLocaleTest {

    /** A locale whose decimal separator is a comma and whose grouping differs. */
    private static final Locale COMMA = Locale.GERMANY;

    @Test
    public void everyPrintedNumberIsLocaleIndependent() {
        List<String> names = new ArrayList<String>();
        List<Supplier<String>> printers = new ArrayList<Supplier<String>>();

        final DMatrix a = new DMatrix(2, 2);
        a.set(0, 0, 1.5).set(0, 1, -2.25).set(1, 0, 1000.125).set(1, 1, 0.0009);
        add(names, printers, "DMatrix.toString", new Supplier<String>() {
            @Override
            public String get() {
                return a.toString();
            }
        });

        final FMatrix f = new FMatrix(2, 2);
        f.set(0, 0, 1.5f).set(0, 1, -2.25f).set(1, 0, 1000.125f).set(1, 1, 0.0009f);
        add(names, printers, "FMatrix.toString", new Supplier<String>() {
            @Override
            public String get() {
                return f.toString();
            }
        });

        final DoubleArrayList list = new DoubleArrayList(new double[] { 1.5, -2.25, 1000.125 });
        add(names, printers, "DoubleArrayList.toString", new Supplier<String>() {
            @Override
            public String get() {
                return list.toString();
            }
        });

        for (int t = 0; t < 2; ++t) {
            final DoubleStatistics stats = DoubleStatistics.newInstance(t == 1);
            stats.accept(1.5);
            stats.accept(-2.25);
            stats.accept(1000.125);
            add(names, printers, "DoubleStatistics(threadsafe=" + (t == 1) + ").toString", new Supplier<String>() {
                @Override
                public String get() {
                    return stats.toString();
                }
            });
        }

        double[] sample = new double[64];
        long lcg = 99L;
        for (int i = 0; i < sample.length; ++i) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            sample[i] = (lcg >>> 11) * 0x1.0p-53;
        }
        final Bootstrap bootstrap = new Bootstrap(sample, new SampleStatistic() {
            @Override
            public double apply(double[] values) {
                double sum = 0.0;
                for (int i = 0; i < values.length; ++i) {
                    sum += values[i];
                }
                return sum / values.length;
            }
        }, 500, 4242L);
        add(names, printers, "Bootstrap.summary", new Supplier<String>() {
            @Override
            public String get() {
                return bootstrap.summary(0.95);
            }
        });
        add(names, printers, "Bootstrap.toString", new Supplier<String>() {
            @Override
            public String get() {
                return bootstrap.toString();
            }
        });

        final DMatrix x = design(40, 3);
        final DMatrix y = regressand(x);
        add(names, printers, "Lasso.Result.toString", new Supplier<String>() {
            @Override
            public String get() {
                return Lasso.estimate(x, y, 0.05).toString();
            }
        });
        add(names, printers, "Lasso.Path.toString", new Supplier<String>() {
            @Override
            public String get() {
                return Lasso.path(x, y, 1.0, 12, 0.01).toString();
            }
        });
        add(names, printers, "Lasso.CvResult.toString", new Supplier<String>() {
            @Override
            public String get() {
                return Lasso.cv(x, y, 1.0, 4, 7L).toString();
            }
        });
        add(names, printers, "Ridge.Result.toString", new Supplier<String>() {
            @Override
            public String get() {
                return Ridge.estimate(x, y, 0.5).toString();
            }
        });

        final DFunction bell = new DFunction() {
            @Override
            public double apply(double t) {
                return Math.exp(-t * t) + 0.25 * t;
            }
        };
        add(names, printers, "AdaptiveGaussKronrod.IntegralResult.toString", new Supplier<String>() {
            @Override
            public String get() {
                return AdaptiveGaussKronrod.integrate1D(AdaptiveGaussKronrod.G7_K15.POINTS_15, bell, 0.0, 1.0)
                        .toString();
            }
        });
        add(names, printers, "ClenshawCurtis.IntegralResult.toString", new Supplier<String>() {
            @Override
            public String get() {
                return ClenshawCurtis.integrate1D(bell, 0.0, 1.0, 1.0e-10).toString();
            }
        });
        add(names, printers, "DoubleExponential.IntegralResult.toString", new Supplier<String>() {
            @Override
            public String get() {
                return DoubleExponential.integrate1D(bell, 0.0, 1.0, 1.0e-10).toString();
            }
        });
        add(names, printers, "BrentMinimizer.Bracket.toString", new Supplier<String>() {
            @Override
            public String get() {
                return new BrentMinimizer().bracket(bell, 0.0, 1.0).toString();
            }
        });
        add(names, printers, "BrentMinimizer.Result.toString", new Supplier<String>() {
            @Override
            public String get() {
                BrentMinimizer minimizer = new BrentMinimizer();
                return minimizer.minimize(bell, minimizer.bracket(bell, 0.0, 1.0)).toString();
            }
        });

        for (int i = 0; i < printers.size(); ++i) {
            String neutral = under(Locale.ROOT, printers.get(i));
            String german = under(COMMA, printers.get(i));
            assertEquals(names.get(i) + " depends on the default locale", neutral, german);
            assertFalse(names.get(i) + " printed a decimal comma: " + german, german.matches("(?s).*\\d,\\d.*"));
        }
    }

    private static void add(List<String> names, List<Supplier<String>> printers, String name,
            Supplier<String> printer) {
        names.add(name);
        printers.add(printer);
    }

    /** Evaluates {@code printer} with {@code locale} installed as the default. */
    private static String under(Locale locale, Supplier<String> printer) {
        Locale previous = Locale.getDefault();
        try {
            Locale.setDefault(locale);
            return printer.get();
        } finally {
            Locale.setDefault(previous);
        }
    }

    /** A design matrix with no constant column, from an inline LCG. */
    private static DMatrix design(int rows, int columns) {
        DMatrix x = new DMatrix(rows, columns);
        long lcg = 20260821L;
        for (int j = 0; j < columns; ++j) {
            for (int i = 0; i < rows; ++i) {
                lcg = lcg * 6364136223846793005L + 1442695040888963407L;
                x.set(i, j, ((lcg >>> 11) * 0x1.0p-53) * 2.0 - 1.0);
            }
        }
        return x;
    }

    private static DMatrix regressand(DMatrix x) {
        DMatrix y = new DMatrix(x.numRows(), 1);
        for (int i = 0; i < x.numRows(); ++i) {
            y.set(i, 0, 2.0 * x.get(i, 0) - 0.5 * x.get(i, 1) + 0.125);
        }
        return y;
    }
}
