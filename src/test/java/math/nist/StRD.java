package math.nist;

/**
 * The reference datasets of the NIST Statistical Reference Datasets collection,
 * embedded as literals together with the values NIST certifies for them.
 * <p>
 * Source:
 * <a href="https://www.itl.nist.gov/div898/strd/">NIST/ITL StRD</a>, retrieved
 * 2026-08-22 -- the linear sets from {@code strd/lls/data/LINKS/DATA} and the
 * nonlinear ones from {@code strd/nls/data/LINKS/DATA}. A publication of the
 * United States government, released expressly so that numerical software can
 * be tested against it, and in the public domain.
 * <p>
 * Everything else in this repository is checked against itself: an invariant, a
 * second algorithm, a closed form. Both sides of such a check can be wrong
 * together. These values were computed elsewhere, in extended precision, by
 * people whose purpose was to have something to check against -- they are the
 * only external oracle the library has.
 * <p>
 * Each set carries the difficulty NIST assigns it, the certified parameters and
 * their standard deviations, and the residual measure: {@code R^2} and the
 * residual standard deviation for a linear fit, the residual sum of squares for
 * a nonlinear one. The nonlinear sets also carry the two starting points NIST
 * prescribes, a far one and a near one, because a solver that only works from
 * the near one is only half a solver.
 */
public final class StRD {

    /** A linear least squares reference set. */
    public static final class LinearSet {
        /** The name NIST gives the set. */
        public final String name;
        /** The difficulty NIST assigns: lower, average or higher. */
        public final String difficulty;
        /** Number of observations. */
        public final int observations;
        /** Number of parameters, the intercept included. */
        public final int parameters;
        /**
         * The predictor columns as the file prints them, {@code [row][column]}.
         * A set with one column and more than two parameters is a polynomial
         * fit; see {@link #isPolynomial()}.
         */
        private final double[][] predictors;
        /** The response. */
        private final double[] response;
        /** The certified parameter estimates. */
        private final double[] beta;
        /** The certified standard deviations of the estimates. */
        private final double[] standardDeviation;
        /** The certified residual standard deviation. */
        public final double residualStandardDeviation;
        /** The certified coefficient of determination. */
        public final double rSquared;

        LinearSet(String name, String difficulty, int observations, int parameters, double[][] predictors,
                double[] response, double[] beta, double[] standardDeviation, double residualStandardDeviation,
                double rSquared) {
            this.name = name;
            this.difficulty = difficulty;
            this.observations = observations;
            this.parameters = parameters;
            this.predictors = predictors;
            this.response = response;
            this.beta = beta;
            this.standardDeviation = standardDeviation;
            this.residualStandardDeviation = residualStandardDeviation;
            this.rSquared = rSquared;
        }

        /**
         * Whether the model is a polynomial in the single predictor rather than
         * a regression on several.
         *
         * @return {@code true} if one predictor column carries more than two
         *         parameters
         */
        public boolean isPolynomial() {
            return predictors[0].length == 1 && parameters > 2;
        }

        /**
         * The design matrix the model calls for, flat and column-major with the
         * intercept first, ready for {@code FlatParallelJacobiSVD}.
         *
         * @return a fresh {@code observations x parameters} matrix
         */
        public double[] design() {
            double[] a = new double[observations * parameters];
            for (int i = 0; i < observations; ++i) {
                a[i] = 1.0;
                if (isPolynomial()) {
                    double power = 1.0;
                    for (int j = 1; j < parameters; ++j) {
                        power *= predictors[i][0];
                        a[j * observations + i] = power;
                    }
                } else {
                    for (int j = 1; j < parameters; ++j) {
                        a[j * observations + i] = predictors[i][j - 1];
                    }
                }
            }
            return a;
        }

        /**
         * The response.
         *
         * @return a fresh copy
         */
        public double[] response() {
            return response.clone();
        }

        /**
         * The certified parameter estimates.
         *
         * @return a fresh copy
         */
        public double[] certifiedBeta() {
            return beta.clone();
        }

        /**
         * The certified standard deviations of the estimates.
         *
         * @return a fresh copy
         */
        public double[] certifiedStandardDeviation() {
            return standardDeviation.clone();
        }
    }

    /** A nonlinear least squares reference set. */
    public static final class NonlinearSet {
        /** The name NIST gives the set. */
        public final String name;
        /** The difficulty NIST assigns: lower, average or higher. */
        public final String difficulty;
        /** Number of observations. */
        public final int observations;
        /** Number of parameters. */
        public final int parameters;
        private final double[] x;
        private final double[] y;
        private final double[] start1;
        private final double[] start2;
        private final double[] beta;
        private final double[] standardDeviation;
        /** The certified residual sum of squares. */
        public final double residualSumOfSquares;
        /** The certified residual standard deviation. */
        public final double residualStandardDeviation;

        NonlinearSet(String name, String difficulty, int observations, int parameters, double[] x, double[] y,
                double[] start1, double[] start2, double[] beta, double[] standardDeviation,
                double residualSumOfSquares, double residualStandardDeviation) {
            this.name = name;
            this.difficulty = difficulty;
            this.observations = observations;
            this.parameters = parameters;
            this.x = x;
            this.y = y;
            this.start1 = start1;
            this.start2 = start2;
            this.beta = beta;
            this.standardDeviation = standardDeviation;
            this.residualSumOfSquares = residualSumOfSquares;
            this.residualStandardDeviation = residualStandardDeviation;
        }

        /**
         * The predictor.
         *
         * @return a fresh copy
         */
        public double[] x() {
            return x.clone();
        }

        /**
         * The response.
         *
         * @return a fresh copy
         */
        public double[] y() {
            return y.clone();
        }

        /**
         * One of the two starting points NIST prescribes.
         *
         * @param which
         *            {@code 1} for the far start, {@code 2} for the near one
         * @return a fresh copy
         */
        public double[] start(int which) {
            return which == 1 ? start1.clone() : start2.clone();
        }

        /**
         * The certified parameter estimates.
         *
         * @return a fresh copy
         */
        public double[] certifiedBeta() {
            return beta.clone();
        }

        /**
         * The certified standard deviations of the estimates.
         *
         * @return a fresh copy
         */
        public double[] certifiedStandardDeviation() {
            return standardDeviation.clone();
        }
    }

    private static final double[] Norris_Y = {
            0.1, 338.8, 118.1, 888.0, 9.2, 228.1,
            668.5, 998.5, 449.1, 778.9, 559.2, 0.3,
            0.1, 778.1, 668.8, 339.3, 448.9, 10.8,
            557.7, 228.3, 998.0, 888.8, 119.6, 0.3,
            0.6, 557.6, 339.3, 888.0, 998.5, 778.9,
            10.2, 117.6, 228.9, 668.4, 449.2, 0.2
    };

    private static final double[][] Norris_X = {
            { 0.2 },
            { 337.4 },
            { 118.2 },
            { 884.6 },
            { 10.1 },
            { 226.5 },
            { 666.3 },
            { 996.3 },
            { 448.6 },
            { 777.0 },
            { 558.2 },
            { 0.4 },
            { 0.6 },
            { 775.5 },
            { 666.9 },
            { 338.0 },
            { 447.5 },
            { 11.6 },
            { 556.0 },
            { 228.1 },
            { 995.8 },
            { 887.6 },
            { 120.2 },
            { 0.3 },
            { 0.3 },
            { 556.8 },
            { 339.1 },
            { 887.2 },
            { 999.0 },
            { 779.0 },
            { 11.1 },
            { 118.3 },
            { 229.2 },
            { 669.1 },
            { 448.9 },
            { 0.5 },
    };

    private static final double[] Norris_BETA = { -0.262323073774029, 1.00211681802045 };
    private static final double[] Norris_SD = { 0.232818234301152, 0.429796848199937E-03 };

    private static final double[] Pontius_Y = {
            .11019, .21956, .32949, .43899, .54803, .65694,
            .76562, .87487, .98292, 1.09146, 1.20001, 1.30822,
            1.41599, 1.52399, 1.63194, 1.73947, 1.84646, 1.95392,
            2.06128, 2.16844, .11052, .22018, .32939, .43886,
            .54798, .65739, .76596, .87474, .98300, 1.09150,
            1.20004, 1.30818, 1.41613, 1.52408, 1.63159, 1.73965,
            1.84696, 1.95445, 2.06177, 2.16829
    };

    private static final double[][] Pontius_X = {
            { 150000.0 },
            { 300000.0 },
            { 450000.0 },
            { 600000.0 },
            { 750000.0 },
            { 900000.0 },
            { 1050000.0 },
            { 1200000.0 },
            { 1350000.0 },
            { 1500000.0 },
            { 1650000.0 },
            { 1800000.0 },
            { 1950000.0 },
            { 2100000.0 },
            { 2250000.0 },
            { 2400000.0 },
            { 2550000.0 },
            { 2700000.0 },
            { 2850000.0 },
            { 3000000.0 },
            { 150000.0 },
            { 300000.0 },
            { 450000.0 },
            { 600000.0 },
            { 750000.0 },
            { 900000.0 },
            { 1050000.0 },
            { 1200000.0 },
            { 1350000.0 },
            { 1500000.0 },
            { 1650000.0 },
            { 1800000.0 },
            { 1950000.0 },
            { 2100000.0 },
            { 2250000.0 },
            { 2400000.0 },
            { 2550000.0 },
            { 2700000.0 },
            { 2850000.0 },
            { 3000000.0 },
    };

    private static final double[] Pontius_BETA = { 0.673565789473684E-03, 0.732059160401003E-06, -0.316081871345029E-14 };
    private static final double[] Pontius_SD = { 0.107938612033077E-03, 0.157817399981659E-09, 0.486652849992036E-16 };

    private static final double[] Longley_Y = {
            60323.0, 61122.0, 60171.0, 61187.0, 63221.0, 63639.0,
            64989.0, 63761.0, 66019.0, 67857.0, 68169.0, 66513.0,
            68655.0, 69564.0, 69331.0, 70551.0
    };

    private static final double[][] Longley_X = {
            { 83.0, 234289.0, 2356.0, 1590.0, 107608.0, 1947.0 },
            { 88.5, 259426.0, 2325.0, 1456.0, 108632.0, 1948.0 },
            { 88.2, 258054.0, 3682.0, 1616.0, 109773.0, 1949.0 },
            { 89.5, 284599.0, 3351.0, 1650.0, 110929.0, 1950.0 },
            { 96.2, 328975.0, 2099.0, 3099.0, 112075.0, 1951.0 },
            { 98.1, 346999.0, 1932.0, 3594.0, 113270.0, 1952.0 },
            { 99.0, 365385.0, 1870.0, 3547.0, 115094.0, 1953.0 },
            { 100.0, 363112.0, 3578.0, 3350.0, 116219.0, 1954.0 },
            { 101.2, 397469.0, 2904.0, 3048.0, 117388.0, 1955.0 },
            { 104.6, 419180.0, 2822.0, 2857.0, 118734.0, 1956.0 },
            { 108.4, 442769.0, 2936.0, 2798.0, 120445.0, 1957.0 },
            { 110.8, 444546.0, 4681.0, 2637.0, 121950.0, 1958.0 },
            { 112.6, 482704.0, 3813.0, 2552.0, 123366.0, 1959.0 },
            { 114.2, 502601.0, 3931.0, 2514.0, 125368.0, 1960.0 },
            { 115.7, 518173.0, 4806.0, 2572.0, 127852.0, 1961.0 },
            { 116.9, 554894.0, 4007.0, 2827.0, 130081.0, 1962.0 },
    };

    private static final double[] Longley_BETA = { -3482258.63459582, 15.0618722713733, -0.358191792925910E-01, -2.02022980381683, -1.03322686717359, -0.511041056535807E-01, 1829.15146461355 };
    private static final double[] Longley_SD = { 890420.383607373, 84.9149257747669, 0.334910077722432E-01, 0.488399681651699, 0.214274163161675, 0.226073200069370, 455.478499142212 };

    private static final double[] Wampler1_Y = {
            1.0, 6.0, 63.0, 364.0, 1365.0, 3906.0,
            9331.0, 19608.0, 37449.0, 66430.0, 111111.0, 177156.0,
            271453.0, 402234.0, 579195.0, 813616.0, 1118481.0, 1508598.0,
            2000719.0, 2613660.0, 3368421.0
    };

    private static final double[][] Wampler1_X = {
            { 0.0 },
            { 1.0 },
            { 2.0 },
            { 3.0 },
            { 4.0 },
            { 5.0 },
            { 6.0 },
            { 7.0 },
            { 8.0 },
            { 9.0 },
            { 10.0 },
            { 11.0 },
            { 12.0 },
            { 13.0 },
            { 14.0 },
            { 15.0 },
            { 16.0 },
            { 17.0 },
            { 18.0 },
            { 19.0 },
            { 20.0 },
    };

    private static final double[] Wampler1_BETA = { 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000 };
    private static final double[] Wampler1_SD = { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 };

    private static final double[] Wampler2_Y = {
            1.00000, 1.11111, 1.24992, 1.42753, 1.65984, 1.96875,
            2.38336, 2.94117, 3.68928, 4.68559, 6.00000, 7.71561,
            9.92992, 12.75603, 16.32384, 20.78125, 26.29536, 33.05367,
            41.26528, 51.16209, 63.00000
    };

    private static final double[][] Wampler2_X = {
            { 0.0 },
            { 1.0 },
            { 2.0 },
            { 3.0 },
            { 4.0 },
            { 5.0 },
            { 6.0 },
            { 7.0 },
            { 8.0 },
            { 9.0 },
            { 10.0 },
            { 11.0 },
            { 12.0 },
            { 13.0 },
            { 14.0 },
            { 15.0 },
            { 16.0 },
            { 17.0 },
            { 18.0 },
            { 19.0 },
            { 20.0 },
    };

    private static final double[] Wampler2_BETA = { 1.00000000000000, 0.100000000000000, 0.100000000000000E-01, 0.100000000000000E-02, 0.100000000000000E-03, 0.100000000000000E-04 };
    private static final double[] Wampler2_SD = { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 };

    private static final double[] Wampler3_Y = {
            760., -2042., 2111., -1684., 3888., 1858.,
            11379., 17560., 39287., 64382., 113159., 175108.,
            273291., 400186., 581243., 811568., 1121004., 1506550.,
            2002767., 2611612., 3369180.
    };

    private static final double[][] Wampler3_X = {
            { 0.0 },
            { 1.0 },
            { 2.0 },
            { 3.0 },
            { 4.0 },
            { 5.0 },
            { 6.0 },
            { 7.0 },
            { 8.0 },
            { 9.0 },
            { 10.0 },
            { 11.0 },
            { 12.0 },
            { 13.0 },
            { 14.0 },
            { 15.0 },
            { 16.0 },
            { 17.0 },
            { 18.0 },
            { 19.0 },
            { 20.0 },
    };

    private static final double[] Wampler3_BETA = { 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000 };
    private static final double[] Wampler3_SD = { 2152.32624678170, 2363.55173469681, 779.343524331583, 101.475507550350, 5.64566512170752, 0.112324854679312 };

    private static final double[] Wampler4_Y = {
            75901.0, -204794.0, 204863.0, -204436.0, 253665.0, -200894.0,
            214131.0, -185192.0, 221249.0, -138370.0, 315911.0, -27644.0,
            455253.0, 197434.0, 783995.0, 608816.0, 1370781.0, 1303798.0,
            2205519.0, 2408860.0, 3444321.0
    };

    private static final double[][] Wampler4_X = {
            { 0.0 },
            { 1.0 },
            { 2.0 },
            { 3.0 },
            { 4.0 },
            { 5.0 },
            { 6.0 },
            { 7.0 },
            { 8.0 },
            { 9.0 },
            { 10.0 },
            { 11.0 },
            { 12.0 },
            { 13.0 },
            { 14.0 },
            { 15.0 },
            { 16.0 },
            { 17.0 },
            { 18.0 },
            { 19.0 },
            { 20.0 },
    };

    private static final double[] Wampler4_BETA = { 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000 };
    private static final double[] Wampler4_SD = { 215232.624678170, 236355.173469681, 77934.3524331583, 10147.5507550350, 564.566512170752, 11.2324854679312 };

    private static final double[] Wampler5_Y = {
            7590001.0, -20479994.0, 20480063.0, -20479636.0, 25231365.0, -20476094.0,
            20489331.0, -20460392.0, 18417449.0, -20413570.0, 20591111.0, -20302844.0,
            18651453.0, -20077766.0, 21059195.0, -19666384.0, 26348481.0, -18971402.0,
            22480719.0, -17866340.0, 10958421.0
    };

    private static final double[][] Wampler5_X = {
            { 0.0 },
            { 1.0 },
            { 2.0 },
            { 3.0 },
            { 4.0 },
            { 5.0 },
            { 6.0 },
            { 7.0 },
            { 8.0 },
            { 9.0 },
            { 10.0 },
            { 11.0 },
            { 12.0 },
            { 13.0 },
            { 14.0 },
            { 15.0 },
            { 16.0 },
            { 17.0 },
            { 18.0 },
            { 19.0 },
            { 20.0 },
    };

    private static final double[] Wampler5_BETA = { 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000 };
    private static final double[] Wampler5_SD = { 21523262.4678170, 23635517.3469681, 7793435.24331583, 1014755.07550350, 56456.6512170752, 1123.24854679312 };

    private static final double[] Filip_Y = {
            0.8116, 0.9072, 0.9052, 0.9039, 0.8053, 0.8377,
            0.8667, 0.8809, 0.7975, 0.8162, 0.8515, 0.8766,
            0.8885, 0.8859, 0.8959, 0.8913, 0.8959, 0.8971,
            0.9021, 0.909, 0.9139, 0.9199, 0.8692, 0.8872,
            0.89, 0.891, 0.8977, 0.9035, 0.9078, 0.7675,
            0.7705, 0.7713, 0.7736, 0.7775, 0.7841, 0.7971,
            0.8329, 0.8641, 0.8804, 0.7668, 0.7633, 0.7678,
            0.7697, 0.77, 0.7749, 0.7796, 0.7897, 0.8131,
            0.8498, 0.8741, 0.8061, 0.846, 0.8751, 0.8856,
            0.8919, 0.8934, 0.894, 0.8957, 0.9047, 0.9129,
            0.9209, 0.9219, 0.7739, 0.7681, 0.7665, 0.7703,
            0.7702, 0.7761, 0.7809, 0.7961, 0.8253, 0.8602,
            0.8809, 0.8301, 0.8664, 0.8834, 0.8898, 0.8964,
            0.8963, 0.9074, 0.9119, 0.9228
    };

    private static final double[][] Filip_X = {
            { -6.860120914 },
            { -4.324130045 },
            { -4.358625055 },
            { -4.358426747 },
            { -6.955852379 },
            { -6.661145254 },
            { -6.355462942 },
            { -6.118102026 },
            { -7.115148017 },
            { -6.815308569 },
            { -6.519993057 },
            { -6.204119983 },
            { -5.853871964 },
            { -6.109523091 },
            { -5.79832982 },
            { -5.482672118 },
            { -5.171791386 },
            { -4.851705903 },
            { -4.517126416 },
            { -4.143573228 },
            { -3.709075441 },
            { -3.499489089 },
            { -6.300769497 },
            { -5.953504836 },
            { -5.642065153 },
            { -5.031376979 },
            { -4.680685696 },
            { -4.329846955 },
            { -3.928486195 },
            { -8.56735134 },
            { -8.363211311 },
            { -8.107682739 },
            { -7.823908741 },
            { -7.522878745 },
            { -7.218819279 },
            { -6.920818754 },
            { -6.628932138 },
            { -6.323946875 },
            { -5.991399828 },
            { -8.781464495 },
            { -8.663140179 },
            { -8.473531488 },
            { -8.247337057 },
            { -7.971428747 },
            { -7.676129393 },
            { -7.352812702 },
            { -7.072065318 },
            { -6.774174009 },
            { -6.478861916 },
            { -6.159517513 },
            { -6.835647144 },
            { -6.53165267 },
            { -6.224098421 },
            { -5.910094889 },
            { -5.598599459 },
            { -5.290645224 },
            { -4.974284616 },
            { -4.64454848 },
            { -4.290560426 },
            { -3.885055584 },
            { -3.408378962 },
            { -3.13200249 },
            { -8.726767166 },
            { -8.66695597 },
            { -8.511026475 },
            { -8.165388579 },
            { -7.886056648 },
            { -7.588043762 },
            { -7.283412422 },
            { -6.995678626 },
            { -6.691862621 },
            { -6.392544977 },
            { -6.067374056 },
            { -6.684029655 },
            { -6.378719832 },
            { -6.065855188 },
            { -5.752272167 },
            { -5.132414673 },
            { -4.811352704 },
            { -4.098269308 },
            { -3.66174277 },
            { -3.2644011 },
    };

    private static final double[] Filip_BETA = { -1467.48961422980, -2772.17959193342, -2316.37108160893, -1127.97394098372, -354.478233703349, -75.1242017393757, -10.8753180355343, -1.06221498588947, -0.670191154593408E-01, -0.246781078275479E-02, -0.402962525080404E-04 };
    private static final double[] Filip_SD = { 298.084530995537, 559.779865474950, 466.477572127796, 227.204274477751, 71.6478660875927, 15.2897178747400, 2.23691159816033, 0.221624321934227, 0.142363763154724E-01, 0.535617408889821E-03, 0.896632837373868E-05 };

    private static final double[] Misra1a_Y = {
            10.07E0, 14.73E0, 17.94E0, 23.93E0, 29.61E0, 35.18E0,
            40.02E0, 44.82E0, 50.76E0, 55.05E0, 61.01E0, 66.40E0,
            75.47E0, 81.78E0
    };

    private static final double[] Misra1a_X = {
            77.6E0, 114.9E0, 141.1E0, 190.8E0, 239.9E0, 289.0E0,
            332.8E0, 378.4E0, 434.8E0, 477.3E0, 536.8E0, 593.1E0,
            689.1E0, 760.0E0
    };

    private static final double[] Misra1a_START1 = { 500.0, 0.0001 };
    private static final double[] Misra1a_START2 = { 250.0, 0.0005 };
    private static final double[] Misra1a_BETA = { 2.3894212918E+02, 5.5015643181E-04 };
    private static final double[] Misra1a_SD = { 2.7070075241E+00, 7.2668688436E-06 };

    private static final double[] Chwirut1_Y = {
            92.9000E0, 78.7000E0, 64.2000E0, 64.9000E0, 57.1000E0, 43.3000E0,
            31.1000E0, 23.6000E0, 31.0500E0, 23.7750E0, 17.7375E0, 13.8000E0,
            11.5875E0, 9.4125E0, 7.7250E0, 7.3500E0, 8.0250E0, 90.6000E0,
            76.9000E0, 71.6000E0, 63.6000E0, 54.0000E0, 39.2000E0, 29.3000E0,
            21.4000E0, 29.1750E0, 22.1250E0, 17.5125E0, 14.2500E0, 9.4500E0,
            9.1500E0, 7.9125E0, 8.4750E0, 6.1125E0, 80.0000E0, 79.0000E0,
            63.8000E0, 57.2000E0, 53.2000E0, 42.5000E0, 26.8000E0, 20.4000E0,
            26.8500E0, 21.0000E0, 16.4625E0, 12.5250E0, 10.5375E0, 8.5875E0,
            7.1250E0, 6.1125E0, 5.9625E0, 74.1000E0, 67.3000E0, 60.8000E0,
            55.5000E0, 50.3000E0, 41.0000E0, 29.4000E0, 20.4000E0, 29.3625E0,
            21.1500E0, 16.7625E0, 13.2000E0, 10.8750E0, 8.1750E0, 7.3500E0,
            5.9625E0, 5.6250E0, 81.5000E0, 62.4000E0, 32.5000E0, 12.4100E0,
            13.1200E0, 15.5600E0, 5.6300E0, 78.0000E0, 59.9000E0, 33.2000E0,
            13.8400E0, 12.7500E0, 14.6200E0, 3.9400E0, 76.8000E0, 61.0000E0,
            32.9000E0, 13.8700E0, 11.8100E0, 13.3100E0, 5.4400E0, 78.0000E0,
            63.5000E0, 33.8000E0, 12.5600E0, 5.6300E0, 12.7500E0, 13.1200E0,
            5.4400E0, 76.8000E0, 60.0000E0, 47.8000E0, 32.0000E0, 22.2000E0,
            22.5700E0, 18.8200E0, 13.9500E0, 11.2500E0, 9.0000E0, 6.6700E0,
            75.8000E0, 62.0000E0, 48.8000E0, 35.2000E0, 20.0000E0, 20.3200E0,
            19.3100E0, 12.7500E0, 10.4200E0, 7.3100E0, 7.4200E0, 70.5000E0,
            59.5000E0, 48.5000E0, 35.8000E0, 21.0000E0, 21.6700E0, 21.0000E0,
            15.6400E0, 8.1700E0, 8.5500E0, 10.1200E0, 78.0000E0, 66.0000E0,
            62.0000E0, 58.0000E0, 47.7000E0, 37.8000E0, 20.2000E0, 21.0700E0,
            13.8700E0, 9.6700E0, 7.7600E0, 5.4400E0, 4.8700E0, 4.0100E0,
            3.7500E0, 24.1900E0, 25.7600E0, 18.0700E0, 11.8100E0, 12.0700E0,
            16.1200E0, 70.8000E0, 54.7000E0, 48.0000E0, 39.8000E0, 29.8000E0,
            23.7000E0, 29.6200E0, 23.8100E0, 17.7000E0, 11.5500E0, 12.0700E0,
            8.7400E0, 80.7000E0, 61.3000E0, 47.5000E0, 29.0000E0, 24.0000E0,
            17.7000E0, 24.5600E0, 18.6700E0, 16.2400E0, 8.7400E0, 7.8700E0,
            8.5100E0, 66.7000E0, 59.2000E0, 40.8000E0, 30.7000E0, 25.7000E0,
            16.3000E0, 25.9900E0, 16.9500E0, 13.3500E0, 8.6200E0, 7.2000E0,
            6.6400E0, 13.6900E0, 81.0000E0, 64.5000E0, 35.5000E0, 13.3100E0,
            4.8700E0, 12.9400E0, 5.0600E0, 15.1900E0, 14.6200E0, 15.6400E0,
            25.5000E0, 25.9500E0, 81.7000E0, 61.6000E0, 29.8000E0, 29.8100E0,
            17.1700E0, 10.3900E0, 28.4000E0, 28.6900E0, 81.3000E0, 60.9000E0,
            16.6500E0, 10.0500E0, 28.9000E0, 28.9500E0
    };

    private static final double[] Chwirut1_X = {
            0.5000E0, 0.6250E0, 0.7500E0, 0.8750E0, 1.0000E0, 1.2500E0,
            1.7500E0, 2.2500E0, 1.7500E0, 2.2500E0, 2.7500E0, 3.2500E0,
            3.7500E0, 4.2500E0, 4.7500E0, 5.2500E0, 5.7500E0, 0.5000E0,
            0.6250E0, 0.7500E0, 0.8750E0, 1.0000E0, 1.2500E0, 1.7500E0,
            2.2500E0, 1.7500E0, 2.2500E0, 2.7500E0, 3.2500E0, 3.7500E0,
            4.2500E0, 4.7500E0, 5.2500E0, 5.7500E0, 0.5000E0, 0.6250E0,
            0.7500E0, 0.8750E0, 1.0000E0, 1.2500E0, 1.7500E0, 2.2500E0,
            1.7500E0, 2.2500E0, 2.7500E0, 3.2500E0, 3.7500E0, 4.2500E0,
            4.7500E0, 5.2500E0, 5.7500E0, 0.5000E0, 0.6250E0, 0.7500E0,
            0.8750E0, 1.0000E0, 1.2500E0, 1.7500E0, 2.2500E0, 1.7500E0,
            2.2500E0, 2.7500E0, 3.2500E0, 3.7500E0, 4.2500E0, 4.7500E0,
            5.2500E0, 5.7500E0, .5000E0, .7500E0, 1.5000E0, 3.0000E0,
            3.0000E0, 3.0000E0, 6.0000E0, .5000E0, .7500E0, 1.5000E0,
            3.0000E0, 3.0000E0, 3.0000E0, 6.0000E0, .5000E0, .7500E0,
            1.5000E0, 3.0000E0, 3.0000E0, 3.0000E0, 6.0000E0, .5000E0,
            .7500E0, 1.5000E0, 3.0000E0, 6.0000E0, 3.0000E0, 3.0000E0,
            6.0000E0, .5000E0, .7500E0, 1.0000E0, 1.5000E0, 2.0000E0,
            2.0000E0, 2.5000E0, 3.0000E0, 4.0000E0, 5.0000E0, 6.0000E0,
            .5000E0, .7500E0, 1.0000E0, 1.5000E0, 2.0000E0, 2.0000E0,
            2.5000E0, 3.0000E0, 4.0000E0, 5.0000E0, 6.0000E0, .5000E0,
            .7500E0, 1.0000E0, 1.5000E0, 2.0000E0, 2.0000E0, 2.5000E0,
            3.0000E0, 4.0000E0, 5.0000E0, 6.0000E0, .5000E0, .6250E0,
            .7500E0, .8750E0, 1.0000E0, 1.2500E0, 2.2500E0, 2.2500E0,
            2.7500E0, 3.2500E0, 3.7500E0, 4.2500E0, 4.7500E0, 5.2500E0,
            5.7500E0, 3.0000E0, 3.0000E0, 3.0000E0, 3.0000E0, 3.0000E0,
            3.0000E0, .5000E0, .7500E0, 1.0000E0, 1.5000E0, 2.0000E0,
            2.5000E0, 2.0000E0, 2.5000E0, 3.0000E0, 4.0000E0, 5.0000E0,
            6.0000E0, .5000E0, .7500E0, 1.0000E0, 1.5000E0, 2.0000E0,
            2.5000E0, 2.0000E0, 2.5000E0, 3.0000E0, 4.0000E0, 5.0000E0,
            6.0000E0, .5000E0, .7500E0, 1.0000E0, 1.5000E0, 2.0000E0,
            2.5000E0, 2.0000E0, 2.5000E0, 3.0000E0, 4.0000E0, 5.0000E0,
            6.0000E0, 3.0000E0, .5000E0, .7500E0, 1.5000E0, 3.0000E0,
            6.0000E0, 3.0000E0, 6.0000E0, 3.0000E0, 3.0000E0, 3.0000E0,
            1.7500E0, 1.7500E0, .5000E0, .7500E0, 1.7500E0, 1.7500E0,
            2.7500E0, 3.7500E0, 1.7500E0, 1.7500E0, .5000E0, .7500E0,
            2.7500E0, 3.7500E0, 1.7500E0, 1.7500E0
    };

    private static final double[] Chwirut1_START1 = { 0.1, 0.01, 0.02 };
    private static final double[] Chwirut1_START2 = { 0.15, 0.008, 0.010 };
    private static final double[] Chwirut1_BETA = { 1.9027818370E-01, 6.1314004477E-03, 1.0530908399E-02 };
    private static final double[] Chwirut1_SD = { 2.1938557035E-02, 3.4500025051E-04, 7.9281847748E-04 };

    private static final double[] Thurber_Y = {
            80.574E0, 84.248E0, 87.264E0, 87.195E0, 89.076E0, 89.608E0,
            89.868E0, 90.101E0, 92.405E0, 95.854E0, 100.696E0, 101.060E0,
            401.672E0, 390.724E0, 567.534E0, 635.316E0, 733.054E0, 759.087E0,
            894.206E0, 990.785E0, 1090.109E0, 1080.914E0, 1122.643E0, 1178.351E0,
            1260.531E0, 1273.514E0, 1288.339E0, 1327.543E0, 1353.863E0, 1414.509E0,
            1425.208E0, 1421.384E0, 1442.962E0, 1464.350E0, 1468.705E0, 1447.894E0,
            1457.628E0
    };

    private static final double[] Thurber_X = {
            -3.067E0, -2.981E0, -2.921E0, -2.912E0, -2.840E0, -2.797E0,
            -2.702E0, -2.699E0, -2.633E0, -2.481E0, -2.363E0, -2.322E0,
            -1.501E0, -1.460E0, -1.274E0, -1.212E0, -1.100E0, -1.046E0,
            -0.915E0, -0.714E0, -0.566E0, -0.545E0, -0.400E0, -0.309E0,
            -0.109E0, -0.103E0, 0.010E0, 0.119E0, 0.377E0, 0.790E0,
            0.963E0, 1.006E0, 1.115E0, 1.572E0, 1.841E0, 2.047E0,
            2.200E0
    };

    private static final double[] Thurber_START1 = { 1000.0, 1000.0, 400.0, 40.0, 0.7, 0.3, 0.03 };
    private static final double[] Thurber_START2 = { 1300.0, 1500.0, 500.0, 75.0, 1.0, 0.4, 0.05 };
    private static final double[] Thurber_BETA = { 1.2881396800E+03, 1.4910792535E+03, 5.8323836877E+02, 7.5416644291E+01, 9.6629502864E-01, 3.9797285797E-01, 4.9727297349E-02 };
    private static final double[] Thurber_SD = { 4.6647963344E+00, 3.9571156086E+01, 2.8698696102E+01, 5.5675370270E+00, 3.1333340687E-02, 1.4984928198E-02, 6.5842344623E-03 };

    private static final double[] MGH09_Y = {
            1.957000E-01, 1.947000E-01, 1.735000E-01, 1.600000E-01, 8.440000E-02, 6.270000E-02,
            4.560000E-02, 3.420000E-02, 3.230000E-02, 2.350000E-02, 2.460000E-02
    };

    private static final double[] MGH09_X = {
            4.000000E+00, 2.000000E+00, 1.000000E+00, 5.000000E-01, 2.500000E-01, 1.670000E-01,
            1.250000E-01, 1.000000E-01, 8.330000E-02, 7.140000E-02, 6.250000E-02
    };

    private static final double[] MGH09_START1 = { 25.0, 39.0, 41.5, 39.0 };
    private static final double[] MGH09_START2 = { 0.25, 0.39, 0.415, 0.39 };
    private static final double[] MGH09_BETA = { 1.9280693458E-01, 1.9128232873E-01, 1.2305650693E-01, 1.3606233068E-01 };
    private static final double[] MGH09_SD = { 1.1435312227E-02, 1.9633220911E-01, 8.0842031232E-02, 9.0025542308E-02 };

    private static final double[] MGH10_Y = {
            3.478000E+04, 2.861000E+04, 2.365000E+04, 1.963000E+04, 1.637000E+04, 1.372000E+04,
            1.154000E+04, 9.744000E+03, 8.261000E+03, 7.030000E+03, 6.005000E+03, 5.147000E+03,
            4.427000E+03, 3.820000E+03, 3.307000E+03, 2.872000E+03
    };

    private static final double[] MGH10_X = {
            5.000000E+01, 5.500000E+01, 6.000000E+01, 6.500000E+01, 7.000000E+01, 7.500000E+01,
            8.000000E+01, 8.500000E+01, 9.000000E+01, 9.500000E+01, 1.000000E+02, 1.050000E+02,
            1.100000E+02, 1.150000E+02, 1.200000E+02, 1.250000E+02
    };

    private static final double[] MGH10_START1 = { 2.0, 400000.0, 25000.0 };
    private static final double[] MGH10_START2 = { 0.02, 4000.0, 250.0 };
    private static final double[] MGH10_BETA = { 5.6096364710E-03, 6.1813463463E+03, 3.4522363462E+02 };
    private static final double[] MGH10_SD = { 1.5687892471E-04, 2.3309021107E+01, 7.8486103508E-01 };

    /**
     * The nine linear reference sets, easiest first.
     *
     * @return a fresh array of the sets
     */
    public static LinearSet[] linear() {
        return new LinearSet[] {
                new LinearSet("Norris", "lower", 36, 2, Norris_X, Norris_Y, Norris_BETA, Norris_SD, 0.884796396144373, 0.999993745883712),
                new LinearSet("Pontius", "lower", 40, 3, Pontius_X, Pontius_Y, Pontius_BETA, Pontius_SD, 0.205177424076185E-03, 0.999999900178537),
                new LinearSet("Longley", "higher", 16, 7, Longley_X, Longley_Y, Longley_BETA, Longley_SD, 304.854073561965, 0.995479004577296),
                new LinearSet("Wampler1", "higher", 21, 6, Wampler1_X, Wampler1_Y, Wampler1_BETA, Wampler1_SD, 0.000000000000000, 1.00000000000000),
                new LinearSet("Wampler2", "higher", 21, 6, Wampler2_X, Wampler2_Y, Wampler2_BETA, Wampler2_SD, 0.000000000000000, 1.00000000000000),
                new LinearSet("Wampler3", "higher", 21, 6, Wampler3_X, Wampler3_Y, Wampler3_BETA, Wampler3_SD, 2360.14502379268, 0.999995559025820),
                new LinearSet("Wampler4", "higher", 21, 6, Wampler4_X, Wampler4_Y, Wampler4_BETA, Wampler4_SD, 236014.502379268, 0.957478440825662),
                new LinearSet("Wampler5", "higher", 21, 6, Wampler5_X, Wampler5_Y, Wampler5_BETA, Wampler5_SD, 23601450.2379268, 0.224668921574940E-02),
                new LinearSet("Filip", "higher", 82, 11, Filip_X, Filip_Y, Filip_BETA, Filip_SD, 0.334801051324544E-02, 0.996727416185620),
        };
    }

    /**
     * The five nonlinear reference sets, easiest first.
     *
     * @return a fresh array of the sets
     */
    public static NonlinearSet[] nonlinear() {
        return new NonlinearSet[] {
                new NonlinearSet("Misra1a", "lower", 14, 2, Misra1a_X, Misra1a_Y, Misra1a_START1, Misra1a_START2, Misra1a_BETA, Misra1a_SD, 1.2455138894E-01, 1.0187876330E-01),
                new NonlinearSet("Chwirut1", "lower", 214, 3, Chwirut1_X, Chwirut1_Y, Chwirut1_START1, Chwirut1_START2, Chwirut1_BETA, Chwirut1_SD, 2.3844771393E+03, 3.3616721320E+00),
                new NonlinearSet("Thurber", "higher", 37, 7, Thurber_X, Thurber_Y, Thurber_START1, Thurber_START2, Thurber_BETA, Thurber_SD, 5.6427082397E+03, 1.3714600784E+01),
                new NonlinearSet("MGH09", "higher", 11, 4, MGH09_X, MGH09_Y, MGH09_START1, MGH09_START2, MGH09_BETA, MGH09_SD, 3.0750560385E-04, 6.6279236551E-03),
                new NonlinearSet("MGH10", "higher", 16, 3, MGH10_X, MGH10_Y, MGH10_START1, MGH10_START2, MGH10_BETA, MGH10_SD, 8.7945855171E+01, 2.6009740065E+00),
        };
    }

    private StRD() {
        throw new AssertionError();
    }
}
