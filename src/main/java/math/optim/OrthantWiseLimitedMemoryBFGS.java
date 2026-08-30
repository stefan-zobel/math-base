/*
 * Copyright (C) 2005 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

import java.util.logging.Logger;

import math.linalg.VectorOps;


/**
 * Implementation of orthant-wise limited memory quasi Newton method for
 * optimizing convex L1-regularized objectives. See:
 * "Scalable training of l1-regularized log-linear models" by Galen Andrew and
 * Jianfeng Gao in ICML 2007 for details. This code is an adaptation of the
 * freely-available C++ code on Galen's webpage.
 * 
 * <p>
 * Like {@link LimitedMemoryBFGS} this class maximizes, and it adds the L1
 * penalty itself: the objective handed to it is the smooth part alone, as a
 * log-likelihood. {@link #getTermination()} says which rule stopped the
 * search and {@link #getGradientNorm()} gives the norm of the <em>pseudo</em>
 * gradient it fired on, which is the quantity that vanishes at an L1 solution;
 * the ordinary gradient does not, since on the active set it equals the L1
 * weight times the sign of the coordinate.
 *
 * @author Kedar Bellare
 */
public final class OrthantWiseLimitedMemoryBFGS implements Optimizer {
    private static final Logger logger = Logger
            .getLogger(OrthantWiseLimitedMemoryBFGS.class.getName());

    /** Iteration budget if none is given. */
    private static final int DEFAULT_MAX_ITERATIONS = 1000;
    /** Relative change of the objective value that stops the search. */
    private static final double DEFAULT_TOLERANCE = 1.0e-4;
    /** Pseudo gradient norm that stops the search. */
    private static final double DEFAULT_GRADIENT_TOLERANCE = 1.0e-3;
    /** Number of corrections kept for the inverse Hessian. */
    private static final int DEFAULT_M = 4;

    private Termination termination = Termination.NOT_STARTED;
    private double gradientNorm = Double.NaN;
    private Optimizable.ByGradientValue optimizable;
    // name of optimizable for value output
//    private String optName;
    private final int maxIterations;
    private final double tolerance;
    private final double gradientTolerance;
    private final double eps = 1.0e-5;
    private double l1Weight;

    // The number of corrections used in BFGS update
    // ideally 3 <= m <= 7. Larger m means more cpu time, memory.
    private final int m;

    // Upper bound on the backtracking steps of one line search. Halving from
    // 1.0 this often is already far below anything a double can express as a
    // useful step.
    private static final int MAX_BACKTRACKS = 100;

    // State of optimizer search
    // oldValue = value before line search, value = value after line search
    private double oldValue, value, yDotY;
    // grad = gradient
    private double[] grad, oldGrad, direction, steepestDescentDirection,
            parameters, oldParameters;
    // scratch for the parameter and gradient differences of one BFGS update
    private double[] sBuf, yBuf;
    // s = list of m previous "difference in parameters" values
    // y = list of m previous "difference in grad" values
    private SupersedingDoubleArrayQueue s, y;
    // rho = intermediate calculation
    private SupersedingDoubleQueue rhos;
    private double[] alphas;
    private int iterations;

    /**
     * Creates an optimizer without an L1 penalty and with the default
     * tolerances.
     *
     * @param function
     *            the function to optimize
     */
    public OrthantWiseLimitedMemoryBFGS(Optimizable.ByGradientValue function) {
        this(function, 0.0);
    }

    /**
     * Creates an optimizer with the default tolerances.
     *
     * @param function
     *            the function to optimize
     * @param l1wt
     *            weight of the L1 penalty, {@code 0} or greater
     */
    public OrthantWiseLimitedMemoryBFGS(Optimizable.ByGradientValue function,
            double l1wt) {
        this(function, l1wt, DEFAULT_MAX_ITERATIONS, DEFAULT_TOLERANCE,
                DEFAULT_GRADIENT_TOLERANCE, DEFAULT_M);
    }

    /**
     * Creates an optimizer with explicit stopping rules.
     *
     * @param function
     *            the function to optimize
     * @param l1wt
     *            weight of the L1 penalty, {@code 0} or greater
     * @param maxIterations
     *            iteration budget, {@code 1} or greater; exhausting it stops
     *            the search without reporting convergence
     * @param tolerance
     *            relative change of the objective value below which the search
     *            stops, greater than {@code 0}
     * @param gradientTolerance
     *            norm of the pseudo gradient at or below which the search
     *            stops, {@code 0} or greater
     * @param m
     *            number of corrections kept for the inverse Hessian, between
     *            {@code 1} and {@code 100}; ideally between {@code 3} and
     *            {@code 7}
     * @since 1.5.2
     */
    public OrthantWiseLimitedMemoryBFGS(Optimizable.ByGradientValue function,
            double l1wt, int maxIterations, double tolerance,
            double gradientTolerance, int m) {
        if (function == null) {
            throw new IllegalArgumentException("function is null");
        }
        if (!(l1wt >= 0.0) || Double.isInfinite(l1wt)) {
            // a negative weight would be swallowed by the l1Weight > 0 guards
            // below and silently degrade to an unregularized run
            throw new IllegalArgumentException(
                    "l1wt must be finite and non-negative : " + l1wt);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException(
                    "maxIterations must be >= 1 : " + maxIterations);
        }
        if (!(tolerance > 0.0) || Double.isInfinite(tolerance)) {
            throw new IllegalArgumentException(
                    "tolerance must be finite and positive : " + tolerance);
        }
        if (!(gradientTolerance >= 0.0) || Double.isInfinite(gradientTolerance)) {
            throw new IllegalArgumentException(
                    "gradientTolerance must be finite and non-negative : "
                            + gradientTolerance);
        }
        if (m < 1 || m > 100) {
            throw new IllegalArgumentException(
                    "m must be between 1 and 100 : " + m);
        }
        this.maxIterations = maxIterations;
        this.tolerance = tolerance;
        this.gradientTolerance = gradientTolerance;
        this.m = m;
        this.optimizable = function;
        this.l1Weight = l1wt;
//        String parts[] = optimizable.getClass().getName().split("\\.");
//        this.optName = parts[parts.length - 1];

        // initialize optimizer state
        iterations = 0;
        s = new SupersedingDoubleArrayQueue(m);
        y = new SupersedingDoubleArrayQueue(m);
        rhos = new SupersedingDoubleQueue(m);
        alphas = new double[m];
        VectorOps.setAll(alphas, 0.0);
        yDotY = 0;

        int numParameters = optimizable.getNumParameters();

        // get initial parameters
        parameters = new double[numParameters];
        optimizable.getParameters(parameters);

        // get initial value
        value = evalL1();

        // get initial gradient
        grad = new double[numParameters];
        evalGradient();

        // initialize direction
        direction = new double[numParameters];
        steepestDescentDirection = new double[numParameters];

        // initialize backups
        oldParameters = new double[numParameters];
        oldGrad = new double[numParameters];

        // initialize the update scratch buffers
        sBuf = new double[numParameters];
        yBuf = new double[numParameters];
    }

    @Override
    public Optimizable getOptimizable() {
        return optimizable;
    }

    @Override
    public boolean isConverged() {
        return termination.isConvergence();
    }

    /**
     * Why the last search stopped.
     *
     * @return the outcome of the last run, or {@link Termination#NOT_STARTED}
     *         if {@code optimize} was never called
     * @since 1.5.2
     */
    public Termination getTermination() {
        return termination;
    }

    /**
     * Euclidean norm of the <em>pseudo</em> gradient where the last search
     * stopped, that is of the steepest descent direction this method builds
     * from the gradient and the L1 weight. It is the quantity whose vanishing
     * is the optimality condition of the penalized problem.
     *
     * @return the pseudo gradient norm at the exit, or {@code NaN} if
     *         {@code optimize} was never called
     * @since 1.5.2
     */
    public double getGradientNorm() {
        return gradientNorm;
    }

    public int getIteration() {
        return iterations;
    }

    @Override
    public boolean optimize() {
        return optimize(Integer.MAX_VALUE);
    }

    @Override
    public boolean optimize(int numIterations) {
        logger.fine("Entering OWL-BFGS.optimize(). L1 weight=" + l1Weight
                + " Initial Value=" + value);

        for (int iter = 0; iter < numIterations; iter++) {
            // create descent direction
            makeSteepestDescDir();

            // The steepest descent direction is the negative pseudo gradient,
            // and that is the quantity which vanishes at an L1 solution. The
            // gradient of the smooth part does not: on the active set it
            // equals the L1 weight times the sign of the coordinate, so a test
            // on its norm can never fire for the case this class exists for.
            // A pseudo gradient of exactly zero is the optimality condition
            // itself, which is the ordinary situation for an L1 weight above
            // the point where the solution collapses to zero; it must not be
            // mistaken for the non-ascent direction of a wrong gradient.
            double pseudoGradientNorm = VectorOps.twoNorm(steepestDescentDirection);
            gradientNorm = pseudoGradientNorm;
            if (pseudoGradientNorm <= gradientTolerance) {
                logger.fine("Exiting OWL-BFGS: pseudo gradient norm "
                        + pseudoGradientNorm + " <= " + gradientTolerance);
                termination = Termination.GRADIENT_TOLERANCE;
                return true;
            }

            // adjust for curvature
            mapDirByInverseHessian(yDotY);

            // fix direction signs
            fixDirSigns();

            // The sign fix can cancel every coordinate of a quasi-Newton
            // direction that disagrees with the steepest descent everywhere.
            // Fall back to the steepest descent rather than step nowhere.
            if (VectorOps.twoNorm(direction) == 0.0) {
                logger.fine("Curvature-adjusted direction was cancelled; using steepest descent");
                storeSrcInDest(steepestDescentDirection, direction);
            }

            // backup parameters and gradient; then perform line-search
            storeSrcInDest(parameters, oldParameters);
            storeSrcInDest(grad, oldGrad);
            if (!backTrackingLineSearch()) {
                logger.warning("Line search in OWL-BFGS could not find a lower value in "
                        + MAX_BACKTRACKS + " backtracking steps. "
                        + "Stopping at the last accepted parameters.");
                termination = Termination.LINE_SEARCH_STALLED;
                return false;
            }

            // update gradient after line search
            evalGradient();

            // check for termination conditions
            if (checkValueTerminationCondition()) {
                // report the pseudo gradient at the point actually reached rather
                // than the one formed at the top of this iteration
                makeSteepestDescDir();
                gradientNorm = VectorOps.twoNorm(steepestDescentDirection);
//                logger.info("Exiting OWL-BFGS on termination #1:");
//                logger.info("value difference below tolerance (oldValue: "
//                        + oldValue + " newValue: " + value);
                termination = Termination.VALUE_TOLERANCE;
                return true;
            }

            // The pseudo gradient test that used to sit here is now at the top
            // of the loop, where the direction it needs has just been built.

            // update hessian approximation
            shift();

            iterations++;
            if (iterations > maxIterations) {
                // exhausting the budget is not convergence; reporting it as
                // convergence makes the two indistinguishable to the caller
                logger.warning("Too many iterations in OWL-BFGS (" + maxIterations
                        + "). Stopping with the current parameters, which are "
                        + "not known to be optimal.");
                termination = Termination.ITERATION_LIMIT;
                return false;
            }
        }

        // the loop used up the iterations this call was given, which leaves the
        // search unfinished rather than converged
        termination = Termination.PARTIAL_RUN;
        return false;
    }

    /**
     * Evaluate value. Make it a minimization problem.
     */
    private double evalL1() {
        double val = -optimizable.getValue();
        double sumAbsWt = 0;
        if (l1Weight > 0) {
            for (double param : parameters) {
                if (Double.isInfinite(param))
                    continue;
                sumAbsWt += Math.abs(param) * l1Weight;
            }
        }
//        logger.info("getValue() (" + optName + ".getValue() = " + val
//                + " + |w|=" + sumAbsWt + ") = " + (val + sumAbsWt));

        return val + sumAbsWt;
    }

    /**
     * Evaluate gradient, make it a descent direction.
     */
    private void evalGradient() {
        optimizable.getValueGradient(grad);
        adjustGradForInfiniteParams(grad);
        VectorOps.timesEquals(grad, -1.0);
    }

    /**
     * Creates steepest ascent direction from gradient and L1-regularization.
     */
    private void makeSteepestDescDir() {
        if (l1Weight == 0) {
            for (int i = 0; i < grad.length; i++) {
                direction[i] = -grad[i];
            }
        } else {
            for (int i = 0; i < grad.length; i++) {
                if (parameters[i] < 0) {
                    direction[i] = -grad[i] + l1Weight;
                } else if (parameters[i] > 0) {
                    direction[i] = -grad[i] - l1Weight;
                } else {
                    if (grad[i] < -l1Weight) {
                        direction[i] = -grad[i] - l1Weight;
                    } else if (grad[i] > l1Weight) {
                        direction[i] = -grad[i] + l1Weight;
                    } else {
                        direction[i] = 0;
                    }
                }
            }
        }

        storeSrcInDest(direction, steepestDescentDirection);
    }

    private void adjustGradForInfiniteParams(double d[]) {
        for (int i = 0; i < parameters.length; i++) {
            if (Double.isInfinite(parameters[i]))
                d[i] = 0;
        }
    }

    /**
     * Adjusts direction based on approximate hessian inverse.
     * 
     * @param yDotY
     *            y^T * y in BFGS calculation.
     */
    private void mapDirByInverseHessian(double yDotY) {
        if (s.size() == 0)
            return;

        int count = s.size();
        for (int i = count - 1; i >= 0; i--) {
            alphas[i] = -VectorOps.dotProduct(s.get(i), direction)
                    / rhos.get(i);
            VectorOps.plusEquals(direction, y.get(i), alphas[i]);
        }

        double scalar = rhos.get(count - 1) / yDotY;
        logger.fine("Direction multiplier = " + scalar);
        VectorOps.timesEquals(direction, scalar);

        for (int i = 0; i < count; i++) {
            double beta = VectorOps.dotProduct(y.get(i), direction)
                    / rhos.get(i);
            VectorOps.plusEquals(direction, s.get(i), -alphas[i] - beta);
        }
    }

    private void fixDirSigns() {
        if (l1Weight > 0) {
            for (int i = 0; i < direction.length; i++) {
                if (direction[i] * steepestDescentDirection[i] <= 0) {
                    direction[i] = 0;
                }
            }
        }
    }

    private double dirDeriv() {
        if (l1Weight == 0) {
            return VectorOps.dotProduct(direction, grad);
        } else {
            double val = 0.0;
            for (int i = 0; i < direction.length; i++) {
                if (direction[i] != 0) {
                    if (parameters[i] < 0) {
                        val += direction[i] * (grad[i] - l1Weight);
                    } else if (parameters[i] > 0) {
                        val += direction[i] * (grad[i] + l1Weight);
                    } else if (direction[i] < 0) {
                        val += direction[i] * (grad[i] - l1Weight);
                    } else if (direction[i] > 0) {
                        val += direction[i] * (grad[i] + l1Weight);
                    }
                }
            }

            return val;
        }
    }

    private void shift() {
        // The differences go into scratch buffers first. Taking the target
        // array out of the queue before rho is known would destroy the oldest
        // pair even when the update turns out to be unusable.
        double rho = 0.0;
        double yy = 0.0;
        for (int i = 0; i < parameters.length; i++) {
            if (Double.isInfinite(parameters[i])
                    && Double.isInfinite(oldParameters[i])
                    && parameters[i] * oldParameters[i] > 0) {
                sBuf[i] = 0;
            } else {
                sBuf[i] = parameters[i] - oldParameters[i];
            }

            if (Double.isInfinite(grad[i]) && Double.isInfinite(oldGrad[i])
                    && grad[i] * oldGrad[i] > 0) {
                yBuf[i] = 0;
            } else {
                yBuf[i] = grad[i] - oldGrad[i];
            }

            rho += sBuf[i] * yBuf[i];
            yy += yBuf[i] * yBuf[i];
        }

        logger.fine("rho=" + rho);
        // The curvature condition rho > 0 holds for a strictly convex
        // objective, which is the only kind this originally had to serve. It
        // does not hold for one that is merely convex -- a rank deficient
        // least squares, a duplicated column -- where a zero rho is normal
        // rather than a caller error. rho and y.y are both divisors in
        // mapDirByInverseHessian, so a zero or a non-finite value would turn
        // the search direction into an infinity or a NaN. Skip the update and
        // keep the pairs already stored.
        if (!(rho > 0.0) || !(yy > 0.0) || Double.isInfinite(rho)
                || Double.isInfinite(yy)) {
            logger.fine("Skipping the BFGS update: rho = " + rho + ", y.y = "
                    + yy + " gives no usable curvature information");
        } else {
            double[] nextS;
            double[] nextY;
            if (s.size() < m) {
                nextS = new double[parameters.length];
                nextY = new double[parameters.length];
            } else {
                nextS = s.removeFirst();
                nextY = y.removeFirst();
                rhos.removeFirst();
            }
            storeSrcInDest(sBuf, nextS);
            storeSrcInDest(yBuf, nextY);

            s.addLastNoCopy(nextS);
            y.addLastNoCopy(nextY);
            rhos.addLast(rho);
            yDotY = yy;
        }

        // update old params and grad
        storeSrcInDest(parameters, oldParameters);
        storeSrcInDest(grad, oldGrad);
    }

    private static void storeSrcInDest(double src[], double dest[]) {
        System.arraycopy(src, 0, dest, 0, src.length);
    }

    // backtrack line search; returns false if no acceptable step was found
    private boolean backTrackingLineSearch() {
        double origDirDeriv = dirDeriv();
        if (origDirDeriv >= 0) {
            throw new InvalidOptimizableException(
                    "L-BFGS chose a non-ascent direction: check your gradient!");
        }

        double alpha = 1.0;
        double backoff = 0.5;
        if (iterations == 0) {
            // the first step is normalized to unit length, so the norm
            // here must be the safe one -- the naive sum of squares
            // this used to take answers Infinity above about 8.4e152
            // per element and 0.0 below about 1.1e-162, and alpha is
            // one divided by it
            double normDir = VectorOps.twoNorm(direction);
            alpha = 1.0 / normDir;
            backoff = 0.1;
        }

        final double c1 = 1e-4;
        // store old value
        oldValue = value;

        logger.fine("*** Starting line search iter=" + iterations);
        logger.fine("iter[" + iterations + "] Value at start of line search = "
                + value);

        // For a convex objective with a correct gradient the Armijo condition
        // succeeds after a few backtracks. Without a floor, a violated
        // precondition -- or a value that is NaN -- halves alpha all the way
        // to underflow, where the test degenerates to oldValue <= oldValue and
        // passes with zero progress, which the caller then reads as
        // convergence.
        for (int backtrack = 0; backtrack < MAX_BACKTRACKS; backtrack++) {
            // update parameters and gradient
            getNextPoint(alpha);

            // find new value
            value = evalL1();

            logger.fine("iter[" + iterations + "] Using alpha = " + alpha
                    + " new value = " + value + " |grad|="
                    + VectorOps.twoNorm(grad) + " |x|="
                    + VectorOps.twoNorm(parameters));

            if (value <= oldValue + c1 * origDirDeriv * alpha) {
                return true;
            }

            alpha *= backoff;
        }

        // no acceptable step: put the optimizable back where it was
        storeSrcInDest(oldParameters, parameters);
        optimizable.setParameters(parameters);
        value = oldValue;
        return false;
    }

    private void getNextPoint(double alpha) {
        for (int i = 0; i < parameters.length; i++) {
            parameters[i] = oldParameters[i] + direction[i] * alpha;
            if (l1Weight > 0) {
                // do not allow to cross orthant boundaries if using
                // L1-regularization
                if (oldParameters[i] * parameters[i] < 0) {
                    parameters[i] = 0.0;
                }
            }
        }

        optimizable.setParameters(parameters);
    }

    // termination conditions
    private boolean checkValueTerminationCondition() {
        return (2.0 * Math.abs(value - oldValue) <= tolerance
                * (Math.abs(value) + Math.abs(oldValue) + eps));
    }

}
