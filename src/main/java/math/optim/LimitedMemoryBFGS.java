/*
 * Copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

import java.util.logging.Level;
import java.util.logging.Logger;

import math.linalg.VectorOps;


/**
 * Limited Memory BFGS, as described in Byrd, Nocedal, and Schnabel,
 * "Representations of Quasi-Newton Matrices and Their Use in Limited Memory
 * Methods"
 * 
 * @author Aron Culotta <a HREF="mailto:culotta@cs.umass.edu">culotta@cs.umass.edu</a>
 */
public final class LimitedMemoryBFGS implements Optimizer {

    private static final Logger logger = Logger
            .getLogger(LimitedMemoryBFGS.class.getName());

    /** Iteration budget if none is given. */
    private static final int DEFAULT_MAX_ITERATIONS = 1000;
    /** Relative change of the objective value that stops the search. */
    private static final double DEFAULT_TOLERANCE = 1.0e-4;
    /** Gradient norm that stops the search. */
    private static final double DEFAULT_GRADIENT_TOLERANCE = 1.0e-3;
    /** Number of corrections kept for the inverse Hessian. */
    private static final int DEFAULT_M = 4;

    private boolean converged = false;
    private final Optimizable.ByGradientValue optimizable;
    private final int maxIterations;
    // xxx need a more principled stopping point
    private double tolerance;
    private final double gradientTolerance;
    private final double eps = 1.0e-5;

    // The number of corrections used in BFGS update
    // ideally 3 <= m <= 7. Larger m means more cpu time, memory.
    private final int m;

    // Line search function
    private final BackTrackLineSearch lineMaximizer;

    // State of search
    // g = gradient
    // s = list of m previous "parameters" values
    // y = list of m previous "g" values
    // rho = intermediate calculation
    private double[] g, oldg, direction, parameters, oldParameters;
    private SupersedingDoubleArrayQueue s;
    private SupersedingDoubleArrayQueue y;
    private SupersedingDoubleQueue rho;
    private double[] alpha;
    private double step = 1.0;
    private int iterations;

    private OptimizerEvaluator.ByGradient eval = null;

    /**
     * Creates an optimizer with the default stopping rules.
     *
     * @param function
     *            the function to maximize
     */
    public LimitedMemoryBFGS(Optimizable.ByGradientValue function) {
        this(function, DEFAULT_MAX_ITERATIONS, DEFAULT_TOLERANCE,
                DEFAULT_GRADIENT_TOLERANCE, DEFAULT_M);
    }

    /**
     * Creates an optimizer with explicit stopping rules.
     *
     * @param function
     *            the function to maximize
     * @param maxIterations
     *            iteration budget, {@code 1} or greater; exhausting it stops
     *            the search without reporting convergence
     * @param tolerance
     *            relative change of the objective value below which the search
     *            stops, greater than {@code 0}; can still be changed later
     *            through {@link #setTolerance(double)}
     * @param gradientTolerance
     *            gradient norm below which the search stops, {@code 0} or
     *            greater. Loosening it stops the search earlier, but
     *            tightening it beyond what the line search can resolve buys
     *            nothing: the accuracy this class reaches is bounded by the
     *            step tolerances of its internal {@code BackTrackLineSearch}
     *            ({@code 1e-7} relative, {@code 1e-4} absolute), which no
     *            caller can reach.
     * @param m
     *            number of corrections kept for the inverse Hessian, between
     *            {@code 1} and {@code 100}; ideally between {@code 3} and
     *            {@code 7}
     * @since 1.5.2
     */
    public LimitedMemoryBFGS(Optimizable.ByGradientValue function,
            int maxIterations, double tolerance, double gradientTolerance,
            int m) {
        if (function == null) {
            throw new IllegalArgumentException("function is null");
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
        optimizable = function;
        lineMaximizer = new BackTrackLineSearch(function);
    }

    @Override
    public Optimizable getOptimizable() {
        return optimizable;
    }

    @Override
    public boolean isConverged() {
        return converged;
    }

    public void setTolerance(double newtol) {
        tolerance = newtol;
    }

    public void setEvaluator(OptimizerEvaluator.ByGradient eval) {
        this.eval = eval;
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

        double initialValue = optimizable.getValue();
        if (logger.isLoggable(Level.FINE)) {
            logger.fine("Entering L-BFGS.optimize(). Initial Value="
                    + initialValue);
        }

        if (g == null) { // first time through

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("First time through L-BFGS");
            }

            iterations = 0;
            s = new SupersedingDoubleArrayQueue(m);
            y = new SupersedingDoubleArrayQueue(m);
            rho = new SupersedingDoubleQueue(m);
            alpha = new double[m];

            for (int i = 0; i < m; i++) {
                alpha[i] = 0.0;
            }

            parameters = new double[optimizable.getNumParameters()];
            oldParameters = new double[optimizable.getNumParameters()];
            g = new double[optimizable.getNumParameters()];
            oldg = new double[optimizable.getNumParameters()];
            direction = new double[optimizable.getNumParameters()];

            optimizable.getParameters(parameters);
            System.arraycopy(parameters, 0, oldParameters, 0, parameters.length);

            optimizable.getValueGradient(g);
            System.arraycopy(g, 0, oldg, 0, g.length);
            System.arraycopy(g, 0, direction, 0, g.length);

            if (VectorOps.absNormalize(direction) == 0) {
                logger.info("L-BFGS initial gradient is zero; saying converged");
                g = null;
                converged = true;
                return true;
            }

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("direction.2norm: " + VectorOps.twoNorm(direction));
            }

            VectorOps
                    .timesEquals(direction, 1.0 / VectorOps.twoNorm(direction));

            // make initial jump
            if (logger.isLoggable(Level.FINE)) {
                logger.fine("before initial jump: \ndirection.2norm: "
                        + VectorOps.twoNorm(direction) + " \ngradient.2norm: "
                        + VectorOps.twoNorm(g) + "\nparameters.2norm: "
                        + VectorOps.twoNorm(parameters));
            }

            step = lineMaximizer.optimize(direction, step);
            if (step == 0.0) {
                // Could not step in this direction. This is the flat region
                // near the maximum as often as it is a real failure, and the
                // line search has already put the parameters back on the last
                // good point, so report it rather than discard the result.
                g = null; // reset search
                step = 1.0;
                logger.warning("L-BFGS could not step in the current direction "
                        + "on the initial jump. Stopping at the current parameters, "
                        + "which are not known to be optimal.");
                converged = false;
                return false;
            }

            optimizable.getParameters(parameters);
            optimizable.getValueGradient(g);

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("after initial jump: \ndirection.2norm: "
                        + VectorOps.twoNorm(direction) + " \ngradient.2norm: "
                        + VectorOps.twoNorm(g));
            }
        }

        for (int iterationCount = 0; iterationCount < numIterations; iterationCount++) {
            double value = optimizable.getValue();

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("L-BFGS iteration=" + iterationCount + ", value="
                        + value + " g.twoNorm: " + VectorOps.twoNorm(g)
                        + " oldg.twoNorm: " + VectorOps.twoNorm(oldg));
            }

            // get difference between previous 2 gradients and parameters
            double sy = 0.0;
            double yy = 0.0;

            for (int i = 0; i < oldParameters.length; i++) {
                // -inf - (-inf) = 0; inf - inf = 0
                if (Double.isInfinite(parameters[i])
                        && Double.isInfinite(oldParameters[i])
                        && (parameters[i] * oldParameters[i] > 0)) {
                    oldParameters[i] = 0.0;
                } else {
                    oldParameters[i] = parameters[i] - oldParameters[i];
                }

                if (Double.isInfinite(g[i]) && Double.isInfinite(oldg[i])
                        && (g[i] * oldg[i] > 0)) {
                    oldg[i] = 0.0;
                } else {
                    oldg[i] = g[i] - oldg[i];
                }

                sy += oldParameters[i] * oldg[i]; // si * yi
                yy += oldg[i] * oldg[i];
                direction[i] = g[i];
            }

            // This class maximizes, so the curvature condition it needs is
            // sy < 0. It holds for a strictly concave objective, which is the
            // only kind this originally had to serve, and fails routinely on
            // one that is merely concave or not concave at all -- the
            // Rosenbrock function from its classic starting point violates it
            // in the first iterations. sy is a divisor below and yy is one in
            // gamma, so zero or non-finite values would poison the direction
            // with an infinity or a NaN. Skip the update in all these cases
            // and fall back to the plain gradient direction, keeping the pairs
            // already stored.
            double gamma = -1.0;
            if (sy < 0.0 && yy > 0.0 && !Double.isInfinite(sy)
                    && !Double.isInfinite(yy)) {
                gamma = sy / yy; // scaling factor

                rho.addLast(1.0 / sy);
                // These arrays are now the *differences* between parameters
                // and gradient.
                s.addLast(oldParameters);
                y.addLast(oldg);

                assert (s.size() == y.size()) : "s.size: " + s.size()
                        + " y.size: " + y.size();
            } else if (logger.isLoggable(Level.FINE)) {
                logger.fine("Skipping the BFGS update: sy = " + sy + ", yy = "
                        + yy + " gives no usable curvature information");
            }

            //
            // This next section is where we calculate the new direction
            //

            // First work backwards, from the most recent difference vectors
            for (int i = s.size() - 1; i >= 0; i--) {
                alpha[i] = rho.get(i)
                        * VectorOps.dotProduct(s.get(i), direction);
                VectorOps.plusEquals(direction, y.get(i), -1.0 * alpha[i]);
            }

            // Scale the direction by the ratio of s'y and y'y
            VectorOps.timesEquals(direction, gamma);

            // Now work forwards, from the oldest to the newest difference
            // vectors
            for (int i = 0; i < y.size(); i++) {
                double beta = rho.get(i)
                        * VectorOps.dotProduct(y.get(i), direction);
                VectorOps.plusEquals(direction, s.get(i), alpha[i] - beta);
            }

            // Move the current values to the "last iteration" buffers and
            // negate the search direction
            for (int i = 0; i < oldg.length; i++) {
                oldParameters[i] = parameters[i];
                oldg[i] = g[i];
                direction[i] *= -1.0;
            }

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("before linesearch: direction.gradient.dotprod: "
                        + VectorOps.dotProduct(direction, g)
                        + "\ndirection.2norm: " + VectorOps.twoNorm(direction)
                        + "\nparameters.2norm: "
                        + VectorOps.twoNorm(parameters));
            }

            // Do a line search in the current direction
            step = lineMaximizer.optimize(direction, step);

            if (step == 0.0) { // could not step in this direction.
                g = null; // reset search
                step = 1.0;
                logger.warning("L-BFGS could not step in the current direction "
                        + "after " + iterations + " iterations. Stopping at the "
                        + "current parameters, which are not known to be optimal.");
                converged = false;
                return false;
            }
            optimizable.getParameters(parameters);
            optimizable.getValueGradient(g);

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("after linesearch: direction.2norm: "
                        + VectorOps.twoNorm(direction));
            }

            double newValue = optimizable.getValue();

            // Test for terminations
            if (2.0 * Math.abs(newValue - value) <= tolerance
                    * (Math.abs(newValue) + Math.abs(value) + eps)) {
                logger.info("Exiting L-BFGS on termination #1:\nvalue difference below tolerance (oldValue: "
                        + value + " newValue: " + newValue);
                converged = true;
                return true;
            }
            double gg = VectorOps.twoNorm(g);
            if (gg < gradientTolerance) {
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine("Exiting L-BFGS on termination #2: \ngradient="
                            + gg + " < " + gradientTolerance);
                }
                converged = true;
                return true;
            }
            if (gg == 0.0) {
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine("Exiting L-BFGS on termination #3: \ngradient==0.0");
                }
                converged = true;
                return true;
            }
            if (logger.isLoggable(Level.FINE)) {
                logger.fine("Gradient = " + gg);
            }

            iterations++;
            if (iterations > maxIterations) {
                // exhausting the budget is not convergence; reporting it as
                // convergence makes the two indistinguishable to the caller
                logger.warning("Too many iterations in L-BFGS (" + maxIterations
                        + "). Stopping with the current parameters, which are "
                        + "not known to be optimal.");
                converged = false;
                return false;
            }

            // End of iteration. Call evaluator
            if (eval != null && !eval.evaluate(optimizable, iterationCount)) {
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine("Exiting L-BFGS on termination #4: evaluator returned false.");
                }
                // the evaluator aborted the search; that is not convergence,
                // and the returned false already said so
                converged = false;
                return false;
            }
        }
        return false;
    }

    /**
     * Resets the previous gradients and values that are used to approximate the
     * Hessian. NOTE - If the {@link Optimizable} object is modified externally,
     * this method should be called to avoid IllegalStateExceptions.
     */
    public void reset() {
        g = null;
    }
}
