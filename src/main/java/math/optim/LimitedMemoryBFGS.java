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

import math.MatrixOps;


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

    private boolean converged = false;
    private final Optimizable.ByGradientValue optimizable;
    private final int maxIterations = 1000;
    // xxx need a more principled stopping point
    // final double tolerance = .0001;
    private double tolerance = .0001;
    private final double gradientTolerance = .001;
    private final double eps = 1.0e-5;

    // The number of corrections used in BFGS update
    // ideally 3 <= m <= 7. Larger m means more cpu time, memory.
    private final int m = 4;

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

    public LimitedMemoryBFGS(Optimizable.ByGradientValue function) {
        optimizable = function;
        lineMaximizer = new BackTrackLineSearch(function);
    }

    public Optimizable getOptimizable() {
        return optimizable;
    }

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

    public boolean optimize() {
        return optimize(Integer.MAX_VALUE);
    }

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

            if (MatrixOps.absNormalize(direction) == 0) {
                logger.info("L-BFGS initial gradient is zero; saying converged");
                g = null;
                converged = true;
                return true;
            }

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("direction.2norm: " + MatrixOps.twoNorm(direction));
            }

            MatrixOps
                    .timesEquals(direction, 1.0 / MatrixOps.twoNorm(direction));

            // make initial jump
            if (logger.isLoggable(Level.FINE)) {
                logger.fine("before initial jump: \ndirection.2norm: "
                        + MatrixOps.twoNorm(direction) + " \ngradient.2norm: "
                        + MatrixOps.twoNorm(g) + "\nparameters.2norm: "
                        + MatrixOps.twoNorm(parameters));
            }

            step = lineMaximizer.optimize(direction, step);
            if (step == 0.0) {
                // could not step in this direction.
                // give up and say converged.
                g = null; // reset search
                step = 1.0;
                throw new OptimizationException(
                        "Line search could not step in the current direction. "
                                + "(This is not necessarily cause for alarm. Sometimes this happens close to the maximum,"
                                + " where the function may be very flat.)");

                // return false;
            }

            optimizable.getParameters(parameters);
            optimizable.getValueGradient(g);

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("after initial jump: \ndirection.2norm: "
                        + MatrixOps.twoNorm(direction) + " \ngradient.2norm: "
                        + MatrixOps.twoNorm(g));
            }
        }

        for (int iterationCount = 0; iterationCount < numIterations; iterationCount++) {
            double value = optimizable.getValue();

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("L-BFGS iteration=" + iterationCount + ", value="
                        + value + " g.twoNorm: " + MatrixOps.twoNorm(g)
                        + " oldg.twoNorm: " + MatrixOps.twoNorm(oldg));
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

            if (sy > 0) {
                throw new InvalidOptimizableException("sy = " + sy + " > 0");
            }

            double gamma = sy / yy; // scaling factor

            if (gamma > 0) {
                throw new InvalidOptimizableException("gamma = " + gamma
                        + " > 0");
            }

            rho.addLast(1.0 / sy);
            // These arrays are now the *differences* between parameters and
            // gradient.
            s.addLast(oldParameters);
            y.addLast(oldg);

            assert (s.size() == y.size()) : "s.size: " + s.size() + " y.size: "
                    + y.size();

            //
            // This next section is where we calculate the new direction
            //

            // First work backwards, from the most recent difference vectors
            for (int i = s.size() - 1; i >= 0; i--) {
                alpha[i] = rho.get(i)
                        * MatrixOps.dotProduct(s.get(i), direction);
                MatrixOps.plusEquals(direction, y.get(i), -1.0 * alpha[i]);
            }

            // Scale the direction by the ratio of s'y and y'y
            MatrixOps.timesEquals(direction, gamma);

            // Now work forwards, from the oldest to the newest difference
            // vectors
            for (int i = 0; i < y.size(); i++) {
                double beta = rho.get(i)
                        * MatrixOps.dotProduct(y.get(i), direction);
                MatrixOps.plusEquals(direction, s.get(i), alpha[i] - beta);
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
                        + MatrixOps.dotProduct(direction, g)
                        + "\ndirection.2norm: " + MatrixOps.twoNorm(direction)
                        + "\nparameters.2norm: "
                        + MatrixOps.twoNorm(parameters));
            }

            // Do a line search in the current direction
            step = lineMaximizer.optimize(direction, step);

            if (step == 0.0) { // could not step in this direction.
                g = null; // reset search
                step = 1.0;
                // xxx Temporary test; passed OK
                // TestMaximizable.testValueAndGradientInDirection (maxable,
                // direction);
                throw new OptimizationException(
                        "Line search could not step in the current direction. "
                                + "(This is not necessarily cause for alarm. Sometimes this happens close to the maximum,"
                                + " where the function may be very flat.)");
                // return false;
            }
            optimizable.getParameters(parameters);
            optimizable.getValueGradient(g);

            if (logger.isLoggable(Level.FINE)) {
                logger.fine("after linesearch: direction.2norm: "
                        + MatrixOps.twoNorm(direction));
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
            double gg = MatrixOps.twoNorm(g);
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
                System.err
                        .println("Too many iterations in L-BFGS.java. Continuing with current parameters.");
                converged = true;
                return true;
                // throw new IllegalStateException ("Too many iterations.");
            }

            // End of iteration. Call evaluator
            if (eval != null && !eval.evaluate(optimizable, iterationCount)) {
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine("Exiting L-BFGS on termination #4: evaluator returned false.");
                }
                converged = true;
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
