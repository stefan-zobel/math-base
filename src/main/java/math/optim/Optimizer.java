/*
 * Copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

/**
 * @author Andrew McCallum
 *         <a href="mailto:mccallum@cs.umass.edu">mccallum@cs.umass.edu</a>
 */
public interface Optimizer {
    boolean optimize();

    boolean optimize(int numIterations);

    /** Returns true if it has converged */
    boolean isConverged();

    Optimizable getOptimizable();
}
