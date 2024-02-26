/*
 * Copyright (C) 2006 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

/**
 * Exception thrown by optimization algorithms, when the problem is usually due
 * to a problem with the given Maximizable instance.
 * <p>
 * If the optimizer throws this in your code, usually there are two possible
 * causes: (a) you are computing the gradients approximately, (b) your value
 * function and gradient do not match (this can be checking using
 * 
 * 
 * @author <A HREF="mailto:casutton@cs.umass.edu">casutton@cs.umass.edu</A>
 * @version $Id: InvalidOptimizableException.java,v 1.1 2007/10/22 21:37:39
 *          mccallum Exp $
 */
public class InvalidOptimizableException extends OptimizationException {

    private static final long serialVersionUID = -8463567732737341072L;

    public InvalidOptimizableException() {
    }

    public InvalidOptimizableException(String message) {
        super(message);
    }

    public InvalidOptimizableException(String message, Throwable cause) {
        super(message, cause);
    }

    public InvalidOptimizableException(Throwable cause) {
        super(cause);
    }
}
