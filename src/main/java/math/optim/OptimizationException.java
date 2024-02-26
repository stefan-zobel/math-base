/*
 * Copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

/**
 * General exception thrown by optimization algorithms when there is an
 * optimization-specific problem. For example, an exception might be thrown when
 * the gradient is sufficiently large but no step is possible in that direction.
 * 
 * @author Jerod Weinman <a HREF="mailto:weinman@cs.umass.edu">weinman@cs.umass.edu</a>
 */
public class OptimizationException extends RuntimeException {

    private static final long serialVersionUID = 6376656698916305775L;

    public OptimizationException() {
        super();
    }

    public OptimizationException(String message) {
        super(message);
    }

    public OptimizationException(String message, Throwable cause) {
        super(message, cause);
    }

    public OptimizationException(Throwable cause) {
        super(cause);
    }
}
