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
 *         <a HREF="mailto:mccallum@cs.umass.edu">mccallum@cs.umass.edu</a>
 */
public interface Optimizable {

    int getNumParameters();

    void getParameters(double[] buffer);

    double getParameter(int index);

    void setParameters(double[] params);

    void setParameter(int index, double value);

    public interface ByGradientValue extends Optimizable {
        void getValueGradient(double[] buffer);

        double getValue();
    }
}
