/*
 * Copyright (C) 2004 Univ. of Massachusetts Amherst, Computer Science Dept.
 * This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
 * http://mallet.cs.umass.edu/
 * This software is licensed under the terms of the Apache License, Version 2.0
 * or (at your option) any subsequent version.
 */
package math.optim;

import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit Test for class TestOptimizer.java
 * 
 * @author <a href="mailto:casutton@cs.umass.edu">Charles Sutton</a>
 */
public class TestOptimizer extends TestCase {

    public TestOptimizer(String name) {
        super(name);
    }

    // Maximizable for -3x^2 + 5x - 2
    static class SimplePoly implements Optimizable.ByGradientValue {

        double[] params = new double[1];

        public void getParameters(double[] doubleArray) {
            doubleArray[0] = params[0];
        }

        public int getNumParameters() {
            return 1;
        }

        public double getParameter(int n) {
            return params[0];
        };

        public void setParameters(double[] doubleArray) {
            params[0] = doubleArray[0];
        }

        public void setParameter(int n, double d) {
            params[n] = d;
        }

        public double getValue() {
            System.out.println("param = " + params[0] + " value = " + (-3 * params[0] * params[0] + 5 * params[0] - 2));

            return -3 * params[0] * params[0] + 5 * params[0] - 2;
        }

        public void getValueGradient(double[] buffer) {
            buffer[0] = -6 * params[0] + 5;
        }
    }

    public void testLinearLBFGS() {
        SimplePoly poly = new SimplePoly();
        LimitedMemoryBFGS bfgs = new LimitedMemoryBFGS(poly);
        bfgs.optimize();
        assertEquals(5.0 / 6.0, poly.params[0], 1e-16);
    }

    public void testOrthantWiseLBFGSWithoutL1() {
        SimplePoly poly = new SimplePoly();
        OrthantWiseLimitedMemoryBFGS bfgs = new OrthantWiseLimitedMemoryBFGS(poly);
        bfgs.optimize();
        assertEquals(5.0 / 6.0, poly.params[0], 1e-16);
    }

    public void testOrthantWiseLBFGSWithL1() {
        SimplePoly poly = new SimplePoly();
        OrthantWiseLimitedMemoryBFGS bfgs = new OrthantWiseLimitedMemoryBFGS(poly, 3.0);
        bfgs.optimize();
        assertEquals(2.0 / 6.0, poly.params[0], 1e-16);
    }

    /**
     * @return a <code>TestSuite</code>
     */
    public static TestSuite suite() {
        return new TestSuite(TestOptimizer.class);
    }

    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }
}
