package math.linalg;

import static org.junit.Assert.*;
import org.junit.Test;

public class VectorOpsTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testIsVectorized() {
        String versionStr = System.getProperty("java.specification.version");
        double version = Double.parseDouble(versionStr.startsWith("1.") ? versionStr.substring(2) : versionStr);
        if (version >= 25.0) {
            assertTrue("Java 25 should use the Vector API", VectorOps.isVectorized());
        } else {
            assertFalse("Java 8 must not use the Vector API", VectorOps.isVectorized());
        }
    }

    @Test
    public void testTimesEquals() {
        double[] m = {1.0, 2.0, 3.0, 4.0, 5.0};
        VectorOps.timesEquals(m, 2.5);
        
        double[] expected = {2.5, 5.0, 7.5, 10.0, 12.5};
        assertArrayEquals(expected, m, EPSILON);
    }

    @Test
    public void testDotProduct() {
        double[] m1 = {1.0, 2.0, 3.0};
        double[] m2 = {4.0, 5.0, 6.0};
        // (1*4) + (2*5) + (3*6) = 4 + 10 + 18 = 32
        double result = VectorOps.dotProduct(m1, m2);
        assertEquals(32.0, result, EPSILON);
    }

    @Test
    public void testPlusEquals() {
        double[] m1 = {
            10.0, 
            Double.POSITIVE_INFINITY, 
            Double.NEGATIVE_INFINITY, 
            Double.POSITIVE_INFINITY 
        };
        double[] m2 = {
            5.5, 
            Double.NEGATIVE_INFINITY, // Inf + -Inf -> 0.0
            Double.POSITIVE_INFINITY, // -Inf + Inf -> 0.0
            Double.POSITIVE_INFINITY  // Inf + Inf -> Inf
        };
        
        VectorOps.plusEquals(m1, m2);
        
        assertEquals(15.5, m1[0], EPSILON);
        assertEquals("Inf + -Inf muss 0.0 sein", 0.0, m1[1], EPSILON);
        assertEquals("-Inf + Inf muss 0.0 sein", 0.0, m1[2], EPSILON);
        assertEquals("Inf + Inf bleibt Inf", Double.POSITIVE_INFINITY, m1[3], 0.0);
    }

    @Test
    public void testPlusEqualsWithFactor() {
        double[] m1 = {1.0, Double.POSITIVE_INFINITY, 2.0};
        double[] m2 = {3.0, Double.NEGATIVE_INFINITY, 4.0};
        double factor = 2.0;
        
        // Operation: m1[i] += m2[i] * factor
        VectorOps.plusEquals(m1, m2, factor);
        
        assertEquals(1.0 + (3.0 * 2.0), m1[0], EPSILON); // 7.0
        assertEquals("Spezialfall Infinites mit Faktor", 0.0, m1[1], EPSILON);
        assertEquals(2.0 + (4.0 * 2.0), m1[2], EPSILON); // 10.0
    }

    @Test
    public void testTailCleanup() {
        double[] m1 = {1, 1, 1, 1, 1, 1, 1};
        double[] m2 = {2, 2, 2, 2, 2, 2, 2};
        
        VectorOps.plusEquals(m1, m2);
        
        for (double val : m1) {
            assertEquals(3.0, val, EPSILON);
        }
    }
}
