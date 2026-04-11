package math.cern;

import org.junit.Assert;
import org.junit.Test;

public class InterpolationTest {

    private static final double EPS = 1e-12;

    @Test
    public void testInterpolateY() {
        double x0 = 10.0;
        double x1 = 20.0;
        double y0 = 100.0;
        double y1 = 200.0;

        // On the spots
        Assert.assertEquals(100.0, Arithmetic.interpolateY(x0, x1, y0, y1, 10.0), EPS);
        Assert.assertEquals(200.0, Arithmetic.interpolateY(x0, x1, y0, y1, 20.0), EPS);

        // Interpolation
        Assert.assertEquals(150.0, Arithmetic.interpolateY(x0, x1, y0, y1, 15.0), EPS);

        // Extrapolation
        // to the left (x < x0) -> 5.0
        Assert.assertEquals(50.0, Arithmetic.interpolateY(x0, x1, y0, y1, 5.0), EPS);
        // to the right (x > x1) -> 25.0
        Assert.assertEquals(250.0, Arithmetic.interpolateY(x0, x1, y0, y1, 25.0), EPS);
    }

    @Test
    public void testInterpolateX() {
        double x0 = 0.0;
        double x1 = 10.0;
        double y0 = 0.0;
        double y1 = 100.0;

        // On the spots
        Assert.assertEquals(0.0, Arithmetic.interpolateX(x0, x1, y0, y1, 0.0), EPS);
        Assert.assertEquals(10.0, Arithmetic.interpolateX(x0, x1, y0, y1, 100.0), EPS);

        // Interpolation
        Assert.assertEquals(5.0, Arithmetic.interpolateX(x0, x1, y0, y1, 50.0), EPS);

        // Extrapolation
        // below y0
        Assert.assertEquals(-1.0, Arithmetic.interpolateX(x0, x1, y0, y1, -10.0), EPS);
        // above y1
        Assert.assertEquals(11.0, Arithmetic.interpolateX(x0, x1, y0, y1, 110.0), EPS);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testInvalidRange() {
        Arithmetic.interpolateY(20.0, 10.0, 100.0, 200.0, 15.0);
    }
}
