package math.probe;

import java.util.Arrays;

/**
 * Implementation of the P-Square (P2) algorithm for dynamic calculation of the
 * median without storing all samples.
 * 
 * Based on the paper: "The P2 Algorithm for Dynamic Calculation of Quantiles
 * and Histograms without Storing Observations" by Jain and Chlamtac.
 */
public final class PSquaredMedian {

    private final double[] q = new double[5]; // Marker heights (the estimates)
    private final double[] n = new double[5]; // Marker positions (actual)
    private final double[] f = { 0.0, 0.25, 0.5, 0.75, 1.0 }; // Desired marker positions
    private final double[] np = new double[5]; // Desired marker positions (calculated)

    private int count = 0;

    /**
     * Processes a new observation from the data stream.
     *
     * @param x
     *            The new value to be added to the calculation.
     */
    public void accept(double x) {
        count++;

        // Initialization: The first 5 observations are
        // sorted and used as initial markers
        if (count <= 5) {
            q[count - 1] = x;
            if (count == 5) {
                Arrays.sort(q);
                for (int i = 0; i < 5; i++) {
                    n[i] = i + 1;
                }
                updateDesiredPositions();
            }
            return;
        }

        // 1. Find the cell k such that q[k] <= x < q[k+1]
        int k;
        if (x < q[0]) {
            q[0] = x;
            k = 0;
        } else if (x >= q[4]) {
            q[4] = x;
            k = 3;
        } else {
            k = findInterval(x);
        }

        // 2. Increment positions of markers k+1 through 4
        for (int i = k + 1; i < 5; i++) {
            n[i]++;
        }

        // 3. Update desired positions for all markers
        updateDesiredPositions();

        // 4. Adjust heights of markers 1, 2, and 3 (the internal markers)
        for (int i = 1; i <= 3; i++) {
            double d = np[i] - n[i];

            if ((d >= 1 && n[i + 1] - n[i] > 1) || (d <= -1 && n[i - 1] - n[i] < -1)) {
                int sign = (d >= 0) ? 1 : -1;
                double qNext = computeParabolic(i, sign);

                // Ensure the markers remain in non-decreasing order
                if (q[i - 1] < qNext && qNext < q[i + 1]) {
                    q[i] = qNext;
                } else {
                    // Use linear formula if parabolic estimate fails
                    q[i] = q[i] + sign * (q[i + sign] - q[i]) / (n[i + sign] - n[i]);
                }
                n[i] += sign;
            }
        }
    }

    private int findInterval(double x) {
        for (int i = 0; i < 4; i++) {
            if (x < q[i + 1]) {
                return i;
            }
        }
        return 3;
    }

    private void updateDesiredPositions() {
        for (int i = 0; i < 5; i++) {
            np[i] = 1 + (count - 1) * f[i];
        }
    }

    private double computeParabolic(int i, int d) {
        double nl = n[i - 1];
        double nCurr = n[i];
        double nr = n[i + 1];
        double ql = q[i - 1];
        double qCurr = q[i];
        double qr = q[i + 1];

        return qCurr + (double) d / (nr - nl)
                * ((nCurr - nl + d) * (qr - qCurr) / (nr - nCurr) + (nr - nCurr - d) * (qCurr - ql) / (nCurr - nl));
    }

    /**
     * Returns the current estimate of the median.
     *
     * @return Estimated median (0.5 quantile).
     */
    public double getMedian() {
        if (count == 0) {
            return Double.NaN;
        }
        if (count <= 5) {
            double[] temp = Arrays.copyOf(q, count);
            Arrays.sort(temp);
            return temp[count / 2];
        }
        return q[2]; // The 3rd marker (index 2) is the 0.5 quantile
    }

    /**
     * Returns the total number of observations processed so far.
     *
     * @return Number of samples.
     */
    public int getCount() {
        return count;
    }
}
