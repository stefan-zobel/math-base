package math.probe;

import java.util.ArrayList;
import java.util.List;

/**
 * A compact T-Digest implementation. Estimates quantiles (like the median)
 * using centroids and linear interpolation.
 */
public final class SimpleTDigest {

    private static final long COMPRESS_THRESHOLD = 65_536L - 1L;

    private final List<Centroid> centroids = new ArrayList<>();
    // Higher => more centroids (more memory, more centroid looping (slower))
    // => more precision
    private final double compression;
    private long totalWeight = 0L;

    /**
     * Uses compression 50.0
     */
    public SimpleTDigest() {
        this(50.0);
    }

    /**
     * @param compression
     *            Accuracy factor. Higher values (e.g., 100.0) increase
     *            precision but also memory usage.
     */
    public SimpleTDigest(double compression) {
        this.compression = compression;
    }

    private static class Centroid implements Comparable<Centroid> {
        double mean;
        long weight;

        Centroid(double value, long weight) {
            this.mean = value;
            this.weight = weight;
        }

        void add(double value) {
            weight++;
            // Incremental mean calculation
            mean = mean + (value - mean) / weight;
        }

        @Override
        public int compareTo(Centroid other) {
            return Double.compare(this.mean, other.mean);
        }
    }

    /**
     * Adds a new sample to the digest.
     */
    public void accept(double x) {
        if (centroids.isEmpty()) {
            centroids.add(new Centroid(x, 1));
            totalWeight = 1L;
            return;
        }

        // 1. Find the nearest centroid
        Centroid nearest = null;
        double minDist = Double.MAX_VALUE;
        for (Centroid c : centroids) {
            double dist = Math.abs(c.mean - x);
            if (dist < minDist) {
                minDist = dist;
                nearest = c;
            }
        }

        // 2. Calculate the quantile limit for the k-size function
        double q = getApproxQuantile(nearest.mean);
        // This limit ensures more centroids at the edges (0 and 1)
        // and fewer in the middle for better accuracy on skewed data.
        double limit = 4.0 * totalWeight * q * (1.0 - q) / compression;

        if (nearest.weight + 1 <= Math.max(1, limit)) {
            nearest.add(x);
        } else {
            centroids.add(new Centroid(x, 1));
            centroids.sort(null);
        }
        totalWeight++;

        if ((totalWeight & COMPRESS_THRESHOLD) == 0L) {
            compress();
        }
    }

    private double getApproxQuantile(double mean) {
        long cumulativeWeight = 0;
        for (Centroid c : centroids) {
            if (c.mean >= mean) {
                break;
            }
            cumulativeWeight += c.weight;
        }
        return (double) (cumulativeWeight + 0.5) / totalWeight;
    }

    /**
     * Estimates the median using linear interpolation between centroids.
     */
    public double getMedian() {
        return getQuantile(0.5);
    }

    /**
     * Estimates any quantile between 0.0 and 1.0.
     */
    public double getQuantile(double q) {
        if (centroids.isEmpty()) {
            return Double.NaN;
        }
        if (centroids.size() == 1) {
            return centroids.get(0).mean;
        }

        double targetWeight = q * totalWeight;
        double currentWeight = 0.0;

        for (int i = 0; i < centroids.size(); i++) {
            Centroid c = centroids.get(i);
            double nextWeight = currentWeight + c.weight;

            if (nextWeight >= targetWeight) {
                // Special case: If it's the very first centroid and we haven't
                // reached target
                if (i == 0 && targetWeight <= c.weight / 2.0) {
                    return c.mean;
                }
                // If it's the last centroid
                if (i == centroids.size() - 1 && targetWeight >= totalWeight - c.weight / 2.0) {
                    return c.mean;
                }

                // Interpolation logic:
                // We interpolate between the centers of mass of the centroids
                Centroid left = (targetWeight < currentWeight + c.weight / 2.0) ? centroids.get(Math.max(0, i - 1)) : c;
                Centroid right = (left == c) ? centroids.get(Math.min(centroids.size() - 1, i + 1)) : c;

                if (left == right) {
                    return left.mean;
                }

                double leftWeight = (left == c) ? currentWeight + c.weight / 2.0
                        : currentWeight - centroids.get(i - 1).weight / 2.0;
                double rightWeight = (right == c) ? currentWeight + c.weight / 2.0
                        : currentWeight + c.weight + centroids.get(i + 1).weight / 2.0;

                double fraction = (targetWeight - leftWeight) / (rightWeight - leftWeight);
                return left.mean + fraction * (right.mean - left.mean);
            }
            currentWeight = nextWeight;
        }
        return centroids.get(centroids.size() - 1).mean;
    }

    /**
     * Estimates the Cumulative Distribution Function (CDF) for a given value x.
     * Returns the fraction of all samples that are less than or equal to x.
     * 
     * @param x The value to evaluate.
     * @return A value between 0.0 and 1.0.
     */
    public double getCDF(double x) {
        if (centroids.isEmpty()) return Double.NaN;
        if (x < centroids.get(0).mean) return 0.0;
        if (x >= centroids.get(centroids.size() - 1).mean) return 1.0;

        double cumulativeWeight = 0.0;
        
        for (int i = 0; i < centroids.size() - 1; i++) {
            Centroid curr = centroids.get(i);
            Centroid next = centroids.get(i + 1);
            
            // Weight contribution of the current centroid (center of mass)
            double currWeightPos = cumulativeWeight + curr.weight / 2.0;
            double nextWeightPos = cumulativeWeight + curr.weight + next.weight / 2.0;

            if (x >= curr.mean && x < next.mean) {
                // Linear interpolation between the centers of mass of two centroids
                double fraction = (x - curr.mean) / (next.mean - curr.mean);
                double resultWeight = currWeightPos + fraction * (nextWeightPos - currWeightPos);
                return resultWeight / totalWeight;
            }
            cumulativeWeight += curr.weight;
        }
        return 1.0;
    }

    public double getBowleySkewness() {
        if (totalWeight < 4L) return Double.NaN;

        double q1 = getQuantile(0.25);
        double q2 = getQuantile(0.50); // Median
        double q3 = getQuantile(0.75);

        double denominator = q3 - q1;
        return (denominator == 0.0) ? 0.0 : (q3 + q1 - 2 * q2) / denominator;
    }

    /**
     * P99/P50 Ratio
     */
    public double getTailRatio() {
        double p99 = getQuantile(0.99);
        double p50 = getQuantile(0.50);
        double p01 = getQuantile(0.01);
        return (p50 - p01 == 0.0) ? 0.0 : (p99 - p50) / (p50 - p01);
    }

    public void compress() {
        if (centroids.size() <= compression) {
            return;
        }

        centroids.sort(null);
        List<Centroid> compressed = new ArrayList<>();
        Centroid current = centroids.get(0);
        
        double cumulativeWeight = 0.0;
        for (int i = 1; i < centroids.size(); i++) {
            Centroid next = centroids.get(i);
            double q = (cumulativeWeight + current.weight + next.weight / 2.0) / totalWeight;
            double limit = 4.0 * totalWeight * q * (1.0 - q) / compression;

            // Merge if both together are still below the limit
            if (current.weight + next.weight <= Math.max(1, limit)) {
                current.mean = (current.mean * current.weight + next.mean * next.weight) 
                               / (current.weight + next.weight);
                current.weight += next.weight;
            } else {
                compressed.add(current);
                cumulativeWeight += current.weight;
                current = next;
            }
        }
        compressed.add(current);
        centroids.clear();
        centroids.addAll(compressed);
    }

    /**
     * Returns the total number of samples processed.
     */
    public long getTotalWeight() {
        return totalWeight;
    }

    /**
     * Returns the current number of centroids in the digest.
     * This represents the memory footprint of the sketch.
     */
    public int getCentroidCount() {
        return centroids.size();
    }

    public double getCompression() {
        return compression;
    }
}
