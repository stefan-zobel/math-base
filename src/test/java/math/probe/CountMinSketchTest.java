package math.probe;

import org.junit.Test;

import static org.junit.Assert.*;

import java.util.List;

public class CountMinSketchTest {

    @Test
    public void testHeavyHitterDetection() {
        // 5 rows for accuracy, 2000 cols to minimize collisions
        CountMinSketch<String> cms = new CountMinSketch<>(5, 2000);

        // Simulate 100,000 normal users with 1 request each
        for (int i = 0; i < 100_000; i++) {
            cms.add("user_" + i);
        }

        // Simulate 1 "Heavy Hitter" (e.g., a Bot) with 5,000 requests
        String attackerIp = "192.168.1.100";
        for (int i = 0; i < 5_000; i++) {
            cms.add(attackerIp);
        }

        // Estimation
        long estimatedBotCount = cms.estimateCount(attackerIp);
        long estimatedNormalUser = cms.estimateCount("user_500");

        System.out.println("Estimated Bot Requests: " + estimatedBotCount);
        System.out.println("Estimated Normal User: " + estimatedNormalUser);

        // CMS always overestimates, but for the bot it should be very close to 5000
        assertTrue(estimatedBotCount >= 5000 && estimatedBotCount < 5100);
        // Normal users might have collisions, but should be much lower
        assertTrue(estimatedNormalUser < 100);

        List<String> topKList = cms.getTopK();
        System.out.println("\nTop K: " + topKList + "\n");
        for (String user : topKList) {
            System.out.println(user + " -> " + cms.estimateCount(user));
        }
    }

    /**
     * Two sketches built with the same seed have to answer identically. Without
     * a seed they do not: the hash functions are drawn afresh for every
     * instance, so two runs of the same program report different estimates,
     * which is right for an adversarial setting and unusable in a pipeline that
     * has to reproduce its output.
     */
    @Test
    public void testTheSameSeedGivesTheSameAnswers() {
        CountMinSketch<String> first = new CountMinSketch<>(4, 128, 10, 20260822L);
        CountMinSketch<String> second = new CountMinSketch<>(4, 128, 10, 20260822L);
        CountMinSketch<String> other = new CountMinSketch<>(4, 128, 10, 20260823L);

        long lcg = 7L;
        for (int i = 0; i < 20_000; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            String key = "k" + (int) (((lcg >>> 11) * 0x1.0p-53) * 500.0);
            first.add(key);
            second.add(key);
            other.add(key);
        }

        boolean anyDifference = false;
        for (int i = 0; i < 500; i++) {
            String key = "k" + i;
            assertEquals("key " + key, first.estimateCount(key), second.estimateCount(key));
            anyDifference |= first.estimateCount(key) != other.estimateCount(key);
        }
        assertEquals(first.getTotalCount(), second.getTotalCount());
        assertEquals(first.getTopK(), second.getTopK());
        assertTrue("a different seed must give at least one different estimate", anyDifference);
    }

    /** Whatever the seed, the sketch may overcount and may never undercount. */
    @Test
    public void testItNeverUndercounts() {
        CountMinSketch<String> cms = new CountMinSketch<>(3, 64, 5, 99L);
        java.util.HashMap<String, Integer> exact = new java.util.HashMap<>();

        long lcg = 11L;
        for (int i = 0; i < 5_000; i++) {
            lcg = lcg * 6364136223846793005L + 1442695040888963407L;
            String key = "k" + (int) (((lcg >>> 11) * 0x1.0p-53) * 300.0);
            cms.add(key);
            Integer had = exact.get(key);
            exact.put(key, Integer.valueOf(had == null ? 1 : had.intValue() + 1));
        }

        long worstOver = 0L;
        for (java.util.Map.Entry<String, Integer> entry : exact.entrySet()) {
            long estimate = cms.estimateCount(entry.getKey());
            assertTrue(entry.getKey() + " was undercounted", estimate >= entry.getValue().intValue());
            worstOver = Math.max(worstOver, estimate - entry.getValue().intValue());
        }
        assertTrue("a narrow sketch on 300 keys has to collide somewhere", worstOver > 0L);
        assertEquals(5_000L, cms.getTotalCount());
    }
}
