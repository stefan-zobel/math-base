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
}
