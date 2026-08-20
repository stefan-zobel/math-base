package math.lapack;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import org.junit.Test;

import math.gemm.Trans;

/**
 * Tests for {@link Xerbla}, the LAPACK error reporter.
 * <p>
 * It used to write the message to {@code System.err} as well as throwing it. A
 * library without runtime dependencies has no business printing to the console,
 * and the exception already carries the same text, so the write was pure
 * duplication.
 */
public class XerblaTest {

    /** Runs the body with stderr captured, and returns whatever it wrote. */
    private static String stderrOf(Runnable body) {
        PrintStream saved = System.err;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        System.setErr(new PrintStream(captured));
        try {
            body.run();
        } finally {
            System.setErr(saved);
        }
        return captured.toString();
    }

    @Test
    public void testItThrowsAndNamesTheOffendingArgument() {
        String printed = stderrOf(new Runnable() {
            public void run() {
                try {
                    Xerbla.xerbla("DGEQRF", 4);
                    fail("expected an IllegalArgumentException");
                } catch (IllegalArgumentException expected) {
                    assertTrue(expected.getMessage().contains("DGEQRF"));
                    assertTrue(expected.getMessage().contains("'4'"));
                }
            }
        });
        assertEquals("Xerbla must not write to stderr", "", printed);
    }

    @Test
    public void testAnInvalidArgumentReachesTheCallerQuietly() {
        // through a real routine rather than Xerbla directly: an unrecognized
        // uplo makes Dtrtrs report argument 1
        final StringBuilder message = new StringBuilder();
        String printed = stderrOf(new Runnable() {
            public void run() {
                try {
                    Dtrtrs.dtrtrs("X", "N", "N", 2, 1, new double[4], 0, 2, new double[2], 0, 2, new intW(0));
                    fail("expected an IllegalArgumentException");
                } catch (IllegalArgumentException expected) {
                    message.append(expected.getMessage());
                }
            }
        });
        assertTrue(message.toString().contains("DTRTRS"));
        assertEquals("", printed);
    }

    @Test
    public void testAValidCallPrintsNothingEither() {
        // the regression guard on the other side: a well formed call has to
        // stay silent, which it did before and has to keep doing
        String printed = stderrOf(new Runnable() {
            public void run() {
                int m = 20;
                int n = 5;
                double[] a = new double[m * n];
                long state = 4242L;
                for (int i = 0; i < a.length; ++i) {
                    state = state * 6364136223846793005L + 1442695040888963407L;
                    a[i] = 2.0 * ((state >>> 11) * 0x1.0p-53) - 1.0;
                }
                double[] b = new double[m];
                for (int i = 0; i < m; ++i) {
                    b[i] = 1.0;
                }
                double[] work = new double[4096];
                Dgels.dgels(Trans.NO_TRANS, m, n, 1, a, 0, m, b, 0, m, work, 0, work.length);
            }
        });
        assertEquals("", printed);
    }
}
