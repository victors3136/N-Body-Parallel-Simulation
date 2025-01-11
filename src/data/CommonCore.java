package data;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class CommonCore {
    public static final double GravitationalConstant = 6.67300e-11;
    public static final double maxWidth = 1.0e6;
    public static final double maxHeight = 1.0e6;
    public static final double minInitialForce = -10;
    public static final double maxInitialForce = 10;
    public static final double maxBodyRadius = 10;
    public static final double timeIncrement = 0.5;
    public static final double angle = 1.0;
    public static final double defaultMass = 1.0;
    public static final int bodyCount = 10;
    public static final int iterationCount = 5000;
    public static double maxMass = 100;
    public static double maxInitialVelocity = 10;

    public enum Mode {
        DISTRIBUTED, THREADED;

        String toDirName() {
            return switch (this) {
                case DISTRIBUTED -> "distributed";
                case THREADED -> "threaded";
            };
        }
    }

    public static void write(int iterationIndex, List<Point> points, Mode mode) {
        final var dir = mode.toDirName();
        final var filename = String.format("src/visualise/csv/" + dir + "/positions_iter_%d.csv", iterationIndex);
        try (final var fileWriter = new FileWriter(filename);
             final var bufferedWriter = new BufferedWriter(fileWriter);
             final var printWriter = new PrintWriter(bufferedWriter)) {
            printWriter.println("particle_id,x,y");
            for (var index = 0; index < CommonCore.bodyCount; index++) {
                final var point = points.get(index);
                printWriter.printf("%d,%.6f,%.6f%n",
                        index,
                        point.position().horizontal(),
                        point.position().vertical());
            }
        } catch (IOException e) {
            System.err.println("Error writing to output file: " + e.getMessage());
        }
    }
}

