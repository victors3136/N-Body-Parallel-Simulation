package data;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class CommonCore {
    public static final double GravitationalConstant = 1.0e7;
    public static final double maxWidth = 1.0e6;
    public static final double maxHeight = 1.0e6;
    public static final double boundingBoxTolerance = Math.min(maxWidth, maxHeight) / 100;
    public static final double minInitialForce = 0;
    public static final double maxInitialForce = 100;
    public static final double maxBodyRadius = 100;
    public static final double timeIncrement = 2;
    public static final double angle = 1.0;
    public static final int bodyCount = 8;
    public static final int iterationCount = 500;
    public static double maxMass = 1.0e3;
    public static double maxVelocity = 1.0e4;

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
        final var filename = String.format("outputs/" + dir + "/csv/positions_iter_%d.csv", iterationIndex);
        try (final var fileWriter = new FileWriter(filename);
             final var bufferedWriter = new BufferedWriter(fileWriter);
             final var printWriter = new PrintWriter(bufferedWriter)) {
            printWriter.println("particle_id,x,y");
            for (var index = 0; index < CommonCore.bodyCount; index++) {
                final var point = points.get(index);
                printWriter.printf("%4d,%10.2f,%10.2f%n",
                        index,
                        point.position().horizontal(),
                        point.position().vertical());
            }
        } catch (IOException e) {
            System.err.println("Error writing to output file: " + e.getMessage());
        }
    }
}

