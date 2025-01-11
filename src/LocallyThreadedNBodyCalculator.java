import data.*;

import java.util.List;
import java.util.concurrent.Executor;
import java.util.concurrent.Executors;

public class LocallyThreadedNBodyCalculator {

    private final List<Point> points;
    private Quadrant rootQuadrant;
    private Executor taskExecutor;
    private final int partitionSize;

    public LocallyThreadedNBodyCalculator() {
        final var threadCount = Runtime.getRuntime().availableProcessors();
        this.points = RandomGenerator.points(CommonCore.bodyCount);
        taskExecutor = Executors.newFixedThreadPool(threadCount);
        partitionSize = CommonCore.bodyCount / threadCount;

        System.out.println("Initialization:");
        System.out.println("Bodies Count = " + CommonCore.bodyCount);
        System.out.println("Iteration Count = " + CommonCore.iterationCount);
        System.out.println("Thread Count = " + threadCount);
        System.out.println("Particles per thread = " + partitionSize);

    }

    public void run() {
        rootQuadrant = new Quadrant(
                new Dimension(CommonCore.maxWidth, CommonCore.maxHeight),
                new Position()
        );
        points.forEach(point -> rootQuadrant.insert(point));
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
    }
}