import data.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class LocallyThreadedNBodyCalculator {

    private final List<Point> points;
    private Quadrant rootQuadrant;
    private final ExecutorService taskExecutor;
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
        for (var iteration = 0; iteration < CommonCore.iterationCount; iteration++) {
            rootQuadrant = new Quadrant(
                    new Dimension(CommonCore.maxWidth, CommonCore.maxHeight),
                    new Position()
            );

            points.forEach(point -> rootQuadrant.insert(point));

            computeForces();

            points.forEach(Point::updatePosition);
            writeOut(iteration);
        }
    }

    private void writeOut(int iteration) {
        CommonCore.write(iteration, points, CommonCore.Mode.DISTRIBUTED);
    }

    private void computeForces() {
        List<Runnable> tasks = new ArrayList<>();
        for (var index = 0; index < points.size(); index += partitionSize) {
            final var start = index;
            final var end = Math.min(index + partitionSize, points.size());
            tasks.add(() -> {
                for (var jndex = start; jndex < end; jndex++) {
                    rootQuadrant.addForceActingOn(points.get(jndex));
                }
            });
        }

        try {
            taskExecutor.invokeAll(tasks.stream()
                    .map(Executors::callable)
                    .toList());
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Force computation interrupted", e);
        }
    }


    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
        simulation.taskExecutor.shutdown();
    }
}